#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

PILOT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_transition_density_pilot"
)

CUBE_PATH = PILOT_ROOT / "inputs/PYR2/pilot.cistp02.cube"

DIPOLE_TABLE = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_transition_dipole_geometry_audit/"
    "transition_dipole_observations.csv"
)

OUTPUT_ROOT = PILOT_ROOT / "pilot_validation"
METRICS_CSV = OUTPUT_ROOT / "PYR2_S2_TRANSITION_DENSITY_CUBE_METRICS_DAY019.csv"
REPORT_MD = OUTPUT_ROOT / "PYR2_S2_TRANSITION_DENSITY_CUBE_VALIDATION_DAY019.md"

BOHR_TO_ANGSTROM = 0.529177210903

NET_CHARGE_ABS_TARGET = 1.0e-3
DIPOLE_COSINE_TARGET = 0.99
DIPOLE_RELATIVE_ERROR_TARGET = 0.05
BOUNDARY_MAX_RATIO_TARGET = 1.0e-3


def log(message: str = "") -> None:
    print(message, flush=True)


def parse_cube(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise SystemExit(f"Missing cube file: {path}")

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        comment_1 = handle.readline().rstrip("\n")
        comment_2 = handle.readline().rstrip("\n")

        atom_line = handle.readline().split()
        if len(atom_line) < 4:
            raise RuntimeError("Invalid cube atom/origin line.")

        natoms_raw = int(atom_line[0])
        natoms = abs(natoms_raw)
        origin = np.array(
            [float(atom_line[1]), float(atom_line[2]), float(atom_line[3])],
            dtype=np.float64,
        )

        counts: list[int] = []
        axes: list[np.ndarray] = []

        for _ in range(3):
            tokens = handle.readline().split()
            if len(tokens) < 4:
                raise RuntimeError("Invalid cube grid-vector line.")

            counts.append(abs(int(tokens[0])))
            axes.append(
                np.array(
                    [float(tokens[1]), float(tokens[2]), float(tokens[3])],
                    dtype=np.float64,
                )
            )

        symbols_atomic_numbers: list[int] = []
        atom_coordinates: list[np.ndarray] = []

        for _ in range(natoms):
            tokens = handle.readline().split()
            if len(tokens) < 5:
                raise RuntimeError("Invalid cube atom record.")

            symbols_atomic_numbers.append(int(float(tokens[0])))
            atom_coordinates.append(
                np.array(
                    [float(tokens[2]), float(tokens[3]), float(tokens[4])],
                    dtype=np.float64,
                )
            )

        # Negative NATOMS can introduce a dataset-ID record. The current
        # ORCA transition-density cube uses positive NATOMS, but support both.
        if natoms_raw < 0:
            dataset_tokens = handle.readline().split()
            if not dataset_tokens:
                raise RuntimeError("Missing cube dataset-ID record.")

        values: list[float] = []
        for line in handle:
            for token in line.split():
                values.append(float(token))

    nx, ny, nz = counts
    expected_values = nx * ny * nz

    if len(values) != expected_values:
        raise RuntimeError(
            f"Cube value-count mismatch: expected {expected_values}, "
            f"found {len(values)}."
        )

    density = np.asarray(values, dtype=np.float64).reshape(
        (nx, ny, nz),
        order="C",
    )

    axes_matrix = np.vstack(axes)
    voxel_volume = abs(float(np.linalg.det(axes_matrix)))

    if voxel_volume <= 0.0:
        raise RuntimeError("Nonpositive cube voxel volume.")

    atom_coordinates_array = np.vstack(atom_coordinates)
    atomic_numbers = np.asarray(symbols_atomic_numbers, dtype=np.int64)

    return {
        "comment_1": comment_1,
        "comment_2": comment_2,
        "natoms_raw": natoms_raw,
        "natoms": natoms,
        "origin_bohr": origin,
        "counts": np.asarray(counts, dtype=np.int64),
        "axes_bohr": axes_matrix,
        "voxel_volume_bohr3": voxel_volume,
        "density": density,
        "atomic_numbers": atomic_numbers,
        "atom_coordinates_bohr": atom_coordinates_array,
    }


def read_orca_reference() -> dict[str, float]:
    if not DIPOLE_TABLE.is_file():
        raise SystemExit(f"Missing transition-dipole table: {DIPOLE_TABLE}")

    with DIPOLE_TABLE.open(newline="", encoding="utf-8") as handle:
        matches = [
            row
            for row in csv.DictReader(handle)
            if int(row["frame"]) == 0
            and row["site"] == "PYR2"
            and int(row["root"]) == 2
            and row["family"] == "bright_like"
        ]

    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one PYR2 frame000 S2 bright row, "
            f"found {len(matches)}."
        )

    row = matches[0]

    return {
        "DX_au": float(row["DX_au"]),
        "DY_au": float(row["DY_au"]),
        "DZ_au": float(row["DZ_au"]),
        "dipole_magnitude_au": float(row["dipole_magnitude_au"]),
        "fosc": float(row["fosc"]),
        "energy_eV": float(row["state_energy_eV"]),
    }


def vector_cosine(first: np.ndarray, second: np.ndarray) -> float:
    norm_first = float(np.linalg.norm(first))
    norm_second = float(np.linalg.norm(second))

    if norm_first == 0.0 or norm_second == 0.0:
        return float("nan")

    return float(np.dot(first, second) / (norm_first * norm_second))


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    cube = parse_cube(CUBE_PATH)
    reference = read_orca_reference()

    density = np.asarray(cube["density"], dtype=np.float64)
    counts = np.asarray(cube["counts"], dtype=np.int64)
    origin = np.asarray(cube["origin_bohr"], dtype=np.float64)
    axes = np.asarray(cube["axes_bohr"], dtype=np.float64)
    dvol = float(cube["voxel_volume_bohr3"])

    nx, ny, nz = [int(value) for value in counts]

    i, j, k = np.indices((nx, ny, nz), dtype=np.float64)

    x = (
        origin[0]
        + i * axes[0, 0]
        + j * axes[1, 0]
        + k * axes[2, 0]
    )
    y = (
        origin[1]
        + i * axes[0, 1]
        + j * axes[1, 1]
        + k * axes[2, 1]
    )
    z = (
        origin[2]
        + i * axes[0, 2]
        + j * axes[1, 2]
        + k * axes[2, 2]
    )

    net_transition_charge = float(np.sum(density) * dvol)
    absolute_density_integral = float(np.sum(np.abs(density)) * dvol)
    squared_density_integral = float(np.sum(density * density) * dvol)

    atom_coordinates = np.asarray(
        cube["atom_coordinates_bohr"],
        dtype=np.float64,
    )
    atomic_numbers = np.asarray(cube["atomic_numbers"], dtype=np.int64)
    carbon_mask = atomic_numbers == 6

    if int(np.count_nonzero(carbon_mask)) != 16:
        raise RuntimeError(
            f"Expected 16 carbon atoms, found "
            f"{int(np.count_nonzero(carbon_mask))}."
        )

    carbon_centroid = np.mean(atom_coordinates[carbon_mask], axis=0)

    first_moment_global = np.array(
        [
            float(np.sum(density * x) * dvol),
            float(np.sum(density * y) * dvol),
            float(np.sum(density * z) * dvol),
        ],
        dtype=np.float64,
    )

    first_moment_centered = np.array(
        [
            float(np.sum(density * (x - carbon_centroid[0])) * dvol),
            float(np.sum(density * (y - carbon_centroid[1])) * dvol),
            float(np.sum(density * (z - carbon_centroid[2])) * dvol),
        ],
        dtype=np.float64,
    )

    reference_vector = np.array(
        [
            reference["DX_au"],
            reference["DY_au"],
            reference["DZ_au"],
        ],
        dtype=np.float64,
    )

    candidates = {
        "+global_first_moment": first_moment_global,
        "-global_first_moment": -first_moment_global,
        "+centered_first_moment": first_moment_centered,
        "-centered_first_moment": -first_moment_centered,
    }

    comparison_rows: list[dict[str, object]] = []

    for convention, vector in candidates.items():
        magnitude = float(np.linalg.norm(vector))
        cosine = vector_cosine(vector, reference_vector)
        relative_error = abs(
            magnitude - float(reference["dipole_magnitude_au"])
        ) / float(reference["dipole_magnitude_au"])

        comparison_rows.append(
            {
                "convention": convention,
                "mu_x_au": vector[0],
                "mu_y_au": vector[1],
                "mu_z_au": vector[2],
                "magnitude_au": magnitude,
                "cosine_to_ORCA": cosine,
                "absolute_cosine_to_ORCA": abs(cosine),
                "relative_magnitude_error": relative_error,
            }
        )

    best = max(
        comparison_rows,
        key=lambda row: (
            float(row["absolute_cosine_to_ORCA"]),
            -float(row["relative_magnitude_error"]),
        ),
    )

    boundary_mask = np.zeros_like(density, dtype=bool)
    boundary_mask[0, :, :] = True
    boundary_mask[-1, :, :] = True
    boundary_mask[:, 0, :] = True
    boundary_mask[:, -1, :] = True
    boundary_mask[:, :, 0] = True
    boundary_mask[:, :, -1] = True

    max_abs_density = float(np.max(np.abs(density)))
    max_abs_boundary_density = float(
        np.max(np.abs(density[boundary_mask]))
    )
    boundary_max_ratio = (
        max_abs_boundary_density / max_abs_density
        if max_abs_density > 0.0
        else float("nan")
    )

    boundary_abs_integral = float(
        np.sum(np.abs(density[boundary_mask])) * dvol
    )
    boundary_abs_integral_fraction = (
        boundary_abs_integral / absolute_density_integral
        if absolute_density_integral > 0.0
        else float("nan")
    )

    grid_lengths_bohr = np.array(
        [
            float(np.linalg.norm(axes[0])) * (nx - 1),
            float(np.linalg.norm(axes[1])) * (ny - 1),
            float(np.linalg.norm(axes[2])) * (nz - 1),
        ],
        dtype=np.float64,
    )

    net_charge_pass = abs(net_transition_charge) <= NET_CHARGE_ABS_TARGET
    dipole_cosine_pass = (
        float(best["absolute_cosine_to_ORCA"])
        >= DIPOLE_COSINE_TARGET
    )
    dipole_magnitude_pass = (
        float(best["relative_magnitude_error"])
        <= DIPOLE_RELATIVE_ERROR_TARGET
    )
    boundary_pass = boundary_max_ratio <= BOUNDARY_MAX_RATIO_TARGET

    overall_pass = (
        net_charge_pass
        and dipole_cosine_pass
        and dipole_magnitude_pass
        and boundary_pass
    )

    metrics_row = {
        "cube_file": str(CUBE_PATH.relative_to(PROJECT_ROOT)),
        "natoms": int(cube["natoms"]),
        "nx": nx,
        "ny": ny,
        "nz": nz,
        "n_values": int(density.size),
        "voxel_volume_bohr3": dvol,
        "grid_length_x_bohr": grid_lengths_bohr[0],
        "grid_length_y_bohr": grid_lengths_bohr[1],
        "grid_length_z_bohr": grid_lengths_bohr[2],
        "grid_length_x_A": grid_lengths_bohr[0] * BOHR_TO_ANGSTROM,
        "grid_length_y_A": grid_lengths_bohr[1] * BOHR_TO_ANGSTROM,
        "grid_length_z_A": grid_lengths_bohr[2] * BOHR_TO_ANGSTROM,
        "net_transition_charge": net_transition_charge,
        "absolute_density_integral": absolute_density_integral,
        "squared_density_integral": squared_density_integral,
        "maximum_absolute_density": max_abs_density,
        "maximum_absolute_boundary_density": max_abs_boundary_density,
        "boundary_max_ratio": boundary_max_ratio,
        "boundary_absolute_integral_fraction": (
            boundary_abs_integral_fraction
        ),
        "ORCA_mu_x_au": reference_vector[0],
        "ORCA_mu_y_au": reference_vector[1],
        "ORCA_mu_z_au": reference_vector[2],
        "ORCA_mu_magnitude_au": reference["dipole_magnitude_au"],
        "best_cube_dipole_convention": best["convention"],
        "best_cube_mu_x_au": best["mu_x_au"],
        "best_cube_mu_y_au": best["mu_y_au"],
        "best_cube_mu_z_au": best["mu_z_au"],
        "best_cube_mu_magnitude_au": best["magnitude_au"],
        "best_absolute_cosine_to_ORCA": (
            best["absolute_cosine_to_ORCA"]
        ),
        "best_relative_magnitude_error": (
            best["relative_magnitude_error"]
        ),
        "net_charge_pass": net_charge_pass,
        "dipole_cosine_pass": dipole_cosine_pass,
        "dipole_magnitude_pass": dipole_magnitude_pass,
        "boundary_pass": boundary_pass,
        "overall_pilot_pass": overall_pass,
    }

    with METRICS_CSV.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(metrics_row.keys()),
        )
        writer.writeheader()
        writer.writerow(metrics_row)

    with REPORT_MD.open("w", encoding="utf-8") as handle:
        handle.write(
            "# Day019 PYR2 S2 transition-density cube validation\n\n"
        )

        handle.write("## Cube structure\n\n")
        handle.write(f"- Atoms: {int(cube['natoms'])}\n")
        handle.write(f"- Grid: {nx} Ã {ny} Ã {nz}\n")
        handle.write(f"- Values: {density.size}\n")
        handle.write(
            f"- Voxel volume: {dvol:.10e} bohrÂ³\n"
        )
        handle.write(
            "- Grid lengths: "
            f"{grid_lengths_bohr[0] * BOHR_TO_ANGSTROM:.6f} Ã "
            f"{grid_lengths_bohr[1] * BOHR_TO_ANGSTROM:.6f} Ã "
            f"{grid_lengths_bohr[2] * BOHR_TO_ANGSTROM:.6f} Ã\n\n"
        )

        handle.write("## Integral diagnostics\n\n")
        handle.write(
            f"- Net transition charge: {net_transition_charge:.10e}\n"
        )
        handle.write(
            f"- Integral of |rho_tr|: "
            f"{absolute_density_integral:.10e}\n"
        )
        handle.write(
            f"- Integral of rho_trÂ²: "
            f"{squared_density_integral:.10e}\n"
        )
        handle.write(
            f"- Boundary maximum / global maximum: "
            f"{boundary_max_ratio:.10e}\n"
        )
        handle.write(
            f"- Boundary |rho| integral fraction: "
            f"{boundary_abs_integral_fraction:.10e}\n\n"
        )

        handle.write("## Dipole reconstruction\n\n")
        handle.write(
            "- ORCA transition dipole: "
            f"({reference_vector[0]:.10f}, "
            f"{reference_vector[1]:.10f}, "
            f"{reference_vector[2]:.10f}) au\n"
        )
        handle.write(
            f"- ORCA magnitude: "
            f"{float(reference['dipole_magnitude_au']):.10f} au\n"
        )
        handle.write(
            f"- Best cube convention: `{best['convention']}`\n"
        )
        handle.write(
            "- Best cube dipole: "
            f"({float(best['mu_x_au']):.10f}, "
            f"{float(best['mu_y_au']):.10f}, "
            f"{float(best['mu_z_au']):.10f}) au\n"
        )
        handle.write(
            f"- Best cube magnitude: "
            f"{float(best['magnitude_au']):.10f} au\n"
        )
        handle.write(
            f"- Absolute cosine with ORCA: "
            f"{float(best['absolute_cosine_to_ORCA']):.10f}\n"
        )
        handle.write(
            f"- Relative magnitude error: "
            f"{float(best['relative_magnitude_error']):.6%}\n\n"
        )

        handle.write("## Pilot acceptance\n\n")
        handle.write(
            f"- Net charge target |Q| <= {NET_CHARGE_ABS_TARGET:.1e}: "
            f"{'PASS' if net_charge_pass else 'REVIEW'}\n"
        )
        handle.write(
            f"- Dipole cosine target >= {DIPOLE_COSINE_TARGET:.2f}: "
            f"{'PASS' if dipole_cosine_pass else 'REVIEW'}\n"
        )
        handle.write(
            f"- Dipole magnitude error <= "
            f"{DIPOLE_RELATIVE_ERROR_TARGET:.1%}: "
            f"{'PASS' if dipole_magnitude_pass else 'REVIEW'}\n"
        )
        handle.write(
            f"- Boundary maximum ratio <= "
            f"{BOUNDARY_MAX_RATIO_TARGET:.1e}: "
            f"{'PASS' if boundary_pass else 'REVIEW'}\n"
        )
        handle.write(
            f"- Overall pilot status: "
            f"{'PASS' if overall_pass else 'REVIEW'}\n\n"
        )

        handle.write("## Interpretation boundary\n\n")
        handle.write(
            "The generic cube comment `Total electron density` is an ORCA "
            "header string and does not override the transition-density "
            "identity established by the `cistp02` filename and the "
            "`orca_plot` generation path. A production grid should only be "
            "frozen after the integral, dipole, and boundary diagnostics are "
            "reviewed together.\n"
        )

    log("Day019 PYR2 S2 transition-density cube validation completed.")
    log(f"Grid points: {density.size}/64000")
    log(
        f"Net transition charge: {net_transition_charge:.6e}"
    )
    log(
        f"Boundary max ratio: {boundary_max_ratio:.6e}"
    )
    log(
        f"Best dipole convention: {best['convention']}"
    )
    log(
        f"Cube |mu|: {float(best['magnitude_au']):.6f} au"
    )
    log(
        f"ORCA |mu|: "
        f"{float(reference['dipole_magnitude_au']):.6f} au"
    )
    log(
        f"Absolute cosine: "
        f"{float(best['absolute_cosine_to_ORCA']):.6f}"
    )
    log(
        f"Relative magnitude error: "
        f"{float(best['relative_magnitude_error']):.6%}"
    )
    log(
        f"Overall pilot status: "
        f"{'PASS' if overall_pass else 'REVIEW'}"
    )
    log(f"Wrote: {OUTPUT_ROOT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
