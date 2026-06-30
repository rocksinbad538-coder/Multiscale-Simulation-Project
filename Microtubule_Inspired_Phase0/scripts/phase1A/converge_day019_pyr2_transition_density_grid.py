#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

PILOT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_transition_density_pilot"
)

SITE_DIR = PILOT_ROOT / "inputs/PYR2"
GENERATION_ROOT = PILOT_ROOT / "grid_convergence_generation"
ANALYSIS_ROOT = PILOT_ROOT / "grid_convergence_analysis"

GBW_PATH = SITE_DIR / "pilot.gbw"
EXPECTED_CIS_PATH = SITE_DIR / "frame000_PYR2_embedding.cis"
GENERIC_CUBE = SITE_DIR / "pilot.cistp02.cube"

DIPOLE_TABLE = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_transition_dipole_geometry_audit/"
    "transition_dipole_observations.csv"
)

METRICS_CSV = ANALYSIS_ROOT / "PYR2_S2_GRID_CONVERGENCE_METRICS_DAY019.csv"
REPORT_MD = ANALYSIS_ROOT / "PYR2_S2_GRID_CONVERGENCE_DAY019.md"

BOHR_TO_ANGSTROM = 0.529177210903
SQRT2 = math.sqrt(2.0)

NET_CHARGE_ABS_TARGET = 1.0e-3
BOUNDARY_MAX_RATIO_TARGET = 1.0e-3
COSINE_TARGET = 0.99
SCALED_DIPOLE_ERROR_TARGET = 0.02
GRID_CONVERGENCE_TARGET = 0.005


def log(message: str = "") -> None:
    print(message, flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate and audit PYR2 S2 transition-density cubes at "
            "multiple ORCA grid resolutions."
        )
    )
    parser.add_argument(
        "--grids",
        type=int,
        nargs="+",
        default=[40, 80, 120],
        help="Cube grid sizes to analyze. Default: 40 80 120.",
    )
    parser.add_argument(
        "--generate",
        type=int,
        nargs="*",
        default=[80, 120],
        help=(
            "Grid sizes to generate with orca_plot. "
            "Default: 80 120. Pass --generate with no values "
            "to skip generation."
        ),
    )
    return parser.parse_args()


def read_orca_reference() -> np.ndarray:
    if not DIPOLE_TABLE.is_file():
        raise SystemExit(f"Missing transition-dipole table: {DIPOLE_TABLE}")

    with DIPOLE_TABLE.open(newline="", encoding="utf-8") as handle:
        rows = [
            row
            for row in csv.DictReader(handle)
            if int(row["frame"]) == 0
            and row["site"] == "PYR2"
            and int(row["root"]) == 2
            and row["family"] == "bright_like"
        ]

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected one PYR2 frame000 S2 bright row, found {len(rows)}."
        )

    row = rows[0]
    return np.array(
        [
            float(row["DX_au"]),
            float(row["DY_au"]),
            float(row["DZ_au"]),
        ],
        dtype=np.float64,
    )


def target_cube(grid: int) -> Path:
    return SITE_DIR / f"pilot.cistp02.N{grid}.cube"


def preserve_n40() -> None:
    destination = target_cube(40)

    if destination.is_file():
        return

    if not GENERIC_CUBE.is_file():
        raise RuntimeError(
            f"Neither {destination} nor the original {GENERIC_CUBE} exists."
        )

    shutil.copy2(GENERIC_CUBE, destination)
    log(f"Preserved original cube: {destination.relative_to(PROJECT_ROOT)}")


def run_orca_plot(grid: int) -> None:
    executable = shutil.which("orca_plot")

    if executable is None:
        raise SystemExit("orca_plot was not found in PATH.")

    if not GBW_PATH.exists():
        raise SystemExit(f"Missing GBW file: {GBW_PATH}")

    if not EXPECTED_CIS_PATH.exists():
        raise SystemExit(
            f"Missing CIS file with GBW-embedded basename: {EXPECTED_CIS_PATH}"
        )

    GENERATION_ROOT.mkdir(parents=True, exist_ok=True)
    transcript = GENERATION_ROOT / f"orca_plot_PYR2_S2_N{grid}.log"
    destination = target_cube(grid)

    if destination.is_file():
        log(f"[N={grid}] Existing cube retained: {destination.name}")
        return

    if GENERIC_CUBE.exists():
        GENERIC_CUBE.unlink()

    sequence = f"4\n{grid}\n7\ny\n2\n12\n"

    log(f"[N={grid}] Launching orca_plot...")

    process = subprocess.Popen(
        [executable, "pilot.gbw", "-i"],
        cwd=SITE_DIR,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    if process.stdin is None or process.stdout is None:
        raise RuntimeError("Failed to open orca_plot pipes.")

    process.stdin.write(sequence)
    process.stdin.flush()
    process.stdin.close()

    with transcript.open("w", encoding="utf-8") as handle:
        for line in process.stdout:
            print(line, end="", flush=True)
            handle.write(line)

    return_code = process.wait()

    if return_code != 0:
        raise RuntimeError(
            f"orca_plot failed for N={grid} with exit code {return_code}."
        )

    if not GENERIC_CUBE.is_file():
        raise RuntimeError(
            f"orca_plot completed but did not create {GENERIC_CUBE}."
        )

    GENERIC_CUBE.rename(destination)

    if "FINISHED CIS TRANSITION DENSITIES" not in transcript.read_text(
        encoding="utf-8",
        errors="ignore",
    ):
        raise RuntimeError(
            f"Completion marker missing from {transcript}."
        )

    log(
        f"[N={grid}] Generated {destination.name} "
        f"({destination.stat().st_size} bytes)"
    )


def parse_cube(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        comment_1 = handle.readline().rstrip("\n")
        comment_2 = handle.readline().rstrip("\n")

        tokens = handle.readline().split()
        natoms_raw = int(tokens[0])
        natoms = abs(natoms_raw)
        origin = np.array(
            [float(tokens[1]), float(tokens[2]), float(tokens[3])],
            dtype=np.float64,
        )

        counts: list[int] = []
        axes: list[np.ndarray] = []

        for _ in range(3):
            tokens = handle.readline().split()
            counts.append(abs(int(tokens[0])))
            axes.append(
                np.array(
                    [float(tokens[1]), float(tokens[2]), float(tokens[3])],
                    dtype=np.float64,
                )
            )

        atomic_numbers: list[int] = []
        atom_coordinates: list[np.ndarray] = []

        for _ in range(natoms):
            tokens = handle.readline().split()
            atomic_numbers.append(int(float(tokens[0])))
            atom_coordinates.append(
                np.array(
                    [float(tokens[2]), float(tokens[3]), float(tokens[4])],
                    dtype=np.float64,
                )
            )

        if natoms_raw < 0:
            handle.readline()

        values: list[float] = []
        for line in handle:
            values.extend(float(token) for token in line.split())

    nx, ny, nz = counts
    expected = nx * ny * nz

    if len(values) != expected:
        raise RuntimeError(
            f"{path}: expected {expected} values, found {len(values)}."
        )

    density = np.asarray(values, dtype=np.float64).reshape(
        (nx, ny, nz),
        order="C",
    )
    axes_array = np.vstack(axes)
    voxel_volume = abs(float(np.linalg.det(axes_array)))

    return {
        "comment_1": comment_1,
        "comment_2": comment_2,
        "natoms": natoms,
        "origin": origin,
        "counts": np.asarray(counts, dtype=np.int64),
        "axes": axes_array,
        "atomic_numbers": np.asarray(atomic_numbers, dtype=np.int64),
        "atom_coordinates": np.vstack(atom_coordinates),
        "density": density,
        "voxel_volume": voxel_volume,
    }


def analyze_cube(
    grid: int,
    path: Path,
    orca_mu: np.ndarray,
) -> dict[str, object]:
    cube = parse_cube(path)

    density = np.asarray(cube["density"], dtype=np.float64)
    counts = np.asarray(cube["counts"], dtype=np.int64)
    axes = np.asarray(cube["axes"], dtype=np.float64)
    origin = np.asarray(cube["origin"], dtype=np.float64)
    dvol = float(cube["voxel_volume"])

    nx, ny, nz = [int(value) for value in counts]

    if (nx, ny, nz) != (grid, grid, grid):
        raise RuntimeError(
            f"{path}: expected {grid}^3 grid, found {nx}x{ny}x{nz}."
        )

    i, j, k = np.indices((nx, ny, nz), dtype=np.float64)

    x = origin[0] + i * axes[0, 0] + j * axes[1, 0] + k * axes[2, 0]
    y = origin[1] + i * axes[0, 1] + j * axes[1, 1] + k * axes[2, 1]
    z = origin[2] + i * axes[0, 2] + j * axes[1, 2] + k * axes[2, 2]

    atomic_numbers = np.asarray(cube["atomic_numbers"], dtype=np.int64)
    atom_coordinates = np.asarray(
        cube["atom_coordinates"],
        dtype=np.float64,
    )
    carbon_mask = atomic_numbers == 6
    carbon_centroid = np.mean(atom_coordinates[carbon_mask], axis=0)

    net_charge = float(np.sum(density) * dvol)

    raw_mu = np.array(
        [
            float(np.sum(density * (x - carbon_centroid[0])) * dvol),
            float(np.sum(density * (y - carbon_centroid[1])) * dvol),
            float(np.sum(density * (z - carbon_centroid[2])) * dvol),
        ],
        dtype=np.float64,
    )

    dot = float(np.dot(raw_mu, orca_mu))
    gauge_sign = 1.0 if dot >= 0.0 else -1.0
    aligned_mu = gauge_sign * raw_mu

    raw_magnitude = float(np.linalg.norm(aligned_mu))
    orca_magnitude = float(np.linalg.norm(orca_mu))
    cosine = float(
        np.dot(aligned_mu, orca_mu)
        / (raw_magnitude * orca_magnitude)
    )

    sqrt2_mu = SQRT2 * aligned_mu
    sqrt2_magnitude = float(np.linalg.norm(sqrt2_mu))

    raw_ratio = raw_magnitude / orca_magnitude
    sqrt2_ratio = sqrt2_magnitude / orca_magnitude
    empirical_scale = orca_magnitude / raw_magnitude

    boundary = np.zeros_like(density, dtype=bool)
    boundary[0, :, :] = True
    boundary[-1, :, :] = True
    boundary[:, 0, :] = True
    boundary[:, -1, :] = True
    boundary[:, :, 0] = True
    boundary[:, :, -1] = True

    max_abs_density = float(np.max(np.abs(density)))
    max_abs_boundary = float(np.max(np.abs(density[boundary])))
    boundary_ratio = max_abs_boundary / max_abs_density

    grid_lengths_bohr = np.array(
        [
            float(np.linalg.norm(axes[0])) * (nx - 1),
            float(np.linalg.norm(axes[1])) * (ny - 1),
            float(np.linalg.norm(axes[2])) * (nz - 1),
        ],
        dtype=np.float64,
    )

    return {
        "grid": grid,
        "cube_file": str(path.relative_to(PROJECT_ROOT)),
        "file_size_bytes": path.stat().st_size,
        "natoms": int(cube["natoms"]),
        "n_values": int(density.size),
        "voxel_volume_bohr3": dvol,
        "grid_length_x_A": grid_lengths_bohr[0] * BOHR_TO_ANGSTROM,
        "grid_length_y_A": grid_lengths_bohr[1] * BOHR_TO_ANGSTROM,
        "grid_length_z_A": grid_lengths_bohr[2] * BOHR_TO_ANGSTROM,
        "net_transition_charge": net_charge,
        "boundary_max_ratio": boundary_ratio,
        "gauge_sign": int(gauge_sign),
        "aligned_mu_x_au": aligned_mu[0],
        "aligned_mu_y_au": aligned_mu[1],
        "aligned_mu_z_au": aligned_mu[2],
        "raw_mu_magnitude_au": raw_magnitude,
        "ORCA_mu_magnitude_au": orca_magnitude,
        "cosine_to_ORCA": cosine,
        "raw_to_ORCA_ratio": raw_ratio,
        "raw_relative_error": abs(raw_ratio - 1.0),
        "sqrt2_scaled_mu_magnitude_au": sqrt2_magnitude,
        "sqrt2_scaled_to_ORCA_ratio": sqrt2_ratio,
        "sqrt2_scaled_relative_error": abs(sqrt2_ratio - 1.0),
        "empirical_scale_to_ORCA": empirical_scale,
        "empirical_scale_over_sqrt2": empirical_scale / SQRT2,
        "net_charge_pass": abs(net_charge) <= NET_CHARGE_ABS_TARGET,
        "boundary_pass": boundary_ratio <= BOUNDARY_MAX_RATIO_TARGET,
        "cosine_pass": cosine >= COSINE_TARGET,
        "sqrt2_scaled_dipole_pass": (
            abs(sqrt2_ratio - 1.0) <= SCALED_DIPOLE_ERROR_TARGET
        ),
    }


def write_csv(rows: list[dict[str, object]]) -> None:
    ANALYSIS_ROOT.mkdir(parents=True, exist_ok=True)

    with METRICS_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(rows)


def write_report(rows: list[dict[str, object]]) -> None:
    finest = rows[-1]
    previous = rows[-2] if len(rows) >= 2 else None

    if previous is None:
        convergence = float("nan")
        convergence_pass = False
    else:
        convergence = abs(
            float(previous["raw_mu_magnitude_au"])
            - float(finest["raw_mu_magnitude_au"])
        ) / float(finest["raw_mu_magnitude_au"])
        convergence_pass = convergence <= GRID_CONVERGENCE_TARGET

    overall = (
        bool(finest["net_charge_pass"])
        and bool(finest["boundary_pass"])
        and bool(finest["cosine_pass"])
        and bool(finest["sqrt2_scaled_dipole_pass"])
        and convergence_pass
    )

    with REPORT_MD.open("w", encoding="utf-8") as handle:
        handle.write(
            "# Day019 PYR2 S2 transition-density grid convergence\n\n"
        )

        handle.write("## Resolution study\n\n")
        handle.write(
            "| Grid | Values | Raw |mu| (au) | Raw/ORCA | "
            "sqrt(2)-scaled |mu| (au) | Scaled/ORCA | "
            "Empirical scale | Qtr | Boundary ratio |\n"
        )
        handle.write(
            "|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n"
        )

        for row in rows:
            handle.write(
                f"| {row['grid']}^3 "
                f"| {row['n_values']} "
                f"| {float(row['raw_mu_magnitude_au']):.9f} "
                f"| {float(row['raw_to_ORCA_ratio']):.9f} "
                f"| {float(row['sqrt2_scaled_mu_magnitude_au']):.9f} "
                f"| {float(row['sqrt2_scaled_to_ORCA_ratio']):.9f} "
                f"| {float(row['empirical_scale_to_ORCA']):.9f} "
                f"| {float(row['net_transition_charge']):.3e} "
                f"| {float(row['boundary_max_ratio']):.3e} |\n"
            )

        handle.write("\n## Finest-grid validation\n\n")
        handle.write(
            f"- Finest grid: {finest['grid']}^3\n"
        )
        handle.write(
            f"- Directional cosine with ORCA: "
            f"{float(finest['cosine_to_ORCA']):.10f}\n"
        )
        handle.write(
            f"- Raw cube/ORCA dipole ratio: "
            f"{float(finest['raw_to_ORCA_ratio']):.10f}\n"
        )
        handle.write(
            f"- sqrt(2)-scaled cube/ORCA ratio: "
            f"{float(finest['sqrt2_scaled_to_ORCA_ratio']):.10f}\n"
        )
        handle.write(
            f"- sqrt(2)-scaled magnitude error: "
            f"{float(finest['sqrt2_scaled_relative_error']):.6%}\n"
        )
        handle.write(
            f"- Empirical scale factor: "
            f"{float(finest['empirical_scale_to_ORCA']):.10f}\n"
        )
        handle.write(
            f"- Empirical scale / sqrt(2): "
            f"{float(finest['empirical_scale_over_sqrt2']):.10f}\n"
        )
        handle.write(
            f"- Previous-to-finest raw-dipole change: "
            f"{convergence:.6%}\n\n"
        )

        handle.write("## Decision controls\n\n")
        handle.write(
            f"- Net transition charge: "
            f"{'PASS' if finest['net_charge_pass'] else 'REVIEW'}\n"
        )
        handle.write(
            f"- Boundary containment: "
            f"{'PASS' if finest['boundary_pass'] else 'REVIEW'}\n"
        )
        handle.write(
            f"- Dipole direction: "
            f"{'PASS' if finest['cosine_pass'] else 'REVIEW'}\n"
        )
        handle.write(
            f"- sqrt(2)-scaled magnitude within "
            f"{SCALED_DIPOLE_ERROR_TARGET:.1%}: "
            f"{'PASS' if finest['sqrt2_scaled_dipole_pass'] else 'REVIEW'}\n"
        )
        handle.write(
            f"- Grid convergence within "
            f"{GRID_CONVERGENCE_TARGET:.1%}: "
            f"{'PASS' if convergence_pass else 'REVIEW'}\n"
        )
        handle.write(
            f"- Overall convergence status: "
            f"{'PASS' if overall else 'REVIEW'}\n\n"
        )

        handle.write("## Interpretation boundary\n\n")
        handle.write(
            "The near-1/sqrt(2) raw dipole ratio is treated here as a "
            "normalization hypothesis, not as an assumed ORCA convention. "
            "It is accepted only if the ratio is resolution-independent and "
            "the sqrt(2)-scaled moment converges to the independently printed "
            "ORCA transition dipole. Before production transition-density "
            "couplings, the same normalization test must also pass for PYR3, "
            "PYR4, and PYR5. State-specific empirical dipole normalization "
            "remains the fallback because it guarantees the correct far-field "
            "limit of each transition density.\n"
        )

    log("")
    log("PYR2 S2 transition-density grid convergence completed.")
    for row in rows:
        log(
            f"N={row['grid']:3d}: raw |mu|="
            f"{float(row['raw_mu_magnitude_au']):.6f} au, "
            f"sqrt(2)-scaled |mu|="
            f"{float(row['sqrt2_scaled_mu_magnitude_au']):.6f} au, "
            f"scaled error="
            f"{float(row['sqrt2_scaled_relative_error']):.4%}"
        )
    log(
        f"Previous-to-finest change: {convergence:.4%}"
    )
    log(
        f"Overall convergence status: "
        f"{'PASS' if overall else 'REVIEW'}"
    )
    log(f"Wrote: {ANALYSIS_ROOT.relative_to(PROJECT_ROOT)}")


def main() -> None:
    args = parse_args()

    grids = sorted(set(args.grids))
    generate = sorted(set(args.generate))

    if any(grid <= 1 for grid in grids + generate):
        raise SystemExit("All grid sizes must be greater than 1.")

    preserve_n40()

    for grid in generate:
        run_orca_plot(grid)

    missing = [
        target_cube(grid)
        for grid in grids
        if not target_cube(grid).is_file()
    ]

    if missing:
        raise RuntimeError(
            "Missing requested cubes: "
            + ", ".join(str(path) for path in missing)
        )

    orca_mu = read_orca_reference()
    rows = [
        analyze_cube(grid, target_cube(grid), orca_mu)
        for grid in grids
    ]

    write_csv(rows)
    write_report(rows)


if __name__ == "__main__":
    main()
