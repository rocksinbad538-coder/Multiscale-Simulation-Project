#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import os
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

INPUT_ROOT = PILOT_ROOT / "inputs"

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_bright_transition_density_N80"
)

RAW_CUBE_ROOT = OUTPUT_ROOT / "raw_cubes"
NORMALIZED_CUBE_ROOT = OUTPUT_ROOT / "sqrt2_normalized_cubes"
LOG_ROOT = OUTPUT_ROOT / "generation_logs"

METRICS_CSV = OUTPUT_ROOT / "BRIGHT_TRANSITION_DENSITY_N80_METRICS_DAY019.csv"
MANIFEST_CSV = OUTPUT_ROOT / "BRIGHT_TRANSITION_DENSITY_N80_MANIFEST_DAY019.csv"
REPORT_MD = OUTPUT_ROOT / "BRIGHT_TRANSITION_DENSITY_N80_VALIDATION_DAY019.md"

DIPOLE_TABLE = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_transition_dipole_geometry_audit/"
    "transition_dipole_observations.csv"
)

SITES = ("PYR2", "PYR3", "PYR4", "PYR5")
BRIGHT_ROOTS = {
    "PYR2": 2,
    "PYR3": 2,
    "PYR4": 2,
    "PYR5": 1,
}

GRID = 80
SQRT2 = math.sqrt(2.0)
BOHR_TO_ANGSTROM = 0.529177210903

NET_CHARGE_ABS_TARGET = 1.0e-3
BOUNDARY_MAX_RATIO_TARGET = 1.0e-3
DIPOLE_COSINE_TARGET = 0.99
SCALED_DIPOLE_RELATIVE_ERROR_TARGET = 0.02


def log(message: str = "") -> None:
    print(message, flush=True)


def safe_symlink(target: Path, link: Path) -> None:
    if link.exists() or link.is_symlink():
        if link.is_symlink() and link.resolve() == target.resolve():
            return
        link.unlink()

    relative_target = os.path.relpath(target, start=link.parent)
    link.symlink_to(relative_target)


def read_orca_references() -> dict[str, dict[str, object]]:
    if not DIPOLE_TABLE.is_file():
        raise SystemExit(f"Missing dipole table: {DIPOLE_TABLE}")

    references: dict[str, dict[str, object]] = {}

    with DIPOLE_TABLE.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            site = row["site"]
            root = int(row["root"])

            if (
                int(row["frame"]) != 0
                or site not in SITES
                or root != BRIGHT_ROOTS[site]
                or row["family"] != "bright_like"
            ):
                continue

            vector = np.array(
                [
                    float(row["DX_au"]),
                    float(row["DY_au"]),
                    float(row["DZ_au"]),
                ],
                dtype=np.float64,
            )

            references[site] = {
                "root": root,
                "vector_au": vector,
                "magnitude_au": float(np.linalg.norm(vector)),
                "fosc": float(row["fosc"]),
                "energy_eV": float(row["state_energy_eV"]),
            }

    if set(references) != set(SITES):
        raise RuntimeError(
            f"Expected references for {SITES}, found {sorted(references)}."
        )

    return references


def ensure_input_aliases(site: str) -> tuple[Path, Path]:
    site_dir = INPUT_ROOT / site
    gbw = site_dir / "pilot.gbw"
    pilot_cis = site_dir / "pilot.cis"
    expected_cis = site_dir / f"frame000_{site}_embedding.cis"

    if not gbw.exists():
        raise SystemExit(f"Missing GBW file: {gbw}")

    if not pilot_cis.exists():
        raise SystemExit(f"Missing CIS file: {pilot_cis}")

    safe_symlink(pilot_cis.resolve(), expected_cis)

    return gbw, expected_cis


def raw_source_cube(site: str, root: int) -> Path:
    return INPUT_ROOT / site / f"pilot.cistp{root:02d}.N{GRID}.cube"


def generic_cube(site: str, root: int) -> Path:
    return INPUT_ROOT / site / f"pilot.cistp{root:02d}.cube"


def generate_cube(site: str, root: int) -> Path:
    gbw, expected_cis = ensure_input_aliases(site)
    source_cube = raw_source_cube(site, root)

    if source_cube.is_file():
        log(f"[{site}] Reusing existing N{GRID} cube: {source_cube.name}")
        return source_cube

    executable = shutil.which("orca_plot")
    if executable is None:
        raise SystemExit("orca_plot was not found in PATH.")

    site_dir = INPUT_ROOT / site
    generic = generic_cube(site, root)

    if generic.exists():
        generic.unlink()

    LOG_ROOT.mkdir(parents=True, exist_ok=True)
    transcript = LOG_ROOT / f"orca_plot_{site}_bright_S{root}_N{GRID}.log"

    # 4 -> grid intervals, GRID -> N x N x N
    # 7 -> transition density, y -> GBW-embedded CIS filename
    # root -> requested TDA state, 12 -> exit
    sequence = f"4\n{GRID}\n7\ny\n{root}\n12\n"

    log(
        f"[{site}] Generating bright S{root} transition density "
        f"on {GRID}^3 grid..."
    )

    process = subprocess.Popen(
        [executable, gbw.name, "-i"],
        cwd=site_dir,
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
            f"orca_plot failed for {site} S{root} with code {return_code}."
        )

    if not generic.is_file():
        raise RuntimeError(
            f"orca_plot completed but did not create {generic}."
        )

    transcript_text = transcript.read_text(
        encoding="utf-8",
        errors="ignore",
    )

    if "FINISHED CIS TRANSITION DENSITIES" not in transcript_text:
        raise RuntimeError(
            f"Completion marker missing from {transcript}."
        )

    generic.rename(source_cube)

    if not expected_cis.exists():
        raise RuntimeError(
            f"CIS alias disappeared unexpectedly: {expected_cis}"
        )

    log(
        f"[{site}] Generated {source_cube.name} "
        f"({source_cube.stat().st_size} bytes)"
    )

    return source_cube


def parse_cube(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        comment_1 = handle.readline().rstrip("\n")
        comment_2 = handle.readline().rstrip("\n")

        tokens = handle.readline().split()
        if len(tokens) < 4:
            raise RuntimeError(f"Invalid cube origin line: {path}")

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
            if len(tokens) < 4:
                raise RuntimeError(f"Invalid cube axis line: {path}")

            counts.append(abs(int(tokens[0])))
            axes.append(
                np.array(
                    [float(tokens[1]), float(tokens[2]), float(tokens[3])],
                    dtype=np.float64,
                )
            )

        atom_records: list[tuple[int, float, float, float, float]] = []

        for _ in range(natoms):
            tokens = handle.readline().split()
            if len(tokens) < 5:
                raise RuntimeError(f"Invalid cube atom line: {path}")

            atom_records.append(
                (
                    int(float(tokens[0])),
                    float(tokens[1]),
                    float(tokens[2]),
                    float(tokens[3]),
                    float(tokens[4]),
                )
            )

        dataset_line: str | None = None
        if natoms_raw < 0:
            dataset_line = handle.readline().rstrip("\n")

        values: list[float] = []
        for line in handle:
            values.extend(float(token) for token in line.split())

    nx, ny, nz = counts
    expected_values = nx * ny * nz

    if len(values) != expected_values:
        raise RuntimeError(
            f"{path}: expected {expected_values} values, "
            f"found {len(values)}."
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
        "natoms_raw": natoms_raw,
        "natoms": natoms,
        "origin": origin,
        "counts": np.asarray(counts, dtype=np.int64),
        "axes": axes_array,
        "atom_records": atom_records,
        "dataset_line": dataset_line,
        "density": density,
        "voxel_volume": voxel_volume,
    }


def write_scaled_cube(
    source: dict[str, object],
    output_path: Path,
    scale: float,
    site: str,
    root: int,
) -> None:
    counts = np.asarray(source["counts"], dtype=np.int64)
    axes = np.asarray(source["axes"], dtype=np.float64)
    origin = np.asarray(source["origin"], dtype=np.float64)
    density = np.asarray(source["density"], dtype=np.float64) * scale
    atom_records = list(source["atom_records"])

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(
            f"ORCA transition density: frame000 {site} bright S{root}\n"
        )
        handle.write(
            f"Gauge aligned to ORCA transition dipole; density scale={scale:.15g}\n"
        )
        handle.write(
            f"{int(source['natoms_raw']):5d}"
            f"{origin[0]:13.6f}{origin[1]:13.6f}{origin[2]:13.6f}\n"
        )

        for count, axis in zip(counts, axes):
            handle.write(
                f"{int(count):5d}"
                f"{axis[0]:13.6f}{axis[1]:13.6f}{axis[2]:13.6f}\n"
            )

        for atomic_number, charge, x, y, z in atom_records:
            handle.write(
                f"{atomic_number:5d}"
                f"{charge:13.6f}"
                f"{x:13.6f}{y:13.6f}{z:13.6f}\n"
            )

        if source["dataset_line"] is not None:
            handle.write(str(source["dataset_line"]) + "\n")

        flat = density.ravel(order="C")
        for start in range(0, flat.size, 6):
            values = flat[start : start + 6]
            handle.write(
                "".join(f"{value:13.5e}" for value in values) + "\n"
            )


def analyze_and_normalize(
    site: str,
    root: int,
    path: Path,
    reference: dict[str, object],
) -> tuple[dict[str, object], Path]:
    cube = parse_cube(path)

    density = np.asarray(cube["density"], dtype=np.float64)
    counts = np.asarray(cube["counts"], dtype=np.int64)
    axes = np.asarray(cube["axes"], dtype=np.float64)
    origin = np.asarray(cube["origin"], dtype=np.float64)
    dvol = float(cube["voxel_volume"])

    nx, ny, nz = [int(value) for value in counts]

    if (nx, ny, nz) != (GRID, GRID, GRID):
        raise RuntimeError(
            f"{path}: expected {GRID}^3 grid, found "
            f"{nx}x{ny}x{nz}."
        )

    i, j, k = np.indices((nx, ny, nz), dtype=np.float64)

    x = origin[0] + i * axes[0, 0] + j * axes[1, 0] + k * axes[2, 0]
    y = origin[1] + i * axes[0, 1] + j * axes[1, 1] + k * axes[2, 1]
    z = origin[2] + i * axes[0, 2] + j * axes[1, 2] + k * axes[2, 2]

    atom_records = list(cube["atom_records"])
    atomic_numbers = np.array(
        [record[0] for record in atom_records],
        dtype=np.int64,
    )
    atom_coordinates = np.array(
        [[record[2], record[3], record[4]] for record in atom_records],
        dtype=np.float64,
    )

    carbon_mask = atomic_numbers == 6

    if int(np.count_nonzero(carbon_mask)) != 16:
        raise RuntimeError(
            f"{site}: expected 16 carbon atoms, found "
            f"{int(np.count_nonzero(carbon_mask))}."
        )

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

    reference_vector = np.asarray(
        reference["vector_au"],
        dtype=np.float64,
    )

    gauge_sign = 1.0 if float(np.dot(raw_mu, reference_vector)) >= 0.0 else -1.0
    aligned_mu = gauge_sign * raw_mu
    normalized_mu = SQRT2 * aligned_mu

    raw_magnitude = float(np.linalg.norm(aligned_mu))
    normalized_magnitude = float(np.linalg.norm(normalized_mu))
    reference_magnitude = float(np.linalg.norm(reference_vector))

    cosine = float(
        np.dot(aligned_mu, reference_vector)
        / (raw_magnitude * reference_magnitude)
    )

    raw_ratio = raw_magnitude / reference_magnitude
    normalized_ratio = normalized_magnitude / reference_magnitude
    normalized_relative_error = abs(normalized_ratio - 1.0)
    empirical_scale = reference_magnitude / raw_magnitude

    boundary = np.zeros_like(density, dtype=bool)
    boundary[0, :, :] = True
    boundary[-1, :, :] = True
    boundary[:, 0, :] = True
    boundary[:, -1, :] = True
    boundary[:, :, 0] = True
    boundary[:, :, -1] = True

    maximum_absolute_density = float(np.max(np.abs(density)))
    maximum_absolute_boundary_density = float(
        np.max(np.abs(density[boundary]))
    )
    boundary_ratio = (
        maximum_absolute_boundary_density / maximum_absolute_density
    )

    absolute_density_integral = float(
        np.sum(np.abs(density)) * dvol
    )

    normalized_net_charge = gauge_sign * SQRT2 * net_charge

    grid_lengths_bohr = np.array(
        [
            float(np.linalg.norm(axes[0])) * (nx - 1),
            float(np.linalg.norm(axes[1])) * (ny - 1),
            float(np.linalg.norm(axes[2])) * (nz - 1),
        ],
        dtype=np.float64,
    )

    net_charge_pass = abs(net_charge) <= NET_CHARGE_ABS_TARGET
    boundary_pass = boundary_ratio <= BOUNDARY_MAX_RATIO_TARGET
    cosine_pass = cosine >= DIPOLE_COSINE_TARGET
    normalized_dipole_pass = (
        normalized_relative_error
        <= SCALED_DIPOLE_RELATIVE_ERROR_TARGET
    )

    overall_pass = (
        net_charge_pass
        and boundary_pass
        and cosine_pass
        and normalized_dipole_pass
    )

    normalized_path = (
        NORMALIZED_CUBE_ROOT
        / f"frame000_{site}_bright_S{root}_tdens_N{GRID}_sqrt2.cube"
    )

    write_scaled_cube(
        source=cube,
        output_path=normalized_path,
        scale=gauge_sign * SQRT2,
        site=site,
        root=root,
    )

    metrics = {
        "frame": 0,
        "site": site,
        "bright_root": root,
        "grid": GRID,
        "n_values": int(density.size),
        "raw_cube": str(path.relative_to(PROJECT_ROOT)),
        "normalized_cube": str(normalized_path.relative_to(PROJECT_ROOT)),
        "raw_cube_size_bytes": path.stat().st_size,
        "normalized_cube_size_bytes": normalized_path.stat().st_size,
        "voxel_volume_bohr3": dvol,
        "grid_length_x_A": grid_lengths_bohr[0] * BOHR_TO_ANGSTROM,
        "grid_length_y_A": grid_lengths_bohr[1] * BOHR_TO_ANGSTROM,
        "grid_length_z_A": grid_lengths_bohr[2] * BOHR_TO_ANGSTROM,
        "raw_net_transition_charge": net_charge,
        "normalized_net_transition_charge": normalized_net_charge,
        "absolute_density_integral": absolute_density_integral,
        "boundary_max_ratio": boundary_ratio,
        "gauge_sign": int(gauge_sign),
        "raw_mu_x_au_aligned": aligned_mu[0],
        "raw_mu_y_au_aligned": aligned_mu[1],
        "raw_mu_z_au_aligned": aligned_mu[2],
        "raw_mu_magnitude_au": raw_magnitude,
        "ORCA_mu_x_au": reference_vector[0],
        "ORCA_mu_y_au": reference_vector[1],
        "ORCA_mu_z_au": reference_vector[2],
        "ORCA_mu_magnitude_au": reference_magnitude,
        "cosine_to_ORCA": cosine,
        "raw_to_ORCA_ratio": raw_ratio,
        "sqrt2_normalized_mu_magnitude_au": normalized_magnitude,
        "sqrt2_normalized_to_ORCA_ratio": normalized_ratio,
        "sqrt2_normalized_relative_error": normalized_relative_error,
        "empirical_scale_to_ORCA": empirical_scale,
        "empirical_scale_over_sqrt2": empirical_scale / SQRT2,
        "net_charge_pass": net_charge_pass,
        "boundary_pass": boundary_pass,
        "cosine_pass": cosine_pass,
        "sqrt2_normalization_pass": normalized_dipole_pass,
        "overall_site_pass": overall_pass,
    }

    return metrics, normalized_path


def copy_raw_cube(site: str, root: int, source: Path) -> Path:
    RAW_CUBE_ROOT.mkdir(parents=True, exist_ok=True)
    destination = (
        RAW_CUBE_ROOT
        / f"frame000_{site}_bright_S{root}_tdens_N{GRID}_raw.cube"
    )

    if not destination.is_file():
        shutil.copy2(source, destination)

    return destination


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(rows)


def write_report(
    metrics_rows: list[dict[str, object]],
    manifest_rows: list[dict[str, object]],
) -> None:
    all_pass = all(
        bool(row["overall_site_pass"])
        for row in metrics_rows
    )

    maximum_error = max(
        float(row["sqrt2_normalized_relative_error"])
        for row in metrics_rows
    )

    minimum_cosine = min(
        float(row["cosine_to_ORCA"])
        for row in metrics_rows
    )

    maximum_boundary_ratio = max(
        float(row["boundary_max_ratio"])
        for row in metrics_rows
    )

    maximum_abs_raw_charge = max(
        abs(float(row["raw_net_transition_charge"]))
        for row in metrics_rows
    )

    empirical_scale_values = [
        float(row["empirical_scale_to_ORCA"])
        for row in metrics_rows
    ]

    with REPORT_MD.open("w", encoding="utf-8") as handle:
        handle.write(
            "# Day019 bright transition-density N80 validation\n\n"
        )

        handle.write("## Production convention\n\n")
        handle.write("- Frame: 000\n")
        handle.write("- Grid: 80 x 80 x 80\n")
        handle.write(
            "- Bright roots: PYR2=S2, PYR3=S2, PYR4=S2, PYR5=S1\n"
        )
        handle.write(
            "- Each raw transition density is gauge-aligned to the "
            "independently printed ORCA transition dipole.\n"
        )
        handle.write(
            "- The aligned density is multiplied by sqrt(2), following "
            "the normalization established by the PYR2 grid-convergence "
            "study.\n\n"
        )

        handle.write("## Site validation\n\n")
        handle.write(
            "| Site | Root | Gauge sign | Raw |mu| (au) | "
            "ORCA |mu| (au) | sqrt(2)-scaled |mu| (au) | "
            "Scaled error | Cosine | Qtr | Boundary ratio | Status |\n"
        )
        handle.write(
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n"
        )

        for row in metrics_rows:
            handle.write(
                f"| {row['site']} "
                f"| S{row['bright_root']} "
                f"| {row['gauge_sign']} "
                f"| {float(row['raw_mu_magnitude_au']):.9f} "
                f"| {float(row['ORCA_mu_magnitude_au']):.9f} "
                f"| {float(row['sqrt2_normalized_mu_magnitude_au']):.9f} "
                f"| {float(row['sqrt2_normalized_relative_error']):.4%} "
                f"| {float(row['cosine_to_ORCA']):.9f} "
                f"| {float(row['raw_net_transition_charge']):.3e} "
                f"| {float(row['boundary_max_ratio']):.3e} "
                f"| {'PASS' if row['overall_site_pass'] else 'REVIEW'} |\n"
            )

        handle.write("\n## Global controls\n\n")
        handle.write(
            f"- Valid sites: "
            f"{sum(bool(row['overall_site_pass']) for row in metrics_rows)}/4\n"
        )
        handle.write(
            f"- Maximum sqrt(2)-normalized dipole error: "
            f"{maximum_error:.6%}\n"
        )
        handle.write(
            f"- Minimum directional cosine: {minimum_cosine:.10f}\n"
        )
        handle.write(
            f"- Maximum |raw transition charge|: "
            f"{maximum_abs_raw_charge:.6e}\n"
        )
        handle.write(
            f"- Maximum boundary ratio: "
            f"{maximum_boundary_ratio:.6e}\n"
        )
        handle.write(
            f"- Empirical scale range: "
            f"{min(empirical_scale_values):.9f} to "
            f"{max(empirical_scale_values):.9f}\n"
        )
        handle.write(
            f"- Overall production status: "
            f"{'PASS' if all_pass else 'REVIEW'}\n\n"
        )

        handle.write("## Interpretation boundary\n\n")
        handle.write(
            "Passing this validation establishes a common, phase-consistent "
            "and far-field-normalized bright-state transition-density set. "
            "It does not by itself establish dielectric screening. The first "
            "transition-density coupling benchmark should therefore be "
            "computed as an unscreened Coulomb integral and compared directly "
            "with the unscreened point-transition-dipole baseline.\n"
        )


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    RAW_CUBE_ROOT.mkdir(parents=True, exist_ok=True)
    NORMALIZED_CUBE_ROOT.mkdir(parents=True, exist_ok=True)
    LOG_ROOT.mkdir(parents=True, exist_ok=True)

    references = read_orca_references()

    metrics_rows: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []

    log("Day019 bright transition-density N80 production")

    for site in SITES:
        root = BRIGHT_ROOTS[site]

        source = generate_cube(site, root)
        raw_copy = copy_raw_cube(site, root, source)

        metrics, normalized = analyze_and_normalize(
            site=site,
            root=root,
            path=source,
            reference=references[site],
        )

        metrics["raw_cube_copy"] = str(
            raw_copy.relative_to(PROJECT_ROOT)
        )
        metrics_rows.append(metrics)

        manifest_rows.append(
            {
                "frame": 0,
                "site": site,
                "bright_root": root,
                "grid": GRID,
                "raw_source_cube": str(source.relative_to(PROJECT_ROOT)),
                "raw_production_cube": str(
                    raw_copy.relative_to(PROJECT_ROOT)
                ),
                "sqrt2_normalized_cube": str(
                    normalized.relative_to(PROJECT_ROOT)
                ),
                "gauge_sign": metrics["gauge_sign"],
                "density_scale_applied": (
                    int(metrics["gauge_sign"]) * SQRT2
                ),
                "status": (
                    "PASS"
                    if metrics["overall_site_pass"]
                    else "REVIEW"
                ),
            }
        )

        log(
            f"[{site}] S{root}: "
            f"raw |mu|={float(metrics['raw_mu_magnitude_au']):.6f} au, "
            f"sqrt(2)-scaled |mu|="
            f"{float(metrics['sqrt2_normalized_mu_magnitude_au']):.6f} au, "
            f"ORCA |mu|="
            f"{float(metrics['ORCA_mu_magnitude_au']):.6f} au, "
            f"error="
            f"{float(metrics['sqrt2_normalized_relative_error']):.4%}, "
            f"status="
            f"{'PASS' if metrics['overall_site_pass'] else 'REVIEW'}"
        )

    write_csv(METRICS_CSV, metrics_rows)
    write_csv(MANIFEST_CSV, manifest_rows)
    write_report(metrics_rows, manifest_rows)

    valid_count = sum(
        bool(row["overall_site_pass"])
        for row in metrics_rows
    )

    log("")
    log("Day019 bright transition-density N80 production completed.")
    log("Raw cubes: 4/4")
    log("sqrt(2)-normalized cubes: 4/4")
    log(f"Validated sites: {valid_count}/4")
    log(
        f"Maximum normalized dipole error: "
        f"{max(float(row['sqrt2_normalized_relative_error']) for row in metrics_rows):.4%}"
    )
    log(
        f"Minimum directional cosine: "
        f"{min(float(row['cosine_to_ORCA']) for row in metrics_rows):.8f}"
    )
    log(
        f"Overall status: "
        f"{'PASS' if valid_count == 4 else 'REVIEW'}"
    )
    log(f"Wrote: {OUTPUT_ROOT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
