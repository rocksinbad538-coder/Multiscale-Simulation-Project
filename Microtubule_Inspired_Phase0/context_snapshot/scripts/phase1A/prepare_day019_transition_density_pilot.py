#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
from itertools import combinations
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

SOURCE_SELECTION = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_full_two_state_manifold_analysis/"
    "SOURCE_OUTPUT_SELECTION_DAY019.csv"
)

FRAME_SUMMARY = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_full_two_state_manifold_analysis/"
    "two_state_frame_summary.csv"
)

SITE_PAIR_GEOMETRY = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_transition_dipole_geometry_audit/"
    "site_pair_geometry.csv"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_transition_density_pilot"
)

INPUT_ROOT = OUTPUT_ROOT / "inputs"
MANIFEST_CSV = OUTPUT_ROOT / "TRANSITION_DENSITY_PILOT_MANIFEST_DAY019.csv"
GEOMETRY_CSV = OUTPUT_ROOT / "POINT_DIPOLE_GEOMETRY_DIAGNOSTIC_DAY019.csv"
REPORT_MD = OUTPUT_ROOT / "TRANSITION_DENSITY_PILOT_PREPARATION_DAY019.md"
LAUNCHER_SH = OUTPUT_ROOT / "launch_orca_plot_probe.sh"

SITES = ("PYR2", "PYR3", "PYR4", "PYR5")
FRAME = 0


def log(message: str = "") -> None:
    print(message, flush=True)


def as_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows to write: {path}")

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def read_selected_frame000_outputs() -> dict[str, Path]:
    if not SOURCE_SELECTION.is_file():
        raise SystemExit(f"Missing source selection: {SOURCE_SELECTION}")

    selected: dict[str, Path] = {}

    with SOURCE_SELECTION.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if not as_bool(row["selected"]):
                continue

            frame = int(row["frame"])
            site = row["site"]

            if frame != FRAME:
                continue

            path = (PROJECT_ROOT / row["path"]).resolve()

            if site in selected:
                raise RuntimeError(
                    f"Duplicate selected frame000 output for {site}: "
                    f"{selected[site]} and {path}"
                )

            selected[site] = path

    if set(selected) != set(SITES):
        raise RuntimeError(
            f"Expected frame000 outputs for {SITES}, found {sorted(selected)}"
        )

    return selected


def read_bright_roots() -> dict[str, int]:
    if not FRAME_SUMMARY.is_file():
        raise SystemExit(f"Missing frame summary: {FRAME_SUMMARY}")

    roots: dict[str, int] = {}

    with FRAME_SUMMARY.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if int(row["frame"]) != FRAME:
                continue

            roots[row["site"]] = int(row["bright_root"])

    if set(roots) != set(SITES):
        raise RuntimeError(
            f"Expected bright roots for {SITES}, found {sorted(roots)}"
        )

    return roots


def find_companion(output_path: Path, suffix: str) -> Path:
    direct = output_path.with_suffix(suffix)

    if direct.is_file():
        return direct.resolve()

    stem = output_path.stem
    parent = output_path.parent

    exact_candidates = sorted(parent.glob(f"{stem}*{suffix}"))

    if len(exact_candidates) == 1:
        return exact_candidates[0].resolve()

    all_candidates = sorted(parent.glob(f"*{suffix}"))

    if len(all_candidates) == 1:
        return all_candidates[0].resolve()

    raise RuntimeError(
        f"Could not uniquely identify {suffix} companion for {output_path}. "
        f"Exact candidates={exact_candidates}; all candidates={all_candidates}"
    )


def parse_coordinate_line(
    line: str,
) -> tuple[str, float, float, float] | None:
    tokens = line.split()

    def is_float(token: str) -> bool:
        try:
            float(token)
        except ValueError:
            return False
        return True

    if (
        len(tokens) >= 4
        and tokens[0].isalpha()
        and is_float(tokens[1])
        and is_float(tokens[2])
        and is_float(tokens[3])
    ):
        return (
            tokens[0],
            float(tokens[1]),
            float(tokens[2]),
            float(tokens[3]),
        )

    if (
        len(tokens) >= 5
        and tokens[0].lstrip("+-").isdigit()
        and tokens[1].isalpha()
        and is_float(tokens[2])
        and is_float(tokens[3])
        and is_float(tokens[4])
    ):
        return (
            tokens[1],
            float(tokens[2]),
            float(tokens[3]),
            float(tokens[4]),
        )

    return None


def parse_qm_geometry(output_path: Path) -> tuple[list[str], np.ndarray]:
    lines = output_path.read_text(errors="ignore").splitlines()
    blocks: list[tuple[list[str], np.ndarray]] = []

    for index, line in enumerate(lines):
        if "CARTESIAN COORDINATES (ANGSTROEM)" not in line:
            continue

        symbols: list[str] = []
        coordinates: list[list[float]] = []
        started = False

        for candidate in lines[index + 1 :]:
            parsed = parse_coordinate_line(candidate)

            if parsed is None:
                if started:
                    break
                continue

            started = True
            symbol, x, y, z = parsed
            symbols.append(symbol)
            coordinates.append([x, y, z])

        if symbols:
            blocks.append(
                (
                    symbols,
                    np.array(coordinates, dtype=np.float64),
                )
            )

    if not blocks:
        raise RuntimeError(f"No coordinate block found in {output_path}")

    max_atoms = max(len(symbols) for symbols, _ in blocks)
    candidates = [
        block for block in blocks if len(block[0]) == max_atoms
    ]
    symbols, coordinates = candidates[-1]

    if len(symbols) != 26:
        raise RuntimeError(
            f"Expected 26 QM atoms in {output_path}, found {len(symbols)}"
        )

    if symbols.count("C") != 16:
        raise RuntimeError(
            f"Expected 16 carbon atoms in {output_path}, "
            f"found {symbols.count('C')}"
        )

    return symbols, coordinates


def safe_symlink(target: Path, link: Path) -> None:
    if link.exists() or link.is_symlink():
        if link.is_symlink() and link.resolve() == target.resolve():
            return
        link.unlink()

    relative_target = os.path.relpath(target, start=link.parent)
    link.symlink_to(relative_target)


def prepare_inputs(
    selected: dict[str, Path],
    bright_roots: dict[str, int],
) -> tuple[list[dict[str, object]], dict[str, dict[str, float]]]:
    manifest: list[dict[str, object]] = []
    geometry: dict[str, dict[str, float]] = {}

    for site in SITES:
        output_path = selected[site]
        gbw_path = find_companion(output_path, ".gbw")
        cis_path = find_companion(output_path, ".cis")

        site_dir = INPUT_ROOT / site
        site_dir.mkdir(parents=True, exist_ok=True)

        safe_symlink(gbw_path, site_dir / "pilot.gbw")
        safe_symlink(cis_path, site_dir / "pilot.cis")
        safe_symlink(output_path, site_dir / "source.out")

        densities_path = output_path.with_suffix(".densities")
        if densities_path.is_file():
            safe_symlink(
                densities_path.resolve(),
                site_dir / "pilot.densities",
            )

        symbols, coordinates = parse_qm_geometry(output_path)
        symbol_array = np.array(symbols)
        carbon_mask = symbol_array == "C"

        carbon_centroid = np.mean(coordinates[carbon_mask], axis=0)
        all_atom_radius = float(
            np.max(
                np.linalg.norm(
                    coordinates - carbon_centroid,
                    axis=1,
                )
            )
        )
        carbon_radius = float(
            np.max(
                np.linalg.norm(
                    coordinates[carbon_mask] - carbon_centroid,
                    axis=1,
                )
            )
        )

        geometry[site] = {
            "carbon_centroid_x_A": float(carbon_centroid[0]),
            "carbon_centroid_y_A": float(carbon_centroid[1]),
            "carbon_centroid_z_A": float(carbon_centroid[2]),
            "all_atom_radius_A": all_atom_radius,
            "carbon_radius_A": carbon_radius,
        }

        manifest.append(
            {
                "frame": FRAME,
                "site": site,
                "bright_root": bright_roots[site],
                "source_output": str(output_path.relative_to(PROJECT_ROOT)),
                "source_gbw": str(gbw_path.relative_to(PROJECT_ROOT)),
                "source_cis": str(cis_path.relative_to(PROJECT_ROOT)),
                "gbw_size_bytes": gbw_path.stat().st_size,
                "cis_size_bytes": cis_path.stat().st_size,
                "pilot_directory": str(site_dir.relative_to(PROJECT_ROOT)),
                "pilot_gbw": str(
                    (site_dir / "pilot.gbw").relative_to(PROJECT_ROOT)
                ),
                "pilot_cis": str(
                    (site_dir / "pilot.cis").relative_to(PROJECT_ROOT)
                ),
                "all_atom_radius_A": all_atom_radius,
                "carbon_radius_A": carbon_radius,
                "status": "READY",
            }
        )

        log(
            f"[{site}] bright=S{bright_roots[site]} "
            f"GBW={gbw_path.stat().st_size} bytes "
            f"CIS={cis_path.stat().st_size} bytes "
            f"radius={all_atom_radius:.6f} A"
        )

    return manifest, geometry


def build_geometry_diagnostic(
    geometry: dict[str, dict[str, float]],
) -> list[dict[str, object]]:
    if not SITE_PAIR_GEOMETRY.is_file():
        raise SystemExit(
            f"Missing site-pair geometry: {SITE_PAIR_GEOMETRY}"
        )

    pair_rows: dict[tuple[str, str], dict[str, str]] = {}

    with SITE_PAIR_GEOMETRY.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        for row in csv.DictReader(handle):
            pair_rows[(row["site_a"], row["site_b"])] = row

    expected = set(combinations(SITES, 2))

    if set(pair_rows) != expected:
        raise RuntimeError(
            f"Unexpected site-pair rows: {sorted(pair_rows)}"
        )

    diagnostics: list[dict[str, object]] = []

    for site_a, site_b in combinations(SITES, 2):
        pair = pair_rows[(site_a, site_b)]
        distance = float(pair["carbon_centroid_distance_A"])
        radius_sum = (
            geometry[site_a]["all_atom_radius_A"]
            + geometry[site_b]["all_atom_radius_A"]
        )
        carbon_radius_sum = (
            geometry[site_a]["carbon_radius_A"]
            + geometry[site_b]["carbon_radius_A"]
        )

        diagnostics.append(
            {
                "site_a": site_a,
                "site_b": site_b,
                "carbon_centroid_distance_A": distance,
                "minimum_all_atom_distance_A": float(
                    pair["minimum_all_atom_distance_A"]
                ),
                "all_atom_radius_a_A": geometry[site_a][
                    "all_atom_radius_A"
                ],
                "all_atom_radius_b_A": geometry[site_b][
                    "all_atom_radius_A"
                ],
                "sum_all_atom_radii_over_centroid_distance": (
                    radius_sum / distance
                ),
                "sum_carbon_radii_over_centroid_distance": (
                    carbon_radius_sum / distance
                ),
                "inverse_distance_rank_metric_A-1": 1.0 / distance,
            }
        )

    diagnostics.sort(
        key=lambda row: (
            -float(
                row[
                    "sum_all_atom_radii_over_centroid_distance"
                ]
            ),
            float(row["carbon_centroid_distance_A"]),
        )
    )

    for rank, row in enumerate(diagnostics, start=1):
        row["finite_size_benchmark_priority_rank"] = rank

    return diagnostics


def write_launcher() -> None:
    text = r'''#!/usr/bin/env bash
set -euo pipefail

SITE="${1:-PYR2}"
ROOT="runs/phase1A/day016_md_bath_extraction/day019_transition_density_pilot"
SITE_DIR="$ROOT/inputs/$SITE"

case "$SITE" in
  PYR2|PYR3|PYR4|PYR5) ;;
  *)
    echo "ERROR: site must be PYR2, PYR3, PYR4, or PYR5" >&2
    exit 2
    ;;
esac

if [[ ! -e "$SITE_DIR/pilot.gbw" || ! -e "$SITE_DIR/pilot.cis" ]]; then
  echo "ERROR: missing pilot.gbw or pilot.cis under $SITE_DIR" >&2
  exit 3
fi

cd "$SITE_DIR"

echo "Launching ORCA transition-density probe for $SITE"
echo "Transcript: orca_plot_probe_${SITE}.typescript"
echo "Select option 7 (CIS/TD-DFT transition densities)."
echo "This is a menu/protocol discovery run; do not generate production cubes yet."

script -q "orca_plot_probe_${SITE}.typescript" \
  orca_plot pilot.gbw -i
'''
    LAUNCHER_SH.write_text(text, encoding="utf-8")
    LAUNCHER_SH.chmod(0o755)


def write_report(
    manifest: list[dict[str, object]],
    diagnostics: list[dict[str, object]],
) -> None:
    nearest = diagnostics[:3]

    with REPORT_MD.open("w", encoding="utf-8") as handle:
        handle.write(
            "# Day019 transition-density pilot preparation\n\n"
        )

        handle.write("## Scope\n\n")
        handle.write(
            "- Prepared the four frame000 bright-like monomer states.\n"
        )
        handle.write(
            "- Verified matching ORCA `.gbw` and `.cis` files.\n"
        )
        handle.write(
            "- Created standardized per-site symlinks named "
            "`pilot.gbw` and `pilot.cis`.\n"
        )
        handle.write(
            "- Ranked all six site pairs by molecular-size/separation "
            "geometry for transition-density benchmarking.\n\n"
        )

        handle.write("## Pilot states\n\n")
        handle.write(
            "| Site | Bright root | GBW bytes | CIS bytes | "
            "All-atom radius (A) | Status |\n"
        )
        handle.write("|---|---:|---:|---:|---:|---|\n")

        for row in manifest:
            handle.write(
                f"| {row['site']} "
                f"| S{row['bright_root']} "
                f"| {row['gbw_size_bytes']} "
                f"| {row['cis_size_bytes']} "
                f"| {float(row['all_atom_radius_A']):.6f} "
                f"| {row['status']} |\n"
            )

        handle.write("\n## Geometry diagnostic\n\n")
        handle.write(
            "| Rank | Pair | Centroid distance (A) | "
            "Minimum atom distance (A) | "
            "(a_i+a_j)/R |\n"
        )
        handle.write("|---:|---|---:|---:|---:|\n")

        for row in diagnostics:
            handle.write(
                f"| {row['finite_size_benchmark_priority_rank']} "
                f"| {row['site_a']}-{row['site_b']} "
                f"| {float(row['carbon_centroid_distance_A']):.6f} "
                f"| {float(row['minimum_all_atom_distance_A']):.6f} "
                f"| {float(row['sum_all_atom_radii_over_centroid_distance']):.6f} |\n"
            )

        handle.write("\n## Priority pairs\n\n")

        for row in nearest:
            handle.write(
                f"- {row['site_a']}-{row['site_b']}: "
                f"R={float(row['carbon_centroid_distance_A']):.3f} A, "
                f"minimum atom distance="
                f"{float(row['minimum_all_atom_distance_A']):.3f} A.\n"
            )

        handle.write("\n## Next controlled action\n\n")
        handle.write(
            "Run `launch_orca_plot_probe.sh PYR2` and select the "
            "CIS/TD-DFT transition-density option. The purpose of the "
            "first probe is to capture the exact ORCA 6.1.1 interactive "
            "prompt sequence before automating all four cubes. Production "
            "generation must use a single documented grid convention and "
            "must preserve the transition-density sign by matching the "
            "cube-derived dipole to the audited ORCA transition dipole.\n"
        )


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    INPUT_ROOT.mkdir(parents=True, exist_ok=True)

    selected = read_selected_frame000_outputs()
    bright_roots = read_bright_roots()

    log("Day019 transition-density pilot preparation")
    manifest, geometry = prepare_inputs(
        selected=selected,
        bright_roots=bright_roots,
    )

    diagnostics = build_geometry_diagnostic(geometry)

    write_csv(MANIFEST_CSV, manifest)
    write_csv(GEOMETRY_CSV, diagnostics)
    write_launcher()
    write_report(manifest, diagnostics)

    log("")
    log("Day019 transition-density pilot preparation completed.")
    log("Pilot states: 4/4 READY")
    log("GBW companions: 4/4")
    log("CIS companions: 4/4")
    log("Geometry pairs ranked: 6/6")
    log(f"Wrote: {OUTPUT_ROOT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
