#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from statistics import mean, pstdev

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

SOURCE_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_full_two_state_manifold_analysis"
)

FRAME_SUMMARY = SOURCE_ROOT / "two_state_frame_summary.csv"

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_eight_state_diagonal_model"
)

LONG_CSV = OUTPUT_ROOT / "eight_state_diagonal_long.csv"
WIDE_CSV = OUTPUT_ROOT / "eight_state_diagonal_timeseries.csv"
STATE_STATS_CSV = OUTPUT_ROOT / "eight_state_state_statistics.csv"
COV_CSV = OUTPUT_ROOT / "eight_state_covariance_meV2.csv"
CORR_CSV = OUTPUT_ROOT / "eight_state_correlation.csv"
COMMON_DIFF_CSV = OUTPUT_ROOT / "within_site_common_differential_modes.csv"
SNAPSHOT_DIR = OUTPUT_ROOT / "diagonal_snapshots"
BASIS_TXT = OUTPUT_ROOT / "EIGHT_STATE_BASIS_ORDER_DAY019.txt"
REPORT_MD = OUTPUT_ROOT / "EIGHT_STATE_DIAGONAL_MODEL_DAY019.md"

SITES = ("PYR2", "PYR3", "PYR4", "PYR5")
FAMILIES = ("alternate_like", "bright_like")

EXPECTED_BRIGHT_ROOT = {
    "PYR2": 2,
    "PYR3": 2,
    "PYR4": 2,
    "PYR5": 1,
}

EXPECTED_FRAMES = tuple(range(21))


def log(message: str = "") -> None:
    print(message, flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build the Day019 character-indexed eight-state diagonal "
            "Hamiltonian time series from the 84 embedded S1/S2 results."
        )
    )
    parser.add_argument(
        "--frame-spacing-ps",
        type=float,
        default=5.0,
        help="Time separation between consecutive extracted MD frames.",
    )
    parser.add_argument(
        "--energy-zero",
        choices=("absolute", "global_minimum", "per_frame_minimum"),
        default="global_minimum",
        help=(
            "Energy-zero convention used in the exported Hamiltonian "
            "snapshot matrices. Raw absolute energies are always retained "
            "in the CSV outputs."
        ),
    )
    return parser.parse_args()


def read_frame_summary() -> list[dict[str, str]]:
    if not FRAME_SUMMARY.is_file():
        raise SystemExit(f"Missing input table: {FRAME_SUMMARY}")

    with FRAME_SUMMARY.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if len(rows) != 84:
        raise RuntimeError(
            f"Expected 84 frame/site rows, found {len(rows)}."
        )

    return rows


def bool_from_csv(value: str) -> bool:
    return value.strip().lower() == "true"


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    if a.shape != b.shape:
        raise RuntimeError("Pearson vectors have different shapes.")

    a_centered = a - np.mean(a)
    b_centered = b - np.mean(b)

    denominator = math.sqrt(
        float(np.dot(a_centered, a_centered))
        * float(np.dot(b_centered, b_centered))
    )

    if denominator == 0.0:
        return float("nan")

    return float(np.dot(a_centered, b_centered) / denominator)


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
    fieldnames: list[str] | None = None,
) -> None:
    if not rows:
        raise RuntimeError(f"No rows to write: {path}")

    if fieldnames is None:
        fieldnames = list(rows[0].keys())

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_matrix_csv(
    path: Path,
    labels: list[str],
    matrix: np.ndarray,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["state", *labels])

        for label, row in zip(labels, matrix):
            writer.writerow([label, *[f"{value:.12g}" for value in row]])


def main() -> None:
    args = parse_args()

    if args.frame_spacing_ps <= 0.0:
        raise SystemExit("--frame-spacing-ps must be positive.")

    rows = read_frame_summary()
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)

    lookup: dict[tuple[int, str], dict[str, str]] = {}

    for row in rows:
        frame = int(row["frame"])
        site = row["site"]

        key = (frame, site)

        if key in lookup:
            raise RuntimeError(f"Duplicate frame/site row: {key}")

        if frame not in EXPECTED_FRAMES:
            raise RuntimeError(f"Unexpected frame index: {frame}")

        if site not in SITES:
            raise RuntimeError(f"Unexpected site: {site}")

        if not bool_from_csv(row["mapping_preserved"]):
            raise RuntimeError(
                f"Root mapping was not preserved for frame {frame}, {site}."
            )

        expected_root = EXPECTED_BRIGHT_ROOT[site]

        if int(row["bright_root"]) != expected_root:
            raise RuntimeError(
                f"Unexpected bright root for frame {frame}, {site}: "
                f"S{row['bright_root']} instead of S{expected_root}."
            )

        lookup[key] = row

    expected_keys = {
        (frame, site)
        for frame in EXPECTED_FRAMES
        for site in SITES
    }

    if set(lookup) != expected_keys:
        missing = sorted(expected_keys - set(lookup))
        extra = sorted(set(lookup) - expected_keys)
        raise RuntimeError(
            f"Frame/site coverage mismatch. Missing={missing}; extra={extra}"
        )

    basis: list[dict[str, object]] = []

    for site in SITES:
        bright_root = EXPECTED_BRIGHT_ROOT[site]
        alternate_root = 1 if bright_root == 2 else 2

        basis.extend(
            [
                {
                    "state_index": len(basis),
                    "state_label": f"{site}_alternate",
                    "site": site,
                    "family": "alternate_like",
                    "root": alternate_root,
                },
                {
                    "state_index": len(basis) + 1,
                    "state_label": f"{site}_bright",
                    "site": site,
                    "family": "bright_like",
                    "root": bright_root,
                },
            ]
        )

    labels = [str(item["state_label"]) for item in basis]

    with BASIS_TXT.open("w", encoding="utf-8") as handle:
        handle.write("# Day019 eight-state basis order\n")
        handle.write("# index label site family root\n")

        for item in basis:
            handle.write(
                f"{item['state_index']} "
                f"{item['state_label']} "
                f"{item['site']} "
                f"{item['family']} "
                f"S{item['root']}\n"
            )

    absolute = np.zeros((21, 8), dtype=np.float64)
    oscillator_strengths = np.zeros((21, 8), dtype=np.float64)
    long_rows: list[dict[str, object]] = []

    for frame in EXPECTED_FRAMES:
        time_ps = frame * args.frame_spacing_ps

        for item in basis:
            site = str(item["site"])
            family = str(item["family"])
            index = int(item["state_index"])
            row = lookup[(frame, site)]

            if family == "bright_like":
                energy = float(row["bright_energy_eV"])
                fosc = float(row["bright_fosc"])
            else:
                energy = float(row["alternate_energy_eV"])
                fosc = float(row["alternate_fosc"])

            absolute[frame, index] = energy
            oscillator_strengths[frame, index] = fosc

            long_rows.append(
                {
                    "frame": frame,
                    "time_ps": time_ps,
                    "state_index": index,
                    "state_label": item["state_label"],
                    "site": site,
                    "family": family,
                    "root": item["root"],
                    "energy_eV": energy,
                    "oscillator_strength": fosc,
                    "S1_S2_gap_meV": float(row["S1_S2_gap_meV"]),
                }
            )

    means = np.mean(absolute, axis=0)
    fluctuations_meV = (absolute - means) * 1000.0

    covariance = np.cov(
        fluctuations_meV,
        rowvar=False,
        bias=True,
    )

    correlation = np.corrcoef(
        fluctuations_meV,
        rowvar=False,
    )

    if covariance.shape != (8, 8):
        raise RuntimeError("Unexpected covariance-matrix shape.")

    if correlation.shape != (8, 8):
        raise RuntimeError("Unexpected correlation-matrix shape.")

    wide_rows: list[dict[str, object]] = []

    for frame in EXPECTED_FRAMES:
        row: dict[str, object] = {
            "frame": frame,
            "time_ps": frame * args.frame_spacing_ps,
        }

        for index, label in enumerate(labels):
            row[f"{label}_energy_eV"] = absolute[frame, index]
            row[f"{label}_deltaE_meV"] = fluctuations_meV[frame, index]
            row[f"{label}_fosc"] = oscillator_strengths[frame, index]

        wide_rows.append(row)

    write_csv(WIDE_CSV, wide_rows)
    write_csv(LONG_CSV, long_rows)
    write_matrix_csv(COV_CSV, labels, covariance)
    write_matrix_csv(CORR_CSV, labels, correlation)

    stats_rows: list[dict[str, object]] = []

    for item in basis:
        index = int(item["state_index"])
        energies = absolute[:, index]
        fluctuations = fluctuations_meV[:, index]
        foscs = oscillator_strengths[:, index]

        stats_rows.append(
            {
                "state_index": index,
                "state_label": item["state_label"],
                "site": item["site"],
                "family": item["family"],
                "root": item["root"],
                "n_frames": len(energies),
                "mean_energy_eV": float(np.mean(energies)),
                "sd_energy_meV": float(np.std(fluctuations)),
                "minimum_energy_eV": float(np.min(energies)),
                "maximum_energy_eV": float(np.max(energies)),
                "energy_range_meV": float(
                    (np.max(energies) - np.min(energies)) * 1000.0
                ),
                "minimum_deltaE_meV": float(np.min(fluctuations)),
                "maximum_deltaE_meV": float(np.max(fluctuations)),
                "mean_fosc": float(np.mean(foscs)),
                "sd_fosc": float(np.std(foscs)),
                "minimum_fosc": float(np.min(foscs)),
                "maximum_fosc": float(np.max(foscs)),
            }
        )

    expected_sd_meV = np.std(fluctuations_meV, axis=0)
    exported_sd_meV = np.array(
        [float(row["sd_energy_meV"]) for row in stats_rows],
        dtype=np.float64,
    )

    if not np.allclose(
        exported_sd_meV,
        expected_sd_meV,
        atol=1.0e-12,
        rtol=1.0e-12,
    ):
        raise RuntimeError(
            "State-statistics SD unit validation failed: "
            f"exported={exported_sd_meV}, expected={expected_sd_meV}"
        )

    write_csv(STATE_STATS_CSV, stats_rows)

    common_diff_rows: list[dict[str, object]] = []

    for site_index, site in enumerate(SITES):
        alternate_index = site_index * 2
        bright_index = alternate_index + 1

        alternate_energy = absolute[:, alternate_index]
        bright_energy = absolute[:, bright_index]

        common_energy = 0.5 * (
            alternate_energy + bright_energy
        )

        signed_gap = (
            bright_energy - alternate_energy
        )

        common_fluctuation_meV = (
            common_energy - np.mean(common_energy)
        ) * 1000.0

        gap_fluctuation_meV = (
            signed_gap - np.mean(signed_gap)
        ) * 1000.0

        for frame in EXPECTED_FRAMES:
            common_diff_rows.append(
                {
                    "frame": frame,
                    "time_ps": frame * args.frame_spacing_ps,
                    "site": site,
                    "alternate_energy_eV": alternate_energy[frame],
                    "bright_energy_eV": bright_energy[frame],
                    "common_mode_energy_eV": common_energy[frame],
                    "signed_bright_minus_alternate_gap_meV": (
                        signed_gap[frame] * 1000.0
                    ),
                    "common_mode_fluctuation_meV": (
                        common_fluctuation_meV[frame]
                    ),
                    "gap_fluctuation_meV": gap_fluctuation_meV[frame],
                }
            )

    write_csv(COMMON_DIFF_CSV, common_diff_rows)

    if args.energy_zero == "absolute":
        snapshot_values = absolute.copy()
        zero_description = "absolute TDDFT excitation energies"
    elif args.energy_zero == "global_minimum":
        global_minimum = float(np.min(absolute))
        snapshot_values = absolute - global_minimum
        zero_description = (
            f"global minimum energy shifted to zero "
            f"(E0={global_minimum:.9f} eV)"
        )
    else:
        frame_minima = np.min(absolute, axis=1, keepdims=True)
        snapshot_values = absolute - frame_minima
        zero_description = "minimum state energy shifted to zero in each frame"

    for frame in EXPECTED_FRAMES:
        matrix = np.diag(snapshot_values[frame])

        path = SNAPSHOT_DIR / f"Hdiag_frame{frame:03d}.dat"

        with path.open("w", encoding="utf-8") as handle:
            handle.write(
                f"# Day019 diagonal-only eight-state Hamiltonian, "
                f"frame={frame:03d}, "
                f"time_ps={frame * args.frame_spacing_ps:.6f}\n"
            )
            handle.write(f"# energy_zero: {zero_description}\n")
            handle.write("# basis: " + " ".join(labels) + "\n")
            np.savetxt(handle, matrix, fmt="%.12f")

    bright_indices = [
        int(item["state_index"])
        for item in basis
        if item["family"] == "bright_like"
    ]

    alternate_indices = [
        int(item["state_index"])
        for item in basis
        if item["family"] == "alternate_like"
    ]

    bright_sds = [
        float(np.std(fluctuations_meV[:, index]))
        for index in bright_indices
    ]

    alternate_sds = [
        float(np.std(fluctuations_meV[:, index]))
        for index in alternate_indices
    ]

    within_site_correlations: dict[str, float] = {}

    for site_index, site in enumerate(SITES):
        alternate_index = site_index * 2
        bright_index = alternate_index + 1

        within_site_correlations[site] = pearson(
            fluctuations_meV[:, alternate_index],
            fluctuations_meV[:, bright_index],
        )

    maximum_absolute_offdiag_correlation = max(
        abs(float(correlation[i, j]))
        for i in range(8)
        for j in range(8)
        if i != j
    )

    with REPORT_MD.open("w", encoding="utf-8") as handle:
        handle.write(
            "# Day019 character-indexed eight-state diagonal model\n\n"
        )

        handle.write("## Basis definition\n\n")
        handle.write(
            "| Index | Label | Site | Family | TDDFT root |\n"
        )
        handle.write("|---:|---|---|---|---:|\n")

        for item in basis:
            handle.write(
                f"| {item['state_index']} "
                f"| `{item['state_label']}` "
                f"| {item['site']} "
                f"| {item['family']} "
                f"| S{item['root']} |\n"
            )

        handle.write("\n")
        handle.write(
            "The labels follow electronic character, not a globally fixed "
            "root number. PYR5 therefore uses S1 as the bright-like state, "
            "whereas PYR2-PYR4 use S2.\n\n"
        )

        handle.write("## Validation\n\n")
        handle.write("- Frame/site rows consumed: 84/84\n")
        handle.write("- Character-indexed state observations: 168/168\n")
        handle.write("- Frames represented: 21/21\n")
        handle.write(
            f"- Frame spacing: {args.frame_spacing_ps:.6f} ps\n"
        )
        handle.write("- Bright-root mapping preserved: 84/84\n")
        handle.write(
            "- State-statistics energy SD values were validated and exported "
            "in meV.\n"
        )
        handle.write(
            "- Hamiltonian snapshots are diagonal-only; all off-diagonal "
            "couplings are intentionally zero placeholders.\n\n"
        )

        handle.write("## State statistics\n\n")
        handle.write(
            "| State | Mean energy (eV) | SD (meV) | Range (meV) | "
            "Mean fosc |\n"
        )
        handle.write("|---|---:|---:|---:|---:|\n")

        for row in stats_rows:
            handle.write(
                f"| `{row['state_label']}` "
                f"| {float(row['mean_energy_eV']):.6f} "
                f"| {float(row['sd_energy_meV']):.3f} "
                f"| {float(row['energy_range_meV']):.3f} "
                f"| {float(row['mean_fosc']):.6f} |\n"
            )

        handle.write("\n## Aggregate fluctuation structure\n\n")
        handle.write(
            f"- Bright-like energy SD range: "
            f"{min(bright_sds):.3f}-{max(bright_sds):.3f} meV\n"
        )
        handle.write(
            f"- Alternate-like energy SD range: "
            f"{min(alternate_sds):.3f}-{max(alternate_sds):.3f} meV\n"
        )
        handle.write(
            f"- Maximum absolute off-diagonal energy correlation: "
            f"{maximum_absolute_offdiag_correlation:.6f}\n"
        )

        for site in SITES:
            handle.write(
                f"- {site} alternate/bright fluctuation correlation: "
                f"{within_site_correlations[site]:.6f}\n"
            )

        handle.write("\n## Hamiltonian status\n\n")
        handle.write(
            "The resulting matrices provide the complete time-dependent "
            "diagonal component of the eight-state Hamiltonian. They are not "
            "yet the final excitonic Hamiltonian because intersite couplings "
            "and possible same-site interstate couplings have not been "
            "computed. No dynamical propagation should treat the current "
            "zero off-diagonal entries as physical coupling estimates.\n"
        )

    log("Day019 eight-state diagonal model completed.")
    log("Frame/site rows: 84/84")
    log("State observations: 168/168")
    log("Hamiltonian snapshots: 21/21")
    log("State-statistics SD unit validation: PASS")
    log("Basis states: 8")
    log(f"Frame spacing: {args.frame_spacing_ps:.3f} ps")
    log(f"Energy-zero convention: {args.energy_zero}")
    log(
        "Mean-energy range: "
        f"{float(np.min(means)):.6f}-{float(np.max(means)):.6f} eV"
    )
    log(
        "Energy-SD range: "
        f"{float(np.min(np.std(fluctuations_meV, axis=0))):.3f}-"
        f"{float(np.max(np.std(fluctuations_meV, axis=0))):.3f} meV"
    )
    log(f"Wrote: {OUTPUT_ROOT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
