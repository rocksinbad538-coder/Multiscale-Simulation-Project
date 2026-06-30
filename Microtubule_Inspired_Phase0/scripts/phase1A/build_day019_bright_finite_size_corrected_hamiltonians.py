#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
from itertools import combinations, permutations
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

POINT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_point_transition_dipole_couplings"
)

TDCAC_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_bright_tdc_atomic_charge_couplings"
)

POINT_SNAPSHOT_ROOT = POINT_ROOT / "hamiltonian_snapshots"
TDCAC_COUPLINGS = (
    TDCAC_ROOT / "tdc_atomic_charge_couplings_frame000.csv"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_bright_finite_size_corrected_hamiltonians"
)

FULL8_ROOT = OUTPUT_ROOT / "hamiltonian_snapshots_full8_hybrid"
BRIGHT4_ROOT = OUTPUT_ROOT / "hamiltonian_snapshots_bright4"

FACTORS_CSV = OUTPUT_ROOT / "bright_pair_correction_factors.csv"
FRAME_SUMMARY_CSV = OUTPUT_ROOT / "frame_correction_summary.csv"
PAIR_TIME_SERIES_CSV = OUTPUT_ROOT / "bright_pair_corrected_couplings.csv"
FULL8_EIGEN_CSV = OUTPUT_ROOT / "full8_eigenvalue_comparison.csv"
BRIGHT4_EIGEN_CSV = OUTPUT_ROOT / "bright4_eigenvalue_comparison.csv"
BRIGHT4_OVERLAP_CSV = OUTPUT_ROOT / "bright4_eigenvector_overlap_comparison.csv"
MEAN_DELTA_MATRIX_CSV = OUTPUT_ROOT / "mean_bright_correction_matrix_meV.csv"
MAX_DELTA_MATRIX_CSV = OUTPUT_ROOT / "maximum_absolute_bright_correction_matrix_meV.csv"
REPORT_MD = OUTPUT_ROOT / "BRIGHT_FINITE_SIZE_CORRECTED_HAMILTONIANS_DAY019.md"

SITES = ("PYR2", "PYR3", "PYR4", "PYR5")
SITE_ORDER = {site: index for index, site in enumerate(SITES)}

EXPECTED_FRAMES = 21
FRAME_SPACING_PS = 5.0
MAX_DIAGONAL_SD_MEV = 16.034
MIN_LOCAL_GAP_MEV = 53.0

TIME_RE = re.compile(r"time_ps=([0-9.+\-Ee]+)")
FRAME_RE = re.compile(r"frame=(\d+)")


def log(message: str = "") -> None:
    print(message, flush=True)


def canonical_pair(site_a: str, site_b: str) -> tuple[str, str]:
    if SITE_ORDER[site_a] < SITE_ORDER[site_b]:
        return site_a, site_b
    return site_b, site_a


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
    labels: tuple[str, ...],
    matrix: np.ndarray,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["site", *labels])

        for label, row in zip(labels, matrix):
            writer.writerow(
                [label, *[f"{float(value):.12g}" for value in row]]
            )


def read_correction_factors() -> tuple[
    dict[tuple[str, str], float],
    list[dict[str, object]],
]:
    if not TDCAC_COUPLINGS.is_file():
        raise SystemExit(f"Missing TDC-AC couplings: {TDCAC_COUPLINGS}")

    factors: dict[tuple[str, str], float] = {}
    rows: list[dict[str, object]] = []

    with TDCAC_COUPLINGS.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            site_a = row["site_a"]
            site_b = row["site_b"]
            key = canonical_pair(site_a, site_b)

            point_meV = float(row["saved_point_dipole_J_meV"])
            tdcac_meV = float(row["corrected_TDCAC_J_meV"])
            factor = float(row["TDCAC_over_point_ratio"])

            if abs(point_meV) < 1.0e-14:
                raise RuntimeError(
                    f"Cannot define correction factor for zero coupling: {key}"
                )

            reproduced = tdcac_meV / point_meV

            if abs(reproduced - factor) > 1.0e-10:
                raise RuntimeError(
                    f"Inconsistent correction factor for {key}: "
                    f"CSV={factor}, recomputed={reproduced}"
                )

            factors[key] = factor
            rows.append(
                {
                    "site_a": key[0],
                    "site_b": key[1],
                    "frame000_point_J_meV": point_meV,
                    "frame000_TDCAC_J_meV": tdcac_meV,
                    "correction_factor": factor,
                    "percent_change": 100.0 * (factor - 1.0),
                    "sign_preserved": (
                        np.sign(point_meV) == np.sign(tdcac_meV)
                    ),
                }
            )

    expected = set(combinations(SITES, 2))

    if set(factors) != expected:
        raise RuntimeError(
            f"Expected six pair factors, found {sorted(factors)}."
        )

    return factors, rows


def parse_snapshot(path: Path) -> tuple[
    int,
    float,
    tuple[str, ...],
    np.ndarray,
    list[str],
]:
    if not path.is_file():
        raise SystemExit(f"Missing Hamiltonian snapshot: {path}")

    lines = path.read_text(encoding="utf-8").splitlines()
    header = [line for line in lines if line.startswith("#")]

    frame: int | None = None
    time_ps: float | None = None
    basis: tuple[str, ...] | None = None

    for line in header:
        frame_match = FRAME_RE.search(line)
        if frame_match:
            frame = int(frame_match.group(1))

        time_match = TIME_RE.search(line)
        if time_match:
            time_ps = float(time_match.group(1))

        if line.startswith("# basis:"):
            basis = tuple(line.split(":", 1)[1].split())

    if frame is None:
        match = re.search(r"frame(\d{3})", path.name)
        if not match:
            raise RuntimeError(f"Could not determine frame from {path}")
        frame = int(match.group(1))

    if time_ps is None:
        time_ps = frame * FRAME_SPACING_PS

    if basis is None:
        raise RuntimeError(f"Missing basis comment in {path}")

    matrix = np.loadtxt(path, comments="#", dtype=np.float64)

    if matrix.shape != (8, 8):
        raise RuntimeError(
            f"Expected 8x8 Hamiltonian in {path}, found {matrix.shape}."
        )

    if not np.allclose(matrix, matrix.T, atol=1.0e-12, rtol=0.0):
        raise RuntimeError(f"Hamiltonian is not symmetric: {path}")

    return frame, time_ps, basis, matrix, header


def write_snapshot(
    path: Path,
    matrix: np.ndarray,
    frame: int,
    time_ps: float,
    basis: tuple[str, ...],
    model_description: str,
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(
            f"# Day019 {model_description}, "
            f"frame={frame:03d}, time_ps={time_ps:.6f}\n"
        )
        handle.write(
            "# energy units: eV; matrix is real symmetric\n"
        )
        handle.write(
            "# bright-bright finite-size correction factors were "
            "derived from frame000 TDC-AC / point-dipole ratios\n"
        )
        handle.write(
            "# no dielectric screening is included\n"
        )
        handle.write("# basis: " + " ".join(basis) + "\n")

        for row in matrix:
            handle.write(
                " ".join(f"{float(value):.12f}" for value in row)
                + "\n"
            )


def match_bright_eigenvectors(
    vectors_raw: np.ndarray,
    vectors_corrected: np.ndarray,
) -> tuple[tuple[int, ...], np.ndarray]:
    overlap = np.abs(vectors_raw.T @ vectors_corrected)

    best_permutation: tuple[int, ...] | None = None
    best_score = -np.inf

    for permutation in permutations(range(4)):
        score = sum(
            overlap[index, permutation[index]]
            for index in range(4)
        )

        if score > best_score:
            best_score = score
            best_permutation = permutation

    if best_permutation is None:
        raise RuntimeError("Failed to match bright eigenvectors.")

    matched_overlaps = np.array(
        [
            overlap[index, best_permutation[index]]
            for index in range(4)
        ],
        dtype=np.float64,
    )

    return best_permutation, matched_overlaps


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    FULL8_ROOT.mkdir(parents=True, exist_ok=True)
    BRIGHT4_ROOT.mkdir(parents=True, exist_ok=True)

    factors, factor_rows = read_correction_factors()
    write_csv(FACTORS_CSV, factor_rows)

    snapshot_paths = sorted(
        POINT_SNAPSHOT_ROOT.glob("H_point_dipole_frame*.dat")
    )

    if len(snapshot_paths) != EXPECTED_FRAMES:
        raise RuntimeError(
            f"Expected {EXPECTED_FRAMES} point-dipole snapshots, "
            f"found {len(snapshot_paths)}."
        )

    pair_time_rows: list[dict[str, object]] = []
    frame_rows: list[dict[str, object]] = []
    full8_eigen_rows: list[dict[str, object]] = []
    bright4_eigen_rows: list[dict[str, object]] = []
    bright4_overlap_rows: list[dict[str, object]] = []

    delta_matrices_meV: list[np.ndarray] = []

    log("Day019 finite-size-corrected bright Hamiltonian construction")
    log("Primary model: four bright local states")
    log(
        "Sensitivity model: full eight-state Hamiltonian with only "
        "bright-bright couplings corrected"
    )

    reference_basis: tuple[str, ...] | None = None

    for path in snapshot_paths:
        frame, time_ps, basis, raw_full8, _ = parse_snapshot(path)

        if reference_basis is None:
            reference_basis = basis
        elif basis != reference_basis:
            raise RuntimeError(
                f"Basis changed at frame {frame}: {basis}"
            )

        expected_basis = (
            "PYR2_alternate",
            "PYR2_bright",
            "PYR3_alternate",
            "PYR3_bright",
            "PYR4_alternate",
            "PYR4_bright",
            "PYR5_alternate",
            "PYR5_bright",
        )

        if basis != expected_basis:
            raise RuntimeError(
                f"Unexpected eight-state basis: {basis}"
            )

        bright_indices = {
            site: basis.index(f"{site}_bright")
            for site in SITES
        }
        bright_basis = tuple(f"{site}_bright" for site in SITES)

        corrected_full8 = raw_full8.copy()
        frame_delta_meV = np.zeros((4, 4), dtype=np.float64)

        for site_a, site_b in combinations(SITES, 2):
            key = canonical_pair(site_a, site_b)
            factor = factors[key]
            index_a = bright_indices[site_a]
            index_b = bright_indices[site_b]

            point_eV = raw_full8[index_a, index_b]
            corrected_eV = factor * point_eV
            delta_meV = (corrected_eV - point_eV) * 1000.0

            corrected_full8[index_a, index_b] = corrected_eV
            corrected_full8[index_b, index_a] = corrected_eV

            site_index_a = SITE_ORDER[site_a]
            site_index_b = SITE_ORDER[site_b]
            frame_delta_meV[site_index_a, site_index_b] = delta_meV
            frame_delta_meV[site_index_b, site_index_a] = delta_meV

            pair_time_rows.append(
                {
                    "frame": frame,
                    "time_ps": time_ps,
                    "site_a": site_a,
                    "site_b": site_b,
                    "correction_factor": factor,
                    "point_J_meV": point_eV * 1000.0,
                    "corrected_J_meV": corrected_eV * 1000.0,
                    "delta_J_meV": delta_meV,
                    "absolute_delta_J_meV": abs(delta_meV),
                }
            )

        if not np.allclose(
            np.diag(corrected_full8),
            np.diag(raw_full8),
            atol=0.0,
            rtol=0.0,
        ):
            raise RuntimeError(
                f"Diagonal changed unexpectedly at frame {frame}."
            )

        if not np.allclose(
            corrected_full8,
            corrected_full8.T,
            atol=1.0e-12,
            rtol=0.0,
        ):
            raise RuntimeError(
                f"Corrected Hamiltonian lost symmetry at frame {frame}."
            )

        bright_index_list = [bright_indices[site] for site in SITES]
        raw_bright4 = raw_full8[np.ix_(
            bright_index_list,
            bright_index_list,
        )]
        corrected_bright4 = corrected_full8[np.ix_(
            bright_index_list,
            bright_index_list,
        )]

        raw_full8_eigenvalues = np.linalg.eigvalsh(raw_full8)
        corrected_full8_eigenvalues = np.linalg.eigvalsh(
            corrected_full8
        )
        full8_shift_meV = (
            corrected_full8_eigenvalues
            - raw_full8_eigenvalues
        ) * 1000.0

        raw_bright_eigenvalues, raw_bright_vectors = np.linalg.eigh(
            raw_bright4
        )
        corrected_bright_eigenvalues, corrected_bright_vectors = (
            np.linalg.eigh(corrected_bright4)
        )

        bright_permutation, bright_overlaps = (
            match_bright_eigenvectors(
                raw_bright_vectors,
                corrected_bright_vectors,
            )
        )

        matched_corrected_eigenvalues = np.array(
            [
                corrected_bright_eigenvalues[
                    bright_permutation[index]
                ]
                for index in range(4)
            ],
            dtype=np.float64,
        )
        bright_shift_meV = (
            matched_corrected_eigenvalues
            - raw_bright_eigenvalues
        ) * 1000.0

        full8_path = (
            FULL8_ROOT
            / f"H_full8_hybrid_tdcac_frame{frame:03d}.dat"
        )
        bright4_path = (
            BRIGHT4_ROOT
            / f"H_bright4_tdcac_frame{frame:03d}.dat"
        )

        write_snapshot(
            path=full8_path,
            matrix=corrected_full8,
            frame=frame,
            time_ps=time_ps,
            basis=basis,
            model_description=(
                "full-eight-state hybrid TDC-AC-corrected Hamiltonian"
            ),
        )
        write_snapshot(
            path=bright4_path,
            matrix=corrected_bright4,
            frame=frame,
            time_ps=time_ps,
            basis=bright_basis,
            model_description=(
                "four-state bright TDC-AC-corrected Hamiltonian"
            ),
        )

        for eigen_index, (
            raw_value,
            corrected_value,
            shift,
        ) in enumerate(
            zip(
                raw_full8_eigenvalues,
                corrected_full8_eigenvalues,
                full8_shift_meV,
            ),
            start=1,
        ):
            full8_eigen_rows.append(
                {
                    "frame": frame,
                    "time_ps": time_ps,
                    "eigen_index": eigen_index,
                    "point_eigenvalue_eV": raw_value,
                    "corrected_eigenvalue_eV": corrected_value,
                    "shift_meV": shift,
                    "absolute_shift_meV": abs(shift),
                }
            )

        for raw_index in range(4):
            corrected_index = bright_permutation[raw_index]

            bright4_eigen_rows.append(
                {
                    "frame": frame,
                    "time_ps": time_ps,
                    "raw_eigen_index": raw_index + 1,
                    "matched_corrected_eigen_index": (
                        corrected_index + 1
                    ),
                    "point_eigenvalue_eV": (
                        raw_bright_eigenvalues[raw_index]
                    ),
                    "corrected_eigenvalue_eV": (
                        corrected_bright_eigenvalues[
                            corrected_index
                        ]
                    ),
                    "shift_meV": bright_shift_meV[raw_index],
                    "absolute_shift_meV": abs(
                        bright_shift_meV[raw_index]
                    ),
                }
            )

            bright4_overlap_rows.append(
                {
                    "frame": frame,
                    "time_ps": time_ps,
                    "raw_eigen_index": raw_index + 1,
                    "matched_corrected_eigen_index": (
                        corrected_index + 1
                    ),
                    "absolute_eigenvector_overlap": (
                        bright_overlaps[raw_index]
                    ),
                    "one_minus_overlap": (
                        1.0 - bright_overlaps[raw_index]
                    ),
                }
            )

        all_delta_values = np.array(
            [
                row["delta_J_meV"]
                for row in pair_time_rows
                if int(row["frame"]) == frame
            ],
            dtype=np.float64,
        )

        frame_rows.append(
            {
                "frame": frame,
                "time_ps": time_ps,
                "maximum_absolute_delta_J_meV": float(
                    np.max(np.abs(all_delta_values))
                ),
                "rms_delta_J_meV": float(
                    np.sqrt(np.mean(all_delta_values**2))
                ),
                "maximum_full8_eigenvalue_shift_meV": float(
                    np.max(np.abs(full8_shift_meV))
                ),
                "rms_full8_eigenvalue_shift_meV": float(
                    np.sqrt(np.mean(full8_shift_meV**2))
                ),
                "maximum_bright4_eigenvalue_shift_meV": float(
                    np.max(np.abs(bright_shift_meV))
                ),
                "minimum_matched_bright4_eigenvector_overlap": float(
                    np.min(bright_overlaps)
                ),
                "full8_snapshot": str(
                    full8_path.relative_to(PROJECT_ROOT)
                ),
                "bright4_snapshot": str(
                    bright4_path.relative_to(PROJECT_ROOT)
                ),
            }
        )

        delta_matrices_meV.append(frame_delta_meV)

        log(
            f"[frame {frame:03d}/{EXPECTED_FRAMES - 1:03d}] "
            f"max|delta J|="
            f"{float(np.max(np.abs(all_delta_values))):.6f} meV, "
            f"max bright4 eig shift="
            f"{float(np.max(np.abs(bright_shift_meV))):.6f} meV, "
            f"min eigvec overlap="
            f"{float(np.min(bright_overlaps)):.8f}"
        )

    if reference_basis is None:
        raise RuntimeError("No Hamiltonian snapshots were processed.")

    write_csv(PAIR_TIME_SERIES_CSV, pair_time_rows)
    write_csv(FRAME_SUMMARY_CSV, frame_rows)
    write_csv(FULL8_EIGEN_CSV, full8_eigen_rows)
    write_csv(BRIGHT4_EIGEN_CSV, bright4_eigen_rows)
    write_csv(BRIGHT4_OVERLAP_CSV, bright4_overlap_rows)

    delta_stack = np.stack(delta_matrices_meV, axis=0)
    mean_delta_matrix = np.mean(delta_stack, axis=0)
    max_abs_delta_matrix = np.max(np.abs(delta_stack), axis=0)

    write_matrix_csv(
        MEAN_DELTA_MATRIX_CSV,
        SITES,
        mean_delta_matrix,
    )
    write_matrix_csv(
        MAX_DELTA_MATRIX_CSV,
        SITES,
        max_abs_delta_matrix,
    )

    maximum_absolute_delta_J = max(
        float(row["absolute_delta_J_meV"])
        for row in pair_time_rows
    )
    maximum_full8_shift = max(
        float(row["absolute_shift_meV"])
        for row in full8_eigen_rows
    )
    maximum_bright4_shift = max(
        float(row["absolute_shift_meV"])
        for row in bright4_eigen_rows
    )
    minimum_bright4_overlap = min(
        float(row["absolute_eigenvector_overlap"])
        for row in bright4_overlap_rows
    )

    nearest_pairs = {
        ("PYR2", "PYR3"),
        ("PYR3", "PYR4"),
        ("PYR4", "PYR5"),
    }
    distant_pairs = set(combinations(SITES, 2)) - nearest_pairs

    nearest_factor_changes = [
        abs(factors[pair] - 1.0)
        for pair in nearest_pairs
    ]
    distant_factor_changes = [
        abs(factors[pair] - 1.0)
        for pair in distant_pairs
    ]

    with REPORT_MD.open("w", encoding="utf-8") as handle:
        handle.write(
            "# Day019 bright finite-size-corrected Hamiltonians\n\n"
        )

        handle.write("## Construction\n\n")
        handle.write(
            "- Six static correction factors were defined as "
            "`frame000 TDC-AC / frame000 point-dipole`.\n"
        )
        handle.write(
            "- Each factor was applied to the corresponding bright-bright "
            "point-dipole coupling in all 21 solvent frames.\n"
        )
        handle.write(
            "- The four-state bright model is the primary finite-size-"
            "corrected coupling model.\n"
        )
        handle.write(
            "- A full eight-state hybrid sensitivity model was also "
            "written; only its bright-bright elements are corrected. "
            "Mixed and alternate-state couplings remain the original "
            "point-dipole values and are not transition-density "
            "benchmarked.\n\n"
        )

        handle.write("## Pair correction factors\n\n")
        handle.write(
            "| Pair | Point frame000 (meV) | TDC-AC frame000 (meV) | "
            "Factor | Change | Sign preserved |\n"
        )
        handle.write("|---|---:|---:|---:|---:|---|\n")

        for row in factor_rows:
            handle.write(
                f"| {row['site_a']}-{row['site_b']} "
                f"| {float(row['frame000_point_J_meV']):+.6f} "
                f"| {float(row['frame000_TDCAC_J_meV']):+.6f} "
                f"| {float(row['correction_factor']):.6f} "
                f"| {float(row['percent_change']):+.3f}% "
                f"| {'yes' if row['sign_preserved'] else 'no'} |\n"
            )

        handle.write("\n## Dynamic impact across 21 frames\n\n")
        handle.write(
            f"- Hamiltonian snapshots generated: "
            f"{len(frame_rows)}/{EXPECTED_FRAMES}\n"
        )
        handle.write(
            f"- Maximum |coupling correction|: "
            f"{maximum_absolute_delta_J:.6f} meV\n"
        )
        handle.write(
            f"- Maximum full-eight-state eigenvalue shift: "
            f"{maximum_full8_shift:.6f} meV\n"
        )
        handle.write(
            f"- Maximum four-state bright eigenvalue shift: "
            f"{maximum_bright4_shift:.6f} meV\n"
        )
        handle.write(
            f"- Minimum matched bright eigenvector overlap: "
            f"{minimum_bright4_overlap:.10f}\n"
        )
        handle.write(
            f"- Maximum correction / maximum diagonal-energy SD: "
            f"{maximum_absolute_delta_J / MAX_DIAGONAL_SD_MEV:.6f}\n"
        )
        handle.write(
            f"- Maximum correction / minimum local S1-S2 gap: "
            f"{maximum_absolute_delta_J / MIN_LOCAL_GAP_MEV:.6f}\n\n"
        )

        handle.write("## Pair-class summary\n\n")
        handle.write(
            "- Mean absolute factor change, nearest pairs: "
            f"{100.0 * float(np.mean(nearest_factor_changes)):.4f}%\n"
        )
        handle.write(
            "- Mean absolute factor change, distant pairs: "
            f"{100.0 * float(np.mean(distant_factor_changes)):.4f}%\n\n"
        )

        handle.write("## Interpretation boundary\n\n")
        handle.write(
            "The pair-specific factors transfer a static frame000 "
            "finite-size correction to all solvent frames. This is "
            "consistent with the frozen chromophore geometries and the "
            "observed stability of the bright transition-dipole directions, "
            "but it does not constitute an embedded transition-density "
            "calculation for every frame. The four-state bright Hamiltonian "
            "is therefore the defensible corrected control model. The "
            "eight-state hybrid is retained only for sensitivity analysis "
            "until alternate-state transition densities and their gauge are "
            "treated explicitly. No dielectric screening is included.\n"
        )

    log("")
    log(
        "Day019 bright finite-size-corrected Hamiltonians completed."
    )
    log(f"Correction factors: {len(factors)}/6")
    log(f"Full-eight-state snapshots: {len(frame_rows)}/21")
    log(f"Bright-four-state snapshots: {len(frame_rows)}/21")
    log(
        f"Maximum |delta J|: "
        f"{maximum_absolute_delta_J:.6f} meV"
    )
    log(
        f"Maximum full8 eigenvalue shift: "
        f"{maximum_full8_shift:.6f} meV"
    )
    log(
        f"Maximum bright4 eigenvalue shift: "
        f"{maximum_bright4_shift:.6f} meV"
    )
    log(
        f"Minimum matched bright4 eigenvector overlap: "
        f"{minimum_bright4_overlap:.10f}"
    )
    log(f"Wrote: {OUTPUT_ROOT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
