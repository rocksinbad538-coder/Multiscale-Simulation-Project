#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
from itertools import combinations
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_bright_finite_size_corrected_hamiltonians/"
    "hamiltonian_snapshots_bright4"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_bright_model_readiness_audit"
)

DIAGONAL_TS_CSV = OUTPUT_ROOT / "bright_diagonal_time_series.csv"
ENERGY_STATS_CSV = OUTPUT_ROOT / "bright_energy_statistics.csv"
COUPLING_TS_CSV = OUTPUT_ROOT / "bright_coupling_time_series.csv"
COUPLING_STATS_CSV = OUTPUT_ROOT / "bright_coupling_statistics.csv"
DETUNING_STATS_CSV = OUTPUT_ROOT / "bright_pair_detuning_statistics.csv"
AUTOCORRELATION_CSV = OUTPUT_ROOT / "bright_energy_autocorrelation.csv"
EIGENSTATE_CSV = OUTPUT_ROOT / "bright_eigenstate_statistics.csv"
FRAME_SUMMARY_CSV = OUTPUT_ROOT / "bright_frame_readiness_summary.csv"

CORRELATION_MATRIX_CSV = (
    OUTPUT_ROOT / "bright_energy_correlation_matrix.csv"
)
COVARIANCE_MATRIX_CSV = (
    OUTPUT_ROOT / "bright_energy_covariance_meV2.csv"
)
MEAN_HAMILTONIAN = OUTPUT_ROOT / "mean_bright_hamiltonian_eV.dat"
REPORT_MD = OUTPUT_ROOT / "BRIGHT_MODEL_READINESS_AUDIT_DAY019.md"

EXPECTED_FRAMES = 21
EXPECTED_BASIS = (
    "PYR2_bright",
    "PYR3_bright",
    "PYR4_bright",
    "PYR5_bright",
)

FRAME_SPACING_PS = 5.0
HBAR_MEV_PS = 0.6582119569

FRAME_RE = re.compile(r"frame=(\d+)")
TIME_RE = re.compile(r"time_ps=([0-9.+\-Ee]+)")


def log(message: str = "") -> None:
    print(message, flush=True)


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
        writer.writerow(["state", *labels])

        for label, row in zip(labels, matrix):
            writer.writerow(
                [
                    label,
                    *[
                        f"{float(value):.12g}"
                        for value in row
                    ],
                ]
            )


def parse_snapshot(
    path: Path,
) -> tuple[int, float, tuple[str, ...], np.ndarray]:
    lines = path.read_text(encoding="utf-8").splitlines()

    frame: int | None = None
    time_ps: float | None = None
    basis: tuple[str, ...] | None = None

    for line in lines:
        if not line.startswith("#"):
            continue

        frame_match = FRAME_RE.search(line)
        if frame_match is not None:
            frame = int(frame_match.group(1))

        time_match = TIME_RE.search(line)
        if time_match is not None:
            time_ps = float(time_match.group(1))

        if line.startswith("# basis:"):
            basis = tuple(line.split(":", 1)[1].split())

    if frame is None:
        filename_match = re.search(r"frame(\d{3})", path.name)

        if filename_match is None:
            raise RuntimeError(
                f"Could not identify frame in {path}"
            )

        frame = int(filename_match.group(1))

    if time_ps is None:
        time_ps = frame * FRAME_SPACING_PS

    if basis is None:
        raise RuntimeError(f"Missing basis header in {path}")

    matrix = np.loadtxt(
        path,
        comments="#",
        dtype=np.float64,
    )

    if matrix.shape != (4, 4):
        raise RuntimeError(
            f"Expected 4x4 Hamiltonian in {path}, "
            f"found {matrix.shape}"
        )

    if not np.allclose(
        matrix,
        matrix.T,
        atol=1.0e-12,
        rtol=0.0,
    ):
        raise RuntimeError(
            f"Hamiltonian is not symmetric: {path}"
        )

    if not np.all(np.isfinite(matrix)):
        raise RuntimeError(
            f"Hamiltonian contains nonfinite values: {path}"
        )

    return frame, time_ps, basis, matrix


def normalized_autocorrelation(
    values: np.ndarray,
) -> np.ndarray:
    centered = values - np.mean(values)
    variance = float(np.mean(centered * centered))

    if variance <= 0.0:
        return np.full(values.size, np.nan)

    result = np.zeros(values.size, dtype=np.float64)

    for lag in range(values.size):
        result[lag] = float(
            np.mean(
                centered[: values.size - lag]
                * centered[lag:]
            )
            / variance
        )

    return result


def safe_ratio(
    numerator: np.ndarray,
    denominator: np.ndarray,
) -> np.ndarray:
    result = np.full_like(
        numerator,
        np.nan,
        dtype=np.float64,
    )

    mask = np.abs(denominator) > 1.0e-12
    result[mask] = numerator[mask] / denominator[mask]

    return result


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    snapshot_paths = sorted(
        INPUT_ROOT.glob("H_bright4_tdcac_frame*.dat")
    )

    if len(snapshot_paths) != EXPECTED_FRAMES:
        raise RuntimeError(
            f"Expected {EXPECTED_FRAMES} bright snapshots, "
            f"found {len(snapshot_paths)}"
        )

    frames: list[int] = []
    times_ps: list[float] = []
    matrices: list[np.ndarray] = []

    for path in snapshot_paths:
        frame, time_ps, basis, matrix = parse_snapshot(path)

        if basis != EXPECTED_BASIS:
            raise RuntimeError(
                f"Unexpected basis in {path}: {basis}"
            )

        frames.append(frame)
        times_ps.append(time_ps)
        matrices.append(matrix)

    expected_frame_indices = list(range(EXPECTED_FRAMES))

    if frames != expected_frame_indices:
        raise RuntimeError(
            f"Unexpected frame sequence: {frames}"
        )

    times = np.asarray(times_ps, dtype=np.float64)
    hamiltonians = np.stack(matrices, axis=0)

    diagonal_eV = np.diagonal(
        hamiltonians,
        axis1=1,
        axis2=2,
    )
    diagonal_meV = diagonal_eV * 1000.0
    diagonal_fluctuations_meV = (
        diagonal_meV
        - np.mean(diagonal_meV, axis=0, keepdims=True)
    )

    diagonal_rows: list[dict[str, object]] = []
    energy_stat_rows: list[dict[str, object]] = []
    coupling_rows: list[dict[str, object]] = []
    coupling_stat_rows: list[dict[str, object]] = []
    detuning_stat_rows: list[dict[str, object]] = []
    autocorrelation_rows: list[dict[str, object]] = []
    eigenstate_rows: list[dict[str, object]] = []
    frame_summary_rows: list[dict[str, object]] = []

    for site_index, state_label in enumerate(EXPECTED_BASIS):
        energies = diagonal_eV[:, site_index]
        fluctuations = diagonal_fluctuations_meV[:, site_index]

        for frame, time_ps, energy_eV, fluctuation_meV in zip(
            frames,
            times,
            energies,
            fluctuations,
        ):
            diagonal_rows.append(
                {
                    "frame": frame,
                    "time_ps": time_ps,
                    "state": state_label,
                    "energy_eV": energy_eV,
                    "fluctuation_meV": fluctuation_meV,
                }
            )

        energy_stat_rows.append(
            {
                "state": state_label,
                "mean_energy_eV": float(np.mean(energies)),
                "sd_energy_meV": float(
                    np.std(energies, ddof=0) * 1000.0
                ),
                "minimum_energy_eV": float(np.min(energies)),
                "maximum_energy_eV": float(np.max(energies)),
                "range_energy_meV": float(
                    (np.max(energies) - np.min(energies))
                    * 1000.0
                ),
            }
        )

        autocorrelation = normalized_autocorrelation(
            fluctuations
        )

        for lag, value in enumerate(autocorrelation):
            autocorrelation_rows.append(
                {
                    "state": state_label,
                    "lag_frames": lag,
                    "lag_ps": lag * FRAME_SPACING_PS,
                    "normalized_autocorrelation": value,
                    "n_products": EXPECTED_FRAMES - lag,
                }
            )

    pair_data: dict[
        tuple[int, int],
        dict[str, np.ndarray],
    ] = {}

    for index_a, index_b in combinations(range(4), 2):
        state_a = EXPECTED_BASIS[index_a]
        state_b = EXPECTED_BASIS[index_b]

        coupling_meV = (
            hamiltonians[:, index_a, index_b]
            * 1000.0
        )
        detuning_meV = (
            diagonal_eV[:, index_a]
            - diagonal_eV[:, index_b]
        ) * 1000.0

        absolute_detuning_meV = np.abs(detuning_meV)
        absolute_coupling_meV = np.abs(coupling_meV)

        coupling_to_detuning = safe_ratio(
            absolute_coupling_meV,
            absolute_detuning_meV,
        )

        pair_data[(index_a, index_b)] = {
            "coupling_meV": coupling_meV,
            "detuning_meV": detuning_meV,
            "ratio": coupling_to_detuning,
        }

        for (
            frame,
            time_ps,
            coupling,
            detuning,
            ratio,
        ) in zip(
            frames,
            times,
            coupling_meV,
            detuning_meV,
            coupling_to_detuning,
        ):
            coupling_rows.append(
                {
                    "frame": frame,
                    "time_ps": time_ps,
                    "state_a": state_a,
                    "state_b": state_b,
                    "J_meV": coupling,
                    "absolute_J_meV": abs(coupling),
                    "detuning_meV": detuning,
                    "absolute_detuning_meV": abs(detuning),
                    "absolute_J_over_absolute_detuning": ratio,
                }
            )

        coupling_stat_rows.append(
            {
                "state_a": state_a,
                "state_b": state_b,
                "mean_J_meV": float(np.mean(coupling_meV)),
                "sd_J_meV": float(
                    np.std(coupling_meV, ddof=0)
                ),
                "mean_absolute_J_meV": float(
                    np.mean(absolute_coupling_meV)
                ),
                "maximum_absolute_J_meV": float(
                    np.max(absolute_coupling_meV)
                ),
                "minimum_J_meV": float(np.min(coupling_meV)),
                "maximum_J_meV": float(np.max(coupling_meV)),
                "sign_changes": int(
                    np.count_nonzero(
                        np.sign(coupling_meV[1:])
                        != np.sign(coupling_meV[:-1])
                    )
                ),
            }
        )

        finite_ratios = coupling_to_detuning[
            np.isfinite(coupling_to_detuning)
        ]

        detuning_stat_rows.append(
            {
                "state_a": state_a,
                "state_b": state_b,
                "mean_detuning_meV": float(
                    np.mean(detuning_meV)
                ),
                "sd_detuning_meV": float(
                    np.std(detuning_meV, ddof=0)
                ),
                "mean_absolute_detuning_meV": float(
                    np.mean(absolute_detuning_meV)
                ),
                "minimum_absolute_detuning_meV": float(
                    np.min(absolute_detuning_meV)
                ),
                "maximum_absolute_detuning_meV": float(
                    np.max(absolute_detuning_meV)
                ),
                "maximum_absolute_J_over_detuning": float(
                    np.max(finite_ratios)
                ),
                "mean_absolute_J_over_detuning": float(
                    np.mean(finite_ratios)
                ),
                "frames_ratio_ge_0p1": int(
                    np.count_nonzero(finite_ratios >= 0.1)
                ),
                "frames_ratio_ge_0p5": int(
                    np.count_nonzero(finite_ratios >= 0.5)
                ),
                "frames_ratio_ge_1": int(
                    np.count_nonzero(finite_ratios >= 1.0)
                ),
            }
        )

    all_eigenvalues: list[np.ndarray] = []
    all_participation_ratios: list[np.ndarray] = []
    all_minimum_spacings_meV: list[float] = []

    for frame_index, (
        frame,
        time_ps,
        matrix,
    ) in enumerate(
        zip(frames, times, hamiltonians)
    ):
        eigenvalues, eigenvectors = np.linalg.eigh(matrix)
        populations = eigenvectors * eigenvectors
        participation_ratios = 1.0 / np.sum(
            populations * populations,
            axis=0,
        )

        adjacent_spacings_meV = (
            np.diff(eigenvalues) * 1000.0
        )
        minimum_spacing_meV = float(
            np.min(adjacent_spacings_meV)
        )

        all_eigenvalues.append(eigenvalues)
        all_participation_ratios.append(
            participation_ratios
        )
        all_minimum_spacings_meV.append(
            minimum_spacing_meV
        )

        frame_pair_rows = [
            row
            for row in coupling_rows
            if int(row["frame"]) == frame
        ]

        frame_ratios = np.array(
            [
                float(
                    row[
                        "absolute_J_over_absolute_detuning"
                    ]
                )
                for row in frame_pair_rows
                if np.isfinite(
                    float(
                        row[
                            "absolute_J_over_absolute_detuning"
                        ]
                    )
                )
            ],
            dtype=np.float64,
        )

        frame_couplings = np.array(
            [
                float(row["absolute_J_meV"])
                for row in frame_pair_rows
            ],
            dtype=np.float64,
        )

        frame_detunings = np.array(
            [
                float(row["absolute_detuning_meV"])
                for row in frame_pair_rows
            ],
            dtype=np.float64,
        )

        lowest_state_pyr5_population = float(
            populations[3, 0]
        )

        frame_summary_rows.append(
            {
                "frame": frame,
                "time_ps": time_ps,
                "maximum_absolute_J_meV": float(
                    np.max(frame_couplings)
                ),
                "minimum_absolute_site_detuning_meV": float(
                    np.min(frame_detunings)
                ),
                "maximum_absolute_J_over_detuning": float(
                    np.max(frame_ratios)
                ),
                "minimum_adjacent_eigenvalue_spacing_meV": (
                    minimum_spacing_meV
                ),
                "minimum_participation_ratio": float(
                    np.min(participation_ratios)
                ),
                "maximum_participation_ratio": float(
                    np.max(participation_ratios)
                ),
                "lowest_eigenstate_PYR5_population": (
                    lowest_state_pyr5_population
                ),
            }
        )

        for eigen_index in range(4):
            state_populations = populations[:, eigen_index]
            dominant_index = int(
                np.argmax(state_populations)
            )

            eigenstate_rows.append(
                {
                    "frame": frame,
                    "time_ps": time_ps,
                    "eigen_index": eigen_index + 1,
                    "eigenvalue_eV": eigenvalues[eigen_index],
                    "eigenvalue_relative_to_lowest_meV": (
                        eigenvalues[eigen_index]
                        - eigenvalues[0]
                    )
                    * 1000.0,
                    "participation_ratio": (
                        participation_ratios[eigen_index]
                    ),
                    "dominant_local_state": (
                        EXPECTED_BASIS[dominant_index]
                    ),
                    "dominant_population": float(
                        state_populations[dominant_index]
                    ),
                    "PYR2_population": float(
                        state_populations[0]
                    ),
                    "PYR3_population": float(
                        state_populations[1]
                    ),
                    "PYR4_population": float(
                        state_populations[2]
                    ),
                    "PYR5_population": float(
                        state_populations[3]
                    ),
                }
            )

    energy_correlation = np.corrcoef(
        diagonal_fluctuations_meV.T
    )
    energy_covariance = np.cov(
        diagonal_fluctuations_meV.T,
        ddof=0,
    )

    mean_hamiltonian = np.mean(
        hamiltonians,
        axis=0,
    )

    with MEAN_HAMILTONIAN.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day019 mean four-state bright "
            "TDC-AC-corrected Hamiltonian\n"
        )
        handle.write("# energy units: eV\n")
        handle.write(
            "# basis: " + " ".join(EXPECTED_BASIS) + "\n"
        )

        for row in mean_hamiltonian:
            handle.write(
                " ".join(
                    f"{float(value):.12f}"
                    for value in row
                )
                + "\n"
            )

    write_csv(DIAGONAL_TS_CSV, diagonal_rows)
    write_csv(ENERGY_STATS_CSV, energy_stat_rows)
    write_csv(COUPLING_TS_CSV, coupling_rows)
    write_csv(COUPLING_STATS_CSV, coupling_stat_rows)
    write_csv(DETUNING_STATS_CSV, detuning_stat_rows)
    write_csv(AUTOCORRELATION_CSV, autocorrelation_rows)
    write_csv(EIGENSTATE_CSV, eigenstate_rows)
    write_csv(FRAME_SUMMARY_CSV, frame_summary_rows)

    write_matrix_csv(
        CORRELATION_MATRIX_CSV,
        EXPECTED_BASIS,
        energy_correlation,
    )
    write_matrix_csv(
        COVARIANCE_MATRIX_CSV,
        EXPECTED_BASIS,
        energy_covariance,
    )

    all_coupling_values_meV = np.array(
        [
            float(row["J_meV"])
            for row in coupling_rows
        ],
        dtype=np.float64,
    )

    all_absolute_detunings_meV = np.array(
        [
            float(row["absolute_detuning_meV"])
            for row in coupling_rows
        ],
        dtype=np.float64,
    )

    all_ratios = np.array(
        [
            float(
                row[
                    "absolute_J_over_absolute_detuning"
                ]
            )
            for row in coupling_rows
            if np.isfinite(
                float(
                    row[
                        "absolute_J_over_absolute_detuning"
                    ]
                )
            )
        ],
        dtype=np.float64,
    )

    participation_array = np.stack(
        all_participation_ratios,
        axis=0,
    )
    eigenvalue_array = np.stack(
        all_eigenvalues,
        axis=0,
    )

    maximum_absolute_J_meV = float(
        np.max(np.abs(all_coupling_values_meV))
    )
    coupling_timescale_ps = (
        HBAR_MEV_PS / maximum_absolute_J_meV
    )
    sampling_to_coupling_timescale_ratio = (
        FRAME_SPACING_PS / coupling_timescale_ps
    )

    energy_sd_values = np.array(
        [
            float(row["sd_energy_meV"])
            for row in energy_stat_rows
        ],
        dtype=np.float64,
    )

    minimum_detuning_index = int(
        np.argmin(all_absolute_detunings_meV)
    )
    minimum_detuning_row = coupling_rows[
        minimum_detuning_index
    ]

    maximum_ratio_index = int(np.argmax(all_ratios))
    finite_ratio_rows = [
        row
        for row in coupling_rows
        if np.isfinite(
            float(
                row[
                    "absolute_J_over_absolute_detuning"
                ]
            )
        )
    ]
    maximum_ratio_row = finite_ratio_rows[
        maximum_ratio_index
    ]

    lowest_state_pyr5_populations = np.array(
        [
            float(
                row[
                    "lowest_eigenstate_PYR5_population"
                ]
            )
            for row in frame_summary_rows
        ],
        dtype=np.float64,
    )

    first_lag_autocorrelations: dict[str, float] = {}

    for state in EXPECTED_BASIS:
        matches = [
            row
            for row in autocorrelation_rows
            if row["state"] == state
            and int(row["lag_frames"]) == 1
        ]

        if len(matches) != 1:
            raise RuntimeError(
                f"Missing lag-1 autocorrelation for {state}"
            )

        first_lag_autocorrelations[state] = float(
            matches[0]["normalized_autocorrelation"]
        )

    static_ensemble_ready = (
        len(snapshot_paths) == EXPECTED_FRAMES
        and np.all(np.isfinite(hamiltonians))
    )

    continuous_stochastic_ready = (
        FRAME_SPACING_PS
        <= coupling_timescale_ps / 5.0
    )

    spectral_density_ready = (
        EXPECTED_FRAMES >= 100
        and FRAME_SPACING_PS
        <= coupling_timescale_ps / 5.0
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day019 bright-model readiness audit\n\n"
        )

        handle.write("## Input model\n\n")
        handle.write(
            "- Four bright local states: PYR2, PYR3, "
            "PYR4, and PYR5.\n"
        )
        handle.write(
            "- TDC-AC finite-size-corrected "
            "bright-bright couplings.\n"
        )
        handle.write(
            f"- Snapshots: {len(snapshot_paths)}/"
            f"{EXPECTED_FRAMES}.\n"
        )
        handle.write(
            f"- Sampling interval: "
            f"{FRAME_SPACING_PS:.3f} ps.\n"
        )
        handle.write(
            "- Chromophore geometries are frozen; "
            "the time dependence originates from "
            "the solvent embedding.\n\n"
        )

        handle.write("## Diagonal disorder\n\n")
        handle.write(
            "| State | Mean energy (eV) | SD (meV) | "
            "Range (meV) |\n"
        )
        handle.write("|---|---:|---:|---:|\n")

        for row in energy_stat_rows:
            handle.write(
                f"| {row['state']} "
                f"| {float(row['mean_energy_eV']):.6f} "
                f"| {float(row['sd_energy_meV']):.3f} "
                f"| {float(row['range_energy_meV']):.3f} |\n"
            )

        handle.write("\n## Coupling and detuning scales\n\n")
        handle.write(
            f"- Overall maximum |J|: "
            f"{maximum_absolute_J_meV:.6f} meV.\n"
        )
        handle.write(
            f"- Diagonal-energy SD range: "
            f"{float(np.min(energy_sd_values)):.3f}-"
            f"{float(np.max(energy_sd_values)):.3f} meV.\n"
        )
        handle.write(
            f"- Minimum sampled absolute site detuning: "
            f"{float(minimum_detuning_row['absolute_detuning_meV']):.6f} "
            f"meV at frame "
            f"{int(minimum_detuning_row['frame']):03d} "
            f"for {minimum_detuning_row['state_a']}-"
            f"{minimum_detuning_row['state_b']}.\n"
        )
        handle.write(
            f"- Maximum sampled |J|/|Delta|: "
            f"{float(maximum_ratio_row['absolute_J_over_absolute_detuning']):.6f} "
            f"at frame {int(maximum_ratio_row['frame']):03d} "
            f"for {maximum_ratio_row['state_a']}-"
            f"{maximum_ratio_row['state_b']}.\n\n"
        )

        handle.write("## Eigenstate localization\n\n")
        handle.write(
            f"- Participation-ratio range: "
            f"{float(np.min(participation_array)):.6f}-"
            f"{float(np.max(participation_array)):.6f}.\n"
        )
        handle.write(
            f"- Mean participation ratio: "
            f"{float(np.mean(participation_array)):.6f}.\n"
        )
        handle.write(
            f"- Minimum PYR5 population in the lowest "
            f"eigenstate: "
            f"{float(np.min(lowest_state_pyr5_populations)):.8f}.\n"
        )
        handle.write(
            f"- Minimum adjacent eigenvalue spacing: "
            f"{float(np.min(all_minimum_spacings_meV)):.6f} meV.\n"
        )
        handle.write(
            f"- Full sampled eigenvalue range: "
            f"{float((np.max(eigenvalue_array) - np.min(eigenvalue_array)) * 1000.0):.6f} "
            f"meV.\n\n"
        )

        handle.write("## Time-resolution audit\n\n")
        handle.write(
            f"- Electronic coupling timescale hbar/max|J|: "
            f"{coupling_timescale_ps:.6f} ps.\n"
        )
        handle.write(
            f"- Snapshot interval / coupling timescale: "
            f"{sampling_to_coupling_timescale_ratio:.3f}.\n"
        )

        for state in EXPECTED_BASIS:
            handle.write(
                f"- Lag-1 ({FRAME_SPACING_PS:.1f} ps) "
                f"autocorrelation for {state}: "
                f"{first_lag_autocorrelations[state]:+.6f}.\n"
            )

        handle.write("\n")
        handle.write(
            "- The autocorrelations are descriptive only: "
            "21 points are insufficient for a converged "
            "spectral density.\n"
        )
        handle.write(
            "- The 5 ps sampling interval is substantially "
            "longer than the electronic coupling timescale.\n\n"
        )

        handle.write("## Readiness decision\n\n")
        handle.write(
            f"- Static-disorder ensemble: "
            f"{'READY' if static_ensemble_ready else 'NOT READY'}.\n"
        )
        handle.write(
            f"- Continuous stochastic Hamiltonian propagation: "
            f"{'READY' if continuous_stochastic_ready else 'NOT READY'}.\n"
        )
        handle.write(
            f"- Bath autocorrelation/spectral-density extraction: "
            f"{'READY' if spectral_density_ready else 'NOT READY'}.\n\n"
        )

        handle.write("## Accepted use of the present dataset\n\n")
        handle.write(
            "The 21 Hamiltonians can be used as a "
            "quasi-static disorder ensemble: each snapshot "
            "defines one frozen Hamiltonian realization for "
            "independent coherent propagation, ensemble "
            "averaging, eigenstate analysis, or comparison "
            "with phenomenological dephasing models. They "
            "should not be connected sequentially as a "
            "continuous 5 ps-resolved stochastic trajectory.\n\n"
        )

        handle.write("## Required data for dynamical bath models\n\n")
        handle.write(
            "A defensible time-dependent bath treatment "
            "requires a substantially finer embedding "
            "sampling interval and a longer trajectory. "
            "The required interval must resolve both the "
            "sub-picosecond electronic coupling scale and "
            "the relevant solvent fluctuations. The current "
            "21-point series cannot determine a reliable "
            "spectral density, memory kernel, or microscopic "
            "pure-dephasing rate.\n"
        )

    log("Day019 bright-model readiness audit completed.")
    log(f"Snapshots: {len(snapshot_paths)}/{EXPECTED_FRAMES}")
    log(
        f"Energy-SD range: "
        f"{float(np.min(energy_sd_values)):.3f}-"
        f"{float(np.max(energy_sd_values)):.3f} meV"
    )
    log(
        f"Maximum |J|: "
        f"{maximum_absolute_J_meV:.6f} meV"
    )
    log(
        f"Coupling timescale hbar/max|J|: "
        f"{coupling_timescale_ps:.6f} ps"
    )
    log(
        f"Sampling/coupling-timescale ratio: "
        f"{sampling_to_coupling_timescale_ratio:.3f}"
    )
    log(
        f"Maximum |J|/|detuning|: "
        f"{float(np.max(all_ratios)):.6f}"
    )
    log(
        f"Participation-ratio range: "
        f"{float(np.min(participation_array)):.6f}-"
        f"{float(np.max(participation_array)):.6f}"
    )
    log(
        f"Static-disorder ensemble: "
        f"{'READY' if static_ensemble_ready else 'NOT READY'}"
    )
    log(
        f"Continuous stochastic propagation: "
        f"{'READY' if continuous_stochastic_ready else 'NOT READY'}"
    )
    log(
        f"Spectral-density extraction: "
        f"{'READY' if spectral_density_ready else 'NOT READY'}"
    )
    log(f"Wrote: {OUTPUT_ROOT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
