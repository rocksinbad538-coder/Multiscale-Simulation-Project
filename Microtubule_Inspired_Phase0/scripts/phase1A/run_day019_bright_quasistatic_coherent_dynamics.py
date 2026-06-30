#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

POINT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_point_transition_dipole_couplings/"
    "hamiltonian_snapshots"
)

CORRECTED_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_bright_finite_size_corrected_hamiltonians/"
    "hamiltonian_snapshots_bright4"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_bright_quasistatic_coherent_dynamics"
)

ENSEMBLE_CSV = OUTPUT_ROOT / "ensemble_population_statistics.csv"
MODEL_DIFFERENCE_CSV = (
    OUTPUT_ROOT / "point_vs_tdcac_population_difference.csv"
)
FRAME_METRICS_CSV = (
    OUTPUT_ROOT / "frame_coherent_dynamics_metrics.csv"
)
SUMMARY_CSV = OUTPUT_ROOT / "initial_state_dynamics_summary.csv"
DIAGONAL_ENSEMBLE_CSV = (
    OUTPUT_ROOT / "diagonal_ensemble_populations.csv"
)
NPZ_PATH = OUTPUT_ROOT / "population_trajectories.npz"
REPORT_MD = (
    OUTPUT_ROOT
    / "BRIGHT_QUASISTATIC_COHERENT_DYNAMICS_DAY019.md"
)

EXPECTED_FRAMES = 21

SITES = (
    "PYR2_bright",
    "PYR3_bright",
    "PYR4_bright",
    "PYR5_bright",
)

MODELS = (
    "point_dipole",
    "tdcac_corrected",
)

T_MAX_PS = 20.0
DT_PS = 0.005
HBAR_EV_PS = 6.582119569e-4

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
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_snapshot(
    path: Path,
    expected_shape: tuple[int, int],
) -> tuple[int, float, tuple[str, ...], np.ndarray]:
    lines = path.read_text(
        encoding="utf-8",
    ).splitlines()

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
            basis = tuple(
                line.split(":", 1)[1].split()
            )

    if frame is None:
        filename_match = re.search(
            r"frame(\d{3})",
            path.name,
        )

        if filename_match is None:
            raise RuntimeError(
                f"Could not determine frame from {path}"
            )

        frame = int(filename_match.group(1))

    if time_ps is None:
        time_ps = frame * 5.0

    if basis is None:
        raise RuntimeError(
            f"Missing basis header in {path}"
        )

    matrix = np.loadtxt(
        path,
        comments="#",
        dtype=np.float64,
    )

    if matrix.shape != expected_shape:
        raise RuntimeError(
            f"Expected matrix shape {expected_shape} "
            f"in {path}, found {matrix.shape}"
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
            f"Nonfinite Hamiltonian values: {path}"
        )

    return frame, time_ps, basis, matrix


def read_hamiltonians(
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    point_paths = sorted(
        POINT_ROOT.glob(
            "H_point_dipole_frame*.dat"
        )
    )

    corrected_paths = sorted(
        CORRECTED_ROOT.glob(
            "H_bright4_tdcac_frame*.dat"
        )
    )

    if len(point_paths) != EXPECTED_FRAMES:
        raise RuntimeError(
            f"Expected {EXPECTED_FRAMES} point snapshots, "
            f"found {len(point_paths)}"
        )

    if len(corrected_paths) != EXPECTED_FRAMES:
        raise RuntimeError(
            f"Expected {EXPECTED_FRAMES} corrected snapshots, "
            f"found {len(corrected_paths)}"
        )

    point_by_frame: dict[
        int,
        tuple[float, np.ndarray],
    ] = {}

    corrected_by_frame: dict[
        int,
        tuple[float, np.ndarray],
    ] = {}

    expected_full_basis = (
        "PYR2_alternate",
        "PYR2_bright",
        "PYR3_alternate",
        "PYR3_bright",
        "PYR4_alternate",
        "PYR4_bright",
        "PYR5_alternate",
        "PYR5_bright",
    )

    for path in point_paths:
        frame, time_ps, basis, matrix = (
            parse_snapshot(
                path,
                expected_shape=(8, 8),
            )
        )

        if basis != expected_full_basis:
            raise RuntimeError(
                "Unexpected point-dipole basis in "
                f"{path}: {basis}"
            )

        bright_indices = [
            basis.index(site)
            for site in SITES
        ]

        bright_matrix = matrix[
            np.ix_(
                bright_indices,
                bright_indices,
            )
        ]

        point_by_frame[frame] = (
            time_ps,
            bright_matrix,
        )

    for path in corrected_paths:
        frame, time_ps, basis, matrix = (
            parse_snapshot(
                path,
                expected_shape=(4, 4),
            )
        )

        if basis != SITES:
            raise RuntimeError(
                f"Unexpected corrected basis "
                f"in {path}: {basis}"
            )

        corrected_by_frame[frame] = (
            time_ps,
            matrix,
        )

    expected_frames = set(
        range(EXPECTED_FRAMES)
    )

    if set(point_by_frame) != expected_frames:
        raise RuntimeError(
            "Unexpected point frame set: "
            f"{sorted(point_by_frame)}"
        )

    if set(corrected_by_frame) != expected_frames:
        raise RuntimeError(
            "Unexpected corrected frame set: "
            f"{sorted(corrected_by_frame)}"
        )

    times_ps: list[float] = []
    point_matrices: list[np.ndarray] = []
    corrected_matrices: list[np.ndarray] = []

    for frame in range(EXPECTED_FRAMES):
        point_time, point_matrix = (
            point_by_frame[frame]
        )

        corrected_time, corrected_matrix = (
            corrected_by_frame[frame]
        )

        if abs(
            point_time - corrected_time
        ) > 1.0e-10:
            raise RuntimeError(
                f"Time mismatch at frame {frame}: "
                f"{point_time} vs {corrected_time}"
            )

        if not np.allclose(
            np.diag(point_matrix),
            np.diag(corrected_matrix),
            atol=1.0e-12,
            rtol=0.0,
        ):
            raise RuntimeError(
                "Diagonal mismatch between models "
                f"at frame {frame}"
            )

        times_ps.append(point_time)
        point_matrices.append(point_matrix)
        corrected_matrices.append(
            corrected_matrix
        )

    return (
        np.asarray(
            times_ps,
            dtype=np.float64,
        ),
        np.stack(
            point_matrices,
            axis=0,
        ),
        np.stack(
            corrected_matrices,
            axis=0,
        ),
    )


def propagate_exact(
    hamiltonian_eV: np.ndarray,
    initial_index: int,
    times_ps: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    # A scalar energy shift affects only the global phase.
    shifted = (
        hamiltonian_eV
        - np.eye(4)
        * float(
            np.trace(hamiltonian_eV) / 4.0
        )
    )

    eigenvalues, eigenvectors = np.linalg.eigh(
        shifted
    )

    initial = np.zeros(
        4,
        dtype=np.complex128,
    )
    initial[initial_index] = 1.0 + 0.0j

    coefficients = (
        eigenvectors.T.conj()
        @ initial
    )

    phases = np.exp(
        -1j
        * times_ps[:, None]
        * eigenvalues[None, :]
        / HBAR_EV_PS
    )

    amplitudes = (
        phases
        * coefficients[None, :]
    ) @ eigenvectors.T

    populations = np.abs(amplitudes) ** 2
    norms = np.sum(
        populations,
        axis=1,
    )

    initial_eigenstate_weights = (
        np.abs(
            eigenvectors[initial_index, :]
        ) ** 2
    )

    diagonal_ensemble = np.sum(
        initial_eigenstate_weights[None, :]
        * np.abs(eigenvectors) ** 2,
        axis=1,
    )

    return (
        populations,
        norms,
        diagonal_ensemble,
    )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    (
        snapshot_times_ps,
        point_hamiltonians,
        corrected_hamiltonians,
    ) = read_hamiltonians()

    times_ps = np.arange(
        0.0,
        T_MAX_PS + 0.5 * DT_PS,
        DT_PS,
        dtype=np.float64,
    )

    if abs(
        times_ps[-1] - T_MAX_PS
    ) > 1.0e-12:
        raise RuntimeError(
            "Time grid does not terminate "
            "at T_MAX_PS"
        )

    model_hamiltonians = {
        "point_dipole": point_hamiltonians,
        "tdcac_corrected": (
            corrected_hamiltonians
        ),
    }

    trajectories: dict[
        tuple[str, int],
        np.ndarray,
    ] = {}

    diagonal_ensembles: dict[
        tuple[str, int],
        np.ndarray,
    ] = {}

    frame_metric_rows: list[
        dict[str, object]
    ] = []

    diagonal_ensemble_rows: list[
        dict[str, object]
    ] = []

    maximum_norm_error = 0.0

    log(
        "Day019 quasi-static coherent "
        "bright-state dynamics"
    )
    log(
        f"Exact propagation: {times_ps.size} "
        f"points from 0 to {T_MAX_PS:.3f} ps, "
        f"dt={DT_PS:.3f} ps"
    )
    log(
        "Each of the 21 solvent frames "
        "is propagated independently."
    )

    for model_name in MODELS:
        hamiltonians = (
            model_hamiltonians[model_name]
        )

        for initial_index, initial_state in enumerate(
            SITES
        ):
            population_stack = np.empty(
                (
                    EXPECTED_FRAMES,
                    times_ps.size,
                    4,
                ),
                dtype=np.float64,
            )

            diagonal_stack = np.empty(
                (
                    EXPECTED_FRAMES,
                    4,
                ),
                dtype=np.float64,
            )

            for frame in range(
                EXPECTED_FRAMES
            ):
                (
                    populations,
                    norms,
                    diagonal_ensemble,
                ) = propagate_exact(
                    hamiltonian_eV=(
                        hamiltonians[frame]
                    ),
                    initial_index=initial_index,
                    times_ps=times_ps,
                )

                norm_error = float(
                    np.max(
                        np.abs(norms - 1.0)
                    )
                )

                maximum_norm_error = max(
                    maximum_norm_error,
                    norm_error,
                )

                if norm_error > 1.0e-10:
                    raise RuntimeError(
                        "Norm failure for "
                        f"{model_name}, "
                        f"{initial_state}, "
                        f"frame {frame}: "
                        f"{norm_error:.3e}"
                    )

                population_stack[frame] = (
                    populations
                )
                diagonal_stack[frame] = (
                    diagonal_ensemble
                )

                survival = populations[
                    :,
                    initial_index,
                ]

                minimum_survival_index = int(
                    np.argmin(survival)
                )

                metric_row: dict[
                    str,
                    object,
                ] = {
                    "model": model_name,
                    "frame": frame,
                    "snapshot_time_ps": (
                        snapshot_times_ps[frame]
                    ),
                    "initial_state": (
                        initial_state
                    ),
                    "minimum_survival_probability": (
                        float(
                            survival[
                                minimum_survival_index
                            ]
                        )
                    ),
                    "time_of_minimum_survival_ps": (
                        float(
                            times_ps[
                                minimum_survival_index
                            ]
                        )
                    ),
                    "maximum_norm_error": (
                        norm_error
                    ),
                }

                for (
                    target_index,
                    target_state,
                ) in enumerate(SITES):
                    target_population = (
                        populations[
                            :,
                            target_index,
                        ]
                    )

                    maximum_index = int(
                        np.argmax(
                            target_population
                        )
                    )

                    metric_row[
                        "maximum_population_"
                        f"{target_state}"
                    ] = float(
                        target_population[
                            maximum_index
                        ]
                    )

                    metric_row[
                        "time_of_maximum_"
                        f"{target_state}_ps"
                    ] = float(
                        times_ps[
                            maximum_index
                        ]
                    )

                    metric_row[
                        "time_average_"
                        f"{target_state}"
                    ] = float(
                        np.mean(
                            target_population
                        )
                    )

                    metric_row[
                        "diagonal_ensemble_"
                        f"{target_state}"
                    ] = float(
                        diagonal_ensemble[
                            target_index
                        ]
                    )

                    diagonal_ensemble_rows.append(
                        {
                            "model": model_name,
                            "frame": frame,
                            "snapshot_time_ps": (
                                snapshot_times_ps[
                                    frame
                                ]
                            ),
                            "initial_state": (
                                initial_state
                            ),
                            "target_state": (
                                target_state
                            ),
                            "population": float(
                                diagonal_ensemble[
                                    target_index
                                ]
                            ),
                        }
                    )

                frame_metric_rows.append(
                    metric_row
                )

            trajectories[
                (
                    model_name,
                    initial_index,
                )
            ] = population_stack

            diagonal_ensembles[
                (
                    model_name,
                    initial_index,
                )
            ] = diagonal_stack

            log(
                f"[{model_name}] "
                f"initial={initial_state}: "
                f"completed "
                f"{EXPECTED_FRAMES}/"
                f"{EXPECTED_FRAMES}"
            )

    ensemble_rows: list[
        dict[str, object]
    ] = []

    summary_rows: list[
        dict[str, object]
    ] = []

    difference_rows: list[
        dict[str, object]
    ] = []

    for model_name in MODELS:
        for initial_index, initial_state in enumerate(
            SITES
        ):
            stack = trajectories[
                (
                    model_name,
                    initial_index,
                )
            ]

            mean = np.mean(
                stack,
                axis=0,
            )
            sd = np.std(
                stack,
                axis=0,
                ddof=0,
            )
            q05 = np.quantile(
                stack,
                0.05,
                axis=0,
            )
            q50 = np.quantile(
                stack,
                0.50,
                axis=0,
            )
            q95 = np.quantile(
                stack,
                0.95,
                axis=0,
            )

            for (
                time_index,
                time_ps,
            ) in enumerate(times_ps):
                row: dict[str, object] = {
                    "model": model_name,
                    "initial_state": (
                        initial_state
                    ),
                    "time_ps": float(
                        time_ps
                    ),
                }

                for (
                    target_index,
                    target_state,
                ) in enumerate(SITES):
                    row[
                        f"mean_{target_state}"
                    ] = float(
                        mean[
                            time_index,
                            target_index,
                        ]
                    )

                    row[
                        f"sd_{target_state}"
                    ] = float(
                        sd[
                            time_index,
                            target_index,
                        ]
                    )

                    row[
                        f"q05_{target_state}"
                    ] = float(
                        q05[
                            time_index,
                            target_index,
                        ]
                    )

                    row[
                        f"q50_{target_state}"
                    ] = float(
                        q50[
                            time_index,
                            target_index,
                        ]
                    )

                    row[
                        f"q95_{target_state}"
                    ] = float(
                        q95[
                            time_index,
                            target_index,
                        ]
                    )

                ensemble_rows.append(row)

            survival_mean = mean[
                :,
                initial_index,
            ]

            minimum_survival_index = int(
                np.argmin(
                    survival_mean
                )
            )

            summary_row: dict[
                str,
                object,
            ] = {
                "model": model_name,
                "initial_state": (
                    initial_state
                ),
                "ensemble_minimum_survival_probability": (
                    float(
                        survival_mean[
                            minimum_survival_index
                        ]
                    )
                ),
                "time_of_ensemble_minimum_survival_ps": (
                    float(
                        times_ps[
                            minimum_survival_index
                        ]
                    )
                ),
                "mean_frame_minimum_survival_probability": (
                    float(
                        np.mean(
                            np.min(
                                stack[
                                    :,
                                    :,
                                    initial_index,
                                ],
                                axis=1,
                            )
                        )
                    )
                ),
                "maximum_norm_error": (
                    maximum_norm_error
                ),
            }

            diagonal_mean = np.mean(
                diagonal_ensembles[
                    (
                        model_name,
                        initial_index,
                    )
                ],
                axis=0,
            )

            for (
                target_index,
                target_state,
            ) in enumerate(SITES):
                target_mean = mean[
                    :,
                    target_index,
                ]

                maximum_index = int(
                    np.argmax(
                        target_mean
                    )
                )

                summary_row[
                    "maximum_ensemble_mean_"
                    f"{target_state}"
                ] = float(
                    target_mean[
                        maximum_index
                    ]
                )

                summary_row[
                    "time_of_maximum_ensemble_mean_"
                    f"{target_state}_ps"
                ] = float(
                    times_ps[
                        maximum_index
                    ]
                )

                summary_row[
                    "time_average_ensemble_mean_"
                    f"{target_state}"
                ] = float(
                    np.mean(
                        target_mean
                    )
                )

                summary_row[
                    "mean_diagonal_ensemble_"
                    f"{target_state}"
                ] = float(
                    diagonal_mean[
                        target_index
                    ]
                )

            summary_rows.append(
                summary_row
            )

    global_max_population_difference = 0.0

    global_difference_context: (
        tuple[str, str, float] | None
    ) = None

    for initial_index, initial_state in enumerate(
        SITES
    ):
        point_mean = np.mean(
            trajectories[
                (
                    "point_dipole",
                    initial_index,
                )
            ],
            axis=0,
        )

        corrected_mean = np.mean(
            trajectories[
                (
                    "tdcac_corrected",
                    initial_index,
                )
            ],
            axis=0,
        )

        difference = (
            corrected_mean - point_mean
        )

        for (
            time_index,
            time_ps,
        ) in enumerate(times_ps):
            row: dict[str, object] = {
                "initial_state": (
                    initial_state
                ),
                "time_ps": float(
                    time_ps
                ),
            }

            for (
                target_index,
                target_state,
            ) in enumerate(SITES):
                value = float(
                    difference[
                        time_index,
                        target_index,
                    ]
                )

                row[
                    f"delta_{target_state}"
                ] = value

                row[
                    "absolute_delta_"
                    f"{target_state}"
                ] = abs(value)

                if (
                    abs(value)
                    > global_max_population_difference
                ):
                    global_max_population_difference = (
                        abs(value)
                    )

                    global_difference_context = (
                        initial_state,
                        target_state,
                        float(time_ps),
                    )

            difference_rows.append(row)

    npz_payload: dict[
        str,
        np.ndarray,
    ] = {
        "times_ps": times_ps,
        "snapshot_times_ps": (
            snapshot_times_ps
        ),
        "site_labels": np.asarray(
            SITES
        ),
    }

    for model_name in MODELS:
        for initial_index, initial_state in enumerate(
            SITES
        ):
            key = (
                f"{model_name}__"
                f"initial_{initial_state}"
            )

            npz_payload[key] = trajectories[
                (
                    model_name,
                    initial_index,
                )
            ]

    np.savez_compressed(
        NPZ_PATH,
        **npz_payload,
    )

    write_csv(
        ENSEMBLE_CSV,
        ensemble_rows,
    )

    write_csv(
        MODEL_DIFFERENCE_CSV,
        difference_rows,
    )

    write_csv(
        FRAME_METRICS_CSV,
        frame_metric_rows,
    )

    write_csv(
        SUMMARY_CSV,
        summary_rows,
    )

    write_csv(
        DIAGONAL_ENSEMBLE_CSV,
        diagonal_ensemble_rows,
    )

    def summary_for(
        model: str,
        initial: str,
    ) -> dict[str, object]:
        matches = [
            row
            for row in summary_rows
            if row["model"] == model
            and row["initial_state"] == initial
        ]

        if len(matches) != 1:
            raise RuntimeError(
                "Could not locate summary for "
                f"{model}/{initial}"
            )

        return matches[0]

    corrected_high_energy_pyr5_maxima = []
    corrected_high_energy_pyr5_diagonal = []

    for initial_state in SITES[:3]:
        row = summary_for(
            "tdcac_corrected",
            initial_state,
        )

        corrected_high_energy_pyr5_maxima.append(
            float(
                row[
                    "maximum_ensemble_mean_"
                    "PYR5_bright"
                ]
            )
        )

        corrected_high_energy_pyr5_diagonal.append(
            float(
                row[
                    "mean_diagonal_ensemble_"
                    "PYR5_bright"
                ]
            )
        )

    maximum_single_frame_pyr5_transfer = 0.0

    maximum_single_frame_pyr5_context: (
        tuple[str, int] | None
    ) = None

    for initial_index, initial_state in enumerate(
        SITES[:3]
    ):
        stack = trajectories[
            (
                "tdcac_corrected",
                initial_index,
            )
        ]

        frame_maxima = np.max(
            stack[
                :,
                :,
                3,
            ],
            axis=1,
        )

        frame_index = int(
            np.argmax(
                frame_maxima
            )
        )

        value = float(
            frame_maxima[
                frame_index
            ]
        )

        if (
            value
            > maximum_single_frame_pyr5_transfer
        ):
            maximum_single_frame_pyr5_transfer = (
                value
            )

            maximum_single_frame_pyr5_context = (
                initial_state,
                frame_index,
            )

    if global_difference_context is None:
        raise RuntimeError(
            "Population difference context "
            "was not set"
        )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day019 quasi-static coherent "
            "bright-state dynamics\n\n"
        )

        handle.write("## Protocol\n\n")

        handle.write(
            "- Twenty-one solvent snapshots were "
            "propagated independently as frozen "
            "Hamiltonian realizations.\n"
        )

        handle.write(
            "- Exact unitary propagation was used; "
            "no numerical time integrator or "
            "interpolation between solvent frames "
            "was introduced.\n"
        )

        handle.write(
            f"- Propagation interval: "
            f"0-{T_MAX_PS:.3f} ps.\n"
        )

        handle.write(
            f"- Output interval: "
            f"{DT_PS:.3f} ps.\n"
        )

        handle.write(
            "- Initial conditions: one localized "
            "excitation on each of the four bright "
            "local states.\n"
        )

        handle.write(
            "- Compared models: point-transition-"
            "dipole and TDC-AC finite-size-corrected "
            "bright Hamiltonians.\n\n"
        )

        handle.write(
            "## Numerical validation\n\n"
        )

        handle.write(
            f"- Maximum population-norm error: "
            f"{maximum_norm_error:.3e}.\n"
        )

        handle.write(
            "- Frames per model and initial "
            f"condition: {EXPECTED_FRAMES}/"
            f"{EXPECTED_FRAMES}.\n"
        )

        handle.write(
            "- Maximum absolute difference between "
            "ensemble-mean point and corrected "
            f"populations: "
            f"{global_max_population_difference:.6e}.\n"
        )

        handle.write(
            "- Difference maximum context: initial "
            f"{global_difference_context[0]}, target "
            f"{global_difference_context[1]}, time "
            f"{global_difference_context[2]:.3f} ps.\n\n"
        )

        handle.write(
            "## Corrected-model ensemble summary\n\n"
        )

        handle.write(
            "| Initial state | Minimum ensemble "
            "survival | Time (ps) | Mean diagonal-"
            "ensemble survival |\n"
        )

        handle.write(
            "|---|---:|---:|---:|\n"
        )

        for initial_state in SITES:
            row = summary_for(
                "tdcac_corrected",
                initial_state,
            )

            handle.write(
                f"| {initial_state} "
                f"| {float(row['ensemble_minimum_survival_probability']):.6f} "
                f"| {float(row['time_of_ensemble_minimum_survival_ps']):.3f} "
                f"| {float(row[f'mean_diagonal_ensemble_{initial_state}']):.6f} |\n"
            )

        handle.write(
            "\n## PYR5 coherent-accessibility audit\n\n"
        )

        handle.write(
            "- Maximum ensemble-mean PYR5 "
            "population from a PYR2/PYR3/PYR4 "
            "localized initial excitation: "
            f"{max(corrected_high_energy_pyr5_maxima):.6e}.\n"
        )

        handle.write(
            "- Maximum mean diagonal-ensemble PYR5 "
            "population from those three initial "
            "conditions: "
            f"{max(corrected_high_energy_pyr5_diagonal):.6e}.\n"
        )

        if (
            maximum_single_frame_pyr5_context
            is not None
        ):
            handle.write(
                "- Maximum PYR5 population in any "
                "individual snapshot from a high-"
                "energy initial state: "
                f"{maximum_single_frame_pyr5_transfer:.6e} "
                f"(initial "
                f"{maximum_single_frame_pyr5_context[0]}, "
                f"frame "
                f"{maximum_single_frame_pyr5_context[1]:03d}).\n"
            )

        handle.write(
            "\nThe large PYR5 energy offset suppresses "
            "direct coherent transfer into PYR5 in the "
            "present closed-system model. This result "
            "does not exclude bath-assisted downhill "
            "relaxation, which is absent from unitary "
            "propagation.\n\n"
        )

        handle.write(
            "## Interpretation boundary\n\n"
        )

        handle.write(
            "Damping of ensemble-averaged oscillations "
            "in this analysis is inhomogeneous dephasing "
            "caused by averaging over static Hamiltonian "
            "realizations. It is not a microscopic "
            "decoherence rate. The calculation contains "
            "no population relaxation, pure-dephasing "
            "operator, spectral density, or sequential "
            "5 ps solvent dynamics. Near-resonant "
            "PYR2-PYR3-PYR4 frames can support transient "
            "coherent mixing, whereas PYR5 remains "
            "energetically isolated in the closed "
            "bright-state model.\n"
        )

    log("")
    log(
        "Day019 quasi-static coherent "
        "dynamics completed."
    )
    log(
        f"Models: {len(MODELS)}/2"
    )
    log(
        f"Initial states: {len(SITES)}/4"
    )
    log(
        f"Frames per condition: "
        f"{EXPECTED_FRAMES}/21"
    )
    log(
        f"Time points: {times_ps.size}"
    )
    log(
        f"Maximum norm error: "
        f"{maximum_norm_error:.3e}"
    )
    log(
        "Maximum point-vs-corrected ensemble "
        "population difference: "
        f"{global_max_population_difference:.6e}"
    )
    log(
        "Maximum ensemble-mean PYR5 population "
        "from PYR2/PYR3/PYR4: "
        f"{max(corrected_high_energy_pyr5_maxima):.6e}"
    )
    log(
        f"Wrote: "
        f"{OUTPUT_ROOT.relative_to(PROJECT_ROOT)}"
    )


if __name__ == "__main__":
    main()
