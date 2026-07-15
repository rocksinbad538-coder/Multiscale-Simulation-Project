#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np

try:
    import scipy
    from scipy.sparse import csr_matrix
    from scipy.sparse.linalg import expm_multiply
except ImportError as exc:
    raise SystemExit(
        "SciPy is required for the exact Liouvillian propagation. "
        "Activate the project environment and install scipy if necessary."
    ) from exc


PROJECT_ROOT = Path(__file__).resolve().parents[2]

HAMILTONIAN_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_bright_finite_size_corrected_hamiltonians/"
    "hamiltonian_snapshots_bright4"
)

COHERENT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_bright_quasistatic_coherent_dynamics"
)

COHERENT_NPZ = COHERENT_ROOT / "population_trajectories.npz"

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day020_bright_haken_strobl_high_dephasing"
)

ENSEMBLE_CSV = OUTPUT_ROOT / "ensemble_population_dephasing_statistics.csv"
SUMMARY_CSV = OUTPUT_ROOT / "initial_state_dephasing_summary.csv"
FRAME_METRICS_CSV = OUTPUT_ROOT / "frame_dephasing_metrics.csv"
GAMMA_AGGREGATE_CSV = OUTPUT_ROOT / "gamma_aggregate_summary.csv"
CROSSOVER_CSV = OUTPUT_ROOT / "dephasing_crossover_summary.csv"
VALIDATION_CSV = OUTPUT_ROOT / "numerical_validation.csv"
NPZ_PATH = OUTPUT_ROOT / "dephasing_population_trajectories.npz"
REPORT_MD = OUTPUT_ROOT / "HAKEN_STROBL_HIGH_DEPHASING_DAY020.md"

EXPECTED_FRAMES = 21
N_STATES = 4

SITES = (
    "PYR2_bright",
    "PYR3_bright",
    "PYR4_bright",
    "PYR5_bright",
)

HIGH_ENERGY_INDICES = (0, 1, 2)
PYR5_INDEX = 3

GAMMA_PHI_PS_INV = np.array(
    [0.0, 0.05, 0.10, 0.20, 0.50, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0],
    dtype=np.float64,
)

T_MAX_PS = 20.0
DT_PS = 0.01
HBAR_EV_PS = 6.582119569e-4

TRACE_TOL = 1.0e-10
HERMITICITY_TOL = 1.0e-10
POSITIVITY_TOL = 1.0e-10
GAMMA0_MATCH_TOL = 1.0e-9

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
            raise RuntimeError(f"Could not determine frame from {path}")

        frame = int(filename_match.group(1))

    if time_ps is None:
        time_ps = frame * 5.0

    if basis is None:
        raise RuntimeError(f"Missing basis header in {path}")

    matrix = np.loadtxt(
        path,
        comments="#",
        dtype=np.float64,
    )

    if matrix.shape != (N_STATES, N_STATES):
        raise RuntimeError(
            f"Expected {N_STATES}x{N_STATES} matrix in {path}, "
            f"found {matrix.shape}"
        )

    if not np.allclose(
        matrix,
        matrix.T,
        atol=1.0e-12,
        rtol=0.0,
    ):
        raise RuntimeError(f"Hamiltonian is not symmetric: {path}")

    if not np.all(np.isfinite(matrix)):
        raise RuntimeError(f"Nonfinite Hamiltonian values: {path}")

    return frame, time_ps, basis, matrix


def read_hamiltonians() -> tuple[np.ndarray, np.ndarray]:
    paths = sorted(
        HAMILTONIAN_ROOT.glob("H_bright4_tdcac_frame*.dat")
    )

    if len(paths) != EXPECTED_FRAMES:
        raise RuntimeError(
            f"Expected {EXPECTED_FRAMES} Hamiltonians, found {len(paths)}"
        )

    frame_records: dict[int, tuple[float, np.ndarray]] = {}

    for path in paths:
        frame, time_ps, basis, matrix = parse_snapshot(path)

        if basis != SITES:
            raise RuntimeError(
                f"Unexpected bright-state basis in {path}: {basis}"
            )

        frame_records[frame] = (time_ps, matrix)

    expected_frames = set(range(EXPECTED_FRAMES))

    if set(frame_records) != expected_frames:
        raise RuntimeError(
            f"Unexpected frame set: {sorted(frame_records)}"
        )

    snapshot_times: list[float] = []
    matrices: list[np.ndarray] = []

    for frame in range(EXPECTED_FRAMES):
        time_ps, matrix = frame_records[frame]
        snapshot_times.append(time_ps)
        matrices.append(matrix)

    return (
        np.asarray(snapshot_times, dtype=np.float64),
        np.stack(matrices, axis=0),
    )


def build_dephasing_superoperator() -> np.ndarray:
    identity = np.eye(N_STATES, dtype=np.complex128)
    dissipator = np.zeros(
        (N_STATES**2, N_STATES**2),
        dtype=np.complex128,
    )

    for site_index in range(N_STATES):
        projector = np.zeros(
            (N_STATES, N_STATES),
            dtype=np.complex128,
        )
        projector[site_index, site_index] = 1.0

        dissipator += (
            np.kron(projector.conj(), projector)
            - 0.5 * np.kron(identity, projector)
            - 0.5 * np.kron(projector.T, identity)
        )

    return dissipator


def validate_dephasing_superoperator(
    dissipator: np.ndarray,
) -> None:
    test_matrix = np.zeros(
        (N_STATES, N_STATES),
        dtype=np.complex128,
    )
    test_matrix[0, 1] = 1.0

    test_vector = test_matrix.reshape(
        N_STATES**2,
        order="F",
    )

    result = dissipator @ test_vector
    expected = -test_vector

    if not np.allclose(
        result,
        expected,
        atol=1.0e-12,
        rtol=0.0,
    ):
        raise RuntimeError(
            "The pure-dephasing superoperator failed the analytical "
            "off-diagonal decay test."
        )


def build_liouvillian(
    hamiltonian_eV: np.ndarray,
    gamma_phi_ps_inv: float,
    dephasing_superoperator: np.ndarray,
) -> np.ndarray:
    identity = np.eye(N_STATES, dtype=np.complex128)

    coherent = (
        -1j
        / HBAR_EV_PS
        * (
            np.kron(identity, hamiltonian_eV)
            - np.kron(hamiltonian_eV.T, identity)
        )
    )

    return (
        coherent
        + gamma_phi_ps_inv * dephasing_superoperator
    )


def build_initial_density_vectors() -> np.ndarray:
    vectors = np.zeros(
        (N_STATES**2, N_STATES),
        dtype=np.complex128,
    )

    for initial_index in range(N_STATES):
        diagonal_vector_index = (
            initial_index
            + initial_index * N_STATES
        )
        vectors[diagonal_vector_index, initial_index] = 1.0

    return vectors


def matrix_from_vector(vector: np.ndarray) -> np.ndarray:
    return vector.reshape(
        (N_STATES, N_STATES),
        order="F",
    )


def time_slug(value: float) -> str:
    text = f"{value:g}"
    return text.replace(".", "p")


def nearest_time_index(
    times_ps: np.ndarray,
    requested_time_ps: float,
) -> int:
    return int(np.argmin(np.abs(times_ps - requested_time_ps)))


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    snapshot_times_ps, hamiltonians = read_hamiltonians()

    times_ps = np.arange(
        0.0,
        T_MAX_PS + 0.5 * DT_PS,
        DT_PS,
        dtype=np.float64,
    )

    if abs(times_ps[-1] - T_MAX_PS) > 1.0e-12:
        raise RuntimeError(
            "The propagation grid does not terminate at T_MAX_PS."
        )

    n_gamma = GAMMA_PHI_PS_INV.size
    n_times = times_ps.size

    populations = np.empty(
        (
            n_gamma,
            N_STATES,
            EXPECTED_FRAMES,
            n_times,
            N_STATES,
        ),
        dtype=np.float64,
    )

    coherence_l1 = np.empty(
        (
            n_gamma,
            N_STATES,
            EXPECTED_FRAMES,
            n_times,
        ),
        dtype=np.float64,
    )

    diagonal_indices = np.array(
        [
            index + index * N_STATES
            for index in range(N_STATES)
        ],
        dtype=np.int64,
    )

    offdiagonal_indices = np.array(
        [
            row + column * N_STATES
            for column in range(N_STATES)
            for row in range(N_STATES)
            if row != column
        ],
        dtype=np.int64,
    )

    validation_time_indices = np.unique(
        np.linspace(
            0,
            n_times - 1,
            21,
            dtype=np.int64,
        )
    )

    dephasing_superoperator = build_dephasing_superoperator()
    validate_dephasing_superoperator(dephasing_superoperator)

    initial_vectors = build_initial_density_vectors()

    maximum_trace_error = 0.0
    maximum_hermiticity_error = 0.0
    minimum_sampled_density_eigenvalue = np.inf
    minimum_raw_population = np.inf
    maximum_raw_population = -np.inf

    frame_metric_rows: list[dict[str, object]] = []

    log("Day020 Haken-Strobl high-dephasing sensitivity")
    log(
        f"Gamma values: {', '.join(f'{value:g}' for value in GAMMA_PHI_PS_INV)} ps^-1"
    )
    log(
        f"Exact Liouvillian propagation: {n_times} points, "
        f"0-{T_MAX_PS:.3f} ps, dt={DT_PS:.3f} ps"
    )
    log(
        "All 21 solvent Hamiltonians are propagated independently."
    )

    for gamma_index, gamma_phi in enumerate(GAMMA_PHI_PS_INV):
        for frame in range(EXPECTED_FRAMES):
            liouvillian = build_liouvillian(
                hamiltonian_eV=hamiltonians[frame],
                gamma_phi_ps_inv=float(gamma_phi),
                dephasing_superoperator=dephasing_superoperator,
            )

            evolution = expm_multiply(
                csr_matrix(liouvillian),
                initial_vectors,
                start=0.0,
                stop=T_MAX_PS,
                num=n_times,
                endpoint=True,
            )

            if evolution.shape != (
                n_times,
                N_STATES**2,
                N_STATES,
            ):
                raise RuntimeError(
                    f"Unexpected expm_multiply shape: {evolution.shape}"
                )

            for initial_index, initial_state in enumerate(SITES):
                states = evolution[:, :, initial_index]

                raw_populations = np.real(
                    states[:, diagonal_indices]
                )

                minimum_raw_population = min(
                    minimum_raw_population,
                    float(np.min(raw_populations)),
                )
                maximum_raw_population = max(
                    maximum_raw_population,
                    float(np.max(raw_populations)),
                )

                if np.min(raw_populations) < -POSITIVITY_TOL:
                    raise RuntimeError(
                        f"Negative population below tolerance for "
                        f"gamma={gamma_phi}, frame={frame}, "
                        f"initial={initial_state}: "
                        f"{np.min(raw_populations):.3e}"
                    )

                if np.max(raw_populations) > 1.0 + POSITIVITY_TOL:
                    raise RuntimeError(
                        f"Population above unity for "
                        f"gamma={gamma_phi}, frame={frame}, "
                        f"initial={initial_state}: "
                        f"{np.max(raw_populations):.3e}"
                    )

                clean_populations = np.clip(
                    raw_populations,
                    0.0,
                    1.0,
                )

                trace_values = np.sum(
                    states[:, diagonal_indices],
                    axis=1,
                )

                trace_error = float(
                    np.max(
                        np.abs(trace_values - 1.0)
                    )
                )

                maximum_trace_error = max(
                    maximum_trace_error,
                    trace_error,
                )

                if trace_error > TRACE_TOL:
                    raise RuntimeError(
                        f"Trace failure for gamma={gamma_phi}, "
                        f"frame={frame}, initial={initial_state}: "
                        f"{trace_error:.3e}"
                    )

                coherence = np.sum(
                    np.abs(
                        states[:, offdiagonal_indices]
                    ),
                    axis=1,
                )

                populations[
                    gamma_index,
                    initial_index,
                    frame,
                ] = clean_populations.astype(np.float64)

                coherence_l1[
                    gamma_index,
                    initial_index,
                    frame,
                ] = coherence.astype(np.float64)

                minimum_eigenvalue_this_trajectory = np.inf
                maximum_hermiticity_this_trajectory = 0.0

                for time_index in validation_time_indices:
                    density = matrix_from_vector(
                        states[time_index]
                    )

                    hermiticity_error = float(
                        np.max(
                            np.abs(
                                density
                                - density.conj().T
                            )
                        )
                    )

                    maximum_hermiticity_error = max(
                        maximum_hermiticity_error,
                        hermiticity_error,
                    )
                    maximum_hermiticity_this_trajectory = max(
                        maximum_hermiticity_this_trajectory,
                        hermiticity_error,
                    )

                    hermitian_density = 0.5 * (
                        density + density.conj().T
                    )

                    minimum_eigenvalue = float(
                        np.min(
                            np.linalg.eigvalsh(
                                hermitian_density
                            )
                        )
                    )

                    minimum_sampled_density_eigenvalue = min(
                        minimum_sampled_density_eigenvalue,
                        minimum_eigenvalue,
                    )
                    minimum_eigenvalue_this_trajectory = min(
                        minimum_eigenvalue_this_trajectory,
                        minimum_eigenvalue,
                    )

                if (
                    maximum_hermiticity_this_trajectory
                    > HERMITICITY_TOL
                ):
                    raise RuntimeError(
                        f"Hermiticity failure for gamma={gamma_phi}, "
                        f"frame={frame}, initial={initial_state}: "
                        f"{maximum_hermiticity_this_trajectory:.3e}"
                    )

                if (
                    minimum_eigenvalue_this_trajectory
                    < -POSITIVITY_TOL
                ):
                    raise RuntimeError(
                        f"Density-matrix positivity failure for "
                        f"gamma={gamma_phi}, frame={frame}, "
                        f"initial={initial_state}: "
                        f"{minimum_eigenvalue_this_trajectory:.3e}"
                    )

                survival = clean_populations[
                    :,
                    initial_index,
                ]

                minimum_survival_index = int(
                    np.argmin(survival)
                )

                frame_row: dict[str, object] = {
                    "gamma_phi_ps_inv": float(gamma_phi),
                    "dephasing_time_ps": (
                        np.inf
                        if gamma_phi == 0.0
                        else 1.0 / float(gamma_phi)
                    ),
                    "frame": frame,
                    "snapshot_time_ps": snapshot_times_ps[frame],
                    "initial_state": initial_state,
                    "minimum_survival_probability": float(
                        survival[minimum_survival_index]
                    ),
                    "time_of_minimum_survival_ps": float(
                        times_ps[minimum_survival_index]
                    ),
                    "time_average_survival_probability": float(
                        np.mean(survival)
                    ),
                    "final_survival_probability": float(
                        survival[-1]
                    ),
                    "time_integrated_l1_coherence": float(
                        np.trapezoid(
                            coherence,
                            times_ps,
                        )
                    ),
                    "maximum_trace_error": trace_error,
                    "maximum_sampled_hermiticity_error": (
                        maximum_hermiticity_this_trajectory
                    ),
                    "minimum_sampled_density_eigenvalue": (
                        minimum_eigenvalue_this_trajectory
                    ),
                }

                for target_index, target_state in enumerate(SITES):
                    target_population = clean_populations[
                        :,
                        target_index,
                    ]

                    maximum_index = int(
                        np.argmax(target_population)
                    )

                    frame_row[
                        f"maximum_population_{target_state}"
                    ] = float(
                        target_population[maximum_index]
                    )
                    frame_row[
                        f"time_of_maximum_{target_state}_ps"
                    ] = float(
                        times_ps[maximum_index]
                    )
                    frame_row[
                        f"time_average_{target_state}"
                    ] = float(
                        np.mean(target_population)
                    )
                    frame_row[
                        f"final_population_{target_state}"
                    ] = float(
                        target_population[-1]
                    )

                frame_metric_rows.append(frame_row)

        log(
            f"[gamma={gamma_phi:g} ps^-1] "
            f"completed {EXPECTED_FRAMES}/{EXPECTED_FRAMES} frames"
        )

    gamma0_index = int(
        np.where(
            np.isclose(
                GAMMA_PHI_PS_INV,
                0.0,
                atol=0.0,
                rtol=0.0,
            )
        )[0][0]
    )

    if not COHERENT_NPZ.is_file():
        raise RuntimeError(
            f"Missing coherent reference archive: {COHERENT_NPZ}"
        )

    coherent_reference = np.load(
        COHERENT_NPZ,
        allow_pickle=False,
    )

    previous_times = coherent_reference["times_ps"]

    previous_time_indices = np.rint(
        times_ps / float(previous_times[1] - previous_times[0])
    ).astype(np.int64)

    if np.max(previous_time_indices) >= previous_times.size:
        raise RuntimeError(
            "Current time grid extends beyond the coherent reference."
        )

    if not np.allclose(
        previous_times[previous_time_indices],
        times_ps,
        atol=1.0e-12,
        rtol=0.0,
    ):
        raise RuntimeError(
            "The current time grid cannot be aligned exactly with "
            "the coherent reference."
        )

    gamma0_maximum_population_error = 0.0

    for initial_index, initial_state in enumerate(SITES):
        key = (
            "tdcac_corrected__initial_"
            f"{initial_state}"
        )

        if key not in coherent_reference:
            raise RuntimeError(
                f"Missing coherent reference key: {key}"
            )

        reference_population = coherent_reference[key][
            :,
            previous_time_indices,
            :,
        ]

        current_population = populations[
            gamma0_index,
            initial_index,
        ].astype(np.float64)

        error = float(
            np.max(
                np.abs(
                    current_population
                    - reference_population
                )
            )
        )

        gamma0_maximum_population_error = max(
            gamma0_maximum_population_error,
            error,
        )

    if gamma0_maximum_population_error > GAMMA0_MATCH_TOL:
        raise RuntimeError(
            "The gamma=0 Liouvillian propagation does not reproduce "
            "the accepted coherent trajectories. "
            f"Maximum error={gamma0_maximum_population_error:.3e}"
        )

    ensemble_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    gamma_aggregate_rows: list[dict[str, object]] = []

    selected_coherence_times = (
        0.5,
        1.0,
        2.0,
        5.0,
        10.0,
        20.0,
    )

    for gamma_index, gamma_phi in enumerate(GAMMA_PHI_PS_INV):
        for initial_index, initial_state in enumerate(SITES):
            population_stack = populations[
                gamma_index,
                initial_index,
            ].astype(np.float64)

            coherence_stack = coherence_l1[
                gamma_index,
                initial_index,
            ].astype(np.float64)

            mean_population = np.mean(
                population_stack,
                axis=0,
            )
            sd_population = np.std(
                population_stack,
                axis=0,
                ddof=0,
            )
            q05_population = np.quantile(
                population_stack,
                0.05,
                axis=0,
            )
            q50_population = np.quantile(
                population_stack,
                0.50,
                axis=0,
            )
            q95_population = np.quantile(
                population_stack,
                0.95,
                axis=0,
            )

            mean_coherence = np.mean(
                coherence_stack,
                axis=0,
            )
            sd_coherence = np.std(
                coherence_stack,
                axis=0,
                ddof=0,
            )

            for time_index, time_ps in enumerate(times_ps):
                row: dict[str, object] = {
                    "gamma_phi_ps_inv": float(gamma_phi),
                    "dephasing_time_ps": (
                        np.inf
                        if gamma_phi == 0.0
                        else 1.0 / float(gamma_phi)
                    ),
                    "initial_state": initial_state,
                    "time_ps": float(time_ps),
                    "mean_l1_coherence": float(
                        mean_coherence[time_index]
                    ),
                    "sd_l1_coherence": float(
                        sd_coherence[time_index]
                    ),
                }

                for target_index, target_state in enumerate(SITES):
                    row[
                        f"mean_{target_state}"
                    ] = float(
                        mean_population[
                            time_index,
                            target_index,
                        ]
                    )
                    row[
                        f"sd_{target_state}"
                    ] = float(
                        sd_population[
                            time_index,
                            target_index,
                        ]
                    )
                    row[
                        f"q05_{target_state}"
                    ] = float(
                        q05_population[
                            time_index,
                            target_index,
                        ]
                    )
                    row[
                        f"q50_{target_state}"
                    ] = float(
                        q50_population[
                            time_index,
                            target_index,
                        ]
                    )
                    row[
                        f"q95_{target_state}"
                    ] = float(
                        q95_population[
                            time_index,
                            target_index,
                        ]
                    )

                ensemble_rows.append(row)

            survival = mean_population[
                :,
                initial_index,
            ]

            minimum_survival_index = int(
                np.argmin(survival)
            )

            if initial_index in HIGH_ENERGY_INDICES:
                other_high_indices = [
                    index
                    for index in HIGH_ENERGY_INDICES
                    if index != initial_index
                ]

                other_high_population = np.sum(
                    mean_population[
                        :,
                        other_high_indices,
                    ],
                    axis=1,
                )
            else:
                other_high_population = np.sum(
                    mean_population[
                        :,
                        HIGH_ENERGY_INDICES,
                    ],
                    axis=1,
                )

            maximum_other_high_index = int(
                np.argmax(other_high_population)
            )

            pyr5_population = mean_population[
                :,
                PYR5_INDEX,
            ]
            maximum_pyr5_index = int(
                np.argmax(pyr5_population)
            )

            summary_row: dict[str, object] = {
                "gamma_phi_ps_inv": float(gamma_phi),
                "dephasing_time_ps": (
                    np.inf
                    if gamma_phi == 0.0
                    else 1.0 / float(gamma_phi)
                ),
                "initial_state": initial_state,
                "ensemble_minimum_survival_probability": float(
                    survival[minimum_survival_index]
                ),
                "time_of_ensemble_minimum_survival_ps": float(
                    times_ps[minimum_survival_index]
                ),
                "time_average_ensemble_survival": float(
                    np.mean(survival)
                ),
                "final_ensemble_survival": float(
                    survival[-1]
                ),
                "maximum_ensemble_transfer_out": float(
                    np.max(1.0 - survival)
                ),
                "maximum_ensemble_other_high_population": float(
                    other_high_population[
                        maximum_other_high_index
                    ]
                ),
                "time_of_maximum_other_high_population_ps": float(
                    times_ps[
                        maximum_other_high_index
                    ]
                ),
                "maximum_ensemble_PYR5_population": float(
                    pyr5_population[
                        maximum_pyr5_index
                    ]
                ),
                "time_of_maximum_PYR5_population_ps": float(
                    times_ps[maximum_pyr5_index]
                ),
                "time_average_l1_coherence": float(
                    np.mean(mean_coherence)
                ),
                "time_integrated_l1_coherence": float(
                    np.trapezoid(
                        mean_coherence,
                        times_ps,
                    )
                ),
                "mean_frame_minimum_survival": float(
                    np.mean(
                        np.min(
                            population_stack[
                                :,
                                :,
                                initial_index,
                            ],
                            axis=1,
                        )
                    )
                ),
                "maximum_single_frame_transfer_out": float(
                    np.max(
                        1.0
                        - population_stack[
                            :,
                            :,
                            initial_index,
                        ]
                    )
                ),
            }

            for requested_time in selected_coherence_times:
                time_index = nearest_time_index(
                    times_ps,
                    requested_time,
                )

                summary_row[
                    "mean_l1_coherence_at_"
                    f"{time_slug(requested_time)}ps"
                ] = float(
                    mean_coherence[time_index]
                )

            for target_index, target_state in enumerate(SITES):
                target_population = mean_population[
                    :,
                    target_index,
                ]
                maximum_index = int(
                    np.argmax(target_population)
                )

                summary_row[
                    f"maximum_ensemble_{target_state}"
                ] = float(
                    target_population[maximum_index]
                )
                summary_row[
                    f"time_of_maximum_{target_state}_ps"
                ] = float(
                    times_ps[maximum_index]
                )
                summary_row[
                    f"time_average_{target_state}"
                ] = float(
                    np.mean(target_population)
                )
                summary_row[
                    f"final_{target_state}"
                ] = float(
                    target_population[-1]
                )

            summary_rows.append(summary_row)

        high_rows = [
            row
            for row in summary_rows
            if float(row["gamma_phi_ps_inv"]) == float(gamma_phi)
            and row["initial_state"] in SITES[:3]
        ]

        gamma_aggregate_rows.append(
            {
                "gamma_phi_ps_inv": float(gamma_phi),
                "dephasing_time_ps": (
                    np.inf
                    if gamma_phi == 0.0
                    else 1.0 / float(gamma_phi)
                ),
                "mean_minimum_survival_PYR2_to_PYR4": float(
                    np.mean(
                        [
                            float(
                                row[
                                    "ensemble_minimum_survival_probability"
                                ]
                            )
                            for row in high_rows
                        ]
                    )
                ),
                "mean_time_average_survival_PYR2_to_PYR4": float(
                    np.mean(
                        [
                            float(
                                row[
                                    "time_average_ensemble_survival"
                                ]
                            )
                            for row in high_rows
                        ]
                    )
                ),
                "mean_maximum_other_high_population": float(
                    np.mean(
                        [
                            float(
                                row[
                                    "maximum_ensemble_other_high_population"
                                ]
                            )
                            for row in high_rows
                        ]
                    )
                ),
                "maximum_PYR5_population_from_high_states": float(
                    np.max(
                        [
                            float(
                                row[
                                    "maximum_ensemble_PYR5_population"
                                ]
                            )
                            for row in high_rows
                        ]
                    )
                ),
                "mean_time_integrated_l1_coherence": float(
                    np.mean(
                        [
                            float(
                                row[
                                    "time_integrated_l1_coherence"
                                ]
                            )
                            for row in high_rows
                        ]
                    )
                ),
                "mean_l1_coherence_at_20ps": float(
                    np.mean(
                        [
                            float(
                                row[
                                    "mean_l1_coherence_at_20ps"
                                ]
                            )
                            for row in high_rows
                        ]
                    )
                ),
            }
        )

    crossover_rows: list[dict[str, object]] = []

    for initial_state in SITES:
        rows = [
            row
            for row in summary_rows
            if row["initial_state"] == initial_state
        ]

        gamma0_row = [
            row
            for row in rows
            if float(row["gamma_phi_ps_inv"]) == 0.0
        ][0]

        maximum_transfer_row = max(
            rows,
            key=lambda row: float(
                row["maximum_ensemble_transfer_out"]
            ),
        )

        minimum_time_average_survival_row = min(
            rows,
            key=lambda row: float(
                row["time_average_ensemble_survival"]
            ),
        )

        maximum_other_high_row = max(
            rows,
            key=lambda row: float(
                row["maximum_ensemble_other_high_population"]
            ),
        )

        maximum_pyr5_row = max(
            rows,
            key=lambda row: float(
                row["maximum_ensemble_PYR5_population"]
            ),
        )

        crossover_rows.append(
            {
                "initial_state": initial_state,
                "gamma0_maximum_transfer_out": float(
                    gamma0_row["maximum_ensemble_transfer_out"]
                ),
                "gamma_maximizing_transfer_out_ps_inv": float(
                    maximum_transfer_row["gamma_phi_ps_inv"]
                ),
                "maximum_transfer_out": float(
                    maximum_transfer_row[
                        "maximum_ensemble_transfer_out"
                    ]
                ),
                "transfer_out_change_vs_gamma0": float(
                    maximum_transfer_row[
                        "maximum_ensemble_transfer_out"
                    ]
                )
                - float(
                    gamma0_row[
                        "maximum_ensemble_transfer_out"
                    ]
                ),
                "gamma_minimizing_time_average_survival_ps_inv": float(
                    minimum_time_average_survival_row[
                        "gamma_phi_ps_inv"
                    ]
                ),
                "minimum_time_average_survival": float(
                    minimum_time_average_survival_row[
                        "time_average_ensemble_survival"
                    ]
                ),
                "gamma_maximizing_other_high_population_ps_inv": float(
                    maximum_other_high_row["gamma_phi_ps_inv"]
                ),
                "maximum_other_high_population": float(
                    maximum_other_high_row[
                        "maximum_ensemble_other_high_population"
                    ]
                ),
                "gamma_maximizing_PYR5_population_ps_inv": float(
                    maximum_pyr5_row["gamma_phi_ps_inv"]
                ),
                "maximum_PYR5_population": float(
                    maximum_pyr5_row[
                        "maximum_ensemble_PYR5_population"
                    ]
                ),
            }
        )

    validation_rows = [
        {
            "scipy_version": scipy.__version__,
            "n_gamma_values": int(n_gamma),
            "n_frames": EXPECTED_FRAMES,
            "n_initial_states": N_STATES,
            "n_time_points": int(n_times),
            "maximum_trace_error": maximum_trace_error,
            "maximum_hermiticity_error": maximum_hermiticity_error,
            "minimum_sampled_density_eigenvalue": (
                minimum_sampled_density_eigenvalue
            ),
            "minimum_raw_population": minimum_raw_population,
            "maximum_raw_population": maximum_raw_population,
            "gamma0_maximum_population_error_vs_coherent_reference": (
                gamma0_maximum_population_error
            ),
            "trace_validation_pass": (
                maximum_trace_error <= TRACE_TOL
            ),
            "hermiticity_validation_pass": (
                maximum_hermiticity_error
                <= HERMITICITY_TOL
            ),
            "positivity_validation_pass": (
                minimum_sampled_density_eigenvalue
                >= -POSITIVITY_TOL
                and minimum_raw_population
                >= -POSITIVITY_TOL
                and maximum_raw_population
                <= 1.0 + POSITIVITY_TOL
            ),
            "gamma0_reference_validation_pass": (
                gamma0_maximum_population_error
                <= GAMMA0_MATCH_TOL
            ),
        }
    ]

    write_csv(ENSEMBLE_CSV, ensemble_rows)
    write_csv(SUMMARY_CSV, summary_rows)
    write_csv(FRAME_METRICS_CSV, frame_metric_rows)
    write_csv(GAMMA_AGGREGATE_CSV, gamma_aggregate_rows)
    write_csv(CROSSOVER_CSV, crossover_rows)
    write_csv(VALIDATION_CSV, validation_rows)

    np.savez_compressed(
        NPZ_PATH,
        gamma_phi_ps_inv=GAMMA_PHI_PS_INV,
        times_ps=times_ps,
        snapshot_times_ps=snapshot_times_ps,
        site_labels=np.asarray(SITES),
        populations=populations,
        coherence_l1=coherence_l1,
    )

    gamma0_aggregate = [
        row
        for row in gamma_aggregate_rows
        if float(row["gamma_phi_ps_inv"]) == 0.0
    ][0]

    maximum_high_transfer_row = max(
        gamma_aggregate_rows,
        key=lambda row: float(
            row["mean_maximum_other_high_population"]
        ),
    )

    maximum_pyr5_aggregate_row = max(
        gamma_aggregate_rows,
        key=lambda row: float(
            row["maximum_PYR5_population_from_high_states"]
        ),
    )

    minimum_survival_aggregate_row = min(
        gamma_aggregate_rows,
        key=lambda row: float(
            row["mean_time_average_survival_PYR2_to_PYR4"]
        ),
    )

    all_validation_pass = all(
        bool(value)
        for key, value in validation_rows[0].items()
        if key.endswith("_pass")
    )

    with REPORT_MD.open("w", encoding="utf-8") as handle:
        handle.write(
            "# Day020 Haken-Strobl high-dephasing sensitivity\n\n"
        )

        handle.write("## Protocol\n\n")
        handle.write(
            "- Primary Hamiltonian: four-state bright TDC-AC-corrected model.\n"
        )
        handle.write(
            "- Twenty-one solvent snapshots were propagated independently.\n"
        )
        handle.write(
            "- Lindblad operators: "
            "`L_i = sqrt(gamma_phi) |i><i|` for each local bright state.\n"
        )
        handle.write(
            "- Therefore each off-diagonal density-matrix element decays "
            "directly at `gamma_phi` in the absence of the Hamiltonian.\n"
        )
        handle.write(
            f"- Gamma sweep: "
            f"{', '.join(f'{value:g}' for value in GAMMA_PHI_PS_INV)} ps^-1.\n"
        )
        handle.write(
            f"- Propagation interval: 0-{T_MAX_PS:.3f} ps; "
            f"output interval: {DT_PS:.3f} ps.\n"
        )
        handle.write(
            "- No population-relaxation operators or detailed-balance "
            "constraints were included.\n\n"
        )

        handle.write("## Numerical validation\n\n")
        handle.write(
            f"- Maximum trace error: {maximum_trace_error:.3e}.\n"
        )
        handle.write(
            f"- Maximum Hermiticity error: "
            f"{maximum_hermiticity_error:.3e}.\n"
        )
        handle.write(
            f"- Minimum sampled density-matrix eigenvalue: "
            f"{minimum_sampled_density_eigenvalue:.3e}.\n"
        )
        handle.write(
            f"- Gamma=0 maximum population error versus the accepted "
            f"coherent trajectories: "
            f"{gamma0_maximum_population_error:.3e}.\n"
        )
        handle.write(
            f"- Overall numerical validation: "
            f"{'PASS' if all_validation_pass else 'REVIEW'}.\n\n"
        )

        handle.write("## Aggregate high-energy-manifold response\n\n")
        handle.write(
            "| gamma_phi (ps^-1) | T_phi (ps) | "
            "Mean minimum survival | Mean time-averaged survival | "
            "Mean maximum population on another PYR2-PYR4 site | "
            "Maximum PYR5 population | Integrated l1 coherence |\n"
        )
        handle.write(
            "|---:|---:|---:|---:|---:|---:|---:|\n"
        )

        for row in gamma_aggregate_rows:
            handle.write(
                f"| {float(row['gamma_phi_ps_inv']):g} "
                f"| {float(row['dephasing_time_ps']):g} "
                f"| {float(row['mean_minimum_survival_PYR2_to_PYR4']):.6f} "
                f"| {float(row['mean_time_average_survival_PYR2_to_PYR4']):.6f} "
                f"| {float(row['mean_maximum_other_high_population']):.6f} "
                f"| {float(row['maximum_PYR5_population_from_high_states']):.6e} "
                f"| {float(row['mean_time_integrated_l1_coherence']):.6f} |\n"
            )

        handle.write("\n## Sensitivity extrema\n\n")
        handle.write(
            f"- Gamma=0 mean maximum high-manifold transfer: "
            f"{float(gamma0_aggregate['mean_maximum_other_high_population']):.6f}.\n"
        )
        handle.write(
            f"- Gamma maximizing mean high-manifold transfer: "
            f"{float(maximum_high_transfer_row['gamma_phi_ps_inv']):g} ps^-1, "
            f"with value "
            f"{float(maximum_high_transfer_row['mean_maximum_other_high_population']):.6f}.\n"
        )
        handle.write(
            f"- Gamma minimizing mean time-averaged survival: "
            f"{float(minimum_survival_aggregate_row['gamma_phi_ps_inv']):g} ps^-1, "
            f"with value "
            f"{float(minimum_survival_aggregate_row['mean_time_average_survival_PYR2_to_PYR4']):.6f}.\n"
        )
        handle.write(
            f"- Largest PYR5 population anywhere in the phenomenological "
            f"gamma sweep: "
            f"{float(maximum_pyr5_aggregate_row['maximum_PYR5_population_from_high_states']):.6e} "
            f"at gamma="
            f"{float(maximum_pyr5_aggregate_row['gamma_phi_ps_inv']):g} ps^-1.\n\n"
        )

        handle.write("## Interpretation boundary\n\n")
        handle.write(
            "This calculation is a phenomenological robustness analysis. "
            "The gamma values are not extracted from the 21 solvent snapshots "
            "and are not microscopic dephasing rates. Pure dephasing can "
            "convert coherent coupling into population redistribution, but "
            "it contains no energy-selective downhill relaxation and does not "
            "enforce thermal detailed balance. Any increase of PYR5 population "
            "must therefore not be interpreted as physical bath-assisted "
            "relaxation. A subsequent relaxation model requires explicitly "
            "declared rates or a separately justified bath spectral density.\n"
        )

    log("")
    log("Day020 Haken-Strobl high-dephasing sensitivity completed.")
    log(f"Gamma values: {n_gamma}/{n_gamma}")
    log(f"Frames per gamma: {EXPECTED_FRAMES}/{EXPECTED_FRAMES}")
    log(f"Initial states: {N_STATES}/{N_STATES}")
    log(f"Time points: {n_times}")
    log(f"Maximum trace error: {maximum_trace_error:.3e}")
    log(
        f"Minimum sampled density eigenvalue: "
        f"{minimum_sampled_density_eigenvalue:.3e}"
    )
    log(
        f"Gamma=0 coherent-reference error: "
        f"{gamma0_maximum_population_error:.3e}"
    )
    log(
        "Gamma maximizing mean high-manifold transfer: "
        f"{float(maximum_high_transfer_row['gamma_phi_ps_inv']):g} ps^-1"
    )
    log(
        "Maximum PYR5 population across sweep: "
        f"{float(maximum_pyr5_aggregate_row['maximum_PYR5_population_from_high_states']):.6e}"
    )
    log(
        f"Overall numerical validation: "
        f"{'PASS' if all_validation_pass else 'REVIEW'}"
    )
    log(f"Wrote: {OUTPUT_ROOT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
