#!/usr/bin/env python3

from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np

try:
    from scipy.sparse import csr_matrix
    from scipy.sparse.linalg import expm_multiply
except ImportError as exc:
    raise SystemExit(
        "SciPy is required. Activate the project environment."
    ) from exc


PROJECT_ROOT = Path(__file__).resolve().parents[2]

HAMILTONIAN_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_bright_finite_size_corrected_hamiltonians/"
    "hamiltonian_snapshots_bright4"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day020_bright_detailed_balance_relaxation"
)

FRAME_METRICS_CSV = OUTPUT_ROOT / "frame_relaxation_metrics.csv"
CONDITION_SUMMARY_CSV = OUTPUT_ROOT / "condition_summary.csv"
EQUILIBRIUM_CSV = OUTPUT_ROOT / "equilibrium_populations.csv"
VALIDATION_CSV = OUTPUT_ROOT / "numerical_validation.csv"
NPZ_PATH = OUTPUT_ROOT / "ensemble_relaxation_trajectories.npz"
REPORT_MD = OUTPUT_ROOT / "DETAILED_BALANCE_RELAXATION_DAY020.md"

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

TEMPERATURES_K = np.array(
    [150.0, 200.0, 250.0, 300.0],
    dtype=np.float64,
)

# Phenomenological system-bath relaxation amplitudes.
# These are sensitivity parameters, not microscopic rates.
KAPPA_REF_PS_INV = np.array(
    [0.1, 1.0, 10.0],
    dtype=np.float64,
)

T_MAX_PS = 100.0
DT_PS = 0.05

HBAR_EV_PS = 6.582119569e-4
KB_EV_K = 8.617333262145e-5

TRACE_TOL = 1.0e-10
HERMITICITY_TOL = 1.0e-10
POSITIVITY_TOL = 1.0e-10
DETAILED_BALANCE_TOL = 1.0e-10
GIBBS_STATIONARITY_TOL = 1.0e-10
TRACE_PRESERVATION_TOL = 1.0e-10

VALIDATION_STRIDE = 20

FRAME_RE = re.compile(r"frame=(\d+)")
TIME_RE = re.compile(r"time_ps=([0-9.+\-Ee]+)")


def log(message: str = "") -> None:
    print(message, flush=True)


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
) -> None:
    if not rows:
        raise RuntimeError(f"No rows available for {path}")

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0].keys()),
        )
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
        raise RuntimeError(f"Missing frame number in {path}")

    if time_ps is None:
        raise RuntimeError(f"Missing snapshot time in {path}")

    if basis != SITES:
        raise RuntimeError(
            f"Unexpected basis in {path}: {basis}"
        )

    matrix = np.loadtxt(path, comments="#", dtype=np.float64)

    if matrix.shape != (N_STATES, N_STATES):
        raise RuntimeError(
            f"Unexpected matrix shape in {path}: {matrix.shape}"
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

    return frame, time_ps, basis, matrix


def read_hamiltonians() -> tuple[np.ndarray, np.ndarray]:
    files = sorted(
        HAMILTONIAN_ROOT.glob(
            "H_bright4_tdcac_frame*.dat"
        )
    )

    if len(files) != EXPECTED_FRAMES:
        raise RuntimeError(
            f"Expected {EXPECTED_FRAMES} Hamiltonians, "
            f"found {len(files)}"
        )

    snapshot_times: list[float] = []
    matrices: list[np.ndarray] = []

    for expected_frame, path in enumerate(files):
        frame, time_ps, _, matrix = parse_snapshot(path)

        if frame != expected_frame:
            raise RuntimeError(
                f"Frame mismatch in {path}: "
                f"expected {expected_frame}, found {frame}"
            )

        snapshot_times.append(time_ps)
        matrices.append(matrix)

    return (
        np.asarray(snapshot_times, dtype=np.float64),
        np.stack(matrices),
    )


def vectorize(matrix: np.ndarray) -> np.ndarray:
    return matrix.reshape(N_STATES**2, order="F")


def matrix_from_vector(vector: np.ndarray) -> np.ndarray:
    return vector.reshape(
        (N_STATES, N_STATES),
        order="F",
    )


def build_initial_density_vectors() -> np.ndarray:
    vectors = np.zeros(
        (N_STATES**2, N_STATES),
        dtype=np.complex128,
    )

    for initial_index in range(N_STATES):
        rho = np.zeros(
            (N_STATES, N_STATES),
            dtype=np.complex128,
        )
        rho[initial_index, initial_index] = 1.0
        vectors[:, initial_index] = vectorize(rho)

    return vectors


def jump_dissipator(
    jump_operator: np.ndarray,
) -> np.ndarray:
    identity = np.eye(
        N_STATES,
        dtype=np.complex128,
    )

    product = (
        jump_operator.conj().T
        @ jump_operator
    )

    return (
        np.kron(
            jump_operator.conj(),
            jump_operator,
        )
        - 0.5 * np.kron(identity, product)
        - 0.5 * np.kron(product.T, identity)
    )


def bose_population(
    gap_eV: float,
    temperature_K: float,
) -> float:
    x = gap_eV / (KB_EV_K * temperature_K)

    if x > 700.0:
        return 0.0

    return float(1.0 / np.expm1(x))


def build_global_thermal_liouvillian(
    hamiltonian_eV: np.ndarray,
    temperature_K: float,
    kappa_ref_ps_inv: float,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    float,
]:
    identity = np.eye(
        N_STATES,
        dtype=np.complex128,
    )

    coherent = (
        -1j
        / HBAR_EV_PS
        * (
            np.kron(identity, hamiltonian_eV)
            - np.kron(
                hamiltonian_eV.T,
                identity,
            )
        )
    )

    eigenvalues, eigenvectors = np.linalg.eigh(
        hamiltonian_eV
    )

    shifted = eigenvalues - np.min(eigenvalues)
    gibbs_weights = np.exp(
        -shifted / (KB_EV_K * temperature_K)
    )
    gibbs_weights /= np.sum(gibbs_weights)

    gibbs_density = (
        eigenvectors
        @ np.diag(gibbs_weights)
        @ eigenvectors.conj().T
    )

    liouvillian = coherent.copy()

    rate_matrix = np.zeros(
        (N_STATES, N_STATES),
        dtype=np.float64,
    )

    overlap_weights = np.zeros(
        (N_STATES, N_STATES),
        dtype=np.float64,
    )

    maximum_detailed_balance_error = 0.0

    site_probabilities = np.abs(eigenvectors) ** 2

    for low in range(N_STATES):
        for high in range(low + 1, N_STATES):
            gap_eV = float(
                eigenvalues[high]
                - eigenvalues[low]
            )

            overlap_weight = float(
                np.sum(
                    site_probabilities[:, low]
                    * site_probabilities[:, high]
                )
            )

            overlap_weights[low, high] = (
                overlap_weight
            )
            overlap_weights[high, low] = (
                overlap_weight
            )

            n_bose = bose_population(
                gap_eV,
                temperature_K,
            )

            downward_rate = (
                kappa_ref_ps_inv
                * overlap_weight
                * (n_bose + 1.0)
            )

            upward_rate = (
                kappa_ref_ps_inv
                * overlap_weight
                * n_bose
            )

            rate_matrix[high, low] = downward_rate
            rate_matrix[low, high] = upward_rate

            ket_low = eigenvectors[:, low]
            ket_high = eigenvectors[:, high]

            if downward_rate > 0.0:
                jump_down = (
                    np.sqrt(downward_rate)
                    * np.outer(
                        ket_low,
                        ket_high.conj(),
                    )
                )

                liouvillian += jump_dissipator(
                    jump_down
                )

            if upward_rate > 0.0:
                jump_up = (
                    np.sqrt(upward_rate)
                    * np.outer(
                        ket_high,
                        ket_low.conj(),
                    )
                )

                liouvillian += jump_dissipator(
                    jump_up
                )

            if (
                upward_rate > 0.0
                and downward_rate > 0.0
            ):
                beta_gap = (
                    gap_eV
                    / (KB_EV_K * temperature_K)
                )

                detailed_balance_error = abs(
                    np.log(
                        upward_rate
                        / downward_rate
                    )
                    + beta_gap
                )

                maximum_detailed_balance_error = max(
                    maximum_detailed_balance_error,
                    detailed_balance_error,
                )

    trace_vector = np.zeros(
        N_STATES**2,
        dtype=np.complex128,
    )

    for index in range(N_STATES):
        trace_vector[
            index + index * N_STATES
        ] = 1.0

    trace_preservation_error = float(
        np.linalg.norm(
            trace_vector @ liouvillian
        )
    )

    gibbs_stationarity_error = float(
        np.linalg.norm(
            liouvillian
            @ vectorize(gibbs_density)
        )
    )

    return (
        liouvillian,
        gibbs_density,
        rate_matrix,
        overlap_weights,
        max(
            maximum_detailed_balance_error,
            trace_preservation_error,
            gibbs_stationarity_error,
        ),
    )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    snapshot_times_ps, hamiltonians = (
        read_hamiltonians()
    )

    times_ps = np.arange(
        0.0,
        T_MAX_PS + 0.5 * DT_PS,
        DT_PS,
        dtype=np.float64,
    )

    if abs(times_ps[-1] - T_MAX_PS) > 1.0e-12:
        raise RuntimeError(
            "Time grid does not end at T_MAX_PS"
        )

    initial_vectors = (
        build_initial_density_vectors()
    )

    n_temperature = TEMPERATURES_K.size
    n_kappa = KAPPA_REF_PS_INV.size
    n_times = times_ps.size

    ensemble_sum = np.zeros(
        (
            n_temperature,
            n_kappa,
            N_STATES,
            n_times,
            N_STATES,
        ),
        dtype=np.float64,
    )

    ensemble_sum_sq = np.zeros_like(
        ensemble_sum
    )

    equilibrium_sum = np.zeros(
        (
            n_temperature,
            N_STATES,
        ),
        dtype=np.float64,
    )

    frame_rows: list[dict[str, object]] = []
    equilibrium_rows: list[dict[str, object]] = []

    maximum_trace_error = 0.0
    maximum_hermiticity_error = 0.0
    minimum_density_eigenvalue = np.inf
    maximum_detailed_balance_error = 0.0
    maximum_gibbs_stationarity_error = 0.0
    maximum_trace_preservation_error = 0.0

    log(
        "Day020 detailed-balance relaxation sensitivity"
    )
    log(
        "Temperatures: "
        + ", ".join(
            f"{value:g}"
            for value in TEMPERATURES_K
        )
        + " K"
    )
    log(
        "Kappa values: "
        + ", ".join(
            f"{value:g}"
            for value in KAPPA_REF_PS_INV
        )
        + " ps^-1"
    )
    log(
        f"Propagation: 0-{T_MAX_PS:g} ps, "
        f"dt={DT_PS:g} ps, "
        f"{n_times} points"
    )

    for temperature_index, temperature_K in enumerate(
        TEMPERATURES_K
    ):
        for frame_index, hamiltonian in enumerate(
            hamiltonians
        ):
            (
                unit_liouvillian,
                gibbs_density,
                _,
                _,
                unit_validation,
            ) = build_global_thermal_liouvillian(
                hamiltonian,
                float(temperature_K),
                1.0,
            )

            gibbs_site_population = np.real(
                np.diag(gibbs_density)
            )

            equilibrium_sum[
                temperature_index
            ] += gibbs_site_population

            equilibrium_rows.append(
                {
                    "temperature_K": float(
                        temperature_K
                    ),
                    "frame": frame_index,
                    "snapshot_time_ps": float(
                        snapshot_times_ps[
                            frame_index
                        ]
                    ),
                    "PYR2_equilibrium_population": float(
                        gibbs_site_population[0]
                    ),
                    "PYR3_equilibrium_population": float(
                        gibbs_site_population[1]
                    ),
                    "PYR4_equilibrium_population": float(
                        gibbs_site_population[2]
                    ),
                    "PYR5_equilibrium_population": float(
                        gibbs_site_population[3]
                    ),
                }
            )

            del unit_liouvillian
            del unit_validation

            for kappa_index, kappa_ref in enumerate(
                KAPPA_REF_PS_INV
            ):
                (
                    liouvillian,
                    gibbs_density,
                    rate_matrix,
                    overlap_weights,
                    _,
                ) = build_global_thermal_liouvillian(
                    hamiltonian,
                    float(temperature_K),
                    float(kappa_ref),
                )

                trace_vector = np.zeros(
                    N_STATES**2,
                    dtype=np.complex128,
                )

                for index in range(N_STATES):
                    trace_vector[
                        index + index * N_STATES
                    ] = 1.0

                trace_preservation_error = float(
                    np.linalg.norm(
                        trace_vector
                        @ liouvillian
                    )
                )

                gibbs_stationarity_error = float(
                    np.linalg.norm(
                        liouvillian
                        @ vectorize(gibbs_density)
                    )
                )

                maximum_trace_preservation_error = max(
                    maximum_trace_preservation_error,
                    trace_preservation_error,
                )

                maximum_gibbs_stationarity_error = max(
                    maximum_gibbs_stationarity_error,
                    gibbs_stationarity_error,
                )

                eigenvalues, eigenvectors = np.linalg.eigh(
                    hamiltonian
                )

                for low in range(N_STATES):
                    for high in range(
                        low + 1,
                        N_STATES,
                    ):
                        down = rate_matrix[
                            high,
                            low,
                        ]
                        up = rate_matrix[
                            low,
                            high,
                        ]

                        if up > 0.0 and down > 0.0:
                            gap = (
                                eigenvalues[high]
                                - eigenvalues[low]
                            )

                            error = abs(
                                np.log(up / down)
                                + gap
                                / (
                                    KB_EV_K
                                    * temperature_K
                                )
                            )

                            maximum_detailed_balance_error = max(
                                maximum_detailed_balance_error,
                                float(error),
                            )

                propagated = expm_multiply(
                    csr_matrix(liouvillian),
                    initial_vectors,
                    start=0.0,
                    stop=T_MAX_PS,
                    num=n_times,
                    endpoint=True,
                )

                populations = np.empty(
                    (
                        N_STATES,
                        n_times,
                        N_STATES,
                    ),
                    dtype=np.float64,
                )

                for initial_index in range(N_STATES):
                    for site_index in range(N_STATES):
                        vector_index = (
                            site_index
                            + site_index * N_STATES
                        )

                        populations[
                            initial_index,
                            :,
                            site_index,
                        ] = np.real(
                            propagated[
                                :,
                                vector_index,
                                initial_index,
                            ]
                        )

                ensemble_sum[
                    temperature_index,
                    kappa_index,
                ] += populations

                ensemble_sum_sq[
                    temperature_index,
                    kappa_index,
                ] += populations**2

                validation_indices = list(
                    range(
                        0,
                        n_times,
                        VALIDATION_STRIDE,
                    )
                )

                if validation_indices[-1] != n_times - 1:
                    validation_indices.append(
                        n_times - 1
                    )

                for initial_index in range(N_STATES):
                    for time_index in validation_indices:
                        rho = matrix_from_vector(
                            propagated[
                                time_index,
                                :,
                                initial_index,
                            ]
                        )

                        trace_error = abs(
                            np.trace(rho) - 1.0
                        )

                        hermiticity_error = float(
                            np.max(
                                np.abs(
                                    rho
                                    - rho.conj().T
                                )
                            )
                        )

                        minimum_eigenvalue = float(
                            np.min(
                                np.linalg.eigvalsh(
                                    0.5
                                    * (
                                        rho
                                        + rho.conj().T
                                    )
                                )
                            )
                        )

                        maximum_trace_error = max(
                            maximum_trace_error,
                            float(trace_error),
                        )

                        maximum_hermiticity_error = max(
                            maximum_hermiticity_error,
                            hermiticity_error,
                        )

                        minimum_density_eigenvalue = min(
                            minimum_density_eigenvalue,
                            minimum_eigenvalue,
                        )

                    trajectory = populations[
                        initial_index
                    ]

                    pyr5_trajectory = trajectory[
                        :,
                        PYR5_INDEX,
                    ]

                    maximum_pyr5_index = int(
                        np.argmax(pyr5_trajectory)
                    )

                    final_population = trajectory[-1]
                    gibbs_site_population = np.real(
                        np.diag(gibbs_density)
                    )

                    frame_rows.append(
                        {
                            "temperature_K": float(
                                temperature_K
                            ),
                            "kappa_ref_ps_inv": float(
                                kappa_ref
                            ),
                            "frame": frame_index,
                            "snapshot_time_ps": float(
                                snapshot_times_ps[
                                    frame_index
                                ]
                            ),
                            "initial_state": SITES[
                                initial_index
                            ],
                            "final_survival": float(
                                final_population[
                                    initial_index
                                ]
                            ),
                            "final_PYR5_population": float(
                                final_population[
                                    PYR5_INDEX
                                ]
                            ),
                            "maximum_PYR5_population": float(
                                pyr5_trajectory[
                                    maximum_pyr5_index
                                ]
                            ),
                            "time_of_maximum_PYR5_ps": float(
                                times_ps[
                                    maximum_pyr5_index
                                ]
                            ),
                            "final_l1_distance_to_Gibbs": float(
                                np.sum(
                                    np.abs(
                                        final_population
                                        - gibbs_site_population
                                    )
                                )
                            ),
                            "Gibbs_PYR5_population": float(
                                gibbs_site_population[
                                    PYR5_INDEX
                                ]
                            ),
                            "maximum_energy_basis_rate_ps_inv": float(
                                np.max(rate_matrix)
                            ),
                            "maximum_site_projector_overlap_weight": float(
                                np.max(
                                    overlap_weights
                                )
                            ),
                        }
                    )

        log(
            f"[T={temperature_K:g} K] "
            f"completed {EXPECTED_FRAMES}/"
            f"{EXPECTED_FRAMES} frames"
        )

    ensemble_mean = (
        ensemble_sum / EXPECTED_FRAMES
    )

    ensemble_variance = (
        ensemble_sum_sq / EXPECTED_FRAMES
        - ensemble_mean**2
    )

    ensemble_variance = np.maximum(
        ensemble_variance,
        0.0,
    )

    ensemble_sd = np.sqrt(
        ensemble_variance
    )

    equilibrium_mean = (
        equilibrium_sum / EXPECTED_FRAMES
    )

    condition_rows: list[dict[str, object]] = []

    for temperature_index, temperature_K in enumerate(
        TEMPERATURES_K
    ):
        for kappa_index, kappa_ref in enumerate(
            KAPPA_REF_PS_INV
        ):
            high_initial_pyr5_maxima: list[float] = []
            high_initial_pyr5_final: list[float] = []
            high_initial_final_distance: list[float] = []

            for initial_index in HIGH_ENERGY_INDICES:
                trajectory = ensemble_mean[
                    temperature_index,
                    kappa_index,
                    initial_index,
                ]

                high_initial_pyr5_maxima.append(
                    float(
                        np.max(
                            trajectory[
                                :,
                                PYR5_INDEX,
                            ]
                        )
                    )
                )

                high_initial_pyr5_final.append(
                    float(
                        trajectory[
                            -1,
                            PYR5_INDEX,
                        ]
                    )
                )

                high_initial_final_distance.append(
                    float(
                        np.sum(
                            np.abs(
                                trajectory[-1]
                                - equilibrium_mean[
                                    temperature_index
                                ]
                            )
                        )
                    )
                )

            condition_rows.append(
                {
                    "temperature_K": float(
                        temperature_K
                    ),
                    "kappa_ref_ps_inv": float(
                        kappa_ref
                    ),
                    "mean_maximum_PYR5_population_from_high_states": float(
                        np.mean(
                            high_initial_pyr5_maxima
                        )
                    ),
                    "mean_final_PYR5_population_from_high_states": float(
                        np.mean(
                            high_initial_pyr5_final
                        )
                    ),
                    "mean_final_l1_distance_to_Gibbs_from_high_states": float(
                        np.mean(
                            high_initial_final_distance
                        )
                    ),
                    "ensemble_Gibbs_PYR5_population": float(
                        equilibrium_mean[
                            temperature_index,
                            PYR5_INDEX,
                        ]
                    ),
                }
            )

    validation_rows = [
        {
            "n_temperatures": int(
                TEMPERATURES_K.size
            ),
            "n_kappa_values": int(
                KAPPA_REF_PS_INV.size
            ),
            "n_frames": EXPECTED_FRAMES,
            "n_initial_states": N_STATES,
            "n_time_points": n_times,
            "maximum_trace_error": maximum_trace_error,
            "maximum_hermiticity_error": (
                maximum_hermiticity_error
            ),
            "minimum_sampled_density_eigenvalue": (
                minimum_density_eigenvalue
            ),
            "maximum_detailed_balance_error": (
                maximum_detailed_balance_error
            ),
            "maximum_Gibbs_stationarity_error": (
                maximum_gibbs_stationarity_error
            ),
            "maximum_trace_preservation_error": (
                maximum_trace_preservation_error
            ),
            "trace_validation_pass": (
                maximum_trace_error <= TRACE_TOL
            ),
            "hermiticity_validation_pass": (
                maximum_hermiticity_error
                <= HERMITICITY_TOL
            ),
            "positivity_validation_pass": (
                minimum_density_eigenvalue
                >= -POSITIVITY_TOL
            ),
            "detailed_balance_validation_pass": (
                maximum_detailed_balance_error
                <= DETAILED_BALANCE_TOL
            ),
            "Gibbs_stationarity_validation_pass": (
                maximum_gibbs_stationarity_error
                <= GIBBS_STATIONARITY_TOL
            ),
            "trace_preservation_validation_pass": (
                maximum_trace_preservation_error
                <= TRACE_PRESERVATION_TOL
            ),
        }
    ]

    write_csv(
        FRAME_METRICS_CSV,
        frame_rows,
    )

    write_csv(
        CONDITION_SUMMARY_CSV,
        condition_rows,
    )

    write_csv(
        EQUILIBRIUM_CSV,
        equilibrium_rows,
    )

    write_csv(
        VALIDATION_CSV,
        validation_rows,
    )

    np.savez_compressed(
        NPZ_PATH,
        temperatures_K=TEMPERATURES_K,
        kappa_ref_ps_inv=KAPPA_REF_PS_INV,
        times_ps=times_ps,
        snapshot_times_ps=snapshot_times_ps,
        site_labels=np.asarray(SITES),
        ensemble_mean_populations=ensemble_mean,
        ensemble_sd_populations=ensemble_sd,
        ensemble_gibbs_site_populations=(
            equilibrium_mean
        ),
    )

    overall_pass = all(
        bool(value)
        for key, value in validation_rows[0].items()
        if key.endswith("_pass")
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Detailed-Balance Relaxation Sensitivity\n\n"
        )

        handle.write("## Scope\n\n")
        handle.write(
            "This calculation introduces population relaxation "
            "separately from the previous pure-dephasing model.\n\n"
        )

        handle.write(
            "The jump operators are constructed in the instantaneous "
            "energy eigenbasis of each bright-state Hamiltonian. "
            "Upward and downward rates satisfy detailed balance "
            "at the selected temperature.\n\n"
        )

        handle.write(
            "The relaxation amplitude `kappa_ref` is phenomenological "
            "and is not derived from a microscopic spectral density.\n\n"
        )

        handle.write("## Parameter sweep\n\n")
        handle.write(
            "- Temperatures: "
            + ", ".join(
                f"{value:g}"
                for value in TEMPERATURES_K
            )
            + " K.\n"
        )
        handle.write(
            "- Kappa values: "
            + ", ".join(
                f"{value:g}"
                for value in KAPPA_REF_PS_INV
            )
            + " ps^-1.\n"
        )
        handle.write(
            f"- Propagation interval: 0-{T_MAX_PS:g} ps.\n"
        )
        handle.write(
            f"- Output interval: {DT_PS:g} ps.\n\n"
        )

        handle.write("## Numerical validation\n\n")
        handle.write(
            f"- Maximum trace error: "
            f"{maximum_trace_error:.3e}.\n"
        )
        handle.write(
            f"- Maximum Hermiticity error: "
            f"{maximum_hermiticity_error:.3e}.\n"
        )
        handle.write(
            f"- Minimum sampled density eigenvalue: "
            f"{minimum_density_eigenvalue:.3e}.\n"
        )
        handle.write(
            f"- Maximum detailed-balance error: "
            f"{maximum_detailed_balance_error:.3e}.\n"
        )
        handle.write(
            f"- Maximum Gibbs-stationarity error: "
            f"{maximum_gibbs_stationarity_error:.3e}.\n"
        )
        handle.write(
            f"- Overall validation: "
            f"{'PASS' if overall_pass else 'FAIL'}.\n\n"
        )

        handle.write("## Interpretation limits\n\n")
        handle.write(
            "This model tests whether a thermally consistent "
            "population-relaxation bath can transfer excitation "
            "toward the low-energy PYR5 state. It does not provide "
            "a microscopic relaxation time because no bath spectral "
            "density has been derived from the current trajectory.\n"
        )

    log("")
    log(
        "Day020 detailed-balance relaxation sensitivity completed."
    )
    log(
        f"Temperatures: "
        f"{TEMPERATURES_K.size}/"
        f"{TEMPERATURES_K.size}"
    )
    log(
        f"Kappa values: "
        f"{KAPPA_REF_PS_INV.size}/"
        f"{KAPPA_REF_PS_INV.size}"
    )
    log(
        f"Frames per condition: "
        f"{EXPECTED_FRAMES}/"
        f"{EXPECTED_FRAMES}"
    )
    log(
        f"Maximum trace error: "
        f"{maximum_trace_error:.3e}"
    )
    log(
        f"Maximum Hermiticity error: "
        f"{maximum_hermiticity_error:.3e}"
    )
    log(
        f"Minimum sampled density eigenvalue: "
        f"{minimum_density_eigenvalue:.3e}"
    )
    log(
        f"Maximum detailed-balance error: "
        f"{maximum_detailed_balance_error:.3e}"
    )
    log(
        f"Maximum Gibbs-stationarity error: "
        f"{maximum_gibbs_stationarity_error:.3e}"
    )
    log(
        f"Overall numerical validation: "
        f"{'PASS' if overall_pass else 'FAIL'}"
    )
    log(
        f"Wrote: "
        f"{OUTPUT_ROOT.relative_to(PROJECT_ROOT)}"
    )


if __name__ == "__main__":
    main()
