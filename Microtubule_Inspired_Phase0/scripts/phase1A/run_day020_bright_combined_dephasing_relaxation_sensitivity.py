#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

try:
    from scipy.optimize import linear_sum_assignment
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

RELAXATION_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day020_bright_detailed_balance_relaxation"
)

RELAXATION_REFERENCE_NPZ = (
    RELAXATION_ROOT / "ensemble_relaxation_trajectories.npz"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day020_bright_combined_dephasing_relaxation"
)

FRAME_METRICS_CSV = (
    OUTPUT_ROOT / "frame_combined_metrics.csv"
)
CONDITION_SUMMARY_CSV = (
    OUTPUT_ROOT / "condition_summary.csv"
)
STEADY_STATE_CSV = (
    OUTPUT_ROOT / "steady_state_metrics.csv"
)
VALIDATION_CSV = (
    OUTPUT_ROOT / "numerical_validation.csv"
)
NPZ_PATH = (
    OUTPUT_ROOT / "combined_population_trajectories.npz"
)
REPORT_MD = (
    OUTPUT_ROOT
    / "COMBINED_DEPHASING_RELAXATION_DAY020.md"
)

EXPECTED_FRAMES = 21
N_STATES = 4

SITES = (
    "PYR2_bright",
    "PYR3_bright",
    "PYR4_bright",
    "PYR5_bright",
)

HIGH_SITE_INDICES = (0, 1, 2)
PYR5_SITE_INDEX = 3

TEMPERATURES_K = np.array(
    [150.0, 300.0],
    dtype=np.float64,
)

KAPPA_REF_PS_INV = np.array(
    [0.1, 1.0, 10.0],
    dtype=np.float64,
)

GAMMA_PHI_PS_INV = np.array(
    [0.0, 1.0, 20.0, 100.0],
    dtype=np.float64,
)

T_MAX_PS = 100.0
DT_PS = 0.05

HBAR_EV_PS = 6.582119569e-4
KB_EV_K = 8.617333262145e-5

TRACE_TOL = 1.0e-10
HERMITICITY_TOL = 1.0e-10
POSITIVITY_TOL = 1.0e-10
TRACE_PRESERVATION_TOL = 1.0e-10
DETAILED_BALANCE_TOL = 1.0e-10
RELAXATION_GIBBS_TOL = 1.0e-10
STEADY_STATE_RESIDUAL_TOL = 1.0e-9
GAMMA0_REFERENCE_TOL = 1.0e-9

VALIDATION_STRIDE = 20

FRAME_RE = re.compile(r"frame=(\d+)")
TIME_RE = re.compile(
    r"time_ps=([0-9.+\-Ee]+)"
)


def log(message: str = "") -> None:
    print(message, flush=True)


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
) -> None:
    if not rows:
        raise RuntimeError(
            f"No rows available for {path}"
        )

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(rows)


def vectorize(
    matrix: np.ndarray,
) -> np.ndarray:
    return matrix.reshape(
        N_STATES**2,
        order="F",
    )


def matrix_from_vector(
    vector: np.ndarray,
) -> np.ndarray:
    return vector.reshape(
        (N_STATES, N_STATES),
        order="F",
    )


def parse_snapshot(
    path: Path,
) -> tuple[
    int,
    float,
    tuple[str, ...],
    np.ndarray,
]:
    lines = path.read_text(
        encoding="utf-8"
    ).splitlines()

    frame: int | None = None
    time_ps: float | None = None
    basis: tuple[str, ...] | None = None

    for line in lines:
        if not line.startswith("#"):
            continue

        frame_match = FRAME_RE.search(line)

        if frame_match is not None:
            frame = int(
                frame_match.group(1)
            )

        time_match = TIME_RE.search(line)

        if time_match is not None:
            time_ps = float(
                time_match.group(1)
            )

        if line.startswith("# basis:"):
            basis = tuple(
                line.split(":", 1)[1].split()
            )

    if frame is None or time_ps is None:
        raise RuntimeError(
            f"Missing frame metadata in {path}"
        )

    if basis != SITES:
        raise RuntimeError(
            f"Unexpected basis in {path}: {basis}"
        )

    matrix = np.loadtxt(
        path,
        comments="#",
        dtype=np.float64,
    )

    if matrix.shape != (
        N_STATES,
        N_STATES,
    ):
        raise RuntimeError(
            f"Unexpected matrix shape in "
            f"{path}: {matrix.shape}"
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


def read_hamiltonians(
) -> tuple[np.ndarray, np.ndarray]:
    files = sorted(
        HAMILTONIAN_ROOT.glob(
            "H_bright4_tdcac_frame*.dat"
        )
    )

    if len(files) != EXPECTED_FRAMES:
        raise RuntimeError(
            f"Expected {EXPECTED_FRAMES} "
            f"Hamiltonians, found {len(files)}"
        )

    snapshot_times: list[float] = []
    matrices: list[np.ndarray] = []

    for expected_frame, path in enumerate(
        files
    ):
        frame, time_ps, _, matrix = (
            parse_snapshot(path)
        )

        if frame != expected_frame:
            raise RuntimeError(
                f"Frame mismatch in {path}: "
                f"expected {expected_frame}, "
                f"found {frame}"
            )

        snapshot_times.append(time_ps)
        matrices.append(matrix)

    return (
        np.asarray(
            snapshot_times,
            dtype=np.float64,
        ),
        np.stack(matrices),
    )


def build_initial_density_vectors(
) -> np.ndarray:
    vectors = np.zeros(
        (N_STATES**2, N_STATES),
        dtype=np.complex128,
    )

    for initial_index in range(
        N_STATES
    ):
        rho = np.zeros(
            (N_STATES, N_STATES),
            dtype=np.complex128,
        )
        rho[
            initial_index,
            initial_index,
        ] = 1.0

        vectors[:, initial_index] = (
            vectorize(rho)
        )

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
        - 0.5
        * np.kron(
            identity,
            product,
        )
        - 0.5
        * np.kron(
            product.T,
            identity,
        )
    )


def build_pure_dephasing_superoperator(
) -> np.ndarray:
    superoperator = np.zeros(
        (
            N_STATES**2,
            N_STATES**2,
        ),
        dtype=np.complex128,
    )

    for site_index in range(
        N_STATES
    ):
        projector = np.zeros(
            (N_STATES, N_STATES),
            dtype=np.complex128,
        )

        projector[
            site_index,
            site_index,
        ] = 1.0

        superoperator += (
            jump_dissipator(projector)
        )

    test_matrix = np.zeros(
        (N_STATES, N_STATES),
        dtype=np.complex128,
    )

    test_matrix[0, 1] = 1.0

    test_vector = vectorize(
        test_matrix
    )

    if not np.allclose(
        superoperator @ test_vector,
        -test_vector,
        atol=1.0e-12,
        rtol=0.0,
    ):
        raise RuntimeError(
            "Pure-dephasing superoperator "
            "failed the analytical test"
        )

    return superoperator


def bose_population(
    gap_eV: float,
    temperature_K: float,
) -> float:
    exponent = gap_eV / (
        KB_EV_K * temperature_K
    )

    if exponent > 700.0:
        return 0.0

    return float(
        1.0 / np.expm1(exponent)
    )


def coherent_superoperator(
    hamiltonian_eV: np.ndarray,
) -> np.ndarray:
    identity = np.eye(
        N_STATES,
        dtype=np.complex128,
    )

    return (
        -1j
        / HBAR_EV_PS
        * (
            np.kron(
                identity,
                hamiltonian_eV,
            )
            - np.kron(
                hamiltonian_eV.T,
                identity,
            )
        )
    )


def build_thermal_relaxation_superoperator(
    hamiltonian_eV: np.ndarray,
    temperature_K: float,
    kappa_ref_ps_inv: float,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    dict[int, int],
    int,
    float,
]:
    eigenvalues, eigenvectors = (
        np.linalg.eigh(
            hamiltonian_eV
        )
    )

    site_probabilities = (
        np.abs(eigenvectors) ** 2
    )

    (
        site_indices,
        eigenstate_indices,
    ) = linear_sum_assignment(
        -site_probabilities
    )

    site_to_eigenstate = {
        int(site): int(eigenstate)
        for site, eigenstate in zip(
            site_indices,
            eigenstate_indices,
        )
    }

    eigenstate_to_site = {
        eigenstate: site
        for site, eigenstate
        in site_to_eigenstate.items()
    }

    sink_eigenstate = (
        site_to_eigenstate[
            PYR5_SITE_INDEX
        ]
    )

    if sink_eigenstate != 0:
        raise RuntimeError(
            "The PYR5-like state is not the "
            "lowest-energy eigenstate"
        )

    shifted = (
        eigenvalues
        - np.min(eigenvalues)
    )

    gibbs_weights = np.exp(
        -shifted
        / (
            KB_EV_K
            * temperature_K
        )
    )

    gibbs_weights /= np.sum(
        gibbs_weights
    )

    gibbs_density = (
        eigenvectors
        @ np.diag(gibbs_weights)
        @ eigenvectors.conj().T
    )

    relaxation = np.zeros(
        (
            N_STATES**2,
            N_STATES**2,
        ),
        dtype=np.complex128,
    )

    maximum_detailed_balance_error = 0.0

    for low in range(N_STATES):
        for high in range(
            low + 1,
            N_STATES,
        ):
            gap_eV = float(
                eigenvalues[high]
                - eigenvalues[low]
            )

            overlap_weight = float(
                np.sum(
                    site_probabilities[
                        :,
                        low,
                    ]
                    * site_probabilities[
                        :,
                        high,
                    ]
                )
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

            ket_low = eigenvectors[:, low]
            ket_high = eigenvectors[:, high]

            if downward_rate > 0.0:
                jump_down = (
                    np.sqrt(
                        downward_rate
                    )
                    * np.outer(
                        ket_low,
                        ket_high.conj(),
                    )
                )

                relaxation += (
                    jump_dissipator(
                        jump_down
                    )
                )

            if upward_rate > 0.0:
                jump_up = (
                    np.sqrt(
                        upward_rate
                    )
                    * np.outer(
                        ket_high,
                        ket_low.conj(),
                    )
                )

                relaxation += (
                    jump_dissipator(
                        jump_up
                    )
                )

            if (
                upward_rate > 0.0
                and downward_rate > 0.0
            ):
                error = abs(
                    np.log(
                        upward_rate
                        / downward_rate
                    )
                    + gap_eV
                    / (
                        KB_EV_K
                        * temperature_K
                    )
                )

                maximum_detailed_balance_error = max(
                    maximum_detailed_balance_error,
                    float(error),
                )

    return (
        relaxation,
        gibbs_density,
        eigenvalues,
        eigenvectors,
        eigenstate_to_site,
        sink_eigenstate,
        maximum_detailed_balance_error,
    )


def trace_vector() -> np.ndarray:
    vector = np.zeros(
        N_STATES**2,
        dtype=np.complex128,
    )

    for index in range(N_STATES):
        vector[
            index
            + index * N_STATES
        ] = 1.0

    return vector


def solve_stationary_density(
    liouvillian: np.ndarray,
) -> tuple[np.ndarray, float]:
    trace_row = trace_vector()

    augmented = liouvillian.copy()

    right_hand_side = np.zeros(
        N_STATES**2,
        dtype=np.complex128,
    )

    augmented[-1, :] = trace_row
    right_hand_side[-1] = 1.0

    stationary_vector, *_ = (
        np.linalg.lstsq(
            augmented,
            right_hand_side,
            rcond=None,
        )
    )

    rho = matrix_from_vector(
        stationary_vector
    )

    rho = 0.5 * (
        rho + rho.conj().T
    )

    rho /= np.trace(rho)

    residual = float(
        np.linalg.norm(
            liouvillian
            @ vectorize(rho)
        )
    )

    return rho, residual


def load_relaxation_reference(
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    with np.load(
        RELAXATION_REFERENCE_NPZ
    ) as data:
        temperatures = data[
            "temperatures_K"
        ].astype(np.float64)

        kappas = data[
            "kappa_ref_ps_inv"
        ].astype(np.float64)

        times = data[
            "times_ps"
        ].astype(np.float64)

        populations = data[
            "ensemble_mean_populations"
        ].astype(np.float64)

    return (
        temperatures,
        kappas,
        times,
        populations,
    )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    (
        snapshot_times_ps,
        hamiltonians,
    ) = read_hamiltonians()

    initial_vectors = (
        build_initial_density_vectors()
    )

    dephasing_superoperator = (
        build_pure_dephasing_superoperator()
    )

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
            "Time grid does not end "
            "at T_MAX_PS"
        )

    (
        reference_temperatures,
        reference_kappas,
        reference_times,
        reference_populations,
    ) = load_relaxation_reference()

    if not np.allclose(
        reference_times,
        times_ps,
        atol=1.0e-12,
        rtol=0.0,
    ):
        raise RuntimeError(
            "Reference and combined "
            "time grids differ"
        )

    n_temperature = (
        TEMPERATURES_K.size
    )
    n_kappa = (
        KAPPA_REF_PS_INV.size
    )
    n_gamma = (
        GAMMA_PHI_PS_INV.size
    )
    n_times = times_ps.size

    ensemble_sum = np.zeros(
        (
            n_temperature,
            n_kappa,
            n_gamma,
            N_STATES,
            n_times,
            N_STATES,
        ),
        dtype=np.float64,
    )

    ensemble_sum_sq = np.zeros_like(
        ensemble_sum
    )

    coherence_sum = np.zeros(
        (
            n_temperature,
            n_kappa,
            n_gamma,
            N_STATES,
            n_times,
        ),
        dtype=np.float64,
    )

    frame_rows: list[
        dict[str, object]
    ] = []

    steady_rows: list[
        dict[str, object]
    ] = []

    maximum_trace_error = 0.0
    maximum_hermiticity_error = 0.0
    minimum_density_eigenvalue = np.inf

    maximum_trace_preservation_error = 0.0
    maximum_detailed_balance_error = 0.0
    maximum_relaxation_gibbs_error = 0.0
    maximum_stationary_residual = 0.0
    minimum_stationary_eigenvalue = np.inf
    maximum_gamma0_reference_error = 0.0

    trace_row = trace_vector()

    log(
        "Day020 combined dephasing "
        "+ detailed-balance relaxation"
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
            for value
            in KAPPA_REF_PS_INV
        )
        + " ps^-1"
    )

    log(
        "Gamma values: "
        + ", ".join(
            f"{value:g}"
            for value
            in GAMMA_PHI_PS_INV
        )
        + " ps^-1"
    )

    log(
        f"Propagation: 0-{T_MAX_PS:g} ps, "
        f"dt={DT_PS:g} ps, "
        f"{n_times} points"
    )

    for (
        temperature_index,
        temperature_K,
    ) in enumerate(TEMPERATURES_K):

        for (
            kappa_index,
            kappa_ref,
        ) in enumerate(
            KAPPA_REF_PS_INV
        ):

            for (
                gamma_index,
                gamma_phi,
            ) in enumerate(
                GAMMA_PHI_PS_INV
            ):

                for (
                    frame_index,
                    hamiltonian,
                ) in enumerate(
                    hamiltonians
                ):

                    (
                        relaxation_superoperator,
                        gibbs_density,
                        eigenvalues,
                        eigenvectors,
                        eigenstate_to_site,
                        sink_eigenstate,
                        detailed_balance_error,
                    ) = (
                        build_thermal_relaxation_superoperator(
                            hamiltonian,
                            float(
                                temperature_K
                            ),
                            float(
                                kappa_ref
                            ),
                        )
                    )

                    maximum_detailed_balance_error = max(
                        maximum_detailed_balance_error,
                        detailed_balance_error,
                    )

                    relaxation_gibbs_error = float(
                        np.linalg.norm(
                            relaxation_superoperator
                            @ vectorize(
                                gibbs_density
                            )
                        )
                    )

                    maximum_relaxation_gibbs_error = max(
                        maximum_relaxation_gibbs_error,
                        relaxation_gibbs_error,
                    )

                    liouvillian = (
                        coherent_superoperator(
                            hamiltonian
                        )
                        + relaxation_superoperator
                        + float(
                            gamma_phi
                        )
                        * dephasing_superoperator
                    )

                    trace_preservation_error = float(
                        np.linalg.norm(
                            trace_row
                            @ liouvillian
                        )
                    )

                    maximum_trace_preservation_error = max(
                        maximum_trace_preservation_error,
                        trace_preservation_error,
                    )

                    # For gamma_phi=0, the thermal generator has the
                    # analytically known Gibbs stationary state. Using it
                    # directly avoids loss of positivity from an ill-conditioned
                    # least-squares nullspace solve when Gibbs is nearly pure.
                    if np.isclose(
                        gamma_phi,
                        0.0,
                        atol=1.0e-15,
                        rtol=0.0,
                    ):
                        stationary_density = (
                            gibbs_density.copy()
                        )

                        stationary_residual = float(
                            np.linalg.norm(
                                liouvillian
                                @ vectorize(
                                    stationary_density
                                )
                            )
                        )
                    else:
                        (
                            stationary_density,
                            stationary_residual,
                        ) = solve_stationary_density(
                            liouvillian
                        )

                    maximum_stationary_residual = max(
                        maximum_stationary_residual,
                        stationary_residual,
                    )

                    stationary_eigenvalue = float(
                        np.min(
                            np.linalg.eigvalsh(
                                stationary_density
                            )
                        )
                    )

                    minimum_stationary_eigenvalue = min(
                        minimum_stationary_eigenvalue,
                        stationary_eigenvalue,
                    )

                    stationary_site_population = np.real(
                        np.diag(
                            stationary_density
                        )
                    )

                    gibbs_site_population = np.real(
                        np.diag(
                            gibbs_density
                        )
                    )

                    steady_rows.append(
                        {
                            "temperature_K": float(
                                temperature_K
                            ),
                            "kappa_ref_ps_inv": float(
                                kappa_ref
                            ),
                            "gamma_phi_ps_inv": float(
                                gamma_phi
                            ),
                            "frame": frame_index,
                            "snapshot_time_ps": float(
                                snapshot_times_ps[
                                    frame_index
                                ]
                            ),
                            "steady_PYR2_population": float(
                                stationary_site_population[
                                    0
                                ]
                            ),
                            "steady_PYR3_population": float(
                                stationary_site_population[
                                    1
                                ]
                            ),
                            "steady_PYR4_population": float(
                                stationary_site_population[
                                    2
                                ]
                            ),
                            "steady_PYR5_population": float(
                                stationary_site_population[
                                    3
                                ]
                            ),
                            "Gibbs_PYR5_population": float(
                                gibbs_site_population[
                                    3
                                ]
                            ),
                            "steady_l1_distance_to_Gibbs": float(
                                np.sum(
                                    np.abs(
                                        stationary_site_population
                                        - gibbs_site_population
                                    )
                                )
                            ),
                            "stationary_residual": (
                                stationary_residual
                            ),
                            "minimum_stationary_density_eigenvalue": (
                                stationary_eigenvalue
                            ),
                        }
                    )

                    propagated = expm_multiply(
                        csr_matrix(
                            liouvillian
                        ),
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

                    coherence_l1 = np.empty(
                        (
                            N_STATES,
                            n_times,
                        ),
                        dtype=np.float64,
                    )

                    validation_indices = list(
                        range(
                            0,
                            n_times,
                            VALIDATION_STRIDE,
                        )
                    )

                    if (
                        validation_indices[-1]
                        != n_times - 1
                    ):
                        validation_indices.append(
                            n_times - 1
                        )

                    energy_projectors = [
                        np.outer(
                            eigenvectors[
                                :,
                                state,
                            ],
                            eigenvectors[
                                :,
                                state,
                            ].conj(),
                        )
                        for state
                        in range(N_STATES)
                    ]

                    high_eigenstates = [
                        state
                        for state
                        in range(N_STATES)
                        if state
                        != sink_eigenstate
                    ]

                    for initial_index in range(
                        N_STATES
                    ):
                        energy_populations = np.empty(
                            (
                                n_times,
                                N_STATES,
                            ),
                            dtype=np.float64,
                        )

                        for time_index in range(
                            n_times
                        ):
                            rho = matrix_from_vector(
                                propagated[
                                    time_index,
                                    :,
                                    initial_index,
                                ]
                            )

                            populations[
                                initial_index,
                                time_index,
                            ] = np.real(
                                np.diag(rho)
                            )

                            off_diagonal = (
                                rho
                                - np.diag(
                                    np.diag(rho)
                                )
                            )

                            coherence_l1[
                                initial_index,
                                time_index,
                            ] = float(
                                np.sum(
                                    np.abs(
                                        off_diagonal
                                    )
                                )
                            )

                            for (
                                state_index,
                                projector,
                            ) in enumerate(
                                energy_projectors
                            ):
                                energy_populations[
                                    time_index,
                                    state_index,
                                ] = float(
                                    np.real(
                                        np.trace(
                                            projector
                                            @ rho
                                        )
                                    )
                                )

                        for time_index in (
                            validation_indices
                        ):
                            rho = matrix_from_vector(
                                propagated[
                                    time_index,
                                    :,
                                    initial_index,
                                ]
                            )

                            trace_error = abs(
                                np.trace(rho)
                                - 1.0
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
                                float(
                                    trace_error
                                ),
                            )

                            maximum_hermiticity_error = max(
                                maximum_hermiticity_error,
                                hermiticity_error,
                            )

                            minimum_density_eigenvalue = min(
                                minimum_density_eigenvalue,
                                minimum_eigenvalue,
                            )

                        downward_flux_by_site = {
                            site_index: 0.0
                            for site_index
                            in HIGH_SITE_INDICES
                        }

                        upward_flux_by_site = {
                            site_index: 0.0
                            for site_index
                            in HIGH_SITE_INDICES
                        }

                        site_probabilities = (
                            np.abs(
                                eigenvectors
                            ) ** 2
                        )

                        for high_state in (
                            high_eigenstates
                        ):
                            high_site = (
                                eigenstate_to_site[
                                    high_state
                                ]
                            )

                            gap_eV = float(
                                eigenvalues[
                                    high_state
                                ]
                                - eigenvalues[
                                    sink_eigenstate
                                ]
                            )

                            overlap_weight = float(
                                np.sum(
                                    site_probabilities[
                                        :,
                                        sink_eigenstate,
                                    ]
                                    * site_probabilities[
                                        :,
                                        high_state,
                                    ]
                                )
                            )

                            n_bose = bose_population(
                                gap_eV,
                                float(
                                    temperature_K
                                ),
                            )

                            downward_rate = (
                                float(
                                    kappa_ref
                                )
                                * overlap_weight
                                * (
                                    n_bose
                                    + 1.0
                                )
                            )

                            upward_rate = (
                                float(
                                    kappa_ref
                                )
                                * overlap_weight
                                * n_bose
                            )

                            downward_flux_by_site[
                                high_site
                            ] += float(
                                np.trapezoid(
                                    downward_rate
                                    * energy_populations[
                                        :,
                                        high_state,
                                    ],
                                    times_ps,
                                )
                            )

                            upward_flux_by_site[
                                high_site
                            ] += float(
                                np.trapezoid(
                                    upward_rate
                                    * energy_populations[
                                        :,
                                        sink_eigenstate,
                                    ],
                                    times_ps,
                                )
                            )

                        total_downward_flux = sum(
                            downward_flux_by_site.values()
                        )

                        total_upward_flux = sum(
                            upward_flux_by_site.values()
                        )

                        gateway_fractions = {
                            site_index: (
                                downward_flux_by_site[
                                    site_index
                                ]
                                / total_downward_flux
                                if total_downward_flux
                                > 0.0
                                else 0.0
                            )
                            for site_index
                            in HIGH_SITE_INDICES
                        }

                        trajectory = populations[
                            initial_index
                        ]

                        pyr5_trajectory = trajectory[
                            :,
                            PYR5_SITE_INDEX,
                        ]

                        maximum_pyr5_index = int(
                            np.argmax(
                                pyr5_trajectory
                            )
                        )

                        final_population = trajectory[
                            -1
                        ]

                        frame_rows.append(
                            {
                                "temperature_K": float(
                                    temperature_K
                                ),
                                "kappa_ref_ps_inv": float(
                                    kappa_ref
                                ),
                                "gamma_phi_ps_inv": float(
                                    gamma_phi
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
                                        PYR5_SITE_INDEX
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
                                "final_l1_distance_to_steady": float(
                                    np.sum(
                                        np.abs(
                                            final_population
                                            - stationary_site_population
                                        )
                                    )
                                ),
                                "final_l1_distance_to_Gibbs": float(
                                    np.sum(
                                        np.abs(
                                            final_population
                                            - gibbs_site_population
                                        )
                                    )
                                ),
                                "steady_PYR5_population": float(
                                    stationary_site_population[
                                        PYR5_SITE_INDEX
                                    ]
                                ),
                                "integrated_downward_flux_to_PYR5": (
                                    total_downward_flux
                                ),
                                "integrated_upward_flux_from_PYR5": (
                                    total_upward_flux
                                ),
                                "net_thermal_flux_to_PYR5": (
                                    total_downward_flux
                                    - total_upward_flux
                                ),
                                "PYR2_gateway_fraction": (
                                    gateway_fractions[
                                        0
                                    ]
                                ),
                                "PYR3_gateway_fraction": (
                                    gateway_fractions[
                                        1
                                    ]
                                ),
                                "PYR4_gateway_fraction": (
                                    gateway_fractions[
                                        2
                                    ]
                                ),
                            }
                        )

                    ensemble_sum[
                        temperature_index,
                        kappa_index,
                        gamma_index,
                    ] += populations

                    ensemble_sum_sq[
                        temperature_index,
                        kappa_index,
                        gamma_index,
                    ] += populations**2

                    coherence_sum[
                        temperature_index,
                        kappa_index,
                        gamma_index,
                    ] += coherence_l1

                log(
                    f"[T={temperature_K:g} K, "
                    f"kappa={kappa_ref:g}, "
                    f"gamma={gamma_phi:g}] "
                    f"completed "
                    f"{EXPECTED_FRAMES}/"
                    f"{EXPECTED_FRAMES} frames"
                )

    ensemble_mean = (
        ensemble_sum
        / EXPECTED_FRAMES
    )

    ensemble_variance = (
        ensemble_sum_sq
        / EXPECTED_FRAMES
        - ensemble_mean**2
    )

    ensemble_variance = np.maximum(
        ensemble_variance,
        0.0,
    )

    ensemble_sd = np.sqrt(
        ensemble_variance
    )

    coherence_mean = (
        coherence_sum
        / EXPECTED_FRAMES
    )

    for (
        temperature_index,
        temperature_K,
    ) in enumerate(TEMPERATURES_K):

        reference_temperature_matches = (
            np.where(
                np.isclose(
                    reference_temperatures,
                    temperature_K,
                )
            )[0]
        )

        if (
            reference_temperature_matches.size
            != 1
        ):
            raise RuntimeError(
                f"Could not locate reference "
                f"T={temperature_K:g} K"
            )

        reference_temperature_index = int(
            reference_temperature_matches[
                0
            ]
        )

        for (
            kappa_index,
            kappa_ref,
        ) in enumerate(
            KAPPA_REF_PS_INV
        ):
            reference_kappa_matches = (
                np.where(
                    np.isclose(
                        reference_kappas,
                        kappa_ref,
                    )
                )[0]
            )

            if (
                reference_kappa_matches.size
                != 1
            ):
                raise RuntimeError(
                    f"Could not locate reference "
                    f"kappa={kappa_ref:g}"
                )

            reference_kappa_index = int(
                reference_kappa_matches[
                    0
                ]
            )

            gamma0_index = int(
                np.where(
                    np.isclose(
                        GAMMA_PHI_PS_INV,
                        0.0,
                    )
                )[0][0]
            )

            gamma0_error = float(
                np.max(
                    np.abs(
                        ensemble_mean[
                            temperature_index,
                            kappa_index,
                            gamma0_index,
                        ]
                        - reference_populations[
                            reference_temperature_index,
                            reference_kappa_index,
                        ]
                    )
                )
            )

            maximum_gamma0_reference_error = max(
                maximum_gamma0_reference_error,
                gamma0_error,
            )

    grouped_frame_rows: dict[
        tuple[float, float, float],
        list[dict[str, object]],
    ] = defaultdict(list)

    grouped_steady_rows: dict[
        tuple[float, float, float],
        list[dict[str, object]],
    ] = defaultdict(list)

    for row in frame_rows:
        key = (
            float(
                row["temperature_K"]
            ),
            float(
                row["kappa_ref_ps_inv"]
            ),
            float(
                row["gamma_phi_ps_inv"]
            ),
        )

        grouped_frame_rows[
            key
        ].append(row)

    for row in steady_rows:
        key = (
            float(
                row["temperature_K"]
            ),
            float(
                row["kappa_ref_ps_inv"]
            ),
            float(
                row["gamma_phi_ps_inv"]
            ),
        )

        grouped_steady_rows[
            key
        ].append(row)

    condition_rows: list[
        dict[str, object]
    ] = []

    for temperature_K in TEMPERATURES_K:
        for kappa_ref in (
            KAPPA_REF_PS_INV
        ):
            for gamma_phi in (
                GAMMA_PHI_PS_INV
            ):
                key = (
                    float(
                        temperature_K
                    ),
                    float(
                        kappa_ref
                    ),
                    float(
                        gamma_phi
                    ),
                )

                rows = (
                    grouped_frame_rows[
                        key
                    ]
                )

                steady_condition_rows = (
                    grouped_steady_rows[
                        key
                    ]
                )

                high_rows = [
                    row
                    for row in rows
                    if row[
                        "initial_state"
                    ] in SITES[:3]
                ]

                downward_fluxes = (
                    np.asarray(
                        [
                            float(
                                row[
                                    "integrated_downward_flux_to_PYR5"
                                ]
                            )
                            for row
                            in high_rows
                        ],
                        dtype=np.float64,
                    )
                )

                weighted_gateway_totals = {
                    site: float(
                        np.sum(
                            [
                                float(
                                    row[
                                        f"{site}_gateway_fraction"
                                    ]
                                )
                                * float(
                                    row[
                                        "integrated_downward_flux_to_PYR5"
                                    ]
                                )
                                for row
                                in high_rows
                            ]
                        )
                    )
                    for site in (
                        "PYR2",
                        "PYR3",
                        "PYR4",
                    )
                }

                total_gateway_flux = float(
                    np.sum(
                        downward_fluxes
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
                        "gamma_phi_ps_inv": float(
                            gamma_phi
                        ),
                        "mean_maximum_PYR5_population_from_high_states": float(
                            np.mean(
                                [
                                    float(
                                        row[
                                            "maximum_PYR5_population"
                                        ]
                                    )
                                    for row
                                    in high_rows
                                ]
                            )
                        ),
                        "mean_final_PYR5_population_from_high_states": float(
                            np.mean(
                                [
                                    float(
                                        row[
                                            "final_PYR5_population"
                                        ]
                                    )
                                    for row
                                    in high_rows
                                ]
                            )
                        ),
                        "mean_steady_PYR5_population": float(
                            np.mean(
                                [
                                    float(
                                        row[
                                            "steady_PYR5_population"
                                        ]
                                    )
                                    for row
                                    in steady_condition_rows
                                ]
                            )
                        ),
                        "mean_final_l1_distance_to_steady_from_high_states": float(
                            np.mean(
                                [
                                    float(
                                        row[
                                            "final_l1_distance_to_steady"
                                        ]
                                    )
                                    for row
                                    in high_rows
                                ]
                            )
                        ),
                        "mean_final_l1_distance_to_Gibbs_from_high_states": float(
                            np.mean(
                                [
                                    float(
                                        row[
                                            "final_l1_distance_to_Gibbs"
                                        ]
                                    )
                                    for row
                                    in high_rows
                                ]
                            )
                        ),
                        "mean_integrated_downward_flux_to_PYR5": float(
                            np.mean(
                                downward_fluxes
                            )
                        ),
                        "mean_net_thermal_flux_to_PYR5": float(
                            np.mean(
                                [
                                    float(
                                        row[
                                            "net_thermal_flux_to_PYR5"
                                        ]
                                    )
                                    for row
                                    in high_rows
                                ]
                            )
                        ),
                        "PYR2_gateway_fraction": (
                            weighted_gateway_totals[
                                "PYR2"
                            ]
                            / total_gateway_flux
                            if total_gateway_flux
                            > 0.0
                            else 0.0
                        ),
                        "PYR3_gateway_fraction": (
                            weighted_gateway_totals[
                                "PYR3"
                            ]
                            / total_gateway_flux
                            if total_gateway_flux
                            > 0.0
                            else 0.0
                        ),
                        "PYR4_gateway_fraction": (
                            weighted_gateway_totals[
                                "PYR4"
                            ]
                            / total_gateway_flux
                            if total_gateway_flux
                            > 0.0
                            else 0.0
                        ),
                        "mean_steady_l1_distance_to_Gibbs": float(
                            np.mean(
                                [
                                    float(
                                        row[
                                            "steady_l1_distance_to_Gibbs"
                                        ]
                                    )
                                    for row
                                    in steady_condition_rows
                                ]
                            )
                        ),
                    }
                )

    write_csv(
        FRAME_METRICS_CSV,
        frame_rows,
    )

    write_csv(
        CONDITION_SUMMARY_CSV,
        condition_rows,
    )

    write_csv(
        STEADY_STATE_CSV,
        steady_rows,
    )

    validation_rows = [
        {
            "n_temperatures": int(
                n_temperature
            ),
            "n_kappa_values": int(
                n_kappa
            ),
            "n_gamma_values": int(
                n_gamma
            ),
            "n_frames": EXPECTED_FRAMES,
            "n_initial_states": N_STATES,
            "n_time_points": n_times,
            "maximum_trace_error": (
                maximum_trace_error
            ),
            "maximum_hermiticity_error": (
                maximum_hermiticity_error
            ),
            "minimum_sampled_density_eigenvalue": (
                minimum_density_eigenvalue
            ),
            "maximum_trace_preservation_error": (
                maximum_trace_preservation_error
            ),
            "maximum_detailed_balance_error": (
                maximum_detailed_balance_error
            ),
            "maximum_relaxation_Gibbs_stationarity_error": (
                maximum_relaxation_gibbs_error
            ),
            "maximum_combined_stationary_residual": (
                maximum_stationary_residual
            ),
            "minimum_combined_stationary_density_eigenvalue": (
                minimum_stationary_eigenvalue
            ),
            "maximum_gamma0_reference_error": (
                maximum_gamma0_reference_error
            ),
            "trace_validation_pass": (
                maximum_trace_error
                <= TRACE_TOL
            ),
            "hermiticity_validation_pass": (
                maximum_hermiticity_error
                <= HERMITICITY_TOL
            ),
            "positivity_validation_pass": (
                minimum_density_eigenvalue
                >= -POSITIVITY_TOL
            ),
            "trace_preservation_validation_pass": (
                maximum_trace_preservation_error
                <= TRACE_PRESERVATION_TOL
            ),
            "detailed_balance_validation_pass": (
                maximum_detailed_balance_error
                <= DETAILED_BALANCE_TOL
            ),
            "relaxation_Gibbs_stationarity_validation_pass": (
                maximum_relaxation_gibbs_error
                <= RELAXATION_GIBBS_TOL
            ),
            "combined_stationary_residual_validation_pass": (
                maximum_stationary_residual
                <= STEADY_STATE_RESIDUAL_TOL
            ),
            "combined_stationary_positivity_validation_pass": (
                minimum_stationary_eigenvalue
                >= -POSITIVITY_TOL
            ),
            "gamma0_reference_validation_pass": (
                maximum_gamma0_reference_error
                <= GAMMA0_REFERENCE_TOL
            ),
        }
    ]

    write_csv(
        VALIDATION_CSV,
        validation_rows,
    )

    np.savez_compressed(
        NPZ_PATH,
        temperatures_K=TEMPERATURES_K,
        kappa_ref_ps_inv=(
            KAPPA_REF_PS_INV
        ),
        gamma_phi_ps_inv=(
            GAMMA_PHI_PS_INV
        ),
        times_ps=times_ps,
        snapshot_times_ps=(
            snapshot_times_ps
        ),
        site_labels=np.asarray(
            SITES
        ),
        ensemble_mean_populations=(
            ensemble_mean
        ),
        ensemble_sd_populations=(
            ensemble_sd
        ),
        ensemble_mean_coherence_l1=(
            coherence_mean
        ),
    )

    overall_pass = all(
        bool(value)
        for key, value
        in validation_rows[0].items()
        if key.endswith("_pass")
    )

    pyr4_gateway_values = [
        float(
            row[
                "PYR4_gateway_fraction"
            ]
        )
        for row in condition_rows
    ]

    steady_pyr5_values = [
        float(
            row[
                "mean_steady_PYR5_population"
            ]
        )
        for row in condition_rows
    ]

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Combined Dephasing "
            "and Relaxation Sensitivity\n\n"
        )

        handle.write("## Scope\n\n")

        handle.write(
            "This calculation combines coherent "
            "dynamics, local site-basis pure "
            "dephasing, and global energy-basis "
            "detailed-balance relaxation.\n\n"
        )

        handle.write(
            "The dephasing and relaxation terms "
            "represent separate phenomenological "
            "environments. For gamma_phi > 0, the "
            "combined stationary state is not "
            "assumed to equal the Gibbs state and "
            "is solved explicitly.\n\n"
        )

        handle.write(
            "## Parameter sweep\n\n"
        )

        handle.write(
            "- Temperatures: "
            + ", ".join(
                f"{value:g}"
                for value
                in TEMPERATURES_K
            )
            + " K.\n"
        )

        handle.write(
            "- Kappa values: "
            + ", ".join(
                f"{value:g}"
                for value
                in KAPPA_REF_PS_INV
            )
            + " ps^-1.\n"
        )

        handle.write(
            "- Gamma values: "
            + ", ".join(
                f"{value:g}"
                for value
                in GAMMA_PHI_PS_INV
            )
            + " ps^-1.\n"
        )

        handle.write(
            f"- Propagation interval: "
            f"0-{T_MAX_PS:g} ps.\n\n"
        )

        handle.write(
            "## Numerical validation\n\n"
        )

        handle.write(
            f"- Maximum trace error: "
            f"{maximum_trace_error:.3e}.\n"
        )

        handle.write(
            f"- Maximum Hermiticity error: "
            f"{maximum_hermiticity_error:.3e}.\n"
        )

        handle.write(
            f"- Minimum sampled density "
            f"eigenvalue: "
            f"{minimum_density_eigenvalue:.3e}.\n"
        )

        handle.write(
            f"- Maximum combined stationary "
            f"residual: "
            f"{maximum_stationary_residual:.3e}.\n"
        )

        handle.write(
            f"- Gamma=0 reference error: "
            f"{maximum_gamma0_reference_error:.3e}.\n"
        )

        handle.write(
            f"- Overall validation: "
            f"{'PASS' if overall_pass else 'FAIL'}."
            "\n\n"
        )

        handle.write(
            "## Aggregate diagnostic ranges\n\n"
        )

        handle.write(
            f"- PYR4 thermal gateway fraction: "
            f"{min(pyr4_gateway_values):.6f} "
            f"to "
            f"{max(pyr4_gateway_values):.6f}.\n"
        )

        handle.write(
            f"- Mean stationary PYR5 population: "
            f"{min(steady_pyr5_values):.6f} "
            f"to "
            f"{max(steady_pyr5_values):.6f}.\n\n"
        )

        handle.write(
            "## Interpretation limit\n\n"
        )

        handle.write(
            "Absolute relaxation times remain "
            "proportional to the phenomenological "
            "kappa scale. Local pure dephasing is "
            "not a thermally balanced population "
            "bath, so deviations of the combined "
            "stationary state from Gibbs are "
            "expected and must be interpreted as "
            "model sensitivity rather than a "
            "microscopic equilibrium prediction.\n"
        )

    log("")
    log(
        "Day020 combined dephasing "
        "+ relaxation completed."
    )

    log(
        f"Temperatures: "
        f"{n_temperature}/"
        f"{n_temperature}"
    )

    log(
        f"Kappa values: "
        f"{n_kappa}/{n_kappa}"
    )

    log(
        f"Gamma values: "
        f"{n_gamma}/{n_gamma}"
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
        "Minimum sampled density eigenvalue: "
        f"{minimum_density_eigenvalue:.3e}"
    )

    log(
        "Maximum combined stationary residual: "
        f"{maximum_stationary_residual:.3e}"
    )

    log(
        "Maximum gamma=0 reference error: "
        f"{maximum_gamma0_reference_error:.3e}"
    )

    log(
        "PYR4 gateway fraction range: "
        f"{min(pyr4_gateway_values):.6f} "
        f"to "
        f"{max(pyr4_gateway_values):.6f}"
    )

    log(
        "Mean stationary PYR5 population range: "
        f"{min(steady_pyr5_values):.6f} "
        f"to "
        f"{max(steady_pyr5_values):.6f}"
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
