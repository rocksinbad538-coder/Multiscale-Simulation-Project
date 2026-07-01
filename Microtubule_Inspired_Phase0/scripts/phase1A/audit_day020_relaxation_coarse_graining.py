#!/usr/bin/env python3

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

try:
    from scipy.optimize import linear_sum_assignment
except ImportError as exc:
    raise SystemExit(
        "SciPy is required for the state-character assignment."
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

PATHWAY_ROOT = (
    RELAXATION_ROOT
    / "relaxation_pathway_audit"
)

APPARENT_CAPTURE_CSV = (
    PATHWAY_ROOT
    / "apparent_capture_times.csv"
)

OUTPUT_ROOT = (
    RELAXATION_ROOT
    / "relaxation_coarse_graining_audit"
)

FRAME_CSV = OUTPUT_ROOT / "frame_coarse_graining_metrics.csv"
SUMMARY_CSV = OUTPUT_ROOT / "coarse_grained_capture_summary.csv"
REPORT_MD = OUTPUT_ROOT / "RELAXATION_COARSE_GRAINING_DAY020.md"

EXPECTED_FRAMES = 21
N_STATES = 4

SITES = (
    "PYR2_bright",
    "PYR3_bright",
    "PYR4_bright",
    "PYR5_bright",
)

PYR5_SITE_INDEX = 3

TEMPERATURES_K = (
    150.0,
    200.0,
    250.0,
    300.0,
)

KAPPA_VALUES_PS_INV = (
    0.1,
    1.0,
    10.0,
)

KB_EV_K = 8.617333262145e-5
FRAME_RE = re.compile(r"frame=(\d+)")


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


def bose_population(
    gap_eV: float,
    temperature_K: float,
) -> float:
    exponent = gap_eV / (
        KB_EV_K * temperature_K
    )

    if exponent > 700.0:
        return 0.0

    return float(1.0 / np.expm1(exponent))


def read_hamiltonians() -> np.ndarray:
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

    matrices: list[np.ndarray] = []

    for expected_frame, path in enumerate(files):
        frame: int | None = None

        for line in path.read_text(
            encoding="utf-8"
        ).splitlines():
            match = FRAME_RE.search(line)

            if match is not None:
                frame = int(match.group(1))
                break

        if frame != expected_frame:
            raise RuntimeError(
                f"Unexpected frame in {path}: {frame}"
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
                f"Unexpected shape in {path}: "
                f"{matrix.shape}"
            )

        if not np.allclose(
            matrix,
            matrix.T,
            atol=1.0e-12,
            rtol=0.0,
        ):
            raise RuntimeError(
                f"Non-symmetric Hamiltonian: {path}"
            )

        matrices.append(matrix)

    return np.stack(matrices)


def build_rate_generator(
    energies_eV: np.ndarray,
    eigenvectors: np.ndarray,
    temperature_K: float,
) -> tuple[np.ndarray, float]:
    site_weights = np.abs(eigenvectors) ** 2

    generator = np.zeros(
        (N_STATES, N_STATES),
        dtype=np.float64,
    )

    maximum_detailed_balance_error = 0.0

    for low in range(N_STATES):
        for high in range(
            low + 1,
            N_STATES,
        ):
            gap_eV = float(
                energies_eV[high]
                - energies_eV[low]
            )

            overlap_weight = float(
                np.sum(
                    site_weights[:, low]
                    * site_weights[:, high]
                )
            )

            n_bose = bose_population(
                gap_eV,
                temperature_K,
            )

            downward_rate = (
                overlap_weight
                * (n_bose + 1.0)
            )

            upward_rate = (
                overlap_weight
                * n_bose
            )

            # Column convention:
            # generator[destination, source].
            generator[low, high] += downward_rate
            generator[high, low] += upward_rate

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

    for source in range(N_STATES):
        generator[source, source] = -float(
            np.sum(generator[:, source])
        )

    return (
        generator,
        maximum_detailed_balance_error,
    )


def spectral_gap(
    generator: np.ndarray,
) -> float:
    eigenvalues = np.linalg.eigvals(generator)

    rates = [
        -float(value.real)
        for value in eigenvalues
        if value.real < -1.0e-14
    ]

    if not rates:
        raise RuntimeError(
            "No nonzero decay rate found"
        )

    return min(rates)


def internal_generator(
    full_generator: np.ndarray,
    high_indices: list[int],
) -> np.ndarray:
    size = len(high_indices)

    reduced = np.zeros(
        (size, size),
        dtype=np.float64,
    )

    for destination_position, destination in enumerate(
        high_indices
    ):
        for source_position, source in enumerate(
            high_indices
        ):
            if destination == source:
                continue

            reduced[
                destination_position,
                source_position,
            ] = full_generator[
                destination,
                source,
            ]

    for source_position in range(size):
        reduced[
            source_position,
            source_position,
        ] = -float(
            np.sum(
                reduced[:, source_position]
            )
        )

    return reduced


def read_apparent_capture_times() -> dict[
    tuple[float, float],
    float,
]:
    values: dict[
        tuple[float, float],
        float,
    ] = {}

    with APPARENT_CAPTURE_CSV.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(handle)

        for row in reader:
            key = (
                float(row["temperature_K"]),
                float(row["kappa_ref_ps_inv"]),
            )

            values[key] = float(
                row[
                    "apparent_capture_timescale_ns"
                ]
            )

    return values


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    hamiltonians = read_hamiltonians()
    apparent_times = read_apparent_capture_times()

    frame_rows: list[dict[str, object]] = []

    maximum_column_sum_error = 0.0
    maximum_detailed_balance_error = 0.0
    minimum_assigned_weight = 1.0

    for frame, hamiltonian in enumerate(
        hamiltonians
    ):
        energies, eigenvectors = np.linalg.eigh(
            hamiltonian
        )

        site_weights = np.abs(
            eigenvectors
        ) ** 2

        site_indices, eigenstate_indices = (
            linear_sum_assignment(
                -site_weights
            )
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

        assigned_weights = [
            float(
                site_weights[
                    site,
                    eigenstate,
                ]
            )
            for site, eigenstate
            in site_to_eigenstate.items()
        ]

        minimum_assigned_weight = min(
            minimum_assigned_weight,
            min(assigned_weights),
        )

        sink_index = site_to_eigenstate[
            PYR5_SITE_INDEX
        ]

        high_indices = [
            index
            for index in range(N_STATES)
            if index != sink_index
        ]

        for temperature_K in TEMPERATURES_K:
            generator, detailed_balance_error = (
                build_rate_generator(
                    energies,
                    eigenvectors,
                    temperature_K,
                )
            )

            maximum_detailed_balance_error = max(
                maximum_detailed_balance_error,
                detailed_balance_error,
            )

            column_sum_error = float(
                np.max(
                    np.abs(
                        np.sum(
                            generator,
                            axis=0,
                        )
                    )
                )
            )

            maximum_column_sum_error = max(
                maximum_column_sum_error,
                column_sum_error,
            )

            shifted_high_energies = (
                energies[high_indices]
                - np.min(
                    energies[high_indices]
                )
            )

            conditional_weights = np.exp(
                -shifted_high_energies
                / (
                    KB_EV_K
                    * temperature_K
                )
            )

            conditional_weights /= np.sum(
                conditional_weights
            )

            capture_rates = np.asarray(
                [
                    generator[
                        sink_index,
                        high_index,
                    ]
                    for high_index in high_indices
                ],
                dtype=np.float64,
            )

            contributions = (
                conditional_weights
                * capture_rates
            )

            effective_capture_rate = float(
                np.sum(contributions)
            )

            if effective_capture_rate <= 0.0:
                raise RuntimeError(
                    "Non-positive effective capture rate"
                )

            contribution_fractions = (
                contributions
                / effective_capture_rate
            )

            contribution_by_site = {
                site: 0.0
                for site in SITES[:3]
            }

            conditional_weight_by_site = {
                site: 0.0
                for site in SITES[:3]
            }

            direct_rate_by_site = {
                site: 0.0
                for site in SITES[:3]
            }

            for position, high_index in enumerate(
                high_indices
            ):
                site_index = eigenstate_to_site[
                    high_index
                ]
                site_label = SITES[site_index]

                contribution_by_site[
                    site_label
                ] = float(
                    contribution_fractions[
                        position
                    ]
                )

                conditional_weight_by_site[
                    site_label
                ] = float(
                    conditional_weights[
                        position
                    ]
                )

                direct_rate_by_site[
                    site_label
                ] = float(
                    capture_rates[position]
                )

            high_generator = internal_generator(
                generator,
                high_indices,
            )

            internal_relaxation_rate = spectral_gap(
                high_generator
            )

            exact_full_slow_rate = spectral_gap(
                generator
            )

            relative_slow_mode_error = abs(
                effective_capture_rate
                - exact_full_slow_rate
            ) / exact_full_slow_rate

            frame_rows.append(
                {
                    "temperature_K": temperature_K,
                    "frame": frame,
                    "effective_capture_rate_at_kappa1_ps_inv": (
                        effective_capture_rate
                    ),
                    "effective_capture_timescale_at_kappa1_ps": (
                        1.0 / effective_capture_rate
                    ),
                    "exact_full_slow_rate_at_kappa1_ps_inv": (
                        exact_full_slow_rate
                    ),
                    "exact_full_slow_timescale_at_kappa1_ps": (
                        1.0 / exact_full_slow_rate
                    ),
                    "internal_relaxation_rate_at_kappa1_ps_inv": (
                        internal_relaxation_rate
                    ),
                    "internal_relaxation_timescale_at_kappa1_ps": (
                        1.0 / internal_relaxation_rate
                    ),
                    "timescale_separation_ratio": (
                        internal_relaxation_rate
                        / effective_capture_rate
                    ),
                    "relative_effective_vs_exact_slow_mode_error": (
                        relative_slow_mode_error
                    ),
                    "PYR2_conditional_high_population": (
                        conditional_weight_by_site[
                            "PYR2_bright"
                        ]
                    ),
                    "PYR3_conditional_high_population": (
                        conditional_weight_by_site[
                            "PYR3_bright"
                        ]
                    ),
                    "PYR4_conditional_high_population": (
                        conditional_weight_by_site[
                            "PYR4_bright"
                        ]
                    ),
                    "PYR2_direct_capture_rate_at_kappa1_ps_inv": (
                        direct_rate_by_site[
                            "PYR2_bright"
                        ]
                    ),
                    "PYR3_direct_capture_rate_at_kappa1_ps_inv": (
                        direct_rate_by_site[
                            "PYR3_bright"
                        ]
                    ),
                    "PYR4_direct_capture_rate_at_kappa1_ps_inv": (
                        direct_rate_by_site[
                            "PYR4_bright"
                        ]
                    ),
                    "PYR2_capture_contribution_fraction": (
                        contribution_by_site[
                            "PYR2_bright"
                        ]
                    ),
                    "PYR3_capture_contribution_fraction": (
                        contribution_by_site[
                            "PYR3_bright"
                        ]
                    ),
                    "PYR4_capture_contribution_fraction": (
                        contribution_by_site[
                            "PYR4_bright"
                        ]
                    ),
                }
            )

    grouped: dict[
        float,
        list[dict[str, object]],
    ] = defaultdict(list)

    for row in frame_rows:
        grouped[
            float(row["temperature_K"])
        ].append(row)

    summary_rows: list[dict[str, object]] = []

    for temperature_K in TEMPERATURES_K:
        rows = grouped[temperature_K]

        unit_effective_rates = np.asarray(
            [
                float(
                    row[
                        "effective_capture_rate_at_kappa1_ps_inv"
                    ]
                )
                for row in rows
            ],
            dtype=np.float64,
        )

        unit_exact_rates = np.asarray(
            [
                float(
                    row[
                        "exact_full_slow_rate_at_kappa1_ps_inv"
                    ]
                )
                for row in rows
            ],
            dtype=np.float64,
        )

        unit_internal_rates = np.asarray(
            [
                float(
                    row[
                        "internal_relaxation_rate_at_kappa1_ps_inv"
                    ]
                )
                for row in rows
            ],
            dtype=np.float64,
        )

        contribution_totals = {
            site: float(
                np.sum(
                    [
                        float(
                            row[
                                f"{site}_capture_contribution_fraction"
                            ]
                        )
                        * float(
                            row[
                                "effective_capture_rate_at_kappa1_ps_inv"
                            ]
                        )
                        for row in rows
                    ]
                )
            )
            for site in (
                "PYR2",
                "PYR3",
                "PYR4",
            )
        }

        total_effective_rate = float(
            np.sum(unit_effective_rates)
        )

        for kappa in KAPPA_VALUES_PS_INV:
            mean_effective_rate = float(
                np.mean(unit_effective_rates)
                * kappa
            )

            mean_exact_rate = float(
                np.mean(unit_exact_rates)
                * kappa
            )

            mean_internal_rate = float(
                np.mean(unit_internal_rates)
                * kappa
            )

            summary_rows.append(
                {
                    "temperature_K": temperature_K,
                    "kappa_ref_ps_inv": kappa,
                    "mean_effective_capture_rate_ps_inv": (
                        mean_effective_rate
                    ),
                    "inverse_mean_effective_rate_ns": (
                        1.0
                        / mean_effective_rate
                        / 1000.0
                    ),
                    "mean_frame_effective_capture_time_ns": float(
                        np.mean(
                            1.0
                            / (
                                unit_effective_rates
                                * kappa
                            )
                        )
                        / 1000.0
                    ),
                    "mean_exact_full_slow_rate_ps_inv": (
                        mean_exact_rate
                    ),
                    "inverse_mean_exact_full_slow_rate_ns": (
                        1.0
                        / mean_exact_rate
                        / 1000.0
                    ),
                    "mean_internal_relaxation_rate_ps_inv": (
                        mean_internal_rate
                    ),
                    "mean_internal_relaxation_time_ps": (
                        1.0 / mean_internal_rate
                    ),
                    "mean_timescale_separation_ratio": float(
                        np.mean(
                            unit_internal_rates
                            / unit_effective_rates
                        )
                    ),
                    "maximum_relative_effective_vs_exact_slow_mode_error": float(
                        np.max(
                            np.abs(
                                unit_effective_rates
                                - unit_exact_rates
                            )
                            / unit_exact_rates
                        )
                    ),
                    "PYR2_gateway_contribution_fraction": (
                        contribution_totals["PYR2"]
                        / total_effective_rate
                    ),
                    "PYR3_gateway_contribution_fraction": (
                        contribution_totals["PYR3"]
                        / total_effective_rate
                    ),
                    "PYR4_gateway_contribution_fraction": (
                        contribution_totals["PYR4"]
                        / total_effective_rate
                    ),
                    "apparent_100ps_capture_time_ns": (
                        apparent_times[
                            (
                                temperature_K,
                                kappa,
                            )
                        ]
                    ),
                }
            )

    write_csv(
        FRAME_CSV,
        frame_rows,
    )

    write_csv(
        SUMMARY_CSV,
        summary_rows,
    )

    minimum_separation = min(
        float(
            row[
                "timescale_separation_ratio"
            ]
        )
        for row in frame_rows
    )

    maximum_slow_mode_error = max(
        float(
            row[
                "relative_effective_vs_exact_slow_mode_error"
            ]
        )
        for row in frame_rows
    )

    pyr4_contributions = [
        float(
            row[
                "PYR4_gateway_contribution_fraction"
            ]
        )
        for row in summary_rows
    ]

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Relaxation Coarse-Graining Audit\n\n"
        )

        handle.write("## Purpose\n\n")
        handle.write(
            "This audit tests whether the rapidly mixing "
            "PYR2-PYR4 manifold can be represented as a "
            "single upper reservoir coupled to the low-energy "
            "PYR5 sink.\n\n"
        )

        handle.write("## Validation\n\n")
        handle.write(
            f"- Hamiltonians: {len(hamiltonians)}/"
            f"{EXPECTED_FRAMES}.\n"
        )
        handle.write(
            f"- Minimum assigned site weight: "
            f"{minimum_assigned_weight:.8f}.\n"
        )
        handle.write(
            f"- Maximum generator column-sum error: "
            f"{maximum_column_sum_error:.3e}.\n"
        )
        handle.write(
            f"- Maximum detailed-balance error: "
            f"{maximum_detailed_balance_error:.3e}.\n\n"
        )

        handle.write("## Coarse-graining result\n\n")
        handle.write(
            f"- Minimum internal-mixing/capture rate ratio: "
            f"{minimum_separation:.3e}.\n"
        )
        handle.write(
            f"- Maximum relative error between the "
            f"coarse-grained capture rate and the exact slow "
            f"mode: {maximum_slow_mode_error:.3e}.\n"
        )
        handle.write(
            f"- PYR4 gateway contribution range: "
            f"{min(pyr4_contributions):.6f} to "
            f"{max(pyr4_contributions):.6f}.\n\n"
        )

        handle.write("## Interpretation\n\n")
        handle.write(
            "A large separation ratio supports rapid "
            "pre-equilibration within PYR2-PYR4 followed by "
            "slow leakage toward PYR5. Agreement with the "
            "exact slow mode validates the reduced "
            "upper-reservoir plus PYR5-sink description.\n\n"
        )

        handle.write(
            "All absolute capture times remain proportional "
            "to the phenomenological kappa parameter and are "
            "not microscopic predictions.\n"
        )

    log(
        "Day020 relaxation coarse-graining audit completed."
    )
    log(
        f"Hamiltonians: "
        f"{len(hamiltonians)}/{EXPECTED_FRAMES}"
    )
    log(
        f"Minimum assigned site weight: "
        f"{minimum_assigned_weight:.8f}"
    )
    log(
        f"Maximum generator column-sum error: "
        f"{maximum_column_sum_error:.3e}"
    )
    log(
        f"Maximum detailed-balance error: "
        f"{maximum_detailed_balance_error:.3e}"
    )
    log(
        f"Minimum internal-mixing/capture rate ratio: "
        f"{minimum_separation:.3e}"
    )
    log(
        f"Maximum coarse-grained slow-mode error: "
        f"{maximum_slow_mode_error:.3e}"
    )
    log(
        f"PYR4 gateway contribution range: "
        f"{min(pyr4_contributions):.6f} to "
        f"{max(pyr4_contributions):.6f}"
    )
    log(
        f"Wrote: "
        f"{OUTPUT_ROOT.relative_to(PROJECT_ROOT)}"
    )


if __name__ == "__main__":
    main()
