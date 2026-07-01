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
        "SciPy is required for the energy-state character assignment."
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

OUTPUT_ROOT = RELAXATION_ROOT / "relaxation_pathway_audit"

CHARACTER_RAW_CSV = OUTPUT_ROOT / "energy_state_character_raw.csv"
CHARACTER_SUMMARY_CSV = OUTPUT_ROOT / "energy_state_character_summary.csv"
PAIR_RAW_CSV = OUTPUT_ROOT / "thermal_pair_rates_raw.csv"
PAIR_SUMMARY_CSV = OUTPUT_ROOT / "thermal_pair_rate_summary.csv"
PYR5_PATHWAYS_CSV = OUTPUT_ROOT / "PYR5_pathway_summary.csv"
CAPTURE_CSV = OUTPUT_ROOT / "apparent_capture_times.csv"
REPORT_MD = OUTPUT_ROOT / "RELAXATION_PATHWAY_AUDIT_DAY020.md"

CONDITION_SUMMARY_CSV = RELAXATION_ROOT / "condition_summary.csv"

EXPECTED_FRAMES = 21
N_STATES = 4

SITES = (
    "PYR2_bright",
    "PYR3_bright",
    "PYR4_bright",
    "PYR5_bright",
)

SITE_ORDER = {
    label: index
    for index, label in enumerate(SITES)
}

TEMPERATURES_K = (
    150.0,
    200.0,
    250.0,
    300.0,
)

KB_EV_K = 8.617333262145e-5
KAPPA_REFERENCE_PS_INV = 1.0
PROPAGATION_TIME_PS = 100.0

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
                f"Unexpected matrix shape in {path}: "
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


def canonical_pair(
    label_a: str,
    label_b: str,
) -> tuple[str, str]:
    if SITE_ORDER[label_a] < SITE_ORDER[label_b]:
        return label_a, label_b

    return label_b, label_a


def mean(values: list[float]) -> float:
    return float(
        np.mean(
            np.asarray(values, dtype=np.float64)
        )
    )


def sd(values: list[float]) -> float:
    return float(
        np.std(
            np.asarray(values, dtype=np.float64)
        )
    )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    hamiltonians = read_hamiltonians()

    character_rows: list[dict[str, object]] = []
    pair_rows: list[dict[str, object]] = []

    maximum_assignment_error = 0.0
    maximum_detailed_balance_error = 0.0
    minimum_assigned_weight = 1.0

    for frame_index, hamiltonian in enumerate(
        hamiltonians
    ):
        energies, eigenvectors = np.linalg.eigh(
            hamiltonian
        )

        site_weights = (
            np.abs(eigenvectors) ** 2
        )

        row_indices, column_indices = (
            linear_sum_assignment(
                -site_weights
            )
        )

        site_to_eigenstate = {
            int(site): int(eigenstate)
            for site, eigenstate in zip(
                row_indices,
                column_indices,
            )
        }

        if len(site_to_eigenstate) != N_STATES:
            raise RuntimeError(
                "Energy-state assignment is incomplete"
            )

        eigenstate_to_site = {
            eigenstate: site
            for site, eigenstate
            in site_to_eigenstate.items()
        }

        assignment_total = sum(
            site_weights[
                site,
                eigenstate,
            ]
            for site, eigenstate
            in site_to_eigenstate.items()
        )

        maximum_assignment_error = max(
            maximum_assignment_error,
            abs(
                assignment_total
                - sum(
                    max(site_weights[:, eigenstate])
                    for eigenstate
                    in range(N_STATES)
                )
            ),
        )

        for site_index, site_label in enumerate(
            SITES
        ):
            eigenstate_index = (
                site_to_eigenstate[site_index]
            )

            weights = site_weights[
                :,
                eigenstate_index,
            ]

            assigned_weight = float(
                weights[site_index]
            )

            minimum_assigned_weight = min(
                minimum_assigned_weight,
                assigned_weight,
            )

            participation_ratio = float(
                1.0 / np.sum(weights**2)
            )

            character_rows.append(
                {
                    "frame": frame_index,
                    "site_character": site_label,
                    "energy_eigenstate_index": (
                        eigenstate_index
                    ),
                    "energy_eV": float(
                        energies[
                            eigenstate_index
                        ]
                    ),
                    "energy_relative_to_minimum_meV": float(
                        (
                            energies[
                                eigenstate_index
                            ]
                            - energies[0]
                        )
                        * 1000.0
                    ),
                    "assigned_site_weight": (
                        assigned_weight
                    ),
                    "participation_ratio": (
                        participation_ratio
                    ),
                    "PYR2_weight": float(
                        weights[0]
                    ),
                    "PYR3_weight": float(
                        weights[1]
                    ),
                    "PYR4_weight": float(
                        weights[2]
                    ),
                    "PYR5_weight": float(
                        weights[3]
                    ),
                }
            )

        for temperature_K in TEMPERATURES_K:
            for low_index in range(N_STATES):
                for high_index in range(
                    low_index + 1,
                    N_STATES,
                ):
                    low_site_index = (
                        eigenstate_to_site[
                            low_index
                        ]
                    )
                    high_site_index = (
                        eigenstate_to_site[
                            high_index
                        ]
                    )

                    low_label = SITES[
                        low_site_index
                    ]
                    high_label = SITES[
                        high_site_index
                    ]

                    site_i, site_j = canonical_pair(
                        low_label,
                        high_label,
                    )

                    gap_eV = float(
                        energies[high_index]
                        - energies[low_index]
                    )

                    overlap_weight = float(
                        np.sum(
                            site_weights[
                                :,
                                low_index,
                            ]
                            * site_weights[
                                :,
                                high_index,
                            ]
                        )
                    )

                    n_bose = bose_population(
                        gap_eV,
                        temperature_K,
                    )

                    downward_rate = (
                        KAPPA_REFERENCE_PS_INV
                        * overlap_weight
                        * (n_bose + 1.0)
                    )

                    upward_rate = (
                        KAPPA_REFERENCE_PS_INV
                        * overlap_weight
                        * n_bose
                    )

                    if (
                        upward_rate > 0.0
                        and downward_rate > 0.0
                    ):
                        detailed_balance_error = abs(
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
                            float(
                                detailed_balance_error
                            ),
                        )

                    pair_rows.append(
                        {
                            "temperature_K": (
                                temperature_K
                            ),
                            "frame": frame_index,
                            "site_i": site_i,
                            "site_j": site_j,
                            "lower_energy_character": (
                                low_label
                            ),
                            "higher_energy_character": (
                                high_label
                            ),
                            "gap_eV": gap_eV,
                            "gap_meV": (
                                1000.0 * gap_eV
                            ),
                            "site_projector_overlap_weight": (
                                overlap_weight
                            ),
                            "bose_population": n_bose,
                            "downward_rate_at_kappa1_ps_inv": (
                                downward_rate
                            ),
                            "upward_rate_at_kappa1_ps_inv": (
                                upward_rate
                            ),
                        }
                    )

    character_groups: dict[
        str,
        list[dict[str, object]],
    ] = defaultdict(list)

    for row in character_rows:
        character_groups[
            str(row["site_character"])
        ].append(row)

    character_summary_rows: list[
        dict[str, object]
    ] = []

    for site_label in SITES:
        rows = character_groups[site_label]

        assigned_weights = [
            float(row["assigned_site_weight"])
            for row in rows
        ]

        participation_ratios = [
            float(row["participation_ratio"])
            for row in rows
        ]

        relative_energies = [
            float(
                row[
                    "energy_relative_to_minimum_meV"
                ]
            )
            for row in rows
        ]

        character_summary_rows.append(
            {
                "site_character": site_label,
                "n_frames": len(rows),
                "mean_assigned_site_weight": mean(
                    assigned_weights
                ),
                "minimum_assigned_site_weight": min(
                    assigned_weights
                ),
                "maximum_assigned_site_weight": max(
                    assigned_weights
                ),
                "mean_participation_ratio": mean(
                    participation_ratios
                ),
                "maximum_participation_ratio": max(
                    participation_ratios
                ),
                "mean_energy_relative_to_minimum_meV": mean(
                    relative_energies
                ),
                "sd_energy_relative_to_minimum_meV": sd(
                    relative_energies
                ),
            }
        )

    pair_groups: dict[
        tuple[float, str, str],
        list[dict[str, object]],
    ] = defaultdict(list)

    for row in pair_rows:
        key = (
            float(row["temperature_K"]),
            str(row["site_i"]),
            str(row["site_j"]),
        )
        pair_groups[key].append(row)

    pair_summary_rows: list[
        dict[str, object]
    ] = []

    for key in sorted(pair_groups):
        temperature_K, site_i, site_j = key
        rows = pair_groups[key]

        gaps = [
            float(row["gap_meV"])
            for row in rows
        ]

        overlaps = [
            float(
                row[
                    "site_projector_overlap_weight"
                ]
            )
            for row in rows
        ]

        downward_rates = [
            float(
                row[
                    "downward_rate_at_kappa1_ps_inv"
                ]
            )
            for row in rows
        ]

        upward_rates = [
            float(
                row[
                    "upward_rate_at_kappa1_ps_inv"
                ]
            )
            for row in rows
        ]

        mean_downward_rate = mean(
            downward_rates
        )

        mean_upward_rate = mean(
            upward_rates
        )

        pair_summary_rows.append(
            {
                "temperature_K": temperature_K,
                "site_i": site_i,
                "site_j": site_j,
                "n_frames": len(rows),
                "fraction_site_i_is_lower_energy": mean(
                    [
                        1.0
                        if row[
                            "lower_energy_character"
                        ]
                        == site_i
                        else 0.0
                        for row in rows
                    ]
                ),
                "mean_gap_meV": mean(gaps),
                "sd_gap_meV": sd(gaps),
                "mean_site_projector_overlap_weight": mean(
                    overlaps
                ),
                "minimum_overlap_weight": min(
                    overlaps
                ),
                "maximum_overlap_weight": max(
                    overlaps
                ),
                "mean_downward_rate_at_kappa1_ps_inv": (
                    mean_downward_rate
                ),
                "mean_upward_rate_at_kappa1_ps_inv": (
                    mean_upward_rate
                ),
                "downward_timescale_at_kappa1_ps": (
                    1.0 / mean_downward_rate
                    if mean_downward_rate > 0.0
                    else float("inf")
                ),
                "downward_timescale_at_kappa10_ps": (
                    1.0
                    / (
                        10.0
                        * mean_downward_rate
                    )
                    if mean_downward_rate > 0.0
                    else float("inf")
                ),
            }
        )

    pyr5_pathway_rows = [
        row
        for row in pair_summary_rows
        if (
            row["site_i"] == "PYR5_bright"
            or row["site_j"] == "PYR5_bright"
        )
    ]

    capture_rows: list[dict[str, object]] = []

    with CONDITION_SUMMARY_CSV.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(handle)

        for row in reader:
            temperature_K = float(
                row["temperature_K"]
            )
            kappa = float(
                row["kappa_ref_ps_inv"]
            )
            final_pyr5 = float(
                row[
                    "mean_final_PYR5_population_from_high_states"
                ]
            )
            equilibrium_pyr5 = float(
                row[
                    "ensemble_Gibbs_PYR5_population"
                ]
            )

            normalized_fraction = (
                final_pyr5 / equilibrium_pyr5
            )

            normalized_fraction = min(
                max(normalized_fraction, 0.0),
                1.0 - 1.0e-15,
            )

            apparent_rate = (
                -np.log(
                    1.0 - normalized_fraction
                )
                / PROPAGATION_TIME_PS
            )

            capture_rows.append(
                {
                    "temperature_K": temperature_K,
                    "kappa_ref_ps_inv": kappa,
                    "final_PYR5_population": (
                        final_pyr5
                    ),
                    "equilibrium_PYR5_population": (
                        equilibrium_pyr5
                    ),
                    "fraction_of_equilibrium_reached": (
                        normalized_fraction
                    ),
                    "apparent_capture_rate_ps_inv": (
                        apparent_rate
                    ),
                    "apparent_capture_timescale_ps": (
                        1.0 / apparent_rate
                        if apparent_rate > 0.0
                        else float("inf")
                    ),
                    "apparent_capture_timescale_ns": (
                        1.0
                        / apparent_rate
                        / 1000.0
                        if apparent_rate > 0.0
                        else float("inf")
                    ),
                }
            )

    write_csv(
        CHARACTER_RAW_CSV,
        character_rows,
    )
    write_csv(
        CHARACTER_SUMMARY_CSV,
        character_summary_rows,
    )
    write_csv(
        PAIR_RAW_CSV,
        pair_rows,
    )
    write_csv(
        PAIR_SUMMARY_CSV,
        pair_summary_rows,
    )
    write_csv(
        PYR5_PATHWAYS_CSV,
        pyr5_pathway_rows,
    )
    write_csv(
        CAPTURE_CSV,
        capture_rows,
    )

    pyr5_overlaps = [
        float(
            row[
                "mean_site_projector_overlap_weight"
            ]
        )
        for row in pyr5_pathway_rows
    ]

    pyr5_times_kappa10 = [
        float(
            row[
                "downward_timescale_at_kappa10_ps"
            ]
        )
        for row in pyr5_pathway_rows
    ]

    capture_times_ns = [
        float(
            row[
                "apparent_capture_timescale_ns"
            ]
        )
        for row in capture_rows
    ]

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Relaxation-Pathway Audit\n\n"
        )

        handle.write("## Purpose\n\n")
        handle.write(
            "This audit identifies the kinetic bottleneck "
            "in the phenomenological detailed-balance "
            "relaxation model.\n\n"
        )

        handle.write("## State character\n\n")
        handle.write(
            f"- Minimum assigned site weight: "
            f"{minimum_assigned_weight:.8f}.\n"
        )
        handle.write(
            "- Energy eigenstates were assigned one-to-one "
            "to the four bright-site characters.\n\n"
        )

        handle.write("## PYR5 pathways\n\n")
        handle.write(
            f"- Mean overlap-weight range: "
            f"{min(pyr5_overlaps):.6e} to "
            f"{max(pyr5_overlaps):.6e}.\n"
        )
        handle.write(
            f"- Direct downward timescale range at "
            f"kappa=10 ps^-1: "
            f"{min(pyr5_times_kappa10):.3f} to "
            f"{max(pyr5_times_kappa10):.3f} ps.\n\n"
        )

        handle.write("## Apparent ensemble capture\n\n")
        handle.write(
            f"- Apparent capture-time range across the "
            f"tested conditions: "
            f"{min(capture_times_ns):.3f} to "
            f"{max(capture_times_ns):.3f} ns.\n\n"
        )

        handle.write("## Interpretation limit\n\n")
        handle.write(
            "The reported times depend linearly on the "
            "phenomenological kappa scale and on the assumed "
            "energy-independent bath amplitude. They are not "
            "microscopic relaxation times.\n"
        )

    log(
        "Day020 relaxation-pathway audit completed."
    )
    log(
        f"Hamiltonians: "
        f"{hamiltonians.shape[0]}/"
        f"{EXPECTED_FRAMES}"
    )
    log(
        f"Minimum assigned site weight: "
        f"{minimum_assigned_weight:.8f}"
    )
    log(
        f"Maximum detailed-balance error: "
        f"{maximum_detailed_balance_error:.3e}"
    )
    log(
        "PYR5 overlap-weight range: "
        f"{min(pyr5_overlaps):.6e} to "
        f"{max(pyr5_overlaps):.6e}"
    )
    log(
        "PYR5 downward-timescale range at "
        "kappa=10 ps^-1: "
        f"{min(pyr5_times_kappa10):.3f} to "
        f"{max(pyr5_times_kappa10):.3f} ps"
    )
    log(
        "Apparent ensemble capture-time range: "
        f"{min(capture_times_ns):.3f} to "
        f"{max(capture_times_ns):.3f} ns"
    )
    log(
        f"Wrote: "
        f"{OUTPUT_ROOT.relative_to(PROJECT_ROOT)}"
    )


if __name__ == "__main__":
    main()
