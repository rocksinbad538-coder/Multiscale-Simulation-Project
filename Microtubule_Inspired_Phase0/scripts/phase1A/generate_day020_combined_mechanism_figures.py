#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day020_bright_combined_dephasing_relaxation"
)

OUTPUT_ROOT = (
    INPUT_ROOT
    / "combined_mechanism_figures"
)

CONDITION_CSV = (
    INPUT_ROOT
    / "condition_summary.csv"
)

TRAJECTORY_NPZ = (
    INPUT_ROOT
    / "combined_population_trajectories.npz"
)

TABLE_CSV = (
    OUTPUT_ROOT
    / "table_day020_combined_mechanism.csv"
)

MANIFEST_MD = (
    OUTPUT_ROOT
    / "FIGURE_MANIFEST_DAY020.md"
)

TEMPERATURES_K = (
    150.0,
    300.0,
)

KAPPA_VALUES = (
    0.1,
    1.0,
    10.0,
)

GAMMA_VALUES = (
    0.0,
    1.0,
    20.0,
    100.0,
)

SITE_LABELS = (
    "PYR2_bright",
    "PYR3_bright",
    "PYR4_bright",
    "PYR5_bright",
)

PYR5_INDEX = 3
HIGH_INITIAL_INDICES = (0, 1, 2)

REPRESENTATIVE_TEMPERATURE_K = 300.0
REPRESENTATIVE_KAPPA_PS_INV = 1.0

DPI = 300


def read_csv(
    path: Path,
) -> list[dict[str, str]]:
    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        return list(
            csv.DictReader(handle)
        )


def number(
    row: dict[str, str],
    key: str,
) -> float:
    return float(row[key])


def find_condition(
    rows: list[dict[str, str]],
    temperature_K: float,
    kappa_ps_inv: float,
    gamma_ps_inv: float,
) -> dict[str, str]:
    matches = [
        row
        for row in rows
        if np.isclose(
            number(
                row,
                "temperature_K",
            ),
            temperature_K,
        )
        and np.isclose(
            number(
                row,
                "kappa_ref_ps_inv",
            ),
            kappa_ps_inv,
        )
        and np.isclose(
            number(
                row,
                "gamma_phi_ps_inv",
            ),
            gamma_ps_inv,
        )
    ]

    if len(matches) != 1:
        raise RuntimeError(
            "Expected exactly one condition for "
            f"T={temperature_K}, "
            f"kappa={kappa_ps_inv}, "
            f"gamma={gamma_ps_inv}; "
            f"found {len(matches)}"
        )

    return matches[0]


def save_figure(
    figure: plt.Figure,
    stem: str,
) -> None:
    png_path = OUTPUT_ROOT / f"{stem}.png"
    pdf_path = OUTPUT_ROOT / f"{stem}.pdf"

    figure.savefig(
        png_path,
        dpi=DPI,
        bbox_inches="tight",
    )

    figure.savefig(
        pdf_path,
        bbox_inches="tight",
    )

    plt.close(figure)


def plot_final_pyr5_vs_gamma(
    rows: list[dict[str, str]],
) -> None:
    figure, axis = plt.subplots(
        figsize=(7.4, 5.2)
    )

    markers = (
        "o",
        "s",
        "^",
    )

    line_styles = {
        150.0: "-",
        300.0: "--",
    }

    gamma_positions = np.arange(
        len(GAMMA_VALUES)
    )

    for temperature_K in TEMPERATURES_K:
        for (
            kappa_index,
            kappa,
        ) in enumerate(KAPPA_VALUES):
            values = []

            for gamma in GAMMA_VALUES:
                row = find_condition(
                    rows,
                    temperature_K,
                    kappa,
                    gamma,
                )

                values.append(
                    number(
                        row,
                        "mean_final_PYR5_population_from_high_states",
                    )
                )

            axis.plot(
                gamma_positions,
                values,
                marker=markers[
                    kappa_index
                ],
                linestyle=line_styles[
                    temperature_K
                ],
                linewidth=1.6,
                markersize=5.5,
                label=(
                    f"{temperature_K:.0f} K, "
                    rf"$\kappa={kappa:g}$ ps$^{{-1}}$"
                ),
            )

    axis.set_xticks(
        gamma_positions,
        [
            f"{gamma:g}"
            for gamma in GAMMA_VALUES
        ],
    )

    axis.set_xlabel(
        r"Pure-dephasing rate "
        r"$\gamma_\phi$ (ps$^{-1}$)"
    )

    axis.set_ylabel(
        r"Mean PYR5 population at 100 ps"
    )

    axis.set_title(
        "Transient access to the low-energy PYR5 state"
    )

    axis.set_ylim(
        bottom=0.0
    )

    axis.grid(
        True,
        alpha=0.25,
    )

    axis.legend(
        fontsize=8,
        frameon=False,
        ncol=2,
    )

    figure.tight_layout()

    save_figure(
        figure,
        "fig_day020_final_PYR5_vs_dephasing",
    )


def plot_steady_pyr5_vs_gamma(
    rows: list[dict[str, str]],
) -> None:
    figure, axis = plt.subplots(
        figsize=(7.4, 5.2)
    )

    markers = (
        "o",
        "s",
        "^",
    )

    line_styles = {
        150.0: "-",
        300.0: "--",
    }

    gamma_positions = np.arange(
        len(GAMMA_VALUES)
    )

    for temperature_K in TEMPERATURES_K:
        for (
            kappa_index,
            kappa,
        ) in enumerate(KAPPA_VALUES):
            values = []

            for gamma in GAMMA_VALUES:
                row = find_condition(
                    rows,
                    temperature_K,
                    kappa,
                    gamma,
                )

                values.append(
                    number(
                        row,
                        "mean_steady_PYR5_population",
                    )
                )

            axis.plot(
                gamma_positions,
                values,
                marker=markers[
                    kappa_index
                ],
                linestyle=line_styles[
                    temperature_K
                ],
                linewidth=1.6,
                markersize=5.5,
                label=(
                    f"{temperature_K:.0f} K, "
                    rf"$\kappa={kappa:g}$ ps$^{{-1}}$"
                ),
            )

    axis.axhline(
        0.25,
        linestyle=":",
        linewidth=1.4,
        label="Uniform four-state limit",
    )

    axis.set_xticks(
        gamma_positions,
        [
            f"{gamma:g}"
            for gamma in GAMMA_VALUES
        ],
    )

    axis.set_xlabel(
        r"Pure-dephasing rate "
        r"$\gamma_\phi$ (ps$^{-1}$)"
    )

    axis.set_ylabel(
        r"Mean stationary PYR5 population"
    )

    axis.set_title(
        "Competition between thermal relaxation and local dephasing"
    )

    axis.set_ylim(
        0.2,
        1.03,
    )

    axis.grid(
        True,
        alpha=0.25,
    )

    axis.legend(
        fontsize=8,
        frameon=False,
        ncol=2,
    )

    figure.tight_layout()

    save_figure(
        figure,
        "fig_day020_steady_PYR5_vs_dephasing",
    )


def plot_gateway_fraction(
    rows: list[dict[str, str]],
) -> None:
    figure, axis = plt.subplots(
        figsize=(7.4, 5.2)
    )

    markers = (
        "o",
        "s",
        "^",
    )

    line_styles = {
        150.0: "-",
        300.0: "--",
    }

    gamma_positions = np.arange(
        len(GAMMA_VALUES)
    )

    for temperature_K in TEMPERATURES_K:
        for (
            kappa_index,
            kappa,
        ) in enumerate(KAPPA_VALUES):
            values = []

            for gamma in GAMMA_VALUES:
                row = find_condition(
                    rows,
                    temperature_K,
                    kappa,
                    gamma,
                )

                values.append(
                    number(
                        row,
                        "PYR4_gateway_fraction",
                    )
                )

            axis.plot(
                gamma_positions,
                values,
                marker=markers[
                    kappa_index
                ],
                linestyle=line_styles[
                    temperature_K
                ],
                linewidth=1.6,
                markersize=5.5,
                label=(
                    f"{temperature_K:.0f} K, "
                    rf"$\kappa={kappa:g}$ ps$^{{-1}}$"
                ),
            )

    axis.set_xticks(
        gamma_positions,
        [
            f"{gamma:g}"
            for gamma in GAMMA_VALUES
        ],
    )

    axis.set_xlabel(
        r"Pure-dephasing rate "
        r"$\gamma_\phi$ (ps$^{-1}$)"
    )

    axis.set_ylabel(
        "PYR4 contribution to downward PYR5 flux"
    )

    axis.set_title(
        "Robustness of PYR4 as the dominant gateway"
    )

    axis.set_ylim(
        0.82,
        0.93,
    )

    axis.grid(
        True,
        alpha=0.25,
    )

    axis.legend(
        fontsize=8,
        frameon=False,
        ncol=2,
    )

    figure.tight_layout()

    save_figure(
        figure,
        "fig_day020_PYR4_gateway_vs_dephasing",
    )


def locate_index(
    values: np.ndarray,
    target: float,
    label: str,
) -> int:
    indices = np.where(
        np.isclose(
            values,
            target,
        )
    )[0]

    if indices.size != 1:
        raise RuntimeError(
            f"Could not uniquely locate "
            f"{label}={target}"
        )

    return int(indices[0])


def plot_representative_time_traces() -> None:
    with np.load(
        TRAJECTORY_NPZ,
        allow_pickle=False,
    ) as data:
        temperatures = data[
            "temperatures_K"
        ].astype(np.float64)

        kappas = data[
            "kappa_ref_ps_inv"
        ].astype(np.float64)

        gammas = data[
            "gamma_phi_ps_inv"
        ].astype(np.float64)

        times_ps = data[
            "times_ps"
        ].astype(np.float64)

        populations = data[
            "ensemble_mean_populations"
        ].astype(np.float64)

        labels = tuple(
            str(value)
            for value in data[
                "site_labels"
            ]
        )

    if labels != SITE_LABELS:
        raise RuntimeError(
            f"Unexpected site labels: {labels}"
        )

    expected_shape = (
        temperatures.size,
        kappas.size,
        gammas.size,
        len(SITE_LABELS),
        times_ps.size,
        len(SITE_LABELS),
    )

    if populations.shape != expected_shape:
        raise RuntimeError(
            f"Unexpected population-array shape: "
            f"{populations.shape}; "
            f"expected {expected_shape}"
        )

    temperature_index = locate_index(
        temperatures,
        REPRESENTATIVE_TEMPERATURE_K,
        "temperature",
    )

    kappa_index = locate_index(
        kappas,
        REPRESENTATIVE_KAPPA_PS_INV,
        "kappa",
    )

    figure, axis = plt.subplots(
        figsize=(7.4, 5.2)
    )

    for gamma in GAMMA_VALUES:
        gamma_index = locate_index(
            gammas,
            gamma,
            "gamma",
        )

        pyr5_trace = np.mean(
            populations[
                temperature_index,
                kappa_index,
                gamma_index,
                HIGH_INITIAL_INDICES,
                :,
                PYR5_INDEX,
            ],
            axis=0,
        )

        axis.plot(
            times_ps,
            pyr5_trace,
            linewidth=1.7,
            label=(
                rf"$\gamma_\phi={gamma:g}$ "
                r"ps$^{-1}$"
            ),
        )

    axis.set_xlabel(
        "Time (ps)"
    )

    axis.set_ylabel(
        "Mean PYR5 population"
    )

    axis.set_title(
        (
            "Representative PYR5 population dynamics\n"
            rf"$T={REPRESENTATIVE_TEMPERATURE_K:.0f}$ K, "
            rf"$\kappa={REPRESENTATIVE_KAPPA_PS_INV:g}$ "
            r"ps$^{-1}$"
        )
    )

    axis.set_xlim(
        0.0,
        float(times_ps[-1]),
    )

    axis.set_ylim(
        bottom=0.0,
    )

    axis.grid(
        True,
        alpha=0.25,
    )

    axis.legend(
        frameon=False,
    )

    figure.tight_layout()

    save_figure(
        figure,
        "fig_day020_representative_PYR5_time_traces",
    )


def write_summary_table(
    rows: list[dict[str, str]],
) -> None:
    table_rows: list[
        dict[str, object]
    ] = []

    for row in rows:
        table_rows.append(
            {
                "temperature_K": number(
                    row,
                    "temperature_K",
                ),
                "kappa_ps_inv": number(
                    row,
                    "kappa_ref_ps_inv",
                ),
                "gamma_phi_ps_inv": number(
                    row,
                    "gamma_phi_ps_inv",
                ),
                "PYR5_at_100ps": number(
                    row,
                    "mean_final_PYR5_population_from_high_states",
                ),
                "stationary_PYR5": number(
                    row,
                    "mean_steady_PYR5_population",
                ),
                "PYR4_gateway_fraction": number(
                    row,
                    "PYR4_gateway_fraction",
                ),
                "steady_L1_distance_to_Gibbs": number(
                    row,
                    "mean_steady_l1_distance_to_Gibbs",
                ),
            }
        )

    with TABLE_CSV.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                table_rows[0].keys()
            ),
        )

        writer.writeheader()
        writer.writerows(
            table_rows
        )


def write_manifest(
    rows: list[dict[str, str]],
) -> None:
    gateway_values = [
        number(
            row,
            "PYR4_gateway_fraction",
        )
        for row in rows
    ]

    final_values = [
        number(
            row,
            "mean_final_PYR5_population_from_high_states",
        )
        for row in rows
    ]

    steady_values = [
        number(
            row,
            "mean_steady_PYR5_population",
        )
        for row in rows
    ]

    gamma100_values = [
        number(
            row,
            "mean_final_PYR5_population_from_high_states",
        )
        for row in rows
        if np.isclose(
            number(
                row,
                "gamma_phi_ps_inv",
            ),
            100.0,
        )
    ]

    with MANIFEST_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Combined-Mechanism Figures\n\n"
        )

        handle.write(
            "## Accepted quantitative ranges\n\n"
        )

        handle.write(
            f"- PYR4 gateway fraction: "
            f"{min(gateway_values):.6f} to "
            f"{max(gateway_values):.6f}.\n"
        )

        handle.write(
            f"- Final PYR5 population at 100 ps: "
            f"{min(final_values):.6f} to "
            f"{max(final_values):.6f}.\n"
        )

        handle.write(
            f"- Final PYR5 population for "
            f"gamma_phi=100 ps^-1: "
            f"{min(gamma100_values):.6f} to "
            f"{max(gamma100_values):.6f}.\n"
        )

        handle.write(
            f"- Stationary PYR5 population: "
            f"{min(steady_values):.6f} to "
            f"{max(steady_values):.6f}.\n\n"
        )

        handle.write(
            "## Figure files\n\n"
        )

        handle.write(
            "1. `fig_day020_final_PYR5_vs_dephasing`: "
            "transient PYR5 population at 100 ps.\n"
        )

        handle.write(
            "2. `fig_day020_steady_PYR5_vs_dephasing`: "
            "competition between thermal relaxation and "
            "local dephasing in the stationary state.\n"
        )

        handle.write(
            "3. `fig_day020_PYR4_gateway_vs_dephasing`: "
            "robustness of PYR4 as the dominant gateway.\n"
        )

        handle.write(
            "4. `fig_day020_representative_PYR5_time_traces`: "
            "representative transient dynamics for "
            "T=300 K and kappa=1 ps^-1.\n\n"
        )

        handle.write(
            "Each figure is written in PNG and PDF formats.\n\n"
        )

        handle.write(
            "## Interpretation boundary\n\n"
        )

        handle.write(
            "The absolute dissipative rates are "
            "phenomenological sensitivity parameters. "
            "The figures support relative mechanistic "
            "conclusions but not microscopic relaxation "
            "times derived from a spectral density.\n"
        )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    rows = read_csv(
        CONDITION_CSV
    )

    expected_conditions = (
        len(TEMPERATURES_K)
        * len(KAPPA_VALUES)
        * len(GAMMA_VALUES)
    )

    if len(rows) != expected_conditions:
        raise RuntimeError(
            f"Expected {expected_conditions} "
            f"conditions, found {len(rows)}"
        )

    plot_final_pyr5_vs_gamma(
        rows
    )

    plot_steady_pyr5_vs_gamma(
        rows
    )

    plot_gateway_fraction(
        rows
    )

    plot_representative_time_traces()

    write_summary_table(
        rows
    )

    write_manifest(
        rows
    )

    expected_stems = (
        "fig_day020_final_PYR5_vs_dephasing",
        "fig_day020_steady_PYR5_vs_dephasing",
        "fig_day020_PYR4_gateway_vs_dephasing",
        "fig_day020_representative_PYR5_time_traces",
    )

    missing_files = []

    for stem in expected_stems:
        for suffix in (
            ".png",
            ".pdf",
        ):
            path = (
                OUTPUT_ROOT
                / f"{stem}{suffix}"
            )

            if (
                not path.exists()
                or path.stat().st_size == 0
            ):
                missing_files.append(
                    str(path)
                )

    if missing_files:
        raise RuntimeError(
            "Missing or empty figure files:\n"
            + "\n".join(
                missing_files
            )
        )

    print(
        "Day020 combined-mechanism "
        "figures completed."
    )

    print(
        f"Conditions: "
        f"{len(rows)}/{expected_conditions}"
    )

    print(
        "Figures: 4 PNG + 4 PDF"
    )

    print(
        f"Table: "
        f"{TABLE_CSV.relative_to(PROJECT_ROOT)}"
    )

    print(
        f"Manifest: "
        f"{MANIFEST_MD.relative_to(PROJECT_ROOT)}"
    )

    print(
        f"Wrote: "
        f"{OUTPUT_ROOT.relative_to(PROJECT_ROOT)}"
    )


if __name__ == "__main__":
    main()
