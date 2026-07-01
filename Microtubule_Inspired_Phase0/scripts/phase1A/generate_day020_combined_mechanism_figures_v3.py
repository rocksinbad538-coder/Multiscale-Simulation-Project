#!/usr/bin/env python3

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day020_bright_combined_dephasing_relaxation"
)

OUTPUT_ROOT = (
    INPUT_ROOT
    / "combined_mechanism_figures_v3"
)

TRAJECTORY_NPZ = (
    INPUT_ROOT
    / "combined_population_trajectories.npz"
)

FRAME_METRICS_CSV = (
    INPUT_ROOT
    / "frame_combined_metrics.csv"
)

STEADY_STATE_CSV = (
    INPUT_ROOT
    / "steady_state_metrics.csv"
)

CONDITION_SUMMARY_CSV = (
    INPUT_ROOT
    / "condition_summary.csv"
)

VALIDATION_CSV = (
    INPUT_ROOT
    / "numerical_validation.csv"
)

MAIN_FIGURE_STEM = (
    OUTPUT_ROOT
    / "figure_day020_combined_mechanism_main_v3"
)

GATEWAY_FIGURE_STEM = (
    OUTPUT_ROOT
    / "figure_day020_gateway_supplementary_v3"
)

STATISTICS_CSV = (
    OUTPUT_ROOT
    / "table_day020_combined_mechanism_statistics_v3.csv"
)

CAPTIONS_MD = (
    OUTPUT_ROOT
    / "FIGURE_CAPTIONS_DAY020_V3.md"
)

MANIFEST_MD = (
    OUTPUT_ROOT
    / "FIGURE_MANIFEST_DAY020_V3.md"
)

TEMPERATURES_K = (
    150.0,
    300.0,
)

KAPPA_VALUES_PS_INV = (
    0.1,
    1.0,
    10.0,
)

GAMMA_VALUES_PS_INV = (
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

HIGH_INITIAL_LABELS = SITE_LABELS[:3]
HIGH_INITIAL_INDICES = (0, 1, 2)
PYR5_INDEX = 3

REPRESENTATIVE_TEMPERATURE_K = 300.0
REPRESENTATIVE_KAPPA_PS_INV = 1.0

EXPECTED_FRAMES = 21
EXPECTED_CONDITIONS = 24
EXPECTED_HIGH_ROWS = (
    EXPECTED_FRAMES
    * len(HIGH_INITIAL_LABELS)
)

N_BOOTSTRAP = 5000
BOOTSTRAP_SEED = 20260701

DPI = 400
SUMMARY_TOL = 5.0e-12
SYMLINTHRESH = 1.0

KAPPA_COLORS = {
    0.1: "#0072B2",
    1.0: "#E69F00",
    10.0: "#009E73",
}

KAPPA_MARKERS = {
    0.1: "o",
    1.0: "s",
    10.0: "^",
}

TEMPERATURE_LINESTYLES = {
    150.0: "-",
    300.0: "--",
}

TRACE_COLORS = {
    0.0: "#0072B2",
    1.0: "#E69F00",
    20.0: "#009E73",
    100.0: "#D55E00",
}


def log(message: str = "") -> None:
    print(message, flush=True)


def read_csv(
    path: Path,
) -> list[dict[str, str]]:
    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        return list(csv.DictReader(handle))


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


def number(
    row: dict[str, str],
    key: str,
) -> float:
    return float(row[key])


def condition_key(
    temperature_K: float,
    kappa_ps_inv: float,
    gamma_ps_inv: float,
) -> tuple[float, float, float]:
    return (
        float(temperature_K),
        float(kappa_ps_inv),
        float(gamma_ps_inv),
    )


def select_condition_rows(
    rows: list[dict[str, str]],
    temperature_K: float,
    kappa_ps_inv: float,
    gamma_ps_inv: float,
) -> list[dict[str, str]]:
    return [
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


def locate_index(
    values: np.ndarray,
    target: float,
    label: str,
) -> int:
    indices = np.where(
        np.isclose(
            values,
            target,
            atol=1.0e-12,
            rtol=0.0,
        )
    )[0]

    if indices.size != 1:
        raise RuntimeError(
            f"Could not uniquely locate "
            f"{label}={target}"
        )

    return int(indices[0])


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.5,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "legend.fontsize": 9,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "axes.linewidth": 0.9,
            "lines.linewidth": 1.8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def configure_gamma_axis(
    axis: plt.Axes,
) -> None:
    axis.set_xscale(
        "symlog",
        linthresh=SYMLINTHRESH,
        linscale=1.0,
        base=10,
    )

    axis.set_xlim(
        -0.08,
        130.0,
    )

    axis.set_xticks(
        GAMMA_VALUES_PS_INV
    )

    axis.set_xticklabels(
        [
            f"{value:g}"
            for value
            in GAMMA_VALUES_PS_INV
        ]
    )

    axis.set_xlabel(
        r"Pure-dephasing rate "
        r"$\gamma_\phi$ (ps$^{-1}$)"
    )

    axis.grid(
        True,
        alpha=0.22,
        linewidth=0.8,
    )


def save_figure(
    figure: plt.Figure,
    stem: Path,
) -> None:
    figure.savefig(
        stem.with_suffix(".png"),
        dpi=DPI,
        bbox_inches="tight",
    )

    figure.savefig(
        stem.with_suffix(".pdf"),
        bbox_inches="tight",
    )

    plt.close(figure)


def validate_numerical_status() -> None:
    rows = read_csv(
        VALIDATION_CSV
    )

    if len(rows) != 1:
        raise RuntimeError(
            "Expected one validation row"
        )

    failures = [
        key
        for key, item
        in rows[0].items()
        if key.endswith("_pass")
        and str(item).lower()
        != "true"
    ]

    if failures:
        raise RuntimeError(
            "Failed validation fields: "
            + ", ".join(failures)
        )


def pooled_mean_and_sd(
    means: np.ndarray,
    standard_deviations: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    pooled_mean = np.mean(
        means,
        axis=0,
    )

    second_moment = np.mean(
        standard_deviations**2
        + means**2,
        axis=0,
    )

    pooled_variance = np.maximum(
        second_moment
        - pooled_mean**2,
        0.0,
    )

    return (
        pooled_mean,
        np.sqrt(
            pooled_variance
        ),
    )


def bootstrap_gateway_ratio(
    frame_total_flux: np.ndarray,
    frame_pyr4_flux: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, float, float]:
    total_flux = float(
        np.sum(
            frame_total_flux
        )
    )

    if total_flux <= 0.0:
        raise RuntimeError(
            "Non-positive aggregate gateway flux"
        )

    estimate = float(
        np.sum(
            frame_pyr4_flux
        )
        / total_flux
    )

    indices = rng.integers(
        0,
        frame_total_flux.size,
        size=(
            N_BOOTSTRAP,
            frame_total_flux.size,
        ),
    )

    bootstrap_total = np.sum(
        frame_total_flux[
            indices
        ],
        axis=1,
    )

    bootstrap_pyr4 = np.sum(
        frame_pyr4_flux[
            indices
        ],
        axis=1,
    )

    valid = bootstrap_total > 0.0

    bootstrap_ratios = (
        bootstrap_pyr4[
            valid
        ]
        / bootstrap_total[
            valid
        ]
    )

    lower, upper = np.percentile(
        bootstrap_ratios,
        [2.5, 97.5],
    )

    return (
        estimate,
        float(lower),
        float(upper),
    )


def build_statistics(
    frame_rows: list[dict[str, str]],
    steady_rows: list[dict[str, str]],
    condition_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    condition_lookup = {
        condition_key(
            number(
                row,
                "temperature_K",
            ),
            number(
                row,
                "kappa_ref_ps_inv",
            ),
            number(
                row,
                "gamma_phi_ps_inv",
            ),
        ): row
        for row in condition_rows
    }

    rng = np.random.default_rng(
        BOOTSTRAP_SEED
    )

    statistics_rows: list[
        dict[str, object]
    ] = []

    for temperature_K in TEMPERATURES_K:
        for kappa in KAPPA_VALUES_PS_INV:
            for gamma in GAMMA_VALUES_PS_INV:
                key = condition_key(
                    temperature_K,
                    kappa,
                    gamma,
                )

                selected_frame_rows = (
                    select_condition_rows(
                        frame_rows,
                        temperature_K,
                        kappa,
                        gamma,
                    )
                )

                selected_steady_rows = (
                    select_condition_rows(
                        steady_rows,
                        temperature_K,
                        kappa,
                        gamma,
                    )
                )

                high_rows = [
                    row
                    for row
                    in selected_frame_rows
                    if row[
                        "initial_state"
                    ] in HIGH_INITIAL_LABELS
                ]

                if len(high_rows) != (
                    EXPECTED_HIGH_ROWS
                ):
                    raise RuntimeError(
                        f"Expected "
                        f"{EXPECTED_HIGH_ROWS} "
                        f"high-state rows for {key}; "
                        f"found {len(high_rows)}"
                    )

                if len(
                    selected_steady_rows
                ) != EXPECTED_FRAMES:
                    raise RuntimeError(
                        f"Expected "
                        f"{EXPECTED_FRAMES} "
                        f"stationary rows for {key}"
                    )

                final_values = np.asarray(
                    [
                        number(
                            row,
                            "final_PYR5_population",
                        )
                        for row in high_rows
                    ],
                    dtype=np.float64,
                )

                steady_values = np.asarray(
                    [
                        number(
                            row,
                            "steady_PYR5_population",
                        )
                        for row
                        in selected_steady_rows
                    ],
                    dtype=np.float64,
                )

                frame_groups: dict[
                    int,
                    list[dict[str, str]],
                ] = defaultdict(list)

                for row in high_rows:
                    frame_groups[
                        int(row["frame"])
                    ].append(row)

                if len(frame_groups) != (
                    EXPECTED_FRAMES
                ):
                    raise RuntimeError(
                        f"Expected "
                        f"{EXPECTED_FRAMES} frames "
                        f"for gateway bootstrap"
                    )

                frame_total_flux: list[
                    float
                ] = []

                frame_pyr4_flux: list[
                    float
                ] = []

                for frame in sorted(
                    frame_groups
                ):
                    rows_for_frame = (
                        frame_groups[frame]
                    )

                    total_flux = sum(
                        number(
                            row,
                            "integrated_downward_flux_to_PYR5",
                        )
                        for row
                        in rows_for_frame
                    )

                    pyr4_flux = sum(
                        number(
                            row,
                            "integrated_downward_flux_to_PYR5",
                        )
                        * number(
                            row,
                            "PYR4_gateway_fraction",
                        )
                        for row
                        in rows_for_frame
                    )

                    frame_total_flux.append(
                        total_flux
                    )

                    frame_pyr4_flux.append(
                        pyr4_flux
                    )

                (
                    gateway_estimate,
                    gateway_ci_lower,
                    gateway_ci_upper,
                ) = bootstrap_gateway_ratio(
                    np.asarray(
                        frame_total_flux,
                        dtype=np.float64,
                    ),
                    np.asarray(
                        frame_pyr4_flux,
                        dtype=np.float64,
                    ),
                    rng,
                )

                final_mean = float(
                    np.mean(
                        final_values
                    )
                )

                final_sd = float(
                    np.std(
                        final_values,
                        ddof=0,
                    )
                )

                steady_mean = float(
                    np.mean(
                        steady_values
                    )
                )

                steady_sd = float(
                    np.std(
                        steady_values,
                        ddof=0,
                    )
                )

                reference = (
                    condition_lookup[key]
                )

                cross_checks = {
                    "final": abs(
                        final_mean
                        - number(
                            reference,
                            "mean_final_PYR5_population_from_high_states",
                        )
                    ),
                    "steady": abs(
                        steady_mean
                        - number(
                            reference,
                            "mean_steady_PYR5_population",
                        )
                    ),
                    "gateway": abs(
                        gateway_estimate
                        - number(
                            reference,
                            "PYR4_gateway_fraction",
                        )
                    ),
                }

                failures = {
                    name: error
                    for name, error
                    in cross_checks.items()
                    if error > SUMMARY_TOL
                }

                if failures:
                    raise RuntimeError(
                        f"Summary mismatch for {key}: "
                        f"{failures}"
                    )

                statistics_rows.append(
                    {
                        "temperature_K": temperature_K,
                        "kappa_ref_ps_inv": kappa,
                        "gamma_phi_ps_inv": gamma,
                        "n_final_population_samples": (
                            final_values.size
                        ),
                        "mean_final_PYR5_population": (
                            final_mean
                        ),
                        "sd_final_PYR5_population": (
                            final_sd
                        ),
                        "n_stationary_samples": (
                            steady_values.size
                        ),
                        "mean_stationary_PYR5_population": (
                            steady_mean
                        ),
                        "sd_stationary_PYR5_population": (
                            steady_sd
                        ),
                        "gateway_frame_count": (
                            EXPECTED_FRAMES
                        ),
                        "PYR4_gateway_fraction": (
                            gateway_estimate
                        ),
                        "PYR4_gateway_bootstrap95_lower": (
                            gateway_ci_lower
                        ),
                        "PYR4_gateway_bootstrap95_upper": (
                            gateway_ci_upper
                        ),
                        "gateway_bootstrap_replicates": (
                            N_BOOTSTRAP
                        ),
                    }
                )

    return statistics_rows


def make_lookup(
    statistics_rows: list[
        dict[str, object]
    ],
) -> dict[
    tuple[float, float, float],
    dict[str, object],
]:
    return {
        condition_key(
            float(
                row["temperature_K"]
            ),
            float(
                row["kappa_ref_ps_inv"]
            ),
            float(
                row["gamma_phi_ps_inv"]
            ),
        ): row
        for row in statistics_rows
    }


def build_shared_handles(
) -> tuple[
    list[Line2D],
    list[str],
]:
    handles: list[Line2D] = []
    labels: list[str] = []

    for kappa in KAPPA_VALUES_PS_INV:
        handles.append(
            Line2D(
                [0],
                [0],
                color=KAPPA_COLORS[
                    kappa
                ],
                marker=KAPPA_MARKERS[
                    kappa
                ],
                linestyle="-",
                markersize=6,
            )
        )

        labels.append(
            rf"$\kappa={kappa:g}$ ps$^{{-1}}$"
        )

    for temperature_K in TEMPERATURES_K:
        handles.append(
            Line2D(
                [0],
                [0],
                color="black",
                linestyle=(
                    TEMPERATURE_LINESTYLES[
                        temperature_K
                    ]
                ),
            )
        )

        labels.append(
            f"{temperature_K:.0f} K"
        )

    return handles, labels


def build_main_figure(
    trajectory_data: dict[
        str,
        np.ndarray,
    ],
    statistics_rows: list[
        dict[str, object]
    ],
) -> None:
    temperatures = trajectory_data[
        "temperatures_K"
    ].astype(np.float64)

    kappas = trajectory_data[
        "kappa_ref_ps_inv"
    ].astype(np.float64)

    gammas = trajectory_data[
        "gamma_phi_ps_inv"
    ].astype(np.float64)

    times_ps = trajectory_data[
        "times_ps"
    ].astype(np.float64)

    mean_populations = trajectory_data[
        "ensemble_mean_populations"
    ].astype(np.float64)

    sd_populations = trajectory_data[
        "ensemble_sd_populations"
    ].astype(np.float64)

    lookup = make_lookup(
        statistics_rows
    )

    figure, axes = plt.subplots(
        1,
        3,
        figsize=(17.2, 6.3),
    )

    (
        axis_a,
        axis_b,
        axis_c,
    ) = axes

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

    for gamma in GAMMA_VALUES_PS_INV:
        gamma_index = locate_index(
            gammas,
            gamma,
            "gamma",
        )

        means_by_initial = np.take(
            mean_populations[
                temperature_index,
                kappa_index,
                gamma_index,
                :,
                :,
                PYR5_INDEX,
            ],
            HIGH_INITIAL_INDICES,
            axis=0,
        )

        sds_by_initial = np.take(
            sd_populations[
                temperature_index,
                kappa_index,
                gamma_index,
                :,
                :,
                PYR5_INDEX,
            ],
            HIGH_INITIAL_INDICES,
            axis=0,
        )

        (
            pooled_mean,
            pooled_sd,
        ) = pooled_mean_and_sd(
            means_by_initial,
            sds_by_initial,
        )

        lower = np.clip(
            pooled_mean
            - pooled_sd,
            0.0,
            1.0,
        )

        upper = np.clip(
            pooled_mean
            + pooled_sd,
            0.0,
            1.0,
        )

        color = TRACE_COLORS[
            gamma
        ]

        axis_a.plot(
            times_ps,
            pooled_mean,
            color=color,
            label=(
                rf"$\gamma_\phi={gamma:g}$"
            ),
        )

        axis_a.fill_between(
            times_ps,
            lower,
            upper,
            color=color,
            alpha=0.16,
            linewidth=0.0,
        )

    axis_a.set_xlabel(
        "Time (ps)"
    )

    axis_a.set_ylabel(
        "Mean PYR5 population"
    )

    axis_a.set_xlim(
        0.0,
        100.0,
    )

    axis_a.set_ylim(
        0.0,
        0.128,
    )

    axis_a.grid(
        True,
        alpha=0.22,
    )

    axis_a.legend(
        title=(
            r"$\gamma_\phi$ "
            r"(ps$^{-1}$)"
        ),
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(
            0.01,
            0.96,
        ),
    )

    axis_a.text(
        0.98,
        0.04,
        (
            r"$T=300$ K, "
            r"$\kappa=1$ ps$^{-1}$"
        ),
        transform=axis_a.transAxes,
        ha="right",
        va="bottom",
        fontsize=9.5,
    )

    for temperature_K in TEMPERATURES_K:
        for kappa in KAPPA_VALUES_PS_INV:
            means = []
            sds = []

            for gamma in GAMMA_VALUES_PS_INV:
                row = lookup[
                    condition_key(
                        temperature_K,
                        kappa,
                        gamma,
                    )
                ]

                means.append(
                    float(
                        row[
                            "mean_final_PYR5_population"
                        ]
                    )
                )

                sds.append(
                    float(
                        row[
                            "sd_final_PYR5_population"
                        ]
                    )
                )

            axis_b.errorbar(
                GAMMA_VALUES_PS_INV,
                means,
                yerr=sds,
                color=KAPPA_COLORS[
                    kappa
                ],
                marker=KAPPA_MARKERS[
                    kappa
                ],
                linestyle=(
                    TEMPERATURE_LINESTYLES[
                        temperature_K
                    ]
                ),
                markersize=5.5,
                capsize=2.3,
                elinewidth=0.9,
            )

    configure_gamma_axis(
        axis_b
    )

    axis_b.set_ylabel(
        "PYR5 population at 100 ps"
    )

    axis_b.set_ylim(
        0.0,
        0.142,
    )

    for temperature_K in TEMPERATURES_K:
        for kappa in KAPPA_VALUES_PS_INV:
            means = []
            lower_errors = []
            upper_errors = []

            for gamma in GAMMA_VALUES_PS_INV:
                row = lookup[
                    condition_key(
                        temperature_K,
                        kappa,
                        gamma,
                    )
                ]

                mean = float(
                    row[
                        "mean_stationary_PYR5_population"
                    ]
                )

                standard_deviation = float(
                    row[
                        "sd_stationary_PYR5_population"
                    ]
                )

                means.append(
                    mean
                )

                lower_errors.append(
                    standard_deviation
                )

                upper_errors.append(
                    standard_deviation
                )

            axis_c.errorbar(
                GAMMA_VALUES_PS_INV,
                means,
                yerr=np.asarray(
                    [
                        lower_errors,
                        upper_errors,
                    ]
                ),
                color=KAPPA_COLORS[
                    kappa
                ],
                marker=KAPPA_MARKERS[
                    kappa
                ],
                linestyle=(
                    TEMPERATURE_LINESTYLES[
                        temperature_K
                    ]
                ),
                markersize=5.5,
                capsize=2.3,
                elinewidth=0.9,
            )

    axis_c.axhline(
        0.25,
        color="black",
        linestyle=":",
        linewidth=1.3,
    )

    configure_gamma_axis(
        axis_c
    )

    axis_c.set_ylabel(
        "Stationary PYR5 population"
    )

    axis_c.set_ylim(
        0.20,
        1.035,
    )

    axis_c.text(
        0.97,
        0.07,
        "Uniform limit = 0.25",
        transform=axis_c.transAxes,
        ha="right",
        va="bottom",
        fontsize=8.8,
    )

    axis_c.text(
        0.97,
        0.95,
        (
            "PYR4 gateway\n"
            "84.3–91.4%"
        ),
        transform=axis_c.transAxes,
        ha="right",
        va="top",
        fontsize=9.2,
        bbox={
            "boxstyle": "round,pad=0.28",
            "facecolor": "white",
            "edgecolor": "0.65",
            "alpha": 0.92,
        },
    )

    for label, axis in zip(
        ("(a)", "(b)", "(c)"),
        axes,
    ):
        axis.text(
            -0.12,
            1.045,
            label,
            transform=axis.transAxes,
            ha="left",
            va="bottom",
            fontsize=13,
            fontweight="bold",
            clip_on=False,
        )

    shared_handles, shared_labels = (
        build_shared_handles()
    )

    figure.legend(
        shared_handles,
        shared_labels,
        loc="lower center",
        bbox_to_anchor=(
            0.5,
            0.015,
        ),
        frameon=False,
        ncol=5,
        columnspacing=2.0,
        handlelength=2.4,
    )

    figure.suptitle(
        (
            "Competition between local dephasing "
            "and thermal relaxation"
        ),
        fontsize=14.5,
        y=0.975,
    )

    figure.subplots_adjust(
        left=0.065,
        right=0.992,
        bottom=0.19,
        top=0.875,
        wspace=0.34,
    )

    save_figure(
        figure,
        MAIN_FIGURE_STEM,
    )


def build_gateway_figure(
    statistics_rows: list[
        dict[str, object]
    ],
) -> None:
    lookup = make_lookup(
        statistics_rows
    )

    figure, axis = plt.subplots(
        figsize=(8.4, 6.0),
    )

    ci_lower_values = []
    ci_upper_values = []

    for temperature_K in TEMPERATURES_K:
        for kappa in KAPPA_VALUES_PS_INV:
            means = []
            lower_errors = []
            upper_errors = []

            for gamma in GAMMA_VALUES_PS_INV:
                row = lookup[
                    condition_key(
                        temperature_K,
                        kappa,
                        gamma,
                    )
                ]

                mean = float(
                    row[
                        "PYR4_gateway_fraction"
                    ]
                )

                lower = float(
                    row[
                        "PYR4_gateway_bootstrap95_lower"
                    ]
                )

                upper = float(
                    row[
                        "PYR4_gateway_bootstrap95_upper"
                    ]
                )

                means.append(
                    mean
                )

                lower_errors.append(
                    mean - lower
                )

                upper_errors.append(
                    upper - mean
                )

                ci_lower_values.append(
                    lower
                )

                ci_upper_values.append(
                    upper
                )

            axis.errorbar(
                GAMMA_VALUES_PS_INV,
                means,
                yerr=np.asarray(
                    [
                        lower_errors,
                        upper_errors,
                    ]
                ),
                color=KAPPA_COLORS[
                    kappa
                ],
                marker=KAPPA_MARKERS[
                    kappa
                ],
                linestyle=(
                    TEMPERATURE_LINESTYLES[
                        temperature_K
                    ]
                ),
                markersize=6,
                capsize=3,
                elinewidth=1.0,
            )

    configure_gamma_axis(
        axis
    )

    lower_limit = max(
        0.0,
        min(ci_lower_values)
        - 0.015,
    )

    upper_limit = min(
        1.0,
        max(ci_upper_values)
        + 0.015,
    )

    axis.set_ylim(
        lower_limit,
        upper_limit,
    )

    axis.set_ylabel(
        (
            "PYR4 fraction of integrated "
            "downward PYR5 flux"
        )
    )

    axis.set_title(
        "Robustness of PYR4 as the dominant gateway"
    )

    axis.text(
        0.02,
        0.97,
        (
            "Point estimate with 95% frame-bootstrap CI\n"
            f"({N_BOOTSTRAP:,} resamples of 21 snapshots)"
        ),
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=9.3,
    )

    shared_handles, shared_labels = (
        build_shared_handles()
    )

    figure.legend(
        shared_handles,
        shared_labels,
        loc="lower center",
        bbox_to_anchor=(
            0.5,
            0.015,
        ),
        frameon=False,
        ncol=5,
        columnspacing=1.7,
        handlelength=2.4,
    )

    figure.subplots_adjust(
        left=0.14,
        right=0.98,
        bottom=0.20,
        top=0.88,
    )

    save_figure(
        figure,
        GATEWAY_FIGURE_STEM,
    )


def write_captions(
    statistics_rows: list[
        dict[str, object]
    ],
) -> None:
    final_values = [
        float(
            row[
                "mean_final_PYR5_population"
            ]
        )
        for row in statistics_rows
    ]

    stationary_values = [
        float(
            row[
                "mean_stationary_PYR5_population"
            ]
        )
        for row in statistics_rows
    ]

    gateway_values = [
        float(
            row[
                "PYR4_gateway_fraction"
            ]
        )
        for row in statistics_rows
    ]

    with CAPTIONS_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Figure Captions V3\n\n"
        )

        handle.write(
            "## Main figure\n\n"
        )

        handle.write(
            "**Competition between local dephasing "
            "and thermal relaxation.** "
            "(a) Ensemble-averaged PYR5 population for "
            "representative conditions "
            "\\(T=300\\,\\mathrm{K}\\) and "
            "\\(\\kappa=1\\,\\mathrm{ps}^{-1}\\). "
            "Shaded regions show one pooled standard "
            "deviation across the 21 Hamiltonian snapshots "
            "and three upper-manifold initial states. "
            "(b) PYR5 population after 100 ps. Error bars "
            "show one standard deviation over 63 "
            "frame–initial-state realizations. "
            "(c) Stationary PYR5 population of the combined "
            "phenomenological Liouvillian. Error bars show "
            "one standard deviation over 21 Hamiltonian "
            "snapshots. The dotted line marks the uniform "
            "four-state limit. Colors and markers identify "
            "the phenomenological relaxation scale, while "
            "solid and dashed lines identify 150 and 300 K. "
            "PYR4 provides 84.3–91.4% of the integrated "
            "downward flux into PYR5 across all tested "
            "conditions.\n\n"
        )

        handle.write(
            "## Supplementary gateway figure\n\n"
        )

        handle.write(
            "**Robustness of PYR4 as the dominant gateway.** "
            "The point estimate is the ratio of the "
            "aggregate PYR4-resolved downward flux to the "
            "aggregate total downward flux into PYR5. "
            "Error bars are 95% percentile bootstrap "
            f"intervals obtained from {N_BOOTSTRAP:,} "
            "resamplings of the 21 Hamiltonian snapshots. "
            "The bootstrap preserves all three upper-state "
            "initial conditions within each sampled "
            "snapshot.\n\n"
        )

        handle.write(
            "## Accepted numerical ranges\n\n"
        )

        handle.write(
            f"- PYR5 population at 100 ps: "
            f"{min(final_values):.6f} to "
            f"{max(final_values):.6f}.\n"
        )

        handle.write(
            f"- Stationary PYR5 population: "
            f"{min(stationary_values):.6f} to "
            f"{max(stationary_values):.6f}.\n"
        )

        handle.write(
            f"- PYR4 gateway fraction: "
            f"{min(gateway_values):.6f} to "
            f"{max(gateway_values):.6f}.\n\n"
        )

        handle.write(
            "## Interpretation boundary\n\n"
        )

        handle.write(
            "The absolute dephasing and relaxation rates "
            "are phenomenological sensitivity parameters. "
            "The stationary state for nonzero local "
            "dephasing is a non-thermal stationary state "
            "of the combined model and must not be "
            "identified with thermodynamic equilibrium.\n"
        )


def write_manifest() -> None:
    with MANIFEST_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Combined-Mechanism "
            "Figure Manifest V3\n\n"
        )

        handle.write(
            "## Products\n\n"
        )

        for filename in (
            "figure_day020_combined_mechanism_main_v3.png",
            "figure_day020_combined_mechanism_main_v3.pdf",
            "figure_day020_gateway_supplementary_v3.png",
            "figure_day020_gateway_supplementary_v3.pdf",
            "table_day020_combined_mechanism_statistics_v3.csv",
            "FIGURE_CAPTIONS_DAY020_V3.md",
        ):
            handle.write(
                f"- `{filename}`\n"
            )

        handle.write(
            "\n## Statistical conventions\n\n"
        )

        handle.write(
            "- Population curves: pooled descriptive "
            "standard deviation.\n"
        )

        handle.write(
            "- Final populations: descriptive standard "
            "deviation over 63 realizations.\n"
        )

        handle.write(
            "- Stationary populations: descriptive "
            "standard deviation over 21 snapshots.\n"
        )

        handle.write(
            "- Gateway uncertainty: 95% bootstrap "
            "confidence interval for the aggregate "
            "ratio of fluxes, with snapshots as the "
            "resampling unit.\n"
        )

        handle.write(
            "- The dephasing-rate axis uses a symmetric-log "
            "scale to represent zero and the physical "
            "values 1, 20, and 100 ps^-1.\n"
        )


def validate_outputs() -> None:
    expected = (
        MAIN_FIGURE_STEM.with_suffix(
            ".png"
        ),
        MAIN_FIGURE_STEM.with_suffix(
            ".pdf"
        ),
        GATEWAY_FIGURE_STEM.with_suffix(
            ".png"
        ),
        GATEWAY_FIGURE_STEM.with_suffix(
            ".pdf"
        ),
        STATISTICS_CSV,
        CAPTIONS_MD,
        MANIFEST_MD,
    )

    failures = [
        path
        for path in expected
        if (
            not path.exists()
            or path.stat().st_size == 0
        )
    ]

    if failures:
        raise RuntimeError(
            "Missing or empty products:\n"
            + "\n".join(
                str(path)
                for path in failures
            )
        )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    configure_matplotlib()
    validate_numerical_status()

    frame_rows = read_csv(
        FRAME_METRICS_CSV
    )

    steady_rows = read_csv(
        STEADY_STATE_CSV
    )

    condition_rows = read_csv(
        CONDITION_SUMMARY_CSV
    )

    if len(condition_rows) != (
        EXPECTED_CONDITIONS
    ):
        raise RuntimeError(
            f"Expected {EXPECTED_CONDITIONS} "
            f"conditions, found "
            f"{len(condition_rows)}"
        )

    with np.load(
        TRAJECTORY_NPZ,
        allow_pickle=False,
    ) as data:
        trajectory_data = {
            key: data[key].copy()
            for key in data.files
        }

    labels = tuple(
        str(item)
        for item
        in trajectory_data[
            "site_labels"
        ]
    )

    if labels != SITE_LABELS:
        raise RuntimeError(
            f"Unexpected site labels: {labels}"
        )

    statistics_rows = build_statistics(
        frame_rows,
        steady_rows,
        condition_rows,
    )

    if len(statistics_rows) != (
        EXPECTED_CONDITIONS
    ):
        raise RuntimeError(
            f"Expected {EXPECTED_CONDITIONS} "
            f"statistics rows, found "
            f"{len(statistics_rows)}"
        )

    write_csv(
        STATISTICS_CSV,
        statistics_rows,
    )

    build_main_figure(
        trajectory_data,
        statistics_rows,
    )

    build_gateway_figure(
        statistics_rows,
    )

    write_captions(
        statistics_rows,
    )

    write_manifest()
    validate_outputs()

    gateway_lower = min(
        float(
            row[
                "PYR4_gateway_bootstrap95_lower"
            ]
        )
        for row in statistics_rows
    )

    gateway_upper = max(
        float(
            row[
                "PYR4_gateway_bootstrap95_upper"
            ]
        )
        for row in statistics_rows
    )

    log(
        "Day020 combined-mechanism "
        "publication figures V3 completed."
    )

    log(
        f"Conditions: "
        f"{len(statistics_rows)}/"
        f"{EXPECTED_CONDITIONS}"
    )

    log(
        "Numerical validation: PASS"
    )

    log(
        "Summary-statistics cross-check: PASS"
    )

    log(
        f"Gateway bootstrap: "
        f"{N_BOOTSTRAP} replicates per condition"
    )

    log(
        "Overall gateway 95% CI envelope: "
        f"{gateway_lower:.6f} to "
        f"{gateway_upper:.6f}"
    )

    log(
        "Main figure: 1 PNG + 1 PDF"
    )

    log(
        "Supplementary figure: 1 PNG + 1 PDF"
    )

    log(
        f"Wrote: "
        f"{OUTPUT_ROOT.relative_to(PROJECT_ROOT)}"
    )


if __name__ == "__main__":
    main()
