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
    / "combined_mechanism_figures_v2"
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
    / "figure_day020_combined_mechanism_main"
)

GATEWAY_FIGURE_STEM = (
    OUTPUT_ROOT
    / "figure_day020_gateway_supplementary"
)

STATISTICS_CSV = (
    OUTPUT_ROOT
    / "table_day020_combined_mechanism_statistics.csv"
)

CAPTIONS_MD = (
    OUTPUT_ROOT
    / "FIGURE_CAPTIONS_DAY020.md"
)

MANIFEST_MD = (
    OUTPUT_ROOT
    / "FIGURE_MANIFEST_DAY020_V2.md"
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

EXPECTED_CONDITIONS = 24
EXPECTED_FRAMES = 21
EXPECTED_HIGH_INITIAL_ROWS = (
    EXPECTED_FRAMES
    * len(HIGH_INITIAL_LABELS)
)

DPI = 400
SYMLINTHRESH = 1.0
SUMMARY_TOL = 5.0e-12

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


def value(
    row: dict[str, str],
    key: str,
) -> float:
    return float(row[key])


def locate_index(
    array: np.ndarray,
    target: float,
    label: str,
) -> int:
    matches = np.where(
        np.isclose(
            array,
            target,
            atol=1.0e-12,
            rtol=0.0,
        )
    )[0]

    if matches.size != 1:
        raise RuntimeError(
            f"Could not uniquely locate "
            f"{label}={target}"
        )

    return int(matches[0])


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


def select_rows(
    rows: list[dict[str, str]],
    temperature_K: float,
    kappa_ps_inv: float,
    gamma_ps_inv: float,
) -> list[dict[str, str]]:
    return [
        row
        for row in rows
        if np.isclose(
            value(
                row,
                "temperature_K",
            ),
            temperature_K,
        )
        and np.isclose(
            value(
                row,
                "kappa_ref_ps_inv",
            ),
            kappa_ps_inv,
        )
        and np.isclose(
            value(
                row,
                "gamma_phi_ps_inv",
            ),
            gamma_ps_inv,
        )
    ]


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
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
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
            f"{gamma:g}"
            for gamma
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


def pooled_mean_and_sd(
    means: np.ndarray,
    standard_deviations: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    if means.shape != standard_deviations.shape:
        raise RuntimeError(
            "Mean and SD arrays have "
            "different shapes"
        )

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
        np.sqrt(pooled_variance),
    )


def weighted_mean_and_sd(
    observations: np.ndarray,
    weights: np.ndarray,
) -> tuple[float, float]:
    positive = weights > 0.0

    observations = observations[
        positive
    ]

    weights = weights[
        positive
    ]

    if observations.size == 0:
        raise RuntimeError(
            "No positive-weight observations"
        )

    weight_sum = float(
        np.sum(weights)
    )

    mean = float(
        np.sum(
            weights * observations
        )
        / weight_sum
    )

    variance = float(
        np.sum(
            weights
            * (
                observations
                - mean
            )**2
        )
        / weight_sum
    )

    return (
        mean,
        float(
            np.sqrt(
                max(
                    variance,
                    0.0,
                )
            )
        ),
    )


def load_trajectory_data() -> dict[str, np.ndarray]:
    with np.load(
        TRAJECTORY_NPZ,
        allow_pickle=False,
    ) as data:
        result = {
            key: data[key].copy()
            for key in data.files
        }

    labels = tuple(
        str(item)
        for item
        in result["site_labels"]
    )

    if labels != SITE_LABELS:
        raise RuntimeError(
            f"Unexpected site labels: {labels}"
        )

    return result


def validate_numerical_status() -> None:
    rows = read_csv(
        VALIDATION_CSV
    )

    if len(rows) != 1:
        raise RuntimeError(
            "Expected one numerical-validation row"
        )

    failures = [
        key
        for key, item
        in rows[0].items()
        if key.endswith("_pass")
        and str(item).lower() != "true"
    ]

    if failures:
        raise RuntimeError(
            "Numerical validation failed: "
            + ", ".join(failures)
        )


def build_statistics(
    frame_rows: list[dict[str, str]],
    steady_rows: list[dict[str, str]],
    condition_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    condition_lookup: dict[
        tuple[float, float, float],
        dict[str, str],
    ] = {}

    for row in condition_rows:
        key = condition_key(
            value(
                row,
                "temperature_K",
            ),
            value(
                row,
                "kappa_ref_ps_inv",
            ),
            value(
                row,
                "gamma_phi_ps_inv",
            ),
        )

        condition_lookup[key] = row

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

                current_frame_rows = select_rows(
                    frame_rows,
                    temperature_K,
                    kappa,
                    gamma,
                )

                current_steady_rows = select_rows(
                    steady_rows,
                    temperature_K,
                    kappa,
                    gamma,
                )

                high_rows = [
                    row
                    for row
                    in current_frame_rows
                    if row[
                        "initial_state"
                    ] in HIGH_INITIAL_LABELS
                ]

                if len(high_rows) != (
                    EXPECTED_HIGH_INITIAL_ROWS
                ):
                    raise RuntimeError(
                        f"Expected "
                        f"{EXPECTED_HIGH_INITIAL_ROWS} "
                        f"high-state rows for {key}, "
                        f"found {len(high_rows)}"
                    )

                if len(current_steady_rows) != (
                    EXPECTED_FRAMES
                ):
                    raise RuntimeError(
                        f"Expected {EXPECTED_FRAMES} "
                        f"steady-state rows for {key}, "
                        f"found "
                        f"{len(current_steady_rows)}"
                    )

                final_values = np.asarray(
                    [
                        value(
                            row,
                            "final_PYR5_population",
                        )
                        for row in high_rows
                    ],
                    dtype=np.float64,
                )

                steady_values = np.asarray(
                    [
                        value(
                            row,
                            "steady_PYR5_population",
                        )
                        for row
                        in current_steady_rows
                    ],
                    dtype=np.float64,
                )

                gateway_values = np.asarray(
                    [
                        value(
                            row,
                            "PYR4_gateway_fraction",
                        )
                        for row in high_rows
                    ],
                    dtype=np.float64,
                )

                gateway_weights = np.asarray(
                    [
                        value(
                            row,
                            "integrated_downward_flux_to_PYR5",
                        )
                        for row in high_rows
                    ],
                    dtype=np.float64,
                )

                final_mean = float(
                    np.mean(final_values)
                )

                final_sd = float(
                    np.std(
                        final_values,
                        ddof=0,
                    )
                )

                steady_mean = float(
                    np.mean(steady_values)
                )

                steady_sd = float(
                    np.std(
                        steady_values,
                        ddof=0,
                    )
                )

                (
                    gateway_mean,
                    gateway_sd,
                ) = weighted_mean_and_sd(
                    gateway_values,
                    gateway_weights,
                )

                reference = (
                    condition_lookup[key]
                )

                reference_final = value(
                    reference,
                    "mean_final_PYR5_population_from_high_states",
                )

                reference_steady = value(
                    reference,
                    "mean_steady_PYR5_population",
                )

                reference_gateway = value(
                    reference,
                    "PYR4_gateway_fraction",
                )

                final_error = abs(
                    final_mean
                    - reference_final
                )

                steady_error = abs(
                    steady_mean
                    - reference_steady
                )

                gateway_error = abs(
                    gateway_mean
                    - reference_gateway
                )

                if final_error > SUMMARY_TOL:
                    raise RuntimeError(
                        f"Final-PYR5 summary mismatch "
                        f"for {key}: {final_error:.3e}"
                    )

                if steady_error > SUMMARY_TOL:
                    raise RuntimeError(
                        f"Steady-PYR5 summary mismatch "
                        f"for {key}: {steady_error:.3e}"
                    )

                if gateway_error > SUMMARY_TOL:
                    raise RuntimeError(
                        f"Gateway summary mismatch "
                        f"for {key}: {gateway_error:.3e}"
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
                        "n_gateway_samples": (
                            gateway_values.size
                        ),
                        "weighted_mean_PYR4_gateway_fraction": (
                            gateway_mean
                        ),
                        "weighted_sd_PYR4_gateway_fraction": (
                            gateway_sd
                        ),
                        "total_gateway_weight": float(
                            np.sum(
                                gateway_weights
                            )
                        ),
                    }
                )

    return statistics_rows


def statistics_lookup(
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

    expected_shape = (
        len(TEMPERATURES_K),
        len(KAPPA_VALUES_PS_INV),
        len(GAMMA_VALUES_PS_INV),
        len(SITE_LABELS),
        times_ps.size,
        len(SITE_LABELS),
    )

    if mean_populations.shape != (
        expected_shape
    ):
        raise RuntimeError(
            f"Unexpected trajectory shape: "
            f"{mean_populations.shape}"
        )

    lookup = statistics_lookup(
        statistics_rows
    )

    figure, axes = plt.subplots(
        1,
        3,
        figsize=(16.2, 5.1),
        constrained_layout=True,
    )

    axis_a, axis_b, axis_c = axes

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
            pooled_mean - pooled_sd,
            0.0,
            1.0,
        )

        upper = np.clip(
            pooled_mean + pooled_sd,
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
        float(times_ps[-1]),
    )

    axis_a.set_ylim(
        0.0,
        0.125,
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
    )

    axis_a.text(
        0.02,
        0.97,
        "(a)",
        transform=axis_a.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
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
                markersize=5.2,
                capsize=2.2,
                elinewidth=0.9,
                alpha=0.96,
            )

    configure_gamma_axis(
        axis_b
    )

    axis_b.set_ylabel(
        "PYR5 population at 100 ps"
    )

    axis_b.set_ylim(
        0.0,
        0.125,
    )

    axis_b.text(
        0.02,
        0.97,
        "(b)",
        transform=axis_b.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
    )

    kappa_handles = [
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
            label=(
                rf"$\kappa={kappa:g}$"
            ),
        )
        for kappa
        in KAPPA_VALUES_PS_INV
    ]

    temperature_handles = [
        Line2D(
            [0],
            [0],
            color="black",
            linestyle=(
                TEMPERATURE_LINESTYLES[
                    temperature_K
                ]
            ),
            label=(
                f"{temperature_K:.0f} K"
            ),
        )
        for temperature_K
        in TEMPERATURES_K
    ]

    first_legend = axis_b.legend(
        handles=kappa_handles,
        title=(
            r"Relaxation scale "
            r"(ps$^{-1}$)"
        ),
        frameon=False,
        loc="upper left",
    )

    axis_b.add_artist(
        first_legend
    )

    axis_b.legend(
        handles=temperature_handles,
        title="Temperature",
        frameon=False,
        loc="center left",
        bbox_to_anchor=(
            0.0,
            0.54,
        ),
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
                            "mean_stationary_PYR5_population"
                        ]
                    )
                )

                sds.append(
                    float(
                        row[
                            "sd_stationary_PYR5_population"
                        ]
                    )
                )

            axis_c.errorbar(
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
                markersize=5.2,
                capsize=2.2,
                elinewidth=0.9,
                alpha=0.96,
            )

    axis_c.axhline(
        0.25,
        color="black",
        linestyle=":",
        linewidth=1.35,
    )

    axis_c.text(
        55.0,
        0.268,
        "Uniform four-state limit",
        ha="center",
        va="bottom",
        fontsize=8.8,
    )

    configure_gamma_axis(
        axis_c
    )

    axis_c.set_ylabel(
        (
            "Stationary PYR5 population\n"
            "(combined phenomenological Liouvillian)"
        )
    )

    axis_c.set_ylim(
        0.20,
        1.035,
    )

    axis_c.text(
        0.02,
        0.97,
        "(c)",
        transform=axis_c.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
    )

    axis_c.text(
        0.98,
        0.96,
        (
            "PYR4 gateway:\n"
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
            "alpha": 0.90,
        },
    )

    figure.suptitle(
        (
            "Combined coherent, dephasing, and "
            "thermal-relaxation dynamics"
        ),
        fontsize=14,
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
    lookup = statistics_lookup(
        statistics_rows
    )

    figure, axis = plt.subplots(
        figsize=(7.5, 5.4),
        constrained_layout=True,
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
                            "weighted_mean_PYR4_gateway_fraction"
                        ]
                    )
                )

                sds.append(
                    float(
                        row[
                            "weighted_sd_PYR4_gateway_fraction"
                        ]
                    )
                )

            axis.errorbar(
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
                markersize=5.6,
                capsize=2.4,
                elinewidth=0.9,
                alpha=0.96,
            )

    configure_gamma_axis(
        axis
    )

    axis.set_ylabel(
        (
            "PYR4 fraction of integrated "
            "downward PYR5 flux"
        )
    )

    axis.set_ylim(
        0.76,
        0.96,
    )

    axis.set_title(
        "Robustness of PYR4 as the dominant gateway"
    )

    kappa_handles = [
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
            label=(
                rf"$\kappa={kappa:g}$ "
                r"ps$^{-1}$"
            ),
        )
        for kappa
        in KAPPA_VALUES_PS_INV
    ]

    temperature_handles = [
        Line2D(
            [0],
            [0],
            color="black",
            linestyle=(
                TEMPERATURE_LINESTYLES[
                    temperature_K
                ]
            ),
            label=(
                f"{temperature_K:.0f} K"
            ),
        )
        for temperature_K
        in TEMPERATURES_K
    ]

    first_legend = axis.legend(
        handles=kappa_handles,
        title="Relaxation scale",
        frameon=False,
        loc="lower right",
    )

    axis.add_artist(
        first_legend
    )

    axis.legend(
        handles=temperature_handles,
        title="Temperature",
        frameon=False,
        loc="lower left",
    )

    axis.text(
        0.02,
        0.97,
        (
            "Weighted mean ± weighted SD\n"
            "across frame–initial-state contributions"
        ),
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=9,
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
    final_values = np.asarray(
        [
            float(
                row[
                    "mean_final_PYR5_population"
                ]
            )
            for row
            in statistics_rows
        ],
        dtype=np.float64,
    )

    stationary_values = np.asarray(
        [
            float(
                row[
                    "mean_stationary_PYR5_population"
                ]
            )
            for row
            in statistics_rows
        ],
        dtype=np.float64,
    )

    gateway_values = np.asarray(
        [
            float(
                row[
                    "weighted_mean_PYR4_gateway_fraction"
                ]
            )
            for row
            in statistics_rows
        ],
        dtype=np.float64,
    )

    gamma100_final = np.asarray(
        [
            float(
                row[
                    "mean_final_PYR5_population"
                ]
            )
            for row
            in statistics_rows
            if np.isclose(
                float(
                    row[
                        "gamma_phi_ps_inv"
                    ]
                ),
                100.0,
            )
        ],
        dtype=np.float64,
    )

    with CAPTIONS_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Figure Captions\n\n"
        )

        handle.write(
            "## Main figure\n\n"
        )

        handle.write(
            "**Combined coherent, dephasing, and "
            "thermal-relaxation dynamics.** "
            "(a) Ensemble-averaged PYR5 population "
            "for representative conditions "
            "\\(T=300\\,\\mathrm{K}\\) and "
            "\\(\\kappa=1\\,\\mathrm{ps}^{-1}\\). "
            "Solid lines are means over the three "
            "upper-manifold initial states and 21 "
            "Hamiltonian snapshots; shaded regions "
            "show one pooled ensemble standard deviation. "
            "(b) PYR5 population after 100 ps as a "
            "function of the local pure-dephasing rate. "
            "Symbols and colors identify the relaxation "
            "scale \\(\\kappa\\), while solid and dashed "
            "lines identify 150 and 300 K, respectively. "
            "Error bars show one standard deviation over "
            "63 frame–initial-state realizations. "
            "(c) Stationary PYR5 population of the "
            "combined phenomenological Liouvillian. "
            "Error bars show one standard deviation over "
            "21 Hamiltonian snapshots. The dotted line "
            "marks the uniform four-state limit "
            "\\(P_{\\mathrm{PYR5}}=1/4\\). Thermal "
            "relaxation alone drives PYR5 population "
            "toward unity, whereas sufficiently strong "
            "local dephasing drives the stationary state "
            "toward a non-thermal, nearly uniform "
            "distribution. Across all tested conditions, "
            "PYR4 supplies 84.3–91.4% of the integrated "
            "downward flux into PYR5.\n\n"
        )

        handle.write(
            "## Supplementary gateway figure\n\n"
        )

        handle.write(
            "**Robustness of PYR4 as the dominant "
            "gateway into PYR5.** Weighted mean PYR4 "
            "contribution to the integrated downward "
            "thermal flux into PYR5 as a function of "
            "the local pure-dephasing rate. Error bars "
            "show the weighted standard deviation across "
            "the 21 Hamiltonian snapshots and three "
            "upper-manifold initial states. PYR4 remains "
            "the dominant gateway for every tested "
            "temperature, relaxation scale, and "
            "dephasing rate.\n\n"
        )

        handle.write(
            "## Accepted quantitative ranges\n\n"
        )

        handle.write(
            f"- PYR5 population at 100 ps: "
            f"{np.min(final_values):.6f} to "
            f"{np.max(final_values):.6f}.\n"
        )

        handle.write(
            f"- PYR5 population at 100 ps for "
            f"`gamma_phi = 100 ps^-1`: "
            f"{np.min(gamma100_final):.6f} to "
            f"{np.max(gamma100_final):.6f}.\n"
        )

        handle.write(
            f"- Stationary PYR5 population: "
            f"{np.min(stationary_values):.6f} to "
            f"{np.max(stationary_values):.6f}.\n"
        )

        handle.write(
            f"- PYR4 gateway contribution: "
            f"{np.min(gateway_values):.6f} to "
            f"{np.max(gateway_values):.6f}.\n\n"
        )

        handle.write(
            "## Interpretation boundary\n\n"
        )

        handle.write(
            "The absolute values of "
            "\\(\\gamma_\\phi\\) and \\(\\kappa\\) are "
            "phenomenological sensitivity parameters. "
            "The calculations establish relative "
            "mechanistic robustness but do not provide "
            "microscopic relaxation times derived from "
            "a molecular spectral density.\n"
        )


def write_manifest() -> None:
    with MANIFEST_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Combined-Mechanism "
            "Figure Manifest V2\n\n"
        )

        handle.write(
            "## Main products\n\n"
        )

        handle.write(
            "- `figure_day020_combined_mechanism_main.png`\n"
        )

        handle.write(
            "- `figure_day020_combined_mechanism_main.pdf`\n"
        )

        handle.write(
            "- `figure_day020_gateway_supplementary.png`\n"
        )

        handle.write(
            "- `figure_day020_gateway_supplementary.pdf`\n"
        )

        handle.write(
            "- `table_day020_combined_mechanism_statistics.csv`\n"
        )

        handle.write(
            "- `FIGURE_CAPTIONS_DAY020.md`\n\n"
        )

        handle.write(
            "## Statistical conventions\n\n"
        )

        handle.write(
            "- Main panel (a): pooled ensemble mean and "
            "population SD across 21 Hamiltonian snapshots "
            "and three high-state initial conditions.\n"
        )

        handle.write(
            "- Main panel (b): mean and population SD "
            "across 63 frame–initial-state realizations.\n"
        )

        handle.write(
            "- Main panel (c): mean and population SD "
            "across 21 stationary states.\n"
        )

        handle.write(
            "- Supplementary gateway figure: flux-weighted "
            "mean and weighted SD across 63 "
            "frame–initial-state contributions.\n"
        )

        handle.write(
            "- Dephasing axes use the physical values "
            "`0, 1, 20, 100 ps^-1` on a symmetric-log "
            "scale with a linear region around zero.\n"
        )


def validate_outputs() -> None:
    expected_paths = [
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
    ]

    missing = [
        path
        for path in expected_paths
        if (
            not path.exists()
            or path.stat().st_size == 0
        )
    ]

    if missing:
        raise RuntimeError(
            "Missing or empty V2 products:\n"
            + "\n".join(
                str(path)
                for path in missing
            )
        )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    configure_matplotlib()
    validate_numerical_status()

    trajectory_data = (
        load_trajectory_data()
    )

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
            f"condition-summary rows, found "
            f"{len(condition_rows)}"
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
        statistics_rows
    )

    write_captions(
        statistics_rows
    )

    write_manifest()
    validate_outputs()

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
                "weighted_mean_PYR4_gateway_fraction"
            ]
        )
        for row in statistics_rows
    ]

    log(
        "Day020 combined-mechanism "
        "publication figures V2 completed."
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
        "Main figure: 1 PNG + 1 PDF"
    )

    log(
        "Supplementary figure: 1 PNG + 1 PDF"
    )

    log(
        "Final PYR5 range: "
        f"{min(final_values):.6f} to "
        f"{max(final_values):.6f}"
    )

    log(
        "Stationary PYR5 range: "
        f"{min(stationary_values):.6f} to "
        f"{max(stationary_values):.6f}"
    )

    log(
        "PYR4 gateway range: "
        f"{min(gateway_values):.6f} to "
        f"{max(gateway_values):.6f}"
    )

    log(
        f"Wrote: "
        f"{OUTPUT_ROOT.relative_to(PROJECT_ROOT)}"
    )


if __name__ == "__main__":
    main()
