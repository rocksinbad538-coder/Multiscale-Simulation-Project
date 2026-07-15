#!/usr/bin/env python3

from __future__ import annotations

import csv
import itertools
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

ANALYSIS_ROOT = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol/"
    "execution/08_nvt_mobile_100ps/"
    "matched_mobile_vs_frozen_water"
)

PAIRED_CSV = (
    ANALYSIS_ROOT
    / "matched_framewise_water_differences.csv"
)

SOURCE_SUMMARY_CSV = (
    ANALYSIS_ROOT
    / "matched_mobile_frozen_water_comparison_summary.csv"
)

BLOCK_EFFECT_CSV = (
    ANALYSIS_ROOT
    / "matched_water_10ps_block_effects.csv"
)

METRIC_SUMMARY_CSV = (
    ANALYSIS_ROOT
    / "matched_water_block_inference_summary.csv"
)

EVENT_KINETICS_CSV = (
    ANALYSIS_ROOT
    / "matched_water_reentry_event_kinetics.csv"
)

FINAL_SUMMARY_CSV = (
    ANALYSIS_ROOT
    / "matched_water_scientific_review_summary.csv"
)

REPORT_MD = (
    ANALYSIS_ROOT
    / "MATCHED_WATER_SCIENTIFIC_REVIEW_DAY022.md"
)

EXPECTED_FRAMES = 201
FRAME_INTERVAL_PS = 0.5
BLOCK_COUNT = 10
BOOTSTRAP_REPLICATES = 50000
RNG_SEED = 20260707

INFERENTIAL_METRICS = [
    "lumen_water_count",
    "lumen_number_density_nm-3",
]

for pyrene_index in range(1, 5):
    INFERENTIAL_METRICS.extend(
        [
            f"PYR{pyrene_index}_waterO_within_35pm",
            f"PYR{pyrene_index}_waterO_within_50pm",
        ]
    )


def require_file(path: Path) -> None:
    if (
        not path.exists()
        or not path.is_file()
        or path.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing or empty file: {path}"
        )


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

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

    fields: list[str] = []
    seen: set[str] = set()

    for row in rows:
        for field in row:
            if field not in seen:
                seen.add(field)
                fields.append(field)

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(
                {
                    field: row.get(field, "")
                    for field in fields
                }
            )


def as_float(
    row: dict[str, str],
    field: str,
) -> float:
    try:
        value = float(row[field])
    except (KeyError, ValueError) as error:
        raise RuntimeError(
            f"Invalid numeric field {field!r}"
        ) from error

    if not math.isfinite(value):
        raise RuntimeError(
            f"Non-finite field {field!r}: {value}"
        )

    return value


def exact_sign_flip_pvalue(
    block_differences: np.ndarray,
) -> float:
    values = np.asarray(
        block_differences,
        dtype=float,
    )

    if len(values) == 0:
        return math.nan

    observed = abs(
        float(values.mean())
    )

    extreme = 0
    total = 0

    for signs_tuple in itertools.product(
        (-1.0, 1.0),
        repeat=len(values),
    ):
        signs = np.asarray(
            signs_tuple,
            dtype=float,
        )

        statistic = abs(
            float(
                np.mean(
                    signs * values
                )
            )
        )

        if statistic >= observed - 1.0e-14:
            extreme += 1

        total += 1

    return extreme / total


def bootstrap_block_mean_ci(
    block_differences: np.ndarray,
    rng: np.random.Generator,
) -> tuple[float, float]:
    values = np.asarray(
        block_differences,
        dtype=float,
    )

    indices = rng.integers(
        0,
        len(values),
        size=(
            BOOTSTRAP_REPLICATES,
            len(values),
        ),
    )

    bootstrap_means = np.mean(
        values[indices],
        axis=1,
    )

    lower, upper = np.quantile(
        bootstrap_means,
        [0.025, 0.975],
    )

    return (
        float(lower),
        float(upper),
    )


def lag1_autocorrelation(
    values: np.ndarray,
) -> float:
    values = np.asarray(
        values,
        dtype=float,
    )

    if (
        len(values) < 3
        or float(values.std()) < 1.0e-15
    ):
        return 0.0

    return float(
        np.corrcoef(
            values[:-1],
            values[1:],
        )[0, 1]
    )


def positive_episodes(
    values: np.ndarray,
) -> list[tuple[int, int]]:
    positive = np.asarray(
        values,
        dtype=float,
    ) > 0.0

    episodes: list[tuple[int, int]] = []
    start: int | None = None

    for index, is_positive in enumerate(
        positive
    ):
        if is_positive and start is None:
            start = index

        if (
            not is_positive
            and start is not None
        ):
            episodes.append(
                (
                    start,
                    index - 1,
                )
            )

            start = None

    if start is not None:
        episodes.append(
            (
                start,
                len(positive) - 1,
            )
        )

    return episodes


def event_kinetics(
    dataset: str,
    values: np.ndarray,
) -> dict[str, object]:
    values = np.asarray(
        values,
        dtype=float,
    )

    episodes = positive_episodes(
        values
    )

    episode_durations = np.asarray(
        [
            (
                stop
                - start
                + 1
            )
            * FRAME_INTERVAL_PS
            for start, stop in episodes
        ],
        dtype=float,
    )

    positive_values = values[
        values > 0.0
    ]

    return {
        "dataset": dataset,
        "frames": len(values),
        "zero_occupancy_fraction": float(
            np.mean(
                values == 0.0
            )
        ),
        "positive_occupancy_fraction": float(
            np.mean(
                values > 0.0
            )
        ),
        "mean_occupancy": float(
            values.mean()
        ),
        "conditional_mean_when_positive": (
            float(
                positive_values.mean()
            )
            if len(positive_values)
            else 0.0
        ),
        "maximum_occupancy": float(
            values.max()
        ),
        "positive_episode_count": len(
            episodes
        ),
        "mean_positive_episode_duration_ps": (
            float(
                episode_durations.mean()
            )
            if len(episode_durations)
            else 0.0
        ),
        "maximum_positive_episode_duration_ps": (
            float(
                episode_durations.max()
            )
            if len(episode_durations)
            else 0.0
        ),
        "sampled_positive_time_ps": float(
            np.count_nonzero(
                values > 0.0
            )
            * FRAME_INTERVAL_PS
        ),
        "lag1_autocorrelation": (
            lag1_autocorrelation(
                values
            )
        ),
    }


def main() -> None:
    paired_rows = read_csv(
        PAIRED_CSV
    )

    source_summary_rows = read_csv(
        SOURCE_SUMMARY_CSV
    )

    if len(paired_rows) != EXPECTED_FRAMES:
        raise RuntimeError(
            f"Expected {EXPECTED_FRAMES} paired rows; "
            f"found {len(paired_rows)}"
        )

    if len(source_summary_rows) != 1:
        raise RuntimeError(
            "Expected one matched-comparison summary row"
        )

    paired_rows.sort(
        key=lambda row: as_float(
            row,
            "relative_time_ps",
        )
    )

    times = np.asarray(
        [
            as_float(
                row,
                "relative_time_ps",
            )
            for row in paired_rows
        ],
        dtype=float,
    )

    expected_times = np.arange(
        EXPECTED_FRAMES,
        dtype=float,
    ) * FRAME_INTERVAL_PS

    if not np.allclose(
        times,
        expected_times,
        atol=1.0e-9,
        rtol=0.0,
    ):
        raise RuntimeError(
            "Paired frame times are not the expected "
            "0-100 ps sequence at 0.5 ps intervals"
        )

    block_masks = []

    for block_index in range(
        BLOCK_COUNT
    ):
        lower = 10.0 * block_index
        upper = lower + 10.0

        if block_index == BLOCK_COUNT - 1:
            mask = (
                (times >= lower)
                & (times <= upper)
            )
        else:
            mask = (
                (times >= lower)
                & (times < upper)
            )

        if not np.any(mask):
            raise RuntimeError(
                f"No frames in block {block_index + 1}"
            )

        block_masks.append(mask)

    rng = np.random.default_rng(
        RNG_SEED
    )

    block_rows: list[
        dict[str, object]
    ] = []

    metric_summary_rows: list[
        dict[str, object]
    ] = []

    for metric in INFERENTIAL_METRICS:
        mobile_column = (
            f"mobile_{metric}"
        )

        frozen_column = (
            f"matched_frozen_{metric}"
        )

        difference_column = (
            f"mobile_minus_frozen_{metric}"
        )

        mobile_values = np.asarray(
            [
                as_float(
                    row,
                    mobile_column,
                )
                for row in paired_rows
            ],
            dtype=float,
        )

        frozen_values = np.asarray(
            [
                as_float(
                    row,
                    frozen_column,
                )
                for row in paired_rows
            ],
            dtype=float,
        )

        differences = np.asarray(
            [
                as_float(
                    row,
                    difference_column,
                )
                for row in paired_rows
            ],
            dtype=float,
        )

        block_differences = []

        for block_index, mask in enumerate(
            block_masks,
            start=1,
        ):
            mobile_mean = float(
                mobile_values[mask].mean()
            )

            frozen_mean = float(
                frozen_values[mask].mean()
            )

            difference_mean = float(
                differences[mask].mean()
            )

            block_differences.append(
                difference_mean
            )

            block_rows.append(
                {
                    "metric": metric,
                    "block": block_index,
                    "relative_start_ps": (
                        10.0
                        * (
                            block_index - 1
                        )
                    ),
                    "relative_end_ps": (
                        10.0
                        * block_index
                    ),
                    "frame_count": int(
                        np.count_nonzero(
                            mask
                        )
                    ),
                    "mobile_mean": mobile_mean,
                    "matched_frozen_mean": frozen_mean,
                    "mobile_minus_frozen_mean": (
                        difference_mean
                    ),
                }
            )

        block_array = np.asarray(
            block_differences,
            dtype=float,
        )

        ci_lower, ci_upper = (
            bootstrap_block_mean_ci(
                block_array,
                rng,
            )
        )

        sign_flip_p = (
            exact_sign_flip_pvalue(
                block_array
            )
        )

        block_std = float(
            block_array.std(
                ddof=1
            )
        )

        standardized_paired_effect = (
            float(
                block_array.mean()
                / block_std
            )
            if block_std > 1.0e-15
            else math.inf
        )

        frozen_mean = float(
            frozen_values.mean()
        )

        relative_difference_pct = (
            100.0
            * float(
                differences.mean()
            )
            / frozen_mean
            if abs(frozen_mean) > 1.0e-15
            else math.nan
        )

        screening_supported = (
            ci_lower > 0.0
            and sign_flip_p <= 0.05
            and float(
                np.mean(
                    block_array > 0.0
                )
            )
            >= 0.70
        )

        metric_summary_rows.append(
            {
                "metric": metric,
                "mobile_frame_mean": float(
                    mobile_values.mean()
                ),
                "matched_frozen_frame_mean": (
                    frozen_mean
                ),
                "frame_level_mean_difference": float(
                    differences.mean()
                ),
                "relative_difference_pct_secondary": (
                    relative_difference_pct
                ),
                "block_mean_difference": float(
                    block_array.mean()
                ),
                "block_difference_std": (
                    block_std
                ),
                "block_bootstrap_CI95_lower": (
                    ci_lower
                ),
                "block_bootstrap_CI95_upper": (
                    ci_upper
                ),
                "exact_sign_flip_pvalue": (
                    sign_flip_p
                ),
                "positive_block_fraction": float(
                    np.mean(
                        block_array > 0.0
                    )
                ),
                "zero_block_fraction": float(
                    np.mean(
                        block_array == 0.0
                    )
                ),
                "negative_block_fraction": float(
                    np.mean(
                        block_array < 0.0
                    )
                ),
                "standardized_paired_block_effect": (
                    standardized_paired_effect
                ),
                "screening_effect_supported": (
                    screening_supported
                ),
            }
        )

    write_csv(
        BLOCK_EFFECT_CSV,
        block_rows,
    )

    write_csv(
        METRIC_SUMMARY_CSV,
        metric_summary_rows,
    )

    occupancy_mobile = np.asarray(
        [
            as_float(
                row,
                "mobile_lumen_water_count",
            )
            for row in paired_rows
        ],
        dtype=float,
    )

    occupancy_frozen = np.asarray(
        [
            as_float(
                row,
                "matched_frozen_lumen_water_count",
            )
            for row in paired_rows
        ],
        dtype=float,
    )

    event_rows = [
        event_kinetics(
            "mobile",
            occupancy_mobile,
        ),
        event_kinetics(
            "matched_frozen",
            occupancy_frozen,
        ),
    ]

    write_csv(
        EVENT_KINETICS_CSV,
        event_rows,
    )

    metric_by_name = {
        row["metric"]: row
        for row in metric_summary_rows
    }

    occupancy_inference = metric_by_name[
        "lumen_water_count"
    ]

    occupancy_supported = bool(
        occupancy_inference[
            "screening_effect_supported"
        ]
    )

    hydration_supported_metrics = [
        str(row["metric"])
        for row in metric_summary_rows
        if (
            str(row["metric"]).startswith(
                "PYR"
            )
            and bool(
                row[
                    "screening_effect_supported"
                ]
            )
        )
    ]

    hydration_suggestive_metrics = [
        str(row["metric"])
        for row in metric_summary_rows
        if (
            str(row["metric"]).startswith(
                "PYR"
            )
            and float(
                row[
                    "block_bootstrap_CI95_lower"
                ]
            )
            > 0.0
        )
    ]

    source_summary = (
        source_summary_rows[0]
    )

    radial_js = float(
        source_summary[
            "radial_JS_divergence"
        ]
    )

    axial_js = float(
        source_summary[
            "axial_JS_divergence"
        ]
    )

    if occupancy_supported:
        screening_interpretation = (
            "MOBILITY_ASSOCIATED_TRANSIENT_"
            "REHYDRATION_SUPPORTED"
        )

        next_step = (
            "REVIEW_21_MATCHED_SNAPSHOT_PAIRS_"
            "FOR_LIMITED_MOBILE_QM_PILOT"
        )

    elif hydration_supported_metrics:
        screening_interpretation = (
            "LOCAL_PYRENE_HYDRATION_"
            "REORGANIZATION_SUPPORTED"
        )

        next_step = (
            "REVIEW_LOCAL_HYDRATION_SNAPSHOT_"
            "PAIRS_BEFORE_QM"
        )

    else:
        screening_interpretation = (
            "INCONCLUSIVE_AT_CURRENT_"
            "SINGLE_TRAJECTORY_SAMPLING"
        )

        next_step = (
            "CONSIDER_REPLICATE_TRAJECTORIES_"
            "BEFORE_ELECTRONIC_RECALCULATION"
        )

    mobile_events = event_rows[0]
    frozen_events = event_rows[1]

    final_summary = {
        "frames_per_dataset": (
            EXPECTED_FRAMES
        ),
        "matched_time_window_ps": (
            "44_to_144"
        ),
        "block_count": (
            BLOCK_COUNT
        ),
        "block_duration_ps": 10.0,
        "bootstrap_replicates": (
            BOOTSTRAP_REPLICATES
        ),
        "occupancy_mobile_mean": float(
            occupancy_mobile.mean()
        ),
        "occupancy_matched_frozen_mean": float(
            occupancy_frozen.mean()
        ),
        "occupancy_absolute_difference": float(
            occupancy_mobile.mean()
            - occupancy_frozen.mean()
        ),
        "occupancy_block_CI95_lower": (
            occupancy_inference[
                "block_bootstrap_CI95_lower"
            ]
        ),
        "occupancy_block_CI95_upper": (
            occupancy_inference[
                "block_bootstrap_CI95_upper"
            ]
        ),
        "occupancy_exact_sign_flip_pvalue": (
            occupancy_inference[
                "exact_sign_flip_pvalue"
            ]
        ),
        "occupancy_positive_block_fraction": (
            occupancy_inference[
                "positive_block_fraction"
            ]
        ),
        "mobile_zero_occupancy_fraction": (
            mobile_events[
                "zero_occupancy_fraction"
            ]
        ),
        "matched_frozen_zero_occupancy_fraction": (
            frozen_events[
                "zero_occupancy_fraction"
            ]
        ),
        "mobile_positive_episode_count": (
            mobile_events[
                "positive_episode_count"
            ]
        ),
        "matched_frozen_positive_episode_count": (
            frozen_events[
                "positive_episode_count"
            ]
        ),
        "mobile_longest_positive_episode_ps": (
            mobile_events[
                "maximum_positive_episode_duration_ps"
            ]
        ),
        "matched_frozen_longest_positive_episode_ps": (
            frozen_events[
                "maximum_positive_episode_duration_ps"
            ]
        ),
        "radial_JS_divergence_descriptive": (
            radial_js
        ),
        "axial_JS_divergence_descriptive": (
            axial_js
        ),
        "occupancy_screening_effect_supported": (
            occupancy_supported
        ),
        "supported_PYR_hydration_metric_count": (
            len(
                hydration_supported_metrics
            )
        ),
        "supported_PYR_hydration_metrics": (
            " | ".join(
                hydration_supported_metrics
            )
        ),
        "suggestive_PYR_hydration_metrics": (
            " | ".join(
                hydration_suggestive_metrics
            )
        ),
        "screening_interpretation": (
            screening_interpretation
        ),
        "publication_level_causal_claim_authorized": (
            False
        ),
        "electronic_recalculation_authorized": (
            False
        ),
        "longer_mobile_production_authorized": (
            False
        ),
        "authorized_next_step": (
            next_step
        ),
    }

    write_csv(
        FINAL_SUMMARY_CSV,
        [final_summary],
    )

    REPORT_MD.write_text(
        f"""# Matched Water Scientific Review

## Scope

- Matched branch interval: **44-144 ps**
- Frames per trajectory: **201**
- Coordinate interval: **0.5 ps**
- Temporal blocks: **10 × 10 ps**
- Block-bootstrap replicates: **{BOOTSTRAP_REPLICATES}**

This is a screening-level paired time-series analysis. The ten
temporal blocks reduce frame-level autocorrelation but do not replace
independent trajectory replicas.

## Lumen occupancy

- Mobile mean: {occupancy_mobile.mean():.6f}
- Matched frozen mean: {occupancy_frozen.mean():.6f}
- Absolute mobile-minus-frozen difference:
  {final_summary['occupancy_absolute_difference']:.6f} waters
- Block-bootstrap 95% interval:
  [{final_summary['occupancy_block_CI95_lower']:.6f},
  {final_summary['occupancy_block_CI95_upper']:.6f}]
- Exact block sign-flip p value:
  {final_summary['occupancy_exact_sign_flip_pvalue']:.8f}
- Fraction of blocks with positive mobile-minus-frozen difference:
  {final_summary['occupancy_positive_block_fraction']:.6f}
- Screening effect supported:
  **{'YES' if occupancy_supported else 'NO'}**

The previously reported 316% relative difference is secondary because
the matched-frozen mean is below one water molecule.

## Reentry-event kinetics

- Mobile zero-occupancy fraction:
  {mobile_events['zero_occupancy_fraction']:.6f}
- Matched-frozen zero-occupancy fraction:
  {frozen_events['zero_occupancy_fraction']:.6f}
- Mobile/matched-frozen positive episodes:
  {mobile_events['positive_episode_count']}/
  {frozen_events['positive_episode_count']}
- Longest mobile/matched-frozen positive episode:
  {mobile_events['maximum_positive_episode_duration_ps']:.3f}/
  {frozen_events['maximum_positive_episode_duration_ps']:.3f} ps

## Spatial distributions

- Radial Jensen-Shannon divergence:
  {radial_js:.8f}
- Axial Jensen-Shannon divergence:
  {axial_js:.8f}

These divergences remain descriptive because the distributions contain
few lumen-water observations and no independent trajectory replicas.

## Local pyrene hydration

- Supported PYR/cutoff metrics:
  {', '.join(hydration_supported_metrics) if hydration_supported_metrics else 'NONE'}
- Suggestive PYR/cutoff metrics:
  {', '.join(hydration_suggestive_metrics) if hydration_suggestive_metrics else 'NONE'}

## Decision

- Screening interpretation:
  **{screening_interpretation}**
- Publication-level causal claim authorized: **NO**
- Electronic recalculation authorized: **NO**
- Longer mobile production authorized: **NO**
- Authorized next step:
  `{next_step}`

The scientifically appropriate interpretation remains transient
rehydration or local hydration reorganization, not stable confined
water retention.
""",
        encoding="utf-8",
    )

    print(
        "Day022 matched-water scientific review completed."
    )

    print(
        "Matched frames / blocks: "
        f"{EXPECTED_FRAMES} / {BLOCK_COUNT}"
    )

    print(
        "Mobile / matched-frozen occupancy mean: "
        f"{occupancy_mobile.mean():.6f}/"
        f"{occupancy_frozen.mean():.6f}"
    )

    print(
        "Absolute occupancy difference: "
        f"{final_summary['occupancy_absolute_difference']:.6f} waters"
    )

    print(
        "Occupancy block-bootstrap 95% CI: "
        f"[{final_summary['occupancy_block_CI95_lower']:.6f}, "
        f"{final_summary['occupancy_block_CI95_upper']:.6f}]"
    )

    print(
        "Occupancy exact sign-flip p-value: "
        f"{final_summary['occupancy_exact_sign_flip_pvalue']:.8f}"
    )

    print(
        "Positive occupancy-difference block fraction: "
        f"{final_summary['occupancy_positive_block_fraction']:.6f}"
    )

    print(
        "Mobile/frozen zero-occupancy fractions: "
        f"{mobile_events['zero_occupancy_fraction']:.6f}/"
        f"{frozen_events['zero_occupancy_fraction']:.6f}"
    )

    print(
        "Mobile/frozen positive episode counts: "
        f"{mobile_events['positive_episode_count']}/"
        f"{frozen_events['positive_episode_count']}"
    )

    print(
        "Mobile/frozen longest positive episode: "
        f"{mobile_events['maximum_positive_episode_duration_ps']:.3f}/"
        f"{frozen_events['maximum_positive_episode_duration_ps']:.3f} ps"
    )

    print(
        "Supported PYR hydration metrics: "
        + (
            " | ".join(
                hydration_supported_metrics
            )
            if hydration_supported_metrics
            else "NONE"
        )
    )

    print(
        "Screening interpretation: "
        f"{screening_interpretation}"
    )

    print(
        "Publication-level causal claim authorized: NO"
    )

    print(
        "Electronic recalculation authorized: NO"
    )

    print(
        "Longer mobile production authorized: NO"
    )

    print(
        f"Authorized next step: {next_step}"
    )

    print(
        "Wrote: "
        + str(
            REPORT_MD.relative_to(ROOT)
        )
    )


if __name__ == "__main__":
    main()
