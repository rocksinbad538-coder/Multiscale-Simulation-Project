#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

ANALYSIS_ROOT = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol/"
    "execution/08_nvt_mobile_100ps/"
    "matched_mobile_vs_frozen_water"
)

TIMESERIES_CSV = (
    ANALYSIS_ROOT
    / "matched_mobile_frozen_water_timeseries.csv"
)

SNAPSHOT_CSV = (
    ANALYSIS_ROOT
    / "matched_representative_snapshot_pairs.csv"
)

REVIEW_CSV = (
    ANALYSIS_ROOT
    / "matched_snapshot_pair_review.csv"
)

PILOT_SELECTION_CSV = (
    ANALYSIS_ROOT
    / "limited_mobile_frozen_qm_pilot_candidates.csv"
)

SUMMARY_CSV = (
    ANALYSIS_ROOT
    / "matched_snapshot_pair_review_summary.csv"
)

REPORT_MD = (
    ANALYSIS_ROOT
    / "MATCHED_SNAPSHOT_PAIR_REVIEW_DAY022.md"
)

EXPECTED_PAIRS = 21

LOCAL_PATTERN = re.compile(
    r"^PYR(\d+)_waterO_within_(35|50)pm$"
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


def read_csv(
    path: Path,
) -> list[dict[str, str]]:
    require_file(path)

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        return list(
            csv.DictReader(handle)
        )


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
                fields.append(field)
                seen.add(field)

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
                    field: row.get(
                        field,
                        "",
                    )
                    for field in fields
                }
            )


def parse_float(
    row: dict[str, str],
    field: str,
) -> float:
    try:
        value = float(
            row[field]
        )
    except (
        KeyError,
        ValueError,
        TypeError,
    ) as error:
        raise RuntimeError(
            f"Invalid numeric field {field!r}"
        ) from error

    if not math.isfinite(value):
        raise RuntimeError(
            f"Non-finite value in {field!r}: "
            f"{value}"
        )

    return value


def rounded_time(value: float) -> float:
    return round(
        float(value),
        6,
    )


def find_snapshot_time_column(
    rows: list[dict[str, str]],
    valid_mobile_times: set[float],
) -> str:
    if not rows:
        raise RuntimeError(
            "Snapshot-pair CSV is empty"
        )

    fields = list(
        rows[0].keys()
    )

    preferred = [
        "relative_time_ps",
        "matched_relative_time_ps",
        "mobile_relative_time_ps",
        "mobile_time_ps",
        "target_time_ps",
        "time_ps",
    ]

    candidates: list[str] = []

    for field in preferred + fields:
        if field not in fields:
            continue

        if field in candidates:
            continue

        lowered = field.lower()

        if (
            "time" not in lowered
            or "ps" not in lowered
        ):
            continue

        try:
            values = [
                rounded_time(
                    float(
                        row[field]
                    )
                )
                for row in rows
            ]
        except (
            ValueError,
            TypeError,
            KeyError,
        ):
            continue

        if len(set(values)) != len(values):
            continue

        if not all(
            value in valid_mobile_times
            for value in values
        ):
            continue

        candidates.append(
            field
        )

    if not candidates:
        raise RuntimeError(
            "Could not identify the mobile-relative "
            "time column in the representative-pair CSV. "
            f"Available fields: {' | '.join(fields)}"
        )

    for field in preferred:
        if field in candidates:
            return field

    return candidates[0]


def temporal_quartile(
    time_ps: float,
) -> int:
    if time_ps < 25.0:
        return 1

    if time_ps < 50.0:
        return 2

    if time_ps < 75.0:
        return 3

    return 4


def add_selection_reason(
    selected: dict[float, dict[str, object]],
    row: dict[str, object],
    reason: str,
) -> None:
    time_ps = float(
        row["relative_time_ps"]
    )

    if time_ps not in selected:
        selected[time_ps] = {
            **row,
            "selection_reasons": reason,
        }

        return

    existing = str(
        selected[time_ps][
            "selection_reasons"
        ]
    )

    reasons = existing.split(
        " | "
    )

    if reason not in reasons:
        selected[time_ps][
            "selection_reasons"
        ] = (
            existing
            + " | "
            + reason
        )


def greedily_select(
    candidates: list[dict[str, object]],
    selected: dict[float, dict[str, object]],
    limit: int,
    reason: str,
    minimum_separation_ps: float,
) -> int:
    added = 0

    for row in candidates:
        if added >= limit:
            break

        time_ps = float(
            row["relative_time_ps"]
        )

        if time_ps in selected:
            add_selection_reason(
                selected,
                row,
                reason,
            )

            continue

        if any(
            abs(
                time_ps
                - existing_time
            )
            < minimum_separation_ps
            for existing_time in selected
        ):
            continue

        add_selection_reason(
            selected,
            row,
            reason,
        )

        added += 1

    return added


def main() -> None:
    timeseries_rows = read_csv(
        TIMESERIES_CSV
    )

    snapshot_rows = read_csv(
        SNAPSHOT_CSV
    )

    if len(snapshot_rows) != EXPECTED_PAIRS:
        raise RuntimeError(
            f"Expected {EXPECTED_PAIRS} representative "
            f"snapshot pairs; found {len(snapshot_rows)}"
        )

    fields = set(
        timeseries_rows[0].keys()
    )

    required_fields = {
        "dataset",
        "relative_time_ps",
        "lumen_water_count",
    }

    missing_fields = (
        required_fields
        - fields
    )

    if missing_fields:
        raise RuntimeError(
            "Missing required time-series fields: "
            + " | ".join(
                sorted(
                    missing_fields
                )
            )
        )

    local_fields = sorted(
        field
        for field in fields
        if LOCAL_PATTERN.match(
            field
        )
    )

    if not local_fields:
        raise RuntimeError(
            "No PYR local-hydration fields were found"
        )

    local_35_fields = [
        field
        for field in local_fields
        if field.endswith(
            "_35pm"
        )
    ]

    local_50_fields = [
        field
        for field in local_fields
        if field.endswith(
            "_50pm"
        )
    ]

    mobile_rows = [
        row
        for row in timeseries_rows
        if row[
            "dataset"
        ].strip().lower()
        == "mobile"
    ]

    frozen_rows = [
        row
        for row in timeseries_rows
        if row[
            "dataset"
        ].strip().lower()
        == "frozen"
    ]

    if (
        len(mobile_rows) != 201
        or len(frozen_rows) != 201
    ):
        raise RuntimeError(
            "Expected 201 mobile and 201 frozen "
            "time-series rows"
        )

    mobile_map = {
        rounded_time(
            parse_float(
                row,
                "relative_time_ps",
            )
        ): row
        for row in mobile_rows
    }

    frozen_map = {
        rounded_time(
            parse_float(
                row,
                "relative_time_ps",
            )
        ): row
        for row in frozen_rows
    }

    if set(
        mobile_map
    ) != set(
        frozen_map
    ):
        raise RuntimeError(
            "Mobile and frozen relative time grids differ"
        )

    snapshot_time_field = (
        find_snapshot_time_column(
            snapshot_rows,
            set(
                mobile_map
            ),
        )
    )

    selected_times = [
        rounded_time(
            parse_float(
                row,
                snapshot_time_field,
            )
        )
        for row in snapshot_rows
    ]

    if len(
        set(
            selected_times
        )
    ) != EXPECTED_PAIRS:
        raise RuntimeError(
            "Representative snapshot times are not unique"
        )

    review_rows: list[
        dict[str, object]
    ] = []

    for pair_index, time_ps in enumerate(
        sorted(
            selected_times
        ),
        start=1,
    ):
        if (
            time_ps not in mobile_map
            or time_ps not in frozen_map
        ):
            raise RuntimeError(
                "Representative time is absent from "
                f"the matched time series: {time_ps}"
            )

        mobile = mobile_map[
            time_ps
        ]

        frozen = frozen_map[
            time_ps
        ]

        mobile_lumen = parse_float(
            mobile,
            "lumen_water_count",
        )

        frozen_lumen = parse_float(
            frozen,
            "lumen_water_count",
        )

        lumen_difference = (
            mobile_lumen
            - frozen_lumen
        )

        local_35_differences: list[
            float
        ] = []

        local_50_differences: list[
            float
        ] = []

        output: dict[str, object] = {
            "pair_index": pair_index,
            "relative_time_ps": time_ps,
            "matched_frozen_absolute_time_ps": (
                44.0 + time_ps
            ),
            "temporal_quartile": (
                temporal_quartile(
                    time_ps
                )
            ),
            "mobile_lumen_water_count": (
                mobile_lumen
            ),
            "matched_frozen_lumen_water_count": (
                frozen_lumen
            ),
            "mobile_minus_frozen_lumen_count": (
                lumen_difference
            ),
            "absolute_lumen_difference": abs(
                lumen_difference
            ),
        }

        for field in local_fields:
            mobile_value = parse_float(
                mobile,
                field,
            )

            frozen_value = parse_float(
                frozen,
                field,
            )

            difference = (
                mobile_value
                - frozen_value
            )

            output[
                f"mobile_{field}"
            ] = mobile_value

            output[
                f"matched_frozen_{field}"
            ] = frozen_value

            output[
                f"mobile_minus_frozen_{field}"
            ] = difference

            if field in local_35_fields:
                local_35_differences.append(
                    difference
                )

            if field in local_50_fields:
                local_50_differences.append(
                    difference
                )

        max_abs_local_35 = (
            max(
                abs(value)
                for value
                in local_35_differences
            )
            if local_35_differences
            else 0.0
        )

        max_abs_local_50 = (
            max(
                abs(value)
                for value
                in local_50_differences
            )
            if local_50_differences
            else 0.0
        )

        sum_abs_local_50 = float(
            sum(
                abs(value)
                for value
                in local_50_differences
            )
        )

        local_50_contrast_count = int(
            sum(
                abs(value) >= 2.0
                for value
                in local_50_differences
            )
        )

        strong_local_contrast = bool(
            max_abs_local_50 >= 2.0
            or max_abs_local_35 >= 1.0
        )

        strong_lumen_contrast = bool(
            abs(
                lumen_difference
            )
            >= 2.0
        )

        control_like = bool(
            abs(
                lumen_difference
            )
            <= 1.0
            and max_abs_local_35 <= 1.0
            and max_abs_local_50 <= 1.0
        )

        output.update(
            {
                "max_absolute_local_35pm_difference": (
                    max_abs_local_35
                ),
                "max_absolute_local_50pm_difference": (
                    max_abs_local_50
                ),
                "sum_absolute_local_50pm_differences": (
                    sum_abs_local_50
                ),
                "local_50pm_contrast_count_ge_2": (
                    local_50_contrast_count
                ),
                "strong_local_contrast": (
                    strong_local_contrast
                ),
                "strong_lumen_contrast": (
                    strong_lumen_contrast
                ),
                "control_like_pair": (
                    control_like
                ),
                "occupancy_direction": (
                    "MOBILE_HIGHER"
                    if lumen_difference > 0.0
                    else (
                        "FROZEN_HIGHER"
                        if lumen_difference < 0.0
                        else "EQUAL"
                    )
                ),
                "informative_screening_score": (
                    3.0
                    * max_abs_local_50
                    + 2.0
                    * max_abs_local_35
                    + abs(
                        lumen_difference
                    )
                    + 0.5
                    * sum_abs_local_50
                ),
            }
        )

        review_rows.append(
            output
        )

    write_csv(
        REVIEW_CSV,
        review_rows,
    )

    strong_local_rows = sorted(
        [
            row
            for row in review_rows
            if bool(
                row[
                    "strong_local_contrast"
                ]
            )
        ],
        key=lambda row: (
            float(
                row[
                    "max_absolute_local_50pm_difference"
                ]
            ),
            float(
                row[
                    "max_absolute_local_35pm_difference"
                ]
            ),
            float(
                row[
                    "absolute_lumen_difference"
                ]
            ),
        ),
        reverse=True,
    )

    strong_lumen_rows = sorted(
        [
            row
            for row in review_rows
            if bool(
                row[
                    "strong_lumen_contrast"
                ]
            )
        ],
        key=lambda row: (
            float(
                row[
                    "absolute_lumen_difference"
                ]
            ),
            float(
                row[
                    "max_absolute_local_50pm_difference"
                ]
            ),
        ),
        reverse=True,
    )

    control_rows = sorted(
        [
            row
            for row in review_rows
            if bool(
                row[
                    "control_like_pair"
                ]
            )
        ],
        key=lambda row: (
            float(
                row[
                    "absolute_lumen_difference"
                ]
            )
            + float(
                row[
                    "max_absolute_local_50pm_difference"
                ]
            )
            + float(
                row[
                    "max_absolute_local_35pm_difference"
                ]
            ),
            abs(
                float(
                    row[
                        "relative_time_ps"
                    ]
                )
                - 50.0
            ),
        ),
    )

    selected: dict[
        float,
        dict[str, object]
    ] = {}

    greedily_select(
        strong_local_rows,
        selected,
        limit=3,
        reason="LOCAL_PYRENE_HYDRATION_CONTRAST",
        minimum_separation_ps=10.0,
    )

    greedily_select(
        strong_lumen_rows,
        selected,
        limit=2,
        reason="LUMEN_REHYDRATION_CONTRAST",
        minimum_separation_ps=7.5,
    )

    greedily_select(
        control_rows,
        selected,
        limit=2,
        reason="MATCHED_LOW_CONTRAST_CONTROL",
        minimum_separation_ps=15.0,
    )

    if len(selected) < 5:
        remaining = sorted(
            [
                row
                for row in review_rows
                if float(
                    row[
                        "relative_time_ps"
                    ]
                ) not in selected
            ],
            key=lambda row: float(
                row[
                    "informative_screening_score"
                ]
            ),
            reverse=True,
        )

        for row in remaining:
            if len(selected) >= 5:
                break

            add_selection_reason(
                selected,
                row,
                "TEMPORAL_AND_INFORMATION_COVERAGE",
            )

    selected_rows = sorted(
        selected.values(),
        key=lambda row: float(
            row[
                "relative_time_ps"
            ]
        ),
    )

    for selection_index, row in enumerate(
        selected_rows,
        start=1,
    ):
        row[
            "pilot_selection_index"
        ] = selection_index

    write_csv(
        PILOT_SELECTION_CSV,
        selected_rows,
    )

    selected_times_array = np.asarray(
        [
            float(
                row[
                    "relative_time_ps"
                ]
            )
            for row in selected_rows
        ],
        dtype=float,
    )

    selected_quartiles = {
        int(
            row[
                "temporal_quartile"
            ]
        )
        for row in selected_rows
    }

    selected_local_count = sum(
        bool(
            row[
                "strong_local_contrast"
            ]
        )
        for row in selected_rows
    )

    selected_control_count = sum(
        bool(
            row[
                "control_like_pair"
            ]
        )
        for row in selected_rows
    )

    selected_informative_count = sum(
        bool(
            row[
                "strong_local_contrast"
            ]
        )
        or bool(
            row[
                "strong_lumen_contrast"
            ]
        )
        for row in selected_rows
    )

    temporal_span_ps = (
        float(
            selected_times_array.max()
            - selected_times_array.min()
        )
        if len(
            selected_times_array
        )
        else 0.0
    )

    gates = {
        "all_21_pairs_valid": (
            len(
                review_rows
            )
            == EXPECTED_PAIRS
        ),
        "at_least_2_local_contrast_pairs": (
            len(
                strong_local_rows
            )
            >= 2
        ),
        "at_least_2_control_like_pairs": (
            len(
                control_rows
            )
            >= 2
        ),
        "selected_5_to_7_pairs": (
            5
            <= len(
                selected_rows
            )
            <= 7
        ),
        "at_least_2_selected_local_contrasts": (
            selected_local_count
            >= 2
        ),
        "at_least_3_selected_informative_pairs": (
            selected_informative_count
            >= 3
        ),
        "at_least_3_temporal_quartiles": (
            len(
                selected_quartiles
            )
            >= 3
        ),
        "selected_temporal_span_at_least_50ps": (
            temporal_span_ps
            >= 50.0
        ),
        "at_least_1_selected_control": (
            selected_control_count
            >= 1
        ),
    }

    failed_gates = [
        name
        for name, passed
        in gates.items()
        if not passed
    ]

    limited_qm_pilot_authorized = (
        len(
            failed_gates
        )
        == 0
    )

    if limited_qm_pilot_authorized:
        decision = (
            "LIMITED_PAIRED_QM_PILOT_JUSTIFIED"
        )

        authorized_next_step = (
            "PREPARE_SELECTED_PAIRED_QM_INPUTS_"
            "WITHOUT_STARTING_FULL_RECALCULATION"
        )
    else:
        decision = (
            "LIMITED_PAIRED_QM_PILOT_NOT_JUSTIFIED"
        )

        authorized_next_step = (
            "RETAIN_FROZEN_QM_BASELINE_AND_DOCUMENT_"
            "MOBILITY_ASSOCIATED_SOLVENT_KINETICS"
        )

    summary = {
        "representative_pairs_audited": (
            len(
                review_rows
            )
        ),
        "snapshot_time_field": (
            snapshot_time_field
        ),
        "strong_local_contrast_pair_count": (
            len(
                strong_local_rows
            )
        ),
        "strong_lumen_contrast_pair_count": (
            len(
                strong_lumen_rows
            )
        ),
        "control_like_pair_count": (
            len(
                control_rows
            )
        ),
        "selected_pair_count": (
            len(
                selected_rows
            )
        ),
        "selected_relative_times_ps": (
            " | ".join(
                f"{float(row['relative_time_ps']):.3f}"
                for row in selected_rows
            )
        ),
        "selected_temporal_quartile_count": (
            len(
                selected_quartiles
            )
        ),
        "selected_temporal_span_ps": (
            temporal_span_ps
        ),
        "selected_local_contrast_count": (
            selected_local_count
        ),
        "selected_informative_pair_count": (
            selected_informative_count
        ),
        "selected_control_pair_count": (
            selected_control_count
        ),
        "failed_gates": (
            " | ".join(
                failed_gates
            )
        ),
        "decision": decision,
        "limited_mobile_frozen_QM_pilot_authorized": (
            limited_qm_pilot_authorized
        ),
        "full_electronic_recalculation_authorized": (
            False
        ),
        "longer_mobile_production_authorized": (
            False
        ),
        "publication_level_causal_claim_authorized": (
            False
        ),
        "authorized_next_step": (
            authorized_next_step
        ),
    }

    write_csv(
        SUMMARY_CSV,
        [summary],
    )

    gate_lines = "\n".join(
        (
            f"- {name}: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed
        in gates.items()
    )

    selection_lines = "\n".join(
        (
            f"- {float(row['relative_time_ps']):.3f} ps: "
            f"{row['selection_reasons']}; "
            f"ΔN_lumen="
            f"{float(row['mobile_minus_frozen_lumen_count']):.3f}; "
            f"max |ΔN_PYR,0.35|="
            f"{float(row['max_absolute_local_35pm_difference']):.3f}; "
            f"max |ΔN_PYR,0.50|="
            f"{float(row['max_absolute_local_50pm_difference']):.3f}"
        )
        for row in selected_rows
    )

    REPORT_MD.write_text(
        f"""# Matched Snapshot-Pair Review for a Limited QM Pilot

## Scope

- Representative matched pairs audited: **{len(review_rows)}**
- Snapshot time field: `{snapshot_time_field}`
- Strong local PYR-hydration contrast pairs:
  **{len(strong_local_rows)}**
- Strong lumen-occupancy contrast pairs:
  **{len(strong_lumen_rows)}**
- Low-contrast control pairs:
  **{len(control_rows)}**

The purpose of this review is screening only. Selection from the same
trajectory cannot support a publication-level causal claim and cannot
replace independent trajectory replicas.

## Screening definitions

A pair is considered locally informative when either:

- the maximum absolute PYR hydration-shell difference within the
  nominal 0.50 nm shell is at least 2 waters; or
- the maximum absolute PYR hydration-shell difference within the
  nominal 0.35 nm shell is at least 1 water.

A lumen contrast requires an absolute mobile-minus-frozen difference
of at least 2 waters.

A control-like pair requires:

- absolute lumen difference no greater than 1 water;
- maximum local 0.35 nm difference no greater than 1 water; and
- maximum local 0.50 nm difference no greater than 1 water.

## Selected candidate pairs

{selection_lines if selection_lines else '- NONE'}

## Pilot gates

{gate_lines}

## Decision

- Decision: **{decision}**
- Limited paired mobile-versus-frozen QM pilot authorized:
  **{'YES' if limited_qm_pilot_authorized else 'NO'}**
- Full electronic recalculation authorized: **NO**
- Longer mobile production authorized: **NO**
- Publication-level causal claim authorized: **NO**
- Authorized next step:
  `{authorized_next_step}`

The matched-water result remains a finding of mobility-associated
transient rehydration. It is not evidence of persistent confined-water
stabilization.
""",
        encoding="utf-8",
    )

    print(
        "Day022 matched snapshot-pair review completed."
    )

    print(
        "Representative pairs audited: "
        f"{len(review_rows)}/{EXPECTED_PAIRS}"
    )

    print(
        "Detected local hydration metrics: "
        f"{len(local_fields)}"
    )

    print(
        "Strong local / strong lumen / control-like pairs: "
        f"{len(strong_local_rows)}/"
        f"{len(strong_lumen_rows)}/"
        f"{len(control_rows)}"
    )

    print(
        "Selected pilot-pair count: "
        f"{len(selected_rows)}"
    )

    print(
        "Selected relative times: "
        + (
            " | ".join(
                f"{float(row['relative_time_ps']):.3f} ps"
                for row in selected_rows
            )
            if selected_rows
            else "NONE"
        )
    )

    print(
        "Selected local/informative/control counts: "
        f"{selected_local_count}/"
        f"{selected_informative_count}/"
        f"{selected_control_count}"
    )

    print(
        "Selected temporal quartiles / span: "
        f"{len(selected_quartiles)}/4 / "
        f"{temporal_span_ps:.3f} ps"
    )

    print(
        "Failed pilot gates: "
        + (
            " | ".join(
                failed_gates
            )
            if failed_gates
            else "NONE"
        )
    )

    print(
        f"Decision: {decision}"
    )

    print(
        "Limited mobile-vs-frozen QM pilot authorized: "
        f"{'YES' if limited_qm_pilot_authorized else 'NO'}"
    )

    print(
        "Full electronic recalculation authorized: NO"
    )

    print(
        "Longer mobile production authorized: NO"
    )

    print(
        "Publication-level causal claim authorized: NO"
    )

    print(
        "Authorized next step: "
        f"{authorized_next_step}"
    )

    print(
        "Wrote: "
        + str(
            REPORT_MD.relative_to(
                ROOT
            )
        )
    )


if __name__ == "__main__":
    main()
