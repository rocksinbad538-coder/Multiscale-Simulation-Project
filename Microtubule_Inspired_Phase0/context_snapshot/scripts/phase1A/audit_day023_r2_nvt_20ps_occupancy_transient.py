#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

DAY023_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design"
)

NVT_ROOT = (
    DAY023_ROOT
    / "16_r2_frozen_solute_nvt_20ps"
)

OUTPUT_ROOT = (
    DAY023_ROOT
    / "17_r2_nvt_20ps_occupancy_transient_audit"
)

SOURCE_SUMMARY = (
    NVT_ROOT
    / "r2_frozen_solute_nvt_20ps_summary.csv"
)

SOURCE_GATES = (
    NVT_ROOT
    / "r2_frozen_solute_nvt_20ps_gates.csv"
)

SOURCE_OCCUPANCY = (
    NVT_ROOT
    / "r2_frozen_solute_nvt_20ps_lumen_occupancy.csv"
)

WINDOWS_CSV = (
    OUTPUT_ROOT
    / "r2_nvt_20ps_occupancy_window_metrics.csv"
)

BLOCKS_CSV = (
    OUTPUT_ROOT
    / "r2_nvt_20ps_occupancy_block_metrics.csv"
)

AUDIT_SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_nvt_20ps_occupancy_transient_audit_summary.csv"
)

AUDIT_GATES_CSV = (
    OUTPUT_ROOT
    / "r2_nvt_20ps_occupancy_transient_audit_gates.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_NVT_20PS_OCCUPANCY_TRANSIENT_AUDIT_DAY023.md"
)

EXPECTED_SOURCE_DECISION = (
    "R2_FROZEN_SOLUTE_NVT_20PS_REQUIRES_REVIEW"
)

EXPECTED_FAILED_GATE = (
    "second_half_occupancy_slope_is_acceptable"
)

INITIAL_OCCUPANCY = 428

ORIGINAL_SLOPE_LIMIT_WATER_PS = 0.50

MINIMUM_OCCUPANCY_FRACTION = 0.80
MINIMUM_ENDPOINT_OCCUPANCY_FRACTION = 0.90
MINIMUM_ENDPOINT_IDENTITY_FRACTION = 0.90

MIN_CAP_OW_DISTANCE_NM = 0.15

FROZEN_TOLERANCE_NM = 1.0e-6

TEMPERATURE_TARGET_K = 300.0
TEMPERATURE_MEAN_TOLERANCE_K = 5.0
MAX_TEMPERATURE_SLOPE_K_PS = 0.10

MAX_CAP_SOL_LJ_KJ_MOL = 100.0

WINDOW_STARTS_PS = (
    0.0,
    5.0,
    10.0,
    12.5,
    15.0,
)

BLOCKS_PS = (
    (0.0, 5.0),
    (5.0, 10.0),
    (10.0, 15.0),
    (15.0, 20.0),
)


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
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


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


def read_csv_rows(
    path: Path,
) -> list[dict[str, str]]:
    require_file(path)

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(
            csv.DictReader(handle)
        )

    if not rows:
        raise RuntimeError(
            f"No rows found in {path}"
        )

    return rows


def read_single_csv_row(
    path: Path,
) -> dict[str, str]:
    rows = read_csv_rows(path)

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected one row in {path}; found {len(rows)}"
        )

    return rows[0]


def write_csv(
    path: Path,
    rows: list[dict[str, Any]],
) -> None:
    if not rows:
        raise RuntimeError(
            f"No rows available for {path}"
        )

    fields: list[str] = []

    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

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
                    key: row.get(key, "")
                    for key in fields
                }
            )


def as_float(
    row: dict[str, str],
    key: str,
) -> float:
    try:
        value = float(row[key])
    except (
        KeyError,
        TypeError,
        ValueError,
    ) as exc:
        raise RuntimeError(
            f"Could not parse numeric field {key!r}"
        ) from exc

    if not math.isfinite(value):
        raise RuntimeError(
            f"Non-finite value in field {key!r}"
        )

    return value


def as_int(
    row: dict[str, str],
    key: str,
) -> int:
    return int(
        round(
            as_float(
                row,
                key,
            )
        )
    )


def linear_slope(
    time: np.ndarray,
    values: np.ndarray,
) -> float:
    if len(time) < 2:
        return math.nan

    return float(
        np.polyfit(
            time,
            values,
            1,
        )[0]
    )


def describe_window(
    time: np.ndarray,
    occupancy: np.ndarray,
    retained: np.ndarray,
    noninitial: np.ndarray,
    start_ps: float,
) -> dict[str, Any]:
    mask = (
        time
        >= start_ps - 1.0e-9
    )

    window_time = time[mask]
    window_occupancy = occupancy[mask]
    window_retained = retained[mask]
    window_noninitial = noninitial[mask]

    if len(window_time) < 2:
        raise RuntimeError(
            f"Insufficient data for window {start_ps}-20 ps"
        )

    return {
        "window_start_ps": start_ps,
        "window_end_ps": float(
            window_time[-1]
        ),
        "points": len(window_time),
        "occupancy_start": int(
            window_occupancy[0]
        ),
        "occupancy_end": int(
            window_occupancy[-1]
        ),
        "occupancy_change": int(
            window_occupancy[-1]
            - window_occupancy[0]
        ),
        "occupancy_mean": float(
            np.mean(
                window_occupancy
            )
        ),
        "occupancy_minimum": int(
            np.min(
                window_occupancy
            )
        ),
        "occupancy_maximum": int(
            np.max(
                window_occupancy
            )
        ),
        "occupancy_slope_water_ps": (
            linear_slope(
                window_time,
                window_occupancy,
            )
        ),
        "retained_initial_start": int(
            window_retained[0]
        ),
        "retained_initial_end": int(
            window_retained[-1]
        ),
        "retained_initial_change": int(
            window_retained[-1]
            - window_retained[0]
        ),
        "retained_initial_slope_water_ps": (
            linear_slope(
                window_time,
                window_retained,
            )
        ),
        "noninitial_occupancy_start": int(
            window_noninitial[0]
        ),
        "noninitial_occupancy_end": int(
            window_noninitial[-1]
        ),
        "noninitial_occupancy_maximum": int(
            np.max(
                window_noninitial
            )
        ),
    }


def describe_block(
    time: np.ndarray,
    occupancy: np.ndarray,
    retained: np.ndarray,
    noninitial: np.ndarray,
    start_ps: float,
    end_ps: float,
    *,
    include_end: bool,
) -> dict[str, Any]:
    if include_end:
        mask = (
            (time >= start_ps - 1.0e-9)
            & (time <= end_ps + 1.0e-9)
        )
    else:
        mask = (
            (time >= start_ps - 1.0e-9)
            & (time < end_ps - 1.0e-9)
        )

    block_time = time[mask]
    block_occupancy = occupancy[mask]
    block_retained = retained[mask]
    block_noninitial = noninitial[mask]

    if len(block_time) < 2:
        raise RuntimeError(
            f"Insufficient data for block {start_ps}-{end_ps} ps"
        )

    return {
        "block_start_ps": start_ps,
        "block_end_ps": end_ps,
        "points": len(block_time),
        "occupancy_start": int(
            block_occupancy[0]
        ),
        "occupancy_end": int(
            block_occupancy[-1]
        ),
        "occupancy_change": int(
            block_occupancy[-1]
            - block_occupancy[0]
        ),
        "occupancy_mean": float(
            np.mean(
                block_occupancy
            )
        ),
        "occupancy_minimum": int(
            np.min(
                block_occupancy
            )
        ),
        "occupancy_maximum": int(
            np.max(
                block_occupancy
            )
        ),
        "occupancy_slope_water_ps": (
            linear_slope(
                block_time,
                block_occupancy,
            )
        ),
        "retained_initial_start": int(
            block_retained[0]
        ),
        "retained_initial_end": int(
            block_retained[-1]
        ),
        "retained_initial_change": int(
            block_retained[-1]
            - block_retained[0]
        ),
        "noninitial_occupancy_mean": float(
            np.mean(
                block_noninitial
            )
        ),
        "noninitial_occupancy_maximum": int(
            np.max(
                block_noninitial
            )
        ),
    }


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        SOURCE_SUMMARY,
        SOURCE_GATES,
        SOURCE_OCCUPANCY,
    ):
        require_file(required)

    summary = read_single_csv_row(
        SOURCE_SUMMARY
    )

    source_gate_rows = read_csv_rows(
        SOURCE_GATES
    )

    occupancy_rows = read_csv_rows(
        SOURCE_OCCUPANCY
    )

    time = np.asarray(
        [
            as_float(
                row,
                "time_ps",
            )
            for row in occupancy_rows
        ],
        dtype=float,
    )

    occupancy = np.asarray(
        [
            as_int(
                row,
                "lumen_occupancy",
            )
            for row in occupancy_rows
        ],
        dtype=int,
    )

    retained = np.asarray(
        [
            as_int(
                row,
                "initial_lumen_waters_retained",
            )
            for row in occupancy_rows
        ],
        dtype=int,
    )

    cap_distance = np.asarray(
        [
            as_float(
                row,
                "minimum_CAP_OW_distance_nm",
            )
            for row in occupancy_rows
        ],
        dtype=float,
    )

    if not (
        len(time)
        == len(occupancy)
        == len(retained)
        == len(cap_distance)
        == 41
    ):
        raise RuntimeError(
            "Expected 41 matched trajectory records."
        )

    if not np.allclose(
        np.diff(time),
        0.5,
        rtol=0.0,
        atol=1.0e-6,
    ):
        raise RuntimeError(
            "Unexpected trajectory time grid."
        )

    if not math.isclose(
        float(time[0]),
        0.0,
        rel_tol=0.0,
        abs_tol=1.0e-8,
    ):
        raise RuntimeError(
            "Trajectory does not begin at 0 ps."
        )

    if not math.isclose(
        float(time[-1]),
        20.0,
        rel_tol=0.0,
        abs_tol=1.0e-8,
    ):
        raise RuntimeError(
            "Trajectory does not end at 20 ps."
        )

    noninitial = (
        occupancy
        - retained
    )

    if np.any(noninitial < 0):
        raise RuntimeError(
            "Retained initial identities exceed occupancy."
        )

    failed_source_gates = [
        row.get(
            "gate",
            "",
        )
        for row in source_gate_rows
        if not parse_bool(
            row.get(
                "pass",
                "false",
            )
        )
    ]

    passing_source_gates = [
        row.get(
            "gate",
            "",
        )
        for row in source_gate_rows
        if parse_bool(
            row.get(
                "pass",
                "false",
            )
        )
    ]

    window_rows = [
        describe_window(
            time,
            occupancy,
            retained,
            noninitial,
            start_ps,
        )
        for start_ps in WINDOW_STARTS_PS
    ]

    block_rows = []

    for block_index, (
        start_ps,
        end_ps,
    ) in enumerate(BLOCKS_PS):
        block_rows.append(
            describe_block(
                time,
                occupancy,
                retained,
                noninitial,
                start_ps,
                end_ps,
                include_end=(
                    block_index
                    == len(BLOCKS_PS) - 1
                ),
            )
        )

    write_csv(
        WINDOWS_CSV,
        window_rows,
    )

    write_csv(
        BLOCKS_CSV,
        block_rows,
    )

    initial = int(
        occupancy[0]
    )

    endpoint = int(
        occupancy[-1]
    )

    minimum = int(
        np.min(
            occupancy
        )
    )

    endpoint_retained = int(
        retained[-1]
    )

    endpoint_noninitial = int(
        noninitial[-1]
    )

    maximum_noninitial = int(
        np.max(
            noninitial
        )
    )

    second_half_slope = as_float(
        summary,
        "second_half_occupancy_slope_water_ps",
    )

    slope_excess = max(
        0.0,
        abs(
            second_half_slope
        )
        - ORIGINAL_SLOPE_LIMIT_WATER_PS,
    )

    projected_50ps_occupancy = (
        endpoint
        + second_half_slope
        * 30.0
    )

    temperature_post5_mean = as_float(
        summary,
        "temperature_post5_mean_K",
    )

    temperature_post5_std = as_float(
        summary,
        "temperature_post5_std_K",
    )

    temperature_post5_slope = as_float(
        summary,
        "temperature_post5_slope_K_ps",
    )

    theoretical_temperature_std = as_float(
        summary,
        "theoretical_temperature_std_K",
    )

    minimum_cap_distance = float(
        np.min(
            cap_distance
        )
    )

    cap_sol_lj_maximum = as_float(
        summary,
        "CAP_SOL_LJ_maximum_kJ_mol",
    )

    hbn_displacement = as_float(
        summary,
        "HBN_maximum_displacement_nm",
    )

    pyr_displacement = as_float(
        summary,
        "PYR_maximum_displacement_nm",
    )

    caps_displacement = as_float(
        summary,
        "CAPS_maximum_displacement_nm",
    )

    audit_gates = {
        "source_decision_requires_review": (
            summary.get(
                "decision"
            )
            == EXPECTED_SOURCE_DECISION
        ),
        "exactly_one_source_gate_failed": (
            len(
                failed_source_gates
            )
            == 1
        ),
        "only_stationarity_slope_gate_failed": (
            failed_source_gates
            == [
                EXPECTED_FAILED_GATE
            ]
        ),
        "all_other_source_gates_passed": (
            len(
                passing_source_gates
            )
            == len(
                source_gate_rows
            )
            - 1
        ),
        "mdrun_completed_successfully": (
            as_int(
                summary,
                "mdrun_return_code",
            )
            == 0
        ),
        "trajectory_check_completed_successfully": (
            as_int(
                summary,
                "trajectory_check_return_code",
            )
            == 0
        ),
        "checkpoint_was_written": (
            parse_bool(
                summary.get(
                    "checkpoint_written",
                    "false",
                )
            )
        ),
        "no_instability_signatures": (
            as_int(
                summary,
                "instability_signature_count",
            )
            == 0
        ),
        "trajectory_contains_41_frames": (
            len(time) == 41
        ),
        "initial_occupancy_is_428": (
            initial
            == INITIAL_OCCUPANCY
        ),
        "minimum_occupancy_remains_above_80_percent": (
            minimum
            >= (
                MINIMUM_OCCUPANCY_FRACTION
                * INITIAL_OCCUPANCY
            )
        ),
        "endpoint_occupancy_remains_above_90_percent": (
            endpoint
            >= (
                MINIMUM_ENDPOINT_OCCUPANCY_FRACTION
                * INITIAL_OCCUPANCY
            )
        ),
        "endpoint_identity_retention_remains_above_90_percent": (
            endpoint_retained
            >= (
                MINIMUM_ENDPOINT_IDENTITY_FRACTION
                * INITIAL_OCCUPANCY
            )
        ),
        "noninitial_lumen_waters_demonstrate_exchange": (
            maximum_noninitial
            > 0
        ),
        "minimum_CAP_OW_distance_is_safe": (
            minimum_cap_distance
            >= MIN_CAP_OW_DISTANCE_NM
        ),
        "post5_temperature_mean_is_stable": (
            abs(
                temperature_post5_mean
                - TEMPERATURE_TARGET_K
            )
            <= TEMPERATURE_MEAN_TOLERANCE_K
        ),
        "post5_temperature_std_is_canonical": (
            temperature_post5_std
            <= (
                3.0
                * theoretical_temperature_std
            )
        ),
        "post5_temperature_slope_is_small": (
            abs(
                temperature_post5_slope
            )
            <= MAX_TEMPERATURE_SLOPE_K_PS
        ),
        "HBN_remained_frozen": (
            hbn_displacement
            <= FROZEN_TOLERANCE_NM
        ),
        "PYR_remained_frozen": (
            pyr_displacement
            <= FROZEN_TOLERANCE_NM
        ),
        "CAPS_remained_frozen": (
            caps_displacement
            <= FROZEN_TOLERANCE_NM
        ),
        "CAP_SOL_LJ_remained_below_100kJmol": (
            cap_sol_lj_maximum
            <= MAX_CAP_SOL_LJ_KJ_MOL
        ),
    }

    failed_audit_gates = [
        name
        for name, passed
        in audit_gates.items()
        if not passed
    ]

    extension_justified = (
        len(
            failed_audit_gates
        )
        == 0
    )

    decision = (
        "R2_OCCUPANCY_TRANSIENT_CHECKPOINT_EXTENSION_JUSTIFIED"
        if extension_justified
        else
        "R2_OCCUPANCY_TRANSIENT_AUDIT_REQUIRES_REVIEW"
    )

    required_next_step = (
        "PREPARE_R2_30PS_CHECKPOINT_CONTINUATION_TO_50PS"
        if extension_justified
        else
        "REVIEW_R2_OCCUPANCY_TRANSIENT_AUDIT_FAILURES"
    )

    audit_summary = {
        "decision": decision,
        "source_decision": (
            summary.get(
                "decision",
                "",
            )
        ),
        "source_gate_count": (
            len(
                source_gate_rows
            )
        ),
        "source_failed_gate_count": (
            len(
                failed_source_gates
            )
        ),
        "source_failed_gates": (
            " | ".join(
                failed_source_gates
            )
        ),
        "initial_occupancy": (
            initial
        ),
        "minimum_occupancy": (
            minimum
        ),
        "endpoint_occupancy": (
            endpoint
        ),
        "endpoint_occupancy_fraction": (
            endpoint
            / initial
        ),
        "endpoint_initial_identities_retained": (
            endpoint_retained
        ),
        "endpoint_initial_identity_fraction": (
            endpoint_retained
            / initial
        ),
        "endpoint_noninitial_lumen_waters": (
            endpoint_noninitial
        ),
        "maximum_noninitial_lumen_waters": (
            maximum_noninitial
        ),
        "second_half_occupancy_slope_water_ps": (
            second_half_slope
        ),
        "original_slope_limit_water_ps": (
            ORIGINAL_SLOPE_LIMIT_WATER_PS
        ),
        "slope_limit_excess_water_ps": (
            slope_excess
        ),
        "linear_50ps_projection_from_second_half_slope": (
            projected_50ps_occupancy
        ),
        "minimum_CAP_OW_distance_nm": (
            minimum_cap_distance
        ),
        "temperature_post5_mean_K": (
            temperature_post5_mean
        ),
        "temperature_post5_std_K": (
            temperature_post5_std
        ),
        "temperature_post5_slope_K_ps": (
            temperature_post5_slope
        ),
        "failed_audit_gates": (
            " | ".join(
                failed_audit_gates
            )
        ),
        "20ps_R2_validation_status": (
            "NOT_VALIDATED_NONSTATIONARY"
        ),
        "checkpoint_continuation_preparation_authorized": (
            extension_justified
        ),
        "checkpoint_continuation_execution_authorized": False,
        "long_mobile_MD_authorized": False,
        "multitemperature_MD_authorized": False,
        "QM_recalculation_authorized": False,
        "required_next_step": (
            required_next_step
        ),
    }

    write_csv(
        AUDIT_SUMMARY_CSV,
        [
            audit_summary
        ],
    )

    write_csv(
        AUDIT_GATES_CSV,
        [
            {
                "gate": name,
                "pass": passed,
            }
            for name, passed
            in audit_gates.items()
        ],
    )

    window_lines = "\n".join(
        (
            f"- {row['window_start_ps']:.1f}–"
            f"{row['window_end_ps']:.1f} ps: "
            f"mean={row['occupancy_mean']:.4f}; "
            f"min/max={row['occupancy_minimum']}/"
            f"{row['occupancy_maximum']}; "
            f"change={row['occupancy_change']:+d}; "
            f"slope="
            f"{row['occupancy_slope_water_ps']:.6f} waters/ps"
        )
        for row in window_rows
    )

    block_lines = "\n".join(
        (
            f"- {row['block_start_ps']:.1f}–"
            f"{row['block_end_ps']:.1f} ps: "
            f"mean={row['occupancy_mean']:.4f}; "
            f"change={row['occupancy_change']:+d}; "
            f"slope="
            f"{row['occupancy_slope_water_ps']:.6f} waters/ps; "
            f"maximum noninitial occupancy="
            f"{row['noninitial_occupancy_maximum']}"
        )
        for row in block_rows
    )

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed
        in audit_gates.items()
    )

    REPORT_MD.write_text(
        f"""# R2 20 ps Occupancy-Transient Audit

## Scope

This audit evaluates the completed R2 20 ps frozen-solute trajectory.
It does not repeat molecular dynamics and does not change the original
stationarity threshold.

The original 20 ps gate remains unpassed because the second-half
occupancy slope exceeded the predefined absolute limit of
{ORIGINAL_SLOPE_LIMIT_WATER_PS:.2f} waters/ps.

## Source result

- Source decision:
  **{summary.get('decision', '')}**
- Failed source gates:
  **{'NONE' if not failed_source_gates else ' | '.join(failed_source_gates)}**
- Initial/minimum/endpoint occupancy:
  **{initial}/{minimum}/{endpoint}**
- Endpoint occupancy fraction:
  **{endpoint / initial:.6f}**
- Endpoint retained initial identities:
  **{endpoint_retained}/{initial}**
- Endpoint identity-retention fraction:
  **{endpoint_retained / initial:.6f}**
- Endpoint noninitial luminal waters:
  **{endpoint_noninitial}**
- Maximum noninitial luminal waters:
  **{maximum_noninitial}**

## Stationarity metric

- Second-half occupancy slope:
  **{second_half_slope:.6f} waters/ps**
- Original absolute limit:
  **{ORIGINAL_SLOPE_LIMIT_WATER_PS:.6f} waters/ps**
- Excess beyond limit:
  **{slope_excess:.6f} waters/ps**

The slope is not reclassified as passing. It is used only to determine
whether a checkpoint extension is scientifically justified as a
stationarity diagnostic.

## Cumulative windows

{window_lines}

## Five-picosecond blocks

{block_lines}

## Safety and execution gates

{gate_lines}

## Decision

- Audit decision:
  **{decision}**
- R2 20 ps validation status:
  **NOT VALIDATED — NONSTATIONARY**
- MD rerun required:
  **NO**
- Checkpoint-continuation preparation authorized:
  **{'YES' if extension_justified else 'NO'}**
- Checkpoint-continuation execution authorized:
  **NO**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `{required_next_step}`

The proposed extension must use the exact 20 ps checkpoint and must not
regenerate velocities. Its purpose is to determine whether occupancy
approaches a plateau or continues to decline.
""",
        encoding="utf-8",
    )

    print(
        "Day023 R2 20 ps occupancy-transient "
        "audit completed."
    )

    print(
        "Source decision: "
        f"{summary.get('decision', '')}"
    )

    print(
        "Source gates / failed: "
        f"{len(source_gate_rows)}/"
        f"{len(failed_source_gates)}"
    )

    print(
        "Failed source gates: "
        + (
            "NONE"
            if not failed_source_gates
            else " | ".join(
                failed_source_gates
            )
        )
    )

    print(
        "Initial / minimum / endpoint occupancy: "
        f"{initial}/{minimum}/{endpoint}"
    )

    print(
        "Endpoint occupancy fraction: "
        f"{endpoint / initial:.6f}"
    )

    print(
        "Endpoint retained initial identities: "
        f"{endpoint_retained}/{initial} "
        f"({endpoint_retained / initial:.6f})"
    )

    print(
        "Endpoint / maximum noninitial luminal waters: "
        f"{endpoint_noninitial}/"
        f"{maximum_noninitial}"
    )

    print(
        "Second-half slope / limit / excess: "
        f"{second_half_slope:.6f}/"
        f"{ORIGINAL_SLOPE_LIMIT_WATER_PS:.6f}/"
        f"{slope_excess:.6f} waters/ps"
    )

    for row in window_rows:
        print(
            "Window "
            f"{row['window_start_ps']:.1f}-"
            f"{row['window_end_ps']:.1f} ps | "
            f"mean/min/max="
            f"{row['occupancy_mean']:.4f}/"
            f"{row['occupancy_minimum']}/"
            f"{row['occupancy_maximum']} | "
            f"change="
            f"{row['occupancy_change']:+d} | "
            f"slope="
            f"{row['occupancy_slope_water_ps']:.6f} waters/ps"
        )

    for row in block_rows:
        print(
            "Block "
            f"{row['block_start_ps']:.1f}-"
            f"{row['block_end_ps']:.1f} ps | "
            f"mean="
            f"{row['occupancy_mean']:.4f} | "
            f"change="
            f"{row['occupancy_change']:+d} | "
            f"slope="
            f"{row['occupancy_slope_water_ps']:.6f} waters/ps | "
            f"noninitial max="
            f"{row['noninitial_occupancy_maximum']}"
        )

    print(
        "Minimum CAP-OW distance: "
        f"{minimum_cap_distance:.6f} nm"
    )

    print(
        "Post-5ps temperature mean/std/slope: "
        f"{temperature_post5_mean:.4f}/"
        f"{temperature_post5_std:.4f}/"
        f"{temperature_post5_slope:.6f} "
        "K / K / K ps^-1"
    )

    print(
        f"Decision: {decision}"
    )

    print(
        "Failed audit gates: "
        + (
            "NONE"
            if not failed_audit_gates
            else " | ".join(
                failed_audit_gates
            )
        )
    )

    print(
        "R2 20 ps validation status: "
        "NOT VALIDATED - NONSTATIONARY"
    )

    print(
        "MD rerun required: NO"
    )

    print(
        "Checkpoint-continuation preparation authorized: "
        f"{'YES' if extension_justified else 'NO'}"
    )

    print(
        "Checkpoint-continuation execution authorized: NO"
    )

    print(
        "Long mobile MD authorized: NO"
    )

    print(
        "Multitemperature MD authorized: NO"
    )

    print(
        "QM recalculation authorized: NO"
    )

    print(
        f"Required next step: {required_next_step}"
    )

    for path in (
        WINDOWS_CSV,
        BLOCKS_CSV,
        AUDIT_SUMMARY_CSV,
        AUDIT_GATES_CSV,
        REPORT_MD,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not extension_justified:
        raise RuntimeError(
            "R2 occupancy-transient audit requires review."
        )


if __name__ == "__main__":
    main()
