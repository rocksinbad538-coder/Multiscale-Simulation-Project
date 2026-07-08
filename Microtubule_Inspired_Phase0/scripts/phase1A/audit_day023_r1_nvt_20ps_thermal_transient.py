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
    / "08_r1_frozen_solute_nvt_20ps"
)

PREPARATION_ROOT = (
    DAY023_ROOT
    / "07_r1_frozen_solute_nvt_20ps_preparation"
)

OUTPUT_ROOT = (
    DAY023_ROOT
    / "09_r1_nvt_20ps_thermal_review"
)

ENERGY_CSV = (
    NVT_ROOT
    / "r1_frozen_solute_nvt_20ps_energy_series.csv"
)

OCCUPANCY_CSV = (
    NVT_ROOT
    / "r1_frozen_solute_nvt_20ps_lumen_occupancy.csv"
)

NVT_SUMMARY_CSV = (
    NVT_ROOT
    / "r1_frozen_solute_nvt_20ps_summary.csv"
)

WARNING_SUMMARY_CSV = (
    PREPARATION_ROOT
    / "r1_nvt_grompp_warning_authorization_summary.csv"
)

WINDOW_CSV = (
    OUTPUT_ROOT
    / "r1_nvt_20ps_thermal_window_statistics.csv"
)

INITIAL_SAMPLES_CSV = (
    OUTPUT_ROOT
    / "r1_nvt_20ps_initial_temperature_samples.csv"
)

GATE_CSV = (
    OUTPUT_ROOT
    / "r1_nvt_20ps_thermal_review_gates.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r1_nvt_20ps_thermal_review_summary.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R1_NVT_20PS_THERMAL_TRANSIENT_REVIEW_DAY023.md"
)

WINDOW_STARTS_PS = (
    0.0,
    0.05,
    0.10,
    0.25,
    0.50,
    1.0,
    2.0,
    5.0,
    10.0,
)

TARGET_TEMPERATURE_K = 300.0

MIN_STABLE_MEAN_TEMPERATURE_K = 295.0
MAX_STABLE_MEAN_TEMPERATURE_K = 305.0

MAX_POST_5PS_TEMPERATURE_STD_K = 5.0
MIN_POST_5PS_TEMPERATURE_K = 280.0
MAX_POST_5PS_TEMPERATURE_K = 320.0

MAX_LAST_10PS_TEMPERATURE_SLOPE_K_PER_PS = 0.20

MIN_LUMEN_OCCUPANCY = 386
MIN_CAP_OW_DISTANCE_NM = 0.15

EXPECTED_INITIAL_LUMEN_WATERS = 428
EXPECTED_DURATION_PS = 20.0


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
            f"No data rows in {path}"
        )

    return rows


def read_single_csv_row(
    path: Path,
) -> dict[str, str]:
    rows = read_csv_rows(path)

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected one row in {path}; "
            f"found {len(rows)}"
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


def numeric_column(
    rows: list[dict[str, str]],
    name: str,
) -> np.ndarray:
    if name not in rows[0]:
        raise RuntimeError(
            f"Column {name!r} is absent."
        )

    values = np.asarray(
        [
            float(row[name])
            for row in rows
        ],
        dtype=float,
    )

    if not np.all(
        np.isfinite(values)
    ):
        raise RuntimeError(
            f"Column {name!r} contains non-finite values."
        )

    return values


def linear_slope(
    times: np.ndarray,
    values: np.ndarray,
) -> float:
    if len(times) < 2:
        return math.nan

    coefficients = np.polyfit(
        times,
        values,
        1,
    )

    return float(
        coefficients[0]
    )


def window_statistics(
    times: np.ndarray,
    temperature: np.ndarray,
    potential: np.ndarray,
    total_energy: np.ndarray,
    cap_sol_lj: np.ndarray,
    start_ps: float,
) -> dict[str, Any]:
    mask = times >= (
        start_ps - 1.0e-12
    )

    if np.count_nonzero(mask) < 2:
        raise RuntimeError(
            f"Insufficient data for window {start_ps} ps."
        )

    selected_times = times[mask]
    selected_temperature = temperature[mask]
    selected_potential = potential[mask]
    selected_total_energy = total_energy[mask]
    selected_cap_sol_lj = cap_sol_lj[mask]

    return {
        "window_start_ps": start_ps,
        "window_end_ps": float(
            selected_times[-1]
        ),
        "points": int(
            len(selected_times)
        ),
        "temperature_mean_K": float(
            np.mean(
                selected_temperature
            )
        ),
        "temperature_std_K": float(
            np.std(
                selected_temperature,
                ddof=1,
            )
        ),
        "temperature_min_K": float(
            np.min(
                selected_temperature
            )
        ),
        "temperature_max_K": float(
            np.max(
                selected_temperature
            )
        ),
        "temperature_slope_K_per_ps": (
            linear_slope(
                selected_times,
                selected_temperature,
            )
        ),
        "potential_initial_kJ_mol": float(
            selected_potential[0]
        ),
        "potential_final_kJ_mol": float(
            selected_potential[-1]
        ),
        "potential_change_kJ_mol": float(
            selected_potential[-1]
            - selected_potential[0]
        ),
        "potential_slope_kJ_mol_ps": (
            linear_slope(
                selected_times,
                selected_potential,
            )
        ),
        "total_energy_initial_kJ_mol": float(
            selected_total_energy[0]
        ),
        "total_energy_final_kJ_mol": float(
            selected_total_energy[-1]
        ),
        "total_energy_change_kJ_mol": float(
            selected_total_energy[-1]
            - selected_total_energy[0]
        ),
        "total_energy_slope_kJ_mol_ps": (
            linear_slope(
                selected_times,
                selected_total_energy,
            )
        ),
        "CAP_SOL_LJ_mean_kJ_mol": float(
            np.mean(
                selected_cap_sol_lj
            )
        ),
        "CAP_SOL_LJ_min_kJ_mol": float(
            np.min(
                selected_cap_sol_lj
            )
        ),
        "CAP_SOL_LJ_max_kJ_mol": float(
            np.max(
                selected_cap_sol_lj
            )
        ),
        "CAP_SOL_LJ_final_kJ_mol": float(
            selected_cap_sol_lj[-1]
        ),
    }


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        ENERGY_CSV,
        OCCUPANCY_CSV,
        NVT_SUMMARY_CSV,
        WARNING_SUMMARY_CSV,
    ):
        require_file(required)

    energy_rows = read_csv_rows(
        ENERGY_CSV
    )

    occupancy_rows = read_csv_rows(
        OCCUPANCY_CSV
    )

    nvt_summary = read_single_csv_row(
        NVT_SUMMARY_CSV
    )

    warning_summary = read_single_csv_row(
        WARNING_SUMMARY_CSV
    )

    times = numeric_column(
        energy_rows,
        "time_ps",
    )

    temperature = numeric_column(
        energy_rows,
        "temperature_K",
    )

    potential = numeric_column(
        energy_rows,
        "potential_kJ_mol",
    )

    total_energy = numeric_column(
        energy_rows,
        "total_energy_kJ_mol",
    )

    cap_sol_lj = numeric_column(
        energy_rows,
        "CAP_SOL_LJ_kJ_mol",
    )

    if len(times) < 10:
        raise RuntimeError(
            "The energy series is unexpectedly short."
        )

    if not np.all(
        np.diff(times) > 0.0
    ):
        raise RuntimeError(
            "Energy times are not strictly increasing."
        )

    if not math.isclose(
        float(times[0]),
        0.0,
        rel_tol=0.0,
        abs_tol=1.0e-10,
    ):
        raise RuntimeError(
            "The energy series does not start at 0 ps."
        )

    if not math.isclose(
        float(times[-1]),
        EXPECTED_DURATION_PS,
        rel_tol=0.0,
        abs_tol=1.0e-8,
    ):
        raise RuntimeError(
            "The energy series does not end at 20 ps."
        )

    window_rows = [
        window_statistics(
            times,
            temperature,
            potential,
            total_energy,
            cap_sol_lj,
            start_ps,
        )
        for start_ps in WINDOW_STARTS_PS
    ]

    write_csv(
        WINDOW_CSV,
        window_rows,
    )

    initial_sample_count = min(
        41,
        len(times),
    )

    initial_rows = [
        {
            "index": index,
            "time_ps": float(
                times[index]
            ),
            "temperature_K": float(
                temperature[index]
            ),
            "potential_kJ_mol": float(
                potential[index]
            ),
            "total_energy_kJ_mol": float(
                total_energy[index]
            ),
            "CAP_SOL_LJ_kJ_mol": float(
                cap_sol_lj[index]
            ),
        }
        for index in range(
            initial_sample_count
        )
    ]

    write_csv(
        INITIAL_SAMPLES_CSV,
        initial_rows,
    )

    window_by_start = {
        float(
            row[
                "window_start_ps"
            ]
        ): row
        for row in window_rows
    }

    post_5ps = window_by_start[
        5.0
    ]

    last_10ps = window_by_start[
        10.0
    ]

    degrees_of_freedom = float(
        warning_summary[
            "SOL_degrees_of_freedom"
        ]
    )

    theoretical_temperature_std_K = (
        TARGET_TEMPERATURE_K
        * math.sqrt(
            2.0
            / degrees_of_freedom
        )
    )

    all_temperature_std_K = float(
        np.std(
            temperature,
            ddof=1,
        )
    )

    without_first_temperature_std_K = float(
        np.std(
            temperature[1:],
            ddof=1,
        )
    )

    first_temperature_K = float(
        temperature[0]
    )

    second_temperature_K = float(
        temperature[1]
    )

    first_sample_deviation_K = (
        first_temperature_K
        - TARGET_TEMPERATURE_K
    )

    occupancy_times = numeric_column(
        occupancy_rows,
        "time_ps",
    )

    occupancy = numeric_column(
        occupancy_rows,
        "lumen_occupancy",
    )

    retained_initial = numeric_column(
        occupancy_rows,
        "initial_lumen_waters_retained",
    )

    cap_distance = numeric_column(
        occupancy_rows,
        "minimum_CAP_OW_distance_nm",
    )

    occupancy_post_5_mask = (
        occupancy_times
        >= 5.0 - 1.0e-12
    )

    occupancy_post_5 = occupancy[
        occupancy_post_5_mask
    ]

    retained_post_5 = retained_initial[
        occupancy_post_5_mask
    ]

    cap_distance_post_5 = cap_distance[
        occupancy_post_5_mask
    ]

    original_failed_gates = [
        item.strip()
        for item in nvt_summary.get(
            "failed_gates",
            "",
        ).split("|")
        if item.strip()
    ]

    gates = {
        "original_mdrun_return_code_zero": (
            int(
                nvt_summary[
                    "mdrun_return_code"
                ]
            )
            == 0
        ),
        "original_mdrun_finished": (
            parse_bool(
                nvt_summary[
                    "mdrun_finished"
                ]
            )
        ),
        "original_failure_is_temperature_std_only": (
            original_failed_gates
            == [
                "temperature_standard_deviation_is_acceptable"
            ]
        ),
        "post_5ps_temperature_mean_is_295_to_305K": (
            MIN_STABLE_MEAN_TEMPERATURE_K
            <= float(
                post_5ps[
                    "temperature_mean_K"
                ]
            )
            <= MAX_STABLE_MEAN_TEMPERATURE_K
        ),
        "post_5ps_temperature_std_is_at_most_5K": (
            float(
                post_5ps[
                    "temperature_std_K"
                ]
            )
            <= MAX_POST_5PS_TEMPERATURE_STD_K
        ),
        "post_5ps_temperature_min_is_at_least_280K": (
            float(
                post_5ps[
                    "temperature_min_K"
                ]
            )
            >= MIN_POST_5PS_TEMPERATURE_K
        ),
        "post_5ps_temperature_max_is_at_most_320K": (
            float(
                post_5ps[
                    "temperature_max_K"
                ]
            )
            <= MAX_POST_5PS_TEMPERATURE_K
        ),
        "last_10ps_temperature_mean_is_295_to_305K": (
            MIN_STABLE_MEAN_TEMPERATURE_K
            <= float(
                last_10ps[
                    "temperature_mean_K"
                ]
            )
            <= MAX_STABLE_MEAN_TEMPERATURE_K
        ),
        "last_10ps_temperature_std_is_at_most_5K": (
            float(
                last_10ps[
                    "temperature_std_K"
                ]
            )
            <= MAX_POST_5PS_TEMPERATURE_STD_K
        ),
        "last_10ps_temperature_slope_is_small": (
            abs(
                float(
                    last_10ps[
                        "temperature_slope_K_per_ps"
                    ]
                )
            )
            <= MAX_LAST_10PS_TEMPERATURE_SLOPE_K_PER_PS
        ),
        "post_5ps_lumen_occupancy_remains_above_90_percent": (
            int(
                np.min(
                    occupancy_post_5
                )
            )
            >= MIN_LUMEN_OCCUPANCY
        ),
        "post_5ps_initial_lumen_retention_remains_above_90_percent": (
            int(
                np.min(
                    retained_post_5
                )
            )
            >= MIN_LUMEN_OCCUPANCY
        ),
        "post_5ps_CAP_OW_distance_remains_above_limit": (
            float(
                np.min(
                    cap_distance_post_5
                )
            )
            >= MIN_CAP_OW_DISTANCE_NM
        ),
        "endpoint_lumen_occupancy_is_428": (
            int(
                occupancy[-1]
            )
            == EXPECTED_INITIAL_LUMEN_WATERS
        ),
        "endpoint_initial_lumen_retention_is_428": (
            int(
                retained_initial[-1]
            )
            == EXPECTED_INITIAL_LUMEN_WATERS
        ),
    }

    failed_gates = [
        name
        for name, passed
        in gates.items()
        if not passed
    ]

    transient_only = (
        len(failed_gates) == 0
    )

    decision = (
        "R1_INITIAL_THERMALIZATION_TRANSIENT_CONFIRMED"
        if transient_only
        else
        "R1_THERMAL_STABILITY_REQUIRES_CORRECTION"
    )

    required_next_step = (
        "PREPARE_R1_30PS_CHECKPOINT_CONTINUATION_TO_50PS"
        if transient_only
        else
        "REVIEW_R1_POST_EQUILIBRATION_THERMAL_METRICS"
    )

    gate_rows = [
        {
            "gate": name,
            "pass": passed,
        }
        for name, passed in gates.items()
    ]

    write_csv(
        GATE_CSV,
        gate_rows,
    )

    summary = {
        "decision": decision,
        "energy_points": len(times),
        "first_temperature_K": (
            first_temperature_K
        ),
        "second_temperature_K": (
            second_temperature_K
        ),
        "first_sample_deviation_from_300K": (
            first_sample_deviation_K
        ),
        "all_window_temperature_mean_K": float(
            np.mean(
                temperature
            )
        ),
        "all_window_temperature_std_K": (
            all_temperature_std_K
        ),
        "temperature_std_without_first_sample_K": (
            without_first_temperature_std_K
        ),
        "theoretical_canonical_temperature_std_K": (
            theoretical_temperature_std_K
        ),
        "post_5ps_temperature_mean_K": (
            post_5ps[
                "temperature_mean_K"
            ]
        ),
        "post_5ps_temperature_std_K": (
            post_5ps[
                "temperature_std_K"
            ]
        ),
        "post_5ps_temperature_min_K": (
            post_5ps[
                "temperature_min_K"
            ]
        ),
        "post_5ps_temperature_max_K": (
            post_5ps[
                "temperature_max_K"
            ]
        ),
        "post_5ps_temperature_slope_K_per_ps": (
            post_5ps[
                "temperature_slope_K_per_ps"
            ]
        ),
        "last_10ps_temperature_mean_K": (
            last_10ps[
                "temperature_mean_K"
            ]
        ),
        "last_10ps_temperature_std_K": (
            last_10ps[
                "temperature_std_K"
            ]
        ),
        "last_10ps_temperature_min_K": (
            last_10ps[
                "temperature_min_K"
            ]
        ),
        "last_10ps_temperature_max_K": (
            last_10ps[
                "temperature_max_K"
            ]
        ),
        "last_10ps_temperature_slope_K_per_ps": (
            last_10ps[
                "temperature_slope_K_per_ps"
            ]
        ),
        "post_5ps_potential_change_kJ_mol": (
            post_5ps[
                "potential_change_kJ_mol"
            ]
        ),
        "post_5ps_potential_slope_kJ_mol_ps": (
            post_5ps[
                "potential_slope_kJ_mol_ps"
            ]
        ),
        "last_10ps_potential_change_kJ_mol": (
            last_10ps[
                "potential_change_kJ_mol"
            ]
        ),
        "last_10ps_potential_slope_kJ_mol_ps": (
            last_10ps[
                "potential_slope_kJ_mol_ps"
            ]
        ),
        "post_5ps_total_energy_change_kJ_mol": (
            post_5ps[
                "total_energy_change_kJ_mol"
            ]
        ),
        "last_10ps_total_energy_change_kJ_mol": (
            last_10ps[
                "total_energy_change_kJ_mol"
            ]
        ),
        "post_5ps_CAP_SOL_LJ_mean_kJ_mol": (
            post_5ps[
                "CAP_SOL_LJ_mean_kJ_mol"
            ]
        ),
        "post_5ps_CAP_SOL_LJ_max_kJ_mol": (
            post_5ps[
                "CAP_SOL_LJ_max_kJ_mol"
            ]
        ),
        "post_5ps_minimum_lumen_occupancy": int(
            np.min(
                occupancy_post_5
            )
        ),
        "post_5ps_minimum_initial_lumen_retention": int(
            np.min(
                retained_post_5
            )
        ),
        "post_5ps_minimum_CAP_OW_distance_nm": float(
            np.min(
                cap_distance_post_5
            )
        ),
        "failed_gates": (
            " | ".join(
                failed_gates
            )
        ),
        "checkpoint_continuation_preparation_authorized": (
            transient_only
        ),
        "long_mobile_MD_authorized": False,
        "multitemperature_MD_authorized": False,
        "QM_recalculation_authorized": False,
        "required_next_step": (
            required_next_step
        ),
    }

    write_csv(
        SUMMARY_CSV,
        [summary],
    )

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed in gates.items()
    )

    window_lines = "\n".join(
        (
            f"- {float(row['window_start_ps']):.2f}–"
            f"{float(row['window_end_ps']):.2f} ps: "
            f"T={float(row['temperature_mean_K']):.4f} ± "
            f"{float(row['temperature_std_K']):.4f} K; "
            f"range={float(row['temperature_min_K']):.4f}–"
            f"{float(row['temperature_max_K']):.4f} K; "
            f"slope={float(row['temperature_slope_K_per_ps']):.6f} K/ps; "
            f"potential slope="
            f"{float(row['potential_slope_kJ_mol_ps']):.3f} "
            f"kJ mol^-1 ps^-1"
        )
        for row in window_rows
    )

    REPORT_MD.write_text(
        f"""# R1 NVT 20 ps Thermal-Transient Review

## Purpose

The original 20 ps screening completed without dynamical instability
and failed only the full-trajectory temperature-standard-deviation
gate.

This review determines whether the failure is confined to initial
thermalization or persists after equilibration.

## Whole-trajectory temperature

- First temperature:
  **{first_temperature_K:.6f} K**
- Second temperature:
  **{second_temperature_K:.6f} K**
- Whole-trajectory mean:
  **{np.mean(temperature):.6f} K**
- Whole-trajectory standard deviation:
  **{all_temperature_std_K:.6f} K**
- Standard deviation after removing only the first point:
  **{without_first_temperature_std_K:.6f} K**
- Canonical temperature-fluctuation estimate from
  {degrees_of_freedom:.0f} solvent degrees of freedom:
  **{theoretical_temperature_std_K:.6f} K**

## Temporal windows

{window_lines}

## Confinement after 5 ps

- Minimum lumen occupancy:
  **{int(np.min(occupancy_post_5))} waters**
- Minimum retained initially luminal waters:
  **{int(np.min(retained_post_5))} waters**
- Endpoint occupancy:
  **{int(occupancy[-1])} waters**
- Endpoint initially luminal waters retained:
  **{int(retained_initial[-1])} waters**
- Minimum CAP–OW distance:
  **{np.min(cap_distance_post_5):.6f} nm**

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Preparation of a checkpoint continuation to 50 ps authorized:
  **{'YES' if transient_only else 'NO'}**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `{required_next_step}`

A passing result authorizes only preparation of a 30 ps continuation
from the existing 20 ps checkpoint. It does not authorize velocity
regeneration, a new independent trajectory, mobile-solute MD, or
multitemperature production.
""",
        encoding="utf-8",
    )

    print(
        "Day023 R1 NVT 20 ps thermal-transient "
        "review completed."
    )

    print(
        "Energy points / time range: "
        f"{len(times)} / "
        f"{times[0]:.3f}-"
        f"{times[-1]:.3f} ps"
    )

    print(
        "First / second temperature: "
        f"{first_temperature_K:.4f}/"
        f"{second_temperature_K:.4f} K"
    )

    print(
        "Temperature full mean/std: "
        f"{np.mean(temperature):.4f}/"
        f"{all_temperature_std_K:.4f} K"
    )

    print(
        "Temperature std without first point / "
        "theoretical canonical std: "
        f"{without_first_temperature_std_K:.4f}/"
        f"{theoretical_temperature_std_K:.4f} K"
    )

    print(
        "Post-5ps temperature mean/std/min/max/slope: "
        f"{float(post_5ps['temperature_mean_K']):.4f}/"
        f"{float(post_5ps['temperature_std_K']):.4f}/"
        f"{float(post_5ps['temperature_min_K']):.4f}/"
        f"{float(post_5ps['temperature_max_K']):.4f}/"
        f"{float(post_5ps['temperature_slope_K_per_ps']):.6f} "
        "K/ps"
    )

    print(
        "Last-10ps temperature mean/std/min/max/slope: "
        f"{float(last_10ps['temperature_mean_K']):.4f}/"
        f"{float(last_10ps['temperature_std_K']):.4f}/"
        f"{float(last_10ps['temperature_min_K']):.4f}/"
        f"{float(last_10ps['temperature_max_K']):.4f}/"
        f"{float(last_10ps['temperature_slope_K_per_ps']):.6f} "
        "K/ps"
    )

    print(
        "Post-5ps potential change/slope: "
        f"{float(post_5ps['potential_change_kJ_mol']):.3f}/"
        f"{float(post_5ps['potential_slope_kJ_mol_ps']):.3f} "
        "kJ mol^-1 / kJ mol^-1 ps^-1"
    )

    print(
        "Last-10ps potential change/slope: "
        f"{float(last_10ps['potential_change_kJ_mol']):.3f}/"
        f"{float(last_10ps['potential_slope_kJ_mol_ps']):.3f} "
        "kJ mol^-1 / kJ mol^-1 ps^-1"
    )

    print(
        "Post-5ps occupancy min / "
        "initial-water retention min / endpoint: "
        f"{int(np.min(occupancy_post_5))}/"
        f"{int(np.min(retained_post_5))}/"
        f"{int(occupancy[-1])}"
    )

    print(
        "Post-5ps minimum CAP-OW distance: "
        f"{float(np.min(cap_distance_post_5)):.6f} nm"
    )

    print(
        f"Decision: {decision}"
    )

    print(
        "Failed gates: "
        + (
            "NONE"
            if not failed_gates
            else " | ".join(
                failed_gates
            )
        )
    )

    print(
        "Checkpoint continuation preparation authorized: "
        f"{'YES' if transient_only else 'NO'}"
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

    print(
        f"Wrote: {relative(WINDOW_CSV)}"
    )

    print(
        f"Wrote: {relative(INITIAL_SAMPLES_CSV)}"
    )

    print(
        f"Wrote: {relative(GATE_CSV)}"
    )

    print(
        f"Wrote: {relative(SUMMARY_CSV)}"
    )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )

    if not transient_only:
        raise RuntimeError(
            "R1 post-equilibration thermal stability "
            "requires review."
        )


if __name__ == "__main__":
    main()
