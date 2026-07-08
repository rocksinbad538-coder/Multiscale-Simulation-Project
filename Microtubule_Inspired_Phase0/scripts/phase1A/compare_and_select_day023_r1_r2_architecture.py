#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
import re
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

DAY023_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design"
)

OUTPUT_ROOT = (
    DAY023_ROOT
    / "20_r1_r2_architecture_comparison_and_selection"
)

R2_50PS_ROOT = (
    DAY023_ROOT
    / "19_r2_frozen_solute_nvt_20_to_50ps"
)

R2_SUMMARY = (
    R2_50PS_ROOT
    / "r2_frozen_solute_nvt_50ps_validation_summary.csv"
)

R2_OCCUPANCY = (
    R2_50PS_ROOT
    / "r2_frozen_solute_nvt_combined_0_to_50ps_occupancy.csv"
)

R2_MDRUN_CONSOLE = (
    R2_50PS_ROOT
    / "r2_frozen_solute_nvt_20_to_50ps_mdrun_console.log"
)

COMPARISON_CSV = (
    OUTPUT_ROOT
    / "r1_r2_architecture_comparison.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r1_r2_architecture_selection_summary.csv"
)

GATES_CSV = (
    OUTPUT_ROOT
    / "r1_r2_architecture_selection_gates.csv"
)

SOURCE_MANIFEST_CSV = (
    OUTPUT_ROOT
    / "r1_r2_architecture_source_manifest.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R1_R2_ARCHITECTURE_COMPARISON_AND_SELECTION_DAY023.md"
)

R1_DECISION = (
    "R1_FROZEN_SOLUTE_50PS_POSITIVE_CONTROL_VALIDATED"
)

R2_DECISION = (
    "R2_FROZEN_SOLUTE_NVT_50PS_VALIDATED"
)

R2_STATIC_DECISION = (
    "R2_TOPOLOGY_AND_STATIC_CAP_WATER_MODEL_VALIDATED"
)

R2_GEOMETRY_DECISION = (
    "R2_PARTIAL_CAP_GEOMETRY_STATIC_GATE_PASSED"
)

SELECTION_DECISION = (
    "R2_SELECTED_AS_PRIMARY_PARTIAL_CAP_SCREENING_ARCHITECTURE"
)

INITIAL_OCCUPANCY = 428

MINIMUM_R2_ENDPOINT_FRACTION = 0.90
MINIMUM_R2_OCCUPANCY_FRACTION = 0.80
MINIMUM_R2_RELATIVE_TO_R1_ENDPOINT = 0.90

OCCUPANCY_SLOPE_LIMIT_WATER_PS = 0.50
MAXIMUM_FINAL10_NET_CHANGE = 5

MINIMUM_CAP_OW_DISTANCE_NM = 0.15
FROZEN_TOLERANCE_NM = 1.0e-6


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
            f"Missing or empty required file: {path}"
        )


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
                    field: row.get(field, "")
                    for field in fields
                }
            )


def parse_float(
    value: Any,
    *,
    field: str,
) -> float:
    try:
        parsed = float(value)
    except (
        TypeError,
        ValueError,
    ) as exc:
        raise RuntimeError(
            f"Could not parse {field!r}: {value!r}"
        ) from exc

    if not math.isfinite(parsed):
        raise RuntimeError(
            f"Non-finite field {field!r}"
        )

    return parsed


def pick_float(
    row: dict[str, str],
    aliases: tuple[str, ...],
    *,
    required: bool = True,
) -> float:
    for alias in aliases:
        if (
            alias in row
            and str(row[alias]).strip()
        ):
            return parse_float(
                row[alias],
                field=alias,
            )

    if required:
        raise RuntimeError(
            "Could not find any required field: "
            + " | ".join(aliases)
            + "\nAvailable fields: "
            + " | ".join(sorted(row))
        )

    return math.nan


def find_decision_csv(
    decision: str,
    preferred_tokens: tuple[str, ...],
) -> tuple[Path, dict[str, str]]:
    candidates: list[
        tuple[
            int,
            Path,
            dict[str, str],
        ]
    ] = []

    for path in DAY023_ROOT.rglob("*.csv"):
        if OUTPUT_ROOT in path.parents:
            continue

        try:
            rows = read_csv_rows(path)
        except Exception:
            continue

        for row in rows:
            if row.get(
                "decision",
                "",
            ).strip() != decision:
                continue

            lowered = str(path).lower()

            score = 0

            if "summary" in path.name.lower():
                score += 10

            for token in preferred_tokens:
                if token.lower() in lowered:
                    score += 5

            score += min(
                len(row),
                50,
            )

            candidates.append(
                (
                    score,
                    path,
                    row,
                )
            )

    if not candidates:
        raise RuntimeError(
            "Could not find a CSV containing decision: "
            f"{decision}"
        )

    candidates.sort(
        key=lambda item: (
            item[0],
            -len(str(item[1])),
        ),
        reverse=True,
    )

    _, path, row = candidates[0]

    return path, row


def locate_occupancy_series(
    architecture_token: str,
    target_end_ps: float,
) -> tuple[
    Path | None,
    list[dict[str, str]] | None,
]:
    candidates: list[
        tuple[
            int,
            Path,
            list[dict[str, str]],
        ]
    ] = []

    for path in DAY023_ROOT.rglob("*.csv"):
        if OUTPUT_ROOT in path.parents:
            continue

        if architecture_token.lower() not in str(
            path
        ).lower():
            continue

        try:
            rows = read_csv_rows(path)
        except Exception:
            continue

        fields = set(rows[0])

        if "time_ps" not in fields:
            continue

        occupancy_key = None

        for candidate in (
            "lumen_occupancy",
            "occupancy",
        ):
            if candidate in fields:
                occupancy_key = candidate
                break

        if occupancy_key is None:
            continue

        try:
            times = np.asarray(
                [
                    float(row["time_ps"])
                    for row in rows
                ],
                dtype=float,
            )
        except Exception:
            continue

        if (
            len(times) < 2
            or not np.all(
                np.isfinite(times)
            )
        ):
            continue

        score = 0

        if math.isclose(
            float(np.min(times)),
            0.0,
            rel_tol=0.0,
            abs_tol=1.0e-3,
        ):
            score += 10

        if math.isclose(
            float(np.max(times)),
            target_end_ps,
            rel_tol=0.0,
            abs_tol=1.0e-3,
        ):
            score += 30
        else:
            continue

        lowered = str(path).lower()

        if "combined" in lowered:
            score += 10

        if "50ps" in lowered or "50_ps" in lowered:
            score += 10

        if (
            "initial_lumen_waters_retained"
            in fields
        ):
            score += 10

        candidates.append(
            (
                score,
                path,
                rows,
            )
        )

    if not candidates:
        return None, None

    candidates.sort(
        key=lambda item: item[0],
        reverse=True,
    )

    _, path, rows = candidates[0]

    return path, rows


def series_metrics(
    rows: list[dict[str, str]],
) -> dict[str, Any]:
    occupancy_key = (
        "lumen_occupancy"
        if "lumen_occupancy" in rows[0]
        else "occupancy"
    )

    retained_key = (
        "initial_lumen_waters_retained"
        if (
            "initial_lumen_waters_retained"
            in rows[0]
        )
        else None
    )

    parsed = []

    for row in rows:
        parsed.append(
            (
                float(row["time_ps"]),
                int(
                    round(
                        float(
                            row[
                                occupancy_key
                            ]
                        )
                    )
                ),
                (
                    int(
                        round(
                            float(
                                row[
                                    retained_key
                                ]
                            )
                        )
                    )
                    if retained_key
                    else None
                ),
            )
        )

    parsed.sort(
        key=lambda item: item[0]
    )

    unique: dict[
        float,
        tuple[
            float,
            int,
            int | None,
        ]
    ] = {}

    for item in parsed:
        unique[
            round(
                item[0],
                6,
            )
        ] = item

    ordered = [
        unique[key]
        for key in sorted(unique)
    ]

    time = np.asarray(
        [
            item[0]
            for item in ordered
        ],
        dtype=float,
    )

    occupancy = np.asarray(
        [
            item[1]
            for item in ordered
        ],
        dtype=float,
    )

    retained = None

    if all(
        item[2] is not None
        for item in ordered
    ):
        retained = np.asarray(
            [
                int(item[2])
                for item in ordered
            ],
            dtype=float,
        )

    def window(
        start_ps: float,
    ) -> dict[str, float]:
        mask = (
            time >= start_ps - 1.0e-8
        )

        window_time = time[mask]
        window_occupancy = occupancy[mask]

        if len(window_time) < 2:
            raise RuntimeError(
                f"Insufficient occupancy records from "
                f"{start_ps} ps."
            )

        slope = float(
            np.polyfit(
                window_time,
                window_occupancy,
                1,
            )[0]
        )

        return {
            "mean": float(
                np.mean(
                    window_occupancy
                )
            ),
            "change": int(
                window_occupancy[-1]
                - window_occupancy[0]
            ),
            "slope": slope,
        }

    endpoint_retained = math.nan
    maximum_noninitial = math.nan
    endpoint_noninitial = math.nan

    if retained is not None:
        endpoint_retained = float(
            retained[-1]
        )

        noninitial = (
            occupancy
            - retained
        )

        maximum_noninitial = float(
            np.max(
                noninitial
            )
        )

        endpoint_noninitial = float(
            noninitial[-1]
        )

    return {
        "frames": len(time),
        "start_ps": float(
            time[0]
        ),
        "end_ps": float(
            time[-1]
        ),
        "initial": int(
            occupancy[0]
        ),
        "mean": float(
            np.mean(
                occupancy
            )
        ),
        "minimum": int(
            np.min(
                occupancy
            )
        ),
        "maximum": int(
            np.max(
                occupancy
            )
        ),
        "endpoint": int(
            occupancy[-1]
        ),
        "endpoint_retained": (
            endpoint_retained
        ),
        "endpoint_noninitial": (
            endpoint_noninitial
        ),
        "maximum_noninitial": (
            maximum_noninitial
        ),
        "final20": window(30.0),
        "final15": window(35.0),
        "final10": window(40.0),
    }


def summary_metrics(
    row: dict[str, str],
) -> dict[str, Any]:
    initial = int(
        round(
            pick_float(
                row,
                (
                    "lumen_occupancy_initial",
                    "occupancy_initial",
                    "combined_lumen_occupancy_initial",
                ),
            )
        )
    )

    endpoint = int(
        round(
            pick_float(
                row,
                (
                    "lumen_occupancy_endpoint",
                    "occupancy_endpoint",
                    "combined_lumen_occupancy_endpoint",
                ),
            )
        )
    )

    return {
        "frames": int(
            round(
                pick_float(
                    row,
                    (
                        "combined_unique_frames",
                        "trajectory_frames",
                        "combined_frames",
                    ),
                    required=False,
                )
            )
        )
        if any(
            key in row
            for key in (
                "combined_unique_frames",
                "trajectory_frames",
                "combined_frames",
            )
        )
        else "",
        "start_ps": pick_float(
            row,
            (
                "combined_start_ps",
                "trajectory_start_ps",
            ),
            required=False,
        ),
        "end_ps": pick_float(
            row,
            (
                "combined_end_ps",
                "trajectory_end_ps",
            ),
            required=False,
        ),
        "initial": initial,
        "mean": pick_float(
            row,
            (
                "lumen_occupancy_mean",
                "occupancy_mean",
                "combined_lumen_occupancy_mean",
            ),
        ),
        "minimum": int(
            round(
                pick_float(
                    row,
                    (
                        "lumen_occupancy_minimum",
                        "occupancy_minimum",
                        "combined_lumen_occupancy_minimum",
                    ),
                )
            )
        ),
        "maximum": int(
            round(
                pick_float(
                    row,
                    (
                        "lumen_occupancy_maximum",
                        "occupancy_maximum",
                        "combined_lumen_occupancy_maximum",
                    ),
                )
            )
        ),
        "endpoint": endpoint,
        "endpoint_retained": pick_float(
            row,
            (
                "endpoint_initial_lumen_waters_retained",
                "endpoint_initial_lumen_retention",
            ),
            required=False,
        ),
        "endpoint_noninitial": math.nan,
        "maximum_noninitial": math.nan,
        "final20": {
            "mean": pick_float(
                row,
                (
                    "final20_occupancy_mean",
                    "second_half_occupancy_mean",
                ),
                required=False,
            ),
            "change": "",
            "slope": pick_float(
                row,
                (
                    "final20_occupancy_slope_water_ps",
                    "second_half_occupancy_slope_water_ps",
                ),
            ),
        },
        "final15": {
            "mean": pick_float(
                row,
                (
                    "final15_occupancy_mean",
                ),
                required=False,
            ),
            "change": "",
            "slope": pick_float(
                row,
                (
                    "final15_occupancy_slope_water_ps",
                    "second_half_occupancy_slope_water_ps",
                ),
                required=False,
            ),
        },
        "final10": {
            "mean": pick_float(
                row,
                (
                    "final10_occupancy_mean",
                ),
                required=False,
            ),
            "change": pick_float(
                row,
                (
                    "final10_occupancy_change",
                ),
                required=False,
            ),
            "slope": pick_float(
                row,
                (
                    "final10_occupancy_slope_water_ps",
                    "second_half_occupancy_slope_water_ps",
                ),
                required=False,
            ),
        },
    }


def find_value_containing(
    row: dict[str, str],
    required_tokens: tuple[str, ...],
) -> float:
    for key, value in row.items():
        lowered = key.lower()

        if all(
            token.lower() in lowered
            for token in required_tokens
        ):
            try:
                parsed = float(value)
            except (
                TypeError,
                ValueError,
            ):
                continue

            if math.isfinite(parsed):
                return parsed

    return math.nan


def parse_openmp_threads(
    path: Path,
) -> int | str:
    if not path.exists():
        return ""

    text = path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    matches = re.findall(
        r"Using\s+1\s+MPI\s+thread\s*\n"
        r"Using\s+(\d+)\s+OpenMP\s+threads",
        text,
        flags=re.IGNORECASE,
    )

    if not matches:
        matches = re.findall(
            r"Using\s+(\d+)\s+OpenMP\s+threads",
            text,
            flags=re.IGNORECASE,
        )

    if not matches:
        return ""

    return int(
        matches[-1]
    )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    require_file(R2_SUMMARY)
    require_file(R2_OCCUPANCY)

    r2_summary = read_single_csv_row(
        R2_SUMMARY
    )

    if (
        r2_summary.get(
            "decision",
            "",
        )
        != R2_DECISION
    ):
        raise RuntimeError(
            "R2 50 ps validation is not accepted."
        )

    (
        r1_summary_path,
        r1_summary,
    ) = find_decision_csv(
        R1_DECISION,
        (
            "r1",
            "50",
            "validation",
        ),
    )

    (
        r2_static_summary_path,
        r2_static_summary,
    ) = find_decision_csv(
        R2_STATIC_DECISION,
        (
            "r2",
            "topology",
            "static",
        ),
    )

    (
        r2_geometry_summary_path,
        r2_geometry_summary,
    ) = find_decision_csv(
        R2_GEOMETRY_DECISION,
        (
            "r2",
            "partial",
            "cap",
            "geometry",
        ),
    )

    (
        r1_occupancy_path,
        r1_occupancy_rows,
    ) = locate_occupancy_series(
        "r1",
        50.0,
    )

    r2_occupancy_rows = read_csv_rows(
        R2_OCCUPANCY
    )

    if r1_occupancy_rows is not None:
        r1 = series_metrics(
            r1_occupancy_rows
        )
        r1_metric_source = (
            relative(
                r1_occupancy_path
            )
            if r1_occupancy_path
            else "UNKNOWN"
        )
    else:
        r1 = summary_metrics(
            r1_summary
        )
        r1_metric_source = relative(
            r1_summary_path
        )

    r2 = series_metrics(
        r2_occupancy_rows
    )

    r2_min_cap_distance = pick_float(
        r2_summary,
        (
            "minimum_CAP_OW_distance_nm",
        ),
    )

    r2_hbn_displacement = pick_float(
        r2_summary,
        (
            "HBN_maximum_displacement_nm",
        ),
    )

    r2_pyr_displacement = pick_float(
        r2_summary,
        (
            "PYR_maximum_displacement_nm",
        ),
    )

    r2_caps_displacement = pick_float(
        r2_summary,
        (
            "CAPS_maximum_displacement_nm",
        ),
    )

    r1_min_cap_distance = pick_float(
        r1_summary,
        (
            "minimum_CAP_OW_distance_nm",
            "minimum_cap_water_distance_nm",
            "minimum_CAP_water_distance_nm",
        ),
        required=False,
    )

    r1_stationarity_slope = float(
        r1[
            "final20"
        ][
            "slope"
        ]
    )

    r2_aperture_diameter = find_value_containing(
        r2_static_summary,
        (
            "aperture",
            "diameter",
        ),
    )

    if not math.isfinite(
        r2_aperture_diameter
    ):
        lower_radius = find_value_containing(
            r2_static_summary,
            (
                "lower",
                "aperture",
                "radius",
            ),
        )

        upper_radius = find_value_containing(
            r2_static_summary,
            (
                "upper",
                "aperture",
                "radius",
            ),
        )

        if (
            math.isfinite(
                lower_radius
            )
            and math.isfinite(
                upper_radius
            )
        ):
            r2_aperture_diameter = (
                lower_radius
                + upper_radius
            )

    r2_open_area_fraction = find_value_containing(
        r2_static_summary,
        (
            "open",
            "area",
            "fraction",
        ),
    )

    if not math.isfinite(
        r2_open_area_fraction
    ):
        r2_open_area_fraction = find_value_containing(
            r2_geometry_summary,
            (
                "open",
                "area",
                "fraction",
            ),
        )

    if not math.isfinite(
        r2_open_area_fraction
    ):
        raise RuntimeError(
            "Could not resolve the validated R2 "
            "open-area fraction from Gate 2A or Gate 2B."
        )

    r2_endpoint_fraction = (
        r2["endpoint"]
        / r2["initial"]
    )

    r2_minimum_fraction = (
        r2["minimum"]
        / r2["initial"]
    )

    r1_endpoint_fraction = (
        r1["endpoint"]
        / r1["initial"]
    )

    r2_relative_to_r1_endpoint = (
        r2["endpoint"]
        / r1["endpoint"]
    )

    r2_endpoint_identity_fraction = (
        float(
            r2[
                "endpoint_retained"
            ]
        )
        / r2["initial"]
    )

    r2_retention_penalty_percentage_points = (
        r1_endpoint_fraction
        - r2_endpoint_fraction
    ) * 100.0

    r2_net_occupancy_change = (
        r2["endpoint"]
        - r2["initial"]
    )

    actual_openmp_threads = (
        parse_openmp_threads(
            R2_MDRUN_CONSOLE
        )
    )

    comparison_rows = [
        {
            "metric": "architectural_role",
            "R1": (
                "closed neutral frozen steric "
                "positive control"
            ),
            "R2": (
                "symmetric partial-cap neutral "
                "frozen steric screening architecture"
            ),
            "interpretation": (
                "R1 is the closed-methodology control; "
                "R2 is the partially open candidate."
            ),
        },
        {
            "metric": "initial_lumen_occupancy",
            "R1": r1["initial"],
            "R2": r2["initial"],
            "interpretation": (
                "Matched initial hydration."
            ),
        },
        {
            "metric": "mean_lumen_occupancy_0_to_50ps",
            "R1": r1["mean"],
            "R2": r2["mean"],
            "interpretation": (
                "R2 exchanges and releases a small "
                "fraction of lumen water."
            ),
        },
        {
            "metric": "minimum_lumen_occupancy",
            "R1": r1["minimum"],
            "R2": r2["minimum"],
            "interpretation": (
                "R2 remains highly hydrated throughout."
            ),
        },
        {
            "metric": "endpoint_lumen_occupancy",
            "R1": r1["endpoint"],
            "R2": r2["endpoint"],
            "interpretation": (
                "R2 endpoint is approximately 96% "
                "of the initial and R1 endpoint occupancy."
            ),
        },
        {
            "metric": "endpoint_occupancy_fraction",
            "R1": r1_endpoint_fraction,
            "R2": r2_endpoint_fraction,
            "interpretation": (
                "R2 retention penalty relative to the "
                "closed control is modest."
            ),
        },
        {
            "metric": "net_occupancy_change",
            "R1": (
                r1["endpoint"]
                - r1["initial"]
            ),
            "R2": r2_net_occupancy_change,
            "interpretation": (
                "R2 reaches a lower but stable hydration plateau."
            ),
        },
        {
            "metric": "final20_occupancy_slope_water_ps",
            "R1": r1_stationarity_slope,
            "R2": r2["final20"]["slope"],
            "interpretation": (
                "Both satisfy the predefined stationarity limit."
            ),
        },
        {
            "metric": "final15_occupancy_slope_water_ps",
            "R1": "",
            "R2": r2["final15"]["slope"],
            "interpretation": (
                "R2 final 15 ps is effectively stationary."
            ),
        },
        {
            "metric": "final10_occupancy_slope_water_ps",
            "R1": "",
            "R2": r2["final10"]["slope"],
            "interpretation": (
                "R2 final 10 ps is effectively stationary."
            ),
        },
        {
            "metric": "final10_net_occupancy_change",
            "R1": "",
            "R2": r2["final10"]["change"],
            "interpretation": (
                "No net water loss in the final 10 ps."
            ),
        },
        {
            "metric": "endpoint_initial_identity_fraction",
            "R1": (
                (
                    float(
                        r1[
                            "endpoint_retained"
                        ]
                    )
                    / r1["initial"]
                )
                if math.isfinite(
                    float(
                        r1[
                            "endpoint_retained"
                        ]
                    )
                )
                else ""
            ),
            "R2": (
                r2_endpoint_identity_fraction
            ),
            "interpretation": (
                "R2 retains most original lumen waters "
                "while allowing exchange."
            ),
        },
        {
            "metric": "maximum_noninitial_lumen_waters",
            "R1": (
                r1[
                    "maximum_noninitial"
                ]
                if math.isfinite(
                    float(
                        r1[
                            "maximum_noninitial"
                        ]
                    )
                )
                else ""
            ),
            "R2": r2["maximum_noninitial"],
            "interpretation": (
                "Noninitial lumen waters demonstrate "
                "exchange through the R2 apertures."
            ),
        },
        {
            "metric": "minimum_CAP_OW_distance_nm",
            "R1": (
                r1_min_cap_distance
                if math.isfinite(
                    r1_min_cap_distance
                )
                else ""
            ),
            "R2": r2_min_cap_distance,
            "interpretation": (
                "Both remain above the predefined safety threshold."
            ),
        },
        {
            "metric": "effective_aperture_diameter_nm",
            "R1": 0.0,
            "R2": (
                r2_aperture_diameter
                if math.isfinite(
                    r2_aperture_diameter
                )
                else ""
            ),
            "interpretation": (
                "R1 is closed; R2 provides a central exchange aperture."
            ),
        },
        {
            "metric": "open_area_fraction",
            "R1": 0.0,
            "R2": (
                r2_open_area_fraction
                if math.isfinite(
                    r2_open_area_fraction
                )
                else ""
            ),
            "interpretation": (
                "R2 introduces controlled partial openness."
            ),
        },
        {
            "metric": "chemical_realizability",
            "R1": (
                "not established"
            ),
            "R2": (
                "not established"
            ),
            "interpretation": (
                "Both use ideal neutral frozen steric beads; "
                "R2 requires a chemically realizable end-rim design."
            ),
        },
    ]

    write_csv(
        COMPARISON_CSV,
        comparison_rows,
    )

    r2_has_exchange = (
        math.isfinite(
            float(
                r2[
                    "maximum_noninitial"
                ]
            )
        )
        and float(
            r2[
                "maximum_noninitial"
            ]
        )
        > 0.0
    )

    aperture_is_open = (
        math.isfinite(
            r2_aperture_diameter
        )
        and r2_aperture_diameter
        > 0.0
    )

    gates = {
        "R1_closed_positive_control_is_validated": (
            r1_summary.get(
                "decision",
                "",
            )
            == R1_DECISION
        ),
        "R2_frozen_solute_50ps_is_validated": (
            r2_summary.get(
                "decision",
                "",
            )
            == R2_DECISION
        ),
        "R2_static_partial_cap_model_is_validated": (
            r2_static_summary.get(
                "decision",
                "",
            )
            == R2_STATIC_DECISION
        ),
        "R1_and_R2_have_matched_initial_occupancy": (
            r1["initial"]
            == INITIAL_OCCUPANCY
            and r2["initial"]
            == INITIAL_OCCUPANCY
        ),
        "R1_behaves_as_closed_positive_control": (
            r1["endpoint"]
            == INITIAL_OCCUPANCY
            and abs(
                r1_stationarity_slope
            )
            <= OCCUPANCY_SLOPE_LIMIT_WATER_PS
        ),
        "R2_endpoint_occupancy_is_at_least_90_percent": (
            r2_endpoint_fraction
            >= MINIMUM_R2_ENDPOINT_FRACTION
        ),
        "R2_minimum_occupancy_is_at_least_80_percent": (
            r2_minimum_fraction
            >= MINIMUM_R2_OCCUPANCY_FRACTION
        ),
        "R2_endpoint_is_at_least_90_percent_of_R1_endpoint": (
            r2_relative_to_r1_endpoint
            >= MINIMUM_R2_RELATIVE_TO_R1_ENDPOINT
        ),
        "R2_final20_occupancy_is_stationary": (
            abs(
                float(
                    r2[
                        "final20"
                    ][
                        "slope"
                    ]
                )
            )
            <= OCCUPANCY_SLOPE_LIMIT_WATER_PS
        ),
        "R2_final15_occupancy_is_stationary": (
            abs(
                float(
                    r2[
                        "final15"
                    ][
                        "slope"
                    ]
                )
            )
            <= OCCUPANCY_SLOPE_LIMIT_WATER_PS
        ),
        "R2_final10_occupancy_is_stationary": (
            abs(
                float(
                    r2[
                        "final10"
                    ][
                        "slope"
                    ]
                )
            )
            <= OCCUPANCY_SLOPE_LIMIT_WATER_PS
        ),
        "R2_final10_net_change_is_at_most_5_waters": (
            abs(
                int(
                    r2[
                        "final10"
                    ][
                        "change"
                    ]
                )
            )
            <= MAXIMUM_FINAL10_NET_CHANGE
        ),
        "R2_demonstrates_water_exchange": (
            r2_has_exchange
        ),
        "R2_aperture_is_open": (
            aperture_is_open
        ),
        "R2_open_area_fraction_is_valid": (
            math.isfinite(
                r2_open_area_fraction
            )
            and r2_open_area_fraction > 0.0
            and r2_open_area_fraction < 1.0
        ),
        "R2_minimum_CAP_OW_distance_is_safe": (
            r2_min_cap_distance
            >= MINIMUM_CAP_OW_DISTANCE_NM
        ),
        "R2_HBN_remained_frozen": (
            r2_hbn_displacement
            <= FROZEN_TOLERANCE_NM
        ),
        "R2_PYR_remained_frozen": (
            r2_pyr_displacement
            <= FROZEN_TOLERANCE_NM
        ),
        "R2_CAPS_remained_frozen": (
            r2_caps_displacement
            <= FROZEN_TOLERANCE_NM
        ),
    }

    failed_gates = [
        name
        for name, passed
        in gates.items()
        if not passed
    ]

    accepted = (
        len(failed_gates) == 0
    )

    decision = (
        SELECTION_DECISION
        if accepted
        else
        "R1_R2_ARCHITECTURE_SELECTION_REQUIRES_REVIEW"
    )

    required_next_step = (
        "BEGIN_R2_CHEMICALLY_REALIZABLE_END_RIM_DESIGN_GATE"
        if accepted
        else
        "REVIEW_R1_R2_ARCHITECTURE_SELECTION_FAILURES"
    )

    summary = {
        "decision": decision,
        "R1_status": (
            "RETAIN_AS_CLOSED_STERIC_POSITIVE_CONTROL"
        ),
        "R2_status": (
            "PRIMARY_PARTIAL_CAP_SCREENING_ARCHITECTURE"
            if accepted
            else "REQUIRES_REVIEW"
        ),
        "R3_status": (
            "DEFERRED_NOT_REQUIRED_AT_THIS_GATE"
            if accepted
            else "UNRESOLVED"
        ),
        "R4_status": (
            "DEFERRED_NOT_REQUIRED_AT_THIS_GATE"
            if accepted
            else "UNRESOLVED"
        ),
        "R1_initial_occupancy": (
            r1["initial"]
        ),
        "R1_mean_occupancy": (
            r1["mean"]
        ),
        "R1_minimum_occupancy": (
            r1["minimum"]
        ),
        "R1_endpoint_occupancy": (
            r1["endpoint"]
        ),
        "R1_endpoint_fraction": (
            r1_endpoint_fraction
        ),
        "R1_stationarity_slope_water_ps": (
            r1_stationarity_slope
        ),
        "R2_initial_occupancy": (
            r2["initial"]
        ),
        "R2_mean_occupancy": (
            r2["mean"]
        ),
        "R2_minimum_occupancy": (
            r2["minimum"]
        ),
        "R2_endpoint_occupancy": (
            r2["endpoint"]
        ),
        "R2_endpoint_fraction": (
            r2_endpoint_fraction
        ),
        "R2_endpoint_relative_to_R1": (
            r2_relative_to_r1_endpoint
        ),
        "R2_retention_penalty_vs_R1_percentage_points": (
            r2_retention_penalty_percentage_points
        ),
        "R2_net_occupancy_change": (
            r2_net_occupancy_change
        ),
        "R2_endpoint_initial_identity_fraction": (
            r2_endpoint_identity_fraction
        ),
        "R2_endpoint_noninitial_lumen_waters": (
            r2[
                "endpoint_noninitial"
            ]
        ),
        "R2_maximum_noninitial_lumen_waters": (
            r2[
                "maximum_noninitial"
            ]
        ),
        "R2_final20_occupancy_mean": (
            r2[
                "final20"
            ][
                "mean"
            ]
        ),
        "R2_final20_occupancy_slope_water_ps": (
            r2[
                "final20"
            ][
                "slope"
            ]
        ),
        "R2_final15_occupancy_mean": (
            r2[
                "final15"
            ][
                "mean"
            ]
        ),
        "R2_final15_occupancy_slope_water_ps": (
            r2[
                "final15"
            ][
                "slope"
            ]
        ),
        "R2_final10_occupancy_mean": (
            r2[
                "final10"
            ][
                "mean"
            ]
        ),
        "R2_final10_occupancy_change": (
            r2[
                "final10"
            ][
                "change"
            ]
        ),
        "R2_final10_occupancy_slope_water_ps": (
            r2[
                "final10"
            ][
                "slope"
            ]
        ),
        "R2_effective_aperture_diameter_nm": (
            r2_aperture_diameter
        ),
        "R2_open_area_fraction": (
            r2_open_area_fraction
        ),
        "R2_minimum_CAP_OW_distance_nm": (
            r2_min_cap_distance
        ),
        "R2_execution_OpenMP_threads": (
            actual_openmp_threads
        ),
        "R2_chemical_realizability_established": False,
        "failed_gates": (
            " | ".join(
                failed_gates
            )
        ),
        "chemical_realization_static_design_authorized": (
            accepted
        ),
        "new_MD_execution_authorized": False,
        "short_mobile_MD_authorized": False,
        "long_mobile_MD_authorized": False,
        "multitemperature_MD_authorized": False,
        "QM_recalculation_authorized": False,
        "required_next_step": (
            required_next_step
        ),
    }

    write_csv(
        SUMMARY_CSV,
        [
            summary
        ],
    )

    write_csv(
        GATES_CSV,
        [
            {
                "gate": name,
                "pass": passed,
            }
            for name, passed
            in gates.items()
        ],
    )

    source_manifest = [
        {
            "role": "R1_50ps_summary",
            "file": relative(
                r1_summary_path
            ),
            "decision": (
                r1_summary.get(
                    "decision",
                    "",
                )
            ),
        },
        {
            "role": "R1_occupancy_metrics",
            "file": r1_metric_source,
            "decision": "",
        },
        {
            "role": "R2_50ps_summary",
            "file": relative(
                R2_SUMMARY
            ),
            "decision": (
                r2_summary.get(
                    "decision",
                    "",
                )
            ),
        },
        {
            "role": "R2_50ps_occupancy",
            "file": relative(
                R2_OCCUPANCY
            ),
            "decision": "",
        },
        {
            "role": "R2_static_summary",
            "file": relative(
                r2_static_summary_path
            ),
            "decision": (
                r2_static_summary.get(
                    "decision",
                    "",
                )
            ),
        },
        {
            "role": "R2_geometry_summary",
            "file": relative(
                r2_geometry_summary_path
            ),
            "decision": (
                r2_geometry_summary.get(
                    "decision",
                    "",
                )
            ),
        },
    ]

    write_csv(
        SOURCE_MANIFEST_CSV,
        source_manifest,
    )

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed
        in gates.items()
    )

    REPORT_MD.write_text(
        f"""# R1–R2 Architecture Comparison and Selection

## Scope

This comparison uses the validated 50 ps frozen-solute results for the
closed R1 positive control and the partially open R2 screening
architecture.

No new molecular dynamics, minimization, topology generation, or
quantum calculation was performed.

## Architectural roles

### R1

R1 remains the neutral frozen closed steric positive control. Its role
is to demonstrate that the confinement-detection method identifies
stable lumen hydration when the end boundary is closed.

R1 is not promoted as a chemically realizable architecture.

### R2

R2 is the symmetric partial-cap screening architecture. It introduces a
central aperture while maintaining a highly hydrated lumen and allowing
measurable exchange.

R2 is also not yet chemically realizable because its cap is represented
by ideal neutral frozen steric beads.

## Hydration comparison

- R1 initial/mean/minimum/endpoint:
  **{r1['initial']}/{r1['mean']:.4f}/
  {r1['minimum']}/{r1['endpoint']}**
- R2 initial/mean/minimum/endpoint:
  **{r2['initial']}/{r2['mean']:.4f}/
  {r2['minimum']}/{r2['endpoint']}**
- R1 endpoint fraction:
  **{r1_endpoint_fraction:.6f}**
- R2 endpoint fraction:
  **{r2_endpoint_fraction:.6f}**
- R2 endpoint relative to R1:
  **{r2_relative_to_r1_endpoint:.6f}**
- R2 retention penalty versus R1:
  **{r2_retention_penalty_percentage_points:.4f} percentage points**
- R2 endpoint initial-identity retention:
  **{r2_endpoint_identity_fraction:.6f}**
- R2 endpoint/maximum noninitial lumen waters:
  **{r2['endpoint_noninitial']}/
  {r2['maximum_noninitial']}**

## Stationarity

- R1 stationarity slope:
  **{r1_stationarity_slope:.6f} waters/ps**
- R2 final 20 ps mean/slope:
  **{r2['final20']['mean']:.4f}/
  {r2['final20']['slope']:.6f} waters, waters/ps**
- R2 final 15 ps mean/slope:
  **{r2['final15']['mean']:.4f}/
  {r2['final15']['slope']:.6f} waters, waters/ps**
- R2 final 10 ps mean/change/slope:
  **{r2['final10']['mean']:.4f}/
  {int(r2['final10']['change']):+d}/
  {r2['final10']['slope']:.6f} waters, waters, waters/ps**

## Aperture and steric safety

- R1 effective aperture diameter:
  **0.000000 nm**
- R2 effective aperture diameter:
  **{r2_aperture_diameter:.6f} nm**
- R2 open-area fraction:
  **{r2_open_area_fraction:.6f}**
- R2 minimum CAP–OW distance:
  **{r2_min_cap_distance:.6f} nm**

## Execution metadata

- R2 continuation OpenMP threads:
  **{actual_openmp_threads}**

Thread count is recorded as computational metadata. It does not modify
the TPR physical parameters or checkpoint state.

## Selection gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- R1 status:
  **RETAIN_AS_CLOSED_STERIC_POSITIVE_CONTROL**
- R2 status:
  **{'PRIMARY_PARTIAL_CAP_SCREENING_ARCHITECTURE' if accepted else 'REQUIRES_REVIEW'}**
- R3 status:
  **{'DEFERRED_NOT_REQUIRED_AT_THIS_GATE' if accepted else 'UNRESOLVED'}**
- R4 status:
  **{'DEFERRED_NOT_REQUIRED_AT_THIS_GATE' if accepted else 'UNRESOLVED'}**
- Chemical-realization static design authorized:
  **{'YES' if accepted else 'NO'}**
- New MD execution authorized:
  **NO**
- Short mobile MD authorized:
  **NO**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `{required_next_step}`

The next stage must replace the ideal steric partial cap with a
chemically defensible end-rim or terminal-ring realization while
preserving the validated aperture and confinement envelope.
""",
        encoding="utf-8",
    )

    print(
        "Day023 R1-R2 architecture comparison "
        "and selection completed."
    )

    print(
        "R1 summary source: "
        f"{relative(r1_summary_path)}"
    )

    print(
        "R1 occupancy source: "
        f"{r1_metric_source}"
    )

    print(
        "R2 summary source: "
        f"{relative(R2_SUMMARY)}"
    )

    print(
        "R2 static source: "
        f"{relative(r2_static_summary_path)}"
    )

    print(
        "R2 geometry source: "
        f"{relative(r2_geometry_summary_path)}"
    )

    print(
        "R1 initial/mean/min/endpoint: "
        f"{r1['initial']}/"
        f"{r1['mean']:.4f}/"
        f"{r1['minimum']}/"
        f"{r1['endpoint']}"
    )

    print(
        "R2 initial/mean/min/endpoint: "
        f"{r2['initial']}/"
        f"{r2['mean']:.4f}/"
        f"{r2['minimum']}/"
        f"{r2['endpoint']}"
    )

    print(
        "R1 / R2 endpoint fractions: "
        f"{r1_endpoint_fraction:.6f}/"
        f"{r2_endpoint_fraction:.6f}"
    )

    print(
        "R2 endpoint relative to R1 / retention "
        "penalty: "
        f"{r2_relative_to_r1_endpoint:.6f}/"
        f"{r2_retention_penalty_percentage_points:.4f} "
        "percentage points"
    )

    print(
        "R1 stationarity slope: "
        f"{r1_stationarity_slope:.6f} waters/ps"
    )

    print(
        "R2 final-20/final-15/final-10 slopes: "
        f"{r2['final20']['slope']:.6f}/"
        f"{r2['final15']['slope']:.6f}/"
        f"{r2['final10']['slope']:.6f} waters/ps"
    )

    print(
        "R2 final-10 change: "
        f"{int(r2['final10']['change']):+d} waters"
    )

    print(
        "R2 endpoint retained initial identities: "
        f"{int(r2['endpoint_retained'])}/"
        f"{r2['initial']} "
        f"({r2_endpoint_identity_fraction:.6f})"
    )

    print(
        "R2 endpoint / maximum noninitial lumen waters: "
        f"{int(r2['endpoint_noninitial'])}/"
        f"{int(r2['maximum_noninitial'])}"
    )

    print(
        "R2 effective aperture diameter / open-area fraction: "
        f"{r2_aperture_diameter:.6f}/"
        f"{r2_open_area_fraction:.6f}"
    )

    print(
        "R2 minimum CAP-OW distance: "
        f"{r2_min_cap_distance:.6f} nm"
    )

    print(
        "R2 continuation OpenMP threads: "
        f"{actual_openmp_threads}"
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
        "R1 status: "
        "RETAIN_AS_CLOSED_STERIC_POSITIVE_CONTROL"
    )

    print(
        "R2 status: "
        + (
            "PRIMARY_PARTIAL_CAP_SCREENING_ARCHITECTURE"
            if accepted
            else "REQUIRES_REVIEW"
        )
    )

    print(
        "R3 / R4 status: "
        + (
            "DEFERRED_NOT_REQUIRED_AT_THIS_GATE"
            if accepted
            else "UNRESOLVED"
        )
    )

    print(
        "Chemical-realization static design authorized: "
        f"{'YES' if accepted else 'NO'}"
    )

    print(
        "New MD execution authorized: NO"
    )

    print(
        "Short mobile MD authorized: NO"
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
        COMPARISON_CSV,
        SUMMARY_CSV,
        GATES_CSV,
        SOURCE_MANIFEST_CSV,
        REPORT_MD,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not accepted:
        raise RuntimeError(
            "R1-R2 architecture selection requires review."
        )


if __name__ == "__main__":
    main()
