#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/phase1A/day024_chemical_end_rim_design"

GATE3J = BASE / "12_r2_alternating_bn_trimer_bridge_static_coordinate_embedding"
GATE3K = BASE / "13_r2_trimer_bridge_conformer_and_h_refinement"
GATE3O = BASE / "19_r2_selected_four_atom_heavy_coordinate_embedding"

OUT = (
    BASE
    / "20_r2_aperture_metric_consistency_audit"
)

SCRIPT3J = (
    ROOT
    / "scripts/phase1A/"
    "build_and_validate_day024_r2_alternating_bn_trimer_bridge_static_coordinate_embedding.py"
)

SCRIPT3K = (
    ROOT
    / "scripts/phase1A/"
    "refine_day024_r2_trimer_bridge_conformers_and_h_orientations.py"
)

SCRIPT3O = (
    ROOT
    / "scripts/phase1A/"
    "build_and_validate_day024_r2_selected_four_atom_heavy_coordinate_embedding.py"
)

COORD3K = (
    GATE3K
    / "r2_trimer_bridge_refined_coordinates.csv"
)

COORD3O = (
    GATE3O
    / "r2_selected_four_atom_heavy_coordinates.csv"
)

SUMMARY3J = (
    GATE3J
    / "r2_trimer_bridge_static_embedding_summary.csv"
)

SUMMARY3K = (
    GATE3K
    / "r2_trimer_bridge_refinement_summary.csv"
)

SUMMARY3O = (
    GATE3O
    / "r2_selected_four_atom_heavy_embedding_summary.csv"
)

METRICS = (
    OUT
    / "r2_aperture_metric_comparison.csv"
)

PAIR_DETAILS = (
    OUT
    / "r2_aperture_opposite_pair_details.csv"
)

SOURCE_DEFINITIONS = (
    OUT
    / "r2_aperture_source_definition_extracts.txt"
)

SUMMARY = (
    OUT
    / "r2_aperture_metric_consistency_summary.csv"
)

GATES = (
    OUT
    / "r2_aperture_metric_consistency_gates.csv"
)

JSON_OUT = (
    OUT
    / "r2_aperture_metric_consistency.json"
)

MANIFEST = (
    OUT
    / "r2_aperture_metric_consistency_manifest.csv"
)

REPORT = (
    OUT
    / "R2_APERTURE_METRIC_CONSISTENCY_AUDIT_DAY024.md"
)

EXPECTED_GATE3O_DECISION = (
    "R2_SELECTED_FOUR_ATOM_BN_BRIDGE_"
    "HEAVY_COORDINATE_EMBEDDING_REQUIRES_REVIEW"
)

PASS_DECISION = (
    "R2_APERTURE_METRIC_INCONSISTENCY_IDENTIFIED_"
    "HEAVY_EMBEDDING_GEOMETRY_RETAINED"
)

FAIL_DECISION = (
    "R2_APERTURE_METRIC_AUDIT_REQUIRES_FURTHER_REVIEW"
)

TARGET_APERTURE_NM = 0.839406
MAX_ACCEPTED_RELATIVE_ERROR = 0.10
COORDINATE_TOLERANCE_NM = 1.0e-12


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty required file: {path}"
        )


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(
            f"No rows found in {path}"
        )

    return rows


def read_one(path: Path) -> dict[str, str]:
    rows = read_rows(path)

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected one row in {path}; found {len(rows)}"
        )

    return rows[0]


def write_rows(
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
            f"Non-finite field {key!r}"
        )

    return value


def normalized(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))

    if norm <= 1.0e-12:
        raise RuntimeError(
            "Could not normalize zero vector."
        )

    return vector / norm


def principal_axis(
    positions: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    center = np.mean(
        positions,
        axis=0,
    )

    centered = positions - center

    values, vectors = np.linalg.eigh(
        centered.T @ centered
    )

    axis = normalized(
        vectors[
            :,
            int(np.argmax(values)),
        ]
    )

    if axis[2] < 0.0:
        axis = -axis

    return center, axis


def plane_basis(
    axis: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    reference = np.asarray(
        [1.0, 0.0, 0.0],
        dtype=float,
    )

    if abs(float(np.dot(axis, reference))) > 0.9:
        reference = np.asarray(
            [0.0, 1.0, 0.0],
            dtype=float,
        )

    first = normalized(
        np.cross(axis, reference)
    )

    second = normalized(
        np.cross(axis, first)
    )

    return first, second


def projected_coordinates(
    positions: np.ndarray,
    center: np.ndarray,
    axis: np.ndarray,
) -> np.ndarray:
    basis_1, basis_2 = plane_basis(axis)

    displaced = positions - center

    return np.column_stack(
        (
            displaced @ basis_1,
            displaced @ basis_2,
        )
    )


def minimum_width(points_2d: np.ndarray) -> float:
    angles = np.linspace(
        0.0,
        math.pi,
        7200,
        endpoint=False,
    )

    minimum = math.inf

    for angle in angles:
        direction = np.asarray(
            [
                math.cos(angle),
                math.sin(angle),
            ],
            dtype=float,
        )

        projections = points_2d @ direction

        width = float(
            np.max(projections)
            - np.min(projections)
        )

        minimum = min(
            minimum,
            width,
        )

    return minimum


def opposite_pair_metrics(
    node_ids: list[str],
    points_2d: np.ndarray,
) -> tuple[
    float,
    float,
    list[dict[str, Any]],
]:
    angles = np.arctan2(
        points_2d[:, 1],
        points_2d[:, 0],
    )

    radii = np.linalg.norm(
        points_2d,
        axis=1,
    )

    details = []

    for first_index in range(
        len(node_ids)
    ):
        best_second = None
        best_angular_error = math.inf

        for second_index in range(
            len(node_ids)
        ):
            if first_index == second_index:
                continue

            delta = abs(
                math.atan2(
                    math.sin(
                        angles[second_index]
                        - angles[first_index]
                        - math.pi
                    ),
                    math.cos(
                        angles[second_index]
                        - angles[first_index]
                        - math.pi
                    ),
                )
            )

            if delta < best_angular_error:
                best_angular_error = delta
                best_second = second_index

        if best_second is None:
            raise RuntimeError(
                "Could not resolve opposite node."
            )

        distance = float(
            np.linalg.norm(
                points_2d[first_index]
                - points_2d[best_second]
            )
        )

        details.append(
            {
                "first_node": node_ids[first_index],
                "second_node": node_ids[best_second],
                "first_radius_nm": float(
                    radii[first_index]
                ),
                "second_radius_nm": float(
                    radii[best_second]
                ),
                "opposition_error_deg": math.degrees(
                    best_angular_error
                ),
                "pair_distance_nm": distance,
            }
        )

    unique_pairs = {}

    for row in details:
        key = tuple(
            sorted(
                (
                    row["first_node"],
                    row["second_node"],
                )
            )
        )

        if key not in unique_pairs:
            unique_pairs[key] = row

    pair_distances = [
        float(row["pair_distance_nm"])
        for row in unique_pairs.values()
    ]

    return (
        min(pair_distances),
        float(np.mean(pair_distances)),
        list(unique_pairs.values()),
    )


def extract_aperture_code(path: Path) -> str:
    require_file(path)

    lines = path.read_text(
        encoding="utf-8",
    ).splitlines()

    hits = [
        index
        for index, line in enumerate(lines)
        if re.search(
            r"aperture|inner_radii|inner_radius|diameter",
            line,
            re.IGNORECASE,
        )
    ]

    blocks = []

    for hit in hits:
        start = max(
            0,
            hit - 8,
        )

        stop = min(
            len(lines),
            hit + 14,
        )

        block = "\n".join(
            f"{index + 1:06d}: {lines[index]}"
            for index in range(
                start,
                stop,
            )
        )

        blocks.append(
            block
        )

    return "\n\n".join(
        blocks
    )


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        SCRIPT3J,
        SCRIPT3K,
        SCRIPT3O,
        COORD3K,
        COORD3O,
        SUMMARY3J,
        SUMMARY3K,
        SUMMARY3O,
    ):
        require_file(required)

    summary3j = read_one(
        SUMMARY3J
    )

    summary3k = read_one(
        SUMMARY3K
    )

    summary3o = read_one(
        SUMMARY3O
    )

    coordinates3k_rows = read_rows(
        COORD3K
    )

    coordinates3o_rows = read_rows(
        COORD3O
    )

    coordinates3k = {
        row["node_id"]: np.asarray(
            [
                parse_float(row, "x_nm"),
                parse_float(row, "y_nm"),
                parse_float(row, "z_nm"),
            ],
            dtype=float,
        )
        for row in coordinates3k_rows
        if row["element"] != "H"
    }

    coordinates3o = {
        row["node_id"]: np.asarray(
            [
                parse_float(row, "x_nm"),
                parse_float(row, "y_nm"),
                parse_float(row, "z_nm"),
            ],
            dtype=float,
        )
        for row in coordinates3o_rows
    }

    fixed_common = (
        set(coordinates3k)
        & set(coordinates3o)
    )

    maximum_fixed_coordinate_difference = max(
        float(
            np.linalg.norm(
                coordinates3k[node_id]
                - coordinates3o[node_id]
            )
        )
        for node_id in fixed_common
    )

    parent_ids = [
        row["node_id"]
        for row in coordinates3o_rows
        if row["node_type"] == "PARENT_HBN"
    ]

    parent_positions = np.asarray(
        [
            coordinates3o[node_id]
            for node_id in parent_ids
        ],
        dtype=float,
    )

    tube_center, tube_axis = principal_axis(
        parent_positions
    )

    metric_rows = []
    pair_rows = []

    for end in (
        "LOWER",
        "UPPER",
    ):
        inner_ids = [
            row["node_id"]
            for row in coordinates3o_rows
            if (
                row["end"] == end
                and row["node_type"]
                == "ANNULUS_INNER_BOUNDARY"
            )
        ]

        annulus_ids = [
            row["node_id"]
            for row in coordinates3o_rows
            if (
                row["end"] == end
                and row["node_type"]
                in {
                    "ANNULUS_INTERIOR",
                    "ANNULUS_OUTER_BOUNDARY",
                    "ANNULUS_INNER_BOUNDARY",
                }
            )
        ]

        inner_positions = np.asarray(
            [
                coordinates3o[node_id]
                for node_id in inner_ids
            ],
            dtype=float,
        )

        annulus_positions = np.asarray(
            [
                coordinates3o[node_id]
                for node_id in annulus_ids
            ],
            dtype=float,
        )

        annulus_center = np.mean(
            annulus_positions,
            axis=0,
        )

        projected = projected_coordinates(
            inner_positions,
            annulus_center,
            tube_axis,
        )

        radii = np.linalg.norm(
            projected,
            axis=1,
        )

        radial_diameter = (
            2.0
            * float(
                np.min(radii)
            )
        )

        minimum_pair, mean_pair, details = (
            opposite_pair_metrics(
                inner_ids,
                projected,
            )
        )

        width = minimum_width(
            projected
        )

        for row in details:
            pair_rows.append(
                {
                    "end": end,
                    **row,
                }
            )

        previous_value = (
            parse_float(
                summary3k,
                (
                    "lower_nuclear_aperture_diameter_nm"
                    if end == "LOWER"
                    else "upper_nuclear_aperture_diameter_nm"
                ),
            )
            if (
                (
                    "lower_nuclear_aperture_diameter_nm"
                    if end == "LOWER"
                    else "upper_nuclear_aperture_diameter_nm"
                )
                in summary3k
            )
            else math.nan
        )

        for metric_name, value in (
            (
                "TWO_TIMES_MINIMUM_RADIAL_DISTANCE",
                radial_diameter,
            ),
            (
                "MINIMUM_OPPOSITE_PAIR_DISTANCE",
                minimum_pair,
            ),
            (
                "MEAN_OPPOSITE_PAIR_DISTANCE",
                mean_pair,
            ),
            (
                "MINIMUM_PROJECTED_WIDTH",
                width,
            ),
        ):
            metric_rows.append(
                {
                    "end": end,
                    "metric": metric_name,
                    "diameter_nm": value,
                    "target_nm": TARGET_APERTURE_NM,
                    "relative_error": abs(
                        value
                        - TARGET_APERTURE_NM
                    ) / TARGET_APERTURE_NM,
                    "within10percent": abs(
                        value
                        - TARGET_APERTURE_NM
                    ) / TARGET_APERTURE_NM
                    <= MAX_ACCEPTED_RELATIVE_ERROR,
                }
            )

    write_rows(
        METRICS,
        metric_rows,
    )

    write_rows(
        PAIR_DETAILS,
        pair_rows,
    )

    source_text = "\n\n".join(
        (
            "======================================================================\n"
            f"SOURCE={relative(path)}\n"
            "======================================================================\n"
            f"{extract_aperture_code(path)}"
        )
        for path in (
            SCRIPT3J,
            SCRIPT3K,
            SCRIPT3O,
        )
    )

    SOURCE_DEFINITIONS.write_text(
        source_text
        + "\n",
        encoding="utf-8",
    )

    grouped = {}

    for row in metric_rows:
        grouped[
            (
                row["end"],
                row["metric"],
            )
        ] = float(
            row["diameter_nm"]
        )

    lower_upper_differences = {
        metric: abs(
            grouped[
                (
                    "LOWER",
                    metric,
                )
            ]
            - grouped[
                (
                    "UPPER",
                    metric,
                )
            ]
        )
        for metric in {
            row["metric"]
            for row in metric_rows
        }
    }

    physically_relevant_pass_metrics = [
        row
        for row in metric_rows
        if (
            row["metric"]
            in {
                "MINIMUM_OPPOSITE_PAIR_DISTANCE",
                "MINIMUM_PROJECTED_WIDTH",
            }
            and bool(
                row["within10percent"]
            )
        )
    ]

    radial_metric_rows = [
        row
        for row in metric_rows
        if row["metric"]
        == "TWO_TIMES_MINIMUM_RADIAL_DISTANCE"
    ]

    radial_metric_fails = all(
        not bool(row["within10percent"])
        for row in radial_metric_rows
    )

    alternative_metric_passes_both_ends = any(
        all(
            bool(row["within10percent"])
            for row in metric_rows
            if row["metric"] == metric
        )
        for metric in (
            "MINIMUM_OPPOSITE_PAIR_DISTANCE",
            "MINIMUM_PROJECTED_WIDTH",
        )
    )

    gates = {
        "Gate3O_has_expected_aperture_only_review_decision": (
            summary3o.get(
                "decision"
            )
            == EXPECTED_GATE3O_DECISION
            and summary3o.get(
                "failed_gates"
            )
            == "aperture_errors_are_within10percent"
        ),
        "fixed_annulus_coordinates_are_unchanged_from_Gate3K": (
            maximum_fixed_coordinate_difference
            <= COORDINATE_TOLERANCE_NM
        ),
        "Gate3O_radial_metric_fails_both_ends": (
            radial_metric_fails
        ),
        "at_least_one_geometric_free_aperture_metric_passes_both_ends": (
            alternative_metric_passes_both_ends
        ),
        "all_aperture_metrics_are_lower_upper_symmetric": all(
            difference
            <= 1.0e-9
            for difference
            in lower_upper_differences.values()
        ),
        "heavy_geometry_other_than_aperture_remains_accepted": (
            parse_float(
                summary3o,
                "heavy_heavy_clash_count",
            )
            == 0.0
            and parse_float(
                summary3o,
                "critical_angle_minimum_deg",
            )
            >= 70.0
            and parse_float(
                summary3o,
                "maximum_BN_bond_deviation_nm",
            )
            <= 0.003
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
        PASS_DECISION
        if accepted
        else FAIL_DECISION
    )

    preferred_metric = "NONE"

    for metric in (
        "MINIMUM_OPPOSITE_PAIR_DISTANCE",
        "MINIMUM_PROJECTED_WIDTH",
    ):
        rows = [
            row
            for row in metric_rows
            if row["metric"] == metric
        ]

        if rows and all(
            bool(row["within10percent"])
            for row in rows
        ):
            preferred_metric = metric
            break

    required_next_step = (
        "PATCH_AND_REVALIDATE_R2_HEAVY_EMBEDDING_"
        "WITH_CONSISTENT_APERTURE_METRIC"
        if accepted
        else
        "REVIEW_R2_APERTURE_METRIC_DEFINITION_AND_TARGET"
    )

    summary = {
        "decision": decision,
        "maximum_fixed_coordinate_difference_nm": (
            maximum_fixed_coordinate_difference
        ),
        "preferred_aperture_metric": (
            preferred_metric
        ),
        "lower_preferred_aperture_nm": (
            ""
            if preferred_metric == "NONE"
            else grouped[
                (
                    "LOWER",
                    preferred_metric,
                )
            ]
        ),
        "upper_preferred_aperture_nm": (
            ""
            if preferred_metric == "NONE"
            else grouped[
                (
                    "UPPER",
                    preferred_metric,
                )
            ]
        ),
        "Gate3O_radial_aperture_lower_nm": grouped[
            (
                "LOWER",
                "TWO_TIMES_MINIMUM_RADIAL_DISTANCE",
            )
        ],
        "Gate3O_radial_aperture_upper_nm": grouped[
            (
                "UPPER",
                "TWO_TIMES_MINIMUM_RADIAL_DISTANCE",
            )
        ],
        "aperture_metric_inconsistency_identified": (
            accepted
        ),
        "heavy_coordinate_embedding_retained": (
            accepted
        ),
        "hydrogen_coordinate_generation_authorized": False,
        "molecular_topology_generation_authorized": False,
        "formal_charge_assignment_authorized": False,
        "force_field_parameterization_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "QM_authorized": False,
        "failed_gates": (
            " | ".join(failed_gates)
        ),
        "required_next_step": (
            required_next_step
        ),
    }

    write_rows(
        SUMMARY,
        [summary],
    )

    write_rows(
        GATES,
        [
            {
                "gate": name,
                "pass": passed,
            }
            for name, passed in gates.items()
        ],
    )

    JSON_OUT.write_text(
        json.dumps(
            {
                "summary": summary,
                "metric_rows": metric_rows,
                "gates": gates,
                "limitations": [
                    (
                        "This gate audits aperture definitions only."
                    ),
                    (
                        "No coordinates are changed."
                    ),
                    (
                        "No hydrogen coordinates, topology, charges, "
                        "force-field parameters, minimization, MD or QM "
                        "calculation are generated."
                    ),
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    manifest_rows = [
        {
            "role": "Gate3J_script",
            "file": relative(SCRIPT3J),
            "sha256": sha256(SCRIPT3J),
        },
        {
            "role": "Gate3K_script",
            "file": relative(SCRIPT3K),
            "sha256": sha256(SCRIPT3K),
        },
        {
            "role": "Gate3O_script",
            "file": relative(SCRIPT3O),
            "sha256": sha256(SCRIPT3O),
        },
        {
            "role": "Gate3K_coordinates",
            "file": relative(COORD3K),
            "sha256": sha256(COORD3K),
        },
        {
            "role": "Gate3O_coordinates",
            "file": relative(COORD3O),
            "sha256": sha256(COORD3O),
        },
        {
            "role": "Gate3J_summary",
            "file": relative(SUMMARY3J),
            "sha256": sha256(SUMMARY3J),
        },
        {
            "role": "Gate3K_summary",
            "file": relative(SUMMARY3K),
            "sha256": sha256(SUMMARY3K),
        },
        {
            "role": "Gate3O_summary",
            "file": relative(SUMMARY3O),
            "sha256": sha256(SUMMARY3O),
        },
    ]

    write_rows(
        MANIFEST,
        manifest_rows,
    )

    metric_lines = "\n".join(
        (
            f"- {row['end']} / {row['metric']}: "
            f"{float(row['diameter_nm']):.9f} nm; "
            f"relative error="
            f"{float(row['relative_error']):.6f}; "
            f"pass={row['within10percent']}"
        )
        for row in metric_rows
    )

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed
        in gates.items()
    )

    REPORT.write_text(
        f"""# R2 Aperture Metric Consistency Audit

## Scope

This gate investigates the isolated aperture failure reported by Gate 3O.
No coordinates are modified.

## Coordinate consistency

- Maximum fixed-coordinate difference between Gate 3K and Gate 3O:
  **{maximum_fixed_coordinate_difference:.12e} nm**

## Compared aperture metrics

{metric_lines}

## Preferred metric

- Preferred metric:
  **{preferred_metric}**
- Lower value:
  **{summary['lower_preferred_aperture_nm']} nm**
- Upper value:
  **{summary['upper_preferred_aperture_nm']} nm**

## Gates

{gate_lines}

## Decision

- Decision: **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Heavy coordinate embedding retained:
  **{'YES' if accepted else 'NO'}**
- Hydrogen coordinate generation authorized:
  **NO**
- Molecular topology generation authorized:
  **NO**
- Energy minimization authorized:
  **NO**
- MD authorized:
  **NO**
- QM authorized:
  **NO**
- Required next step:
  `{required_next_step}`
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 aperture metric consistency "
        "audit completed."
    )

    print(
        "Maximum fixed-coordinate difference Gate3K/Gate3O: "
        f"{maximum_fixed_coordinate_difference:.12e} nm"
    )

    for row in metric_rows:
        print(
            f"{row['end']} {row['metric']} / "
            "diameter / relative error / pass: "
            f"{float(row['diameter_nm']):.9f}/"
            f"{float(row['relative_error']):.9f}/"
            f"{row['within10percent']}"
        )

    print(
        f"Preferred aperture metric: {preferred_metric}"
    )

    print(
        "Aperture metric inconsistency identified: "
        f"{accepted}"
    )

    print(
        "Heavy coordinate embedding retained: "
        f"{accepted}"
    )

    print(
        f"Decision: {decision}"
    )

    print(
        "Failed gates: "
        + (
            "NONE"
            if not failed_gates
            else " | ".join(failed_gates)
        )
    )

    print(
        "Hydrogen coordinate generation authorized: NO"
    )

    print(
        "Molecular topology generation authorized: NO"
    )

    print(
        "Formal charge assignment authorized: NO"
    )

    print(
        "Force-field parameterization authorized: NO"
    )

    print(
        "Energy minimization authorized: NO"
    )

    print(
        "MD authorized: NO"
    )

    print(
        "QM authorized: NO"
    )

    print(
        f"Required next step: {required_next_step}"
    )

    for path in (
        METRICS,
        PAIR_DETAILS,
        SOURCE_DEFINITIONS,
        SUMMARY,
        GATES,
        JSON_OUT,
        MANIFEST,
        REPORT,
    ):
        print(
            f"Wrote: {relative(path)}"
        )


if __name__ == "__main__":
    main()
