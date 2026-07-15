#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/phase1A/day024_chemical_end_rim_design"

GATE3M = BASE / "16_r2_selected_full_density_longer_bn_bridge_graph"
GATE3P = BASE / "23_r2_four_atom_hydrogen_coordinate_embedding"
GATE3P1B = BASE / "26_r2_complete_end_symmetry_correspondence"

OUT = BASE / "27_r2_inner_h_reflected_direction_diagnostic"

GRAPH_NODES = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

GRAPH_EDGES = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_edges.csv"
)

FULL_COORDINATES = (
    GATE3P
    / "r2_selected_four_atom_full_coordinates.csv"
)

H_PAIRS = (
    GATE3P1B
    / "r2_complete_lower_upper_hydrogen_pairs.csv"
)

TRANSFORM = (
    GATE3P1B
    / "r2_selected_annulus_symmetry_transform.csv"
)

CORRESPONDENCE_SUMMARY = (
    GATE3P1B
    / "r2_complete_end_symmetry_correspondence_summary.csv"
)

CANDIDATES = (
    OUT
    / "r2_inner_h_reflected_direction_candidates.csv"
)

SCENARIO_SUMMARY = (
    OUT
    / "r2_inner_h_reflected_direction_scenario_summary.csv"
)

SUMMARY = (
    OUT
    / "r2_inner_h_reflected_direction_diagnostic_summary.csv"
)

GATES = (
    OUT
    / "r2_inner_h_reflected_direction_diagnostic_gates.csv"
)

JSON_OUT = (
    OUT
    / "r2_inner_h_reflected_direction_diagnostic.json"
)

REPORT = (
    OUT
    / "R2_INNER_H_REFLECTED_DIRECTION_DIAGNOSTIC_DAY024.md"
)

EXPECTED_CORRESPONDENCE_DECISION = (
    "R2_COMPLETE_END_SYMMETRY_CORRESPONDENCE_VALIDATED"
)

PASS_DECISION = (
    "R2_INNER_H_REFLECTED_DIRECTION_REFINEMENT_PATH_IDENTIFIED"
)

REVIEW_DECISION = (
    "R2_INNER_H_REQUIRES_PAIRED_ORIENTATION_OPTIMIZATION"
)

BH_NM = 0.119
NH_NM = 0.101

MIN_H_HEAVY_NM = 0.070
MIN_H_H_NM = 0.060

MIN_H_ANGLE_DEG = 70.0
MAX_H_ANGLE_DEG = 175.0

MAX_APERTURE_ASYMMETRY_NM = 0.010


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty required file: {path}"
        )


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
    value = float(row[key])

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


def angle_degrees(
    first: np.ndarray,
    center: np.ndarray,
    second: np.ndarray,
) -> float:
    first_vector = normalized(
        first - center
    )

    second_vector = normalized(
        second - center
    )

    cosine = float(
        np.clip(
            np.dot(
                first_vector,
                second_vector,
            ),
            -1.0,
            1.0,
        )
    )

    return math.degrees(
        math.acos(cosine)
    )


def target_XH(element: str) -> float:
    if element == "B":
        return BH_NM

    if element == "N":
        return NH_NM

    raise RuntimeError(
        f"Unsupported attached-heavy element: {element}"
    )


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        GRAPH_NODES,
        GRAPH_EDGES,
        FULL_COORDINATES,
        H_PAIRS,
        TRANSFORM,
        CORRESPONDENCE_SUMMARY,
    ):
        require_file(required)

    correspondence_summary = read_one(
        CORRESPONDENCE_SUMMARY
    )

    if correspondence_summary.get(
        "decision"
    ) != EXPECTED_CORRESPONDENCE_DECISION:
        raise RuntimeError(
            "Complete end correspondence is not accepted."
        )

    node_rows = read_rows(
        GRAPH_NODES
    )

    edge_rows = read_rows(
        GRAPH_EDGES
    )

    coordinate_rows = read_rows(
        FULL_COORDINATES
    )

    H_pair_rows = read_rows(
        H_PAIRS
    )

    transform_rows = read_rows(
        TRANSFORM
    )

    if len(transform_rows) != 3:
        raise RuntimeError(
            "Expected three transform rows."
        )

    transform_rows.sort(
        key=lambda row: int(row["row"])
    )

    rotation = np.asarray(
        [
            [
                parse_float(row, "R0"),
                parse_float(row, "R1"),
                parse_float(row, "R2"),
            ]
            for row in transform_rows
        ],
        dtype=float,
    )

    translation = np.asarray(
        [
            parse_float(
                row,
                "translation_nm",
            )
            for row in transform_rows
        ],
        dtype=float,
    )

    inverse_rotation = rotation.T

    if abs(
        float(
            np.linalg.det(rotation)
        )
        + 1.0
    ) > 1.0e-9:
        raise RuntimeError(
            "Selected transformation is not the expected reflection."
        )

    nodes = {
        row["node_id"]: row
        for row in node_rows
    }

    positions = {
        row["node_id"]: np.asarray(
            [
                parse_float(row, "x_nm"),
                parse_float(row, "y_nm"),
                parse_float(row, "z_nm"),
            ],
            dtype=float,
        )
        for row in coordinate_rows
    }

    adjacency = {
        node_id: set()
        for node_id in nodes
    }

    for row in edge_rows:
        first = row["source_node"]
        second = row["target_node"]

        adjacency[first].add(second)
        adjacency[second].add(first)

    inner_pairs = [
        row
        for row in H_pair_rows
        if row["hydrogen_role"]
        == "ANNULUS_INNER_PASSIVANT_H"
    ]

    if len(inner_pairs) != 12:
        raise RuntimeError(
            f"Expected 12 inner-H pairs; found {len(inner_pairs)}"
        )

    heavy_ids = sorted(
        node_id
        for node_id, row in nodes.items()
        if row["element"] != "H"
    )

    H_ids = sorted(
        node_id
        for node_id, row in nodes.items()
        if row["element"] == "H"
    )

    annulus_center_by_end = {}

    parent_ids = [
        node_id
        for node_id in heavy_ids
        if nodes[node_id]["node_type"]
        == "PARENT_HBN"
    ]

    parent_positions = np.asarray(
        [
            positions[node_id]
            for node_id in parent_ids
        ],
        dtype=float,
    )

    parent_center = np.mean(
        parent_positions,
        axis=0,
    )

    centered = (
        parent_positions
        - parent_center
    )

    values, vectors = np.linalg.eigh(
        centered.T
        @ centered
    )

    tube_axis = normalized(
        vectors[
            :,
            int(
                np.argmax(values)
            ),
        ]
    )

    if tube_axis[2] < 0.0:
        tube_axis = -tube_axis

    for end in (
        "LOWER",
        "UPPER",
    ):
        annulus_ids = [
            node_id
            for node_id in heavy_ids
            if (
                nodes[node_id]["end"] == end
                and nodes[node_id]["node_type"]
                in {
                    "ANNULUS_INTERIOR",
                    "ANNULUS_OUTER_BOUNDARY",
                    "ANNULUS_INNER_BOUNDARY",
                }
            )
        ]

        annulus_center_by_end[end] = np.mean(
            np.asarray(
                [
                    positions[node_id]
                    for node_id in annulus_ids
                ],
                dtype=float,
            ),
            axis=0,
        )

    scenario_positions = {
        "LOWER_DRIVES_UPPER": {
            node_id: np.array(
                point,
                copy=True,
            )
            for node_id, point in positions.items()
        },
        "UPPER_DRIVES_LOWER": {
            node_id: np.array(
                point,
                copy=True,
            )
            for node_id, point in positions.items()
        },
    }

    candidate_rows = []

    for pair in inner_pairs:
        lower_H = pair["lower_H"]
        upper_H = pair["upper_H"]

        lower_heavy = pair[
            "lower_attached_heavy"
        ]

        upper_heavy = pair[
            "upper_attached_heavy"
        ]

        lower_direction = normalized(
            positions[lower_H]
            - positions[lower_heavy]
        )

        mapped_upper_direction = normalized(
            rotation
            @ lower_direction
        )

        reconstructed_upper = (
            positions[upper_heavy]
            + target_XH(
                nodes[upper_heavy]["element"]
            )
            * mapped_upper_direction
        )

        scenario_positions[
            "LOWER_DRIVES_UPPER"
        ][upper_H] = reconstructed_upper

        upper_direction = normalized(
            positions[upper_H]
            - positions[upper_heavy]
        )

        mapped_lower_direction = normalized(
            inverse_rotation
            @ upper_direction
        )

        reconstructed_lower = (
            positions[lower_heavy]
            + target_XH(
                nodes[lower_heavy]["element"]
            )
            * mapped_lower_direction
        )

        scenario_positions[
            "UPPER_DRIVES_LOWER"
        ][lower_H] = reconstructed_lower

        candidate_rows.extend(
            [
                {
                    "scenario": "LOWER_DRIVES_UPPER",
                    "source_H": lower_H,
                    "rebuilt_H": upper_H,
                    "source_heavy": lower_heavy,
                    "rebuilt_heavy": upper_heavy,
                    "source_heavy_element": nodes[
                        lower_heavy
                    ]["element"],
                    "rebuilt_heavy_element": nodes[
                        upper_heavy
                    ]["element"],
                    "rebuilt_target_XH_nm": target_XH(
                        nodes[upper_heavy]["element"]
                    ),
                    "x_nm": reconstructed_upper[0],
                    "y_nm": reconstructed_upper[1],
                    "z_nm": reconstructed_upper[2],
                },
                {
                    "scenario": "UPPER_DRIVES_LOWER",
                    "source_H": upper_H,
                    "rebuilt_H": lower_H,
                    "source_heavy": upper_heavy,
                    "rebuilt_heavy": lower_heavy,
                    "source_heavy_element": nodes[
                        upper_heavy
                    ]["element"],
                    "rebuilt_heavy_element": nodes[
                        lower_heavy
                    ]["element"],
                    "rebuilt_target_XH_nm": target_XH(
                        nodes[lower_heavy]["element"]
                    ),
                    "x_nm": reconstructed_lower[0],
                    "y_nm": reconstructed_lower[1],
                    "z_nm": reconstructed_lower[2],
                },
            ]
        )

    write_rows(
        CANDIDATES,
        candidate_rows,
    )

    scenario_rows = []

    for scenario, trial_positions in scenario_positions.items():
        XH_deviations = []
        angle_violations = 0
        minimum_angle = math.inf
        maximum_angle = -math.inf

        for pair in inner_pairs:
            for hydrogen_id, heavy_id in (
                (
                    pair["lower_H"],
                    pair["lower_attached_heavy"],
                ),
                (
                    pair["upper_H"],
                    pair["upper_attached_heavy"],
                ),
            ):
                distance = float(
                    np.linalg.norm(
                        trial_positions[hydrogen_id]
                        - trial_positions[heavy_id]
                    )
                )

                target = target_XH(
                    nodes[heavy_id]["element"]
                )

                XH_deviations.append(
                    abs(
                        distance
                        - target
                    )
                )

                heavy_neighbors = [
                    node_id
                    for node_id in adjacency[heavy_id]
                    if nodes[node_id]["element"] != "H"
                ]

                if len(heavy_neighbors) != 2:
                    raise RuntimeError(
                        f"{heavy_id}: expected two heavy neighbors."
                    )

                for neighbor_id in heavy_neighbors:
                    value = angle_degrees(
                        trial_positions[neighbor_id],
                        trial_positions[heavy_id],
                        trial_positions[hydrogen_id],
                    )

                    minimum_angle = min(
                        minimum_angle,
                        value,
                    )

                    maximum_angle = max(
                        maximum_angle,
                        value,
                    )

                    if (
                        value < MIN_H_ANGLE_DEG
                        or value > MAX_H_ANGLE_DEG
                    ):
                        angle_violations += 1

        H_heavy_minimum = math.inf
        H_heavy_clashes = 0

        for hydrogen_id in H_ids:
            attached_heavy_neighbors = [
                neighbor
                for neighbor in adjacency[hydrogen_id]
                if nodes[neighbor]["element"] != "H"
            ]

            if len(attached_heavy_neighbors) != 1:
                raise RuntimeError(
                    f"{hydrogen_id}: invalid attachment."
                )

            attached = attached_heavy_neighbors[0]

            excluded = {
                attached,
                *[
                    neighbor
                    for neighbor in adjacency[attached]
                    if nodes[neighbor]["element"] != "H"
                ],
            }

            for heavy_id in heavy_ids:
                if heavy_id in excluded:
                    continue

                distance = float(
                    np.linalg.norm(
                        trial_positions[hydrogen_id]
                        - trial_positions[heavy_id]
                    )
                )

                H_heavy_minimum = min(
                    H_heavy_minimum,
                    distance,
                )

                if distance < MIN_H_HEAVY_NM:
                    H_heavy_clashes += 1

        H_H_minimum = math.inf
        H_H_clashes = 0

        for first_index in range(
            len(H_ids)
        ):
            for second_index in range(
                first_index + 1,
                len(H_ids),
            ):
                distance = float(
                    np.linalg.norm(
                        trial_positions[
                            H_ids[first_index]
                        ]
                        - trial_positions[
                            H_ids[second_index]
                        ]
                    )
                )

                H_H_minimum = min(
                    H_H_minimum,
                    distance,
                )

                if distance < MIN_H_H_NM:
                    H_H_clashes += 1

        aperture_by_end = {}

        for end in (
            "LOWER",
            "UPPER",
        ):
            inner_H_ids = [
                pair[
                    "lower_H"
                    if end == "LOWER"
                    else "upper_H"
                ]
                for pair in inner_pairs
            ]

            center = annulus_center_by_end[
                end
            ]

            radii = []

            for hydrogen_id in inner_H_ids:
                displacement = (
                    trial_positions[hydrogen_id]
                    - center
                )

                displacement -= (
                    np.dot(
                        displacement,
                        tube_axis,
                    )
                    * tube_axis
                )

                radii.append(
                    float(
                        np.linalg.norm(
                            displacement
                        )
                    )
                )

            aperture_by_end[end] = (
                2.0
                * min(radii)
            )

        asymmetry = abs(
            aperture_by_end["LOWER"]
            - aperture_by_end["UPPER"]
        )

        scenario_rows.append(
            {
                "scenario": scenario,
                "maximum_XH_deviation_nm": max(
                    XH_deviations
                ),
                "minimum_H_angle_deg": (
                    minimum_angle
                ),
                "maximum_H_angle_deg": (
                    maximum_angle
                ),
                "H_angle_violation_count": (
                    angle_violations
                ),
                "minimum_H_heavy_nm": (
                    H_heavy_minimum
                ),
                "H_heavy_clash_count": (
                    H_heavy_clashes
                ),
                "minimum_H_H_nm": (
                    H_H_minimum
                ),
                "H_H_clash_count": (
                    H_H_clashes
                ),
                "lower_inner_H_aperture_nm": (
                    aperture_by_end["LOWER"]
                ),
                "upper_inner_H_aperture_nm": (
                    aperture_by_end["UPPER"]
                ),
                "aperture_asymmetry_nm": (
                    asymmetry
                ),
                "passes_all_static_gates": (
                    max(XH_deviations) <= 1.0e-12
                    and angle_violations == 0
                    and H_heavy_clashes == 0
                    and H_H_clashes == 0
                    and asymmetry
                    <= MAX_APERTURE_ASYMMETRY_NM
                ),
            }
        )

    write_rows(
        SCENARIO_SUMMARY,
        scenario_rows,
    )

    passing_scenarios = [
        row
        for row in scenario_rows
        if bool(
            row[
                "passes_all_static_gates"
            ]
        )
    ]

    selected_scenario = (
        min(
            passing_scenarios,
            key=lambda row: (
                float(
                    row[
                        "aperture_asymmetry_nm"
                    ]
                ),
                -float(
                    row[
                        "minimum_H_heavy_nm"
                    ]
                ),
                -float(
                    row[
                        "minimum_H_H_nm"
                    ]
                ),
            ),
        )["scenario"]
        if passing_scenarios
        else "NONE"
    )

    gates = {
        "validated_reflection_was_loaded": (
            abs(
                float(
                    np.linalg.det(rotation)
                )
                + 1.0
            )
            <= 1.0e-9
        ),
        "12_inner_H_pairs_were_tested": (
            len(inner_pairs) == 12
        ),
        "both_directional_reconstruction_scenarios_were_tested": (
            len(scenario_rows) == 2
        ),
        "all_rebuilt_XH_lengths_use_receiver_element_targets": all(
            float(
                row[
                    "maximum_XH_deviation_nm"
                ]
            )
            <= 1.0e-12
            for row in scenario_rows
        ),
        "at_least_one_direct_reflected_direction_scenario_passes": (
            len(
                passing_scenarios
            )
            >= 1
        ),
        "coordinates_are_not_applied_to_Gate3P_structure": True,
        "no_topology_minimization_MD_or_QM_is_generated": True,
    }

    failed_gates = [
        name
        for name, passed in gates.items()
        if not passed
    ]

    direct_path_available = (
        len(
            passing_scenarios
        )
        >= 1
    )

    decision = (
        PASS_DECISION
        if direct_path_available
        else REVIEW_DECISION
    )

    required_next_step = (
        "APPLY_AND_VALIDATE_SELECTED_INNER_H_"
        "REFLECTED_DIRECTION_SCENARIO"
        if direct_path_available
        else
        "OPTIMIZE_INNER_H_ORIENTATIONS_AS_12_"
        "COUPLED_LOWER_UPPER_PAIRS"
    )

    summary = {
        "decision": decision,
        "inner_H_pairs_tested": len(
            inner_pairs
        ),
        "selected_scenario": (
            selected_scenario
        ),
        "passing_scenario_count": len(
            passing_scenarios
        ),
        "coordinates_modified": False,
        "molecular_topology_generated": False,
        "formal_charges_assigned": False,
        "force_field_parameters_assigned": False,
        "energy_minimized": False,
        "MD_performed": False,
        "QM_performed": False,
        "failed_gates": (
            " | ".join(
                failed_gates
            )
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
                "scenario_summaries": scenario_rows,
                "gates": gates,
                "limitations": [
                    (
                        "This is a diagnostic reconstruction of the "
                        "24 inner-rim hydrogens only."
                    ),
                    (
                        "No coordinates are applied to the accepted "
                        "Gate 3P structure."
                    ),
                    (
                        "No topology, charges, parameterization, "
                        "minimization, MD or QM calculation is generated."
                    ),
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    scenario_lines = "\n".join(
        (
            f"- {row['scenario']}: "
            f"aperture L/U/asymmetry="
            f"{float(row['lower_inner_H_aperture_nm']):.9f}/"
            f"{float(row['upper_inner_H_aperture_nm']):.9f}/"
            f"{float(row['aperture_asymmetry_nm']):.9f} nm; "
            f"H-heavy clashes={row['H_heavy_clash_count']}; "
            f"H-H clashes={row['H_H_clash_count']}; "
            f"angle violations={row['H_angle_violation_count']}; "
            f"pass={row['passes_all_static_gates']}"
        )
        for row in scenario_rows
    )

    REPORT.write_text(
        f"""# R2 Inner-H Reflected-Direction Diagnostic

## Scope

The validated inner-boundary reflection is applied to X-H directions,
not directly to hydrogen coordinates. Each rebuilt hydrogen uses the
B-H or N-H length required by its receiving heavy atom.

## Scenarios

{scenario_lines}

## Decision

- Decision: **{decision}**
- Passing scenarios: **{len(passing_scenarios)}**
- Selected scenario: **{selected_scenario}**
- Coordinates modified: **NO**
- Molecular topology generated: **NO**
- Energy minimization performed: **NO**
- MD performed: **NO**
- QM performed: **NO**
- Required next step:
  `{required_next_step}`
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 inner-H reflected-direction "
        "diagnostic completed."
    )

    print(
        "Inner-H pairs tested: "
        f"{len(inner_pairs)}"
    )

    for row in scenario_rows:
        print(
            f"{row['scenario']} "
            "XHdev/angles/H-heavy/H-H/aperture-L/U/asym/pass: "
            f"{float(row['maximum_XH_deviation_nm']):.12e}/"
            f"{row['H_angle_violation_count']}/"
            f"{row['H_heavy_clash_count']}/"
            f"{row['H_H_clash_count']}/"
            f"{float(row['lower_inner_H_aperture_nm']):.9f}/"
            f"{float(row['upper_inner_H_aperture_nm']):.9f}/"
            f"{float(row['aperture_asymmetry_nm']):.9f}/"
            f"{row['passes_all_static_gates']}"
        )

    print(
        "Passing scenarios / selected scenario: "
        f"{len(passing_scenarios)}/"
        f"{selected_scenario}"
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
        "Coordinates modified: NO"
    )

    print(
        "Molecular topology generated: NO"
    )

    print(
        "Energy minimization performed: NO"
    )

    print(
        "MD performed: NO"
    )

    print(
        "QM performed: NO"
    )

    print(
        f"Required next step: {required_next_step}"
    )

    for path in (
        CANDIDATES,
        SCENARIO_SUMMARY,
        SUMMARY,
        GATES,
        JSON_OUT,
        REPORT,
    ):
        print(
            f"Wrote: {path.relative_to(ROOT)}"
        )


if __name__ == "__main__":
    main()
