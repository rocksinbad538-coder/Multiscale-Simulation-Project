#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/phase1A/day024_chemical_end_rim_design"

GATE3M = (
    BASE
    / "16_r2_selected_full_density_longer_bn_bridge_graph"
)

GATE3P = (
    BASE
    / "23_r2_four_atom_hydrogen_coordinate_embedding"
)

GATE3P1B = (
    BASE
    / "26_r2_complete_end_symmetry_correspondence"
)

GATE3P1C = (
    BASE
    / "27_r2_inner_h_reflected_direction_diagnostic"
)

GRAPH_NODES = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

GRAPH_EDGES = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_edges.csv"
)

SOURCE_COORDINATES = (
    GATE3P
    / "r2_selected_four_atom_full_coordinates.csv"
)

SOURCE_SUMMARY = (
    GATE3P
    / "r2_selected_four_atom_hydrogen_embedding_summary.csv"
)

H_PAIRS = (
    GATE3P1B
    / "r2_complete_lower_upper_hydrogen_pairs.csv"
)

CORRESPONDENCE_SUMMARY = (
    GATE3P1B
    / "r2_complete_end_symmetry_correspondence_summary.csv"
)

DIAGNOSTIC_CANDIDATES = (
    GATE3P1C
    / "r2_inner_h_reflected_direction_candidates.csv"
)

DIAGNOSTIC_SUMMARY = (
    GATE3P1C
    / "r2_inner_h_reflected_direction_diagnostic_summary.csv"
)

OUT = (
    BASE
    / "28_r2_inner_h_reflected_direction_refinement"
)

REFINED_COORDINATES = (
    OUT
    / "r2_selected_four_atom_refined_full_coordinates.csv"
)

MODIFIED_COORDINATES = (
    OUT
    / "r2_inner_h_modified_coordinates.csv"
)

BOND_SUMMARY = (
    OUT
    / "r2_refined_hydrogen_bond_summary.csv"
)

ANGLE_SUMMARY = (
    OUT
    / "r2_refined_hydrogen_angle_summary.csv"
)

CONTACT_SUMMARY = (
    OUT
    / "r2_refined_hydrogen_contact_summary.csv"
)

END_SUMMARY = (
    OUT
    / "r2_refined_hydrogen_end_summary.csv"
)

SUMMARY = (
    OUT
    / "r2_inner_h_refinement_summary.csv"
)

GATES = (
    OUT
    / "r2_inner_h_refinement_gates.csv"
)

JSON_OUT = (
    OUT
    / "r2_inner_h_refinement.json"
)

MANIFEST = (
    OUT
    / "r2_inner_h_refinement_manifest.csv"
)

XYZ_OUT = (
    OUT
    / "r2_selected_four_atom_refined_full_embedding.xyz"
)

PDB_OUT = (
    OUT
    / "r2_selected_four_atom_refined_full_embedding.pdb"
)

REPORT = (
    OUT
    / "R2_INNER_H_REFLECTED_DIRECTION_REFINEMENT_DAY024.md"
)

EXPECTED_CORRESPONDENCE_DECISION = (
    "R2_COMPLETE_END_SYMMETRY_CORRESPONDENCE_VALIDATED"
)

EXPECTED_DIAGNOSTIC_DECISION = (
    "R2_INNER_H_REFLECTED_DIRECTION_REFINEMENT_PATH_IDENTIFIED"
)

EXPECTED_SELECTED_SCENARIO = (
    "UPPER_DRIVES_LOWER"
)

PASS_DECISION = (
    "R2_SELECTED_FOUR_ATOM_BN_BRIDGE_"
    "HYDROGEN_COORDINATES_VALIDATED_AFTER_SYMMETRY_REFINEMENT"
)

REVIEW_DECISION = (
    "R2_INNER_H_REFLECTED_DIRECTION_REFINEMENT_REQUIRES_REVIEW"
)

N_HEAVY = 2112
N_H = 204
N_TOTAL = 2316
N_MODIFIED_H = 12

BH_NM = 0.119
NH_NM = 0.101

MAX_XH_DEVIATION_NM = 0.002
MIN_H_ANGLE_DEG = 70.0
MAX_H_ANGLE_DEG = 175.0
MIN_H_HEAVY_NM = 0.070
MIN_H_H_NM = 0.060
MAX_APERTURE_ASYMMETRY_NM = 0.010
COORDINATE_TOLERANCE_NM = 1.0e-12


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty required file: {path}"
        )


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
    )


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
            f"Could not parse field {key!r}"
        ) from exc

    if not math.isfinite(value):
        raise RuntimeError(
            f"Non-finite field {key!r}"
        )

    return value


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


def normalized(vector: np.ndarray) -> np.ndarray:
    norm = float(
        np.linalg.norm(
            vector
        )
    )

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
        math.acos(
            cosine
        )
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
        SOURCE_COORDINATES,
        SOURCE_SUMMARY,
        H_PAIRS,
        CORRESPONDENCE_SUMMARY,
        DIAGNOSTIC_CANDIDATES,
        DIAGNOSTIC_SUMMARY,
    ):
        require_file(required)

    source_summary = read_one(
        SOURCE_SUMMARY
    )

    correspondence_summary = read_one(
        CORRESPONDENCE_SUMMARY
    )

    diagnostic_summary = read_one(
        DIAGNOSTIC_SUMMARY
    )

    if correspondence_summary.get(
        "decision"
    ) != EXPECTED_CORRESPONDENCE_DECISION:
        raise RuntimeError(
            "Complete end correspondence is not accepted."
        )

    if diagnostic_summary.get(
        "decision"
    ) != EXPECTED_DIAGNOSTIC_DECISION:
        raise RuntimeError(
            "Reflected-direction diagnostic is not accepted."
        )

    if diagnostic_summary.get(
        "selected_scenario"
    ) != EXPECTED_SELECTED_SCENARIO:
        raise RuntimeError(
            "Unexpected selected scenario: "
            f"{diagnostic_summary.get('selected_scenario')}"
        )

    node_rows = read_rows(
        GRAPH_NODES
    )

    edge_rows = read_rows(
        GRAPH_EDGES
    )

    source_coordinate_rows = read_rows(
        SOURCE_COORDINATES
    )

    H_pair_rows = read_rows(
        H_PAIRS
    )

    candidate_rows = read_rows(
        DIAGNOSTIC_CANDIDATES
    )

    nodes = {
        row["node_id"]: row
        for row in node_rows
    }

    source_positions = {
        row["node_id"]: np.asarray(
            [
                parse_float(row, "x_nm"),
                parse_float(row, "y_nm"),
                parse_float(row, "z_nm"),
            ],
            dtype=float,
        )
        for row in source_coordinate_rows
    }

    if len(source_positions) != N_TOTAL:
        raise RuntimeError(
            f"Expected {N_TOTAL} source coordinates; "
            f"found {len(source_positions)}."
        )

    positions = {
        node_id: np.array(
            point,
            copy=True,
        )
        for node_id, point in source_positions.items()
    }

    adjacency = {
        node_id: set()
        for node_id in nodes
    }

    for row in edge_rows:
        first = row["source_node"]
        second = row["target_node"]

        adjacency[first].add(
            second
        )

        adjacency[second].add(
            first
        )

    selected_candidates = [
        row
        for row in candidate_rows
        if row["scenario"]
        == EXPECTED_SELECTED_SCENARIO
    ]

    if len(selected_candidates) != N_MODIFIED_H:
        raise RuntimeError(
            f"Expected {N_MODIFIED_H} selected candidate rows; "
            f"found {len(selected_candidates)}."
        )

    modified_rows = []
    modified_H_ids = set()

    for row in selected_candidates:
        hydrogen_id = row[
            "rebuilt_H"
        ]

        if hydrogen_id in modified_H_ids:
            raise RuntimeError(
                f"Duplicate rebuilt H: {hydrogen_id}"
            )

        if nodes[
            hydrogen_id
        ][
            "element"
        ] != "H":
            raise RuntimeError(
                f"Candidate is not H: {hydrogen_id}"
            )

        if nodes[
            hydrogen_id
        ][
            "end"
        ] != "LOWER":
            raise RuntimeError(
                f"Expected LOWER rebuilt H: {hydrogen_id}"
            )

        if nodes[
            hydrogen_id
        ][
            "node_type"
        ] != "ANNULUS_INNER_PASSIVANT_H":
            raise RuntimeError(
                f"Expected inner-rim H: {hydrogen_id}"
            )

        new_position = np.asarray(
            [
                parse_float(row, "x_nm"),
                parse_float(row, "y_nm"),
                parse_float(row, "z_nm"),
            ],
            dtype=float,
        )

        old_position = positions[
            hydrogen_id
        ]

        displacement = float(
            np.linalg.norm(
                new_position
                - old_position
            )
        )

        positions[
            hydrogen_id
        ] = new_position

        modified_H_ids.add(
            hydrogen_id
        )

        modified_rows.append(
            {
                "hydrogen_node": hydrogen_id,
                "hydrogen_role": nodes[
                    hydrogen_id
                ][
                    "node_type"
                ],
                "end": nodes[
                    hydrogen_id
                ][
                    "end"
                ],
                "source_H": row[
                    "source_H"
                ],
                "attached_heavy": row[
                    "rebuilt_heavy"
                ],
                "attached_heavy_element": row[
                    "rebuilt_heavy_element"
                ],
                "old_x_nm": old_position[0],
                "old_y_nm": old_position[1],
                "old_z_nm": old_position[2],
                "new_x_nm": new_position[0],
                "new_y_nm": new_position[1],
                "new_z_nm": new_position[2],
                "displacement_nm": displacement,
                "scenario": EXPECTED_SELECTED_SCENARIO,
            }
        )

    write_rows(
        MODIFIED_COORDINATES,
        modified_rows,
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

    heavy_coordinates_unchanged = all(
        np.linalg.norm(
            positions[node_id]
            - source_positions[node_id]
        )
        <= COORDINATE_TOLERANCE_NM
        for node_id in heavy_ids
    )

    unmodified_H_ids = [
        node_id
        for node_id in H_ids
        if node_id not in modified_H_ids
    ]

    unmodified_H_coordinates_unchanged = all(
        np.linalg.norm(
            positions[node_id]
            - source_positions[node_id]
        )
        <= COORDINATE_TOLERANCE_NM
        for node_id in unmodified_H_ids
    )

    attached_heavy_by_H = {}

    for hydrogen_id in H_ids:
        heavy_neighbors = [
            neighbor
            for neighbor in adjacency[
                hydrogen_id
            ]
            if nodes[
                neighbor
            ][
                "element"
            ]
            != "H"
        ]

        if len(heavy_neighbors) != 1:
            raise RuntimeError(
                f"{hydrogen_id}: expected one attached heavy."
            )

        attached_heavy_by_H[
            hydrogen_id
        ] = heavy_neighbors[0]

    bond_values = {
        "B-H": [],
        "N-H": [],
    }

    angle_rows = []
    angle_violations = 0
    global_minimum_angle = math.inf
    global_maximum_angle = -math.inf

    for hydrogen_id in H_ids:
        heavy_id = attached_heavy_by_H[
            hydrogen_id
        ]

        heavy_element = nodes[
            heavy_id
        ][
            "element"
        ]

        bond_type = (
            f"{heavy_element}-H"
        )

        distance = float(
            np.linalg.norm(
                positions[
                    hydrogen_id
                ]
                - positions[
                    heavy_id
                ]
            )
        )

        bond_values[
            bond_type
        ].append(
            distance
        )

        heavy_heavy_neighbors = [
            neighbor
            for neighbor in adjacency[
                heavy_id
            ]
            if nodes[
                neighbor
            ][
                "element"
            ]
            != "H"
        ]

        if len(heavy_heavy_neighbors) != 2:
            raise RuntimeError(
                f"{heavy_id}: expected two heavy neighbors; "
                f"found {len(heavy_heavy_neighbors)}."
            )

        local_angles = []

        for neighbor_id in heavy_heavy_neighbors:
            value = angle_degrees(
                positions[
                    neighbor_id
                ],
                positions[
                    heavy_id
                ],
                positions[
                    hydrogen_id
                ],
            )

            local_angles.append(
                value
            )

            global_minimum_angle = min(
                global_minimum_angle,
                value,
            )

            global_maximum_angle = max(
                global_maximum_angle,
                value,
            )

            if (
                value
                < MIN_H_ANGLE_DEG
                or value
                > MAX_H_ANGLE_DEG
            ):
                angle_violations += 1

        angle_rows.append(
            {
                "hydrogen_node": hydrogen_id,
                "attached_heavy": heavy_id,
                "attached_heavy_element": heavy_element,
                "hydrogen_role": nodes[
                    hydrogen_id
                ][
                    "node_type"
                ],
                "end": nodes[
                    hydrogen_id
                ][
                    "end"
                ],
                "angle_1_deg": local_angles[0],
                "angle_2_deg": local_angles[1],
                "minimum_angle_deg": min(
                    local_angles
                ),
                "maximum_angle_deg": max(
                    local_angles
                ),
                "violation_count": sum(
                    (
                        value
                        < MIN_H_ANGLE_DEG
                        or value
                        > MAX_H_ANGLE_DEG
                    )
                    for value in local_angles
                ),
                "coordinate_modified_in_Gate3P2": (
                    hydrogen_id
                    in modified_H_ids
                ),
            }
        )

    write_rows(
        ANGLE_SUMMARY,
        angle_rows,
    )

    bond_summary_rows = []

    maximum_XH_deviation = 0.0

    for bond_type, values in sorted(
        bond_values.items()
    ):
        array = np.asarray(
            values,
            dtype=float,
        )

        target = (
            BH_NM
            if bond_type == "B-H"
            else NH_NM
        )

        deviation = float(
            np.max(
                np.abs(
                    array
                    - target
                )
            )
        )

        maximum_XH_deviation = max(
            maximum_XH_deviation,
            deviation,
        )

        bond_summary_rows.append(
            {
                "bond_type": bond_type,
                "count": int(
                    array.size
                ),
                "target_nm": target,
                "minimum_nm": float(
                    np.min(
                        array
                    )
                ),
                "mean_nm": float(
                    np.mean(
                        array
                    )
                ),
                "maximum_nm": float(
                    np.max(
                        array
                    )
                ),
                "maximum_absolute_deviation_nm": (
                    deviation
                ),
            }
        )

    write_rows(
        BOND_SUMMARY,
        bond_summary_rows,
    )

    H_heavy_minimum = math.inf
    H_heavy_clashes = 0
    H_heavy_minimum_pair = ""

    for hydrogen_id in H_ids:
        attached = attached_heavy_by_H[
            hydrogen_id
        ]

        excluded = {
            attached,
            *[
                neighbor
                for neighbor in adjacency[
                    attached
                ]
                if nodes[
                    neighbor
                ][
                    "element"
                ]
                != "H"
            ],
        }

        for heavy_id in heavy_ids:
            if heavy_id in excluded:
                continue

            distance = float(
                np.linalg.norm(
                    positions[
                        hydrogen_id
                    ]
                    - positions[
                        heavy_id
                    ]
                )
            )

            if distance < H_heavy_minimum:
                H_heavy_minimum = distance

                H_heavy_minimum_pair = (
                    f"{hydrogen_id} | {heavy_id}"
                )

            if distance < MIN_H_HEAVY_NM:
                H_heavy_clashes += 1

    H_H_minimum = math.inf
    H_H_clashes = 0
    H_H_minimum_pair = ""

    for first_index in range(
        len(H_ids)
    ):
        for second_index in range(
            first_index + 1,
            len(H_ids),
        ):
            first_id = H_ids[
                first_index
            ]

            second_id = H_ids[
                second_index
            ]

            distance = float(
                np.linalg.norm(
                    positions[
                        first_id
                    ]
                    - positions[
                        second_id
                    ]
                )
            )

            if distance < H_H_minimum:
                H_H_minimum = distance

                H_H_minimum_pair = (
                    f"{first_id} | {second_id}"
                )

            if distance < MIN_H_H_NM:
                H_H_clashes += 1

    write_rows(
        CONTACT_SUMMARY,
        [
            {
                "contact_type": "H_HEAVY",
                "minimum_distance_nm": (
                    H_heavy_minimum
                ),
                "minimum_pair": (
                    H_heavy_minimum_pair
                ),
                "threshold_nm": (
                    MIN_H_HEAVY_NM
                ),
                "clash_count": (
                    H_heavy_clashes
                ),
            },
            {
                "contact_type": "H_H",
                "minimum_distance_nm": (
                    H_H_minimum
                ),
                "minimum_pair": (
                    H_H_minimum_pair
                ),
                "threshold_nm": (
                    MIN_H_H_NM
                ),
                "clash_count": (
                    H_H_clashes
                ),
            },
        ],
    )

    parent_ids = [
        node_id
        for node_id in heavy_ids
        if nodes[
            node_id
        ][
            "node_type"
        ]
        == "PARENT_HBN"
    ]

    parent_positions = np.asarray(
        [
            positions[
                node_id
            ]
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

    eigenvalues, eigenvectors = np.linalg.eigh(
        centered.T
        @ centered
    )

    tube_axis = normalized(
        eigenvectors[
            :,
            int(
                np.argmax(
                    eigenvalues
                )
            ),
        ]
    )

    if tube_axis[2] < 0.0:
        tube_axis = -tube_axis

    end_rows = []

    for end in (
        "LOWER",
        "UPPER",
    ):
        annulus_heavy_ids = [
            node_id
            for node_id in heavy_ids
            if (
                nodes[
                    node_id
                ][
                    "end"
                ]
                == end
                and nodes[
                    node_id
                ][
                    "node_type"
                ]
                in {
                    "ANNULUS_INTERIOR",
                    "ANNULUS_OUTER_BOUNDARY",
                    "ANNULUS_INNER_BOUNDARY",
                }
            )
        ]

        annulus_center = np.mean(
            np.asarray(
                [
                    positions[
                        node_id
                    ]
                    for node_id in annulus_heavy_ids
                ],
                dtype=float,
            ),
            axis=0,
        )

        inner_H_ids = [
            node_id
            for node_id in H_ids
            if (
                nodes[
                    node_id
                ][
                    "end"
                ]
                == end
                and nodes[
                    node_id
                ][
                    "node_type"
                ]
                == "ANNULUS_INNER_PASSIVANT_H"
            )
        ]

        radii = []

        for hydrogen_id in inner_H_ids:
            displacement = (
                positions[
                    hydrogen_id
                ]
                - annulus_center
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

        end_rows.append(
            {
                "end": end,
                "total_H": sum(
                    nodes[
                        node_id
                    ][
                        "end"
                    ]
                    == end
                    for node_id in H_ids
                ),
                "inner_H_count": len(
                    inner_H_ids
                ),
                "modified_inner_H_count": sum(
                    node_id
                    in modified_H_ids
                    for node_id in inner_H_ids
                ),
                "inner_H_nuclear_aperture_diameter_nm": (
                    2.0
                    * min(
                        radii
                    )
                ),
                "inner_H_radius_minimum_nm": min(
                    radii
                ),
                "inner_H_radius_mean_nm": float(
                    np.mean(
                        radii
                    )
                ),
                "inner_H_radius_maximum_nm": max(
                    radii
                ),
            }
        )

    write_rows(
        END_SUMMARY,
        end_rows,
    )

    lower_end = next(
        row
        for row in end_rows
        if row["end"] == "LOWER"
    )

    upper_end = next(
        row
        for row in end_rows
        if row["end"] == "UPPER"
    )

    lower_aperture = float(
        lower_end[
            "inner_H_nuclear_aperture_diameter_nm"
        ]
    )

    upper_aperture = float(
        upper_end[
            "inner_H_nuclear_aperture_diameter_nm"
        ]
    )

    aperture_asymmetry = abs(
        lower_aperture
        - upper_aperture
    )

    refined_coordinate_rows = []

    source_row_by_id = {
        row["node_id"]: row
        for row in source_coordinate_rows
    }

    for node_id in sorted(
        nodes
    ):
        source_row = source_row_by_id[
            node_id
        ]

        point = positions[
            node_id
        ]

        refined_coordinate_rows.append(
            {
                **source_row,
                "x_nm": float(
                    point[0]
                ),
                "y_nm": float(
                    point[1]
                ),
                "z_nm": float(
                    point[2]
                ),
                "coordinate_source": (
                    "GATE3P2_UPPER_DRIVES_LOWER_INNER_H_REFINEMENT"
                    if node_id in modified_H_ids
                    else source_row.get(
                        "coordinate_source",
                        "GATE3P_SOURCE_COORDINATE",
                    )
                ),
                "coordinate_modified_in_Gate3P2": (
                    node_id
                    in modified_H_ids
                ),
                "energy_minimized": False,
                "MD_relaxed": False,
            }
        )

    write_rows(
        REFINED_COORDINATES,
        refined_coordinate_rows,
    )

    gates = {
        "diagnostic_selected_UPPER_DRIVES_LOWER": (
            diagnostic_summary.get(
                "selected_scenario"
            )
            == EXPECTED_SELECTED_SCENARIO
        ),
        "2316_nodes_have_coordinates": (
            len(
                positions
            )
            == N_TOTAL
        ),
        "2112_heavy_coordinates_are_unchanged": (
            heavy_coordinates_unchanged
        ),
        "exactly_12_lower_inner_H_coordinates_were_modified": (
            len(
                modified_H_ids
            )
            == N_MODIFIED_H
            and all(
                nodes[
                    node_id
                ][
                    "end"
                ]
                == "LOWER"
                and nodes[
                    node_id
                ][
                    "node_type"
                ]
                == "ANNULUS_INNER_PASSIVANT_H"
                for node_id in modified_H_ids
            )
        ),
        "remaining_192_H_coordinates_are_unchanged": (
            len(
                unmodified_H_ids
            )
            == N_H
            - N_MODIFIED_H
            and unmodified_H_coordinates_unchanged
        ),
        "102_BH_and_102_NH_bonds_are_present": (
            len(
                bond_values[
                    "B-H"
                ]
            )
            == 102
            and len(
                bond_values[
                    "N-H"
                ]
            )
            == 102
        ),
        "all_XH_bonds_are_within_0p002nm": (
            maximum_XH_deviation
            <= MAX_XH_DEVIATION_NM
        ),
        "all_H_angles_are_within70to175deg": (
            angle_violations
            == 0
        ),
        "no_global_H_heavy_clashes": (
            H_heavy_clashes
            == 0
        ),
        "no_global_H_H_clashes": (
            H_H_clashes
            == 0
        ),
        "lower_upper_inner_H_aperture_asymmetry_is_at_most0p010nm": (
            aperture_asymmetry
            <= MAX_APERTURE_ASYMMETRY_NM
        ),
        "inner_H_aperture_is_recorded_as_geometric_proxy_only": True,
        "no_topology_charges_parameterization_minimization_MD_or_QM": True,
    }

    failed_gates = [
        name
        for name, passed
        in gates.items()
        if not passed
    ]

    accepted = (
        len(
            failed_gates
        )
        == 0
    )

    decision = (
        PASS_DECISION
        if accepted
        else REVIEW_DECISION
    )

    required_next_step = (
        "AUDIT_R2_SELECTED_FOUR_ATOM_BN_BRIDGE_"
        "CHEMICAL_REALIZABILITY_AND_PARAMETERIZATION_SCOPE"
        if accepted
        else
        "REVIEW_R2_INNER_H_REFLECTED_DIRECTION_REFINEMENT"
    )

    summary = {
        "decision": decision,
        "selected_scenario_applied": (
            EXPECTED_SELECTED_SCENARIO
        ),
        "coordinate_nodes": len(
            positions
        ),
        "heavy_coordinates_preserved": (
            heavy_coordinates_unchanged
        ),
        "modified_H_coordinates": len(
            modified_H_ids
        ),
        "unmodified_H_coordinates": len(
            unmodified_H_ids
        ),
        "B_H_bonds": len(
            bond_values[
                "B-H"
            ]
        ),
        "N_H_bonds": len(
            bond_values[
                "N-H"
            ]
        ),
        "maximum_XH_deviation_nm": (
            maximum_XH_deviation
        ),
        "minimum_H_angle_deg": (
            global_minimum_angle
        ),
        "maximum_H_angle_deg": (
            global_maximum_angle
        ),
        "H_angle_violation_count": (
            angle_violations
        ),
        "minimum_H_heavy_nm": (
            H_heavy_minimum
        ),
        "minimum_H_H_nm": (
            H_H_minimum
        ),
        "H_heavy_clash_count": (
            H_heavy_clashes
        ),
        "H_H_clash_count": (
            H_H_clashes
        ),
        "lower_inner_H_nuclear_aperture_nm": (
            lower_aperture
        ),
        "upper_inner_H_nuclear_aperture_nm": (
            upper_aperture
        ),
        "inner_H_aperture_asymmetry_nm": (
            aperture_asymmetry
        ),
        "inner_H_aperture_is_functional_5kBT_result": False,
        "candidate_is_final_parameterized_chemistry": False,
        "molecular_topology_generation_authorized": False,
        "formal_charge_assignment_authorized": False,
        "force_field_parameterization_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "QM_authorized": False,
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
        [
            summary
        ],
    )

    write_rows(
        GATES,
        [
            {
                "gate": name,
                "pass": passed,
            }
            for name, passed
            in gates.items()
        ],
    )

    JSON_OUT.write_text(
        json.dumps(
            {
                "summary": summary,
                "gates": gates,
                "modified_H_nodes": sorted(
                    modified_H_ids
                ),
                "limitations": [
                    (
                        "Only 12 LOWER inner-rim H coordinates were "
                        "replaced using the validated "
                        "UPPER_DRIVES_LOWER reflected-direction scenario."
                    ),
                    (
                        "No heavy coordinates or other hydrogen "
                        "coordinates were modified."
                    ),
                    (
                        "The inner-H aperture remains a nuclear "
                        "geometric proxy, not a functional 5 kBT result."
                    ),
                    (
                        "No topology, formal charges, force-field "
                        "parameters, minimization, MD or QM calculation "
                        "was generated."
                    ),
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    write_rows(
        MANIFEST,
        [
            {
                "role": "Gate3M_graph_nodes",
                "file": relative(
                    GRAPH_NODES
                ),
                "sha256": sha256(
                    GRAPH_NODES
                ),
            },
            {
                "role": "Gate3M_graph_edges",
                "file": relative(
                    GRAPH_EDGES
                ),
                "sha256": sha256(
                    GRAPH_EDGES
                ),
            },
            {
                "role": "Gate3P_source_coordinates",
                "file": relative(
                    SOURCE_COORDINATES
                ),
                "sha256": sha256(
                    SOURCE_COORDINATES
                ),
            },
            {
                "role": "Gate3P1b_H_pairs",
                "file": relative(
                    H_PAIRS
                ),
                "sha256": sha256(
                    H_PAIRS
                ),
            },
            {
                "role": "Gate3P1c_candidates",
                "file": relative(
                    DIAGNOSTIC_CANDIDATES
                ),
                "sha256": sha256(
                    DIAGNOSTIC_CANDIDATES
                ),
            },
            {
                "role": "Gate3P1c_summary",
                "file": relative(
                    DIAGNOSTIC_SUMMARY
                ),
                "sha256": sha256(
                    DIAGNOSTIC_SUMMARY
                ),
            },
        ],
    )

    ordered_ids = sorted(
        nodes
    )

    xyz_lines = [
        str(
            len(
                ordered_ids
            )
        ),
        (
            "R2 four-atom BN bridge static embedding after "
            "inner-H symmetry refinement; not energy minimized"
        ),
    ]

    for node_id in ordered_ids:
        point = (
            positions[
                node_id
            ]
            * 10.0
        )

        xyz_lines.append(
            f"{nodes[node_id]['element']:2s} "
            f"{point[0]: .8f} "
            f"{point[1]: .8f} "
            f"{point[2]: .8f}"
        )

    XYZ_OUT.write_text(
        "\n".join(
            xyz_lines
        )
        + "\n",
        encoding="utf-8",
    )

    pdb_lines = [
        "REMARK R2 FOUR-ATOM BN BRIDGE STATIC EMBEDDING",
        "REMARK INNER-H SYMMETRY REFINEMENT APPLIED",
        "REMARK NOT ENERGY MINIMIZED OR MD RELAXED",
    ]

    for serial, node_id in enumerate(
        ordered_ids,
        start=1,
    ):
        point = (
            positions[
                node_id
            ]
            * 10.0
        )

        row = nodes[
            node_id
        ]

        element = row[
            "element"
        ]

        residue = (
            "HBN"
            if row[
                "node_type"
            ]
            == "PARENT_HBN"
            else (
                "BR4"
                if row[
                    "node_type"
                ]
                == "ALTERNATING_BN_FOUR_ATOM_BRIDGE"
                else (
                    "HYD"
                    if element
                    == "H"
                    else "ANN"
                )
            )
        )

        chain = (
            "L"
            if row[
                "end"
            ]
            == "LOWER"
            else (
                "U"
                if row[
                    "end"
                ]
                == "UPPER"
                else "P"
            )
        )

        atom_name = (
            element
            + str(
                serial
                % 1000
            )
        )[:4]

        pdb_lines.append(
            f"ATOM  {serial:5d} "
            f"{atom_name:>4s} "
            f"{residue:>3s} "
            f"{chain:1s}"
            f"{1:4d}    "
            f"{point[0]:8.3f}"
            f"{point[1]:8.3f}"
            f"{point[2]:8.3f}"
            f"{1.00:6.2f}"
            f"{0.00:6.2f}          "
            f"{element:>2s}"
        )

    pdb_lines.append(
        "END"
    )

    PDB_OUT.write_text(
        "\n".join(
            pdb_lines
        )
        + "\n",
        encoding="utf-8",
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
        f"""# R2 Inner-H Reflected-Direction Refinement

## Applied refinement

- Scenario:
  **{EXPECTED_SELECTED_SCENARIO}**
- Modified coordinates:
  **{len(modified_H_ids)} LOWER inner-rim H**
- Preserved heavy coordinates:
  **{heavy_coordinates_unchanged}**
- Preserved remaining H coordinates:
  **{unmodified_H_coordinates_unchanged}**

## Bond and angle geometry

- B-H/N-H bonds:
  **{len(bond_values['B-H'])}/{len(bond_values['N-H'])}**
- Maximum X-H deviation:
  **{maximum_XH_deviation:.12e} nm**
- H-angle range:
  **{global_minimum_angle:.9f}–{global_maximum_angle:.9f} degrees**
- H-angle violations:
  **{angle_violations}**

## Nonbonded clearance

- Minimum H-heavy / clashes:
  **{H_heavy_minimum:.9f} nm / {H_heavy_clashes}**
- Minimum H-H / clashes:
  **{H_H_minimum:.9f} nm / {H_H_clashes}**

## Inner-H nuclear aperture proxy

- LOWER:
  **{lower_aperture:.9f} nm**
- UPPER:
  **{upper_aperture:.9f} nm**
- Asymmetry:
  **{aperture_asymmetry:.9f} nm**

This remains a nuclear geometric proxy and is not the effective
water-accessible aperture at 5 kBT.

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Molecular topology generation authorized:
  **NO**
- Formal charge assignment authorized:
  **NO**
- Force-field parameterization authorized:
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
        "Day024 R2 inner-H reflected-direction "
        "refinement completed."
    )

    print(
        "Selected scenario applied: "
        f"{EXPECTED_SELECTED_SCENARIO}"
    )

    print(
        "Coordinate nodes heavy/H/total: "
        f"{len(heavy_ids)}/"
        f"{len(H_ids)}/"
        f"{len(positions)}"
    )

    print(
        "Modified lower inner-H / unchanged remaining H: "
        f"{len(modified_H_ids)}/"
        f"{len(unmodified_H_ids)}"
    )

    print(
        "Heavy coordinates unchanged: "
        f"{heavy_coordinates_unchanged}"
    )

    print(
        "Unmodified H coordinates unchanged: "
        f"{unmodified_H_coordinates_unchanged}"
    )

    print(
        "B-H / N-H bond counts: "
        f"{len(bond_values['B-H'])}/"
        f"{len(bond_values['N-H'])}"
    )

    print(
        "Maximum X-H deviation: "
        f"{maximum_XH_deviation:.12e} nm"
    )

    print(
        "H angle minimum / maximum / violations: "
        f"{global_minimum_angle:.9f}/"
        f"{global_maximum_angle:.9f}/"
        f"{angle_violations}"
    )

    print(
        "Minimum H-heavy / H-H distances: "
        f"{H_heavy_minimum:.9f}/"
        f"{H_H_minimum:.9f} nm"
    )

    print(
        "Global H-heavy / H-H clashes: "
        f"{H_heavy_clashes}/"
        f"{H_H_clashes}"
    )

    print(
        "Inner-H aperture proxy lower/upper/asymmetry: "
        f"{lower_aperture:.9f}/"
        f"{upper_aperture:.9f}/"
        f"{aperture_asymmetry:.9f} nm"
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
        REFINED_COORDINATES,
        MODIFIED_COORDINATES,
        BOND_SUMMARY,
        ANGLE_SUMMARY,
        CONTACT_SUMMARY,
        END_SUMMARY,
        SUMMARY,
        GATES,
        JSON_OUT,
        MANIFEST,
        XYZ_OUT,
        PDB_OUT,
        REPORT,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not accepted:
        raise RuntimeError(
            "Inner-H reflected-direction refinement requires review."
        )


if __name__ == "__main__":
    main()
