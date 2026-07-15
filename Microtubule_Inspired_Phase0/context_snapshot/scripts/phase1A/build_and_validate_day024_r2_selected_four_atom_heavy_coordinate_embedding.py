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

GATE3K = BASE / "13_r2_trimer_bridge_conformer_and_h_refinement"
GATE3M = BASE / "16_r2_selected_full_density_longer_bn_bridge_graph"
GATE3N = BASE / "18_r2_selected_four_atom_exact_conformer_replay"
OUT = BASE / "19_r2_selected_four_atom_heavy_coordinate_embedding"

SOURCE_FIXED_COORDINATES = (
    GATE3K
    / "r2_trimer_bridge_refined_coordinates.csv"
)

GRAPH_NODES = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

GRAPH_EDGES = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_edges.csv"
)

GRAPH_PATHS = (
    GATE3M
    / "r2_selected_longer_bn_bridge_paths.csv"
)

GRAPH_SUMMARY = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_summary.csv"
)

EXACT_BRIDGE_COORDINATES = (
    GATE3N
    / "r2_selected_four_atom_exact_bridge_coordinates.csv"
)

REPLAY_SUMMARY = (
    GATE3N
    / "r2_selected_four_atom_exact_conformer_replay_summary.csv"
)

HEAVY_COORDINATES = (
    OUT
    / "r2_selected_four_atom_heavy_coordinates.csv"
)

BOND_LENGTHS = (
    OUT
    / "r2_selected_four_atom_heavy_bond_lengths.csv"
)

ANGLE_SUMMARY = (
    OUT
    / "r2_selected_four_atom_heavy_angle_summary.csv"
)

CONTACT_SUMMARY = (
    OUT
    / "r2_selected_four_atom_heavy_contact_summary.csv"
)

END_SUMMARY = (
    OUT
    / "r2_selected_four_atom_heavy_end_summary.csv"
)

SUMMARY = (
    OUT
    / "r2_selected_four_atom_heavy_embedding_summary.csv"
)

GATES = (
    OUT
    / "r2_selected_four_atom_heavy_embedding_gates.csv"
)

JSON_OUT = (
    OUT
    / "r2_selected_four_atom_heavy_embedding.json"
)

MANIFEST = (
    OUT
    / "r2_selected_four_atom_heavy_embedding_manifest.csv"
)

XYZ_OUT = (
    OUT
    / "r2_selected_four_atom_heavy_embedding.xyz"
)

PDB_OUT = (
    OUT
    / "r2_selected_four_atom_heavy_embedding.pdb"
)

REPORT = (
    OUT
    / "R2_SELECTED_FOUR_ATOM_HEAVY_COORDINATE_EMBEDDING_DAY024.md"
)

EXPECTED_GRAPH_DECISION = (
    "R2_SELECTED_FULL_DENSITY_FOUR_ATOM_BN_BRIDGE_GRAPH_VALIDATED"
)

EXPECTED_REPLAY_DECISION = (
    "R2_SELECTED_FOUR_ATOM_BN_BRIDGE_EXACT_CONFORMERS_RECOVERED"
)

PASS_DECISION = (
    "R2_SELECTED_FOUR_ATOM_BN_BRIDGE_HEAVY_COORDINATE_EMBEDDING_VALIDATED"
)

FAIL_DECISION = (
    "R2_SELECTED_FOUR_ATOM_BN_BRIDGE_HEAVY_COORDINATE_EMBEDDING_REQUIRES_REVIEW"
)

EXPECTED_HEAVY_NODES = 2112
EXPECTED_FIXED_HEAVY_NODES = 1992
EXPECTED_BRIDGE_NODES = 120
EXPECTED_BRIDGE_PATHS = 30
EXPECTED_HEAVY_EDGES = 3066

BN_TARGET_NM = 0.144973
MAX_BN_DEVIATION_NM = 0.003

MIN_ANGLE_DEG = 70.0
MAX_ANGLE_DEG = 175.0
MAX_ANGLE_RMS_DEVIATION_DEG = 30.0

MIN_NONBONDED_HEAVY_DISTANCE_NM = 0.120

TARGET_APERTURE_DIAMETER_NM = 0.839406
TARGET_OUTER_RADIUS_NM = 1.199126

MAX_APERTURE_RELATIVE_ERROR = 0.10
MAX_OUTER_RADIUS_RELATIVE_ERROR = 0.15
MAX_END_ASYMMETRY_NM = 0.010

CRITICAL_NODE_TYPES = {
    "HEXAGONAL_EDGE_COMPLETION_SEED",
    "ALTERNATING_BN_FOUR_ATOM_BRIDGE",
    "ANNULUS_OUTER_BOUNDARY",
    "ANNULUS_INNER_BOUNDARY",
}


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


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


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
            f"Non-finite value in field {key!r}"
        )

    return value


def normalized(
    vector: np.ndarray,
) -> np.ndarray:
    norm = float(np.linalg.norm(vector))

    if norm <= 1.0e-12:
        raise RuntimeError(
            "Could not normalize zero vector."
        )

    return vector / norm


def determine_axis(
    positions: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    center = np.mean(
        positions,
        axis=0,
    )

    centered = positions - center

    eigenvalues, eigenvectors = np.linalg.eigh(
        centered.T @ centered
    )

    axis = normalized(
        eigenvectors[
            :,
            int(np.argmax(eigenvalues)),
        ]
    )

    if axis[2] < 0.0:
        axis = -axis

    return center, axis


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


def radial_distance(
    point: np.ndarray,
    center: np.ndarray,
    axis: np.ndarray,
) -> float:
    displacement = point - center

    radial = (
        displacement
        - float(
            np.dot(
                displacement,
                axis,
            )
        )
        * axis
    )

    return float(
        np.linalg.norm(radial)
    )


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        SOURCE_FIXED_COORDINATES,
        GRAPH_NODES,
        GRAPH_EDGES,
        GRAPH_PATHS,
        GRAPH_SUMMARY,
        EXACT_BRIDGE_COORDINATES,
        REPLAY_SUMMARY,
    ):
        require_file(required)

    graph_summary = read_one(
        GRAPH_SUMMARY
    )

    replay_summary = read_one(
        REPLAY_SUMMARY
    )

    if graph_summary.get(
        "decision"
    ) != EXPECTED_GRAPH_DECISION:
        raise RuntimeError(
            "Gate 3M graph is not accepted."
        )

    if replay_summary.get(
        "decision"
    ) != EXPECTED_REPLAY_DECISION:
        raise RuntimeError(
            "Gate 3N exact replay is not accepted."
        )

    node_rows = read_rows(
        GRAPH_NODES
    )

    edge_rows = read_rows(
        GRAPH_EDGES
    )

    path_rows = read_rows(
        GRAPH_PATHS
    )

    source_coordinate_rows = read_rows(
        SOURCE_FIXED_COORDINATES
    )

    bridge_coordinate_rows = read_rows(
        EXACT_BRIDGE_COORDINATES
    )

    nodes = {
        row["node_id"]: row
        for row in node_rows
    }

    heavy_nodes = {
        node_id: row
        for node_id, row in nodes.items()
        if row["element"] != "H"
    }

    bridge_nodes = {
        node_id
        for node_id, row in heavy_nodes.items()
        if row["node_type"]
        == "ALTERNATING_BN_FOUR_ATOM_BRIDGE"
    }

    fixed_heavy_nodes = (
        set(heavy_nodes)
        - bridge_nodes
    )

    if len(heavy_nodes) != EXPECTED_HEAVY_NODES:
        raise RuntimeError(
            "Unexpected heavy-node population: "
            f"{len(heavy_nodes)}"
        )

    if len(fixed_heavy_nodes) != EXPECTED_FIXED_HEAVY_NODES:
        raise RuntimeError(
            "Unexpected fixed-heavy population: "
            f"{len(fixed_heavy_nodes)}"
        )

    if len(bridge_nodes) != EXPECTED_BRIDGE_NODES:
        raise RuntimeError(
            "Unexpected four-atom bridge population: "
            f"{len(bridge_nodes)}"
        )

    source_coordinates_all = {
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

    missing_fixed = sorted(
        fixed_heavy_nodes
        - set(source_coordinates_all)
    )

    if missing_fixed:
        raise RuntimeError(
            "Fixed heavy nodes missing from Gate 3K coordinates: "
            + " | ".join(missing_fixed[:20])
        )

    bridge_coordinates = {
        row["bridge_node"]: np.asarray(
            [
                parse_float(row, "x_nm"),
                parse_float(row, "y_nm"),
                parse_float(row, "z_nm"),
            ],
            dtype=float,
        )
        for row in bridge_coordinate_rows
    }

    if set(bridge_coordinates) != bridge_nodes:
        missing = sorted(
            bridge_nodes
            - set(bridge_coordinates)
        )

        extra = sorted(
            set(bridge_coordinates)
            - bridge_nodes
        )

        raise RuntimeError(
            "Bridge-coordinate identifiers disagree with Gate 3M. "
            f"Missing={missing[:10]}; extra={extra[:10]}"
        )

    coordinates: dict[str, np.ndarray] = {}

    for node_id in fixed_heavy_nodes:
        coordinates[node_id] = np.array(
            source_coordinates_all[node_id],
            dtype=float,
            copy=True,
        )

    for node_id in bridge_nodes:
        coordinates[node_id] = np.array(
            bridge_coordinates[node_id],
            dtype=float,
            copy=True,
        )

    if set(coordinates) != set(heavy_nodes):
        raise RuntimeError(
            "Heavy coordinate assignment is incomplete."
        )

    if not all(
        np.all(np.isfinite(point))
        for point in coordinates.values()
    ):
        raise RuntimeError(
            "Non-finite heavy coordinate detected."
        )

    heavy_adjacency = {
        node_id: set()
        for node_id in heavy_nodes
    }

    heavy_edge_rows = []

    for row in edge_rows:
        if not parse_bool(
            row["heavy_atom_edge"]
        ):
            continue

        first = row["source_node"]
        second = row["target_node"]

        if (
            first not in heavy_adjacency
            or second not in heavy_adjacency
        ):
            raise RuntimeError(
                "Heavy edge references a non-heavy node."
            )

        heavy_adjacency[first].add(second)
        heavy_adjacency[second].add(first)

        heavy_edge_rows.append(row)

    if len(heavy_edge_rows) != EXPECTED_HEAVY_EDGES:
        raise RuntimeError(
            "Unexpected heavy-edge count: "
            f"{len(heavy_edge_rows)}"
        )

    bond_rows = []

    for row in heavy_edge_rows:
        first = row["source_node"]
        second = row["target_node"]

        length = float(
            np.linalg.norm(
                coordinates[first]
                - coordinates[second]
            )
        )

        bond_rows.append(
            {
                "edge_id": row["edge_id"],
                "source_node": first,
                "target_node": second,
                "edge_type": row["edge_type"],
                "end": row["end"],
                "length_nm": length,
                "target_length_nm": BN_TARGET_NM,
                "deviation_nm": (
                    length
                    - BN_TARGET_NM
                ),
                "absolute_deviation_nm": abs(
                    length
                    - BN_TARGET_NM
                ),
                "contains_four_atom_bridge": (
                    first in bridge_nodes
                    or second in bridge_nodes
                ),
            }
        )

    write_rows(
        BOND_LENGTHS,
        bond_rows,
    )

    maximum_BN_deviation = max(
        float(row["absolute_deviation_nm"])
        for row in bond_rows
    )

    bridge_bond_rows = [
        row
        for row in bond_rows
        if row["contains_four_atom_bridge"]
    ]

    maximum_bridge_BN_deviation = max(
        float(row["absolute_deviation_nm"])
        for row in bridge_bond_rows
    )

    angle_values_by_type: dict[
        str,
        list[float]
    ] = {}

    critical_angles = []

    for center_id, neighbors in heavy_adjacency.items():
        neighbor_list = sorted(neighbors)

        if len(neighbor_list) < 2:
            continue

        node_type = heavy_nodes[
            center_id
        ]["node_type"]

        angle_values_by_type.setdefault(
            node_type,
            [],
        )

        for first_index in range(
            len(neighbor_list)
        ):
            for second_index in range(
                first_index + 1,
                len(neighbor_list),
            ):
                value = angle_degrees(
                    coordinates[
                        neighbor_list[
                            first_index
                        ]
                    ],
                    coordinates[
                        center_id
                    ],
                    coordinates[
                        neighbor_list[
                            second_index
                        ]
                    ],
                )

                angle_values_by_type[
                    node_type
                ].append(value)

                if node_type in CRITICAL_NODE_TYPES:
                    critical_angles.append(value)

    angle_summary_rows = []

    for node_type in sorted(
        angle_values_by_type
    ):
        values = np.asarray(
            angle_values_by_type[node_type],
            dtype=float,
        )

        if values.size == 0:
            continue

        angle_summary_rows.append(
            {
                "center_node_type": node_type,
                "angle_count": int(values.size),
                "minimum_angle_deg": float(
                    np.min(values)
                ),
                "mean_angle_deg": float(
                    np.mean(values)
                ),
                "maximum_angle_deg": float(
                    np.max(values)
                ),
                "RMS_deviation_from_120_deg": float(
                    np.sqrt(
                        np.mean(
                            (
                                values
                                - 120.0
                            )
                            ** 2
                        )
                    )
                ),
            }
        )

    write_rows(
        ANGLE_SUMMARY,
        angle_summary_rows,
    )

    critical_array = np.asarray(
        critical_angles,
        dtype=float,
    )

    critical_minimum = float(
        np.min(critical_array)
    )

    critical_mean = float(
        np.mean(critical_array)
    )

    critical_maximum = float(
        np.max(critical_array)
    )

    critical_rms = float(
        np.sqrt(
            np.mean(
                (
                    critical_array
                    - 120.0
                )
                ** 2
            )
        )
    )

    ordered_heavy_ids = sorted(
        heavy_nodes
    )

    positions = np.asarray(
        [
            coordinates[node_id]
            for node_id in ordered_heavy_ids
        ],
        dtype=float,
    )

    bonded_pairs = {
        tuple(
            sorted(
                (
                    row["source_node"],
                    row["target_node"],
                )
            )
        )
        for row in heavy_edge_rows
    }

    minimum_nonbonded_distance = math.inf
    minimum_nonbonded_pair = ""
    clash_count = 0

    for first_index, first_id in enumerate(
        ordered_heavy_ids
    ):
        if first_index + 1 >= len(
            ordered_heavy_ids
        ):
            break

        distances = np.linalg.norm(
            positions[first_index + 1:]
            - positions[first_index],
            axis=1,
        )

        for offset, distance in enumerate(
            distances,
            start=first_index + 1,
        ):
            second_id = ordered_heavy_ids[
                offset
            ]

            pair = tuple(
                sorted(
                    (
                        first_id,
                        second_id,
                    )
                )
            )

            if pair in bonded_pairs:
                continue

            value = float(distance)

            if value < minimum_nonbonded_distance:
                minimum_nonbonded_distance = value
                minimum_nonbonded_pair = (
                    f"{first_id} | {second_id}"
                )

            if (
                value
                < MIN_NONBONDED_HEAVY_DISTANCE_NM
            ):
                clash_count += 1

    write_rows(
        CONTACT_SUMMARY,
        [
            {
                "category": "HEAVY_HEAVY",
                "minimum_distance_nm": (
                    minimum_nonbonded_distance
                ),
                "minimum_pair": (
                    minimum_nonbonded_pair
                ),
                "threshold_nm": (
                    MIN_NONBONDED_HEAVY_DISTANCE_NM
                ),
                "pairs_below_threshold": (
                    clash_count
                ),
            }
        ],
    )

    parent_ids = [
        node_id
        for node_id, row in heavy_nodes.items()
        if row["node_type"] == "PARENT_HBN"
    ]

    parent_positions = np.asarray(
        [
            coordinates[node_id]
            for node_id in parent_ids
        ],
        dtype=float,
    )

    tube_center, tube_axis = determine_axis(
        parent_positions
    )

    end_rows = []

    for end in (
        "LOWER",
        "UPPER",
    ):
        inner_ids = [
            node_id
            for node_id, row in heavy_nodes.items()
            if (
                row["end"] == end
                and row["node_type"]
                == "ANNULUS_INNER_BOUNDARY"
            )
        ]

        outer_ids = [
            node_id
            for node_id, row in heavy_nodes.items()
            if (
                row["end"] == end
                and row["node_type"]
                == "ANNULUS_OUTER_BOUNDARY"
            )
        ]

        annulus_ids = [
            node_id
            for node_id, row in heavy_nodes.items()
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

        annulus_center = np.mean(
            np.asarray(
                [
                    coordinates[node_id]
                    for node_id in annulus_ids
                ],
                dtype=float,
            ),
            axis=0,
        )

        inner_radii = [
            radial_distance(
                coordinates[node_id],
                annulus_center,
                tube_axis,
            )
            for node_id in inner_ids
        ]

        outer_radii = [
            radial_distance(
                coordinates[node_id],
                annulus_center,
                tube_axis,
            )
            for node_id in outer_ids
        ]

        nuclear_aperture = (
            2.0
            * min(inner_radii)
        )

        mean_outer_radius = float(
            np.mean(outer_radii)
        )

        bridge_ids_end = [
            node_id
            for node_id in bridge_nodes
            if heavy_nodes[node_id]["end"] == end
        ]

        end_rows.append(
            {
                "end": end,
                "bridge_atoms": len(
                    bridge_ids_end
                ),
                "nuclear_aperture_diameter_nm": (
                    nuclear_aperture
                ),
                "nuclear_aperture_relative_error": abs(
                    nuclear_aperture
                    - TARGET_APERTURE_DIAMETER_NM
                ) / TARGET_APERTURE_DIAMETER_NM,
                "outer_annulus_radius_mean_nm": (
                    mean_outer_radius
                ),
                "outer_annulus_radius_relative_error": abs(
                    mean_outer_radius
                    - TARGET_OUTER_RADIUS_NM
                ) / TARGET_OUTER_RADIUS_NM,
            }
        )

    write_rows(
        END_SUMMARY,
        end_rows,
    )

    lower = next(
        row
        for row in end_rows
        if row["end"] == "LOWER"
    )

    upper = next(
        row
        for row in end_rows
        if row["end"] == "UPPER"
    )

    end_asymmetry = max(
        abs(
            float(
                lower[
                    "nuclear_aperture_diameter_nm"
                ]
            )
            - float(
                upper[
                    "nuclear_aperture_diameter_nm"
                ]
            )
        ),
        abs(
            float(
                lower[
                    "outer_annulus_radius_mean_nm"
                ]
            )
            - float(
                upper[
                    "outer_annulus_radius_mean_nm"
                ]
            )
        ),
    )

    fixed_coordinates_unchanged = all(
        np.array_equal(
            coordinates[node_id],
            source_coordinates_all[node_id],
        )
        for node_id in fixed_heavy_nodes
    )

    path_ids_from_graph = {
        row["bridge_path_id"]
        for row in path_rows
    }

    path_ids_from_coordinates = {
        row["bridge_path_id"]
        for row in bridge_coordinate_rows
    }

    coordinate_rows = []

    for node_id in ordered_heavy_ids:
        point = coordinates[node_id]

        coordinate_rows.append(
            {
                "node_id": node_id,
                "element": heavy_nodes[
                    node_id
                ]["element"],
                "node_type": heavy_nodes[
                    node_id
                ]["node_type"],
                "end": heavy_nodes[
                    node_id
                ]["end"],
                "x_nm": float(point[0]),
                "y_nm": float(point[1]),
                "z_nm": float(point[2]),
                "coordinate_source": (
                    "GATE3N_EXACT_SELECTED_CONFORMER"
                    if node_id in bridge_nodes
                    else "GATE3K_FIXED_HEAVY_COORDINATE"
                ),
                "energy_minimized": False,
                "MD_relaxed": False,
            }
        )

    write_rows(
        HEAVY_COORDINATES,
        coordinate_rows,
    )

    gates = {
        "Gate3M_graph_is_accepted": (
            graph_summary.get(
                "decision"
            )
            == EXPECTED_GRAPH_DECISION
        ),
        "Gate3N_exact_replay_is_accepted": (
            replay_summary.get(
                "decision"
            )
            == EXPECTED_REPLAY_DECISION
        ),
        "2112_heavy_nodes_received_coordinates": (
            len(coordinates)
            == EXPECTED_HEAVY_NODES
        ),
        "1992_fixed_heavy_nodes_were_preserved": (
            len(fixed_heavy_nodes)
            == EXPECTED_FIXED_HEAVY_NODES
            and fixed_coordinates_unchanged
        ),
        "120_exact_bridge_coordinates_were_applied": (
            len(bridge_coordinates)
            == EXPECTED_BRIDGE_NODES
        ),
        "30_graph_paths_match_coordinate_paths": (
            len(path_ids_from_graph)
            == EXPECTED_BRIDGE_PATHS
            and path_ids_from_graph
            == path_ids_from_coordinates
        ),
        "3066_heavy_edges_were_audited": (
            len(heavy_edge_rows)
            == EXPECTED_HEAVY_EDGES
        ),
        "all_BN_bonds_are_within_0p003nm": (
            maximum_BN_deviation
            <= MAX_BN_DEVIATION_NM
        ),
        "all_four_atom_bridge_BN_bonds_are_within_0p003nm": (
            maximum_bridge_BN_deviation
            <= MAX_BN_DEVIATION_NM
        ),
        "critical_angle_minimum_is_at_least70deg": (
            critical_minimum
            >= MIN_ANGLE_DEG
        ),
        "critical_angle_maximum_is_at_most175deg": (
            critical_maximum
            <= MAX_ANGLE_DEG
        ),
        "critical_angle_RMS_deviation_is_at_most30deg": (
            critical_rms
            <= MAX_ANGLE_RMS_DEVIATION_DEG
        ),
        "no_nonbonded_heavy_heavy_clashes": (
            clash_count == 0
        ),
        "aperture_errors_are_within10percent": all(
            float(
                row[
                    "nuclear_aperture_relative_error"
                ]
            )
            <= MAX_APERTURE_RELATIVE_ERROR
            for row in end_rows
        ),
        "outer_radius_errors_are_within15percent": all(
            float(
                row[
                    "outer_annulus_radius_relative_error"
                ]
            )
            <= MAX_OUTER_RADIUS_RELATIVE_ERROR
            for row in end_rows
        ),
        "lower_upper_asymmetry_is_within0p010nm": (
            end_asymmetry
            <= MAX_END_ASYMMETRY_NM
        ),
        "no_H_coordinates_were_generated": (
            len(coordinates)
            == EXPECTED_HEAVY_NODES
        ),
    }

    failed_gates = [
        name
        for name, passed in gates.items()
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

    required_next_step = (
        "GENERATE_AND_VALIDATE_R2_SELECTED_FOUR_ATOM_"
        "BN_BRIDGE_HYDROGEN_COORDINATES"
        if accepted
        else
        "REVIEW_R2_SELECTED_FOUR_ATOM_HEAVY_"
        "COORDINATE_EMBEDDING_FAILURES"
    )

    summary = {
        "decision": decision,
        "heavy_nodes": len(coordinates),
        "fixed_heavy_nodes": len(
            fixed_heavy_nodes
        ),
        "exact_bridge_nodes": len(
            bridge_coordinates
        ),
        "bridge_paths": len(
            path_ids_from_graph
        ),
        "heavy_edges": len(
            heavy_edge_rows
        ),
        "maximum_BN_bond_deviation_nm": (
            maximum_BN_deviation
        ),
        "maximum_bridge_BN_bond_deviation_nm": (
            maximum_bridge_BN_deviation
        ),
        "critical_angle_minimum_deg": (
            critical_minimum
        ),
        "critical_angle_mean_deg": (
            critical_mean
        ),
        "critical_angle_maximum_deg": (
            critical_maximum
        ),
        "critical_angle_RMS_deviation_deg": (
            critical_rms
        ),
        "minimum_nonbonded_heavy_heavy_nm": (
            minimum_nonbonded_distance
        ),
        "heavy_heavy_clash_count": (
            clash_count
        ),
        "lower_aperture_diameter_nm": (
            lower[
                "nuclear_aperture_diameter_nm"
            ]
        ),
        "upper_aperture_diameter_nm": (
            upper[
                "nuclear_aperture_diameter_nm"
            ]
        ),
        "maximum_lower_upper_asymmetry_nm": (
            end_asymmetry
        ),
        "fixed_coordinates_unchanged": (
            fixed_coordinates_unchanged
        ),
        "hydrogen_coordinates_generated": False,
        "candidate_is_final_chemistry": False,
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
                "end_summaries": end_rows,
                "gates": gates,
                "limitations": [
                    (
                        "This gate assembles and validates only "
                        "the 2112 heavy-atom coordinates."
                    ),
                    (
                        "The 120 bridge coordinates are the exact "
                        "deterministic Gate 3L selection recovered "
                        "by Gate 3N."
                    ),
                    (
                        "Hydrogen coordinates have not yet been generated."
                    ),
                    (
                        "No molecular topology, charges, force-field "
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

    manifest_rows = [
        {
            "role": "Gate3K_fixed_coordinates",
            "file": relative(
                SOURCE_FIXED_COORDINATES
            ),
            "sha256": sha256(
                SOURCE_FIXED_COORDINATES
            ),
        },
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
            "role": "Gate3M_graph_paths",
            "file": relative(
                GRAPH_PATHS
            ),
            "sha256": sha256(
                GRAPH_PATHS
            ),
        },
        {
            "role": "Gate3M_graph_summary",
            "file": relative(
                GRAPH_SUMMARY
            ),
            "sha256": sha256(
                GRAPH_SUMMARY
            ),
        },
        {
            "role": "Gate3N_exact_bridge_coordinates",
            "file": relative(
                EXACT_BRIDGE_COORDINATES
            ),
            "sha256": sha256(
                EXACT_BRIDGE_COORDINATES
            ),
        },
        {
            "role": "Gate3N_replay_summary",
            "file": relative(
                REPLAY_SUMMARY
            ),
            "sha256": sha256(
                REPLAY_SUMMARY
            ),
        },
    ]

    write_rows(
        MANIFEST,
        manifest_rows,
    )

    xyz_lines = [
        str(len(ordered_heavy_ids)),
        (
            "R2 four-atom BN bridge heavy embedding; "
            "no H; not energy minimized"
        ),
    ]

    for node_id in ordered_heavy_ids:
        point_angstrom = (
            coordinates[node_id]
            * 10.0
        )

        xyz_lines.append(
            f"{heavy_nodes[node_id]['element']:2s} "
            f"{point_angstrom[0]: .8f} "
            f"{point_angstrom[1]: .8f} "
            f"{point_angstrom[2]: .8f}"
        )

    XYZ_OUT.write_text(
        "\n".join(xyz_lines)
        + "\n",
        encoding="utf-8",
    )

    pdb_lines = [
        "REMARK R2 FOUR-ATOM BN BRIDGE HEAVY EMBEDDING",
        "REMARK NO HYDROGEN COORDINATES; NOT ENERGY MINIMIZED",
    ]

    for serial, node_id in enumerate(
        ordered_heavy_ids,
        start=1,
    ):
        point_angstrom = (
            coordinates[node_id]
            * 10.0
        )

        row = heavy_nodes[node_id]

        node_type = row["node_type"]

        if node_type == "PARENT_HBN":
            residue = "HBN"

        elif node_type == "HEXAGONAL_EDGE_COMPLETION_SEED":
            residue = "SED"

        elif node_type == "ALTERNATING_BN_FOUR_ATOM_BRIDGE":
            residue = "BR4"

        else:
            residue = "ANN"

        chain = (
            "L"
            if row["end"] == "LOWER"
            else (
                "U"
                if row["end"] == "UPPER"
                else "P"
            )
        )

        element = row["element"]
        atom_name = (
            element
            + str(serial % 1000)
        )[:4]

        pdb_lines.append(
            f"ATOM  {serial:5d} "
            f"{atom_name:>4s} "
            f"{residue:>3s} "
            f"{chain:1s}"
            f"{1:4d}    "
            f"{point_angstrom[0]:8.3f}"
            f"{point_angstrom[1]:8.3f}"
            f"{point_angstrom[2]:8.3f}"
            f"{1.00:6.2f}"
            f"{0.00:6.2f}          "
            f"{element:>2s}"
        )

    pdb_lines.append("END")

    PDB_OUT.write_text(
        "\n".join(pdb_lines)
        + "\n",
        encoding="utf-8",
    )

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed in gates.items()
    )

    REPORT.write_text(
        f"""# R2 Selected Four-Atom BN Bridge Heavy Coordinate Embedding

## Scope

This gate combines the fixed parent/seed/annulus coordinates with the
exact 120 bridge coordinates recovered from Gate 3L.

Hydrogen coordinates are not generated in this gate.

## Coordinate inventory

- Fixed heavy atoms: **{len(fixed_heavy_nodes)}**
- Four-atom bridge atoms: **{len(bridge_coordinates)}**
- Total heavy atoms: **{len(coordinates)}**
- Bridge paths: **{len(path_ids_from_graph)}**

## Bond geometry

- Maximum B-N deviation:
  **{maximum_BN_deviation:.9f} nm**
- Maximum bridge B-N deviation:
  **{maximum_bridge_BN_deviation:.9f} nm**

## Critical angles

- Minimum/mean/maximum:
  **{critical_minimum:.6f}/
  {critical_mean:.6f}/
  {critical_maximum:.6f} degrees**
- RMS deviation from 120 degrees:
  **{critical_rms:.6f} degrees**

## Heavy-atom clearance

- Minimum nonbonded heavy-heavy distance:
  **{minimum_nonbonded_distance:.9f} nm**
- Heavy-heavy clashes:
  **{clash_count}**

## Aperture and symmetry

- Lower/upper aperture:
  **{float(lower['nuclear_aperture_diameter_nm']):.9f}/
  {float(upper['nuclear_aperture_diameter_nm']):.9f} nm**
- Maximum lower-upper asymmetry:
  **{end_asymmetry:.9f} nm**

## Gates

{gate_lines}

## Decision

- Decision: **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Hydrogen coordinates generated: **NO**
- Candidate is final chemistry: **NO**
- Molecular topology generation authorized: **NO**
- Energy minimization authorized: **NO**
- MD authorized: **NO**
- QM authorized: **NO**
- Required next step:
  `{required_next_step}`
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 selected four-atom heavy "
        "coordinate embedding completed."
    )

    print(
        "Heavy coordinates fixed/bridge/total: "
        f"{len(fixed_heavy_nodes)}/"
        f"{len(bridge_coordinates)}/"
        f"{len(coordinates)}"
    )

    print(
        "Bridge paths / heavy edges: "
        f"{len(path_ids_from_graph)}/"
        f"{len(heavy_edge_rows)}"
    )

    print(
        "Fixed heavy coordinates unchanged: "
        f"{fixed_coordinates_unchanged}"
    )

    print(
        "Maximum BN / bridge-BN bond deviations: "
        f"{maximum_BN_deviation:.9f}/"
        f"{maximum_bridge_BN_deviation:.9f} nm"
    )

    print(
        "Critical angles min/mean/max/RMSdev120: "
        f"{critical_minimum:.6f}/"
        f"{critical_mean:.6f}/"
        f"{critical_maximum:.6f}/"
        f"{critical_rms:.6f} deg"
    )

    print(
        "Minimum nonbonded heavy-heavy / clash count: "
        f"{minimum_nonbonded_distance:.9f}/"
        f"{clash_count}"
    )

    for row in end_rows:
        print(
            f"{row['end']} bridge atoms / aperture/error / "
            "outer-radius/error: "
            f"{row['bridge_atoms']}/"
            f"{float(row['nuclear_aperture_diameter_nm']):.9f}/"
            f"{float(row['nuclear_aperture_relative_error']):.9f}/"
            f"{float(row['outer_annulus_radius_mean_nm']):.9f}/"
            f"{float(row['outer_annulus_radius_relative_error']):.9f}"
        )

    print(
        "Maximum lower-upper asymmetry: "
        f"{end_asymmetry:.9f} nm"
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
        "Hydrogen coordinates generated: NO"
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
        HEAVY_COORDINATES,
        BOND_LENGTHS,
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
            "Heavy coordinate embedding requires review."
        )


if __name__ == "__main__":
    main()
