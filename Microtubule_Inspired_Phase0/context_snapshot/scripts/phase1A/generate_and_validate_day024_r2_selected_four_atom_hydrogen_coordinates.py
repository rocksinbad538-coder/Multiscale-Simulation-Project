#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/phase1A/day024_chemical_end_rim_design"

GATE3K_SCRIPT = (
    ROOT
    / "scripts/phase1A/"
    "refine_day024_r2_trimer_bridge_conformers_and_h_orientations.py"
)

GATE3M = (
    BASE
    / "16_r2_selected_full_density_longer_bn_bridge_graph"
)

GATE3O = (
    BASE
    / "19_r2_selected_four_atom_heavy_coordinate_embedding"
)

GATE3O2 = (
    BASE
    / "21_r2_aperture_target_semantics_resolution"
)

GRAPH_NODES = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

GRAPH_EDGES = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_edges.csv"
)

GRAPH_SUMMARY = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_summary.csv"
)

HEAVY_COORDINATES = (
    GATE3O
    / "r2_selected_four_atom_heavy_coordinates.csv"
)

HEAVY_SUMMARY = (
    GATE3O
    / "r2_selected_four_atom_heavy_embedding_summary.csv"
)

SEMANTICS_SUMMARY = (
    GATE3O2
    / "r2_aperture_target_semantics_resolution_summary.csv"
)

OUT = (
    BASE
    / "23_r2_four_atom_hydrogen_coordinate_embedding"
)

FULL_COORDINATES = (
    OUT
    / "r2_selected_four_atom_full_coordinates.csv"
)

HYDROGEN_ORIENTATIONS = (
    OUT
    / "r2_selected_four_atom_hydrogen_orientations.csv"
)

BOND_SUMMARY = (
    OUT
    / "r2_selected_four_atom_hydrogen_bond_summary.csv"
)

CONTACT_SUMMARY = (
    OUT
    / "r2_selected_four_atom_hydrogen_contact_summary.csv"
)

END_SUMMARY = (
    OUT
    / "r2_selected_four_atom_hydrogen_end_summary.csv"
)

SUMMARY = (
    OUT
    / "r2_selected_four_atom_hydrogen_embedding_summary.csv"
)

GATES = (
    OUT
    / "r2_selected_four_atom_hydrogen_embedding_gates.csv"
)

JSON_OUT = (
    OUT
    / "r2_selected_four_atom_hydrogen_embedding.json"
)

MANIFEST = (
    OUT
    / "r2_selected_four_atom_hydrogen_embedding_manifest.csv"
)

XYZ_OUT = (
    OUT
    / "r2_selected_four_atom_full_embedding.xyz"
)

PDB_OUT = (
    OUT
    / "r2_selected_four_atom_full_embedding.pdb"
)

REPORT = (
    OUT
    / "R2_SELECTED_FOUR_ATOM_HYDROGEN_COORDINATE_EMBEDDING_DAY024.md"
)

EXPECTED_GRAPH_DECISION = (
    "R2_SELECTED_FULL_DENSITY_FOUR_ATOM_BN_BRIDGE_GRAPH_VALIDATED"
)

EXPECTED_SEMANTICS_DECISION = (
    "R2_HEAVY_COORDINATE_EMBEDDING_VALIDATED_"
    "APERTURE_FUNCTIONAL_GATE_DEFERRED"
)

PASS_DECISION = (
    "R2_SELECTED_FOUR_ATOM_BN_BRIDGE_HYDROGEN_COORDINATES_VALIDATED"
)

REVIEW_DECISION = (
    "R2_SELECTED_FOUR_ATOM_BN_BRIDGE_HYDROGEN_COORDINATES_REQUIRE_REVIEW"
)

N_HEAVY = 2112
N_H = 204
N_TOTAL = 2316

BH = 0.119
NH = 0.101

H_DIRECTION_COUNT = 362
H_POOL_SIZE = 48
H_GLOBAL_SWEEPS = 10

MIN_H_HEAVY_NM = 0.070
MIN_H_H_NM = 0.060

SOFT_H_HEAVY_NM = 0.090
SOFT_H_H_NM = 0.080

MIN_H_ANGLE_DEG = 70.0
MAX_H_ANGLE_DEG = 175.0
MAX_XH_DEVIATION_NM = 0.002
MAX_END_ASYMMETRY_NM = 0.010


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


def fibonacci_directions(
    count: int,
) -> np.ndarray:
    if count < 2:
        raise RuntimeError(
            "At least two directions are required."
        )

    golden_angle = (
        math.pi
        * (
            3.0
            - math.sqrt(5.0)
        )
    )

    directions = []

    for index in range(count):
        z_value = (
            1.0
            - 2.0
            * index
            / (
                count - 1
            )
        )

        radius = math.sqrt(
            max(
                0.0,
                1.0
                - z_value
                * z_value,
            )
        )

        theta = (
            golden_angle
            * index
        )

        directions.append(
            [
                radius
                * math.cos(theta),
                radius
                * math.sin(theta),
                z_value,
            ]
        )

    return np.asarray(
        directions,
        dtype=float,
    )


def soft_penalty(
    distances: np.ndarray,
    threshold: float,
) -> float:
    if distances.size == 0:
        return 0.0

    return float(
        np.sum(
            np.maximum(
                threshold
                - distances,
                0.0,
            )
            ** 2
        )
    )


def determine_tube_axis(
    parent_positions: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:
    center = np.mean(
        parent_positions,
        axis=0,
    )

    centered = (
        parent_positions
        - center
    )

    eigenvalues, eigenvectors = np.linalg.eigh(
        centered.T
        @ centered
    )

    axis = normalized(
        eigenvectors[
            :,
            int(
                np.argmax(
                    eigenvalues
                )
            ),
        ]
    )

    if axis[2] < 0.0:
        axis = -axis

    return center, axis


def preferred_direction(
    hydrogen_row: dict[str, str],
    heavy_position: np.ndarray,
    heavy_neighbor_positions: list[np.ndarray],
    tube_center: np.ndarray,
    tube_axis: np.ndarray,
    annulus_center_by_end: dict[str, np.ndarray],
) -> np.ndarray:
    vector_sum = np.zeros(
        3,
        dtype=float,
    )

    for neighbor_position in heavy_neighbor_positions:
        vector_sum += normalized(
            neighbor_position
            - heavy_position
        )

    anti_bond = (
        -vector_sum
    )

    if float(
        np.linalg.norm(
            anti_bond
        )
    ) <= 1.0e-12:
        radial = (
            heavy_position
            - tube_center
        )

        radial -= (
            np.dot(
                radial,
                tube_axis,
            )
            * tube_axis
        )

        anti_bond = radial

    anti_bond = normalized(
        anti_bond
    )

    end = hydrogen_row["end"]
    role = hydrogen_row["node_type"]

    if end in annulus_center_by_end:
        annulus_center = annulus_center_by_end[
            end
        ]

        radial = (
            heavy_position
            - annulus_center
        )

        radial -= (
            np.dot(
                radial,
                tube_axis,
            )
            * tube_axis
        )

        if float(
            np.linalg.norm(
                radial
            )
        ) > 1.0e-12:
            radial = normalized(
                radial
            )

            if role == "ANNULUS_INNER_PASSIVANT_H":
                preferred = -radial

            elif role in {
                "ANNULUS_OUTER_PASSIVANT_H",
                "SEED_PASSIVANT_H",
            }:
                preferred = radial

            else:
                preferred = anti_bond

            combined = (
                0.70
                * anti_bond
                + 0.30
                * preferred
            )

            if float(
                np.linalg.norm(
                    combined
                )
            ) > 1.0e-12:
                return normalized(
                    combined
                )

    return anti_bond


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        GATE3K_SCRIPT,
        GRAPH_NODES,
        GRAPH_EDGES,
        GRAPH_SUMMARY,
        HEAVY_COORDINATES,
        HEAVY_SUMMARY,
        SEMANTICS_SUMMARY,
    ):
        require_file(required)

    graph_summary = read_one(
        GRAPH_SUMMARY
    )

    heavy_summary = read_one(
        HEAVY_SUMMARY
    )

    semantics_summary = read_one(
        SEMANTICS_SUMMARY
    )

    if graph_summary.get(
        "decision"
    ) != EXPECTED_GRAPH_DECISION:
        raise RuntimeError(
            "Gate 3M graph is not accepted."
        )

    if semantics_summary.get(
        "decision"
    ) != EXPECTED_SEMANTICS_DECISION:
        raise RuntimeError(
            "Gate 3O.2 aperture semantics are not accepted."
        )

    if not parse_bool(
        semantics_summary.get(
            "hydrogen_coordinate_generation_authorized",
            "False",
        )
    ):
        raise RuntimeError(
            "Hydrogen-coordinate generation is not authorized."
        )

    node_rows = read_rows(
        GRAPH_NODES
    )

    edge_rows = read_rows(
        GRAPH_EDGES
    )

    heavy_coordinate_rows = read_rows(
        HEAVY_COORDINATES
    )

    nodes = {
        row["node_id"]: row
        for row in node_rows
    }

    heavy_ids = sorted(
        node_id
        for node_id, row
        in nodes.items()
        if row["element"] != "H"
    )

    hydrogen_ids = sorted(
        node_id
        for node_id, row
        in nodes.items()
        if row["element"] == "H"
    )

    if len(heavy_ids) != N_HEAVY:
        raise RuntimeError(
            f"Unexpected heavy population: {len(heavy_ids)}"
        )

    if len(hydrogen_ids) != N_H:
        raise RuntimeError(
            f"Unexpected H population: {len(hydrogen_ids)}"
        )

    positions: dict[str, np.ndarray] = {
        row["node_id"]: np.asarray(
            [
                parse_float(
                    row,
                    "x_nm",
                ),
                parse_float(
                    row,
                    "y_nm",
                ),
                parse_float(
                    row,
                    "z_nm",
                ),
            ],
            dtype=float,
        )
        for row in heavy_coordinate_rows
    }

    if set(positions) != set(
        heavy_ids
    ):
        raise RuntimeError(
            "Heavy coordinates do not match graph heavy nodes."
        )

    original_heavy_positions = {
        node_id: np.array(
            point,
            copy=True,
        )
        for node_id, point
        in positions.items()
    }

    adjacency = {
        node_id: set()
        for node_id in nodes
    }

    for row in edge_rows:
        first = row[
            "source_node"
        ]

        second = row[
            "target_node"
        ]

        adjacency[first].add(
            second
        )

        adjacency[second].add(
            first
        )

    attached_heavy_by_H = {}

    for hydrogen_id in hydrogen_ids:
        heavy_neighbors = [
            node_id
            for node_id
            in adjacency[
                hydrogen_id
            ]
            if nodes[
                node_id
            ][
                "element"
            ]
            != "H"
        ]

        if len(
            heavy_neighbors
        ) != 1:
            raise RuntimeError(
                f"{hydrogen_id}: expected one attached heavy atom."
            )

        attached_heavy_by_H[
            hydrogen_id
        ] = heavy_neighbors[0]

    parent_positions = np.asarray(
        [
            positions[node_id]
            for node_id in heavy_ids
            if nodes[
                node_id
            ][
                "node_type"
            ]
            == "PARENT_HBN"
        ],
        dtype=float,
    )

    tube_center, tube_axis = (
        determine_tube_axis(
            parent_positions
        )
    )

    annulus_center_by_end = {}

    for end in (
        "LOWER",
        "UPPER",
    ):
        annulus_ids = [
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

        annulus_center_by_end[
            end
        ] = np.mean(
            np.asarray(
                [
                    positions[
                        node_id
                    ]
                    for node_id
                    in annulus_ids
                ],
                dtype=float,
            ),
            axis=0,
        )

    directions = fibonacci_directions(
        H_DIRECTION_COUNT
    )

    heavy_position_array = np.asarray(
        [
            positions[
                node_id
            ]
            for node_id
            in heavy_ids
        ],
        dtype=float,
    )

    heavy_index_by_id = {
        node_id: index
        for index, node_id
        in enumerate(
            heavy_ids
        )
    }

    candidate_pools: dict[
        str,
        list[
            dict[str, Any]
        ]
    ] = {}

    print(
        "Generating hydrogen-orientation candidate pools..."
    )

    for hydrogen_index, hydrogen_id in enumerate(
        hydrogen_ids,
        start=1,
    ):
        hydrogen_row = nodes[
            hydrogen_id
        ]

        heavy_id = attached_heavy_by_H[
            hydrogen_id
        ]

        heavy_position = positions[
            heavy_id
        ]

        heavy_heavy_neighbors = [
            node_id
            for node_id
            in adjacency[
                heavy_id
            ]
            if nodes[
                node_id
            ][
                "element"
            ]
            != "H"
        ]

        if len(
            heavy_heavy_neighbors
        ) != 2:
            raise RuntimeError(
                f"{hydrogen_id}: attached heavy atom "
                f"{heavy_id} has "
                f"{len(heavy_heavy_neighbors)} heavy neighbors."
            )

        neighbor_positions = [
            positions[
                node_id
            ]
            for node_id
            in heavy_heavy_neighbors
        ]

        target_direction = preferred_direction(
            hydrogen_row,
            heavy_position,
            neighbor_positions,
            tube_center,
            tube_axis,
            annulus_center_by_end,
        )

        parent_element = nodes[
            heavy_id
        ][
            "element"
        ]

        bond_length = (
            BH
            if parent_element == "B"
            else NH
        )

        excluded_heavy_indices = {
            heavy_index_by_id[
                heavy_id
            ],
            *(
                heavy_index_by_id[
                    node_id
                ]
                for node_id
                in heavy_heavy_neighbors
            ),
        }

        fixed_mask = np.ones(
            len(
                heavy_ids
            ),
            dtype=bool,
        )

        for index in excluded_heavy_indices:
            fixed_mask[
                index
            ] = False

        local_fixed = heavy_position_array[
            fixed_mask
        ]

        pool = []

        for direction in directions:
            direction = normalized(
                direction
            )

            point = (
                heavy_position
                + bond_length
                * direction
            )

            angles = [
                angle_degrees(
                    neighbor_position,
                    heavy_position,
                    point,
                )
                for neighbor_position
                in neighbor_positions
            ]

            angle_violations = sum(
                (
                    value
                    < MIN_H_ANGLE_DEG
                    or value
                    > MAX_H_ANGLE_DEG
                )
                for value in angles
            )

            heavy_distances = np.linalg.norm(
                local_fixed
                - point,
                axis=1,
            )

            heavy_clashes = int(
                np.sum(
                    heavy_distances
                    < MIN_H_HEAVY_NM
                )
            )

            minimum_heavy_clearance = float(
                np.min(
                    heavy_distances
                )
            )

            preferred_alignment = float(
                np.dot(
                    direction,
                    target_direction,
                )
            )

            candidate = {
                "point": point,
                "direction": direction,
                "angles": angles,
                "angle_violations": (
                    angle_violations
                ),
                "heavy_clashes": (
                    heavy_clashes
                ),
                "minimum_heavy_clearance_nm": (
                    minimum_heavy_clearance
                ),
                "heavy_clearance_penalty": (
                    soft_penalty(
                        heavy_distances,
                        SOFT_H_HEAVY_NM,
                    )
                ),
                "preferred_alignment": (
                    preferred_alignment
                ),
            }

            pool.append(
                candidate
            )

        pool.sort(
            key=lambda candidate: (
                int(
                    candidate[
                        "heavy_clashes"
                    ]
                ),
                int(
                    candidate[
                        "angle_violations"
                    ]
                ),
                float(
                    candidate[
                        "heavy_clearance_penalty"
                    ]
                ),
                -float(
                    candidate[
                        "minimum_heavy_clearance_nm"
                    ]
                ),
                -float(
                    candidate[
                        "preferred_alignment"
                    ]
                ),
            )
        )

        candidate_pools[
            hydrogen_id
        ] = pool[
            :H_POOL_SIZE
        ]

        if (
            hydrogen_index
            % 25
            == 0
            or hydrogen_index
            == N_H
        ):
            best = candidate_pools[
                hydrogen_id
            ][0]

            print(
                f"  H pools {hydrogen_index}/{N_H}; "
                "latest best angle violations/"
                "heavy clashes="
                f"{best['angle_violations']}/"
                f"{best['heavy_clashes']}"
            )

    selected_index = {
        hydrogen_id: 0
        for hydrogen_id in hydrogen_ids
    }

    print(
        "Running global hydrogen-orientation coordinate descent..."
    )

    for sweep in range(
        H_GLOBAL_SWEEPS
    ):
        changed = 0

        ordered = (
            hydrogen_ids
            if sweep % 2 == 0
            else list(
                reversed(
                    hydrogen_ids
                )
            )
        )

        for hydrogen_id in ordered:
            other_points = np.asarray(
                [
                    candidate_pools[
                        other_id
                    ][
                        selected_index[
                            other_id
                        ]
                    ][
                        "point"
                    ]
                    for other_id
                    in hydrogen_ids
                    if other_id
                    != hydrogen_id
                ],
                dtype=float,
            )

            scores = []

            for candidate in candidate_pools[
                hydrogen_id
            ]:
                point = candidate[
                    "point"
                ]

                H_distances = np.linalg.norm(
                    other_points
                    - point,
                    axis=1,
                )

                H_clashes = int(
                    np.sum(
                        H_distances
                        < MIN_H_H_NM
                    )
                )

                minimum_H_clearance = float(
                    np.min(
                        H_distances
                    )
                )

                scores.append(
                    (
                        H_clashes,
                        int(
                            candidate[
                                "heavy_clashes"
                            ]
                        ),
                        int(
                            candidate[
                                "angle_violations"
                            ]
                        ),
                        soft_penalty(
                            H_distances,
                            SOFT_H_H_NM,
                        ),
                        float(
                            candidate[
                                "heavy_clearance_penalty"
                            ]
                        ),
                        -minimum_H_clearance,
                        -float(
                            candidate[
                                "minimum_heavy_clearance_nm"
                            ]
                        ),
                        -float(
                            candidate[
                                "preferred_alignment"
                            ]
                        ),
                    )
                )

            best_index = min(
                range(
                    len(
                        scores
                    )
                ),
                key=lambda index: scores[
                    index
                ],
            )

            if best_index != selected_index[
                hydrogen_id
            ]:
                selected_index[
                    hydrogen_id
                ] = best_index

                changed += 1

        print(
            f"  H sweep {sweep + 1}/{H_GLOBAL_SWEEPS}: "
            f"changed={changed}"
        )

        if changed == 0:
            break

    orientation_rows = []

    for hydrogen_id in hydrogen_ids:
        candidate = candidate_pools[
            hydrogen_id
        ][
            selected_index[
                hydrogen_id
            ]
        ]

        positions[
            hydrogen_id
        ] = np.asarray(
            candidate[
                "point"
            ],
            dtype=float,
        )

        heavy_id = attached_heavy_by_H[
            hydrogen_id
        ]

        direction = candidate[
            "direction"
        ]

        orientation_rows.append(
            {
                "hydrogen_node": hydrogen_id,
                "end": nodes[
                    hydrogen_id
                ][
                    "end"
                ],
                "hydrogen_role": nodes[
                    hydrogen_id
                ][
                    "node_type"
                ],
                "attached_heavy_node": (
                    heavy_id
                ),
                "attached_heavy_element": (
                    nodes[
                        heavy_id
                    ][
                        "element"
                    ]
                ),
                "direction_x": float(
                    direction[0]
                ),
                "direction_y": float(
                    direction[1]
                ),
                "direction_z": float(
                    direction[2]
                ),
                "local_angle_violations": (
                    candidate[
                        "angle_violations"
                    ]
                ),
                "local_heavy_clashes": (
                    candidate[
                        "heavy_clashes"
                    ]
                ),
                "minimum_heavy_clearance_nm": (
                    candidate[
                        "minimum_heavy_clearance_nm"
                    ]
                ),
                "heavy_neighbor_angle_1_deg": (
                    candidate[
                        "angles"
                    ][0]
                ),
                "heavy_neighbor_angle_2_deg": (
                    candidate[
                        "angles"
                    ][1]
                ),
                "preferred_alignment": (
                    candidate[
                        "preferred_alignment"
                    ]
                ),
                "selected_pool_index": (
                    selected_index[
                        hydrogen_id
                    ]
                ),
            }
        )

    write_rows(
        HYDROGEN_ORIENTATIONS,
        orientation_rows,
    )

    full_coordinate_rows = []

    for node_id in sorted(
        nodes
    ):
        point = positions[
            node_id
        ]

        full_coordinate_rows.append(
            {
                "node_id": node_id,
                "element": nodes[
                    node_id
                ][
                    "element"
                ],
                "node_type": nodes[
                    node_id
                ][
                    "node_type"
                ],
                "end": nodes[
                    node_id
                ][
                    "end"
                ],
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
                    "GATE3P_HYDROGEN_ORIENTATION_SEARCH"
                    if node_id in hydrogen_ids
                    else "GATE3O_VALIDATED_HEAVY_COORDINATE"
                ),
                "energy_minimized": False,
                "MD_relaxed": False,
            }
        )

    write_rows(
        FULL_COORDINATES,
        full_coordinate_rows,
    )

    heavy_coordinates_unchanged = all(
        np.array_equal(
            positions[
                node_id
            ],
            original_heavy_positions[
                node_id
            ],
        )
        for node_id in heavy_ids
    )

    bond_lengths_by_type = {
        "B-H": [],
        "N-H": [],
    }

    for hydrogen_id in hydrogen_ids:
        heavy_id = attached_heavy_by_H[
            hydrogen_id
        ]

        element = nodes[
            heavy_id
        ][
            "element"
        ]

        bond_type = (
            f"{element}-H"
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

        bond_lengths_by_type[
            bond_type
        ].append(
            distance
        )

    bond_summary_rows = []

    for bond_type, values in sorted(
        bond_lengths_by_type.items()
    ):
        array = np.asarray(
            values,
            dtype=float,
        )

        target = (
            BH
            if bond_type == "B-H"
            else NH
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
                "maximum_absolute_deviation_nm": float(
                    np.max(
                        np.abs(
                            array
                            - target
                        )
                    )
                ),
            }
        )

    write_rows(
        BOND_SUMMARY,
        bond_summary_rows,
    )

    H_positions = np.asarray(
        [
            positions[
                hydrogen_id
            ]
            for hydrogen_id
            in hydrogen_ids
        ],
        dtype=float,
    )

    H_heavy_minimum = math.inf
    H_heavy_clashes = 0
    H_heavy_pair = ""

    for hydrogen_index, hydrogen_id in enumerate(
        hydrogen_ids
    ):
        heavy_id = attached_heavy_by_H[
            hydrogen_id
        ]

        excluded = {
            heavy_id,
            *[
                node_id
                for node_id
                in adjacency[
                    heavy_id
                ]
                if nodes[
                    node_id
                ][
                    "element"
                ]
                != "H"
            ],
        }

        for heavy_node_id in heavy_ids:
            if heavy_node_id in excluded:
                continue

            distance = float(
                np.linalg.norm(
                    positions[
                        hydrogen_id
                    ]
                    - positions[
                        heavy_node_id
                    ]
                )
            )

            if distance < H_heavy_minimum:
                H_heavy_minimum = distance
                H_heavy_pair = (
                    f"{hydrogen_id} | "
                    f"{heavy_node_id}"
                )

            if distance < MIN_H_HEAVY_NM:
                H_heavy_clashes += 1

    H_H_minimum = math.inf
    H_H_clashes = 0
    H_H_pair = ""

    for first_index in range(
        len(
            hydrogen_ids
        )
    ):
        for second_index in range(
            first_index + 1,
            len(
                hydrogen_ids
            ),
        ):
            distance = float(
                np.linalg.norm(
                    H_positions[
                        first_index
                    ]
                    - H_positions[
                        second_index
                    ]
                )
            )

            if distance < H_H_minimum:
                H_H_minimum = distance
                H_H_pair = (
                    f"{hydrogen_ids[first_index]} | "
                    f"{hydrogen_ids[second_index]}"
                )

            if distance < MIN_H_H_NM:
                H_H_clashes += 1

    write_rows(
        CONTACT_SUMMARY,
        [
            {
                "category": "H_HEAVY",
                "minimum_distance_nm": (
                    H_heavy_minimum
                ),
                "minimum_pair": (
                    H_heavy_pair
                ),
                "threshold_nm": (
                    MIN_H_HEAVY_NM
                ),
                "clash_count": (
                    H_heavy_clashes
                ),
            },
            {
                "category": "H_H",
                "minimum_distance_nm": (
                    H_H_minimum
                ),
                "minimum_pair": (
                    H_H_pair
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

    end_rows = []

    for end in (
        "LOWER",
        "UPPER",
    ):
        inner_H_ids = [
            hydrogen_id
            for hydrogen_id
            in hydrogen_ids
            if (
                nodes[
                    hydrogen_id
                ][
                    "end"
                ]
                == end
                and nodes[
                    hydrogen_id
                ][
                    "node_type"
                ]
                == "ANNULUS_INNER_PASSIVANT_H"
            )
        ]

        annulus_center = (
            annulus_center_by_end[
                end
            ]
        )

        inner_radii = []

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

            inner_radii.append(
                float(
                    np.linalg.norm(
                        displacement
                    )
                )
            )

        inner_H_diameter = (
            2.0
            * min(
                inner_radii
            )
        )

        end_rows.append(
            {
                "end": end,
                "total_H": sum(
                    nodes[
                        hydrogen_id
                    ][
                        "end"
                    ]
                    == end
                    for hydrogen_id
                    in hydrogen_ids
                ),
                "inner_H_count": len(
                    inner_H_ids
                ),
                "inner_H_nuclear_aperture_diameter_nm": (
                    inner_H_diameter
                ),
                "inner_H_radius_minimum_nm": min(
                    inner_radii
                ),
                "inner_H_radius_mean_nm": float(
                    np.mean(
                        inner_radii
                    )
                ),
                "inner_H_radius_maximum_nm": max(
                    inner_radii
                ),
            }
        )

    write_rows(
        END_SUMMARY,
        end_rows,
    )

    lower = next(
        row
        for row in end_rows
        if row["end"]
        == "LOWER"
    )

    upper = next(
        row
        for row in end_rows
        if row["end"]
        == "UPPER"
    )

    aperture_asymmetry = abs(
        float(
            lower[
                "inner_H_nuclear_aperture_diameter_nm"
            ]
        )
        - float(
            upper[
                "inner_H_nuclear_aperture_diameter_nm"
            ]
        )
    )

    H_role_counts = Counter(
        nodes[
            hydrogen_id
        ][
            "node_type"
        ]
        for hydrogen_id
        in hydrogen_ids
    )

    H_end_counts = Counter(
        nodes[
            hydrogen_id
        ][
            "end"
        ]
        for hydrogen_id
        in hydrogen_ids
    )

    H_parent_element_counts = Counter(
        nodes[
            attached_heavy_by_H[
                hydrogen_id
            ]
        ][
            "element"
        ]
        for hydrogen_id
        in hydrogen_ids
    )

    maximum_XH_deviation = max(
        float(
            row[
                "maximum_absolute_deviation_nm"
            ]
        )
        for row in bond_summary_rows
    )

    local_angle_violations = sum(
        int(
            row[
                "local_angle_violations"
            ]
        )
        for row in orientation_rows
    )

    local_heavy_clashes = sum(
        int(
            row[
                "local_heavy_clashes"
            ]
        )
        for row in orientation_rows
    )

    gates = {
        "Gate3M_graph_is_accepted": (
            graph_summary.get(
                "decision"
            )
            == EXPECTED_GRAPH_DECISION
        ),
        "Gate3O2_semantics_resolution_is_accepted": (
            semantics_summary.get(
                "decision"
            )
            == EXPECTED_SEMANTICS_DECISION
        ),
        "2316_nodes_received_coordinates": (
            len(
                positions
            )
            == N_TOTAL
        ),
        "2112_heavy_coordinates_are_unchanged": (
            heavy_coordinates_unchanged
        ),
        "204_H_coordinates_were_generated": (
            len(
                hydrogen_ids
            )
            == N_H
        ),
        "102_H_were_generated_per_end": (
            H_end_counts[
                "LOWER"
            ]
            == 102
            and H_end_counts[
                "UPPER"
            ]
            == 102
        ),
        "120_bridge_H_were_generated": (
            H_role_counts[
                "BRIDGE_PASSIVANT_H"
            ]
            == 120
        ),
        "30_seed_H_were_generated": (
            H_role_counts[
                "SEED_PASSIVANT_H"
            ]
            == 30
        ),
        "30_outer_H_were_generated": (
            H_role_counts[
                "ANNULUS_OUTER_PASSIVANT_H"
            ]
            == 30
        ),
        "24_inner_H_were_generated": (
            H_role_counts[
                "ANNULUS_INNER_PASSIVANT_H"
            ]
            == 24
        ),
        "102_BH_and_102_NH_bonds_are_present": (
            H_parent_element_counts[
                "B"
            ]
            == 102
            and H_parent_element_counts[
                "N"
            ]
            == 102
        ),
        "all_XH_bonds_are_within_0p002nm": (
            maximum_XH_deviation
            <= MAX_XH_DEVIATION_NM
        ),
        "all_local_H_angles_are_within70to175deg": (
            local_angle_violations
            == 0
        ),
        "all_local_H_heavy_clearance_checks_pass": (
            local_heavy_clashes
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
        "inner_H_aperture_is_lower_upper_symmetric": (
            aperture_asymmetry
            <= MAX_END_ASYMMETRY_NM
        ),
        "inner_H_aperture_is_recorded_as_geometric_proxy_only": True,
        "no_energy_minimization_or_MD_was_performed": True,
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
        "REFINE_R2_SELECTED_FOUR_ATOM_HYDROGEN_ORIENTATIONS"
    )

    summary = {
        "decision": decision,
        "coordinate_nodes": len(
            positions
        ),
        "heavy_coordinates_preserved": (
            heavy_coordinates_unchanged
        ),
        "H_coordinates_generated": len(
            hydrogen_ids
        ),
        "H_lower": H_end_counts[
            "LOWER"
        ],
        "H_upper": H_end_counts[
            "UPPER"
        ],
        "B_H_bonds": H_parent_element_counts[
            "B"
        ],
        "N_H_bonds": H_parent_element_counts[
            "N"
        ],
        "maximum_XH_deviation_nm": (
            maximum_XH_deviation
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
        "local_H_angle_violations": (
            local_angle_violations
        ),
        "local_H_heavy_clashes": (
            local_heavy_clashes
        ),
        "lower_inner_H_nuclear_aperture_nm": (
            lower[
                "inner_H_nuclear_aperture_diameter_nm"
            ]
        ),
        "upper_inner_H_nuclear_aperture_nm": (
            upper[
                "inner_H_nuclear_aperture_diameter_nm"
            ]
        ),
        "inner_H_aperture_asymmetry_nm": (
            aperture_asymmetry
        ),
        "inner_H_aperture_is_functional_5kBT_result": False,
        "candidate_is_final_chemistry": False,
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
                "end_summaries": end_rows,
                "gates": gates,
                "parameters": {
                    "B_H_target_nm": BH,
                    "N_H_target_nm": NH,
                    "H_direction_count": H_DIRECTION_COUNT,
                    "H_pool_size": H_POOL_SIZE,
                    "H_global_sweeps": H_GLOBAL_SWEEPS,
                    "minimum_H_heavy_nm": MIN_H_HEAVY_NM,
                    "minimum_H_H_nm": MIN_H_H_NM,
                },
                "limitations": [
                    (
                        "Hydrogen placement is a deterministic static "
                        "geometric search, not an energy optimization."
                    ),
                    (
                        "The inner-H nuclear aperture is a geometric "
                        "proxy only and is not the effective 5 kBT "
                        "water-accessible aperture."
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

    write_rows(
        MANIFEST,
        [
            {
                "role": "Gate3K_H_orientation_reference",
                "file": relative(
                    GATE3K_SCRIPT
                ),
                "sha256": sha256(
                    GATE3K_SCRIPT
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
                "role": "Gate3O_heavy_coordinates",
                "file": relative(
                    HEAVY_COORDINATES
                ),
                "sha256": sha256(
                    HEAVY_COORDINATES
                ),
            },
            {
                "role": "Gate3O2_semantics_summary",
                "file": relative(
                    SEMANTICS_SUMMARY
                ),
                "sha256": sha256(
                    SEMANTICS_SUMMARY
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
            "R2 four-atom BN bridge full static embedding; "
            "not energy minimized"
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
        "REMARK R2 FOUR-ATOM BN BRIDGE FULL STATIC EMBEDDING",
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
        f"""# R2 Four-Atom BN Bridge Hydrogen Coordinate Embedding

## Inventory

- Total coordinate nodes: **{len(positions)}**
- Heavy atoms: **{len(heavy_ids)}**
- Hydrogen atoms: **{len(hydrogen_ids)}**
- Lower/upper H: **{H_end_counts['LOWER']}/{H_end_counts['UPPER']}**
- B-H/N-H: **{H_parent_element_counts['B']}/{H_parent_element_counts['N']}**

## X-H geometry

- Maximum X-H deviation:
  **{maximum_XH_deviation:.9f} nm**
- Local H-angle violations:
  **{local_angle_violations}**

## Nonbonded clearance

- Minimum H-heavy / clashes:
  **{H_heavy_minimum:.9f}/{H_heavy_clashes}**
- Minimum H-H / clashes:
  **{H_H_minimum:.9f}/{H_H_clashes}**

## Inner-H aperture proxy

- Lower:
  **{float(lower['inner_H_nuclear_aperture_diameter_nm']):.9f} nm**
- Upper:
  **{float(upper['inner_H_nuclear_aperture_diameter_nm']):.9f} nm**
- Asymmetry:
  **{aperture_asymmetry:.9f} nm**

This is a nuclear geometric proxy, not the effective 5 kBT aperture.

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Candidate is final chemistry:
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
        "Day024 R2 four-atom hydrogen coordinate "
        "embedding completed."
    )

    print(
        "Coordinate nodes heavy/H/total: "
        f"{len(heavy_ids)}/"
        f"{len(hydrogen_ids)}/"
        f"{len(positions)}"
    )

    print(
        "H lower/upper and B-H/N-H: "
        f"{H_end_counts['LOWER']}/"
        f"{H_end_counts['UPPER']} and "
        f"{H_parent_element_counts['B']}/"
        f"{H_parent_element_counts['N']}"
    )

    print(
        "Heavy coordinates unchanged: "
        f"{heavy_coordinates_unchanged}"
    )

    print(
        "Maximum X-H deviation: "
        f"{maximum_XH_deviation:.9f} nm"
    )

    print(
        "Local H angle violations / local heavy clashes: "
        f"{local_angle_violations}/"
        f"{local_heavy_clashes}"
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
        f"{float(lower['inner_H_nuclear_aperture_diameter_nm']):.9f}/"
        f"{float(upper['inner_H_nuclear_aperture_diameter_nm']):.9f}/"
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
        FULL_COORDINATES,
        HYDROGEN_ORIENTATIONS,
        BOND_SUMMARY,
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
            "Hydrogen coordinate embedding requires review."
        )


if __name__ == "__main__":
    main()
