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

G3I = BASE / "11_r2_alternating_bn_trimer_bridge_graph"
G3K = BASE / "13_r2_trimer_bridge_conformer_and_h_refinement"

OUT = (
    BASE
    / "14_r2_trimer_bridge_search_completeness_audit"
)

GRAPH_NODES = (
    G3I
    / "r2_alternating_bn_trimer_bridge_graph_nodes.csv"
)

GRAPH_EDGES = (
    G3I
    / "r2_alternating_bn_trimer_bridge_graph_edges.csv"
)

BRIDGE_PATHS = (
    G3I
    / "r2_alternating_bn_trimer_bridge_paths.csv"
)

GRAPH_SUMMARY = (
    G3I
    / "r2_alternating_bn_trimer_bridge_graph_summary.csv"
)

SOURCE_COORDINATES = (
    G3K
    / "r2_trimer_bridge_refined_coordinates.csv"
)

SOURCE_CONFORMERS = (
    G3K
    / "r2_trimer_bridge_refined_conformers.csv"
)

SOURCE_SUMMARY = (
    G3K
    / "r2_trimer_bridge_refinement_summary.csv"
)

LIBRARY_SUMMARY = (
    OUT
    / "r2_trimer_independent_angle_library_summary.csv"
)

PATH_RESULTS = (
    OUT
    / "r2_trimer_expanded_search_path_results.csv"
)

BEST_CANDIDATES = (
    OUT
    / "r2_trimer_expanded_search_best_candidate_coordinates.csv"
)

SUMMARY = (
    OUT
    / "r2_trimer_search_completeness_summary.csv"
)

GATES = (
    OUT
    / "r2_trimer_search_completeness_gates.csv"
)

JSON_OUT = (
    OUT
    / "r2_trimer_search_completeness.json"
)

MANIFEST = (
    OUT
    / "r2_trimer_search_completeness_source_manifest.csv"
)

REPORT = (
    OUT
    / "R2_TRIMER_BRIDGE_SEARCH_COMPLETENESS_AUDIT_DAY024.md"
)

EXPECTED_GRAPH_DECISION = (
    "R2_ALTERNATING_BN_TRIMER_BRIDGE_GRAPH_VALIDATED"
)

EXPECTED_SOURCE_DECISION = (
    "R2_ALTERNATING_BN_TRIMER_BRIDGE_REFINEMENT_"
    "REQUIRES_BRIDGE_TOPOLOGY_REDESIGN"
)

RETAIN_DECISION = (
    "R2_TRIMER_BRIDGE_TOPOLOGY_RETAINED_AFTER_"
    "EXPANDED_CONFORMER_SEARCH"
)

REDESIGN_DECISION = (
    "R2_TRIMER_BRIDGE_TOPOLOGY_REDESIGN_CONFIRMED_"
    "BY_EXPANDED_CONFORMER_SEARCH"
)

BN_TARGET_NM = 0.144973

ANGLE_VALUES_DEG = tuple(
    float(value)
    for value in range(
        105,
        136,
        5,
    )
)

TORSION_VALUES_DEG = tuple(
    float(value)
    for value in range(
        0,
        360,
        30,
    )
)

AZIMUTH_VALUES_DEG = tuple(
    float(value)
    for value in range(
        0,
        360,
        5,
    )
)

BASE_CONFORMERS_PER_PATH = 600
LOCAL_FIXED_RADIUS_NM = 0.80

MIN_ACCEPTED_ANGLE_DEG = 70.0
MAX_ACCEPTED_ANGLE_DEG = 175.0

MAX_BOND_DEVIATION_NM = 0.003
MIN_NONBONDED_HEAVY_DISTANCE_NM = 0.120

EXPECTED_TOTAL_PATHS = 30
EXPECTED_BRIDGE_ATOMS_PER_PATH = 3


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
    )


def require_file(path: Path) -> None:
    if (
        not path.is_file()
        or path.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing or empty required file: {path}"
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


def read_rows(
    path: Path,
) -> list[dict[str, str]]:
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


def read_one(
    path: Path,
) -> dict[str, str]:
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
            f"Non-finite numeric field {key!r}"
        )

    return value


def parse_int(
    row: dict[str, str],
    key: str,
) -> int:
    try:
        return int(float(row[key]))
    except (
        KeyError,
        TypeError,
        ValueError,
    ) as exc:
        raise RuntimeError(
            f"Could not parse integer field {key!r}"
        ) from exc


def normalized(
    vector: np.ndarray,
    fallback: np.ndarray | None = None,
) -> np.ndarray:
    norm = float(
        np.linalg.norm(vector)
    )

    if norm > 1.0e-12:
        return vector / norm

    if fallback is None:
        raise RuntimeError(
            "Could not normalize a zero vector."
        )

    fallback_norm = float(
        np.linalg.norm(fallback)
    )

    if fallback_norm <= 1.0e-12:
        raise RuntimeError(
            "Fallback vector is also zero."
        )

    return fallback / fallback_norm


def rotation_about_axis(
    axis: np.ndarray,
    angle_rad: float,
) -> np.ndarray:
    x_value, y_value, z_value = normalized(axis)

    cosine = math.cos(angle_rad)
    sine = math.sin(angle_rad)
    factor = 1.0 - cosine

    return np.asarray(
        [
            [
                cosine
                + x_value * x_value * factor,
                x_value * y_value * factor
                - z_value * sine,
                x_value * z_value * factor
                + y_value * sine,
            ],
            [
                y_value * x_value * factor
                + z_value * sine,
                cosine
                + y_value * y_value * factor,
                y_value * z_value * factor
                - x_value * sine,
            ],
            [
                z_value * x_value * factor
                - y_value * sine,
                z_value * y_value * factor
                + x_value * sine,
                cosine
                + z_value * z_value * factor,
            ],
        ],
        dtype=float,
    )


def rotation_from_vectors(
    source: np.ndarray,
    target: np.ndarray,
) -> np.ndarray:
    source_unit = normalized(source)
    target_unit = normalized(target)

    cross = np.cross(
        source_unit,
        target_unit,
    )

    sine = float(
        np.linalg.norm(cross)
    )

    cosine = float(
        np.clip(
            np.dot(
                source_unit,
                target_unit,
            ),
            -1.0,
            1.0,
        )
    )

    if sine <= 1.0e-12:
        if cosine > 0.0:
            return np.eye(
                3,
                dtype=float,
            )

        reference = np.asarray(
            [
                1.0,
                0.0,
                0.0,
            ],
            dtype=float,
        )

        if abs(
            float(
                np.dot(
                    source_unit,
                    reference,
                )
            )
        ) > 0.90:
            reference = np.asarray(
                [
                    0.0,
                    1.0,
                    0.0,
                ],
                dtype=float,
            )

        return rotation_about_axis(
            normalized(
                np.cross(
                    source_unit,
                    reference,
                )
            ),
            math.pi,
        )

    return rotation_about_axis(
        cross / sine,
        math.atan2(
            sine,
            cosine,
        ),
    )


def place_next(
    first: np.ndarray,
    second: np.ndarray,
    third: np.ndarray,
    angle_rad: float,
    dihedral_rad: float,
) -> np.ndarray:
    previous = normalized(
        third - second
    )

    normal = np.cross(
        second - first,
        previous,
    )

    if float(
        np.linalg.norm(normal)
    ) <= 1.0e-12:
        reference = np.asarray(
            [
                0.0,
                0.0,
                1.0,
            ],
            dtype=float,
        )

        if abs(
            float(
                np.dot(
                    previous,
                    reference,
                )
            )
        ) > 0.90:
            reference = np.asarray(
                [
                    0.0,
                    1.0,
                    0.0,
                ],
                dtype=float,
            )

        normal = np.cross(
            previous,
            reference,
        )

    normal = normalized(normal)

    in_plane = normalized(
        np.cross(
            normal,
            previous,
        )
    )

    direction = (
        -math.cos(angle_rad)
        * previous
        + math.sin(angle_rad)
        * (
            math.cos(dihedral_rad)
            * in_plane
            + math.sin(dihedral_rad)
            * normal
        )
    )

    return (
        third
        + BN_TARGET_NM
        * direction
    )


def canonical_chain(
    angle_1_deg: float,
    angle_2_deg: float,
    angle_3_deg: float,
    torsion_1_deg: float,
    torsion_2_deg: float,
) -> np.ndarray:
    angle_1 = math.radians(
        angle_1_deg
    )

    point_0 = np.zeros(
        3,
        dtype=float,
    )

    point_1 = np.asarray(
        [
            BN_TARGET_NM,
            0.0,
            0.0,
        ],
        dtype=float,
    )

    point_2 = (
        point_1
        + BN_TARGET_NM
        * np.asarray(
            [
                -math.cos(angle_1),
                math.sin(angle_1),
                0.0,
            ],
            dtype=float,
        )
    )

    point_3 = place_next(
        point_0,
        point_1,
        point_2,
        math.radians(
            angle_2_deg
        ),
        math.radians(
            torsion_1_deg
        ),
    )

    point_4 = place_next(
        point_1,
        point_2,
        point_3,
        math.radians(
            angle_3_deg
        ),
        math.radians(
            torsion_2_deg
        ),
    )

    return np.asarray(
        [
            point_0,
            point_1,
            point_2,
            point_3,
            point_4,
        ],
        dtype=float,
    )


def build_library() -> tuple[
    np.ndarray,
    list[dict[str, float]],
]:
    chains = []
    metadata = []

    for angle_1 in ANGLE_VALUES_DEG:
        for angle_2 in ANGLE_VALUES_DEG:
            for angle_3 in ANGLE_VALUES_DEG:
                for torsion_1 in TORSION_VALUES_DEG:
                    for torsion_2 in TORSION_VALUES_DEG:
                        chain = canonical_chain(
                            angle_1,
                            angle_2,
                            angle_3,
                            torsion_1,
                            torsion_2,
                        )

                        chains.append(
                            chain
                        )

                        metadata.append(
                            {
                                "angle_1_deg": angle_1,
                                "angle_2_deg": angle_2,
                                "angle_3_deg": angle_3,
                                "torsion_1_deg": torsion_1,
                                "torsion_2_deg": torsion_2,
                                "endpoint_distance_nm": float(
                                    np.linalg.norm(
                                        chain[-1]
                                        - chain[0]
                                    )
                                ),
                            }
                        )

    return (
        np.asarray(
            chains,
            dtype=float,
        ),
        metadata,
    )


def map_chain(
    chain: np.ndarray,
    start: np.ndarray,
    finish: np.ndarray,
    azimuth_deg: float,
    mirror: bool,
) -> np.ndarray:
    working = np.array(
        chain,
        dtype=float,
        copy=True,
    )

    if mirror:
        working[:, 2] *= -1.0

    endpoint_vector = (
        finish - start
    )

    rotation = rotation_from_vectors(
        working[-1]
        - working[0],
        endpoint_vector,
    )

    mapped = (
        working
        - working[0]
    ) @ rotation.T

    mapped = (
        mapped
        @ rotation_about_axis(
            normalized(
                endpoint_vector
            ),
            math.radians(
                azimuth_deg
            ),
        ).T
    )

    correction = (
        endpoint_vector
        - mapped[-1]
    )

    for index in range(
        1,
        5,
    ):
        mapped[index] += (
            index / 4.0
        ) * correction

    mapped += start

    return mapped


def angle_degrees(
    first: np.ndarray,
    center: np.ndarray,
    second: np.ndarray,
) -> float:
    cosine = float(
        np.clip(
            np.dot(
                normalized(
                    first - center
                ),
                normalized(
                    second - center
                ),
            ),
            -1.0,
            1.0,
        )
    )

    return math.degrees(
        math.acos(cosine)
    )


def evaluate_candidate(
    chain: np.ndarray,
    start_neighbors: list[np.ndarray],
    finish_neighbors: list[np.ndarray],
    local_fixed: np.ndarray,
) -> dict[str, Any]:
    bond_lengths = np.linalg.norm(
        np.diff(
            chain,
            axis=0,
        ),
        axis=1,
    )

    bond_deviations = np.abs(
        bond_lengths
        - BN_TARGET_NM
    )

    internal_angles = [
        angle_degrees(
            chain[index - 1],
            chain[index],
            chain[index + 1],
        )
        for index in (
            1,
            2,
            3,
        )
    ]

    seed_angles = [
        angle_degrees(
            neighbor,
            chain[0],
            chain[1],
        )
        for neighbor
        in start_neighbors
    ]

    annulus_angles = [
        angle_degrees(
            neighbor,
            chain[4],
            chain[3],
        )
        for neighbor
        in finish_neighbors
    ]

    all_angles = np.asarray(
        [
            *internal_angles,
            *seed_angles,
            *annulus_angles,
        ],
        dtype=float,
    )

    angle_violations = int(
        np.sum(
            (
                all_angles
                < MIN_ACCEPTED_ANGLE_DEG
            )
            | (
                all_angles
                > MAX_ACCEPTED_ANGLE_DEG
            )
        )
    )

    angle_deficit = float(
        np.sum(
            np.maximum(
                MIN_ACCEPTED_ANGLE_DEG
                - all_angles,
                0.0,
            )
            + np.maximum(
                all_angles
                - MAX_ACCEPTED_ANGLE_DEG,
                0.0,
            )
        )
    )

    own_pairs = (
        (0, 2),
        (0, 3),
        (0, 4),
        (1, 3),
        (1, 4),
        (2, 4),
    )

    own_distances = np.asarray(
        [
            float(
                np.linalg.norm(
                    chain[first]
                    - chain[second]
                )
            )
            for first, second
            in own_pairs
        ],
        dtype=float,
    )

    internal = chain[
        1:4
    ]

    if local_fixed.size:
        fixed_distances = np.linalg.norm(
            internal[:, None, :]
            - local_fixed[
                None,
                :,
                :,
            ],
            axis=2,
        ).reshape(-1)
    else:
        fixed_distances = np.asarray(
            [],
            dtype=float,
        )

    all_nonbonded = np.concatenate(
        (
            own_distances,
            fixed_distances,
        )
    )

    clash_count = int(
        np.sum(
            all_nonbonded
            < MIN_NONBONDED_HEAVY_DISTANCE_NM
        )
    )

    minimum_clearance = float(
        np.min(
            all_nonbonded
        )
    )

    maximum_bond_deviation = float(
        np.max(
            bond_deviations
        )
    )

    bond_failure_count = int(
        np.sum(
            bond_deviations
            > MAX_BOND_DEVIATION_NM
        )
    )

    passes = (
        angle_violations == 0
        and clash_count == 0
        and bond_failure_count == 0
    )

    return {
        "bond_lengths": bond_lengths,
        "maximum_bond_deviation_nm": (
            maximum_bond_deviation
        ),
        "bond_failure_count": (
            bond_failure_count
        ),
        "internal_angles": (
            internal_angles
        ),
        "seed_angles": (
            seed_angles
        ),
        "annulus_angles": (
            annulus_angles
        ),
        "minimum_angle_deg": float(
            np.min(
                all_angles
            )
        ),
        "maximum_angle_deg": float(
            np.max(
                all_angles
            )
        ),
        "angle_violations": (
            angle_violations
        ),
        "angle_deficit_deg": (
            angle_deficit
        ),
        "clash_count": (
            clash_count
        ),
        "minimum_clearance_nm": (
            minimum_clearance
        ),
        "passes": passes,
    }


def candidate_score(
    candidate: dict[str, Any],
) -> tuple[Any, ...]:
    return (
        0
        if bool(
            candidate[
                "passes"
            ]
        )
        else 1,
        int(
            candidate[
                "angle_violations"
            ]
        ),
        int(
            candidate[
                "clash_count"
            ]
        ),
        int(
            candidate[
                "bond_failure_count"
            ]
        ),
        float(
            candidate[
                "angle_deficit_deg"
            ]
        ),
        float(
            candidate[
                "maximum_bond_deviation_nm"
            ]
        ),
        -float(
            candidate[
                "minimum_clearance_nm"
            ]
        ),
        float(
            candidate[
                "library_distance_error_nm"
            ]
        ),
    )


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        GRAPH_NODES,
        GRAPH_EDGES,
        BRIDGE_PATHS,
        GRAPH_SUMMARY,
        SOURCE_COORDINATES,
        SOURCE_CONFORMERS,
        SOURCE_SUMMARY,
    ):
        require_file(required)

    graph_summary = read_one(
        GRAPH_SUMMARY
    )

    source_summary = read_one(
        SOURCE_SUMMARY
    )

    node_rows = read_rows(
        GRAPH_NODES
    )

    edge_rows = read_rows(
        GRAPH_EDGES
    )

    path_rows = read_rows(
        BRIDGE_PATHS
    )

    coordinate_rows = read_rows(
        SOURCE_COORDINATES
    )

    conformer_rows = read_rows(
        SOURCE_CONFORMERS
    )

    if graph_summary.get(
        "decision"
    ) != EXPECTED_GRAPH_DECISION:
        raise RuntimeError(
            "Gate 3I graph is not accepted."
        )

    if source_summary.get(
        "decision"
    ) != EXPECTED_SOURCE_DECISION:
        raise RuntimeError(
            "Gate 3K does not contain the expected "
            "bridge-topology review decision."
        )

    if len(path_rows) != EXPECTED_TOTAL_PATHS:
        raise RuntimeError(
            "Unexpected number of bridge paths: "
            f"{len(path_rows)}/"
            f"{EXPECTED_TOTAL_PATHS}"
        )

    nodes = {
        row[
            "node_id"
        ]: row
        for row in node_rows
    }

    coordinates = {
        row[
            "node_id"
        ]: np.asarray(
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
        for row in coordinate_rows
    }

    if set(nodes) != set(coordinates):
        raise RuntimeError(
            "Graph nodes and refined coordinates disagree."
        )

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

    conformers_by_path = {
        row[
            "bridge_path_id"
        ]: row
        for row in conformer_rows
    }

    paths_by_id = {
        row[
            "bridge_path_id"
        ]: row
        for row in path_rows
    }

    if set(conformers_by_path) != set(
        paths_by_id
    ):
        raise RuntimeError(
            "Bridge-path and refined-conformer identifiers disagree."
        )

    angle_fields = (
        "internal_angle_1_deg",
        "internal_angle_2_deg",
        "internal_angle_3_deg",
        "seed_junction_angle_1_deg",
        "seed_junction_angle_2_deg",
        "annulus_junction_angle_1_deg",
        "annulus_junction_angle_2_deg",
    )

    failed_path_ids = []

    original_failure_rows = []

    for path_id, row in conformers_by_path.items():
        values = [
            parse_float(
                row,
                field,
            )
            for field in angle_fields
        ]

        failures = [
            field
            for field, value
            in zip(
                angle_fields,
                values,
            )
            if (
                value
                < MIN_ACCEPTED_ANGLE_DEG
                or value
                > MAX_ACCEPTED_ANGLE_DEG
            )
        ]

        if failures:
            failed_path_ids.append(
                path_id
            )

            original_failure_rows.append(
                {
                    "bridge_path_id": (
                        path_id
                    ),
                    "end": row[
                        "end"
                    ],
                    "bridge_index": row[
                        "bridge_index"
                    ],
                    "original_failed_angle_fields": (
                        " | ".join(
                            failures
                        )
                    ),
                    "original_minimum_angle_deg": (
                        min(values)
                    ),
                    "original_maximum_angle_deg": (
                        max(values)
                    ),
                }
            )

    failed_path_ids.sort()

    if not failed_path_ids:
        raise RuntimeError(
            "No failed trimer paths were found in Gate 3K."
        )

    print(
        "Failed trimer paths identified: "
        f"{len(failed_path_ids)}"
    )

    library, library_metadata = (
        build_library()
    )

    endpoint_distances = np.linalg.norm(
        library[:, -1, :]
        - library[:, 0, :],
        axis=1,
    )

    library_summary_row = {
        "library_conformers": int(
            library.shape[0]
        ),
        "independent_internal_angle_values": (
            len(
                ANGLE_VALUES_DEG
            )
        ),
        "angle_values_deg": (
            " | ".join(
                f"{value:.1f}"
                for value
                in ANGLE_VALUES_DEG
            )
        ),
        "torsion_values": (
            len(
                TORSION_VALUES_DEG
            )
        ),
        "torsion_step_deg": 30,
        "minimum_endpoint_distance_nm": float(
            np.min(
                endpoint_distances
            )
        ),
        "maximum_endpoint_distance_nm": float(
            np.max(
                endpoint_distances
            )
        ),
        "base_conformers_retained_per_path": (
            BASE_CONFORMERS_PER_PATH
        ),
        "azimuth_values_per_mirror": (
            len(
                AZIMUTH_VALUES_DEG
            )
        ),
        "mirrors_screened": 2,
    }

    write_rows(
        LIBRARY_SUMMARY,
        [
            library_summary_row
        ],
    )

    heavy_ids = [
        node_id
        for node_id, row
        in nodes.items()
        if row[
            "element"
        ]
        != "H"
    ]

    bridge_node_ids = {
        node_id
        for node_id, row
        in nodes.items()
        if row[
            "node_type"
        ]
        == "ALTERNATING_BN_TRIMER_BRIDGE"
    }

    fixed_heavy_ids = [
        node_id
        for node_id in heavy_ids
        if node_id not in bridge_node_ids
    ]

    result_rows = []
    best_coordinate_rows = []

    for audit_index, path_id in enumerate(
        failed_path_ids,
        start=1,
    ):
        path = paths_by_id[
            path_id
        ]

        start_id = path[
            "seed_node"
        ]

        finish_id = path[
            "annulus_node"
        ]

        first_bridge = path[
            "bridge_node_1"
        ]

        third_bridge = path[
            "bridge_node_3"
        ]

        start = coordinates[
            start_id
        ]

        finish = coordinates[
            finish_id
        ]

        start_neighbors = [
            coordinates[
                neighbor
            ]
            for neighbor in adjacency[
                start_id
            ]
            if (
                nodes[
                    neighbor
                ][
                    "element"
                ]
                != "H"
                and neighbor
                != first_bridge
            )
        ]

        finish_neighbors = [
            coordinates[
                neighbor
            ]
            for neighbor in adjacency[
                finish_id
            ]
            if (
                nodes[
                    neighbor
                ][
                    "element"
                ]
                != "H"
                and neighbor
                != third_bridge
            )
        ]

        if (
            len(
                start_neighbors
            )
            != 2
            or len(
                finish_neighbors
            )
            != 2
        ):
            raise RuntimeError(
                f"{path_id}: endpoint heavy-neighbor count "
                "is not 2/2."
            )

        midpoint = (
            start
            + finish
        ) / 2.0

        local_fixed_ids = [
            node_id
            for node_id in fixed_heavy_ids
            if (
                node_id
                not in {
                    start_id,
                    finish_id,
                }
                and float(
                    np.linalg.norm(
                        coordinates[
                            node_id
                        ]
                        - midpoint
                    )
                )
                <= LOCAL_FIXED_RADIUS_NM
            )
        ]

        if local_fixed_ids:
            local_fixed = np.asarray(
                [
                    coordinates[
                        node_id
                    ]
                    for node_id
                    in local_fixed_ids
                ],
                dtype=float,
            )
        else:
            local_fixed = np.empty(
                (
                    0,
                    3,
                ),
                dtype=float,
            )

        target_distance = float(
            np.linalg.norm(
                finish
                - start
            )
        )

        distance_errors = np.abs(
            endpoint_distances
            - target_distance
        )

        retained_count = min(
            BASE_CONFORMERS_PER_PATH,
            distance_errors.size,
        )

        retained_indices = np.argpartition(
            distance_errors,
            retained_count - 1,
        )[
            :retained_count
        ]

        retained_indices = retained_indices[
            np.argsort(
                distance_errors[
                    retained_indices
                ]
            )
        ]

        tested_candidates = 0
        passing_candidates = 0
        best_candidate: dict[str, Any] | None = None
        best_chain: np.ndarray | None = None

        for library_index in retained_indices:
            base_chain = library[
                library_index
            ]

            metadata = library_metadata[
                int(
                    library_index
                )
            ]

            library_error = float(
                distance_errors[
                    library_index
                ]
            )

            for mirror in (
                False,
                True,
            ):
                for azimuth in AZIMUTH_VALUES_DEG:
                    mapped = map_chain(
                        base_chain,
                        start,
                        finish,
                        azimuth,
                        mirror,
                    )

                    metrics = evaluate_candidate(
                        mapped,
                        start_neighbors,
                        finish_neighbors,
                        local_fixed,
                    )

                    candidate = {
                        "library_index": int(
                            library_index
                        ),
                        "library_angle_1_deg": (
                            metadata[
                                "angle_1_deg"
                            ]
                        ),
                        "library_angle_2_deg": (
                            metadata[
                                "angle_2_deg"
                            ]
                        ),
                        "library_angle_3_deg": (
                            metadata[
                                "angle_3_deg"
                            ]
                        ),
                        "library_torsion_1_deg": (
                            metadata[
                                "torsion_1_deg"
                            ]
                        ),
                        "library_torsion_2_deg": (
                            metadata[
                                "torsion_2_deg"
                            ]
                        ),
                        "library_endpoint_distance_nm": (
                            metadata[
                                "endpoint_distance_nm"
                            ]
                        ),
                        "library_distance_error_nm": (
                            library_error
                        ),
                        "mirror": mirror,
                        "azimuth_deg": azimuth,
                        **metrics,
                    }

                    tested_candidates += 1

                    if bool(
                        candidate[
                            "passes"
                        ]
                    ):
                        passing_candidates += 1

                    if (
                        best_candidate is None
                        or candidate_score(
                            candidate
                        )
                        < candidate_score(
                            best_candidate
                        )
                    ):
                        best_candidate = (
                            candidate
                        )

                        best_chain = np.asarray(
                            mapped,
                            dtype=float,
                        )

        if (
            best_candidate is None
            or best_chain is None
        ):
            raise RuntimeError(
                f"No candidates evaluated for {path_id}"
            )

        original = next(
            row
            for row in original_failure_rows
            if row[
                "bridge_path_id"
            ]
            == path_id
        )

        result_rows.append(
            {
                **original,
                "target_endpoint_distance_nm": (
                    target_distance
                ),
                "local_fixed_heavy_atoms": (
                    len(
                        local_fixed_ids
                    )
                ),
                "base_library_conformers_screened": (
                    retained_count
                ),
                "mapped_candidates_tested": (
                    tested_candidates
                ),
                "passing_candidates": (
                    passing_candidates
                ),
                "local_geometry_solution_found": (
                    bool(
                        best_candidate[
                            "passes"
                        ]
                    )
                ),
                "best_library_angle_1_deg": (
                    best_candidate[
                        "library_angle_1_deg"
                    ]
                ),
                "best_library_angle_2_deg": (
                    best_candidate[
                        "library_angle_2_deg"
                    ]
                ),
                "best_library_angle_3_deg": (
                    best_candidate[
                        "library_angle_3_deg"
                    ]
                ),
                "best_library_torsion_1_deg": (
                    best_candidate[
                        "library_torsion_1_deg"
                    ]
                ),
                "best_library_torsion_2_deg": (
                    best_candidate[
                        "library_torsion_2_deg"
                    ]
                ),
                "best_mirror": (
                    best_candidate[
                        "mirror"
                    ]
                ),
                "best_azimuth_deg": (
                    best_candidate[
                        "azimuth_deg"
                    ]
                ),
                "best_library_distance_error_nm": (
                    best_candidate[
                        "library_distance_error_nm"
                    ]
                ),
                "best_maximum_bond_deviation_nm": (
                    best_candidate[
                        "maximum_bond_deviation_nm"
                    ]
                ),
                "best_bond_failure_count": (
                    best_candidate[
                        "bond_failure_count"
                    ]
                ),
                "best_minimum_angle_deg": (
                    best_candidate[
                        "minimum_angle_deg"
                    ]
                ),
                "best_maximum_angle_deg": (
                    best_candidate[
                        "maximum_angle_deg"
                    ]
                ),
                "best_angle_violations": (
                    best_candidate[
                        "angle_violations"
                    ]
                ),
                "best_angle_deficit_deg": (
                    best_candidate[
                        "angle_deficit_deg"
                    ]
                ),
                "best_clash_count": (
                    best_candidate[
                        "clash_count"
                    ]
                ),
                "best_minimum_clearance_nm": (
                    best_candidate[
                        "minimum_clearance_nm"
                    ]
                ),
                "best_internal_angle_1_deg": (
                    best_candidate[
                        "internal_angles"
                    ][0]
                ),
                "best_internal_angle_2_deg": (
                    best_candidate[
                        "internal_angles"
                    ][1]
                ),
                "best_internal_angle_3_deg": (
                    best_candidate[
                        "internal_angles"
                    ][2]
                ),
                "best_seed_junction_angle_1_deg": (
                    best_candidate[
                        "seed_angles"
                    ][0]
                ),
                "best_seed_junction_angle_2_deg": (
                    best_candidate[
                        "seed_angles"
                    ][1]
                ),
                "best_annulus_junction_angle_1_deg": (
                    best_candidate[
                        "annulus_angles"
                    ][0]
                ),
                "best_annulus_junction_angle_2_deg": (
                    best_candidate[
                        "annulus_angles"
                    ][1]
                ),
            }
        )

        bridge_ids = [
            path[
                "bridge_node_1"
            ],
            path[
                "bridge_node_2"
            ],
            path[
                "bridge_node_3"
            ],
        ]

        for bridge_position, (
            bridge_id,
            coordinate,
        ) in enumerate(
            zip(
                bridge_ids,
                best_chain[
                    1:4
                ],
            ),
            start=1,
        ):
            best_coordinate_rows.append(
                {
                    "bridge_path_id": (
                        path_id
                    ),
                    "end": path[
                        "end"
                    ],
                    "bridge_index": path[
                        "bridge_index"
                    ],
                    "bridge_position": (
                        bridge_position
                    ),
                    "bridge_node": (
                        bridge_id
                    ),
                    "candidate_is_local_pass": (
                        bool(
                            best_candidate[
                                "passes"
                            ]
                        )
                    ),
                    "x_nm": float(
                        coordinate[0]
                    ),
                    "y_nm": float(
                        coordinate[1]
                    ),
                    "z_nm": float(
                        coordinate[2]
                    ),
                    "coordinate_status": (
                        "SEARCH_CANDIDATE_ONLY_"
                        "NOT_APPLIED_TO_STRUCTURE"
                    ),
                }
            )

        print(
            f"  {audit_index:02d}/"
            f"{len(failed_path_ids):02d} "
            f"{path_id}: tested="
            f"{tested_candidates}; passing="
            f"{passing_candidates}; best angles="
            f"{best_candidate['minimum_angle_deg']:.3f}-"
            f"{best_candidate['maximum_angle_deg']:.3f}; "
            f"clashes={best_candidate['clash_count']}; "
            f"bond-dev="
            f"{best_candidate['maximum_bond_deviation_nm']:.6f}"
        )

    write_rows(
        PATH_RESULTS,
        result_rows,
    )

    write_rows(
        BEST_CANDIDATES,
        best_coordinate_rows,
    )

    solved_paths = [
        row
        for row in result_rows
        if bool(
            row[
                "local_geometry_solution_found"
            ]
        )
    ]

    unsolved_paths = [
        row
        for row in result_rows
        if not bool(
            row[
                "local_geometry_solution_found"
            ]
        )
    ]

    lower_failed = sum(
        row[
            "end"
        ]
        == "LOWER"
        for row in result_rows
    )

    upper_failed = sum(
        row[
            "end"
        ]
        == "UPPER"
        for row in result_rows
    )

    lower_solved = sum(
        row[
            "end"
        ]
        == "LOWER"
        and bool(
            row[
                "local_geometry_solution_found"
            ]
        )
        for row in result_rows
    )

    upper_solved = sum(
        row[
            "end"
        ]
        == "UPPER"
        and bool(
            row[
                "local_geometry_solution_found"
            ]
        )
        for row in result_rows
    )

    audit_gates = {
        "Gate3I_graph_is_accepted": (
            graph_summary.get(
                "decision"
            )
            == EXPECTED_GRAPH_DECISION
        ),
        "Gate3K_has_expected_topology_review_decision": (
            source_summary.get(
                "decision"
            )
            == EXPECTED_SOURCE_DECISION
        ),
        "all_Gate3K_angle_failed_paths_were_identified": (
            len(
                result_rows
            )
            == len(
                failed_path_ids
            )
            and len(
                failed_path_ids
            )
            > 0
        ),
        "independent_angle_library_is_nonempty": (
            library.shape[0]
            > 0
        ),
        "every_failed_path_received_expanded_search": all(
            int(
                row[
                    "mapped_candidates_tested"
                ]
            )
            > 0
            for row in result_rows
        ),
        "all_reported_search_metrics_are_finite": all(
            all(
                math.isfinite(
                    float(
                        row[field]
                    )
                )
                for field in (
                    "target_endpoint_distance_nm",
                    "best_library_distance_error_nm",
                    "best_maximum_bond_deviation_nm",
                    "best_minimum_angle_deg",
                    "best_maximum_angle_deg",
                    "best_angle_deficit_deg",
                    "best_minimum_clearance_nm",
                )
            )
            for row in result_rows
        ),
        "candidate_coordinates_are_search_only": all(
            row[
                "coordinate_status"
            ]
            == (
                "SEARCH_CANDIDATE_ONLY_"
                "NOT_APPLIED_TO_STRUCTURE"
            )
            for row in best_coordinate_rows
        ),
    }

    failed_audit_gates = [
        name
        for name, passed
        in audit_gates.items()
        if not passed
    ]

    audit_integrity_pass = (
        len(
            failed_audit_gates
        )
        == 0
    )

    all_failed_paths_solved = (
        audit_integrity_pass
        and len(
            unsolved_paths
        )
        == 0
    )

    decision = (
        RETAIN_DECISION
        if all_failed_paths_solved
        else REDESIGN_DECISION
    )

    required_next_step = (
        "RUN_R2_GLOBAL_TRIMER_BRIDGE_REFINEMENT_WITH_"
        "INDEPENDENT_INTERNAL_ANGLE_LIBRARY"
        if all_failed_paths_solved
        else
        "SCREEN_R2_SPARSE_AND_LONGER_BRIDGE_TOPOLOGIES"
    )

    summary = {
        "decision": decision,
        "original_angle_failed_paths": (
            len(
                result_rows
            )
        ),
        "lower_original_angle_failed_paths": (
            lower_failed
        ),
        "upper_original_angle_failed_paths": (
            upper_failed
        ),
        "locally_solved_paths": (
            len(
                solved_paths
            )
        ),
        "locally_unsolved_paths": (
            len(
                unsolved_paths
            )
        ),
        "lower_locally_solved_paths": (
            lower_solved
        ),
        "upper_locally_solved_paths": (
            upper_solved
        ),
        "independent_angle_library_conformers": int(
            library.shape[0]
        ),
        "base_conformers_retained_per_path": (
            BASE_CONFORMERS_PER_PATH
        ),
        "azimuth_values_per_mirror": (
            len(
                AZIMUTH_VALUES_DEG
            )
        ),
        "maximum_passing_candidate_count": max(
            int(
                row[
                    "passing_candidates"
                ]
            )
            for row in result_rows
        ),
        "minimum_passing_candidate_count": min(
            int(
                row[
                    "passing_candidates"
                ]
            )
            for row in result_rows
        ),
        "best_global_minimum_angle_deg": min(
            float(
                row[
                    "best_minimum_angle_deg"
                ]
            )
            for row in result_rows
        ),
        "best_global_maximum_bond_deviation_nm": max(
            float(
                row[
                    "best_maximum_bond_deviation_nm"
                ]
            )
            for row in result_rows
        ),
        "best_global_minimum_clearance_nm": min(
            float(
                row[
                    "best_minimum_clearance_nm"
                ]
            )
            for row in result_rows
        ),
        "audit_integrity_pass": (
            audit_integrity_pass
        ),
        "current_trimer_graph_retained": (
            all_failed_paths_solved
        ),
        "current_trimer_graph_redesign_confirmed": (
            not all_failed_paths_solved
        ),
        "global_coordinate_refinement_authorized": (
            all_failed_paths_solved
        ),
        "candidate_coordinates_applied": False,
        "molecular_topology_generation_authorized": False,
        "formal_charge_assignment_authorized": False,
        "force_field_parameterization_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "QM_authorized": False,
        "failed_audit_gates": (
            " | ".join(
                failed_audit_gates
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
            in audit_gates.items()
        ],
    )

    JSON_OUT.write_text(
        json.dumps(
            {
                "summary": summary,
                "library": (
                    library_summary_row
                ),
                "path_results": (
                    result_rows
                ),
                "audit_gates": (
                    audit_gates
                ),
                "limitations": [
                    (
                        "The expanded search varies the three "
                        "internal trimer angles independently."
                    ),
                    (
                        "The search evaluates local compatibility "
                        "against the fixed parent, seed and annulus "
                        "heavy-atom scaffold."
                    ),
                    (
                        "It does not yet optimize all 30 bridges "
                        "simultaneously."
                    ),
                    (
                        "Candidate coordinates are diagnostic only "
                        "and are not applied to the accepted structure."
                    ),
                    (
                        "No topology, charges, force-field parameters, "
                        "minimization, MD or QM calculation was generated."
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
            "role": (
                "Gate3I_graph_nodes"
            ),
            "file": relative(
                GRAPH_NODES
            ),
            "sha256": sha256(
                GRAPH_NODES
            ),
        },
        {
            "role": (
                "Gate3I_graph_edges"
            ),
            "file": relative(
                GRAPH_EDGES
            ),
            "sha256": sha256(
                GRAPH_EDGES
            ),
        },
        {
            "role": (
                "Gate3I_bridge_paths"
            ),
            "file": relative(
                BRIDGE_PATHS
            ),
            "sha256": sha256(
                BRIDGE_PATHS
            ),
        },
        {
            "role": (
                "Gate3I_graph_summary"
            ),
            "file": relative(
                GRAPH_SUMMARY
            ),
            "sha256": sha256(
                GRAPH_SUMMARY
            ),
        },
        {
            "role": (
                "Gate3K_refined_coordinates"
            ),
            "file": relative(
                SOURCE_COORDINATES
            ),
            "sha256": sha256(
                SOURCE_COORDINATES
            ),
        },
        {
            "role": (
                "Gate3K_refined_conformers"
            ),
            "file": relative(
                SOURCE_CONFORMERS
            ),
            "sha256": sha256(
                SOURCE_CONFORMERS
            ),
        },
        {
            "role": (
                "Gate3K_refinement_summary"
            ),
            "file": relative(
                SOURCE_SUMMARY
            ),
            "sha256": sha256(
                SOURCE_SUMMARY
            ),
        },
    ]

    write_rows(
        MANIFEST,
        manifest_rows,
    )

    path_lines = "\n".join(
        (
            f"- `{row['bridge_path_id']}`: "
            f"passing candidates="
            f"{row['passing_candidates']}; "
            f"best angle range="
            f"{float(row['best_minimum_angle_deg']):.3f}–"
            f"{float(row['best_maximum_angle_deg']):.3f}°; "
            f"best max bond deviation="
            f"{float(row['best_maximum_bond_deviation_nm']):.6f} nm; "
            f"best clashes="
            f"{row['best_clash_count']}; "
            f"solved="
            f"{row['local_geometry_solution_found']}"
        )
        for row in result_rows
    )

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed
        in audit_gates.items()
    )

    REPORT.write_text(
        f"""# R2 Trimer-Bridge Search-Completeness Audit

## Scope

This gate tests whether the Gate 3K angle failures result from the
trimer topology itself or from the restricted equal-angle conformer
library previously used.

The three internal trimer angles are varied independently. Candidate
coordinates are diagnostic only and are not applied to the accepted
structure.

## Search space

- Independent angle values:
  **{len(ANGLE_VALUES_DEG)}**
- Torsion combinations:
  **{len(TORSION_VALUES_DEG) ** 2}**
- Total library conformers:
  **{library.shape[0]}**
- Base conformers retained per failed path:
  **{BASE_CONFORMERS_PER_PATH}**
- Azimuths per mirror:
  **{len(AZIMUTH_VALUES_DEG)}**
- Mirrors:
  **2**

## Original failed paths

- Total:
  **{len(result_rows)}**
- Lower/upper:
  **{lower_failed}/{upper_failed}**

## Expanded-search results

{path_lines}

## Aggregate result

- Locally solved paths:
  **{len(solved_paths)}**
- Locally unsolved paths:
  **{len(unsolved_paths)}**
- Lower solved:
  **{lower_solved}/{lower_failed}**
- Upper solved:
  **{upper_solved}/{upper_failed}**

## Audit gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed audit-integrity gates:
  **{'NONE' if not failed_audit_gates else ' | '.join(failed_audit_gates)}**
- Current trimer graph retained:
  **{'YES' if all_failed_paths_solved else 'NO'}**
- Global coordinate refinement authorized:
  **{'YES' if all_failed_paths_solved else 'NO'}**
- Candidate coordinates applied:
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
        "Day024 R2 trimer-bridge search-completeness "
        "audit completed."
    )

    print(
        "Independent-angle library conformers: "
        f"{library.shape[0]}"
    )

    print(
        "Original failed paths lower/upper/total: "
        f"{lower_failed}/"
        f"{upper_failed}/"
        f"{len(result_rows)}"
    )

    print(
        "Locally solved paths lower/upper/total: "
        f"{lower_solved}/"
        f"{upper_solved}/"
        f"{len(solved_paths)}"
    )

    print(
        "Locally unsolved paths: "
        f"{len(unsolved_paths)}"
    )

    print(
        "Best global minimum angle / "
        "max bond deviation / minimum clearance: "
        f"{summary['best_global_minimum_angle_deg']:.3f}/"
        f"{summary['best_global_maximum_bond_deviation_nm']:.6f}/"
        f"{summary['best_global_minimum_clearance_nm']:.6f}"
    )

    print(
        f"Decision: {decision}"
    )

    print(
        "Failed audit-integrity gates: "
        + (
            "NONE"
            if not failed_audit_gates
            else " | ".join(
                failed_audit_gates
            )
        )
    )

    print(
        "Current trimer graph retained: "
        f"{'YES' if all_failed_paths_solved else 'NO'}"
    )

    print(
        "Global coordinate refinement authorized: "
        f"{'YES' if all_failed_paths_solved else 'NO'}"
    )

    print(
        "Candidate coordinates applied: NO"
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
        LIBRARY_SUMMARY,
        PATH_RESULTS,
        BEST_CANDIDATES,
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
