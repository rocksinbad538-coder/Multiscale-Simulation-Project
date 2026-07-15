#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/phase1A/day024_chemical_end_rim_design"

G3A = BASE / "01_r2_parent_rim_chemical_audit"
G3I = BASE / "11_r2_alternating_bn_trimer_bridge_graph"
G3J = BASE / "12_r2_alternating_bn_trimer_bridge_static_coordinate_embedding"

OUT = (
    BASE
    / "13_r2_trimer_bridge_conformer_and_h_refinement"
)

PARENT_SUMMARY = (
    G3A
    / "r2_parent_rim_chemical_audit_summary.csv"
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

GRAPH_GATES = (
    G3I
    / "r2_alternating_bn_trimer_bridge_graph_gates.csv"
)

SOURCE_COORDINATES = (
    G3J
    / "r2_trimer_bridge_static_coordinates.csv"
)

SOURCE_EMBEDDING_SUMMARY = (
    G3J
    / "r2_trimer_bridge_static_embedding_summary.csv"
)

REFINED_COORDINATES = (
    OUT
    / "r2_trimer_bridge_refined_coordinates.csv"
)

REFINED_CONFORMERS = (
    OUT
    / "r2_trimer_bridge_refined_conformers.csv"
)

REFINED_HYDROGENS = (
    OUT
    / "r2_trimer_bridge_refined_hydrogen_orientations.csv"
)

BOND_LENGTHS = (
    OUT
    / "r2_trimer_bridge_refined_bond_lengths.csv"
)

BOND_SUMMARY = (
    OUT
    / "r2_trimer_bridge_refined_bond_type_summary.csv"
)

ANGLE_SUMMARY = (
    OUT
    / "r2_trimer_bridge_refined_angle_summary.csv"
)

CONTACT_SUMMARY = (
    OUT
    / "r2_trimer_bridge_refined_nonbonded_contact_summary.csv"
)

END_SUMMARY = (
    OUT
    / "r2_trimer_bridge_refined_end_summary.csv"
)

SUMMARY = (
    OUT
    / "r2_trimer_bridge_refinement_summary.csv"
)

GATES = (
    OUT
    / "r2_trimer_bridge_refinement_gates.csv"
)

JSON_OUT = (
    OUT
    / "r2_trimer_bridge_refinement.json"
)

MANIFEST = (
    OUT
    / "r2_trimer_bridge_refinement_source_manifest.csv"
)

XYZ_OUT = (
    OUT
    / "r2_trimer_bridge_refined_embedding.xyz"
)

PDB_OUT = (
    OUT
    / "r2_trimer_bridge_refined_embedding.pdb"
)

REPORT = (
    OUT
    / "R2_TRIMER_BRIDGE_CONFORMER_AND_H_REFINEMENT_DAY024.md"
)

EXPECTED_PARENT_DECISION = (
    "R2_PARENT_RIM_AND_CHEMICAL_CONSTRAINTS_AUDITED"
)

EXPECTED_GRAPH_DECISION = (
    "R2_ALTERNATING_BN_TRIMER_BRIDGE_GRAPH_VALIDATED"
)

EXPECTED_SOURCE_DECISION = (
    "R2_ALTERNATING_BN_TRIMER_BRIDGE_STATIC_"
    "COORDINATE_EMBEDDING_REQUIRES_CONFORMER_REFINEMENT"
)

PASS_DECISION = (
    "R2_ALTERNATING_BN_TRIMER_BRIDGE_CONFORMER_"
    "AND_H_REFINEMENT_VALIDATED"
)

HEAVY_REVIEW_DECISION = (
    "R2_ALTERNATING_BN_TRIMER_BRIDGE_REFINEMENT_"
    "REQUIRES_BRIDGE_TOPOLOGY_REDESIGN"
)

H_REVIEW_DECISION = (
    "R2_ALTERNATING_BN_TRIMER_BRIDGE_REFINEMENT_"
    "REQUIRES_SECONDARY_H_ORIENTATION_OPTIMIZATION"
)

N_PARENT = 1680
N_HEAVY = 2082
N_H = 174
N_TOTAL = 2256
N_BRIDGES = 30

BN = 0.144973
BH = 0.119
NH = 0.101

LIBRARY_ANGLE_MIN = 105.0
LIBRARY_ANGLE_MAX = 135.0
LIBRARY_ANGLE_STEP = 1.0
LIBRARY_TORSION_STEP = 10.0

BASE_CONFORMERS_PER_BRIDGE = 80
AZIMUTH_STEP = 10.0
BRIDGE_POOL_SIZE = 96
BRIDGE_GLOBAL_SWEEPS = 12

H_DIRECTION_COUNT = 362
H_POOL_SIZE = 48
H_GLOBAL_SWEEPS = 10

LOCAL_FIXED_RADIUS_NM = 0.80

MAX_LIBRARY_DISTANCE_ERROR_NM = 0.0005
MAX_BN_DEVIATION_NM = 0.003
MAX_XH_DEVIATION_NM = 0.002

MIN_CRITICAL_ANGLE_DEG = 70.0
MAX_CRITICAL_ANGLE_DEG = 175.0
MAX_CRITICAL_RMS_DEVIATION_DEG = 30.0

MIN_NONBONDED_HEAVY_HEAVY_NM = 0.120
MIN_NONBONDED_H_HEAVY_NM = 0.070
MIN_NONBONDED_H_H_NM = 0.060

SOFT_HEAVY_CLEARANCE_NM = 0.140
SOFT_H_HEAVY_CLEARANCE_NM = 0.090
SOFT_H_H_CLEARANCE_NM = 0.080

MAX_APERTURE_ERROR = 0.10
MAX_OUTER_RADIUS_ERROR = 0.15
MAX_END_ASYMMETRY_NM = 0.010

CRITICAL_NODE_TYPES = {
    "HEXAGONAL_EDGE_COMPLETION_SEED",
    "ALTERNATING_BN_TRIMER_BRIDGE",
    "ANNULUS_OUTER_BOUNDARY",
    "ANNULUS_INNER_BOUNDARY",
}


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

    fieldnames: list[str] = []

    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(
                {
                    key: row.get(key, "")
                    for key in fieldnames
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


def determine_axis(
    positions: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:
    center = np.mean(
        positions,
        axis=0,
    )

    centered = (
        positions
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


def rotation_about_axis(
    axis: np.ndarray,
    angle_rad: float,
) -> np.ndarray:
    x_value, y_value, z_value = normalized(
        axis
    )

    cosine = math.cos(
        angle_rad
    )

    sine = math.sin(
        angle_rad
    )

    factor = (
        1.0
        - cosine
    )

    return np.asarray(
        [
            [
                cosine
                + x_value
                * x_value
                * factor,
                x_value
                * y_value
                * factor
                - z_value
                * sine,
                x_value
                * z_value
                * factor
                + y_value
                * sine,
            ],
            [
                y_value
                * x_value
                * factor
                + z_value
                * sine,
                cosine
                + y_value
                * y_value
                * factor,
                y_value
                * z_value
                * factor
                - x_value
                * sine,
            ],
            [
                z_value
                * x_value
                * factor
                - y_value
                * sine,
                z_value
                * y_value
                * factor
                + x_value
                * sine,
                cosine
                + z_value
                * z_value
                * factor,
            ],
        ],
        dtype=float,
    )


def rotation_from_vectors(
    source: np.ndarray,
    target: np.ndarray,
) -> np.ndarray:
    source_unit = normalized(
        source
    )

    target_unit = normalized(
        target
    )

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
        third
        - second
    )

    normal = np.cross(
        second
        - first,
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

    normal = normalized(
        normal
    )

    in_plane = normalized(
        np.cross(
            normal,
            previous,
        )
    )

    direction = (
        -math.cos(
            angle_rad
        )
        * previous
        + math.sin(
            angle_rad
        )
        * (
            math.cos(
                dihedral_rad
            )
            * in_plane
            + math.sin(
                dihedral_rad
            )
            * normal
        )
    )

    return (
        third
        + BN
        * direction
    )


def canonical_chain(
    angle_deg: float,
    phi1_deg: float,
    phi2_deg: float,
) -> np.ndarray:
    angle_rad = math.radians(
        angle_deg
    )

    point_0 = np.zeros(
        3,
        dtype=float,
    )

    point_1 = np.asarray(
        [
            BN,
            0.0,
            0.0,
        ],
        dtype=float,
    )

    point_2 = (
        point_1
        + BN
        * np.asarray(
            [
                -math.cos(
                    angle_rad
                ),
                math.sin(
                    angle_rad
                ),
                0.0,
            ],
            dtype=float,
        )
    )

    point_3 = place_next(
        point_0,
        point_1,
        point_2,
        angle_rad,
        math.radians(
            phi1_deg
        ),
    )

    point_4 = place_next(
        point_1,
        point_2,
        point_3,
        angle_rad,
        math.radians(
            phi2_deg
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


def build_conformer_library() -> dict[str, np.ndarray]:
    angles = []
    torsions_1 = []
    torsions_2 = []
    distances = []

    for angle in np.arange(
        LIBRARY_ANGLE_MIN,
        LIBRARY_ANGLE_MAX
        + 0.5
        * LIBRARY_ANGLE_STEP,
        LIBRARY_ANGLE_STEP,
    ):
        for phi1 in np.arange(
            0.0,
            360.0,
            LIBRARY_TORSION_STEP,
        ):
            for phi2 in np.arange(
                0.0,
                360.0,
                LIBRARY_TORSION_STEP,
            ):
                chain = canonical_chain(
                    float(angle),
                    float(phi1),
                    float(phi2),
                )

                angles.append(
                    float(angle)
                )

                torsions_1.append(
                    float(phi1)
                )

                torsions_2.append(
                    float(phi2)
                )

                distances.append(
                    float(
                        np.linalg.norm(
                            chain[-1]
                            - chain[0]
                        )
                    )
                )

    return {
        "angle": np.asarray(
            angles,
            dtype=float,
        ),
        "phi1": np.asarray(
            torsions_1,
            dtype=float,
        ),
        "phi2": np.asarray(
            torsions_2,
            dtype=float,
        ),
        "distance": np.asarray(
            distances,
            dtype=float,
        ),
    }


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
        finish
        - start
    )

    alignment = rotation_from_vectors(
        working[-1]
        - working[0],
        endpoint_vector,
    )

    mapped = (
        working
        - working[0]
    ) @ alignment.T

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
                    first
                    - center
                ),
                normalized(
                    second
                    - center
                ),
            ),
            -1.0,
            1.0,
        )
    )

    return math.degrees(
        math.acos(cosine)
    )


def xh_target(
    element: str,
) -> float:
    if element == "B":
        return BH

    if element == "N":
        return NH

    raise RuntimeError(
        f"Unexpected X-H parent element: {element}"
    )


def clearance_penalty(
    distances: np.ndarray,
    soft_threshold: float,
) -> float:
    if distances.size == 0:
        return 0.0

    return float(
        np.sum(
            np.maximum(
                soft_threshold
                - distances,
                0.0,
            )
            ** 2
        )
    )


def candidate_local_metrics(
    chain: np.ndarray,
    start_existing_neighbors: list[np.ndarray],
    finish_existing_neighbors: list[np.ndarray],
    local_fixed_positions: np.ndarray,
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
        - BN
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

    endpoint_angles = []

    for neighbor in start_existing_neighbors:
        endpoint_angles.append(
            angle_degrees(
                neighbor,
                chain[0],
                chain[1],
            )
        )

    for neighbor in finish_existing_neighbors:
        endpoint_angles.append(
            angle_degrees(
                neighbor,
                chain[4],
                chain[3],
            )
        )

    critical_angles = np.asarray(
        [
            *internal_angles,
            *endpoint_angles,
        ],
        dtype=float,
    )

    angle_violations = int(
        np.sum(
            (
                critical_angles
                < MIN_CRITICAL_ANGLE_DEG
            )
            | (
                critical_angles
                > MAX_CRITICAL_ANGLE_DEG
            )
        )
    )

    angle_penalty = float(
        np.sum(
            (
                critical_angles
                - 120.0
            )
            ** 2
        )
    )

    internal = chain[
        1:4
    ]

    if local_fixed_positions.size:
        fixed_distances = np.linalg.norm(
            internal[:, None, :]
            - local_fixed_positions[
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

    own_nonbonded_pairs = (
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
            in own_nonbonded_pairs
        ],
        dtype=float,
    )

    all_clearances = np.concatenate(
        (
            fixed_distances,
            own_distances,
        )
    )

    local_clashes = int(
        np.sum(
            all_clearances
            < MIN_NONBONDED_HEAVY_HEAVY_NM
        )
    )

    minimum_clearance = float(
        np.min(
            all_clearances
        )
    )

    soft_penalty = clearance_penalty(
        all_clearances,
        SOFT_HEAVY_CLEARANCE_NM,
    )

    bond_failures = int(
        np.sum(
            bond_deviations
            > MAX_BN_DEVIATION_NM
        )
    )

    return {
        "bond_lengths": (
            bond_lengths
        ),
        "bond_deviations": (
            bond_deviations
        ),
        "internal_angles": (
            internal_angles
        ),
        "endpoint_angles": (
            endpoint_angles
        ),
        "angle_violations": (
            angle_violations
        ),
        "angle_penalty": (
            angle_penalty
        ),
        "local_clashes": (
            local_clashes
        ),
        "minimum_clearance_nm": (
            minimum_clearance
        ),
        "clearance_penalty": (
            soft_penalty
        ),
        "bond_failures": (
            bond_failures
        ),
    }


def bridge_local_score(
    candidate: dict[str, Any],
) -> tuple[Any, ...]:
    return (
        int(
            candidate[
                "angle_violations"
            ]
        ),
        int(
            candidate[
                "local_clashes"
            ]
        ),
        int(
            candidate[
                "bond_failures"
            ]
        ),
        float(
            candidate[
                "clearance_penalty"
            ]
        ),
        float(
            candidate[
                "angle_penalty"
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


def global_bridge_score(
    candidate: dict[str, Any],
    other_internal_positions: np.ndarray,
) -> tuple[Any, ...]:
    internal = candidate[
        "chain"
    ][
        1:4
    ]

    if other_internal_positions.size:
        distances = np.linalg.norm(
            internal[:, None, :]
            - other_internal_positions[
                None,
                :,
                :,
            ],
            axis=2,
        ).reshape(-1)

        inter_clashes = int(
            np.sum(
                distances
                < MIN_NONBONDED_HEAVY_HEAVY_NM
            )
        )

        inter_penalty = clearance_penalty(
            distances,
            SOFT_HEAVY_CLEARANCE_NM,
        )

        inter_minimum = float(
            np.min(
                distances
            )
        )
    else:
        inter_clashes = 0
        inter_penalty = 0.0
        inter_minimum = math.inf

    total_clashes = (
        int(
            candidate[
                "local_clashes"
            ]
        )
        + inter_clashes
    )

    minimum_clearance = min(
        float(
            candidate[
                "minimum_clearance_nm"
            ]
        ),
        inter_minimum,
    )

    return (
        int(
            candidate[
                "angle_violations"
            ]
        ),
        total_clashes,
        int(
            candidate[
                "bond_failures"
            ]
        ),
        float(
            candidate[
                "clearance_penalty"
            ]
        )
        + inter_penalty,
        float(
            candidate[
                "angle_penalty"
            ]
        ),
        -minimum_clearance,
        float(
            candidate[
                "library_distance_error_nm"
            ]
        ),
    )


def fibonacci_directions(
    count: int,
) -> np.ndarray:
    directions = []

    golden_angle = (
        math.pi
        * (
            3.0
            - math.sqrt(5.0)
        )
    )

    for index in range(count):
        z_value = (
            1.0
            - 2.0
            * (
                index
                + 0.5
            )
            / count
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


def preferred_h_direction(
    hydrogen_row: dict[str, str],
    heavy_position: np.ndarray,
    heavy_neighbor_positions: list[np.ndarray],
    tube_center: np.ndarray,
    tube_axis: np.ndarray,
    annulus_center_by_end: dict[
        str,
        np.ndarray,
    ],
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

    opposite_bisector = normalized(
        -vector_sum,
        np.asarray(
            [
                1.0,
                0.0,
                0.0,
            ],
            dtype=float,
        ),
    )

    end = hydrogen_row[
        "end"
    ]

    node_type = hydrogen_row[
        "node_type"
    ]

    outward = (
        -tube_axis
        if end == "LOWER"
        else tube_axis
    )

    annulus_center = (
        annulus_center_by_end[
            end
        ]
    )

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

    radial = normalized(
        radial,
        np.asarray(
            [
                1.0,
                0.0,
                0.0,
            ],
            dtype=float,
        ),
    )

    if node_type == "ANNULUS_INNER_PASSIVANT_H":
        role_direction = -radial

    elif node_type == "ANNULUS_OUTER_PASSIVANT_H":
        role_direction = radial

    elif node_type == "SEED_PASSIVANT_H":
        role_direction = outward

    else:
        role_direction = radial

    return normalized(
        (
            0.80
            * opposite_bisector
            + 0.20
            * role_direction
        ),
        opposite_bisector,
    )


def hydrogen_local_metrics(
    point: np.ndarray,
    center: np.ndarray,
    heavy_neighbor_positions: list[np.ndarray],
    all_heavy_positions: np.ndarray,
    attached_heavy_index: int,
    preferred_direction: np.ndarray,
) -> dict[str, Any]:
    direction = normalized(
        point
        - center
    )

    angles = np.asarray(
        [
            angle_degrees(
                neighbor,
                center,
                point,
            )
            for neighbor
            in heavy_neighbor_positions
        ],
        dtype=float,
    )

    angle_violations = int(
        np.sum(
            (
                angles
                < MIN_CRITICAL_ANGLE_DEG
            )
            | (
                angles
                > MAX_CRITICAL_ANGLE_DEG
            )
        )
    )

    angle_penalty = float(
        np.sum(
            (
                angles
                - 120.0
            )
            ** 2
        )
    )

    distances = np.linalg.norm(
        all_heavy_positions
        - point,
        axis=1,
    )

    distances = np.delete(
        distances,
        attached_heavy_index,
    )

    heavy_clashes = int(
        np.sum(
            distances
            < MIN_NONBONDED_H_HEAVY_NM
        )
    )

    heavy_penalty = clearance_penalty(
        distances,
        SOFT_H_HEAVY_CLEARANCE_NM,
    )

    minimum_heavy_clearance = float(
        np.min(
            distances
        )
    )

    directional_penalty = (
        1.0
        - float(
            np.clip(
                np.dot(
                    direction,
                    preferred_direction,
                ),
                -1.0,
                1.0,
            )
        )
    )

    return {
        "angles": (
            angles
        ),
        "angle_violations": (
            angle_violations
        ),
        "angle_penalty": (
            angle_penalty
        ),
        "heavy_clashes": (
            heavy_clashes
        ),
        "heavy_penalty": (
            heavy_penalty
        ),
        "minimum_heavy_clearance_nm": (
            minimum_heavy_clearance
        ),
        "directional_penalty": (
            directional_penalty
        ),
    }


def hydrogen_local_score(
    candidate: dict[str, Any],
) -> tuple[Any, ...]:
    return (
        int(
            candidate[
                "angle_violations"
            ]
        ),
        int(
            candidate[
                "heavy_clashes"
            ]
        ),
        float(
            candidate[
                "heavy_penalty"
            ]
        ),
        float(
            candidate[
                "angle_penalty"
            ]
        ),
        float(
            candidate[
                "directional_penalty"
            ]
        ),
        -float(
            candidate[
                "minimum_heavy_clearance_nm"
            ]
        ),
    )


def hydrogen_global_score(
    candidate: dict[str, Any],
    other_h_positions: np.ndarray,
) -> tuple[Any, ...]:
    point = candidate[
        "point"
    ]

    if other_h_positions.size:
        distances = np.linalg.norm(
            other_h_positions
            - point,
            axis=1,
        )

        hh_clashes = int(
            np.sum(
                distances
                < MIN_NONBONDED_H_H_NM
            )
        )

        hh_penalty = clearance_penalty(
            distances,
            SOFT_H_H_CLEARANCE_NM,
        )

        hh_minimum = float(
            np.min(
                distances
            )
        )
    else:
        hh_clashes = 0
        hh_penalty = 0.0
        hh_minimum = math.inf

    total_clashes = (
        int(
            candidate[
                "heavy_clashes"
            ]
        )
        + hh_clashes
    )

    return (
        int(
            candidate[
                "angle_violations"
            ]
        ),
        total_clashes,
        float(
            candidate[
                "heavy_penalty"
            ]
        )
        + hh_penalty,
        float(
            candidate[
                "angle_penalty"
            ]
        ),
        float(
            candidate[
                "directional_penalty"
            ]
        ),
        -min(
            float(
                candidate[
                    "minimum_heavy_clearance_nm"
                ]
            ),
            hh_minimum,
        ),
    )


def summarize_bond_types(
    bond_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    groups: dict[
        str,
        list[dict[str, Any]]
    ] = defaultdict(list)

    for row in bond_rows:
        groups[
            str(
                row[
                    "edge_type"
                ]
            )
        ].append(row)

    summaries = []

    for edge_type in sorted(groups):
        rows = groups[
            edge_type
        ]

        lengths = np.asarray(
            [
                float(
                    row[
                        "length_nm"
                    ]
                )
                for row in rows
            ],
            dtype=float,
        )

        deviations = np.asarray(
            [
                float(
                    row[
                        "deviation_nm"
                    ]
                )
                for row in rows
            ],
            dtype=float,
        )

        summaries.append(
            {
                "edge_type": edge_type,
                "edge_count": len(rows),
                "mean_length_nm": float(
                    np.mean(lengths)
                ),
                "minimum_length_nm": float(
                    np.min(lengths)
                ),
                "maximum_length_nm": float(
                    np.max(lengths)
                ),
                "RMS_deviation_nm": float(
                    np.sqrt(
                        np.mean(
                            deviations
                            ** 2
                        )
                    )
                ),
                "maximum_absolute_deviation_nm": float(
                    np.max(
                        np.abs(
                            deviations
                        )
                    )
                ),
            }
        )

    return summaries


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        PARENT_SUMMARY,
        GRAPH_NODES,
        GRAPH_EDGES,
        BRIDGE_PATHS,
        GRAPH_SUMMARY,
        GRAPH_GATES,
        SOURCE_COORDINATES,
        SOURCE_EMBEDDING_SUMMARY,
    ):
        require_file(required)

    parent_summary = read_one(
        PARENT_SUMMARY
    )

    graph_summary = read_one(
        GRAPH_SUMMARY
    )

    source_summary = read_one(
        SOURCE_EMBEDDING_SUMMARY
    )

    graph_gate_rows = read_rows(
        GRAPH_GATES
    )

    node_rows = read_rows(
        GRAPH_NODES
    )

    edge_rows = read_rows(
        GRAPH_EDGES
    )

    bridge_path_rows = read_rows(
        BRIDGE_PATHS
    )

    source_coordinate_rows = read_rows(
        SOURCE_COORDINATES
    )

    if parent_summary.get(
        "decision"
    ) != EXPECTED_PARENT_DECISION:
        raise RuntimeError(
            "Gate 3A is not accepted."
        )

    if graph_summary.get(
        "decision"
    ) != EXPECTED_GRAPH_DECISION:
        raise RuntimeError(
            "Gate 3I is not accepted."
        )

    if source_summary.get(
        "decision"
    ) != EXPECTED_SOURCE_DECISION:
        raise RuntimeError(
            "Gate 3J does not contain the expected "
            "conformer-refinement decision."
        )

    failed_graph_gates = [
        row.get(
            "gate",
            "",
        )
        for row in graph_gate_rows
        if not parse_bool(
            row.get(
                "pass",
                "false",
            )
        )
    ]

    if failed_graph_gates:
        raise RuntimeError(
            "Gate 3I still contains failed gates: "
            + " | ".join(
                failed_graph_gates
            )
        )

    if (
        len(node_rows) != N_TOTAL
        or len(bridge_path_rows) != N_BRIDGES
    ):
        raise RuntimeError(
            "Unexpected Gate 3I node or bridge count."
        )

    nodes = {
        row[
            "node_id"
        ]: row
        for row in node_rows
    }

    if len(nodes) != len(node_rows):
        raise RuntimeError(
            "Duplicate node identifiers."
        )

    adjacency: dict[
        str,
        set[str]
    ] = {
        node_id: set()
        for node_id in nodes
    }

    edge_pairs: set[
        tuple[str, str]
    ] = set()

    for row in edge_rows:
        first = row[
            "source_node"
        ]

        second = row[
            "target_node"
        ]

        if (
            first not in adjacency
            or second not in adjacency
        ):
            raise RuntimeError(
                "Graph edge references an unknown node."
            )

        pair = tuple(
            sorted(
                (
                    first,
                    second,
                )
            )
        )

        if pair in edge_pairs:
            raise RuntimeError(
                f"Duplicate graph edge: {pair}"
            )

        edge_pairs.add(pair)
        adjacency[first].add(second)
        adjacency[second].add(first)

    source_coordinates = {
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
        for row in source_coordinate_rows
    }

    if set(source_coordinates) != set(nodes):
        raise RuntimeError(
            "Gate 3I nodes and Gate 3J coordinates disagree."
        )

    parent_ids = sorted(
        (
            node_id
            for node_id, row
            in nodes.items()
            if row[
                "node_type"
            ]
            == "PARENT_HBN"
        ),
        key=lambda node_id: int(
            node_id.split(
                ":"
            )[1]
        ),
    )

    if len(parent_ids) != N_PARENT:
        raise RuntimeError(
            "Unexpected parent population."
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

    hydrogen_ids = [
        node_id
        for node_id, row
        in nodes.items()
        if row[
            "element"
        ]
        == "H"
    ]

    bridge_ids = {
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
        if node_id not in bridge_ids
    ]

    if (
        len(heavy_ids) != N_HEAVY
        or len(hydrogen_ids) != N_H
        or len(bridge_ids) != 90
    ):
        raise RuntimeError(
            "Unexpected heavy, H, or bridge populations."
        )

    positions: dict[
        str,
        np.ndarray
    ] = {
        node_id: np.array(
            source_coordinates[
                node_id
            ],
            dtype=float,
            copy=True,
        )
        for node_id in fixed_heavy_ids
    }

    parent_positions = np.asarray(
        [
            positions[
                node_id
            ]
            for node_id in parent_ids
        ],
        dtype=float,
    )

    tube_center, tube_axis = determine_axis(
        parent_positions
    )

    target_aperture_diameter_nm = parse_float(
        parent_summary,
        "target_aperture_diameter_nm",
    )

    target_outer_radius_nm = parse_float(
        parent_summary,
        "parent_rim_mean_radius_nm",
    )

    library = build_conformer_library()

    bridge_paths = sorted(
        bridge_path_rows,
        key=lambda row: (
            row[
                "end"
            ],
            parse_int(
                row,
                "bridge_index",
            ),
        ),
    )

    bridge_candidate_pools: dict[
        str,
        list[dict[str, Any]]
    ] = {}

    print(
        "Generating bridge conformer candidate pools..."
    )

    for path_index, path in enumerate(
        bridge_paths,
        start=1,
    ):
        path_id = path[
            "bridge_path_id"
        ]

        start_id = path[
            "seed_node"
        ]

        finish_id = path[
            "annulus_node"
        ]

        start = positions[
            start_id
        ]

        finish = positions[
            finish_id
        ]

        target_distance = float(
            np.linalg.norm(
                finish
                - start
            )
        )

        distance_errors = np.abs(
            library[
                "distance"
            ]
            - target_distance
        )

        ranking_metric = (
            distance_errors
            + 1.0e-8
            * np.abs(
                library[
                    "angle"
                ]
                - 120.0
            )
        )

        count = min(
            BASE_CONFORMERS_PER_BRIDGE,
            ranking_metric.size,
        )

        top_indices = np.argpartition(
            ranking_metric,
            count - 1,
        )[:count]

        top_indices = top_indices[
            np.argsort(
                ranking_metric[
                    top_indices
                ]
            )
        ]

        first_bridge_id = path[
            "bridge_node_1"
        ]

        third_bridge_id = path[
            "bridge_node_3"
        ]

        start_existing_neighbors = [
            positions[
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
                != first_bridge_id
            )
        ]

        finish_existing_neighbors = [
            positions[
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
                != third_bridge_id
            )
        ]

        if (
            len(
                start_existing_neighbors
            )
            != 2
            or len(
                finish_existing_neighbors
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

        fixed_local_ids = [
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
                        positions[
                            node_id
                        ]
                        - midpoint
                    )
                )
                <= LOCAL_FIXED_RADIUS_NM
            )
        ]

        local_fixed_positions = np.asarray(
            [
                positions[
                    node_id
                ]
                for node_id in fixed_local_ids
            ],
            dtype=float,
        )

        if local_fixed_positions.size == 0:
            local_fixed_positions = np.empty(
                (
                    0,
                    3,
                ),
                dtype=float,
            )

        candidates = []

        for library_index in top_indices:
            angle_value = float(
                library[
                    "angle"
                ][
                    library_index
                ]
            )

            phi1_value = float(
                library[
                    "phi1"
                ][
                    library_index
                ]
            )

            phi2_value = float(
                library[
                    "phi2"
                ][
                    library_index
                ]
            )

            base_chain = canonical_chain(
                angle_value,
                phi1_value,
                phi2_value,
            )

            library_error = float(
                distance_errors[
                    library_index
                ]
            )

            for mirror in (
                False,
                True,
            ):
                for azimuth in np.arange(
                    0.0,
                    360.0,
                    AZIMUTH_STEP,
                ):
                    chain = map_chain(
                        base_chain,
                        start,
                        finish,
                        float(
                            azimuth
                        ),
                        mirror,
                    )

                    metrics = candidate_local_metrics(
                        chain,
                        start_existing_neighbors,
                        finish_existing_neighbors,
                        local_fixed_positions,
                    )

                    candidate = {
                        "chain": chain,
                        "library_angle_deg": (
                            angle_value
                        ),
                        "library_phi1_deg": (
                            phi1_value
                        ),
                        "library_phi2_deg": (
                            phi2_value
                        ),
                        "library_distance_nm": float(
                            library[
                                "distance"
                            ][
                                library_index
                            ]
                        ),
                        "library_distance_error_nm": (
                            library_error
                        ),
                        "mirror": mirror,
                        "azimuth_deg": float(
                            azimuth
                        ),
                        **metrics,
                    }

                    candidates.append(
                        candidate
                    )

        candidates.sort(
            key=bridge_local_score
        )

        bridge_candidate_pools[
            path_id
        ] = candidates[
            :BRIDGE_POOL_SIZE
        ]

        best = bridge_candidate_pools[
            path_id
        ][0]

        print(
            f"  {path_index:02d}/{N_BRIDGES} "
            f"{path_id}: pool="
            f"{len(bridge_candidate_pools[path_id])}; "
            f"best angle violations/clashes="
            f"{best['angle_violations']}/"
            f"{best['local_clashes']}"
        )

    selected_bridge_index = {
        path[
            "bridge_path_id"
        ]: 0
        for path in bridge_paths
    }

    print(
        "Running global bridge conformer coordinate descent..."
    )

    for sweep in range(
        BRIDGE_GLOBAL_SWEEPS
    ):
        changed = 0

        ordered_paths = (
            bridge_paths
            if sweep % 2 == 0
            else list(
                reversed(
                    bridge_paths
                )
            )
        )

        for path in ordered_paths:
            path_id = path[
                "bridge_path_id"
            ]

            other_internal_positions = []

            for other_path in bridge_paths:
                other_id = other_path[
                    "bridge_path_id"
                ]

                if other_id == path_id:
                    continue

                selected_candidate = (
                    bridge_candidate_pools[
                        other_id
                    ][
                        selected_bridge_index[
                            other_id
                        ]
                    ]
                )

                other_internal_positions.extend(
                    selected_candidate[
                        "chain"
                    ][
                        1:4
                    ]
                )

            if other_internal_positions:
                other_array = np.asarray(
                    other_internal_positions,
                    dtype=float,
                )
            else:
                other_array = np.empty(
                    (
                        0,
                        3,
                    ),
                    dtype=float,
                )

            pool = bridge_candidate_pools[
                path_id
            ]

            scores = [
                global_bridge_score(
                    candidate,
                    other_array,
                )
                for candidate in pool
            ]

            best_index = min(
                range(
                    len(pool)
                ),
                key=lambda index: scores[
                    index
                ],
            )

            if (
                best_index
                != selected_bridge_index[
                    path_id
                ]
            ):
                selected_bridge_index[
                    path_id
                ] = best_index

                changed += 1

        print(
            f"  bridge sweep {sweep + 1}/"
            f"{BRIDGE_GLOBAL_SWEEPS}: "
            f"changed={changed}"
        )

        if changed == 0:
            break

    refined_conformer_rows = []

    for path in bridge_paths:
        path_id = path[
            "bridge_path_id"
        ]

        candidate = (
            bridge_candidate_pools[
                path_id
            ][
                selected_bridge_index[
                    path_id
                ]
            ]
        )

        bridge_node_ids = [
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

        for node_id, point in zip(
            bridge_node_ids,
            candidate[
                "chain"
            ][
                1:4
            ],
        ):
            positions[
                node_id
            ] = np.asarray(
                point,
                dtype=float,
            )

        refined_conformer_rows.append(
            {
                "bridge_path_id": path_id,
                "end": path[
                    "end"
                ],
                "bridge_index": path[
                    "bridge_index"
                ],
                "seed_node": path[
                    "seed_node"
                ],
                "annulus_node": path[
                    "annulus_node"
                ],
                "library_angle_deg": (
                    candidate[
                        "library_angle_deg"
                    ]
                ),
                "library_phi1_deg": (
                    candidate[
                        "library_phi1_deg"
                    ]
                ),
                "library_phi2_deg": (
                    candidate[
                        "library_phi2_deg"
                    ]
                ),
                "library_distance_nm": (
                    candidate[
                        "library_distance_nm"
                    ]
                ),
                "library_distance_error_nm": (
                    candidate[
                        "library_distance_error_nm"
                    ]
                ),
                "selected_mirror": (
                    candidate[
                        "mirror"
                    ]
                ),
                "selected_azimuth_deg": (
                    candidate[
                        "azimuth_deg"
                    ]
                ),
                "local_angle_violations": (
                    candidate[
                        "angle_violations"
                    ]
                ),
                "local_heavy_clashes": (
                    candidate[
                        "local_clashes"
                    ]
                ),
                "local_minimum_clearance_nm": (
                    candidate[
                        "minimum_clearance_nm"
                    ]
                ),
                "bond_1_nm": float(
                    candidate[
                        "bond_lengths"
                    ][0]
                ),
                "bond_2_nm": float(
                    candidate[
                        "bond_lengths"
                    ][1]
                ),
                "bond_3_nm": float(
                    candidate[
                        "bond_lengths"
                    ][2]
                ),
                "bond_4_nm": float(
                    candidate[
                        "bond_lengths"
                    ][3]
                ),
                "internal_angle_1_deg": float(
                    candidate[
                        "internal_angles"
                    ][0]
                ),
                "internal_angle_2_deg": float(
                    candidate[
                        "internal_angles"
                    ][1]
                ),
                "internal_angle_3_deg": float(
                    candidate[
                        "internal_angles"
                    ][2]
                ),
                "seed_junction_angle_1_deg": float(
                    candidate[
                        "endpoint_angles"
                    ][0]
                ),
                "seed_junction_angle_2_deg": float(
                    candidate[
                        "endpoint_angles"
                    ][1]
                ),
                "annulus_junction_angle_1_deg": float(
                    candidate[
                        "endpoint_angles"
                    ][2]
                ),
                "annulus_junction_angle_2_deg": float(
                    candidate[
                        "endpoint_angles"
                    ][3]
                ),
            }
        )

    write_rows(
        REFINED_CONFORMERS,
        refined_conformer_rows,
    )

    annulus_center_by_end = {}

    for end in (
        "LOWER",
        "UPPER",
    ):
        annulus_ids = [
            node_id
            for node_id, row
            in nodes.items()
            if (
                row[
                    "end"
                ]
                == end
                and row[
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
                    for node_id in annulus_ids
                ],
                dtype=float,
            ),
            axis=0,
        )

    heavy_ids_sorted = sorted(
        heavy_ids
    )

    heavy_index = {
        node_id: index
        for index, node_id
        in enumerate(
            heavy_ids_sorted
        )
    }

    all_heavy_positions = np.asarray(
        [
            positions[
                node_id
            ]
            for node_id
            in heavy_ids_sorted
        ],
        dtype=float,
    )

    base_directions = fibonacci_directions(
        H_DIRECTION_COUNT
    )

    hydrogen_candidate_pools: dict[
        str,
        list[dict[str, Any]]
    ] = {}

    print(
        "Generating hydrogen-orientation candidate pools..."
    )

    for hydrogen_index, hydrogen_id in enumerate(
        sorted(
            hydrogen_ids
        ),
        start=1,
    ):
        row = nodes[
            hydrogen_id
        ]

        heavy_id = row.get(
            "attached_to",
            "",
        )

        if (
            not heavy_id
            or heavy_id
            not in positions
        ):
            raise RuntimeError(
                f"{hydrogen_id}: invalid attached_to node."
            )

        heavy_neighbor_ids = [
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

        if len(heavy_neighbor_ids) != 2:
            raise RuntimeError(
                f"{hydrogen_id}: attached heavy atom has "
                f"{len(heavy_neighbor_ids)} heavy neighbors."
            )

        center = positions[
            heavy_id
        ]

        heavy_neighbor_positions = [
            positions[
                neighbor
            ]
            for neighbor
            in heavy_neighbor_ids
        ]

        preferred_direction = preferred_h_direction(
            row,
            center,
            heavy_neighbor_positions,
            tube_center,
            tube_axis,
            annulus_center_by_end,
        )

        opposite_bisector = normalized(
            -sum(
                (
                    normalized(
                        neighbor
                        - center
                    )
                    for neighbor
                    in heavy_neighbor_positions
                ),
                np.zeros(
                    3,
                    dtype=float,
                ),
            ),
            preferred_direction,
        )

        candidate_directions = np.vstack(
            (
                base_directions,
                preferred_direction[
                    None,
                    :,
                ],
                opposite_bisector[
                    None,
                    :,
                ],
            )
        )

        candidates = []

        bond_length = xh_target(
            nodes[
                heavy_id
            ][
                "element"
            ]
        )

        for direction in candidate_directions:
            direction = normalized(
                direction
            )

            point = (
                center
                + bond_length
                * direction
            )

            metrics = hydrogen_local_metrics(
                point,
                center,
                heavy_neighbor_positions,
                all_heavy_positions,
                heavy_index[
                    heavy_id
                ],
                preferred_direction,
            )

            candidates.append(
                {
                    "point": point,
                    "direction_x": float(
                        direction[0]
                    ),
                    "direction_y": float(
                        direction[1]
                    ),
                    "direction_z": float(
                        direction[2]
                    ),
                    **metrics,
                }
            )

        candidates.sort(
            key=hydrogen_local_score
        )

        hydrogen_candidate_pools[
            hydrogen_id
        ] = candidates[
            :H_POOL_SIZE
        ]

        if (
            hydrogen_index % 25 == 0
            or hydrogen_index == N_H
        ):
            best = hydrogen_candidate_pools[
                hydrogen_id
            ][0]

            print(
                f"  H pools {hydrogen_index}/{N_H}; "
                f"latest best angle violations/"
                f"heavy clashes="
                f"{best['angle_violations']}/"
                f"{best['heavy_clashes']}"
            )

    selected_h_index = {
        hydrogen_id: 0
        for hydrogen_id in hydrogen_ids
    }

    print(
        "Running global hydrogen-orientation coordinate descent..."
    )

    sorted_hydrogen_ids = sorted(
        hydrogen_ids
    )

    for sweep in range(
        H_GLOBAL_SWEEPS
    ):
        changed = 0

        ordered_hydrogens = (
            sorted_hydrogen_ids
            if sweep % 2 == 0
            else list(
                reversed(
                    sorted_hydrogen_ids
                )
            )
        )

        for hydrogen_id in ordered_hydrogens:
            other_points = [
                hydrogen_candidate_pools[
                    other_id
                ][
                    selected_h_index[
                        other_id
                    ]
                ][
                    "point"
                ]
                for other_id
                in sorted_hydrogen_ids
                if other_id
                != hydrogen_id
            ]

            if other_points:
                other_array = np.asarray(
                    other_points,
                    dtype=float,
                )
            else:
                other_array = np.empty(
                    (
                        0,
                        3,
                    ),
                    dtype=float,
                )

            pool = hydrogen_candidate_pools[
                hydrogen_id
            ]

            scores = [
                hydrogen_global_score(
                    candidate,
                    other_array,
                )
                for candidate in pool
            ]

            best_index = min(
                range(
                    len(pool)
                ),
                key=lambda index: scores[
                    index
                ],
            )

            if (
                best_index
                != selected_h_index[
                    hydrogen_id
                ]
            ):
                selected_h_index[
                    hydrogen_id
                ] = best_index

                changed += 1

        print(
            f"  H sweep {sweep + 1}/"
            f"{H_GLOBAL_SWEEPS}: "
            f"changed={changed}"
        )

        if changed == 0:
            break

    refined_hydrogen_rows = []

    for hydrogen_id in sorted_hydrogen_ids:
        candidate = (
            hydrogen_candidate_pools[
                hydrogen_id
            ][
                selected_h_index[
                    hydrogen_id
                ]
            ]
        )

        positions[
            hydrogen_id
        ] = np.asarray(
            candidate[
                "point"
            ],
            dtype=float,
        )

        refined_hydrogen_rows.append(
            {
                "hydrogen_node": (
                    hydrogen_id
                ),
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
                "attached_heavy_node": nodes[
                    hydrogen_id
                ].get(
                    "attached_to",
                    "",
                ),
                "direction_x": (
                    candidate[
                        "direction_x"
                    ]
                ),
                "direction_y": (
                    candidate[
                        "direction_y"
                    ]
                ),
                "direction_z": (
                    candidate[
                        "direction_z"
                    ]
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
                "heavy_neighbor_angle_1_deg": float(
                    candidate[
                        "angles"
                    ][0]
                ),
                "heavy_neighbor_angle_2_deg": float(
                    candidate[
                        "angles"
                    ][1]
                ),
            }
        )

    write_rows(
        REFINED_HYDROGENS,
        refined_hydrogen_rows,
    )

    missing_positions = [
        node_id
        for node_id in nodes
        if node_id not in positions
    ]

    if (
        missing_positions
        or len(positions) != N_TOTAL
    ):
        raise RuntimeError(
            "Coordinate assignment incomplete: "
            + " | ".join(
                missing_positions[:20]
            )
        )

    coordinate_rows = []

    for node_id in sorted(
        positions
    ):
        point = positions[
            node_id
        ]

        if not np.all(
            np.isfinite(
                point
            )
        ):
            raise RuntimeError(
                f"Non-finite coordinate for {node_id}"
            )

        coordinate_rows.append(
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
                "coordinates_generated_by": (
                    "GLOBAL_DISCRETE_CONFORMER_AND_H_REFINEMENT"
                ),
                "energy_minimized": False,
                "MD_relaxed": False,
            }
        )

    write_rows(
        REFINED_COORDINATES,
        coordinate_rows,
    )

    bond_rows = []

    for row in edge_rows:
        first = row[
            "source_node"
        ]

        second = row[
            "target_node"
        ]

        length = float(
            np.linalg.norm(
                positions[
                    first
                ]
                - positions[
                    second
                ]
            )
        )

        if parse_bool(
            row[
                "heavy_atom_edge"
            ]
        ):
            target = BN

        else:
            heavy_node = (
                first
                if nodes[
                    first
                ][
                    "element"
                ]
                != "H"
                else second
            )

            target = xh_target(
                nodes[
                    heavy_node
                ][
                    "element"
                ]
            )

        bond_rows.append(
            {
                "edge_id": row[
                    "edge_id"
                ],
                "end": row[
                    "end"
                ],
                "edge_type": row[
                    "edge_type"
                ],
                "source_node": first,
                "target_node": second,
                "source_element": nodes[
                    first
                ][
                    "element"
                ],
                "target_element": nodes[
                    second
                ][
                    "element"
                ],
                "length_nm": length,
                "target_length_nm": (
                    target
                ),
                "deviation_nm": (
                    length
                    - target
                ),
                "absolute_deviation_nm": abs(
                    length
                    - target
                ),
            }
        )

    write_rows(
        BOND_LENGTHS,
        bond_rows,
    )

    bond_summary_rows = summarize_bond_types(
        bond_rows
    )

    write_rows(
        BOND_SUMMARY,
        bond_summary_rows,
    )

    heavy_bond_rows = [
        row
        for row in bond_rows
        if (
            nodes[
                row[
                    "source_node"
                ]
            ][
                "element"
            ]
            != "H"
            and nodes[
                row[
                    "target_node"
                ]
            ][
                "element"
            ]
            != "H"
        )
    ]

    bridge_bond_rows = [
        row
        for row in heavy_bond_rows
        if row[
            "edge_type"
        ]
        == "ALTERNATING_BN_TRIMER_BRIDGE"
    ]

    h_bond_rows = [
        row
        for row in bond_rows
        if row not in heavy_bond_rows
    ]

    maximum_BN_deviation = max(
        float(
            row[
                "absolute_deviation_nm"
            ]
        )
        for row in heavy_bond_rows
    )

    maximum_bridge_BN_deviation = max(
        float(
            row[
                "absolute_deviation_nm"
            ]
        )
        for row in bridge_bond_rows
    )

    maximum_XH_deviation = max(
        float(
            row[
                "absolute_deviation_nm"
            ]
        )
        for row in h_bond_rows
    )

    angle_values_by_type: dict[
        str,
        list[float]
    ] = defaultdict(list)

    critical_angles = []

    for center_id, row in nodes.items():
        if row[
            "element"
        ] == "H":
            continue

        neighbors = sorted(
            adjacency[
                center_id
            ]
        )

        if len(neighbors) != 3:
            raise RuntimeError(
                f"Heavy node {center_id} has "
                f"{len(neighbors)} neighbors."
            )

        for first_index in range(3):
            for second_index in range(
                first_index + 1,
                3,
            ):
                value = angle_degrees(
                    positions[
                        neighbors[
                            first_index
                        ]
                    ],
                    positions[
                        center_id
                    ],
                    positions[
                        neighbors[
                            second_index
                        ]
                    ],
                )

                angle_values_by_type[
                    row[
                        "node_type"
                    ]
                ].append(value)

                if row[
                    "node_type"
                ] in CRITICAL_NODE_TYPES:
                    critical_angles.append(
                        value
                    )

    angle_summary_rows = []

    for node_type in sorted(
        angle_values_by_type
    ):
        values = np.asarray(
            angle_values_by_type[
                node_type
            ],
            dtype=float,
        )

        angle_summary_rows.append(
            {
                "center_node_type": node_type,
                "angle_count": int(
                    values.size
                ),
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
        np.min(
            critical_array
        )
    )

    critical_mean = float(
        np.mean(
            critical_array
        )
    )

    critical_maximum = float(
        np.max(
            critical_array
        )
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

    ordered_ids = sorted(
        positions
    )

    index_by_id = {
        node_id: index
        for index, node_id
        in enumerate(
            ordered_ids
        )
    }

    coordinate_array = np.asarray(
        [
            positions[
                node_id
            ]
            for node_id in ordered_ids
        ],
        dtype=float,
    )

    bonded_index_pairs = {
        tuple(
            sorted(
                (
                    index_by_id[
                        first
                    ],
                    index_by_id[
                        second
                    ],
                )
            )
        )
        for first, second
        in edge_pairs
    }

    thresholds = {
        "HEAVY_HEAVY": (
            MIN_NONBONDED_HEAVY_HEAVY_NM
        ),
        "H_HEAVY": (
            MIN_NONBONDED_H_HEAVY_NM
        ),
        "H_H": (
            MIN_NONBONDED_H_H_NM
        ),
    }

    minimum_distances = {
        category: math.inf
        for category in thresholds
    }

    minimum_pairs = {
        category: ""
        for category in thresholds
    }

    clash_counts = {
        category: 0
        for category in thresholds
    }

    for first_index in range(
        len(
            ordered_ids
        )
    ):
        distances = np.linalg.norm(
            coordinate_array[
                first_index + 1:
            ]
            - coordinate_array[
                first_index
            ],
            axis=1,
        )

        first_id = ordered_ids[
            first_index
        ]

        first_is_h = (
            nodes[
                first_id
            ][
                "element"
            ]
            == "H"
        )

        for offset, distance in enumerate(
            distances
        ):
            second_index = (
                first_index
                + 1
                + offset
            )

            if (
                first_index,
                second_index,
            ) in bonded_index_pairs:
                continue

            second_id = ordered_ids[
                second_index
            ]

            second_is_h = (
                nodes[
                    second_id
                ][
                    "element"
                ]
                == "H"
            )

            if first_is_h and second_is_h:
                category = "H_H"

            elif first_is_h or second_is_h:
                category = "H_HEAVY"

            else:
                category = "HEAVY_HEAVY"

            value = float(
                distance
            )

            if value < minimum_distances[
                category
            ]:
                minimum_distances[
                    category
                ] = value

                minimum_pairs[
                    category
                ] = (
                    f"{first_id} | {second_id}"
                )

            if value < thresholds[
                category
            ]:
                clash_counts[
                    category
                ] += 1

    contact_rows = [
        {
            "category": category,
            "minimum_distance_nm": (
                minimum_distances[
                    category
                ]
            ),
            "minimum_pair": (
                minimum_pairs[
                    category
                ]
            ),
            "threshold_nm": (
                thresholds[
                    category
                ]
            ),
            "pairs_below_threshold": (
                clash_counts[
                    category
                ]
            ),
        }
        for category in (
            "HEAVY_HEAVY",
            "H_HEAVY",
            "H_H",
        )
    ]

    write_rows(
        CONTACT_SUMMARY,
        contact_rows,
    )

    end_rows = []

    for end in (
        "LOWER",
        "UPPER",
    ):
        outward = (
            -tube_axis
            if end == "LOWER"
            else tube_axis
        )

        annulus_ids = [
            node_id
            for node_id, row in nodes.items()
            if (
                row[
                    "end"
                ]
                == end
                and row[
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
                    for node_id in annulus_ids
                ],
                dtype=float,
            ),
            axis=0,
        )

        inner_h_ids = [
            node_id
            for node_id, row in nodes.items()
            if (
                row[
                    "end"
                ]
                == end
                and row[
                    "node_type"
                ]
                == "ANNULUS_INNER_PASSIVANT_H"
            )
        ]

        outer_annulus_ids = [
            node_id
            for node_id, row in nodes.items()
            if (
                row[
                    "end"
                ]
                == end
                and row[
                    "node_type"
                ]
                == "ANNULUS_OUTER_BOUNDARY"
            )
        ]

        inner_radii = []

        for node_id in inner_h_ids:
            displacement = (
                positions[
                    node_id
                ]
                - annulus_center
            )

            displacement -= (
                np.dot(
                    displacement,
                    outward,
                )
                * outward
            )

            inner_radii.append(
                float(
                    np.linalg.norm(
                        displacement
                    )
                )
            )

        outer_radii = []

        for node_id in outer_annulus_ids:
            displacement = (
                positions[
                    node_id
                ]
                - annulus_center
            )

            displacement -= (
                np.dot(
                    displacement,
                    outward,
                )
                * outward
            )

            outer_radii.append(
                float(
                    np.linalg.norm(
                        displacement
                    )
                )
            )

        end_conformer_rows = [
            row
            for row in refined_conformer_rows
            if row[
                "end"
            ]
            == end
        ]

        nuclear_aperture = (
            2.0
            * min(
                inner_radii
            )
        )

        outer_radius = float(
            np.mean(
                outer_radii
            )
        )

        end_rows.append(
            {
                "end": end,
                "bridge_conformers": (
                    len(
                        end_conformer_rows
                    )
                ),
                "conformers_with_local_angle_violations": sum(
                    int(
                        row[
                            "local_angle_violations"
                        ]
                    )
                    > 0
                    for row in end_conformer_rows
                ),
                "conformers_with_local_heavy_clashes": sum(
                    int(
                        row[
                            "local_heavy_clashes"
                        ]
                    )
                    > 0
                    for row in end_conformer_rows
                ),
                "minimum_local_bridge_clearance_nm": min(
                    float(
                        row[
                            "local_minimum_clearance_nm"
                        ]
                    )
                    for row in end_conformer_rows
                ),
                "maximum_library_distance_error_nm": max(
                    float(
                        row[
                            "library_distance_error_nm"
                        ]
                    )
                    for row in end_conformer_rows
                ),
                "nuclear_aperture_diameter_nm": (
                    nuclear_aperture
                ),
                "nuclear_aperture_relative_error": abs(
                    nuclear_aperture
                    - target_aperture_diameter_nm
                ) / target_aperture_diameter_nm,
                "outer_annulus_radius_mean_nm": (
                    outer_radius
                ),
                "outer_annulus_radius_relative_error": abs(
                    outer_radius
                    - target_outer_radius_nm
                ) / target_outer_radius_nm,
            }
        )

    write_rows(
        END_SUMMARY,
        end_rows,
    )

    lower = next(
        row
        for row in end_rows
        if row[
            "end"
        ]
        == "LOWER"
    )

    upper = next(
        row
        for row in end_rows
        if row[
            "end"
        ]
        == "UPPER"
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
        abs(
            float(
                lower[
                    "minimum_local_bridge_clearance_nm"
                ]
            )
            - float(
                upper[
                    "minimum_local_bridge_clearance_nm"
                ]
            )
        ),
    )

    parent_coordinates_unchanged = all(
        np.array_equal(
            positions[
                node_id
            ],
            source_coordinates[
                node_id
            ],
        )
        for node_id in parent_ids
    )

    fixed_heavy_coordinates_unchanged = all(
        np.array_equal(
            positions[
                node_id
            ],
            source_coordinates[
                node_id
            ],
        )
        for node_id in fixed_heavy_ids
    )

    maximum_library_error = max(
        float(
            row[
                "library_distance_error_nm"
            ]
        )
        for row in refined_conformer_rows
    )

    heavy_geometry_pass = (
        maximum_BN_deviation
        <= MAX_BN_DEVIATION_NM
        and maximum_bridge_BN_deviation
        <= MAX_BN_DEVIATION_NM
        and critical_minimum
        >= MIN_CRITICAL_ANGLE_DEG
        and critical_maximum
        <= MAX_CRITICAL_ANGLE_DEG
        and critical_rms
        <= MAX_CRITICAL_RMS_DEVIATION_DEG
        and clash_counts[
            "HEAVY_HEAVY"
        ]
        == 0
    )

    hydrogen_geometry_pass = (
        maximum_XH_deviation
        <= MAX_XH_DEVIATION_NM
        and clash_counts[
            "H_HEAVY"
        ]
        == 0
        and clash_counts[
            "H_H"
        ]
        == 0
    )

    gates = {
        "Gate3A_parent_audit_is_accepted": (
            parent_summary.get(
                "decision"
            )
            == EXPECTED_PARENT_DECISION
        ),
        "Gate3I_trimer_bridge_graph_is_accepted": (
            graph_summary.get(
                "decision"
            )
            == EXPECTED_GRAPH_DECISION
        ),
        "Gate3I_has_no_failed_gates": (
            len(
                failed_graph_gates
            )
            == 0
        ),
        "Gate3J_requires_conformer_refinement": (
            source_summary.get(
                "decision"
            )
            == EXPECTED_SOURCE_DECISION
        ),
        "all_2256_nodes_received_finite_coordinates": (
            len(
                positions
            )
            == N_TOTAL
            and all(
                np.all(
                    np.isfinite(
                        point
                    )
                )
                for point in positions.values()
            )
        ),
        "parent_coordinates_are_unchanged": (
            parent_coordinates_unchanged
        ),
        "seed_and_annulus_coordinates_are_unchanged": (
            fixed_heavy_coordinates_unchanged
        ),
        "30_bridge_conformers_were_globally_refined": (
            len(
                refined_conformer_rows
            )
            == N_BRIDGES
        ),
        "174_H_orientations_were_globally_refined": (
            len(
                refined_hydrogen_rows
            )
            == N_H
        ),
        "library_endpoint_errors_are_within_0p0005nm": (
            maximum_library_error
            <= MAX_LIBRARY_DISTANCE_ERROR_NM
        ),
        "all_BN_bonds_are_within_0p003nm": (
            maximum_BN_deviation
            <= MAX_BN_DEVIATION_NM
        ),
        "all_bridge_BN_bonds_are_within_0p003nm": (
            maximum_bridge_BN_deviation
            <= MAX_BN_DEVIATION_NM
        ),
        "all_XH_bonds_are_within_0p002nm": (
            maximum_XH_deviation
            <= MAX_XH_DEVIATION_NM
        ),
        "critical_angle_minimum_is_at_least70deg": (
            critical_minimum
            >= MIN_CRITICAL_ANGLE_DEG
        ),
        "critical_angle_maximum_is_at_most175deg": (
            critical_maximum
            <= MAX_CRITICAL_ANGLE_DEG
        ),
        "critical_angle_RMS_deviation_is_at_most30deg": (
            critical_rms
            <= MAX_CRITICAL_RMS_DEVIATION_DEG
        ),
        "no_nonbonded_heavy_heavy_clashes": (
            clash_counts[
                "HEAVY_HEAVY"
            ]
            == 0
        ),
        "no_nonbonded_H_heavy_clashes": (
            clash_counts[
                "H_HEAVY"
            ]
            == 0
        ),
        "no_nonbonded_H_H_clashes": (
            clash_counts[
                "H_H"
            ]
            == 0
        ),
        "aperture_errors_are_within10percent": all(
            float(
                row[
                    "nuclear_aperture_relative_error"
                ]
            )
            <= MAX_APERTURE_ERROR
            for row in end_rows
        ),
        "outer_radius_errors_are_within15percent": all(
            float(
                row[
                    "outer_annulus_radius_relative_error"
                ]
            )
            <= MAX_OUTER_RADIUS_ERROR
            for row in end_rows
        ),
        "lower_upper_asymmetry_is_within0p010nm": (
            end_asymmetry
            <= MAX_END_ASYMMETRY_NM
        ),
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

    if accepted:
        decision = (
            PASS_DECISION
        )

        required_next_step = (
            "AUDIT_R2_ALTERNATING_BN_TRIMER_BRIDGE_"
            "CHEMICAL_REALIZABILITY_AND_PARAMETERIZATION_SCOPE"
        )

    elif not heavy_geometry_pass:
        decision = (
            HEAVY_REVIEW_DECISION
        )

        required_next_step = (
            "EVALUATE_R2_LONGER_OR_TOPOLOGICALLY_REVISED_"
            "BRIDGE_ARCHITECTURE"
        )

    else:
        decision = (
            H_REVIEW_DECISION
        )

        required_next_step = (
            "RUN_R2_SECONDARY_H_ORIENTATION_AND_"
            "PROTONATION_PATTERN_REFINEMENT"
        )

    summary = {
        "decision": decision,
        "coordinate_nodes": (
            len(
                positions
            )
        ),
        "bridge_conformers_refined": (
            len(
                refined_conformer_rows
            )
        ),
        "H_orientations_refined": (
            len(
                refined_hydrogen_rows
            )
        ),
        "conformer_library_size": int(
            library[
                "distance"
            ].size
        ),
        "maximum_library_distance_error_nm": (
            maximum_library_error
        ),
        "maximum_BN_bond_deviation_nm": (
            maximum_BN_deviation
        ),
        "maximum_bridge_BN_bond_deviation_nm": (
            maximum_bridge_BN_deviation
        ),
        "maximum_XH_bond_deviation_nm": (
            maximum_XH_deviation
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
            minimum_distances[
                "HEAVY_HEAVY"
            ]
        ),
        "minimum_nonbonded_H_heavy_nm": (
            minimum_distances[
                "H_HEAVY"
            ]
        ),
        "minimum_nonbonded_H_H_nm": (
            minimum_distances[
                "H_H"
            ]
        ),
        "heavy_heavy_clash_count": (
            clash_counts[
                "HEAVY_HEAVY"
            ]
        ),
        "H_heavy_clash_count": (
            clash_counts[
                "H_HEAVY"
            ]
        ),
        "H_H_clash_count": (
            clash_counts[
                "H_H"
            ]
        ),
        "lower_nuclear_aperture_diameter_nm": (
            lower[
                "nuclear_aperture_diameter_nm"
            ]
        ),
        "upper_nuclear_aperture_diameter_nm": (
            upper[
                "nuclear_aperture_diameter_nm"
            ]
        ),
        "target_aperture_diameter_nm": (
            target_aperture_diameter_nm
        ),
        "maximum_lower_upper_asymmetry_nm": (
            end_asymmetry
        ),
        "heavy_geometry_pass": (
            heavy_geometry_pass
        ),
        "hydrogen_geometry_pass": (
            hydrogen_geometry_pass
        ),
        "candidate_is_final_chemistry": False,
        "molecular_topology_generation_authorized": False,
        "formal_charge_assignment_authorized": False,
        "force_field_parameterization_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "mobile_MD_authorized": False,
        "multitemperature_MD_authorized": False,
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
                "end_summaries": (
                    end_rows
                ),
                "gates": gates,
                "refinement_parameters": {
                    "base_conformers_per_bridge": (
                        BASE_CONFORMERS_PER_BRIDGE
                    ),
                    "bridge_azimuth_step_deg": (
                        AZIMUTH_STEP
                    ),
                    "bridge_pool_size": (
                        BRIDGE_POOL_SIZE
                    ),
                    "bridge_global_sweeps": (
                        BRIDGE_GLOBAL_SWEEPS
                    ),
                    "H_direction_count": (
                        H_DIRECTION_COUNT
                    ),
                    "H_pool_size": (
                        H_POOL_SIZE
                    ),
                    "H_global_sweeps": (
                        H_GLOBAL_SWEEPS
                    ),
                },
                "limitations": [
                    (
                        "This is a deterministic discrete geometric "
                        "refinement, not an energy minimization."
                    ),
                    (
                        "The graph, parent, seed and annulus coordinates "
                        "were held fixed."
                    ),
                    (
                        "B-H and N-H distances remain provisional "
                        "geometry targets."
                    ),
                    (
                        "Passing this gate would establish geometric "
                        "consistency only, not energetic stability."
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
                "Gate3A_parent_summary"
            ),
            "file": relative(
                PARENT_SUMMARY
            ),
            "sha256": sha256(
                PARENT_SUMMARY
            ),
        },
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
                "Gate3J_source_coordinates"
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
                "Gate3J_source_embedding_summary"
            ),
            "file": relative(
                SOURCE_EMBEDDING_SUMMARY
            ),
            "sha256": sha256(
                SOURCE_EMBEDDING_SUMMARY
            ),
        },
    ]

    write_rows(
        MANIFEST,
        manifest_rows,
    )

    xyz_lines = [
        str(
            len(
                ordered_ids
            )
        ),
        (
            "R2 refined BN trimer-bridge embedding; "
            "not energy minimized"
        ),
    ]

    for node_id in ordered_ids:
        coordinates_angstrom = (
            positions[
                node_id
            ]
            * 10.0
        )

        xyz_lines.append(
            f"{nodes[node_id]['element']:2s} "
            f"{coordinates_angstrom[0]: .8f} "
            f"{coordinates_angstrom[1]: .8f} "
            f"{coordinates_angstrom[2]: .8f}"
        )

    XYZ_OUT.write_text(
        "\n".join(
            xyz_lines
        )
        + "\n",
        encoding="utf-8",
    )

    pdb_lines = [
        "REMARK R2 REFINED BN TRIMER-BRIDGE STATIC EMBEDDING",
        "REMARK NOT ENERGY MINIMIZED; NO FORCE-FIELD TOPOLOGY",
    ]

    for serial, node_id in enumerate(
        ordered_ids,
        start=1,
    ):
        row = nodes[
            node_id
        ]

        coordinates_angstrom = (
            positions[
                node_id
            ]
            * 10.0
        )

        element = row[
            "element"
        ]

        node_type = row[
            "node_type"
        ]

        if node_type == "PARENT_HBN":
            residue = "HBN"

        elif node_type == "HEXAGONAL_EDGE_COMPLETION_SEED":
            residue = "SED"

        elif node_type == "ALTERNATING_BN_TRIMER_BRIDGE":
            residue = "BRG"

        elif node_type.startswith(
            "ANNULUS"
        ):
            residue = "ANN"

        else:
            residue = "PAS"

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
            f"{coordinates_angstrom[0]:8.3f}"
            f"{coordinates_angstrom[1]:8.3f}"
            f"{coordinates_angstrom[2]:8.3f}"
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
        f"""# R2 Trimer-Bridge Conformer and H Refinement

## Scope

This gate globally refines the discrete conformers of all 30 BN trimer
bridges and the orientations of all 174 H passivants.

The validated graph and the parent, seed and annulus coordinates were
not changed.

No topology, formal charges, force-field parameters, minimization, MD
or QM calculation was generated.

## Refined inventory

- Coordinate nodes:
  **{len(positions)}**
- Bridge conformers:
  **{len(refined_conformer_rows)}**
- H orientations:
  **{len(refined_hydrogen_rows)}**
- Parent coordinates unchanged:
  **{parent_coordinates_unchanged}**
- Seed and annulus coordinates unchanged:
  **{fixed_heavy_coordinates_unchanged}**

## Bond geometry

- Maximum B-N deviation:
  **{maximum_BN_deviation:.6f} nm**
- Maximum bridge B-N deviation:
  **{maximum_bridge_BN_deviation:.6f} nm**
- Maximum X-H deviation:
  **{maximum_XH_deviation:.6f} nm**
- Maximum library endpoint mismatch:
  **{maximum_library_error:.6f} nm**

## Critical angles

- Minimum/mean/maximum:
  **{critical_minimum:.3f}/
  {critical_mean:.3f}/
  {critical_maximum:.3f} degrees**
- RMS deviation from 120 degrees:
  **{critical_rms:.3f} degrees**

## Nonbonded contacts

- Heavy-heavy minimum/clashes:
  **{minimum_distances['HEAVY_HEAVY']:.6f}/
  {clash_counts['HEAVY_HEAVY']}**
- H-heavy minimum/clashes:
  **{minimum_distances['H_HEAVY']:.6f}/
  {clash_counts['H_HEAVY']}**
- H-H minimum/clashes:
  **{minimum_distances['H_H']:.6f}/
  {clash_counts['H_H']}**

## Aperture

- Target diameter:
  **{target_aperture_diameter_nm:.6f} nm**
- Lower/upper refined diameter:
  **{float(lower['nuclear_aperture_diameter_nm']):.6f}/
  {float(upper['nuclear_aperture_diameter_nm']):.6f} nm**

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Heavy geometry pass:
  **{heavy_geometry_pass}**
- Hydrogen geometry pass:
  **{hydrogen_geometry_pass}**
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
        "Day024 R2 trimer-bridge conformer and H "
        "refinement completed."
    )

    print(
        "Coordinate nodes / bridge conformers / H orientations: "
        f"{len(positions)}/"
        f"{len(refined_conformer_rows)}/"
        f"{len(refined_hydrogen_rows)}"
    )

    print(
        "Parent / seed-annulus coordinates unchanged: "
        f"{parent_coordinates_unchanged}/"
        f"{fixed_heavy_coordinates_unchanged}"
    )

    print(
        "Maximum library / BN / bridge-BN / XH deviations: "
        f"{maximum_library_error:.6f}/"
        f"{maximum_BN_deviation:.6f}/"
        f"{maximum_bridge_BN_deviation:.6f}/"
        f"{maximum_XH_deviation:.6f} nm"
    )

    print(
        "Critical angles min/mean/max/RMSdev120: "
        f"{critical_minimum:.3f}/"
        f"{critical_mean:.3f}/"
        f"{critical_maximum:.3f}/"
        f"{critical_rms:.3f} deg"
    )

    print(
        "Minimum nonbonded heavy-heavy / H-heavy / H-H: "
        f"{minimum_distances['HEAVY_HEAVY']:.6f}/"
        f"{minimum_distances['H_HEAVY']:.6f}/"
        f"{minimum_distances['H_H']:.6f} nm"
    )

    print(
        "Clash counts heavy-heavy / H-heavy / H-H: "
        f"{clash_counts['HEAVY_HEAVY']}/"
        f"{clash_counts['H_HEAVY']}/"
        f"{clash_counts['H_H']}"
    )

    for row in end_rows:
        print(
            f"{row['end']} local-angle-fail conformers / "
            "local-clash conformers / min-clearance / "
            "aperture/error / outer-radius/error: "
            f"{row['conformers_with_local_angle_violations']}/"
            f"{row['conformers_with_local_heavy_clashes']}/"
            f"{float(row['minimum_local_bridge_clearance_nm']):.6f}/"
            f"{float(row['nuclear_aperture_diameter_nm']):.6f}/"
            f"{float(row['nuclear_aperture_relative_error']):.6f}/"
            f"{float(row['outer_annulus_radius_mean_nm']):.6f}/"
            f"{float(row['outer_annulus_radius_relative_error']):.6f}"
        )

    print(
        "Maximum lower-upper asymmetry: "
        f"{end_asymmetry:.6f} nm"
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
        "Heavy geometry pass: "
        f"{heavy_geometry_pass}"
    )

    print(
        "Hydrogen geometry pass: "
        f"{hydrogen_geometry_pass}"
    )

    print(
        "Candidate is final chemistry: NO"
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
        REFINED_CONFORMERS,
        REFINED_HYDROGENS,
        BOND_LENGTHS,
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


if __name__ == "__main__":
    main()
