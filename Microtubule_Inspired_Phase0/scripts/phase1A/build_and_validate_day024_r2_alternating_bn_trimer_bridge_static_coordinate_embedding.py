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
G3G = BASE / "08_r2_partial_attachment_annulus_static_coordinate_embedding"
G3G1 = BASE / "09_r2_direct_junction_geometric_lower_bound"
G3I = BASE / "11_r2_alternating_bn_trimer_bridge_graph"

OUT = (
    BASE
    / "12_r2_alternating_bn_trimer_bridge_static_coordinate_embedding"
)

PARENT_SUMMARY = (
    G3A
    / "r2_parent_rim_chemical_audit_summary.csv"
)

SOURCE_COORDS = (
    G3G
    / "r2_partial_attachment_static_coordinates.csv"
)

DIRECT_SUMMARY = (
    G3G1
    / "r2_direct_junction_geometric_lower_bound_summary.csv"
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

COORDS_OUT = (
    OUT
    / "r2_trimer_bridge_static_coordinates.csv"
)

CONFORMERS_OUT = (
    OUT
    / "r2_trimer_bridge_selected_conformers.csv"
)

BONDS_OUT = (
    OUT
    / "r2_trimer_bridge_static_bond_lengths.csv"
)

BOND_SUMMARY_OUT = (
    OUT
    / "r2_trimer_bridge_static_bond_type_summary.csv"
)

ANGLE_SUMMARY_OUT = (
    OUT
    / "r2_trimer_bridge_static_angle_summary.csv"
)

CONTACTS_OUT = (
    OUT
    / "r2_trimer_bridge_static_nonbonded_contact_summary.csv"
)

END_SUMMARY_OUT = (
    OUT
    / "r2_trimer_bridge_static_embedding_end_summary.csv"
)

SUMMARY_OUT = (
    OUT
    / "r2_trimer_bridge_static_embedding_summary.csv"
)

GATES_OUT = (
    OUT
    / "r2_trimer_bridge_static_embedding_gates.csv"
)

JSON_OUT = (
    OUT
    / "r2_trimer_bridge_static_embedding.json"
)

MANIFEST_OUT = (
    OUT
    / "r2_trimer_bridge_static_embedding_source_manifest.csv"
)

XYZ_OUT = (
    OUT
    / "r2_trimer_bridge_static_embedding.xyz"
)

PDB_OUT = (
    OUT
    / "r2_trimer_bridge_static_embedding.pdb"
)

REPORT_OUT = (
    OUT
    / "R2_ALTERNATING_BN_TRIMER_BRIDGE_STATIC_EMBEDDING_DAY024.md"
)

EXPECTED_PARENT_DECISION = (
    "R2_PARENT_RIM_AND_CHEMICAL_CONSTRAINTS_AUDITED"
)

EXPECTED_DIRECT_DECISION = (
    "R2_PARTIAL_ATTACHMENT_DIRECT_BN_JUNCTION_"
    "GEOMETRIC_LOWER_BOUND_FAILED"
)

EXPECTED_GRAPH_DECISION = (
    "R2_ALTERNATING_BN_TRIMER_BRIDGE_GRAPH_VALIDATED"
)

PASS_DECISION = (
    "R2_ALTERNATING_BN_TRIMER_BRIDGE_STATIC_"
    "COORDINATE_EMBEDDING_VALIDATED"
)

REVIEW_DECISION = (
    "R2_ALTERNATING_BN_TRIMER_BRIDGE_STATIC_"
    "COORDINATE_EMBEDDING_REQUIRES_CONFORMER_REFINEMENT"
)

N_PARENT = 1680
N_HEAVY = 2082
N_H = 174
N_TOTAL = 2256
N_BRIDGES = 30

BN = 0.144973
BH = 0.119
NH = 0.101

ANGLE_MIN = 105.0
ANGLE_MAX = 135.0
ANGLE_STEP = 1.0
TORSION_STEP = 10.0
AZIMUTH_STEP = 10.0

MAX_BN_DEV = 0.003
MAX_XH_DEV = 0.002
MAX_LIBRARY_DISTANCE_ERROR = 0.0005

MIN_CRITICAL_ANGLE = 70.0
MAX_CRITICAL_ANGLE = 175.0
MAX_CRITICAL_RMS = 30.0

MIN_HH = 0.120
MIN_HX = 0.070
MIN_HYD_HYD = 0.060

MAX_CENTER_OFFSET = 0.050
MAX_APERTURE_ERROR = 0.10
MAX_OUTER_RADIUS_ERROR = 0.15
MAX_END_ASYMMETRY = 0.010

CRITICAL_TYPES = {
    "HEXAGONAL_EDGE_COMPLETION_SEED",
    "ALTERNATING_BN_TRIMER_BRIDGE",
    "ANNULUS_OUTER_BOUNDARY",
    "ANNULUS_INNER_BOUNDARY",
}


def rel(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
    )


def require(path: Path) -> None:
    if (
        not path.is_file()
        or path.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing or empty required file: {path}"
        )


def digest(path: Path) -> str:
    value = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            value.update(chunk)

    return value.hexdigest()


def read_rows(
    path: Path,
) -> list[dict[str, str]]:
    require(path)

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

        writer.writerows(
            {
                field: row.get(field, "")
                for field in fields
            }
            for row in rows
        )


def as_bool(value: Any) -> bool:
    return str(value).strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


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
            f"Non-finite numeric field {key!r}"
        )

    return value


def as_int(
    row: dict[str, str],
    key: str,
) -> int:
    try:
        return int(
            float(row[key])
        )
    except (
        KeyError,
        TypeError,
        ValueError,
    ) as exc:
        raise RuntimeError(
            f"Could not parse integer field {key!r}"
        ) from exc


def unit(
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
            "Could not normalize a zero-length vector."
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
) -> tuple[np.ndarray, np.ndarray]:
    center = np.mean(
        positions,
        axis=0,
    )

    centered = (
        positions
        - center
    )

    eigenvalues, eigenvectors = (
        np.linalg.eigh(
            centered.T
            @ centered
        )
    )

    axis = unit(
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


def local_xy(
    row: dict[str, str],
) -> np.ndarray:
    lattice_x = as_int(
        row,
        "lattice_x",
    )

    lattice_y = as_int(
        row,
        "lattice_y",
    )

    return np.asarray(
        [
            lattice_x
            * BN
            / 2.0,
            lattice_y
            * math.sqrt(3.0)
            * BN
            / 2.0,
        ],
        dtype=float,
    )


def fit_orthogonal(
    source: np.ndarray,
    target: np.ndarray,
    determinant_target: int,
) -> tuple[
    np.ndarray,
    np.ndarray,
    float,
]:
    source_center = np.mean(
        source,
        axis=0,
    )

    target_center = np.mean(
        target,
        axis=0,
    )

    covariance = (
        (
            source
            - source_center
        ).T
        @ (
            target
            - target_center
        )
    )

    left, _, right_t = (
        np.linalg.svd(
            covariance
        )
    )

    rotation = (
        left
        @ right_t
    )

    if int(
        round(
            np.linalg.det(
                rotation
            )
        )
    ) != determinant_target:
        correction = np.eye(
            2,
            dtype=float,
        )

        correction[
            -1,
            -1,
        ] = -1.0

        rotation = (
            left
            @ correction
            @ right_t
        )

    translation = (
        target_center
        - source_center
        @ rotation
    )

    fitted = (
        source
        @ rotation
        + translation
    )

    rms = float(
        np.sqrt(
            np.mean(
                np.sum(
                    (
                        fitted
                        - target
                    )
                    ** 2,
                    axis=1,
                )
            )
        )
    )

    return (
        rotation,
        translation,
        rms,
    )


def rotation_about_axis(
    axis: np.ndarray,
    angle: float,
) -> np.ndarray:
    x_value, y_value, z_value = unit(
        axis
    )

    cosine = math.cos(angle)
    sine = math.sin(angle)
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
    source_unit = unit(source)
    target_unit = unit(target)

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
            unit(
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
    angle: float,
    dihedral: float,
) -> np.ndarray:
    previous = unit(
        third
        - second
    )

    normal = np.cross(
        second
        - first,
        previous,
    )

    if np.linalg.norm(normal) <= 1.0e-12:
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

    normal = unit(normal)

    in_plane = unit(
        np.cross(
            normal,
            previous,
        )
    )

    direction = (
        -math.cos(angle)
        * previous
        + math.sin(angle)
        * (
            math.cos(dihedral)
            * in_plane
            + math.sin(dihedral)
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
    angle = math.radians(
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
                -math.cos(angle),
                math.sin(angle),
                0.0,
            ],
            dtype=float,
        )
    )

    point_3 = place_next(
        point_0,
        point_1,
        point_2,
        angle,
        math.radians(
            phi1_deg
        ),
    )

    point_4 = place_next(
        point_1,
        point_2,
        point_3,
        angle,
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


def conformer_library() -> dict[str, np.ndarray]:
    distances = []
    angles = []
    phi1_values = []
    phi2_values = []

    for angle in np.arange(
        ANGLE_MIN,
        ANGLE_MAX
        + 0.5
        * ANGLE_STEP,
        ANGLE_STEP,
    ):
        for phi1 in np.arange(
            0.0,
            360.0,
            TORSION_STEP,
        ):
            for phi2 in np.arange(
                0.0,
                360.0,
                TORSION_STEP,
            ):
                chain = canonical_chain(
                    float(angle),
                    float(phi1),
                    float(phi2),
                )

                distances.append(
                    float(
                        np.linalg.norm(
                            chain[-1]
                            - chain[0]
                        )
                    )
                )

                angles.append(
                    float(angle)
                )

                phi1_values.append(
                    float(phi1)
                )

                phi2_values.append(
                    float(phi2)
                )

    return {
        "distance": np.asarray(
            distances,
            dtype=float,
        ),
        "angle": np.asarray(
            angles,
            dtype=float,
        ),
        "phi1": np.asarray(
            phi1_values,
            dtype=float,
        ),
        "phi2": np.asarray(
            phi2_values,
            dtype=float,
        ),
    }


def select_base_chain(
    distance: float,
    library: dict[str, np.ndarray],
) -> tuple[
    np.ndarray,
    dict[str, float],
]:
    errors = np.abs(
        library[
            "distance"
        ]
        - distance
    )

    score = (
        errors
        + 1.0e-7
        * np.abs(
            library[
                "angle"
            ]
            - 120.0
        )
    )

    index = int(
        np.argmin(score)
    )

    metadata = {
        "angle_deg": float(
            library[
                "angle"
            ][index]
        ),
        "phi1_deg": float(
            library[
                "phi1"
            ][index]
        ),
        "phi2_deg": float(
            library[
                "phi2"
            ][index]
        ),
        "library_distance_nm": float(
            library[
                "distance"
            ][index]
        ),
        "library_error_nm": float(
            errors[index]
        ),
    }

    return (
        canonical_chain(
            metadata[
                "angle_deg"
            ],
            metadata[
                "phi1_deg"
            ],
            metadata[
                "phi2_deg"
            ],
        ),
        metadata,
    )


def map_chain(
    chain: np.ndarray,
    start: np.ndarray,
    end: np.ndarray,
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

    target = (
        end
        - start
    )

    alignment = rotation_from_vectors(
        working[-1]
        - working[0],
        target,
    )

    mapped = (
        working
        - working[0]
    ) @ alignment.T

    mapped = (
        mapped
        @ rotation_about_axis(
            unit(target),
            math.radians(
                azimuth_deg
            ),
        ).T
    )

    correction = (
        target
        - mapped[-1]
    )

    for index in range(
        1,
        5,
    ):
        mapped[index] += (
            index
            / 4.0
        ) * correction

    mapped += start

    if np.linalg.norm(
        mapped[-1]
        - end
    ) > 1.0e-10:
        raise RuntimeError(
            "Endpoint correction failed."
        )

    return mapped


def clearance_score(
    internal: np.ndarray,
    fixed: np.ndarray,
    radial: np.ndarray,
    midpoint: np.ndarray,
) -> tuple[
    int,
    float,
    float,
    float,
]:
    distances = np.linalg.norm(
        internal[:, None, :]
        - fixed[None, :, :],
        axis=2,
    )

    internal_13 = float(
        np.linalg.norm(
            internal[0]
            - internal[2]
        )
    )

    minimum = min(
        float(
            np.min(distances)
        ),
        internal_13,
    )

    clashes = (
        int(
            np.sum(
                distances
                < MIN_HH
            )
        )
        + int(
            internal_13
            < MIN_HH
        )
    )

    penalty = float(
        np.sum(
            np.maximum(
                0.140
                - distances,
                0.0,
            )
            ** 2
        )
        + max(
            0.140
            - internal_13,
            0.0,
        )
        ** 2
    )

    bulge = float(
        np.dot(
            np.mean(
                internal,
                axis=0,
            )
            - midpoint,
            radial,
        )
    )

    return (
        clashes,
        penalty,
        -bulge,
        -minimum,
    )


def xh_length(
    element: str,
) -> float:
    if element == "B":
        return BH

    if element == "N":
        return NH

    raise RuntimeError(
        f"Unexpected H-attached element: {element}"
    )


def angle_deg(
    first: np.ndarray,
    center: np.ndarray,
    second: np.ndarray,
) -> float:
    cosine = float(
        np.clip(
            np.dot(
                unit(
                    first
                    - center
                ),
                unit(
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


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    required_files = [
        PARENT_SUMMARY,
        SOURCE_COORDS,
        DIRECT_SUMMARY,
        GRAPH_NODES,
        GRAPH_EDGES,
        BRIDGE_PATHS,
        GRAPH_SUMMARY,
        GRAPH_GATES,
    ]

    for path in required_files:
        require(path)

    parent_summary = read_one(
        PARENT_SUMMARY
    )

    direct_summary = read_one(
        DIRECT_SUMMARY
    )

    graph_summary = read_one(
        GRAPH_SUMMARY
    )

    gate_rows = read_rows(
        GRAPH_GATES
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

    source_rows = read_rows(
        SOURCE_COORDS
    )

    if parent_summary.get(
        "decision"
    ) != EXPECTED_PARENT_DECISION:
        raise RuntimeError(
            "Gate 3A is not accepted."
        )

    if direct_summary.get(
        "decision"
    ) != EXPECTED_DIRECT_DECISION:
        raise RuntimeError(
            "Gate 3G.1 does not contain the expected "
            "direct-junction rejection."
        )

    if graph_summary.get(
        "decision"
    ) != EXPECTED_GRAPH_DECISION:
        raise RuntimeError(
            "Gate 3I is not accepted."
        )

    failed_upstream = [
        row.get(
            "gate",
            "",
        )
        for row in gate_rows
        if not as_bool(
            row.get(
                "pass",
                "false",
            )
        )
    ]

    if failed_upstream:
        raise RuntimeError(
            "Gate 3I contains failed gates: "
            + " | ".join(
                failed_upstream
            )
        )

    if (
        len(node_rows)
        != N_TOTAL
        or len(path_rows)
        != N_BRIDGES
    ):
        raise RuntimeError(
            "Unexpected node/path counts: "
            f"{len(node_rows)}/"
            f"{len(path_rows)}"
        )

    nodes = {
        row[
            "node_id"
        ]: row
        for row in node_rows
    }

    if len(nodes) != len(
        node_rows
    ):
        raise RuntimeError(
            "Duplicate graph-node identifiers."
        )

    adjacency = {
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
                "Graph edge references a missing node."
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

    source_coords = {
        row[
            "node_id"
        ]: np.asarray(
            [
                as_float(
                    row,
                    "x_nm",
                ),
                as_float(
                    row,
                    "y_nm",
                ),
                as_float(
                    row,
                    "z_nm",
                ),
            ],
            dtype=float,
        )
        for row in source_rows
    }

    parent_ids = sorted(
        [
            node_id
            for node_id, row
            in nodes.items()
            if row[
                "node_type"
            ]
            == "PARENT_HBN"
        ],
        key=lambda node_id: int(
            node_id.split(
                ":"
            )[1]
        ),
    )

    seed_ids = [
        node_id
        for node_id, row
        in nodes.items()
        if row[
            "node_type"
        ]
        == "HEXAGONAL_EDGE_COMPLETION_SEED"
    ]

    if len(parent_ids) != N_PARENT:
        raise RuntimeError(
            "Unexpected parent-node count."
        )

    positions: dict[
        str,
        np.ndarray
    ] = {}

    for node_id in (
        parent_ids
        + seed_ids
    ):
        if node_id not in source_coords:
            raise RuntimeError(
                f"Missing source coordinate for {node_id}"
            )

        positions[node_id] = np.array(
            source_coords[
                node_id
            ],
            dtype=float,
            copy=True,
        )

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

    target_aperture = as_float(
        parent_summary,
        "target_aperture_diameter_nm",
    )

    target_outer_radius = as_float(
        parent_summary,
        "parent_rim_mean_radius_nm",
    )

    library = conformer_library()

    end_geometry: dict[
        str,
        dict[str, Any]
    ] = {}

    for end in (
        "LOWER",
        "UPPER",
    ):
        end_paths = sorted(
            [
                row
                for row in path_rows
                if row[
                    "end"
                ]
                == end
            ],
            key=lambda row: as_int(
                row,
                "bridge_index",
            ),
        )

        if len(end_paths) != 15:
            raise RuntimeError(
                f"{end}: unexpected bridge-path count."
            )

        gaps = {
            round(
                as_float(
                    row,
                    "selected_gap_nm",
                ),
                12,
            )
            for row in end_paths
        }

        chiralities = {
            as_int(
                row,
                "mapping_chirality",
            )
            for row in end_paths
        }

        if (
            len(gaps) != 1
            or len(chiralities) != 1
        ):
            raise RuntimeError(
                f"{end}: inconsistent gap or chirality metadata."
            )

        gap = float(
            next(
                iter(gaps)
            )
        )

        chirality = int(
            next(
                iter(
                    chiralities
                )
            )
        )

        end_seed_ids = sorted(
            [
                node_id
                for node_id in seed_ids
                if nodes[
                    node_id
                ][
                    "end"
                ]
                == end
            ],
            key=lambda node_id: as_int(
                nodes[
                    node_id
                ],
                "circumferential_index",
            ),
        )

        seed_positions = np.asarray(
            [
                positions[
                    node_id
                ]
                for node_id
                in end_seed_ids
            ],
            dtype=float,
        )

        seed_center = np.mean(
            seed_positions,
            axis=0,
        )

        outward = (
            -tube_axis
            if end == "LOWER"
            else tube_axis
        )

        reference = (
            seed_positions[0]
            - seed_center
        )

        reference -= (
            np.dot(
                reference,
                tube_axis,
            )
            * tube_axis
        )

        basis_x = unit(
            reference,
            np.asarray(
                [
                    1.0,
                    0.0,
                    0.0,
                ],
                dtype=float,
            ),
        )

        basis_y = unit(
            np.cross(
                tube_axis,
                basis_x,
            )
        )

        source_xy = np.asarray(
            [
                local_xy(
                    nodes[
                        row[
                            "annulus_node"
                        ]
                    ]
                )
                for row in end_paths
            ],
            dtype=float,
        )

        target_positions = np.asarray(
            [
                positions[
                    row[
                        "seed_node"
                    ]
                ]
                for row in end_paths
            ],
            dtype=float,
        )

        target_xy = np.asarray(
            [
                [
                    float(
                        np.dot(
                            position
                            - seed_center,
                            basis_x,
                        )
                    ),
                    float(
                        np.dot(
                            position
                            - seed_center,
                            basis_y,
                        )
                    ),
                ]
                for position
                in target_positions
            ],
            dtype=float,
        )

        (
            rotation,
            translation,
            fit_rms,
        ) = fit_orthogonal(
            source_xy,
            target_xy,
            chirality,
        )

        seed_axial = np.asarray(
            [
                float(
                    np.dot(
                        position
                        - seed_center,
                        outward,
                    )
                )
                for position
                in target_positions
            ],
            dtype=float,
        )

        plane_axial = (
            float(
                np.mean(
                    seed_axial
                )
            )
            + gap
        )

        annulus_ids = [
            node_id
            for node_id, row
            in nodes.items()
            if row[
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
        ]

        if len(annulus_ids) != 126:
            raise RuntimeError(
                f"{end}: unexpected annulus population."
            )

        mapped_xy = []

        for node_id in annulus_ids:
            xy_value = (
                local_xy(
                    nodes[
                        node_id
                    ]
                )
                @ rotation
                + translation
            )

            mapped_xy.append(
                xy_value
            )

            positions[node_id] = (
                seed_center
                + basis_x
                * xy_value[0]
                + basis_y
                * xy_value[1]
                + outward
                * plane_axial
            )

        mapped_xy_array = np.asarray(
            mapped_xy,
            dtype=float,
        )

        mean_xy = np.mean(
            mapped_xy_array,
            axis=0,
        )

        annulus_center = (
            seed_center
            + basis_x
            * mean_xy[0]
            + basis_y
            * mean_xy[1]
            + outward
            * plane_axial
        )

        end_geometry[end] = {
            "seed_center": seed_center,
            "outward": outward,
            "basis_x": basis_x,
            "basis_y": basis_y,
            "annulus_center": annulus_center,
            "selected_gap_nm": gap,
            "annulus_center_offset_nm": float(
                np.linalg.norm(
                    mean_xy
                )
            ),
            "endpoint_fit_rms_nm": fit_rms,
        }

    fixed_heavy_ids = [
        node_id
        for node_id, row
        in nodes.items()
        if row[
            "element"
        ]
        != "H"
        and row[
            "node_type"
        ]
        != "ALTERNATING_BN_TRIMER_BRIDGE"
    ]

    missing_fixed = [
        node_id
        for node_id in fixed_heavy_ids
        if node_id not in positions
    ]

    if missing_fixed:
        raise RuntimeError(
            "Fixed-heavy coordinate assignment incomplete: "
            + " | ".join(
                missing_fixed[:20]
            )
        )

    fixed_positions = [
        positions[
            node_id
        ]
        for node_id
        in fixed_heavy_ids
    ]

    conformer_rows: list[
        dict[str, Any]
    ] = []

    sorted_paths = sorted(
        path_rows,
        key=lambda row: (
            row[
                "end"
            ],
            as_int(
                row,
                "bridge_index",
            ),
        ),
    )

    for path in sorted_paths:
        end = path[
            "end"
        ]

        start_id = path[
            "seed_node"
        ]

        finish_id = path[
            "annulus_node"
        ]

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

        base_chain, metadata = select_base_chain(
            target_distance,
            library,
        )

        midpoint = (
            0.5
            * (
                start
                + finish
            )
        )

        radial = (
            midpoint
            - tube_center
        )

        radial -= (
            np.dot(
                radial,
                tube_axis,
            )
            * tube_axis
        )

        radial = unit(
            radial,
            end_geometry[
                end
            ][
                "basis_x"
            ],
        )

        fixed_array = np.asarray(
            fixed_positions,
            dtype=float,
        )

        best_chain: np.ndarray | None = None

        best_metadata: tuple[
            bool,
            float,
            tuple[
                int,
                float,
                float,
                float,
            ],
        ] | None = None

        for mirror in (
            False,
            True,
        ):
            for azimuth in np.arange(
                0.0,
                360.0,
                AZIMUTH_STEP,
            ):
                candidate = map_chain(
                    base_chain,
                    start,
                    finish,
                    float(
                        azimuth
                    ),
                    mirror,
                )

                score = clearance_score(
                    candidate[1:4],
                    fixed_array,
                    radial,
                    midpoint,
                )

                if (
                    best_metadata is None
                    or score
                    < best_metadata[2]
                ):
                    best_chain = candidate

                    best_metadata = (
                        mirror,
                        float(
                            azimuth
                        ),
                        score,
                    )

        if (
            best_chain is None
            or best_metadata is None
        ):
            raise RuntimeError(
                f"No conformer found for {path['bridge_path_id']}"
            )

        for bridge_id, point in zip(
            bridge_ids,
            best_chain[1:4],
        ):
            positions[
                bridge_id
            ] = np.asarray(
                point,
                dtype=float,
            )

            fixed_positions.append(
                positions[
                    bridge_id
                ]
            )

        lengths = np.linalg.norm(
            np.diff(
                best_chain,
                axis=0,
            ),
            axis=1,
        )

        internal_angles = [
            angle_deg(
                best_chain[
                    index - 1
                ],
                best_chain[
                    index
                ],
                best_chain[
                    index + 1
                ],
            )
            for index in (
                1,
                2,
                3,
            )
        ]

        conformer_rows.append(
            {
                "bridge_path_id": (
                    path[
                        "bridge_path_id"
                    ]
                ),
                "end": end,
                "bridge_index": (
                    path[
                        "bridge_index"
                    ]
                ),
                "seed_node": start_id,
                "annulus_node": (
                    finish_id
                ),
                "target_endpoint_distance_nm": (
                    target_distance
                ),
                "library_angle_deg": (
                    metadata[
                        "angle_deg"
                    ]
                ),
                "library_torsion_1_deg": (
                    metadata[
                        "phi1_deg"
                    ]
                ),
                "library_torsion_2_deg": (
                    metadata[
                        "phi2_deg"
                    ]
                ),
                "library_distance_nm": (
                    metadata[
                        "library_distance_nm"
                    ]
                ),
                "library_distance_error_nm": (
                    metadata[
                        "library_error_nm"
                    ]
                ),
                "selected_mirror": (
                    best_metadata[0]
                ),
                "selected_azimuth_deg": (
                    best_metadata[1]
                ),
                "preexisting_heavy_clash_count": (
                    best_metadata[
                        2
                    ][0]
                ),
                "minimum_preexisting_heavy_clearance_nm": (
                    -best_metadata[
                        2
                    ][3]
                ),
                "bond_1_nm": float(
                    lengths[0]
                ),
                "bond_2_nm": float(
                    lengths[1]
                ),
                "bond_3_nm": float(
                    lengths[2]
                ),
                "bond_4_nm": float(
                    lengths[3]
                ),
                "angle_1_deg": float(
                    internal_angles[0]
                ),
                "angle_2_deg": float(
                    internal_angles[1]
                ),
                "angle_3_deg": float(
                    internal_angles[2]
                ),
            }
        )

    if len(conformer_rows) != N_BRIDGES:
        raise RuntimeError(
            "Unexpected selected-conformer count."
        )

    hydrogen_ids = [
        node_id
        for node_id, row
        in nodes.items()
        if row[
            "element"
        ]
        == "H"
    ]

    if len(hydrogen_ids) != N_H:
        raise RuntimeError(
            "Unexpected H-node count."
        )

    for hydrogen_id in hydrogen_ids:
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

        heavy_neighbors = [
            neighbor
            for neighbor
            in adjacency[
                heavy_id
            ]
            if nodes[
                neighbor
            ][
                "element"
            ]
            != "H"
        ]

        if len(heavy_neighbors) != 2:
            raise RuntimeError(
                f"{hydrogen_id}: attached heavy atom "
                "does not have exactly two heavy neighbors."
            )

        heavy_position = positions[
            heavy_id
        ]

        vector_sum = sum(
            (
                unit(
                    positions[
                        neighbor
                    ]
                    - heavy_position
                )
                for neighbor
                in heavy_neighbors
            ),
            np.zeros(
                3,
                dtype=float,
            ),
        )

        end = row[
            "end"
        ]

        geometry = end_geometry[
            end
        ]

        radial = (
            heavy_position
            - geometry[
                "annulus_center"
            ]
        )

        radial -= (
            np.dot(
                radial,
                geometry[
                    "outward"
                ],
            )
            * geometry[
                "outward"
            ]
        )

        if row[
            "node_type"
        ] == "ANNULUS_INNER_PASSIVANT_H":
            fallback = -unit(
                radial,
                geometry[
                    "basis_x"
                ],
            )

        elif row[
            "node_type"
        ] == "ANNULUS_OUTER_PASSIVANT_H":
            fallback = unit(
                radial,
                geometry[
                    "basis_x"
                ],
            )

        elif row[
            "node_type"
        ] == "SEED_PASSIVANT_H":
            fallback = geometry[
                "outward"
            ]

        else:
            fallback = unit(
                radial,
                geometry[
                    "outward"
                ],
            )

        direction = unit(
            -vector_sum,
            fallback,
        )

        positions[
            hydrogen_id
        ] = (
            heavy_position
            + xh_length(
                nodes[
                    heavy_id
                ][
                    "element"
                ]
            )
            * direction
        )

    missing = [
        node_id
        for node_id in nodes
        if node_id not in positions
    ]

    if (
        missing
        or len(positions)
        != N_TOTAL
    ):
        raise RuntimeError(
            "Coordinate assignment incomplete: "
            + " | ".join(
                missing[:20]
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
                "node_id": (
                    node_id
                ),
                "element": (
                    nodes[
                        node_id
                    ][
                        "element"
                    ]
                ),
                "node_type": (
                    nodes[
                        node_id
                    ][
                        "node_type"
                    ]
                ),
                "end": (
                    nodes[
                        node_id
                    ][
                        "end"
                    ]
                ),
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
                    "STATIC_TRIMER_CONFORMER_EMBEDDING"
                ),
                "energy_minimized": False,
                "MD_relaxed": False,
            }
        )

    write_rows(
        COORDS_OUT,
        coordinate_rows,
    )

    write_rows(
        CONFORMERS_OUT,
        conformer_rows,
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

        if as_bool(
            row[
                "heavy_atom_edge"
            ]
        ):
            target = BN

        else:
            heavy = (
                first
                if nodes[
                    first
                ][
                    "element"
                ]
                != "H"
                else second
            )

            target = xh_length(
                nodes[
                    heavy
                ][
                    "element"
                ]
            )

        bond_rows.append(
            {
                "edge_id": (
                    row[
                        "edge_id"
                    ]
                ),
                "end": (
                    row[
                        "end"
                    ]
                ),
                "edge_type": (
                    row[
                        "edge_type"
                    ]
                ),
                "source_node": first,
                "target_node": second,
                "source_element": (
                    nodes[
                        first
                    ][
                        "element"
                    ]
                ),
                "target_element": (
                    nodes[
                        second
                    ][
                        "element"
                    ]
                ),
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
        BONDS_OUT,
        bond_rows,
    )

    grouped_bonds: dict[
        str,
        list[dict[str, Any]]
    ] = defaultdict(list)

    for row in bond_rows:
        grouped_bonds[
            row[
                "edge_type"
            ]
        ].append(row)

    bond_summary_rows = []

    for edge_type, rows in sorted(
        grouped_bonds.items()
    ):
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

        bond_summary_rows.append(
            {
                "edge_type": (
                    edge_type
                ),
                "edge_count": (
                    len(rows)
                ),
                "mean_length_nm": float(
                    np.mean(
                        lengths
                    )
                ),
                "minimum_length_nm": float(
                    np.min(
                        lengths
                    )
                ),
                "maximum_length_nm": float(
                    np.max(
                        lengths
                    )
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

    write_rows(
        BOND_SUMMARY_OUT,
        bond_summary_rows,
    )

    heavy_bonds = [
        row
        for row in bond_rows
        if nodes[
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
    ]

    h_bonds = [
        row
        for row in bond_rows
        if row not in heavy_bonds
    ]

    bridge_bonds = [
        row
        for row in bond_rows
        if row[
            "edge_type"
        ]
        == "ALTERNATING_BN_TRIMER_BRIDGE"
    ]

    max_bn_dev = max(
        float(
            row[
                "absolute_deviation_nm"
            ]
        )
        for row in heavy_bonds
    )

    max_xh_dev = max(
        float(
            row[
                "absolute_deviation_nm"
            ]
        )
        for row in h_bonds
    )

    max_bridge_dev = max(
        float(
            row[
                "absolute_deviation_nm"
            ]
        )
        for row in bridge_bonds
    )

    angle_values_by_type: dict[
        str,
        list[float]
    ] = defaultdict(list)

    critical_values = []

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
                value = angle_deg(
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
                ] in CRITICAL_TYPES:
                    critical_values.append(
                        value
                    )

    angle_summary_rows = []

    for node_type, values_list in sorted(
        angle_values_by_type.items()
    ):
        values = np.asarray(
            values_list,
            dtype=float,
        )

        angle_summary_rows.append(
            {
                "center_node_type": (
                    node_type
                ),
                "angle_count": int(
                    values.size
                ),
                "minimum_angle_deg": float(
                    np.min(
                        values
                    )
                ),
                "mean_angle_deg": float(
                    np.mean(
                        values
                    )
                ),
                "maximum_angle_deg": float(
                    np.max(
                        values
                    )
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
        ANGLE_SUMMARY_OUT,
        angle_summary_rows,
    )

    critical = np.asarray(
        critical_values,
        dtype=float,
    )

    critical_min = float(
        np.min(
            critical
        )
    )

    critical_mean = float(
        np.mean(
            critical
        )
    )

    critical_max = float(
        np.max(
            critical
        )
    )

    critical_rms = float(
        np.sqrt(
            np.mean(
                (
                    critical
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

    coords_array = np.asarray(
        [
            positions[
                node_id
            ]
            for node_id
            in ordered_ids
        ],
        dtype=float,
    )

    bonded_pairs = {
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
            MIN_HH
        ),
        "H_HEAVY": (
            MIN_HX
        ),
        "H_H": (
            MIN_HYD_HYD
        ),
    }

    minima = {
        category: math.inf
        for category
        in thresholds
    }

    minimum_pairs = {
        category: ""
        for category
        in thresholds
    }

    clash_counts = {
        category: 0
        for category
        in thresholds
    }

    for first_index in range(
        len(
            ordered_ids
        )
    ):
        distances = np.linalg.norm(
            coords_array[
                first_index + 1:
            ]
            - coords_array[
                first_index
            ],
            axis=1,
        )

        first_id = ordered_ids[
            first_index
        ]

        first_h = (
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
            ) in bonded_pairs:
                continue

            second_id = ordered_ids[
                second_index
            ]

            second_h = (
                nodes[
                    second_id
                ][
                    "element"
                ]
                == "H"
            )

            if first_h and second_h:
                category = "H_H"

            elif first_h or second_h:
                category = "H_HEAVY"

            else:
                category = "HEAVY_HEAVY"

            value = float(
                distance
            )

            if value < minima[
                category
            ]:
                minima[
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
                minima[
                    category
                ]
            ),
            "minimum_pair": (
                minimum_pairs[
                    category
                ]
            ),
            "clash_threshold_nm": (
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
        CONTACTS_OUT,
        contact_rows,
    )

    end_rows = []

    for end in (
        "LOWER",
        "UPPER",
    ):
        geometry = end_geometry[
            end
        ]

        outward = geometry[
            "outward"
        ]

        center = geometry[
            "annulus_center"
        ]

        inner_h_ids = [
            node_id
            for node_id, row
            in nodes.items()
            if row[
                "end"
            ]
            == end
            and row[
                "node_type"
            ]
            == "ANNULUS_INNER_PASSIVANT_H"
        ]

        outer_ids = [
            node_id
            for node_id, row
            in nodes.items()
            if row[
                "end"
            ]
            == end
            and row[
                "node_type"
            ]
            == "ANNULUS_OUTER_BOUNDARY"
        ]

        inner_radii = []
        outer_radii = []

        for node_id in inner_h_ids:
            displacement = (
                positions[
                    node_id
                ]
                - center
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

        for node_id in outer_ids:
            displacement = (
                positions[
                    node_id
                ]
                - center
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

        end_conformers = [
            row
            for row in conformer_rows
            if row[
                "end"
            ]
            == end
        ]

        endpoint_distances = np.asarray(
            [
                float(
                    row[
                        "target_endpoint_distance_nm"
                    ]
                )
                for row
                in end_conformers
            ],
            dtype=float,
        )

        bridge_deviations = np.asarray(
            [
                max(
                    abs(
                        float(
                            row[field]
                        )
                        - BN
                    )
                    for field in (
                        "bond_1_nm",
                        "bond_2_nm",
                        "bond_3_nm",
                        "bond_4_nm",
                    )
                )
                for row
                in end_conformers
            ],
            dtype=float,
        )

        aperture = (
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
                "bridge_paths": (
                    len(
                        end_conformers
                    )
                ),
                "selected_gap_nm": (
                    geometry[
                        "selected_gap_nm"
                    ]
                ),
                "annulus_center_offset_nm": (
                    geometry[
                        "annulus_center_offset_nm"
                    ]
                ),
                "annulus_endpoint_fit_RMS_nm": (
                    geometry[
                        "endpoint_fit_rms_nm"
                    ]
                ),
                "endpoint_distance_minimum_nm": float(
                    np.min(
                        endpoint_distances
                    )
                ),
                "endpoint_distance_mean_nm": float(
                    np.mean(
                        endpoint_distances
                    )
                ),
                "endpoint_distance_maximum_nm": float(
                    np.max(
                        endpoint_distances
                    )
                ),
                "maximum_bridge_bond_deviation_nm": float(
                    np.max(
                        bridge_deviations
                    )
                ),
                "conformers_with_preexisting_heavy_clashes": sum(
                    int(
                        row[
                            "preexisting_heavy_clash_count"
                        ]
                    )
                    > 0
                    for row
                    in end_conformers
                ),
                "minimum_selected_conformer_clearance_nm": min(
                    float(
                        row[
                            "minimum_preexisting_heavy_clearance_nm"
                        ]
                    )
                    for row
                    in end_conformers
                ),
                "inner_H_atoms": (
                    len(
                        inner_h_ids
                    )
                ),
                "nuclear_aperture_diameter_nm": (
                    aperture
                ),
                "nuclear_aperture_relative_error": abs(
                    aperture
                    - target_aperture
                ) / target_aperture,
                "outer_annulus_radius_mean_nm": (
                    outer_radius
                ),
                "outer_annulus_radius_relative_error": abs(
                    outer_radius
                    - target_outer_radius
                ) / target_outer_radius,
            }
        )

    write_rows(
        END_SUMMARY_OUT,
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

    asymmetry = max(
        abs(
            float(
                lower[
                    "endpoint_distance_mean_nm"
                ]
            )
            - float(
                upper[
                    "endpoint_distance_mean_nm"
                ]
            )
        ),
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
                    "maximum_bridge_bond_deviation_nm"
                ]
            )
            - float(
                upper[
                    "maximum_bridge_bond_deviation_nm"
                ]
            )
        ),
    )

    parent_unchanged = all(
        np.array_equal(
            positions[
                node_id
            ],
            source_coords[
                node_id
            ],
        )
        for node_id
        in parent_ids
    )

    max_library_error = max(
        float(
            row[
                "library_distance_error_nm"
            ]
        )
        for row in conformer_rows
    )

    gates = {
        "Gate3A_parent_audit_is_accepted": (
            parent_summary.get(
                "decision"
            )
            == EXPECTED_PARENT_DECISION
        ),
        "Gate3G1_direct_junction_is_rejected": (
            direct_summary.get(
                "decision"
            )
            == EXPECTED_DIRECT_DECISION
        ),
        "Gate3I_trimer_bridge_graph_is_accepted": (
            graph_summary.get(
                "decision"
            )
            == EXPECTED_GRAPH_DECISION
        ),
        "Gate3I_has_no_failed_gates": (
            len(
                failed_upstream
            )
            == 0
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
                for point
                in positions.values()
            )
        ),
        "parent_coordinates_are_unchanged": (
            parent_unchanged
        ),
        "30_explicit_trimer_conformers_were_embedded": (
            len(
                conformer_rows
            )
            == N_BRIDGES
        ),
        "library_endpoint_distance_error_is_within_0p0005nm": (
            max_library_error
            <= MAX_LIBRARY_DISTANCE_ERROR
        ),
        "all_BN_bonds_are_within_0p003nm_of_target": (
            max_bn_dev
            <= MAX_BN_DEV
        ),
        "all_bridge_BN_bonds_are_within_0p003nm_of_target": (
            max_bridge_dev
            <= MAX_BN_DEV
        ),
        "all_XH_bonds_are_within_0p002nm_of_target": (
            max_xh_dev
            <= MAX_XH_DEV
        ),
        "critical_valence_angle_minimum_is_at_least_70deg": (
            critical_min
            >= MIN_CRITICAL_ANGLE
        ),
        "critical_valence_angle_maximum_is_at_most_175deg": (
            critical_max
            <= MAX_CRITICAL_ANGLE
        ),
        "critical_valence_angle_RMS_deviation_is_at_most_30deg": (
            critical_rms
            <= MAX_CRITICAL_RMS
        ),
        "annulus_center_offsets_are_within_0p050nm": all(
            float(
                row[
                    "annulus_center_offset_nm"
                ]
            )
            <= MAX_CENTER_OFFSET
            for row
            in end_rows
        ),
        "aperture_errors_are_within10percent": all(
            float(
                row[
                    "nuclear_aperture_relative_error"
                ]
            )
            <= MAX_APERTURE_ERROR
            for row
            in end_rows
        ),
        "outer_radius_errors_are_within15percent": all(
            float(
                row[
                    "outer_annulus_radius_relative_error"
                ]
            )
            <= MAX_OUTER_RADIUS_ERROR
            for row
            in end_rows
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
        "lower_and_upper_embeddings_are_symmetric_within_0p010nm": (
            asymmetry
            <= MAX_END_ASYMMETRY
        ),
    }

    failed = [
        name
        for name, passed
        in gates.items()
        if not passed
    ]

    accepted = (
        len(failed)
        == 0
    )

    decision = (
        PASS_DECISION
        if accepted
        else REVIEW_DECISION
    )

    next_step = (
        "AUDIT_R2_ALTERNATING_BN_TRIMER_BRIDGE_"
        "CHEMICAL_REALIZABILITY_AND_PARAMETERIZATION_SCOPE"
        if accepted
        else
        "REFINE_R2_TRIMER_BRIDGE_CONFORMERS_AND_"
        "H_PASSIVANT_ORIENTATIONS"
    )

    summary = {
        "decision": decision,
        "coordinate_nodes": (
            len(positions)
        ),
        "parent_atoms": (
            N_PARENT
        ),
        "total_heavy_atoms": (
            N_HEAVY
        ),
        "total_H_atoms": (
            N_H
        ),
        "bridge_conformers": (
            len(
                conformer_rows
            )
        ),
        "conformer_library_size": int(
            library[
                "distance"
            ].size
        ),
        "maximum_library_distance_error_nm": (
            max_library_error
        ),
        "maximum_BN_bond_deviation_nm": (
            max_bn_dev
        ),
        "maximum_bridge_BN_bond_deviation_nm": (
            max_bridge_dev
        ),
        "maximum_XH_bond_deviation_nm": (
            max_xh_dev
        ),
        "critical_angle_minimum_deg": (
            critical_min
        ),
        "critical_angle_mean_deg": (
            critical_mean
        ),
        "critical_angle_maximum_deg": (
            critical_max
        ),
        "critical_angle_RMS_deviation_deg": (
            critical_rms
        ),
        "minimum_nonbonded_heavy_heavy_nm": (
            minima[
                "HEAVY_HEAVY"
            ]
        ),
        "minimum_nonbonded_H_heavy_nm": (
            minima[
                "H_HEAVY"
            ]
        ),
        "minimum_nonbonded_H_H_nm": (
            minima[
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
            target_aperture
        ),
        "maximum_lower_upper_asymmetry_nm": (
            asymmetry
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
                failed
            )
        ),
        "required_next_step": (
            next_step
        ),
    }

    write_rows(
        SUMMARY_OUT,
        [
            summary
        ],
    )

    write_rows(
        GATES_OUT,
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
                "conformer_model": {
                    "B_N_bond_nm": BN,
                    "internal_angle_grid_deg": [
                        ANGLE_MIN,
                        ANGLE_MAX,
                        ANGLE_STEP,
                    ],
                    "torsion_grid_deg": (
                        TORSION_STEP
                    ),
                    "azimuth_grid_deg": (
                        AZIMUTH_STEP
                    ),
                },
                "limitations": [
                    (
                        "This is a deterministic static conformer "
                        "embedding, not an energy minimization."
                    ),
                    (
                        "The conformer-selection score is based on "
                        "geometric clearance, not an energetic potential."
                    ),
                    (
                        "B-H and N-H distances are provisional geometry "
                        "targets, not force-field parameters."
                    ),
                    (
                        "Passing this gate does not establish chemical "
                        "stability or synthetic feasibility."
                    ),
                    (
                        "No topology, charges, force-field parameters, "
                        "minimization, MD, or QM calculation was generated."
                    ),
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    manifest = [
        {
            "role": (
                "Gate3A_parent_summary"
            ),
            "file": rel(
                PARENT_SUMMARY
            ),
            "sha256": digest(
                PARENT_SUMMARY
            ),
        },
        {
            "role": (
                "Gate3G_parent_and_seed_coordinates"
            ),
            "file": rel(
                SOURCE_COORDS
            ),
            "sha256": digest(
                SOURCE_COORDS
            ),
        },
        {
            "role": (
                "Gate3G1_direct_lower_bound_summary"
            ),
            "file": rel(
                DIRECT_SUMMARY
            ),
            "sha256": digest(
                DIRECT_SUMMARY
            ),
        },
        {
            "role": (
                "Gate3I_graph_nodes"
            ),
            "file": rel(
                GRAPH_NODES
            ),
            "sha256": digest(
                GRAPH_NODES
            ),
        },
        {
            "role": (
                "Gate3I_graph_edges"
            ),
            "file": rel(
                GRAPH_EDGES
            ),
            "sha256": digest(
                GRAPH_EDGES
            ),
        },
        {
            "role": (
                "Gate3I_bridge_paths"
            ),
            "file": rel(
                BRIDGE_PATHS
            ),
            "sha256": digest(
                BRIDGE_PATHS
            ),
        },
        {
            "role": (
                "Gate3I_graph_summary"
            ),
            "file": rel(
                GRAPH_SUMMARY
            ),
            "sha256": digest(
                GRAPH_SUMMARY
            ),
        },
    ]

    write_rows(
        MANIFEST_OUT,
        manifest,
    )

    xyz_lines = [
        str(
            len(
                ordered_ids
            )
        ),
        (
            "R2 alternating BN trimer-bridge static embedding; "
            "not energy minimized"
        ),
    ]

    for node_id in ordered_ids:
        coordinates = (
            positions[
                node_id
            ]
            * 10.0
        )

        xyz_lines.append(
            f"{nodes[node_id]['element']:2s} "
            f"{coordinates[0]: .8f} "
            f"{coordinates[1]: .8f} "
            f"{coordinates[2]: .8f}"
        )

    XYZ_OUT.write_text(
        "\n".join(
            xyz_lines
        )
        + "\n",
        encoding="utf-8",
    )

    pdb_lines = [
        "REMARK R2 ALTERNATING BN TRIMER-BRIDGE STATIC EMBEDDING",
        "REMARK NOT ENERGY MINIMIZED; NO FORCE-FIELD TOPOLOGY",
    ]

    for serial, node_id in enumerate(
        ordered_ids,
        start=1,
    ):
        row = nodes[
            node_id
        ]

        coordinates = (
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
            f"{coordinates[0]:8.3f}"
            f"{coordinates[1]:8.3f}"
            f"{coordinates[2]:8.3f}"
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

    REPORT_OUT.write_text(
        f"""# R2 Alternating BN Trimer-Bridge Static Coordinate Embedding

## Scope

This gate embeds all 30 alternating BN trimer bridges and all 174 H
passivants while retaining the accepted parent coordinates.

No topology, formal charges, force-field parameters, minimization, MD,
or QM calculation was generated.

## Coordinate inventory

- Parent atoms: **{N_PARENT}**
- Total heavy atoms: **{N_HEAVY}**
- H atoms: **{N_H}**
- Total nodes: **{len(positions)}**
- Explicit bridge conformers: **{len(conformer_rows)}**
- Parent coordinates changed:
  **{'NO' if parent_unchanged else 'YES'}**

## Bond geometry

- Maximum B-N deviation:
  **{max_bn_dev:.6f} nm**
- Maximum bridge B-N deviation:
  **{max_bridge_dev:.6f} nm**
- Maximum X-H deviation:
  **{max_xh_dev:.6f} nm**
- Maximum conformer-library distance mismatch:
  **{max_library_error:.6f} nm**

## Critical valence angles

- Minimum/mean/maximum:
  **{critical_min:.3f}/{critical_mean:.3f}/{critical_max:.3f} degrees**
- RMS deviation from 120 degrees:
  **{critical_rms:.3f} degrees**

## Nonbonded contacts

- Heavy-heavy minimum/clashes:
  **{minima['HEAVY_HEAVY']:.6f}/{clash_counts['HEAVY_HEAVY']}**
- H-heavy minimum/clashes:
  **{minima['H_HEAVY']:.6f}/{clash_counts['H_HEAVY']}**
- H-H minimum/clashes:
  **{minima['H_H']:.6f}/{clash_counts['H_H']}**

## Aperture

- Target diameter:
  **{target_aperture:.6f} nm**
- Lower/upper H-defined diameter:
  **{float(lower['nuclear_aperture_diameter_nm']):.6f}/
  {float(upper['nuclear_aperture_diameter_nm']):.6f} nm**

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed else ' | '.join(failed)}**
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
  `{next_step}`

## Interpretation

This is a deterministic static conformer embedding. Passing the gate
would establish geometric consistency only, not energetic stability,
synthetic feasibility, or force-field validity.
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 alternating BN trimer-bridge static "
        "coordinate embedding completed."
    )

    print(
        "Coordinate nodes parent/heavy/H/total: "
        f"{N_PARENT}/"
        f"{N_HEAVY}/"
        f"{N_H}/"
        f"{len(positions)}"
    )

    print(
        "Parent coordinates unchanged: "
        f"{'YES' if parent_unchanged else 'NO'}"
    )

    print(
        "Conformer library / selected bridges / "
        "maximum endpoint mismatch: "
        f"{library['distance'].size}/"
        f"{len(conformer_rows)}/"
        f"{max_library_error:.6f} nm"
    )

    print(
        "Maximum BN / bridge-BN / XH bond deviations: "
        f"{max_bn_dev:.6f}/"
        f"{max_bridge_dev:.6f}/"
        f"{max_xh_dev:.6f} nm"
    )

    print(
        "Critical angles min/mean/max/RMSdev120: "
        f"{critical_min:.3f}/"
        f"{critical_mean:.3f}/"
        f"{critical_max:.3f}/"
        f"{critical_rms:.3f} deg"
    )

    for row in end_rows:
        print(
            f"{row['end']} gap/center-offset/"
            "endpoint-fit-RMS/endpoint-distance-min-max/"
            "aperture/error/outer-radius/error: "
            f"{float(row['selected_gap_nm']):.6f}/"
            f"{float(row['annulus_center_offset_nm']):.6f}/"
            f"{float(row['annulus_endpoint_fit_RMS_nm']):.6f}/"
            f"{float(row['endpoint_distance_minimum_nm']):.6f}-"
            f"{float(row['endpoint_distance_maximum_nm']):.6f}/"
            f"{float(row['nuclear_aperture_diameter_nm']):.6f}/"
            f"{float(row['nuclear_aperture_relative_error']):.6f}/"
            f"{float(row['outer_annulus_radius_mean_nm']):.6f}/"
            f"{float(row['outer_annulus_radius_relative_error']):.6f}"
        )

        print(
            f"{row['end']} conformers with preexisting "
            "heavy clashes / minimum selected clearance: "
            f"{row['conformers_with_preexisting_heavy_clashes']}/"
            f"{float(row['minimum_selected_conformer_clearance_nm']):.6f} nm"
        )

    print(
        "Minimum nonbonded heavy-heavy / H-heavy / H-H: "
        f"{minima['HEAVY_HEAVY']:.6f}/"
        f"{minima['H_HEAVY']:.6f}/"
        f"{minima['H_H']:.6f} nm"
    )

    print(
        "Clash counts heavy-heavy / H-heavy / H-H: "
        f"{clash_counts['HEAVY_HEAVY']}/"
        f"{clash_counts['H_HEAVY']}/"
        f"{clash_counts['H_H']}"
    )

    print(
        "Maximum lower-upper asymmetry: "
        f"{asymmetry:.6f} nm"
    )

    print(
        f"Decision: {decision}"
    )

    print(
        "Failed gates: "
        + (
            "NONE"
            if not failed
            else " | ".join(
                failed
            )
        )
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
        f"Required next step: {next_step}"
    )

    for path in (
        COORDS_OUT,
        CONFORMERS_OUT,
        BONDS_OUT,
        BOND_SUMMARY_OUT,
        ANGLE_SUMMARY_OUT,
        CONTACTS_OUT,
        END_SUMMARY_OUT,
        SUMMARY_OUT,
        GATES_OUT,
        JSON_OUT,
        MANIFEST_OUT,
        XYZ_OUT,
        PDB_OUT,
        REPORT_OUT,
    ):
        print(
            f"Wrote: {rel(path)}"
        )


if __name__ == "__main__":
    main()
