#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

AUDIT_MODULE_PATH = (
    ROOT
    / "scripts/phase1A/"
    "audit_day024_r2_parent_rim_and_chemical_constraints.py"
)

SYSTEM_GRO = (
    ROOT
    / "runs/phase1A/day023_confinement_design/"
    "15_r2_frozen_solute_nvt_20ps_preparation/"
    "r2_frozen_solute_nvt_20ps_input.gro"
)

GATE3A_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "01_r2_parent_rim_chemical_audit"
)

GATE3F_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "07_r2_reconstruction_vs_partial_attachment_contingency"
)

OUTPUT_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "08_r2_partial_attachment_annulus_static_coordinate_embedding"
)

PARENT_SUMMARY_CSV = (
    GATE3A_ROOT
    / "r2_parent_rim_chemical_audit_summary.csv"
)

TERMINAL_ATOMS_CSV = (
    GATE3A_ROOT
    / "r2_parent_terminal_rim_atoms.csv"
)

DESIGN_NODES_CSV = (
    GATE3F_ROOT
    / "r2_partial_attachment_passivated_annulus_nodes.csv"
)

DESIGN_EDGES_CSV = (
    GATE3F_ROOT
    / "r2_partial_attachment_passivated_annulus_edges.csv"
)

DESIGN_SUMMARY_CSV = (
    GATE3F_ROOT
    / "r2_reconstruction_vs_partial_attachment_summary.csv"
)

DESIGN_GATES_CSV = (
    GATE3F_ROOT
    / "r2_reconstruction_vs_partial_attachment_gates.csv"
)

COORDINATES_CSV = (
    OUTPUT_ROOT
    / "r2_partial_attachment_static_coordinates.csv"
)

BOND_LENGTHS_CSV = (
    OUTPUT_ROOT
    / "r2_partial_attachment_static_bond_lengths.csv"
)

BOND_TYPE_SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_partial_attachment_static_bond_type_summary.csv"
)

JUNCTION_ANGLES_CSV = (
    OUTPUT_ROOT
    / "r2_partial_attachment_static_junction_angles.csv"
)

END_SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_partial_attachment_static_embedding_end_summary.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_partial_attachment_static_embedding_summary.csv"
)

GATES_CSV = (
    OUTPUT_ROOT
    / "r2_partial_attachment_static_embedding_gates.csv"
)

EMBEDDING_JSON = (
    OUTPUT_ROOT
    / "r2_partial_attachment_static_embedding.json"
)

SOURCE_MANIFEST_CSV = (
    OUTPUT_ROOT
    / "r2_partial_attachment_static_embedding_source_manifest.csv"
)

XYZ_FILE = (
    OUTPUT_ROOT
    / "r2_partial_attachment_static_embedding.xyz"
)

PDB_FILE = (
    OUTPUT_ROOT
    / "r2_partial_attachment_static_embedding.pdb"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_PARTIAL_ATTACHMENT_STATIC_COORDINATE_EMBEDDING_DAY024.md"
)

EXPECTED_PARENT_DECISION = (
    "R2_PARENT_RIM_AND_CHEMICAL_CONSTRAINTS_AUDITED"
)

EXPECTED_DESIGN_DECISION = (
    "R2_PARTIAL_HETEROPOLAR_ANNULUS_ATTACHMENT_AND_"
    "COMPLEMENTARY_PASSIVATION_GRAPH_VALIDATED"
)

PASS_DECISION = (
    "R2_PARTIAL_ATTACHMENT_ANNULUS_STATIC_COORDINATE_EMBEDDING_VALIDATED"
)

REVIEW_DECISION = (
    "R2_PARTIAL_ATTACHMENT_ANNULUS_STATIC_COORDINATE_EMBEDDING_"
    "REQUIRES_CONSTRAINED_GEOMETRIC_OPTIMIZATION"
)

EXPECTED_PARENT_ATOMS = 1680
EXPECTED_ADDED_HEAVY_ATOMS = 312
EXPECTED_ADDED_H_ATOMS = 84
EXPECTED_TOTAL_COORDINATE_NODES = (
    EXPECTED_PARENT_ATOMS
    + EXPECTED_ADDED_HEAVY_ATOMS
    + EXPECTED_ADDED_H_ATOMS
)

BN_TARGET_NM = 0.144973
BH_TARGET_NM = 0.119
NH_TARGET_NM = 0.101

MAX_PARENT_SEED_BOND_DEVIATION_NM = 0.003
MAX_ANNULUS_BOND_DEVIATION_NM = 0.002

MAX_ATTACHMENT_BOND_RMS_DEVIATION_NM = 0.020
MAX_ATTACHMENT_BOND_MAX_DEVIATION_NM = 0.035

MAX_XH_BOND_DEVIATION_NM = 0.002

MIN_JUNCTION_ANGLE_DEG = 70.0
MAX_JUNCTION_ANGLE_DEG = 170.0
MAX_JUNCTION_ANGLE_RMS_DEVIATION_DEG = 35.0

MAX_ANNULUS_CENTER_OFFSET_NM = 0.050
MAX_ATTACHMENT_PLANE_GAP_NM = 0.250

MAX_NUCLEAR_APERTURE_RELATIVE_ERROR = 0.10
MAX_OUTER_RADIUS_RELATIVE_ERROR = 0.15

MIN_NONBONDED_HEAVY_HEAVY_NM = 0.120
MIN_NONBONDED_H_HEAVY_NM = 0.070
MIN_NONBONDED_H_H_NM = 0.060

MAX_END_ASYMMETRY_NM = 1.0e-5


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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def read_csv_rows(
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


def load_audit_module():
    require_file(
        AUDIT_MODULE_PATH
    )

    spec = importlib.util.spec_from_file_location(
        "day024_parent_rim_audit",
        AUDIT_MODULE_PATH,
    )

    if (
        spec is None
        or spec.loader is None
    ):
        raise RuntimeError(
            "Could not load the Day024 audit module."
        )

    module = importlib.util.module_from_spec(
        spec
    )

    spec.loader.exec_module(
        module
    )

    return module


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


def physical_local_xy_nm(
    lattice_x: int,
    lattice_y: int,
    bond_length_nm: float,
) -> np.ndarray:
    return np.asarray(
        [
            lattice_x
            * bond_length_nm
            / 2.0,
            lattice_y
            * math.sqrt(3.0)
            * bond_length_nm
            / 2.0,
        ],
        dtype=float,
    )


def procrustes_2d(
    source_xy: np.ndarray,
    target_xy: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
    float,
]:
    if (
        source_xy.ndim != 2
        or target_xy.ndim != 2
        or source_xy.shape != target_xy.shape
        or source_xy.shape[1] != 2
    ):
        raise RuntimeError(
            "Invalid 2D Procrustes arrays."
        )

    source_center = np.mean(
        source_xy,
        axis=0,
    )

    target_center = np.mean(
        target_xy,
        axis=0,
    )

    source_centered = (
        source_xy
        - source_center
    )

    target_centered = (
        target_xy
        - target_center
    )

    covariance = (
        source_centered.T
        @ target_centered
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

    translation = (
        target_center
        - source_center
        @ rotation
    )

    fitted = (
        source_xy
        @ rotation
        + translation
    )

    rms = float(
        np.sqrt(
            np.mean(
                np.sum(
                    (
                        fitted
                        - target_xy
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


def optimize_axial_gap(
    lateral_residuals_nm: np.ndarray,
    seed_axial_coordinates_nm: np.ndarray,
    target_bond_nm: float,
) -> tuple[
    float,
    float,
    float,
]:
    mean_seed_axial = float(
        np.mean(
            seed_axial_coordinates_nm
        )
    )

    gaps = np.linspace(
        0.0,
        MAX_ATTACHMENT_PLANE_GAP_NM,
        25001,
    )

    annulus_planes = (
        mean_seed_axial
        + gaps
    )

    axial_differences = (
        annulus_planes[
            None,
            :,
        ]
        - seed_axial_coordinates_nm[
            :,
            None,
        ]
    )

    bond_lengths = np.sqrt(
        lateral_residuals_nm[
            :,
            None,
        ]
        ** 2
        + axial_differences
        ** 2
    )

    deviations = (
        bond_lengths
        - target_bond_nm
    )

    rms_by_gap = np.sqrt(
        np.mean(
            deviations
            ** 2,
            axis=0,
        )
    )

    best_index = int(
        np.argmin(
            rms_by_gap
        )
    )

    best_gap = float(
        gaps[
            best_index
        ]
    )

    best_rms = float(
        rms_by_gap[
            best_index
        ]
    )

    best_max = float(
        np.max(
            np.abs(
                deviations[
                    :,
                    best_index
                ]
            )
        )
    )

    return (
        best_gap,
        best_rms,
        best_max,
    )


def xh_target_nm(
    heavy_element: str,
) -> float:
    if heavy_element == "B":
        return BH_TARGET_NM

    if heavy_element == "N":
        return NH_TARGET_NM

    raise RuntimeError(
        f"Unexpected H-attached element: {heavy_element}"
    )


def angle_degrees(
    first: np.ndarray,
    center: np.ndarray,
    second: np.ndarray,
) -> float:
    vector_a = normalized(
        first - center
    )

    vector_b = normalized(
        second - center
    )

    cosine = float(
        np.clip(
            np.dot(
                vector_a,
                vector_b,
            ),
            -1.0,
            1.0,
        )
    )

    return math.degrees(
        math.acos(cosine)
    )


def edge_type_statistics(
    rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    grouped: dict[
        str,
        list[dict[str, Any]]
    ] = defaultdict(list)

    for row in rows:
        grouped[
            str(
                row[
                    "edge_type"
                ]
            )
        ].append(row)

    summaries = []

    for edge_type in sorted(
        grouped
    ):
        group = grouped[
            edge_type
        ]

        lengths = np.asarray(
            [
                float(
                    row[
                        "length_nm"
                    ]
                )
                for row in group
            ],
            dtype=float,
        )

        targets = np.asarray(
            [
                float(
                    row[
                        "target_length_nm"
                    ]
                )
                for row in group
            ],
            dtype=float,
        )

        deviations = (
            lengths
            - targets
        )

        summaries.append(
            {
                "edge_type": edge_type,
                "edge_count": len(
                    group
                ),
                "mean_length_nm": float(
                    np.mean(lengths)
                ),
                "minimum_length_nm": float(
                    np.min(lengths)
                ),
                "maximum_length_nm": float(
                    np.max(lengths)
                ),
                "mean_target_length_nm": float(
                    np.mean(targets)
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
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        AUDIT_MODULE_PATH,
        SYSTEM_GRO,
        PARENT_SUMMARY_CSV,
        TERMINAL_ATOMS_CSV,
        DESIGN_NODES_CSV,
        DESIGN_EDGES_CSV,
        DESIGN_SUMMARY_CSV,
        DESIGN_GATES_CSV,
    ):
        require_file(required)

    module = load_audit_module()

    gro_atoms, box = module.read_gro(
        SYSTEM_GRO
    )

    if len(gro_atoms) < EXPECTED_PARENT_ATOMS:
        raise RuntimeError(
            "The accepted GRO does not contain "
            "the expected parent HBN atoms."
        )

    parent_summary = read_single_csv_row(
        PARENT_SUMMARY_CSV
    )

    design_summary = read_single_csv_row(
        DESIGN_SUMMARY_CSV
    )

    design_gate_rows = read_csv_rows(
        DESIGN_GATES_CSV
    )

    terminal_rows = read_csv_rows(
        TERMINAL_ATOMS_CSV
    )

    node_rows = read_csv_rows(
        DESIGN_NODES_CSV
    )

    edge_rows = read_csv_rows(
        DESIGN_EDGES_CSV
    )

    if parent_summary.get(
        "decision"
    ) != EXPECTED_PARENT_DECISION:
        raise RuntimeError(
            "Gate 3A is not accepted."
        )

    if design_summary.get(
        "decision"
    ) != EXPECTED_DESIGN_DECISION:
        raise RuntimeError(
            "Gate 3F is not accepted."
        )

    failed_design_gates = [
        row.get(
            "gate",
            "",
        )
        for row in design_gate_rows
        if not parse_bool(
            row.get(
                "pass",
                "false",
            )
        )
    ]

    if failed_design_gates:
        raise RuntimeError(
            "Gate 3F still contains failed gates: "
            + " | ".join(
                failed_design_gates
            )
        )

    target_aperture_diameter_nm = parse_float(
        parent_summary,
        "target_aperture_diameter_nm",
    )

    target_aperture_radius_nm = parse_float(
        parent_summary,
        "target_aperture_radius_nm",
    )

    target_parent_rim_radius_nm = parse_float(
        parent_summary,
        "parent_rim_mean_radius_nm",
    )

    parent_positions = np.asarray(
        [
            gro_atoms[index][
                "position_nm"
            ]
            for index in range(
                EXPECTED_PARENT_ATOMS
            )
        ],
        dtype=float,
    )

    (
        tube_center,
        tube_axis,
        pca_eigenvalues,
    ) = module.determine_tube_axis(
        parent_positions
    )

    tube_axis = normalized(
        np.asarray(
            tube_axis,
            dtype=float,
        )
    )

    nodes_by_id = {
        row["node_id"]: row
        for row in node_rows
    }

    if len(nodes_by_id) != len(
        node_rows
    ):
        raise RuntimeError(
            "Duplicate node identifiers in Gate 3F."
        )

    adjacency: dict[
        str,
        set[str]
    ] = {
        node_id: set()
        for node_id in nodes_by_id
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
            first not in nodes_by_id
            or second not in nodes_by_id
        ):
            raise RuntimeError(
                "Gate 3F edge references a missing node."
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
                f"Duplicate design edge: {pair}"
            )

        edge_pairs.add(
            pair
        )

        adjacency[first].add(
            second
        )

        adjacency[second].add(
            first
        )

    positions: dict[
        str,
        np.ndarray
    ] = {}

    for node_id, row in nodes_by_id.items():
        if row[
            "node_type"
        ] != "PARENT_HBN":
            continue

        index = int(
            float(
                row[
                    "source_index_0based"
                ]
            )
        )

        if not (
            0
            <= index
            < EXPECTED_PARENT_ATOMS
        ):
            raise RuntimeError(
                f"Invalid parent index for {node_id}"
            )

        positions[
            node_id
        ] = np.array(
            parent_positions[
                index
            ],
            dtype=float,
            copy=True,
        )

    terminal_order_by_end: dict[
        str,
        list[int]
    ] = {}

    for end in (
        "LOWER",
        "UPPER",
    ):
        rows = [
            row
            for row in terminal_rows
            if row.get(
                "end"
            )
            == end
        ]

        rows.sort(
            key=lambda row: int(
                float(
                    row[
                        "circumferential_order"
                    ]
                )
            )
        )

        terminal_order_by_end[
            end
        ] = [
            int(
                float(
                    row[
                        "hbn_local_index_0based"
                    ]
                )
            )
            for row in rows
        ]

    embedding_end_rows: list[
        dict[str, Any]
    ] = []

    annulus_centers: dict[
        str,
        np.ndarray
    ] = {}

    annulus_plane_axes: dict[
        str,
        tuple[
            np.ndarray,
            np.ndarray,
            np.ndarray,
        ]
    ] = {}

    for end in (
        "LOWER",
        "UPPER",
    ):
        outward = (
            -tube_axis
            if end == "LOWER"
            else tube_axis
        )

        terminal_parent_ids = [
            f"P:{index}"
            for index
            in terminal_order_by_end[
                end
            ]
        ]

        terminal_positions = np.asarray(
            [
                positions[
                    node_id
                ]
                for node_id
                in terminal_parent_ids
            ],
            dtype=float,
        )

        terminal_center = np.mean(
            terminal_positions,
            axis=0,
        )

        radial_reference = (
            terminal_positions[0]
            - terminal_center
        )

        radial_reference = (
            radial_reference
            - np.dot(
                radial_reference,
                tube_axis,
            )
            * tube_axis
        )

        basis_x = normalized(
            radial_reference,
            fallback=np.asarray(
                [
                    1.0,
                    0.0,
                    0.0,
                ],
                dtype=float,
            ),
        )

        basis_y = normalized(
            np.cross(
                tube_axis,
                basis_x,
            )
        )

        seed_ids = sorted(
            (
                node_id
                for node_id, row
                in nodes_by_id.items()
                if row[
                    "node_type"
                ]
                == "HEXAGONAL_EDGE_COMPLETION_SEED"
                and row[
                    "end"
                ]
                == end
            ),
            key=lambda node_id: int(
                float(
                    nodes_by_id[
                        node_id
                    ][
                        "circumferential_index"
                    ]
                )
            ),
        )

        if len(seed_ids) != 30:
            raise RuntimeError(
                f"{end}: unexpected seed count "
                f"{len(seed_ids)}/30"
            )

        seed_parent_pair_distances = []

        for seed_id in seed_ids:
            parent_neighbors = [
                neighbor
                for neighbor in adjacency[
                    seed_id
                ]
                if nodes_by_id[
                    neighbor
                ][
                    "node_type"
                ]
                == "PARENT_HBN"
            ]

            if len(parent_neighbors) != 2:
                raise RuntimeError(
                    f"{seed_id}: expected two parent neighbors."
                )

            first = positions[
                parent_neighbors[0]
            ]

            second = positions[
                parent_neighbors[1]
            ]

            pair_vector = (
                second
                - first
            )

            pair_distance = float(
                np.linalg.norm(
                    pair_vector
                )
            )

            if pair_distance >= (
                2.0
                * BN_TARGET_NM
            ):
                raise RuntimeError(
                    f"{seed_id}: parent pair separation "
                    f"{pair_distance:.6f} nm is incompatible "
                    "with two BN target bonds."
                )

            pair_unit = normalized(
                pair_vector
            )

            placement_direction = (
                outward
                - np.dot(
                    outward,
                    pair_unit,
                )
                * pair_unit
            )

            placement_direction = normalized(
                placement_direction,
                fallback=outward,
            )

            height_nm = math.sqrt(
                max(
                    0.0,
                    BN_TARGET_NM
                    ** 2
                    - (
                        pair_distance
                        / 2.0
                    )
                    ** 2,
                )
            )

            midpoint = (
                first
                + second
            ) / 2.0

            positions[
                seed_id
            ] = (
                midpoint
                + placement_direction
                * height_nm
            )

            seed_parent_pair_distances.append(
                pair_distance
            )

        attachment_edges = [
            row
            for row in edge_rows
            if row[
                "edge_type"
            ]
            == "PARTIAL_HETEROPOLAR_SEED_TO_ANNULUS"
            and row[
                "end"
            ]
            == end
        ]

        if len(attachment_edges) != 15:
            raise RuntimeError(
                f"{end}: unexpected attachment count "
                f"{len(attachment_edges)}/15"
            )

        attachment_seed_ids = []
        attachment_annulus_ids = []

        for row in attachment_edges:
            first = row[
                "source_node"
            ]

            second = row[
                "target_node"
            ]

            if nodes_by_id[
                first
            ][
                "node_type"
            ] == "HEXAGONAL_EDGE_COMPLETION_SEED":
                seed_id = first
                annulus_id = second
            else:
                seed_id = second
                annulus_id = first

            attachment_seed_ids.append(
                seed_id
            )

            attachment_annulus_ids.append(
                annulus_id
            )

        seed_target_xy = np.asarray(
            [
                [
                    float(
                        np.dot(
                            positions[
                                seed_id
                            ]
                            - terminal_center,
                            basis_x,
                        )
                    ),
                    float(
                        np.dot(
                            positions[
                                seed_id
                            ]
                            - terminal_center,
                            basis_y,
                        )
                    ),
                ]
                for seed_id
                in attachment_seed_ids
            ],
            dtype=float,
        )

        annulus_source_xy = np.asarray(
            [
                physical_local_xy_nm(
                    int(
                        float(
                            nodes_by_id[
                                annulus_id
                            ][
                                "lattice_x"
                            ]
                        )
                    ),
                    int(
                        float(
                            nodes_by_id[
                                annulus_id
                            ][
                                "lattice_y"
                            ]
                        )
                    ),
                    BN_TARGET_NM,
                )
                for annulus_id
                in attachment_annulus_ids
            ],
            dtype=float,
        )

        (
            rotation,
            translation,
            procrustes_rms_nm,
        ) = procrustes_2d(
            annulus_source_xy,
            seed_target_xy,
        )

        fitted_attachment_xy = (
            annulus_source_xy
            @ rotation
            + translation
        )

        lateral_residuals_nm = (
            np.linalg.norm(
                fitted_attachment_xy
                - seed_target_xy,
                axis=1,
            )
        )

        seed_axial_coordinates_nm = np.asarray(
            [
                float(
                    np.dot(
                        positions[
                            seed_id
                        ]
                        - terminal_center,
                        outward,
                    )
                )
                for seed_id
                in attachment_seed_ids
            ],
            dtype=float,
        )

        (
            attachment_plane_gap_nm,
            predicted_attachment_rms_nm,
            predicted_attachment_max_nm,
        ) = optimize_axial_gap(
            lateral_residuals_nm,
            seed_axial_coordinates_nm,
            BN_TARGET_NM,
        )

        mean_seed_axial_nm = float(
            np.mean(
                seed_axial_coordinates_nm
            )
        )

        annulus_axial_coordinate_nm = (
            mean_seed_axial_nm
            + attachment_plane_gap_nm
        )

        annulus_center = (
            terminal_center
            + basis_x
            * translation[0]
            + basis_y
            * translation[1]
            + outward
            * annulus_axial_coordinate_nm
        )

        annulus_centers[
            end
        ] = annulus_center

        annulus_plane_axes[
            end
        ] = (
            basis_x,
            basis_y,
            outward,
        )

        annulus_ids = [
            node_id
            for node_id, row
            in nodes_by_id.items()
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
                f"{end}: unexpected annulus count "
                f"{len(annulus_ids)}/126"
            )

        for annulus_id in annulus_ids:
            row = nodes_by_id[
                annulus_id
            ]

            local_xy = (
                physical_local_xy_nm(
                    int(
                        float(
                            row[
                                "lattice_x"
                            ]
                        )
                    ),
                    int(
                        float(
                            row[
                                "lattice_y"
                            ]
                        )
                    ),
                    BN_TARGET_NM,
                )
            )

            mapped_xy = (
                local_xy
                @ rotation
                + translation
            )

            positions[
                annulus_id
            ] = (
                terminal_center
                + basis_x
                * mapped_xy[0]
                + basis_y
                * mapped_xy[1]
                + outward
                * annulus_axial_coordinate_nm
            )

        embedding_end_rows.append(
            {
                "end": end,
                "seed_atoms": len(
                    seed_ids
                ),
                "annulus_atoms": len(
                    annulus_ids
                ),
                "attachment_edges": len(
                    attachment_edges
                ),
                "mean_parent_pair_distance_nm": float(
                    np.mean(
                        seed_parent_pair_distances
                    )
                ),
                "minimum_parent_pair_distance_nm": float(
                    np.min(
                        seed_parent_pair_distances
                    )
                ),
                "maximum_parent_pair_distance_nm": float(
                    np.max(
                        seed_parent_pair_distances
                    )
                ),
                "Procrustes_RMS_lateral_residual_nm": (
                    procrustes_rms_nm
                ),
                "maximum_lateral_residual_nm": float(
                    np.max(
                        lateral_residuals_nm
                    )
                ),
                "annulus_rotation_determinant": float(
                    np.linalg.det(
                        rotation
                    )
                ),
                "annulus_center_offset_nm": float(
                    np.linalg.norm(
                        translation
                    )
                ),
                "attachment_plane_gap_nm": (
                    attachment_plane_gap_nm
                ),
                "predicted_attachment_RMS_deviation_nm": (
                    predicted_attachment_rms_nm
                ),
                "predicted_attachment_max_deviation_nm": (
                    predicted_attachment_max_nm
                ),
                "annulus_axial_coordinate_from_parent_plane_nm": (
                    annulus_axial_coordinate_nm
                ),
            }
        )

    hydrogen_ids = [
        node_id
        for node_id, row
        in nodes_by_id.items()
        if row[
            "element"
        ]
        == "H"
    ]

    if len(hydrogen_ids) != EXPECTED_ADDED_H_ATOMS:
        raise RuntimeError(
            "Unexpected hydrogen count: "
            f"{len(hydrogen_ids)}/"
            f"{EXPECTED_ADDED_H_ATOMS}"
        )

    for hydrogen_id in hydrogen_ids:
        row = nodes_by_id[
            hydrogen_id
        ]

        attached_to = row.get(
            "attached_to",
            "",
        )

        if (
            not attached_to
            or attached_to
            not in positions
        ):
            raise RuntimeError(
                f"{hydrogen_id}: invalid attached_to node."
            )

        heavy_position = positions[
            attached_to
        ]

        heavy_element = nodes_by_id[
            attached_to
        ][
            "element"
        ]

        bond_length_nm = xh_target_nm(
            heavy_element
        )

        end = row[
            "end"
        ]

        basis_x, basis_y, outward = (
            annulus_plane_axes[
                end
            ]
        )

        node_type = row[
            "node_type"
        ]

        if node_type == "SEED_PASSIVANT_H":
            heavy_neighbors = [
                neighbor
                for neighbor in adjacency[
                    attached_to
                ]
                if nodes_by_id[
                    neighbor
                ][
                    "element"
                ]
                != "H"
            ]

            if len(heavy_neighbors) != 2:
                raise RuntimeError(
                    f"{hydrogen_id}: seed heavy atom does not "
                    "have exactly two heavy neighbors."
                )

            unit_sum = np.zeros(
                3,
                dtype=float,
            )

            for neighbor in heavy_neighbors:
                unit_sum += normalized(
                    positions[
                        neighbor
                    ]
                    - heavy_position
                )

            direction = normalized(
                -unit_sum,
                fallback=outward,
            )

        elif node_type == "ANNULUS_OUTER_PASSIVANT_H":
            radial = (
                heavy_position
                - annulus_centers[
                    end
                ]
            )

            radial = (
                radial
                - np.dot(
                    radial,
                    outward,
                )
                * outward
            )

            direction = normalized(
                radial,
                fallback=basis_x,
            )

        elif node_type == "ANNULUS_INNER_PASSIVANT_H":
            radial = (
                heavy_position
                - annulus_centers[
                    end
                ]
            )

            radial = (
                radial
                - np.dot(
                    radial,
                    outward,
                )
                * outward
            )

            direction = normalized(
                -radial,
                fallback=-basis_x,
            )

        else:
            raise RuntimeError(
                f"Unexpected H node type: {node_type}"
            )

        positions[
            hydrogen_id
        ] = (
            heavy_position
            + direction
            * bond_length_nm
        )

    missing_positions = [
        node_id
        for node_id in nodes_by_id
        if node_id not in positions
    ]

    if missing_positions:
        raise RuntimeError(
            "Coordinate assignment incomplete: "
            + " | ".join(
                missing_positions[:20]
            )
        )

    if len(positions) != EXPECTED_TOTAL_COORDINATE_NODES:
        raise RuntimeError(
            "Unexpected coordinate-node count: "
            f"{len(positions)}/"
            f"{EXPECTED_TOTAL_COORDINATE_NODES}"
        )

    coordinate_rows = []

    for node_id in sorted(
        positions
    ):
        position = positions[
            node_id
        ]

        if not np.all(
            np.isfinite(
                position
            )
        ):
            raise RuntimeError(
                f"Non-finite coordinate for {node_id}"
            )

        source = nodes_by_id[
            node_id
        ]

        coordinate_rows.append(
            {
                "node_id": node_id,
                "element": source[
                    "element"
                ],
                "node_type": source[
                    "node_type"
                ],
                "end": source[
                    "end"
                ],
                "x_nm": float(
                    position[0]
                ),
                "y_nm": float(
                    position[1]
                ),
                "z_nm": float(
                    position[2]
                ),
                "coordinates_generated_by": (
                    "STATIC_ANALYTIC_EMBEDDING"
                ),
                "energy_minimized": False,
                "MD_relaxed": False,
            }
        )

    write_csv(
        COORDINATES_CSV,
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

        length_nm = float(
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
            target_nm = BN_TARGET_NM
        else:
            heavy_node = (
                first
                if nodes_by_id[
                    first
                ][
                    "element"
                ]
                != "H"
                else second
            )

            target_nm = xh_target_nm(
                nodes_by_id[
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
                "source_element": nodes_by_id[
                    first
                ][
                    "element"
                ],
                "target_element": nodes_by_id[
                    second
                ][
                    "element"
                ],
                "length_nm": length_nm,
                "target_length_nm": (
                    target_nm
                ),
                "deviation_nm": (
                    length_nm
                    - target_nm
                ),
                "absolute_deviation_nm": abs(
                    length_nm
                    - target_nm
                ),
            }
        )

    write_csv(
        BOND_LENGTHS_CSV,
        bond_rows,
    )

    bond_type_rows = (
        edge_type_statistics(
            bond_rows
        )
    )

    write_csv(
        BOND_TYPE_SUMMARY_CSV,
        bond_type_rows,
    )

    bond_type_by_name = {
        row[
            "edge_type"
        ]: row
        for row in bond_type_rows
    }

    parent_terminal_nodes = {
        node_id
        for node_id in nodes_by_id
        if node_id.startswith(
            "P:"
        )
        and any(
            neighbor.startswith(
                "S:"
            )
            for neighbor in adjacency[
                node_id
            ]
        )
    }

    junction_nodes = {
        node_id
        for node_id, row
        in nodes_by_id.items()
        if row[
            "node_type"
        ]
        in {
            "HEXAGONAL_EDGE_COMPLETION_SEED",
            "ANNULUS_OUTER_BOUNDARY",
        }
    } | parent_terminal_nodes

    angle_rows = []

    junction_degree_failures = []

    for node_id in sorted(
        junction_nodes
    ):
        neighbors = sorted(
            adjacency[
                node_id
            ]
        )

        if len(neighbors) != 3:
            junction_degree_failures.append(
                node_id
            )
            continue

        angle_index = 0

        for first_index in range(3):
            for second_index in range(
                first_index + 1,
                3,
            ):
                angle_index += 1

                first = neighbors[
                    first_index
                ]

                second = neighbors[
                    second_index
                ]

                angle = angle_degrees(
                    positions[
                        first
                    ],
                    positions[
                        node_id
                    ],
                    positions[
                        second
                    ],
                )

                angle_rows.append(
                    {
                        "center_node": (
                            node_id
                        ),
                        "center_element": (
                            nodes_by_id[
                                node_id
                            ][
                                "element"
                            ]
                        ),
                        "center_node_type": (
                            nodes_by_id[
                                node_id
                            ][
                                "node_type"
                            ]
                        ),
                        "end": nodes_by_id[
                            node_id
                        ][
                            "end"
                        ],
                        "angle_index": (
                            angle_index
                        ),
                        "neighbor_1": (
                            first
                        ),
                        "neighbor_2": (
                            second
                        ),
                        "angle_deg": (
                            angle
                        ),
                        "deviation_from_120_deg": (
                            angle - 120.0
                        ),
                    }
                )

    if not angle_rows:
        raise RuntimeError(
            "No junction angles were calculated."
        )

    write_csv(
        JUNCTION_ANGLES_CSV,
        angle_rows,
    )

    angle_values = np.asarray(
        [
            float(
                row[
                    "angle_deg"
                ]
            )
            for row in angle_rows
        ],
        dtype=float,
    )

    angle_rms_deviation_deg = float(
        np.sqrt(
            np.mean(
                (
                    angle_values
                    - 120.0
                )
                ** 2
            )
        )
    )

    ordered_node_ids = sorted(
        positions
    )

    node_index = {
        node_id: index
        for index, node_id
        in enumerate(
            ordered_node_ids
        )
    }

    coordinate_array = np.asarray(
        [
            positions[
                node_id
            ]
            for node_id
            in ordered_node_ids
        ],
        dtype=float,
    )

    bonded_integer_pairs = {
        tuple(
            sorted(
                (
                    node_index[first],
                    node_index[second],
                )
            )
        )
        for first, second in edge_pairs
    }

    minimum_nonbonded = {
        "HEAVY_HEAVY": math.inf,
        "H_HEAVY": math.inf,
        "H_H": math.inf,
    }

    minimum_nonbonded_pair = {
        "HEAVY_HEAVY": "",
        "H_HEAVY": "",
        "H_H": "",
    }

    clash_counts = {
        "HEAVY_HEAVY": 0,
        "H_HEAVY": 0,
        "H_H": 0,
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

    for first_index in range(
        len(
            ordered_node_ids
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

        first_id = ordered_node_ids[
            first_index
        ]

        first_is_h = (
            nodes_by_id[
                first_id
            ][
                "element"
            ]
            == "H"
        )

        for local_index, distance_nm in enumerate(
            distances
        ):
            second_index = (
                first_index
                + 1
                + local_index
            )

            if (
                first_index,
                second_index,
            ) in bonded_integer_pairs:
                continue

            second_id = ordered_node_ids[
                second_index
            ]

            second_is_h = (
                nodes_by_id[
                    second_id
                ][
                    "element"
                ]
                == "H"
            )

            if (
                first_is_h
                and second_is_h
            ):
                category = "H_H"
            elif (
                first_is_h
                or second_is_h
            ):
                category = "H_HEAVY"
            else:
                category = "HEAVY_HEAVY"

            distance_value = float(
                distance_nm
            )

            if (
                distance_value
                < minimum_nonbonded[
                    category
                ]
            ):
                minimum_nonbonded[
                    category
                ] = distance_value

                minimum_nonbonded_pair[
                    category
                ] = (
                    f"{first_id} | {second_id}"
                )

            if (
                distance_value
                < thresholds[
                    category
                ]
            ):
                clash_counts[
                    category
                ] += 1

    for category in minimum_nonbonded:
        if not math.isfinite(
            minimum_nonbonded[
                category
            ]
        ):
            raise RuntimeError(
                f"No nonbonded pairs found for {category}"
            )

    for end_row in embedding_end_rows:
        end = end_row[
            "end"
        ]

        annulus_center = (
            annulus_centers[
                end
            ]
        )

        basis_x, basis_y, outward = (
            annulus_plane_axes[
                end
            ]
        )

        inner_h_ids = [
            node_id
            for node_id, row
            in nodes_by_id.items()
            if row[
                "end"
            ]
            == end
            and row[
                "node_type"
            ]
            == "ANNULUS_INNER_PASSIVANT_H"
        ]

        outer_annulus_ids = [
            node_id
            for node_id, row
            in nodes_by_id.items()
            if row[
                "end"
            ]
            == end
            and row[
                "node_type"
            ]
            == "ANNULUS_OUTER_BOUNDARY"
        ]

        annulus_ids = [
            node_id
            for node_id, row
            in nodes_by_id.items()
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

        inner_h_radii = []

        for node_id in inner_h_ids:
            displacement = (
                positions[
                    node_id
                ]
                - annulus_center
            )

            displacement = (
                displacement
                - np.dot(
                    displacement,
                    outward,
                )
                * outward
            )

            inner_h_radii.append(
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

            displacement = (
                displacement
                - np.dot(
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

        annulus_axial_values = [
            float(
                np.dot(
                    positions[
                        node_id
                    ]
                    - annulus_center,
                    outward,
                )
            )
            for node_id in annulus_ids
        ]

        nuclear_aperture_diameter_nm = (
            2.0
            * min(
                inner_h_radii
            )
        )

        end_row[
            "inner_H_atoms"
        ] = len(
            inner_h_ids
        )

        end_row[
            "minimum_inner_H_radius_nm"
        ] = min(
            inner_h_radii
        )

        end_row[
            "mean_inner_H_radius_nm"
        ] = float(
            np.mean(
                inner_h_radii
            )
        )

        end_row[
            "nuclear_aperture_diameter_nm"
        ] = (
            nuclear_aperture_diameter_nm
        )

        end_row[
            "nuclear_aperture_relative_error"
        ] = abs(
            nuclear_aperture_diameter_nm
            - target_aperture_diameter_nm
        ) / target_aperture_diameter_nm

        end_row[
            "outer_annulus_radius_mean_nm"
        ] = float(
            np.mean(
                outer_radii
            )
        )

        end_row[
            "outer_annulus_radius_relative_error"
        ] = abs(
            float(
                np.mean(
                    outer_radii
                )
            )
            - target_parent_rim_radius_nm
        ) / target_parent_rim_radius_nm

        end_row[
            "annulus_plane_axial_std_nm"
        ] = float(
            np.std(
                annulus_axial_values
            )
        )

    write_csv(
        END_SUMMARY_CSV,
        embedding_end_rows,
    )

    lower = next(
        row
        for row in embedding_end_rows
        if row[
            "end"
        ]
        == "LOWER"
    )

    upper = next(
        row
        for row in embedding_end_rows
        if row[
            "end"
        ]
        == "UPPER"
    )

    attachment_stats = (
        bond_type_by_name[
            "PARTIAL_HETEROPOLAR_SEED_TO_ANNULUS"
        ]
    )

    parent_seed_stats = (
        bond_type_by_name[
            "PARENT_TO_COMPLETION_SEED"
        ]
    )

    annulus_stats = (
        bond_type_by_name[
            "ANNULUS_BN"
        ]
    )

    h_edge_types = (
        "SEED_H_PASSIVATION",
        "ANNULUS_OUTER_H_PASSIVATION",
        "ANNULUS_INNER_H_PASSIVATION",
    )

    maximum_h_bond_deviation_nm = max(
        float(
            bond_type_by_name[
                edge_type
            ][
                "maximum_absolute_deviation_nm"
            ]
        )
        for edge_type in h_edge_types
    )

    end_asymmetry_metrics = {
        "attachment_plane_gap_nm": abs(
            float(
                lower[
                    "attachment_plane_gap_nm"
                ]
            )
            - float(
                upper[
                    "attachment_plane_gap_nm"
                ]
            )
        ),
        "attachment_RMS_deviation_nm": abs(
            float(
                lower[
                    "predicted_attachment_RMS_deviation_nm"
                ]
            )
            - float(
                upper[
                    "predicted_attachment_RMS_deviation_nm"
                ]
            )
        ),
        "nuclear_aperture_diameter_nm": abs(
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
        "outer_annulus_radius_mean_nm": abs(
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
    }

    maximum_end_asymmetry_nm = max(
        end_asymmetry_metrics.values()
    )

    gates = {
        "Gate3A_parent_audit_is_accepted": (
            parent_summary.get(
                "decision"
            )
            == EXPECTED_PARENT_DECISION
        ),
        "Gate3F_graph_design_is_accepted": (
            design_summary.get(
                "decision"
            )
            == EXPECTED_DESIGN_DECISION
        ),
        "Gate3F_has_no_failed_gates": (
            len(
                failed_design_gates
            )
            == 0
        ),
        "all_2076_nodes_received_finite_coordinates": (
            len(
                positions
            )
            == EXPECTED_TOTAL_COORDINATE_NODES
            and all(
                np.all(
                    np.isfinite(
                        position
                    )
                )
                for position
                in positions.values()
            )
        ),
        "parent_coordinates_are_unchanged": all(
            np.array_equal(
                positions[
                    f"P:{index}"
                ],
                parent_positions[
                    index
                ],
            )
            for index in range(
                EXPECTED_PARENT_ATOMS
            )
        ),
        "parent_seed_bonds_match_target_within_0p003nm": (
            float(
                parent_seed_stats[
                    "maximum_absolute_deviation_nm"
                ]
            )
            <= MAX_PARENT_SEED_BOND_DEVIATION_NM
        ),
        "annulus_internal_bonds_match_target_within_0p002nm": (
            float(
                annulus_stats[
                    "maximum_absolute_deviation_nm"
                ]
            )
            <= MAX_ANNULUS_BOND_DEVIATION_NM
        ),
        "seed_annulus_attachment_RMS_deviation_within_0p020nm": (
            float(
                attachment_stats[
                    "RMS_deviation_nm"
                ]
            )
            <= MAX_ATTACHMENT_BOND_RMS_DEVIATION_NM
        ),
        "seed_annulus_attachment_max_deviation_within_0p035nm": (
            float(
                attachment_stats[
                    "maximum_absolute_deviation_nm"
                ]
            )
            <= MAX_ATTACHMENT_BOND_MAX_DEVIATION_NM
        ),
        "all_XH_bonds_match_provisional_targets_within_0p002nm": (
            maximum_h_bond_deviation_nm
            <= MAX_XH_BOND_DEVIATION_NM
        ),
        "junction_nodes_all_have_three_neighbors": (
            len(
                junction_degree_failures
            )
            == 0
        ),
        "junction_angle_minimum_is_at_least_70deg": (
            float(
                np.min(
                    angle_values
                )
            )
            >= MIN_JUNCTION_ANGLE_DEG
        ),
        "junction_angle_maximum_is_at_most_170deg": (
            float(
                np.max(
                    angle_values
                )
            )
            <= MAX_JUNCTION_ANGLE_DEG
        ),
        "junction_angle_RMS_deviation_is_at_most_35deg": (
            angle_rms_deviation_deg
            <= MAX_JUNCTION_ANGLE_RMS_DEVIATION_DEG
        ),
        "annulus_center_offset_is_at_most_0p050nm": all(
            float(
                row[
                    "annulus_center_offset_nm"
                ]
            )
            <= MAX_ANNULUS_CENTER_OFFSET_NM
            for row in embedding_end_rows
        ),
        "attachment_plane_gap_is_positive_and_at_most_0p250nm": all(
            0.0
            < float(
                row[
                    "attachment_plane_gap_nm"
                ]
            )
            <= MAX_ATTACHMENT_PLANE_GAP_NM
            for row in embedding_end_rows
        ),
        "annulus_is_planar_to_numerical_precision": all(
            float(
                row[
                    "annulus_plane_axial_std_nm"
                ]
            )
            <= 1.0e-10
            for row in embedding_end_rows
        ),
        "nuclear_aperture_error_is_within10percent": all(
            float(
                row[
                    "nuclear_aperture_relative_error"
                ]
            )
            <= MAX_NUCLEAR_APERTURE_RELATIVE_ERROR
            for row in embedding_end_rows
        ),
        "outer_radius_error_is_within15percent": all(
            float(
                row[
                    "outer_annulus_radius_relative_error"
                ]
            )
            <= MAX_OUTER_RADIUS_RELATIVE_ERROR
            for row in embedding_end_rows
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
        "lower_and_upper_embeddings_are_symmetric": (
            maximum_end_asymmetry_nm
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

    decision = (
        PASS_DECISION
        if accepted
        else REVIEW_DECISION
    )

    required_next_step = (
        "AUDIT_R2_PARTIAL_ATTACHMENT_CHEMICAL_REALIZABILITY_"
        "AND_PARAMETERIZATION_SCOPE"
        if accepted
        else
        "OPTIMIZE_R2_PARTIAL_ATTACHMENT_STATIC_EMBEDDING_"
        "WITH_CONSTRAINED_GEOMETRY"
    )

    summary = {
        "decision": decision,
        "coordinate_nodes": len(
            positions
        ),
        "parent_atoms": (
            EXPECTED_PARENT_ATOMS
        ),
        "added_heavy_atoms": (
            EXPECTED_ADDED_HEAVY_ATOMS
        ),
        "added_H_atoms": (
            EXPECTED_ADDED_H_ATOMS
        ),
        "tube_axis_x": float(
            tube_axis[0]
        ),
        "tube_axis_y": float(
            tube_axis[1]
        ),
        "tube_axis_z": float(
            tube_axis[2]
        ),
        "box_x_nm": float(
            box[0]
        ),
        "box_y_nm": float(
            box[1]
        ),
        "box_z_nm": float(
            box[2]
        ),
        "parent_seed_bond_RMS_deviation_nm": (
            parent_seed_stats[
                "RMS_deviation_nm"
            ]
        ),
        "parent_seed_bond_max_deviation_nm": (
            parent_seed_stats[
                "maximum_absolute_deviation_nm"
            ]
        ),
        "attachment_bond_mean_nm": (
            attachment_stats[
                "mean_length_nm"
            ]
        ),
        "attachment_bond_minimum_nm": (
            attachment_stats[
                "minimum_length_nm"
            ]
        ),
        "attachment_bond_maximum_nm": (
            attachment_stats[
                "maximum_length_nm"
            ]
        ),
        "attachment_bond_RMS_deviation_nm": (
            attachment_stats[
                "RMS_deviation_nm"
            ]
        ),
        "attachment_bond_max_deviation_nm": (
            attachment_stats[
                "maximum_absolute_deviation_nm"
            ]
        ),
        "junction_angle_minimum_deg": float(
            np.min(
                angle_values
            )
        ),
        "junction_angle_maximum_deg": float(
            np.max(
                angle_values
            )
        ),
        "junction_angle_mean_deg": float(
            np.mean(
                angle_values
            )
        ),
        "junction_angle_RMS_deviation_deg": (
            angle_rms_deviation_deg
        ),
        "lower_attachment_plane_gap_nm": (
            lower[
                "attachment_plane_gap_nm"
            ]
        ),
        "upper_attachment_plane_gap_nm": (
            upper[
                "attachment_plane_gap_nm"
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
        "lower_nuclear_aperture_relative_error": (
            lower[
                "nuclear_aperture_relative_error"
            ]
        ),
        "upper_nuclear_aperture_relative_error": (
            upper[
                "nuclear_aperture_relative_error"
            ]
        ),
        "minimum_nonbonded_heavy_heavy_nm": (
            minimum_nonbonded[
                "HEAVY_HEAVY"
            ]
        ),
        "minimum_nonbonded_heavy_heavy_pair": (
            minimum_nonbonded_pair[
                "HEAVY_HEAVY"
            ]
        ),
        "minimum_nonbonded_H_heavy_nm": (
            minimum_nonbonded[
                "H_HEAVY"
            ]
        ),
        "minimum_nonbonded_H_heavy_pair": (
            minimum_nonbonded_pair[
                "H_HEAVY"
            ]
        ),
        "minimum_nonbonded_H_H_nm": (
            minimum_nonbonded[
                "H_H"
            ]
        ),
        "minimum_nonbonded_H_H_pair": (
            minimum_nonbonded_pair[
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
        "maximum_lower_upper_asymmetry_nm": (
            maximum_end_asymmetry_nm
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

    EMBEDDING_JSON.write_text(
        json.dumps(
            {
                "summary": summary,
                "end_summaries": (
                    embedding_end_rows
                ),
                "bond_type_summaries": (
                    bond_type_rows
                ),
                "gates": gates,
                "provisional_geometry_targets_nm": {
                    "B_N": (
                        BN_TARGET_NM
                    ),
                    "B_H": (
                        BH_TARGET_NM
                    ),
                    "N_H": (
                        NH_TARGET_NM
                    ),
                },
                "limitations": [
                    (
                        "This is an analytic static coordinate "
                        "embedding, not an energy-minimized structure."
                    ),
                    (
                        "B-H and N-H lengths are provisional "
                        "geometry targets, not force-field parameters."
                    ),
                    (
                        "The nuclear aperture does not equal a "
                        "hydrated free-energy or 5-kBT aperture."
                    ),
                    (
                        "Passing this gate does not prove energetic "
                        "stability or synthetic realizability."
                    ),
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    source_rows = [
        {
            "role": (
                "accepted_parent_GRO"
            ),
            "file": relative(
                SYSTEM_GRO
            ),
            "sha256": sha256(
                SYSTEM_GRO
            ),
        },
        {
            "role": (
                "Gate3A_parent_summary"
            ),
            "file": relative(
                PARENT_SUMMARY_CSV
            ),
            "sha256": sha256(
                PARENT_SUMMARY_CSV
            ),
        },
        {
            "role": (
                "Gate3A_terminal_atoms"
            ),
            "file": relative(
                TERMINAL_ATOMS_CSV
            ),
            "sha256": sha256(
                TERMINAL_ATOMS_CSV
            ),
        },
        {
            "role": (
                "Gate3F_design_nodes"
            ),
            "file": relative(
                DESIGN_NODES_CSV
            ),
            "sha256": sha256(
                DESIGN_NODES_CSV
            ),
        },
        {
            "role": (
                "Gate3F_design_edges"
            ),
            "file": relative(
                DESIGN_EDGES_CSV
            ),
            "sha256": sha256(
                DESIGN_EDGES_CSV
            ),
        },
        {
            "role": (
                "Gate3F_design_summary"
            ),
            "file": relative(
                DESIGN_SUMMARY_CSV
            ),
            "sha256": sha256(
                DESIGN_SUMMARY_CSV
            ),
        },
    ]

    write_csv(
        SOURCE_MANIFEST_CSV,
        source_rows,
    )

    xyz_lines = [
        str(
            len(
                ordered_node_ids
            )
        ),
        (
            "R2 partial-attachment annulus static embedding; "
            "not energy minimized"
        ),
    ]

    for node_id in ordered_node_ids:
        position_angstrom = (
            positions[
                node_id
            ]
            * 10.0
        )

        xyz_lines.append(
            f"{nodes_by_id[node_id]['element']:2s} "
            f"{position_angstrom[0]: .8f} "
            f"{position_angstrom[1]: .8f} "
            f"{position_angstrom[2]: .8f}"
        )

    XYZ_FILE.write_text(
        "\n".join(
            xyz_lines
        )
        + "\n",
        encoding="utf-8",
    )

    pdb_lines = [
        "REMARK R2 PARTIAL-ATTACHMENT STATIC EMBEDDING",
        "REMARK NOT ENERGY MINIMIZED; NO FORCE-FIELD TOPOLOGY",
    ]

    for serial, node_id in enumerate(
        ordered_node_ids,
        start=1,
    ):
        row = nodes_by_id[
            node_id
        ]

        position_angstrom = (
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
                serial % 1000
            )
        )[:4]

        pdb_lines.append(
            f"ATOM  {serial:5d} "
            f"{atom_name:>4s} "
            f"{residue:>3s} "
            f"{chain:1s}"
            f"{1:4d}    "
            f"{position_angstrom[0]:8.3f}"
            f"{position_angstrom[1]:8.3f}"
            f"{position_angstrom[2]:8.3f}"
            f"{1.00:6.2f}"
            f"{0.00:6.2f}          "
            f"{element:>2s}"
        )

    pdb_lines.append(
        "END"
    )

    PDB_FILE.write_text(
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

    REPORT_MD.write_text(
        f"""# R2 Partial-Attachment Annulus Static Coordinate Embedding

## Scope

This stage generated an analytic, non-minimized coordinate embedding
for the Gate 3F graph.

No molecular topology, formal charges, force-field parameters,
minimization, MD, or QM calculation was generated.

## Coordinate inventory

- Parent HBN atoms:
  **{EXPECTED_PARENT_ATOMS}**
- Added B/N atoms:
  **{EXPECTED_ADDED_HEAVY_ATOMS}**
- Added H atoms:
  **{EXPECTED_ADDED_H_ATOMS}**
- Total coordinate nodes:
  **{len(positions)}**
- Parent coordinates modified:
  **NO**

## Attachment geometry

- Seed–annulus mean/minimum/maximum bond:
  **{float(attachment_stats['mean_length_nm']):.6f}/
  {float(attachment_stats['minimum_length_nm']):.6f}/
  {float(attachment_stats['maximum_length_nm']):.6f} nm**
- Seed–annulus RMS/max target deviation:
  **{float(attachment_stats['RMS_deviation_nm']):.6f}/
  {float(attachment_stats['maximum_absolute_deviation_nm']):.6f} nm**
- Lower/upper plane gap:
  **{float(lower['attachment_plane_gap_nm']):.6f}/
  {float(upper['attachment_plane_gap_nm']):.6f} nm**
- Lower/upper annulus-center offset:
  **{float(lower['annulus_center_offset_nm']):.6f}/
  {float(upper['annulus_center_offset_nm']):.6f} nm**

## Junction angles

- Minimum/mean/maximum:
  **{float(np.min(angle_values)):.3f}/
  {float(np.mean(angle_values)):.3f}/
  {float(np.max(angle_values)):.3f} degrees**
- RMS deviation from 120 degrees:
  **{angle_rms_deviation_deg:.3f} degrees**
- Junction degree failures:
  **{len(junction_degree_failures)}**

## Aperture proxy

- Target diameter:
  **{target_aperture_diameter_nm:.6f} nm**
- Lower nuclear H-defined diameter:
  **{float(lower['nuclear_aperture_diameter_nm']):.6f} nm**
- Upper nuclear H-defined diameter:
  **{float(upper['nuclear_aperture_diameter_nm']):.6f} nm**
- Lower/upper relative error:
  **{float(lower['nuclear_aperture_relative_error']):.6f}/
  {float(upper['nuclear_aperture_relative_error']):.6f}**

This is a nucleus-position proxy. It is not a hydrated free-energy
aperture and cannot replace an excluded-volume or PMF calculation.

## Nonbonded contacts

- Heavy–heavy minimum/count below threshold:
  **{minimum_nonbonded['HEAVY_HEAVY']:.6f}/
  {clash_counts['HEAVY_HEAVY']}**
- H–heavy minimum/count below threshold:
  **{minimum_nonbonded['H_HEAVY']:.6f}/
  {clash_counts['H_HEAVY']}**
- H–H minimum/count below threshold:
  **{minimum_nonbonded['H_H']:.6f}/
  {clash_counts['H_H']}**

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
        "Day024 R2 partial-attachment annulus static "
        "coordinate embedding completed."
    )

    print(
        "Coordinate nodes parent/added-heavy/H/total: "
        f"{EXPECTED_PARENT_ATOMS}/"
        f"{EXPECTED_ADDED_HEAVY_ATOMS}/"
        f"{EXPECTED_ADDED_H_ATOMS}/"
        f"{len(positions)}"
    )

    print(
        "Parent coordinates unchanged: YES"
    )

    print(
        "Parent-seed bond RMS/max deviation: "
        f"{float(parent_seed_stats['RMS_deviation_nm']):.6f}/"
        f"{float(parent_seed_stats['maximum_absolute_deviation_nm']):.6f} nm"
    )

    print(
        "Annulus BN bond RMS/max deviation: "
        f"{float(annulus_stats['RMS_deviation_nm']):.6f}/"
        f"{float(annulus_stats['maximum_absolute_deviation_nm']):.6f} nm"
    )

    print(
        "Seed-annulus bond mean/min/max: "
        f"{float(attachment_stats['mean_length_nm']):.6f}/"
        f"{float(attachment_stats['minimum_length_nm']):.6f}/"
        f"{float(attachment_stats['maximum_length_nm']):.6f} nm"
    )

    print(
        "Seed-annulus RMS/max target deviation: "
        f"{float(attachment_stats['RMS_deviation_nm']):.6f}/"
        f"{float(attachment_stats['maximum_absolute_deviation_nm']):.6f} nm"
    )

    print(
        "Junction angles min/mean/max/RMSdev120: "
        f"{float(np.min(angle_values)):.3f}/"
        f"{float(np.mean(angle_values)):.3f}/"
        f"{float(np.max(angle_values)):.3f}/"
        f"{angle_rms_deviation_deg:.3f} deg"
    )

    for row in embedding_end_rows:
        print(
            f"{row['end']} center-offset/gap/"
            "nuclear-aperture/error/outer-radius/error: "
            f"{float(row['annulus_center_offset_nm']):.6f}/"
            f"{float(row['attachment_plane_gap_nm']):.6f}/"
            f"{float(row['nuclear_aperture_diameter_nm']):.6f}/"
            f"{float(row['nuclear_aperture_relative_error']):.6f}/"
            f"{float(row['outer_annulus_radius_mean_nm']):.6f}/"
            f"{float(row['outer_annulus_radius_relative_error']):.6f}"
        )

    print(
        "Minimum nonbonded heavy-heavy / H-heavy / H-H: "
        f"{minimum_nonbonded['HEAVY_HEAVY']:.6f}/"
        f"{minimum_nonbonded['H_HEAVY']:.6f}/"
        f"{minimum_nonbonded['H_H']:.6f} nm"
    )

    print(
        "Clash counts heavy-heavy / H-heavy / H-H: "
        f"{clash_counts['HEAVY_HEAVY']}/"
        f"{clash_counts['H_HEAVY']}/"
        f"{clash_counts['H_H']}"
    )

    print(
        "Maximum lower-upper asymmetry: "
        f"{maximum_end_asymmetry_nm:.12f} nm"
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
        COORDINATES_CSV,
        BOND_LENGTHS_CSV,
        BOND_TYPE_SUMMARY_CSV,
        JUNCTION_ANGLES_CSV,
        END_SUMMARY_CSV,
        SUMMARY_CSV,
        GATES_CSV,
        EMBEDDING_JSON,
        SOURCE_MANIFEST_CSV,
        XYZ_FILE,
        PDB_FILE,
        REPORT_MD,
    ):
        print(
            f"Wrote: {relative(path)}"
        )


if __name__ == "__main__":
    main()
