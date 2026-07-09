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

GATE3F_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "07_r2_reconstruction_vs_partial_attachment_contingency"
)

GATE3G_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "08_r2_partial_attachment_annulus_static_coordinate_embedding"
)

GATE3G1_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "09_r2_direct_junction_geometric_lower_bound"
)

OUTPUT_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "10_r2_alternating_bn_oligomer_bridge_feasibility"
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

STATIC_COORDINATES_CSV = (
    GATE3G_ROOT
    / "r2_partial_attachment_static_coordinates.csv"
)

STATIC_EMBEDDING_SUMMARY_CSV = (
    GATE3G_ROOT
    / "r2_partial_attachment_static_embedding_summary.csv"
)

DIRECT_LOWER_BOUND_SUMMARY_CSV = (
    GATE3G1_ROOT
    / "r2_direct_junction_geometric_lower_bound_summary.csv"
)

CONFORMER_ENVELOPES_CSV = (
    OUTPUT_ROOT
    / "r2_bn_oligomer_bridge_conformer_envelopes.csv"
)

MAPPING_SCREEN_CSV = (
    OUTPUT_ROOT
    / "r2_bn_oligomer_bridge_mapping_screen.csv"
)

CLASS_SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_bn_oligomer_bridge_class_summary.csv"
)

SELECTED_CANDIDATE_CSV = (
    OUTPUT_ROOT
    / "r2_bn_oligomer_bridge_selected_candidate.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_bn_oligomer_bridge_feasibility_summary.csv"
)

GATES_CSV = (
    OUTPUT_ROOT
    / "r2_bn_oligomer_bridge_feasibility_gates.csv"
)

AUDIT_JSON = (
    OUTPUT_ROOT
    / "r2_bn_oligomer_bridge_feasibility.json"
)

SOURCE_MANIFEST_CSV = (
    OUTPUT_ROOT
    / "r2_bn_oligomer_bridge_source_manifest.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_ALTERNATING_BN_OLIGOMER_BRIDGE_FEASIBILITY_DAY024.md"
)

EXPECTED_DESIGN_DECISION = (
    "R2_PARTIAL_HETEROPOLAR_ANNULUS_ATTACHMENT_AND_"
    "COMPLEMENTARY_PASSIVATION_GRAPH_VALIDATED"
)

EXPECTED_STATIC_EMBEDDING_DECISION = (
    "R2_PARTIAL_ATTACHMENT_ANNULUS_STATIC_COORDINATE_EMBEDDING_"
    "REQUIRES_CONSTRAINED_GEOMETRIC_OPTIMIZATION"
)

EXPECTED_DIRECT_LOWER_BOUND_DECISION = (
    "R2_PARTIAL_ATTACHMENT_DIRECT_BN_JUNCTION_GEOMETRIC_"
    "LOWER_BOUND_FAILED"
)

PASS_DECISION = (
    "R2_SHORTEST_ALTERNATING_BN_OLIGOMER_BRIDGE_CLASS_IDENTIFIED"
)

FAIL_DECISION = (
    "R2_ALTERNATING_BN_OLIGOMER_BRIDGE_CLASSES_1_TO_3_"
    "GEOMETRICALLY_INFEASIBLE"
)

BN_TARGET_NM = 0.144973

MIN_INTERNAL_ANGLE_DEG = 105
MAX_INTERNAL_ANGLE_DEG = 135
INTERNAL_ANGLE_STEP_DEG = 5

TORSION_STEP_DEG = 30

MIN_BRIDGE_ATOMS = 1
MAX_BRIDGE_ATOMS = 3

MIN_USEFUL_AXIAL_GAP_NM = 0.050
MAX_USEFUL_AXIAL_GAP_NM = 0.350
AXIAL_GAP_GRID_POINTS = 3001

MAX_ANNULUS_CENTER_OFFSET_NM = 0.050

EXPECTED_SEED_SITES_PER_END = 30
EXPECTED_SELECTED_ATTACHMENTS_PER_END = 15
EXPECTED_ANNULUS_OUTER_SITES_PER_ELEMENT = 15

EXPECTED_FITS_PER_BRIDGE_CLASS_PER_END = 120

EXPECTED_TOTAL_MAPPING_FITS = (
    2
    * (
        MAX_BRIDGE_ATOMS
        - MIN_BRIDGE_ATOMS
        + 1
    )
    * EXPECTED_FITS_PER_BRIDGE_CLASS_PER_END
)

BASE_ADDED_HEAVY_ATOMS_PER_END = 156
BASE_PASSIVANTS_PER_END = 42

TARGET_HEAVY_ATOMS_PER_END = 145.133


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

    covariance = (
        centered.T
        @ centered
    )

    eigenvalues, eigenvectors = np.linalg.eigh(
        covariance
    )

    axis = eigenvectors[
        :,
        int(
            np.argmax(
                eigenvalues
            )
        ),
    ]

    axis = normalized(
        axis
    )

    if axis[2] < 0.0:
        axis = -axis

    return center, axis


def opposite_element(
    element: str,
) -> str:
    if element == "B":
        return "N"

    if element == "N":
        return "B"

    raise RuntimeError(
        f"Unexpected BN element: {element}"
    )


def required_annulus_element(
    seed_element: str,
    bridge_atoms: int,
) -> str:
    total_path_edges = (
        bridge_atoms + 1
    )

    if total_path_edges % 2 == 0:
        return seed_element

    return opposite_element(
        seed_element
    )


def bridge_element_sequence(
    seed_element: str,
    bridge_atoms: int,
) -> list[str]:
    sequence = []
    current = seed_element

    for _ in range(
        bridge_atoms
    ):
        current = opposite_element(
            current
        )

        sequence.append(
            current
        )

    return sequence


def local_xy(
    lattice_x: int,
    lattice_y: int,
) -> np.ndarray:
    return np.asarray(
        [
            lattice_x
            * BN_TARGET_NM
            / 2.0,
            lattice_y
            * math.sqrt(3.0)
            * BN_TARGET_NM
            / 2.0,
        ],
        dtype=float,
    )


def fit_orthogonal(
    source: np.ndarray,
    target: np.ndarray,
    desired_determinant: int,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:
    source_center = np.mean(
        source,
        axis=0,
    )

    target_center = np.mean(
        target,
        axis=0,
    )

    source_centered = (
        source
        - source_center
    )

    target_centered = (
        target
        - target_center
    )

    covariance = (
        source_centered.T
        @ target_centered
    )

    left, _, right_t = np.linalg.svd(
        covariance
    )

    rotation = (
        left
        @ right_t
    )

    observed_determinant = int(
        round(
            np.linalg.det(
                rotation
            )
        )
    )

    if observed_determinant != desired_determinant:
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

    return rotation, translation


def sample_chain_end_to_end_distances(
    bridge_atoms: int,
) -> np.ndarray:
    number_of_bonds = (
        bridge_atoms + 1
    )

    theta_values = np.deg2rad(
        np.arange(
            MIN_INTERNAL_ANGLE_DEG,
            MAX_INTERNAL_ANGLE_DEG + 1,
            INTERNAL_ANGLE_STEP_DEG,
        )
    )

    phi_values = np.deg2rad(
        np.arange(
            0,
            360,
            TORSION_STEP_DEG,
        )
    )

    distances: list[float] = []

    first_bond = np.asarray(
        [
            BN_TARGET_NM,
            0.0,
            0.0,
        ],
        dtype=float,
    )

    def recurse(
        position: np.ndarray,
        previous_bond: np.ndarray,
        completed_bonds: int,
    ) -> None:
        if completed_bonds == number_of_bonds:
            distances.append(
                float(
                    np.linalg.norm(
                        position
                    )
                )
            )

            return

        previous_unit = normalized(
            previous_bond
        )

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
                    previous_unit,
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

        perpendicular_1 = normalized(
            np.cross(
                previous_unit,
                reference,
            )
        )

        perpendicular_2 = normalized(
            np.cross(
                previous_unit,
                perpendicular_1,
            )
        )

        torsions = (
            [0.0]
            if completed_bonds == 1
            else phi_values
        )

        for theta in theta_values:
            forward_angle = (
                math.pi
                - float(theta)
            )

            for phi in torsions:
                direction = (
                    math.cos(
                        forward_angle
                    )
                    * previous_unit
                    + math.sin(
                        forward_angle
                    )
                    * (
                        math.cos(
                            float(phi)
                        )
                        * perpendicular_1
                        + math.sin(
                            float(phi)
                        )
                        * perpendicular_2
                    )
                )

                new_bond = (
                    BN_TARGET_NM
                    * direction
                )

                recurse(
                    position
                    + new_bond,
                    new_bond,
                    completed_bonds + 1,
                )

    recurse(
        first_bond,
        first_bond,
        1,
    )

    values = np.asarray(
        distances,
        dtype=float,
    )

    if (
        values.size == 0
        or not np.all(
            np.isfinite(
                values
            )
        )
    ):
        raise RuntimeError(
            f"Could not sample bridge class {bridge_atoms}."
        )

    return values


def optimize_gap_interval(
    lateral_distances_nm: np.ndarray,
    seed_axial_nm: np.ndarray,
    minimum_span_nm: float,
    maximum_span_nm: float,
) -> dict[str, Any]:
    gap_grid = np.linspace(
        MIN_USEFUL_AXIAL_GAP_NM,
        MAX_USEFUL_AXIAL_GAP_NM,
        AXIAL_GAP_GRID_POINTS,
    )

    mean_seed_axial_nm = float(
        np.mean(
            seed_axial_nm
        )
    )

    annulus_planes_nm = (
        mean_seed_axial_nm
        + gap_grid
    )

    axial_differences_nm = (
        annulus_planes_nm[
            None,
            :,
        ]
        - seed_axial_nm[
            :,
            None,
        ]
    )

    distances_nm = np.sqrt(
        lateral_distances_nm[
            :,
            None,
        ]
        ** 2
        + axial_differences_nm
        ** 2
    )

    feasible_mask = np.all(
        (
            distances_nm
            >= minimum_span_nm
        )
        & (
            distances_nm
            <= maximum_span_nm
        ),
        axis=0,
    )

    feasible_indices = np.flatnonzero(
        feasible_mask
    )

    if feasible_indices.size == 0:
        best_index = int(
            np.argmin(
                np.mean(
                    np.maximum(
                        minimum_span_nm
                        - distances_nm,
                        0.0,
                    )
                    + np.maximum(
                        distances_nm
                        - maximum_span_nm,
                        0.0,
                    ),
                    axis=0,
                )
            )
        )

        selected_distances = (
            distances_nm[
                :,
                best_index
            ]
        )

        violations = (
            np.maximum(
                minimum_span_nm
                - selected_distances,
                0.0,
            )
            + np.maximum(
                selected_distances
                - maximum_span_nm,
                0.0,
            )
        )

        return {
            "gap_feasible": False,
            "feasible_gap_minimum_nm": "",
            "feasible_gap_maximum_nm": "",
            "selected_gap_nm": float(
                gap_grid[
                    best_index
                ]
            ),
            "selected_distance_minimum_nm": float(
                np.min(
                    selected_distances
                )
            ),
            "selected_distance_maximum_nm": float(
                np.max(
                    selected_distances
                )
            ),
            "maximum_envelope_violation_nm": float(
                np.max(
                    violations
                )
            ),
            "mean_envelope_violation_nm": float(
                np.mean(
                    violations
                )
            ),
        }

    first_index = int(
        feasible_indices[0]
    )

    last_index = int(
        feasible_indices[-1]
    )

    selected_index = int(
        feasible_indices[
            feasible_indices.size // 2
        ]
    )

    selected_distances = (
        distances_nm[
            :,
            selected_index
        ]
    )

    return {
        "gap_feasible": True,
        "feasible_gap_minimum_nm": float(
            gap_grid[
                first_index
            ]
        ),
        "feasible_gap_maximum_nm": float(
            gap_grid[
                last_index
            ]
        ),
        "selected_gap_nm": float(
            gap_grid[
                selected_index
            ]
        ),
        "selected_distance_minimum_nm": float(
            np.min(
                selected_distances
            )
        ),
        "selected_distance_maximum_nm": float(
            np.max(
                selected_distances
            )
        ),
        "maximum_envelope_violation_nm": 0.0,
        "mean_envelope_violation_nm": 0.0,
    }


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        DESIGN_NODES_CSV,
        DESIGN_EDGES_CSV,
        DESIGN_SUMMARY_CSV,
        STATIC_COORDINATES_CSV,
        STATIC_EMBEDDING_SUMMARY_CSV,
        DIRECT_LOWER_BOUND_SUMMARY_CSV,
    ):
        require_file(required)

    design_nodes = read_csv_rows(
        DESIGN_NODES_CSV
    )

    design_summary = read_single_csv_row(
        DESIGN_SUMMARY_CSV
    )

    coordinate_rows = read_csv_rows(
        STATIC_COORDINATES_CSV
    )

    static_summary = read_single_csv_row(
        STATIC_EMBEDDING_SUMMARY_CSV
    )

    direct_summary = read_single_csv_row(
        DIRECT_LOWER_BOUND_SUMMARY_CSV
    )

    if design_summary.get(
        "decision"
    ) != EXPECTED_DESIGN_DECISION:
        raise RuntimeError(
            "Gate 3F is not in the accepted state."
        )

    if static_summary.get(
        "decision"
    ) != EXPECTED_STATIC_EMBEDDING_DECISION:
        raise RuntimeError(
            "Gate 3G does not have the expected review decision."
        )

    if direct_summary.get(
        "decision"
    ) != EXPECTED_DIRECT_LOWER_BOUND_DECISION:
        raise RuntimeError(
            "Gate 3G.1 does not contain the expected "
            "direct-junction rejection."
        )

    nodes_by_id = {
        row[
            "node_id"
        ]: row
        for row in design_nodes
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

    if set(nodes_by_id) != set(
        coordinates
    ):
        raise RuntimeError(
            "The Gate 3F nodes and Gate 3G coordinates disagree."
        )

    parent_ids = [
        node_id
        for node_id, row
        in nodes_by_id.items()
        if row[
            "node_type"
        ]
        == "PARENT_HBN"
    ]

    parent_positions = np.asarray(
        [
            coordinates[
                node_id
            ]
            for node_id in parent_ids
        ],
        dtype=float,
    )

    _, tube_axis = determine_axis(
        parent_positions
    )

    conformer_samples: dict[
        int,
        np.ndarray
    ] = {}

    envelope_rows = []

    for bridge_atoms in range(
        MIN_BRIDGE_ATOMS,
        MAX_BRIDGE_ATOMS + 1,
    ):
        samples = (
            sample_chain_end_to_end_distances(
                bridge_atoms
            )
        )

        conformer_samples[
            bridge_atoms
        ] = samples

        percentiles = np.percentile(
            samples,
            [
                1,
                5,
                50,
                95,
                99,
            ],
        )

        envelope_rows.append(
            {
                "bridge_atoms_per_attachment": (
                    bridge_atoms
                ),
                "bond_count_seed_to_annulus": (
                    bridge_atoms + 1
                ),
                "sample_count": int(
                    samples.size
                ),
                "internal_angle_minimum_deg": (
                    MIN_INTERNAL_ANGLE_DEG
                ),
                "internal_angle_maximum_deg": (
                    MAX_INTERNAL_ANGLE_DEG
                ),
                "torsion_step_deg": (
                    TORSION_STEP_DEG
                ),
                "minimum_end_to_end_nm": float(
                    np.min(
                        samples
                    )
                ),
                "percentile_1_nm": float(
                    percentiles[0]
                ),
                "percentile_5_nm": float(
                    percentiles[1]
                ),
                "median_nm": float(
                    percentiles[2]
                ),
                "percentile_95_nm": float(
                    percentiles[3]
                ),
                "percentile_99_nm": float(
                    percentiles[4]
                ),
                "maximum_end_to_end_nm": float(
                    np.max(
                        samples
                    )
                ),
                "maximum_contour_length_nm": (
                    (
                        bridge_atoms + 1
                    )
                    * BN_TARGET_NM
                ),
            }
        )

    write_csv(
        CONFORMER_ENVELOPES_CSV,
        envelope_rows,
    )

    envelope_by_class = {
        int(
            row[
                "bridge_atoms_per_attachment"
            ]
        ): row
        for row in envelope_rows
    }

    mapping_rows = []

    for end in (
        "LOWER",
        "UPPER",
    ):
        outward = (
            -tube_axis
            if end == "LOWER"
            else tube_axis
        )

        seed_rows = [
            row
            for row in design_nodes
            if row[
                "end"
            ]
            == end
            and row[
                "node_type"
            ]
            == "HEXAGONAL_EDGE_COMPLETION_SEED"
        ]

        seed_rows.sort(
            key=lambda row: parse_int(
                row,
                "circumferential_index",
            )
        )

        if len(seed_rows) != EXPECTED_SEED_SITES_PER_END:
            raise RuntimeError(
                f"{end}: unexpected seed count."
            )

        seed_element = seed_rows[0][
            "element"
        ]

        if any(
            row[
                "element"
            ]
            != seed_element
            for row in seed_rows
        ):
            raise RuntimeError(
                f"{end}: seed is not elementally homogeneous."
            )

        seed_positions = np.asarray(
            [
                coordinates[
                    row[
                        "node_id"
                    ]
                ]
                for row in seed_rows
            ],
            dtype=float,
        )

        seed_center = np.mean(
            seed_positions,
            axis=0,
        )

        radial_reference = (
            seed_positions[0]
            - seed_center
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

        all_annulus_rows = [
            row
            for row in design_nodes
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

        all_annulus_xy = np.asarray(
            [
                local_xy(
                    parse_int(
                        row,
                        "lattice_x",
                    ),
                    parse_int(
                        row,
                        "lattice_y",
                    ),
                )
                for row in all_annulus_rows
            ],
            dtype=float,
        )

        for bridge_atoms in range(
            MIN_BRIDGE_ATOMS,
            MAX_BRIDGE_ATOMS + 1,
        ):
            annulus_element = (
                required_annulus_element(
                    seed_element,
                    bridge_atoms,
                )
            )

            bridge_sequence = (
                bridge_element_sequence(
                    seed_element,
                    bridge_atoms,
                )
            )

            envelope = envelope_by_class[
                bridge_atoms
            ]

            minimum_span_nm = float(
                envelope[
                    "minimum_end_to_end_nm"
                ]
            )

            maximum_span_nm = float(
                envelope[
                    "maximum_end_to_end_nm"
                ]
            )

            outer_rows = [
                row
                for row in design_nodes
                if row[
                    "end"
                ]
                == end
                and row[
                    "node_type"
                ]
                == "ANNULUS_OUTER_BOUNDARY"
                and row[
                    "element"
                ]
                == annulus_element
            ]

            outer_rows.sort(
                key=lambda row: parse_float(
                    row,
                    "angle_turns",
                )
            )

            if (
                len(outer_rows)
                != EXPECTED_ANNULUS_OUTER_SITES_PER_ELEMENT
            ):
                raise RuntimeError(
                    f"{end}, bridge={bridge_atoms}: "
                    "unexpected annulus-site count."
                )

            source_base = np.asarray(
                [
                    local_xy(
                        parse_int(
                            row,
                            "lattice_x",
                        ),
                        parse_int(
                            row,
                            "lattice_y",
                        ),
                    )
                    for row in outer_rows
                ],
                dtype=float,
            )

            for seed_parity in (
                0,
                1,
            ):
                selected_seed_rows = [
                    row
                    for row in seed_rows
                    if parse_int(
                        row,
                        "circumferential_index",
                    )
                    % 2
                    == seed_parity
                ]

                if (
                    len(
                        selected_seed_rows
                    )
                    != EXPECTED_SELECTED_ATTACHMENTS_PER_END
                ):
                    raise RuntimeError(
                        f"{end}: parity selection did not "
                        "produce 15 sites."
                    )

                selected_seed_positions = np.asarray(
                    [
                        coordinates[
                            row[
                                "node_id"
                            ]
                        ]
                        for row in selected_seed_rows
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
                        in selected_seed_positions
                    ],
                    dtype=float,
                )

                seed_axial_nm = np.asarray(
                    [
                        float(
                            np.dot(
                                position
                                - seed_center,
                                outward,
                            )
                        )
                        for position
                        in selected_seed_positions
                    ],
                    dtype=float,
                )

                for orientation in (
                    1,
                    -1,
                ):
                    for rotation_index in range(
                        EXPECTED_ANNULUS_OUTER_SITES_PER_ELEMENT
                    ):
                        mapped_indices = [
                            (
                                orientation
                                * index
                                + rotation_index
                            )
                            % EXPECTED_ANNULUS_OUTER_SITES_PER_ELEMENT
                            for index in range(
                                EXPECTED_ANNULUS_OUTER_SITES_PER_ELEMENT
                            )
                        ]

                        source_xy = (
                            source_base[
                                mapped_indices
                            ]
                        )

                        for chirality in (
                            1,
                            -1,
                        ):
                            (
                                rotation,
                                translation,
                            ) = fit_orthogonal(
                                source_xy,
                                target_xy,
                                chirality,
                            )

                            fitted_xy = (
                                source_xy
                                @ rotation
                                + translation
                            )

                            lateral_distances_nm = (
                                np.linalg.norm(
                                    fitted_xy
                                    - target_xy,
                                    axis=1,
                                )
                            )

                            transformed_all_xy = (
                                all_annulus_xy
                                @ rotation
                                + translation
                            )

                            annulus_center_offset_nm = float(
                                np.linalg.norm(
                                    np.mean(
                                        transformed_all_xy,
                                        axis=0,
                                    )
                                )
                            )

                            gap_result = (
                                optimize_gap_interval(
                                    lateral_distances_nm,
                                    seed_axial_nm,
                                    minimum_span_nm,
                                    maximum_span_nm,
                                )
                            )

                            center_pass = (
                                annulus_center_offset_nm
                                <= MAX_ANNULUS_CENTER_OFFSET_NM
                            )

                            mapping_feasible = (
                                bool(
                                    gap_result[
                                        "gap_feasible"
                                    ]
                                )
                                and center_pass
                            )

                            row = {
                                "end": end,
                                "bridge_atoms_per_attachment": (
                                    bridge_atoms
                                ),
                                "seed_element": (
                                    seed_element
                                ),
                                "bridge_element_sequence": (
                                    "-".join(
                                        bridge_sequence
                                    )
                                ),
                                "required_annulus_element": (
                                    annulus_element
                                ),
                                "seed_parity": (
                                    seed_parity
                                ),
                                "orientation": (
                                    orientation
                                ),
                                "rotation_index": (
                                    rotation_index
                                ),
                                "chirality": (
                                    chirality
                                ),
                                "annulus_center_offset_nm": (
                                    annulus_center_offset_nm
                                ),
                                "lateral_distance_minimum_nm": float(
                                    np.min(
                                        lateral_distances_nm
                                    )
                                ),
                                "lateral_distance_mean_nm": float(
                                    np.mean(
                                        lateral_distances_nm
                                    )
                                ),
                                "lateral_distance_maximum_nm": float(
                                    np.max(
                                        lateral_distances_nm
                                    )
                                ),
                                "bridge_envelope_minimum_nm": (
                                    minimum_span_nm
                                ),
                                "bridge_envelope_maximum_nm": (
                                    maximum_span_nm
                                ),
                                **gap_result,
                                "annulus_center_pass": (
                                    center_pass
                                ),
                                "mapping_geometrically_feasible": (
                                    mapping_feasible
                                ),
                            }

                            mapping_rows.append(
                                row
                            )

    if len(mapping_rows) != EXPECTED_TOTAL_MAPPING_FITS:
        raise RuntimeError(
            "Unexpected mapping-fit count: "
            f"{len(mapping_rows)}/"
            f"{EXPECTED_TOTAL_MAPPING_FITS}"
        )

    write_csv(
        MAPPING_SCREEN_CSV,
        mapping_rows,
    )

    class_rows = []

    feasible_bridge_classes = []

    best_mapping_by_class_end: dict[
        tuple[int, str],
        dict[str, Any]
    ] = {}

    for bridge_atoms in range(
        MIN_BRIDGE_ATOMS,
        MAX_BRIDGE_ATOMS + 1,
    ):
        counts_by_end = {}

        for end in (
            "LOWER",
            "UPPER",
        ):
            rows = [
                row
                for row in mapping_rows
                if int(
                    row[
                        "bridge_atoms_per_attachment"
                    ]
                )
                == bridge_atoms
                and row[
                    "end"
                ]
                == end
            ]

            feasible_rows = [
                row
                for row in rows
                if bool(
                    row[
                        "mapping_geometrically_feasible"
                    ]
                )
            ]

            counts_by_end[
                end
            ] = len(
                feasible_rows
            )

            if feasible_rows:
                feasible_rows.sort(
                    key=lambda row: (
                        -(
                            float(
                                row[
                                    "feasible_gap_maximum_nm"
                                ]
                            )
                            - float(
                                row[
                                    "feasible_gap_minimum_nm"
                                ]
                            )
                        ),
                        float(
                            row[
                                "annulus_center_offset_nm"
                            ]
                        ),
                        (
                            float(
                                row[
                                    "selected_distance_maximum_nm"
                                ]
                            )
                            - float(
                                row[
                                    "selected_distance_minimum_nm"
                                ]
                            )
                        ),
                    )
                )

                best_mapping_by_class_end[
                    (
                        bridge_atoms,
                        end,
                    )
                ] = feasible_rows[0]

        uniform_class_feasible = (
            counts_by_end[
                "LOWER"
            ]
            > 0
            and counts_by_end[
                "UPPER"
            ]
            > 0
        )

        if uniform_class_feasible:
            feasible_bridge_classes.append(
                bridge_atoms
            )

        added_bridge_atoms_per_end = (
            EXPECTED_SELECTED_ATTACHMENTS_PER_END
            * bridge_atoms
        )

        total_heavy_atoms_per_end = (
            BASE_ADDED_HEAVY_ATOMS_PER_END
            + added_bridge_atoms_per_end
        )

        total_H_atoms_per_end = (
            BASE_PASSIVANTS_PER_END
            + added_bridge_atoms_per_end
        )

        heavy_atom_relative_error = (
            abs(
                total_heavy_atoms_per_end
                - TARGET_HEAVY_ATOMS_PER_END
            )
            / TARGET_HEAVY_ATOMS_PER_END
        )

        class_rows.append(
            {
                "bridge_atoms_per_attachment": (
                    bridge_atoms
                ),
                "bonds_per_bridge_path": (
                    bridge_atoms + 1
                ),
                "lower_feasible_mapping_count": (
                    counts_by_end[
                        "LOWER"
                    ]
                ),
                "upper_feasible_mapping_count": (
                    counts_by_end[
                        "UPPER"
                    ]
                ),
                "uniform_class_feasible": (
                    uniform_class_feasible
                ),
                "added_bridge_heavy_atoms_per_end": (
                    added_bridge_atoms_per_end
                ),
                "total_added_heavy_atoms_per_end": (
                    total_heavy_atoms_per_end
                ),
                "total_added_H_atoms_per_end": (
                    total_H_atoms_per_end
                ),
                "heavy_atom_relative_error_from_steric_estimate": (
                    heavy_atom_relative_error
                ),
                "candidate_is_final_chemistry": False,
            }
        )

    write_csv(
        CLASS_SUMMARY_CSV,
        class_rows,
    )

    selected_bridge_atoms = (
        min(
            feasible_bridge_classes
        )
        if feasible_bridge_classes
        else None
    )

    selected_rows = []

    if selected_bridge_atoms is not None:
        selected_class = next(
            row
            for row in class_rows
            if int(
                row[
                    "bridge_atoms_per_attachment"
                ]
            )
            == selected_bridge_atoms
        )

        for end in (
            "LOWER",
            "UPPER",
        ):
            selected_mapping = dict(
                best_mapping_by_class_end[
                    (
                        selected_bridge_atoms,
                        end,
                    )
                ]
            )

            selected_mapping[
                "classification"
            ] = (
                "SELECTED_SHORTEST_FEASIBLE_"
                "BRIDGE_CLASS_MAPPING"
            )

            selected_rows.append(
                selected_mapping
            )

        selected_rows.append(
            {
                "classification": (
                    "SELECTED_CLASS_SUMMARY"
                ),
                **selected_class,
            }
        )

        write_csv(
            SELECTED_CANDIDATE_CSV,
            selected_rows,
        )

    audit_gates = {
        "Gate3F_graph_design_is_accepted": (
            design_summary.get(
                "decision"
            )
            == EXPECTED_DESIGN_DECISION
        ),
        "Gate3G_static_embedding_has_expected_review_decision": (
            static_summary.get(
                "decision"
            )
            == EXPECTED_STATIC_EMBEDDING_DECISION
        ),
        "Gate3G1_direct_junction_is_rejected": (
            direct_summary.get(
                "decision"
            )
            == EXPECTED_DIRECT_LOWER_BOUND_DECISION
        ),
        "three_bridge_classes_were_screened": (
            len(
                class_rows
            )
            == 3
        ),
        "720_mapping_fits_were_screened": (
            len(
                mapping_rows
            )
            == EXPECTED_TOTAL_MAPPING_FITS
        ),
        "all_mapping_metrics_are_finite": all(
            all(
                math.isfinite(
                    float(
                        row[field]
                    )
                )
                for field in (
                    "annulus_center_offset_nm",
                    "lateral_distance_minimum_nm",
                    "lateral_distance_mean_nm",
                    "lateral_distance_maximum_nm",
                    "bridge_envelope_minimum_nm",
                    "bridge_envelope_maximum_nm",
                    "selected_gap_nm",
                    "selected_distance_minimum_nm",
                    "selected_distance_maximum_nm",
                    "maximum_envelope_violation_nm",
                )
            )
            for row in mapping_rows
        ),
        "chain_conformer_envelopes_are_ordered_and_nonempty": all(
            int(
                row[
                    "sample_count"
                ]
            )
            > 0
            and float(
                row[
                    "minimum_end_to_end_nm"
                ]
            )
            < float(
                row[
                    "maximum_end_to_end_nm"
                ]
            )
            <= float(
                row[
                    "maximum_contour_length_nm"
                ]
            )
            + 1.0e-12
            for row in envelope_rows
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

    candidate_identified = (
        audit_integrity_pass
        and selected_bridge_atoms is not None
    )

    decision = (
        PASS_DECISION
        if candidate_identified
        else FAIL_DECISION
    )

    required_next_step = (
        "BUILD_AND_VALIDATE_R2_ALTERNATING_BN_OLIGOMER_BRIDGE_GRAPH"
        if candidate_identified
        else
        "EVALUATE_R2_END_SPECIFIC_ORGANIC_OR_HYBRID_LINKER_CLASSES"
    )

    if selected_bridge_atoms is None:
        selected_bridge_sequence_lower = ""
        selected_bridge_sequence_upper = ""
        selected_total_heavy_atoms_per_end = ""
        selected_total_H_atoms_per_end = ""
        selected_heavy_error = ""
        lower_mapping_count = 0
        upper_mapping_count = 0
    else:
        lower_seed_element = next(
            row[
                "seed_element"
            ]
            for row in mapping_rows
            if row[
                "end"
            ]
            == "LOWER"
        )

        upper_seed_element = next(
            row[
                "seed_element"
            ]
            for row in mapping_rows
            if row[
                "end"
            ]
            == "UPPER"
        )

        selected_bridge_sequence_lower = (
            "-".join(
                bridge_element_sequence(
                    lower_seed_element,
                    selected_bridge_atoms,
                )
            )
        )

        selected_bridge_sequence_upper = (
            "-".join(
                bridge_element_sequence(
                    upper_seed_element,
                    selected_bridge_atoms,
                )
            )
        )

        selected_class = next(
            row
            for row in class_rows
            if int(
                row[
                    "bridge_atoms_per_attachment"
                ]
            )
            == selected_bridge_atoms
        )

        selected_total_heavy_atoms_per_end = (
            selected_class[
                "total_added_heavy_atoms_per_end"
            ]
        )

        selected_total_H_atoms_per_end = (
            selected_class[
                "total_added_H_atoms_per_end"
            ]
        )

        selected_heavy_error = (
            selected_class[
                "heavy_atom_relative_error_from_steric_estimate"
            ]
        )

        lower_mapping_count = int(
            selected_class[
                "lower_feasible_mapping_count"
            ]
        )

        upper_mapping_count = int(
            selected_class[
                "upper_feasible_mapping_count"
            ]
        )

    summary = {
        "decision": decision,
        "bridge_classes_screened": (
            len(
                class_rows
            )
        ),
        "mapping_fits_screened": (
            len(
                mapping_rows
            )
        ),
        "feasible_bridge_classes": (
            " | ".join(
                str(value)
                for value
                in feasible_bridge_classes
            )
        ),
        "selected_shortest_bridge_atoms_per_attachment": (
            ""
            if selected_bridge_atoms is None
            else selected_bridge_atoms
        ),
        "selected_bonds_per_bridge_path": (
            ""
            if selected_bridge_atoms is None
            else selected_bridge_atoms + 1
        ),
        "selected_lower_bridge_element_sequence": (
            selected_bridge_sequence_lower
        ),
        "selected_upper_bridge_element_sequence": (
            selected_bridge_sequence_upper
        ),
        "selected_lower_feasible_mapping_count": (
            lower_mapping_count
        ),
        "selected_upper_feasible_mapping_count": (
            upper_mapping_count
        ),
        "selected_total_added_heavy_atoms_per_end": (
            selected_total_heavy_atoms_per_end
        ),
        "selected_total_added_H_atoms_per_end": (
            selected_total_H_atoms_per_end
        ),
        "selected_heavy_atom_relative_error_from_steric_estimate": (
            selected_heavy_error
        ),
        "candidate_is_final_chemistry": False,
        "bridge_graph_generation_authorized": (
            candidate_identified
        ),
        "coordinate_generation_authorized": False,
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
            in audit_gates.items()
        ],
    )

    AUDIT_JSON.write_text(
        json.dumps(
            {
                "summary": summary,
                "conformer_envelopes": (
                    envelope_rows
                ),
                "class_summaries": (
                    class_rows
                ),
                "selected_candidate_rows": (
                    selected_rows
                ),
                "audit_gates": (
                    audit_gates
                ),
                "limitations": [
                    (
                        "The bridge conformer envelopes use fixed "
                        "B-N bond lengths and provisional internal-angle "
                        "ranges. They are not force-field parameters."
                    ),
                    (
                        "A feasible end-to-end envelope does not prove "
                        "that all 15 bridges can be embedded simultaneously "
                        "without steric clashes."
                    ),
                    (
                        "Every internal bridge atom is provisionally "
                        "assumed to require one hydrogen passivant."
                    ),
                    (
                        "The heavy-atom count is recorded but is not used "
                        "as an absolute rejection gate after the failure "
                        "of the direct junction."
                    ),
                    (
                        "No coordinates, topology, charges, minimization, "
                        "MD, or QM calculation were generated."
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
        {
            "role": (
                "Gate3G_static_coordinates"
            ),
            "file": relative(
                STATIC_COORDINATES_CSV
            ),
            "sha256": sha256(
                STATIC_COORDINATES_CSV
            ),
        },
        {
            "role": (
                "Gate3G_static_embedding_summary"
            ),
            "file": relative(
                STATIC_EMBEDDING_SUMMARY_CSV
            ),
            "sha256": sha256(
                STATIC_EMBEDDING_SUMMARY_CSV
            ),
        },
        {
            "role": (
                "Gate3G1_direct_lower_bound_summary"
            ),
            "file": relative(
                DIRECT_LOWER_BOUND_SUMMARY_CSV
            ),
            "sha256": sha256(
                DIRECT_LOWER_BOUND_SUMMARY_CSV
            ),
        },
    ]

    write_csv(
        SOURCE_MANIFEST_CSV,
        source_rows,
    )

    class_lines = "\n".join(
        (
            f"- {row['bridge_atoms_per_attachment']} bridge atoms: "
            f"lower/upper feasible mappings = "
            f"{row['lower_feasible_mapping_count']}/"
            f"{row['upper_feasible_mapping_count']}; "
            f"heavy atoms/end = "
            f"{row['total_added_heavy_atoms_per_end']}; "
            f"H/end = {row['total_added_H_atoms_per_end']}; "
            f"uniformly feasible = "
            f"{row['uniform_class_feasible']}"
        )
        for row in class_rows
    )

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed
        in audit_gates.items()
    )

    REPORT_MD.write_text(
        f"""# R2 Alternating BN Oligomer-Bridge Feasibility Audit

## Scope

This stage screens alternating BN oligomer bridges containing one,
two or three heavy atoms between the validated edge-completion seed
and the n=5, m=2 BN annulus.

No new graph, coordinates, molecular topology, formal charges,
force-field parameters, minimization, MD, or QM calculation were
generated.

## Provisional conformer model

- B-N bond length:
  **{BN_TARGET_NM:.6f} nm**
- Internal-angle range:
  **{MIN_INTERNAL_ANGLE_DEG}–{MAX_INTERNAL_ANGLE_DEG} degrees**
- Torsion step:
  **{TORSION_STEP_DEG} degrees**
- Axial-gap interval:
  **{MIN_USEFUL_AXIAL_GAP_NM:.3f}–
  {MAX_USEFUL_AXIAL_GAP_NM:.3f} nm**

The sampled end-to-end ranges are geometric screening envelopes, not
energetic or force-field predictions.

## Bridge-class results

{class_lines}

## Selection

- Feasible bridge classes:
  **{'NONE' if not feasible_bridge_classes else ', '.join(str(value) for value in feasible_bridge_classes)}**
- Shortest selected class:
  **{'NONE' if selected_bridge_atoms is None else str(selected_bridge_atoms) + ' bridge atoms per attachment'}**
- Lower bridge sequence:
  **{selected_bridge_sequence_lower or 'NONE'}**
- Upper bridge sequence:
  **{selected_bridge_sequence_upper or 'NONE'}**
- Selected total heavy atoms per end:
  **{selected_total_heavy_atoms_per_end or 'N/A'}**
- Selected total H atoms per end:
  **{selected_total_H_atoms_per_end or 'N/A'}**

## Audit gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed audit-integrity gates:
  **{'NONE' if not failed_audit_gates else ' | '.join(failed_audit_gates)}**
- Bridge-graph generation authorized:
  **{'YES' if candidate_identified else 'NO'}**
- Coordinate generation authorized:
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

## Interpretation

The shortest class passing this gate is only a bridge-length and
endpoint-parity candidate. The next gate must construct the complete
graph, assign every bridge atom and H passivant, verify bipartition,
coordination, connectivity, cycle topology and atom counts, and
confirm that the rejected direct seed-annulus edges have been removed.
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 alternating BN oligomer-bridge "
        "feasibility audit completed."
    )

    for row in envelope_rows:
        print(
            "Bridge atoms / bonds / samples / "
            "end-to-end min/max: "
            f"{row['bridge_atoms_per_attachment']}/"
            f"{row['bond_count_seed_to_annulus']}/"
            f"{row['sample_count']}/"
            f"{float(row['minimum_end_to_end_nm']):.6f}/"
            f"{float(row['maximum_end_to_end_nm']):.6f} nm"
        )

    print(
        "Mapping fits screened: "
        f"{len(mapping_rows)}"
    )

    for row in class_rows:
        print(
            "Bridge class atoms / lower feasible / "
            "upper feasible / heavy atoms per end / "
            "H per end / feasible: "
            f"{row['bridge_atoms_per_attachment']}/"
            f"{row['lower_feasible_mapping_count']}/"
            f"{row['upper_feasible_mapping_count']}/"
            f"{row['total_added_heavy_atoms_per_end']}/"
            f"{row['total_added_H_atoms_per_end']}/"
            f"{row['uniform_class_feasible']}"
        )

    print(
        "Feasible bridge classes: "
        + (
            "NONE"
            if not feasible_bridge_classes
            else "/".join(
                str(value)
                for value
                in feasible_bridge_classes
            )
        )
    )

    print(
        "Selected shortest bridge atoms per attachment: "
        + (
            "NONE"
            if selected_bridge_atoms is None
            else str(
                selected_bridge_atoms
            )
        )
    )

    if selected_bridge_atoms is not None:
        print(
            "Selected lower / upper bridge sequences: "
            f"{selected_bridge_sequence_lower}/"
            f"{selected_bridge_sequence_upper}"
        )

        for row in selected_rows:
            if row.get(
                "classification"
            ) != (
                "SELECTED_SHORTEST_FEASIBLE_"
                "BRIDGE_CLASS_MAPPING"
            ):
                continue

            print(
                f"{row['end']} selected mapping "
                "parity/orientation/rotation/chirality/"
                "gap-range/selected-gap: "
                f"{row['seed_parity']}/"
                f"{row['orientation']}/"
                f"{row['rotation_index']}/"
                f"{row['chirality']}/"
                f"{float(row['feasible_gap_minimum_nm']):.6f}-"
                f"{float(row['feasible_gap_maximum_nm']):.6f}/"
                f"{float(row['selected_gap_nm']):.6f} nm"
            )

            print(
                f"{row['end']} selected endpoint "
                "distance min/max and center offset: "
                f"{float(row['selected_distance_minimum_nm']):.6f}/"
                f"{float(row['selected_distance_maximum_nm']):.6f}/"
                f"{float(row['annulus_center_offset_nm']):.6f} nm"
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
        "Bridge-graph generation authorized: "
        f"{'YES' if candidate_identified else 'NO'}"
    )

    print(
        "Coordinate generation authorized: NO"
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
        CONFORMER_ENVELOPES_CSV,
        MAPPING_SCREEN_CSV,
        CLASS_SUMMARY_CSV,
        SUMMARY_CSV,
        GATES_CSV,
        AUDIT_JSON,
        SOURCE_MANIFEST_CSV,
        REPORT_MD,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if selected_bridge_atoms is not None:
        print(
            f"Wrote: {relative(SELECTED_CANDIDATE_CSV)}"
        )


if __name__ == "__main__":
    main()
