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

OUTPUT_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "09_r2_direct_junction_geometric_lower_bound"
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

EMBEDDING_COORDINATES_CSV = (
    GATE3G_ROOT
    / "r2_partial_attachment_static_coordinates.csv"
)

EMBEDDING_SUMMARY_CSV = (
    GATE3G_ROOT
    / "r2_partial_attachment_static_embedding_summary.csv"
)

MAPPING_FITS_CSV = (
    OUTPUT_ROOT
    / "r2_direct_junction_all_mapping_fits.csv"
)

BEST_FITS_CSV = (
    OUTPUT_ROOT
    / "r2_direct_junction_best_fits.csv"
)

END_SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_direct_junction_geometric_lower_bound_end_summary.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_direct_junction_geometric_lower_bound_summary.csv"
)

GATES_CSV = (
    OUTPUT_ROOT
    / "r2_direct_junction_geometric_lower_bound_gates.csv"
)

AUDIT_JSON = (
    OUTPUT_ROOT
    / "r2_direct_junction_geometric_lower_bound.json"
)

SOURCE_MANIFEST_CSV = (
    OUTPUT_ROOT
    / "r2_direct_junction_geometric_lower_bound_source_manifest.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_DIRECT_JUNCTION_GEOMETRIC_LOWER_BOUND_DAY024.md"
)

EXPECTED_DESIGN_DECISION = (
    "R2_PARTIAL_HETEROPOLAR_ANNULUS_ATTACHMENT_AND_"
    "COMPLEMENTARY_PASSIVATION_GRAPH_VALIDATED"
)

EXPECTED_EMBEDDING_DECISION = (
    "R2_PARTIAL_ATTACHMENT_ANNULUS_STATIC_COORDINATE_EMBEDDING_"
    "REQUIRES_CONSTRAINED_GEOMETRIC_OPTIMIZATION"
)

POSITIVE_DECISION = (
    "R2_PARTIAL_ATTACHMENT_LOCAL_CONSTRAINED_OPTIMIZATION_JUSTIFIED"
)

NEGATIVE_DECISION = (
    "R2_PARTIAL_ATTACHMENT_DIRECT_BN_JUNCTION_GEOMETRIC_LOWER_BOUND_FAILED"
)

BN_TARGET_NM = 0.144973

MAX_GAP_NM = 0.250

MAX_ATTACHMENT_RMS_DEVIATION_NM = 0.020
MAX_ATTACHMENT_MAX_DEVIATION_NM = 0.035

MAX_ANNULUS_INTERNAL_BOND_DEVIATION_NM = 0.010
MAX_CENTER_OFFSET_NM = 0.050

MAX_PRINCIPAL_STRAIN = 0.050
MAX_AFFINE_ANISOTROPY = 1.10

EXPECTED_SEED_SITES_PER_END = 30
EXPECTED_SELECTED_SEED_SITES_PER_END = 15
EXPECTED_ATTACHABLE_ANNULUS_SITES_PER_END = 15

EXPECTED_FITS_PER_END = 300
EXPECTED_TOTAL_FITS = 600


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
        rows = list(
            csv.DictReader(handle)
        )

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
            f"Expected one row in {path}; "
            f"found {len(rows)}"
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

    eigenvalues, eigenvectors = (
        np.linalg.eigh(
            covariance
        )
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

    left, _, right_t = (
        np.linalg.svd(
            covariance
        )
    )

    rotation = (
        left
        @ right_t
    )

    observed_det = int(
        round(
            np.linalg.det(
                rotation
            )
        )
    )

    if observed_det != desired_determinant:
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


def fit_similarity(
    source: np.ndarray,
    target: np.ndarray,
    desired_determinant: int,
) -> tuple[
    np.ndarray,
    np.ndarray,
    float,
]:
    rotation, _ = fit_orthogonal(
        source,
        target,
        desired_determinant,
    )

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

    rotated = (
        source_centered
        @ rotation
    )

    denominator = float(
        np.sum(
            source_centered
            ** 2
        )
    )

    if denominator <= 1.0e-16:
        raise RuntimeError(
            "Degenerate similarity source geometry."
        )

    scale = float(
        np.sum(
            rotated
            * target_centered
        )
        / denominator
    )

    if scale <= 0.0:
        raise RuntimeError(
            "Non-positive similarity scale."
        )

    linear = (
        scale
        * rotation
    )

    translation = (
        target_center
        - source_center
        @ linear
    )

    return linear, translation, scale


def fit_affine(
    source: np.ndarray,
    target: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:
    augmented = np.column_stack(
        (
            source,
            np.ones(
                len(source),
                dtype=float,
            ),
        )
    )

    coefficients, _, _, _ = (
        np.linalg.lstsq(
            augmented,
            target,
            rcond=None,
        )
    )

    linear = coefficients[
        :2,
        :,
    ]

    translation = coefficients[
        2,
        :,
    ]

    return linear, translation


def optimize_gap(
    lateral_distances_nm: np.ndarray,
    seed_axial_nm: np.ndarray,
) -> tuple[
    float,
    np.ndarray,
    float,
    float,
]:
    mean_seed_axial_nm = float(
        np.mean(
            seed_axial_nm
        )
    )

    gaps = np.linspace(
        0.0,
        MAX_GAP_NM,
        25001,
    )

    plane_positions = (
        mean_seed_axial_nm
        + gaps
    )

    axial_differences = (
        plane_positions[
            None,
            :,
        ]
        - seed_axial_nm[
            :,
            None,
        ]
    )

    lengths = np.sqrt(
        lateral_distances_nm[
            :,
            None,
        ]
        ** 2
        + axial_differences
        ** 2
    )

    deviations = (
        lengths
        - BN_TARGET_NM
    )

    rms = np.sqrt(
        np.mean(
            deviations
            ** 2,
            axis=0,
        )
    )

    best_index = int(
        np.argmin(
            rms
        )
    )

    selected_lengths = (
        lengths[
            :,
            best_index
        ]
    )

    selected_deviations = (
        selected_lengths
        - BN_TARGET_NM
    )

    return (
        float(
            gaps[
                best_index
            ]
        ),
        selected_lengths,
        float(
            np.sqrt(
                np.mean(
                    selected_deviations
                    ** 2
                )
            )
        ),
        float(
            np.max(
                np.abs(
                    selected_deviations
                )
            )
        ),
    )


def annulus_bond_statistics(
    end: str,
    linear: np.ndarray,
    local_coordinates: dict[
        str,
        np.ndarray,
    ],
    edge_rows: list[
        dict[str, str]
    ],
) -> tuple[
    float,
    float,
    float,
    float,
    int,
]:
    lengths = []

    for row in edge_rows:
        if (
            row.get(
                "end"
            )
            != end
            or row.get(
                "edge_type"
            )
            != "ANNULUS_BN"
        ):
            continue

        first = row[
            "source_node"
        ]

        second = row[
            "target_node"
        ]

        first_xy = (
            local_coordinates[
                first
            ]
            @ linear
        )

        second_xy = (
            local_coordinates[
                second
            ]
            @ linear
        )

        lengths.append(
            float(
                np.linalg.norm(
                    first_xy
                    - second_xy
                )
            )
        )

    if not lengths:
        raise RuntimeError(
            f"No annulus bonds found for {end}."
        )

    values = np.asarray(
        lengths,
        dtype=float,
    )

    deviations = (
        values
        - BN_TARGET_NM
    )

    return (
        float(
            np.mean(values)
        ),
        float(
            np.min(values)
        ),
        float(
            np.max(values)
        ),
        float(
            np.max(
                np.abs(
                    deviations
                )
            )
        ),
        len(values),
    )


def evaluate_transform(
    *,
    end: str,
    mapping_id: str,
    fit_mode: str,
    seed_parity: int,
    orientation: int,
    rotation_index: int,
    chirality: int,
    source_xy: np.ndarray,
    target_xy: np.ndarray,
    target_axial_nm: np.ndarray,
    linear: np.ndarray,
    translation: np.ndarray,
    all_annulus_xy: np.ndarray,
    local_coordinates: dict[str, np.ndarray],
    edge_rows: list[dict[str, str]],
) -> dict[str, Any]:
    fitted = (
        source_xy
        @ linear
        + translation
    )

    lateral_residuals = np.linalg.norm(
        fitted
        - target_xy,
        axis=1,
    )

    (
        optimal_gap_nm,
        bond_lengths_nm,
        attachment_rms_deviation_nm,
        attachment_max_deviation_nm,
    ) = optimize_gap(
        lateral_residuals,
        target_axial_nm,
    )

    (
        annulus_bond_mean_nm,
        annulus_bond_minimum_nm,
        annulus_bond_maximum_nm,
        annulus_bond_max_deviation_nm,
        annulus_bond_count,
    ) = annulus_bond_statistics(
        end,
        linear,
        local_coordinates,
        edge_rows,
    )

    singular_values = np.linalg.svd(
        linear,
        compute_uv=False,
    )

    singular_max = float(
        np.max(
            singular_values
        )
    )

    singular_min = float(
        np.min(
            singular_values
        )
    )

    maximum_principal_strain = max(
        abs(
            singular_max
            - 1.0
        ),
        abs(
            singular_min
            - 1.0
        ),
    )

    anisotropy = (
        singular_max
        / singular_min
        if singular_min > 1.0e-12
        else 1.0e12
    )

    transformed_all = (
        all_annulus_xy
        @ linear
        + translation
    )

    annulus_center_offset_nm = float(
        np.linalg.norm(
            np.mean(
                transformed_all,
                axis=0,
            )
        )
    )

    attachment_pass = (
        attachment_rms_deviation_nm
        <= MAX_ATTACHMENT_RMS_DEVIATION_NM
        and attachment_max_deviation_nm
        <= MAX_ATTACHMENT_MAX_DEVIATION_NM
    )

    annulus_strain_pass = (
        annulus_bond_max_deviation_nm
        <= MAX_ANNULUS_INTERNAL_BOND_DEVIATION_NM
        and maximum_principal_strain
        <= MAX_PRINCIPAL_STRAIN
        and anisotropy
        <= MAX_AFFINE_ANISOTROPY
    )

    center_pass = (
        annulus_center_offset_nm
        <= MAX_CENTER_OFFSET_NM
    )

    positive_gap_pass = (
        0.0
        < optimal_gap_nm
        <= MAX_GAP_NM
    )

    local_preoptimization_pass = (
        attachment_pass
        and annulus_strain_pass
        and center_pass
        and positive_gap_pass
    )

    return {
        "end": end,
        "mapping_id": mapping_id,
        "fit_mode": fit_mode,
        "seed_parity": seed_parity,
        "orientation": orientation,
        "rotation_index": rotation_index,
        "chirality": chirality,
        "linear_determinant": float(
            np.linalg.det(
                linear
            )
        ),
        "singular_value_minimum": (
            singular_min
        ),
        "singular_value_maximum": (
            singular_max
        ),
        "maximum_principal_strain": (
            maximum_principal_strain
        ),
        "affine_anisotropy": (
            anisotropy
        ),
        "annulus_center_offset_nm": (
            annulus_center_offset_nm
        ),
        "optimal_axial_gap_nm": (
            optimal_gap_nm
        ),
        "lateral_residual_mean_nm": float(
            np.mean(
                lateral_residuals
            )
        ),
        "lateral_residual_minimum_nm": float(
            np.min(
                lateral_residuals
            )
        ),
        "lateral_residual_maximum_nm": float(
            np.max(
                lateral_residuals
            )
        ),
        "attachment_bond_mean_nm": float(
            np.mean(
                bond_lengths_nm
            )
        ),
        "attachment_bond_minimum_nm": float(
            np.min(
                bond_lengths_nm
            )
        ),
        "attachment_bond_maximum_nm": float(
            np.max(
                bond_lengths_nm
            )
        ),
        "attachment_RMS_deviation_nm": (
            attachment_rms_deviation_nm
        ),
        "attachment_maximum_deviation_nm": (
            attachment_max_deviation_nm
        ),
        "annulus_bond_count": (
            annulus_bond_count
        ),
        "annulus_bond_mean_nm": (
            annulus_bond_mean_nm
        ),
        "annulus_bond_minimum_nm": (
            annulus_bond_minimum_nm
        ),
        "annulus_bond_maximum_nm": (
            annulus_bond_maximum_nm
        ),
        "annulus_bond_maximum_deviation_nm": (
            annulus_bond_max_deviation_nm
        ),
        "attachment_thresholds_pass": (
            attachment_pass
        ),
        "annulus_strain_thresholds_pass": (
            annulus_strain_pass
        ),
        "annulus_center_threshold_pass": (
            center_pass
        ),
        "positive_axial_gap_pass": (
            positive_gap_pass
        ),
        "local_preoptimization_pass": (
            local_preoptimization_pass
        ),
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
        EMBEDDING_COORDINATES_CSV,
        EMBEDDING_SUMMARY_CSV,
    ):
        require_file(required)

    design_nodes = read_csv_rows(
        DESIGN_NODES_CSV
    )

    design_edges = read_csv_rows(
        DESIGN_EDGES_CSV
    )

    design_summary = read_single_csv_row(
        DESIGN_SUMMARY_CSV
    )

    embedding_rows = read_csv_rows(
        EMBEDDING_COORDINATES_CSV
    )

    embedding_summary = read_single_csv_row(
        EMBEDDING_SUMMARY_CSV
    )

    if design_summary.get(
        "decision"
    ) != EXPECTED_DESIGN_DECISION:
        raise RuntimeError(
            "Gate 3F is not in the accepted state."
        )

    if embedding_summary.get(
        "decision"
    ) != EXPECTED_EMBEDDING_DECISION:
        raise RuntimeError(
            "Gate 3G does not contain the expected "
            "constrained-optimization decision."
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
        for row in embedding_rows
    }

    if set(nodes_by_id) != set(
        coordinates
    ):
        missing_coordinates = (
            set(nodes_by_id)
            - set(coordinates)
        )

        extra_coordinates = (
            set(coordinates)
            - set(nodes_by_id)
        )

        raise RuntimeError(
            "Node/coordinate mismatch. Missing: "
            f"{len(missing_coordinates)}; extra: "
            f"{len(extra_coordinates)}"
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
            for node_id
            in parent_ids
        ],
        dtype=float,
    )

    _, tube_axis = determine_axis(
        parent_positions
    )

    local_coordinates = {}

    for node_id, row in nodes_by_id.items():
        if row[
            "node_type"
        ] not in {
            "ANNULUS_INTERIOR",
            "ANNULUS_OUTER_BOUNDARY",
            "ANNULUS_INNER_BOUNDARY",
        }:
            continue

        local_coordinates[
            node_id
        ] = local_xy(
            parse_int(
                row,
                "lattice_x",
            ),
            parse_int(
                row,
                "lattice_y",
            ),
        )

    fit_rows = []
    end_summary_rows = []

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
                f"{end}: seed count "
                f"{len(seed_rows)}/"
                f"{EXPECTED_SEED_SITES_PER_END}"
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

        complementary_element = (
            "B"
            if seed_element == "N"
            else "N"
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
            == complementary_element
        ]

        outer_rows.sort(
            key=lambda row: parse_float(
                row,
                "angle_turns",
            )
        )

        if (
            len(outer_rows)
            != EXPECTED_ATTACHABLE_ANNULUS_SITES_PER_END
        ):
            raise RuntimeError(
                f"{end}: attachable annulus sites "
                f"{len(outer_rows)}/"
                f"{EXPECTED_ATTACHABLE_ANNULUS_SITES_PER_END}"
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
                local_coordinates[
                    row[
                        "node_id"
                    ]
                ]
                for row in all_annulus_rows
            ],
            dtype=float,
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

        source_base = np.asarray(
            [
                local_coordinates[
                    row[
                        "node_id"
                    ]
                ]
                for row in outer_rows
            ],
            dtype=float,
        )

        end_fit_rows = []

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
                != EXPECTED_SELECTED_SEED_SITES_PER_END
            ):
                raise RuntimeError(
                    f"{end}: parity {seed_parity} did not "
                    "select 15 seed sites."
                )

            target_positions = np.asarray(
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
                    for position in target_positions
                ],
                dtype=float,
            )

            target_axial_nm = np.asarray(
                [
                    float(
                        np.dot(
                            position
                            - seed_center,
                            outward,
                        )
                    )
                    for position in target_positions
                ],
                dtype=float,
            )

            for orientation in (
                1,
                -1,
            ):
                for rotation_index in range(
                    EXPECTED_ATTACHABLE_ANNULUS_SITES_PER_END
                ):
                    mapped_indices = [
                        (
                            orientation
                            * index
                            + rotation_index
                        )
                        % EXPECTED_ATTACHABLE_ANNULUS_SITES_PER_END
                        for index in range(
                            EXPECTED_ATTACHABLE_ANNULUS_SITES_PER_END
                        )
                    ]

                    source_xy = (
                        source_base[
                            mapped_indices
                        ]
                    )

                    mapping_id = (
                        f"{end}:P{seed_parity}:"
                        f"O{orientation}:R{rotation_index}"
                    )

                    for chirality in (
                        1,
                        -1,
                    ):
                        (
                            rigid_linear,
                            rigid_translation,
                        ) = fit_orthogonal(
                            source_xy,
                            target_xy,
                            chirality,
                        )

                        row = evaluate_transform(
                            end=end,
                            mapping_id=mapping_id,
                            fit_mode="RIGID",
                            seed_parity=seed_parity,
                            orientation=orientation,
                            rotation_index=rotation_index,
                            chirality=chirality,
                            source_xy=source_xy,
                            target_xy=target_xy,
                            target_axial_nm=target_axial_nm,
                            linear=rigid_linear,
                            translation=rigid_translation,
                            all_annulus_xy=all_annulus_xy,
                            local_coordinates=local_coordinates,
                            edge_rows=design_edges,
                        )

                        fit_rows.append(row)
                        end_fit_rows.append(row)

                        (
                            similarity_linear,
                            similarity_translation,
                            similarity_scale,
                        ) = fit_similarity(
                            source_xy,
                            target_xy,
                            chirality,
                        )

                        row = evaluate_transform(
                            end=end,
                            mapping_id=mapping_id,
                            fit_mode="SIMILARITY",
                            seed_parity=seed_parity,
                            orientation=orientation,
                            rotation_index=rotation_index,
                            chirality=chirality,
                            source_xy=source_xy,
                            target_xy=target_xy,
                            target_axial_nm=target_axial_nm,
                            linear=similarity_linear,
                            translation=similarity_translation,
                            all_annulus_xy=all_annulus_xy,
                            local_coordinates=local_coordinates,
                            edge_rows=design_edges,
                        )

                        row[
                            "similarity_scale"
                        ] = similarity_scale

                        fit_rows.append(row)
                        end_fit_rows.append(row)

                    (
                        affine_linear,
                        affine_translation,
                    ) = fit_affine(
                        source_xy,
                        target_xy,
                    )

                    row = evaluate_transform(
                        end=end,
                        mapping_id=mapping_id,
                        fit_mode="AFFINE",
                        seed_parity=seed_parity,
                        orientation=orientation,
                        rotation_index=rotation_index,
                        chirality=0,
                        source_xy=source_xy,
                        target_xy=target_xy,
                        target_axial_nm=target_axial_nm,
                        linear=affine_linear,
                        translation=affine_translation,
                        all_annulus_xy=all_annulus_xy,
                        local_coordinates=local_coordinates,
                        edge_rows=design_edges,
                    )

                    fit_rows.append(row)
                    end_fit_rows.append(row)

        if len(end_fit_rows) != EXPECTED_FITS_PER_END:
            raise RuntimeError(
                f"{end}: fit count "
                f"{len(end_fit_rows)}/"
                f"{EXPECTED_FITS_PER_END}"
            )

        best_by_mode = {}

        for fit_mode in (
            "RIGID",
            "SIMILARITY",
            "AFFINE",
        ):
            candidates = [
                row
                for row in end_fit_rows
                if row[
                    "fit_mode"
                ]
                == fit_mode
            ]

            candidates.sort(
                key=lambda row: (
                    float(
                        row[
                            "attachment_RMS_deviation_nm"
                        ]
                    ),
                    float(
                        row[
                            "attachment_maximum_deviation_nm"
                        ]
                    ),
                    float(
                        row[
                            "maximum_principal_strain"
                        ]
                    ),
                    float(
                        row[
                            "annulus_center_offset_nm"
                        ]
                    ),
                )
            )

            best_by_mode[
                fit_mode
            ] = candidates[0]

        feasible_rows = [
            row
            for row in end_fit_rows
            if bool(
                row[
                    "local_preoptimization_pass"
                ]
            )
        ]

        best_any = min(
            end_fit_rows,
            key=lambda row: (
                float(
                    row[
                        "attachment_RMS_deviation_nm"
                    ]
                ),
                float(
                    row[
                        "attachment_maximum_deviation_nm"
                    ]
                ),
                float(
                    row[
                        "maximum_principal_strain"
                    ]
                ),
            ),
        )

        end_summary_rows.append(
            {
                "end": end,
                "fits_screened": len(
                    end_fit_rows
                ),
                "local_preoptimization_pass_count": (
                    len(
                        feasible_rows
                    )
                ),
                "best_overall_fit_mode": (
                    best_any[
                        "fit_mode"
                    ]
                ),
                "best_overall_mapping_id": (
                    best_any[
                        "mapping_id"
                    ]
                ),
                "best_overall_attachment_RMS_deviation_nm": (
                    best_any[
                        "attachment_RMS_deviation_nm"
                    ]
                ),
                "best_overall_attachment_max_deviation_nm": (
                    best_any[
                        "attachment_maximum_deviation_nm"
                    ]
                ),
                "best_overall_gap_nm": (
                    best_any[
                        "optimal_axial_gap_nm"
                    ]
                ),
                "best_overall_maximum_principal_strain": (
                    best_any[
                        "maximum_principal_strain"
                    ]
                ),
                "best_overall_affine_anisotropy": (
                    best_any[
                        "affine_anisotropy"
                    ]
                ),
                "best_rigid_RMS_deviation_nm": (
                    best_by_mode[
                        "RIGID"
                    ][
                        "attachment_RMS_deviation_nm"
                    ]
                ),
                "best_rigid_max_deviation_nm": (
                    best_by_mode[
                        "RIGID"
                    ][
                        "attachment_maximum_deviation_nm"
                    ]
                ),
                "best_rigid_gap_nm": (
                    best_by_mode[
                        "RIGID"
                    ][
                        "optimal_axial_gap_nm"
                    ]
                ),
                "best_similarity_RMS_deviation_nm": (
                    best_by_mode[
                        "SIMILARITY"
                    ][
                        "attachment_RMS_deviation_nm"
                    ]
                ),
                "best_similarity_max_deviation_nm": (
                    best_by_mode[
                        "SIMILARITY"
                    ][
                        "attachment_maximum_deviation_nm"
                    ]
                ),
                "best_similarity_principal_strain": (
                    best_by_mode[
                        "SIMILARITY"
                    ][
                        "maximum_principal_strain"
                    ]
                ),
                "best_similarity_annulus_bond_max_deviation_nm": (
                    best_by_mode[
                        "SIMILARITY"
                    ][
                        "annulus_bond_maximum_deviation_nm"
                    ]
                ),
                "best_affine_RMS_deviation_nm": (
                    best_by_mode[
                        "AFFINE"
                    ][
                        "attachment_RMS_deviation_nm"
                    ]
                ),
                "best_affine_max_deviation_nm": (
                    best_by_mode[
                        "AFFINE"
                    ][
                        "attachment_maximum_deviation_nm"
                    ]
                ),
                "best_affine_principal_strain": (
                    best_by_mode[
                        "AFFINE"
                    ][
                        "maximum_principal_strain"
                    ]
                ),
                "best_affine_anisotropy": (
                    best_by_mode[
                        "AFFINE"
                    ][
                        "affine_anisotropy"
                    ]
                ),
                "local_constrained_optimization_feasible": (
                    len(
                        feasible_rows
                    )
                    > 0
                ),
            }
        )

    if len(fit_rows) != EXPECTED_TOTAL_FITS:
        raise RuntimeError(
            "Unexpected total fit count: "
            f"{len(fit_rows)}/"
            f"{EXPECTED_TOTAL_FITS}"
        )

    write_csv(
        MAPPING_FITS_CSV,
        fit_rows,
    )

    best_fit_rows = []

    for end in (
        "LOWER",
        "UPPER",
    ):
        end_rows = [
            row
            for row in fit_rows
            if row[
                "end"
            ]
            == end
        ]

        for fit_mode in (
            "RIGID",
            "SIMILARITY",
            "AFFINE",
        ):
            mode_rows = [
                row
                for row in end_rows
                if row[
                    "fit_mode"
                ]
                == fit_mode
            ]

            mode_rows.sort(
                key=lambda row: (
                    float(
                        row[
                            "attachment_RMS_deviation_nm"
                        ]
                    ),
                    float(
                        row[
                            "attachment_maximum_deviation_nm"
                        ]
                    ),
                    float(
                        row[
                            "maximum_principal_strain"
                        ]
                    ),
                )
            )

            best_fit_rows.append(
                {
                    "classification": (
                        f"BEST_{fit_mode}"
                    ),
                    **mode_rows[0],
                }
            )

        passing_rows = [
            row
            for row in end_rows
            if bool(
                row[
                    "local_preoptimization_pass"
                ]
            )
        ]

        if passing_rows:
            passing_rows.sort(
                key=lambda row: (
                    float(
                        row[
                            "attachment_RMS_deviation_nm"
                        ]
                    ),
                    float(
                        row[
                            "maximum_principal_strain"
                        ]
                    ),
                )
            )

            best_fit_rows.append(
                {
                    "classification": (
                        "BEST_LOCAL_PREOPTIMIZATION_PASS"
                    ),
                    **passing_rows[0],
                }
            )

    write_csv(
        BEST_FITS_CSV,
        best_fit_rows,
    )

    write_csv(
        END_SUMMARY_CSV,
        end_summary_rows,
    )

    lower = next(
        row
        for row in end_summary_rows
        if row[
            "end"
        ]
        == "LOWER"
    )

    upper = next(
        row
        for row in end_summary_rows
        if row[
            "end"
        ]
        == "UPPER"
    )

    local_optimization_feasible = (
        bool(
            lower[
                "local_constrained_optimization_feasible"
            ]
        )
        and bool(
            upper[
                "local_constrained_optimization_feasible"
            ]
        )
    )

    audit_gates = {
        "Gate3F_graph_design_is_accepted": (
            design_summary.get(
                "decision"
            )
            == EXPECTED_DESIGN_DECISION
        ),
        "Gate3G_embedding_has_expected_review_decision": (
            embedding_summary.get(
                "decision"
            )
            == EXPECTED_EMBEDDING_DECISION
        ),
        "600_mapping_and_transform_fits_were_screened": (
            len(
                fit_rows
            )
            == EXPECTED_TOTAL_FITS
        ),
        "300_fits_were_screened_per_end": all(
            int(
                row[
                    "fits_screened"
                ]
            )
            == EXPECTED_FITS_PER_END
            for row in end_summary_rows
        ),
        "all_required_fit_modes_are_present": all(
            any(
                row[
                    "end"
                ]
                == end
                and row[
                    "fit_mode"
                ]
                == fit_mode
                for row in fit_rows
            )
            for end in (
                "LOWER",
                "UPPER",
            )
            for fit_mode in (
                "RIGID",
                "SIMILARITY",
                "AFFINE",
            )
        ),
        "all_core_fit_metrics_are_finite": all(
            all(
                math.isfinite(
                    float(
                        row[field]
                    )
                )
                for field in (
                    "attachment_RMS_deviation_nm",
                    "attachment_maximum_deviation_nm",
                    "annulus_bond_maximum_deviation_nm",
                    "maximum_principal_strain",
                    "affine_anisotropy",
                    "annulus_center_offset_nm",
                    "optimal_axial_gap_nm",
                )
            )
            for row in fit_rows
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

    optimization_authorized = (
        audit_integrity_pass
        and local_optimization_feasible
    )

    decision = (
        POSITIVE_DECISION
        if optimization_authorized
        else NEGATIVE_DECISION
    )

    required_next_step = (
        "RUN_R2_PARTIAL_ATTACHMENT_LOCAL_CONSTRAINED_GEOMETRIC_OPTIMIZATION"
        if optimization_authorized
        else
        "REDESIGN_R2_ANNULUS_JUNCTION_WITH_EXPLICIT_BRIDGING_LINKER_"
        "OR_REVISED_ATTACHMENT_TOPOLOGY"
    )

    summary = {
        "decision": decision,
        "fits_screened": len(
            fit_rows
        ),
        "lower_local_preoptimization_pass_count": (
            lower[
                "local_preoptimization_pass_count"
            ]
        ),
        "upper_local_preoptimization_pass_count": (
            upper[
                "local_preoptimization_pass_count"
            ]
        ),
        "lower_best_rigid_RMS_deviation_nm": (
            lower[
                "best_rigid_RMS_deviation_nm"
            ]
        ),
        "upper_best_rigid_RMS_deviation_nm": (
            upper[
                "best_rigid_RMS_deviation_nm"
            ]
        ),
        "lower_best_similarity_RMS_deviation_nm": (
            lower[
                "best_similarity_RMS_deviation_nm"
            ]
        ),
        "upper_best_similarity_RMS_deviation_nm": (
            upper[
                "best_similarity_RMS_deviation_nm"
            ]
        ),
        "lower_best_similarity_principal_strain": (
            lower[
                "best_similarity_principal_strain"
            ]
        ),
        "upper_best_similarity_principal_strain": (
            upper[
                "best_similarity_principal_strain"
            ]
        ),
        "lower_best_affine_RMS_deviation_nm": (
            lower[
                "best_affine_RMS_deviation_nm"
            ]
        ),
        "upper_best_affine_RMS_deviation_nm": (
            upper[
                "best_affine_RMS_deviation_nm"
            ]
        ),
        "lower_best_affine_principal_strain": (
            lower[
                "best_affine_principal_strain"
            ]
        ),
        "upper_best_affine_principal_strain": (
            upper[
                "best_affine_principal_strain"
            ]
        ),
        "lower_best_affine_anisotropy": (
            lower[
                "best_affine_anisotropy"
            ]
        ),
        "upper_best_affine_anisotropy": (
            upper[
                "best_affine_anisotropy"
            ]
        ),
        "audit_integrity_pass": (
            audit_integrity_pass
        ),
        "local_constrained_optimization_authorized": (
            optimization_authorized
        ),
        "current_direct_BN_attachment_graph_retained": (
            optimization_authorized
        ),
        "current_direct_BN_attachment_graph_rejected": (
            not optimization_authorized
        ),
        "coordinate_update_authorized": False,
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
                "end_summaries": (
                    end_summary_rows
                ),
                "audit_gates": (
                    audit_gates
                ),
                "thresholds": {
                    "attachment_RMS_deviation_nm": (
                        MAX_ATTACHMENT_RMS_DEVIATION_NM
                    ),
                    "attachment_max_deviation_nm": (
                        MAX_ATTACHMENT_MAX_DEVIATION_NM
                    ),
                    "annulus_internal_bond_deviation_nm": (
                        MAX_ANNULUS_INTERNAL_BOND_DEVIATION_NM
                    ),
                    "maximum_principal_strain": (
                        MAX_PRINCIPAL_STRAIN
                    ),
                    "maximum_affine_anisotropy": (
                        MAX_AFFINE_ANISOTROPY
                    ),
                    "maximum_center_offset_nm": (
                        MAX_CENTER_OFFSET_NM
                    ),
                },
                "interpretation": {
                    "rigid_fit": (
                        "Tests whether any cyclic heteropolar mapping "
                        "works without in-plane annulus deformation."
                    ),
                    "similarity_fit": (
                        "Allows uniform annulus expansion or contraction."
                    ),
                    "affine_fit": (
                        "Provides a permissive lower bound that allows "
                        "global anisotropic deformation."
                    ),
                    "negative_result": (
                        "Failure even under mild affine deformation blocks "
                        "local constrained optimization of the current "
                        "direct B-N attachment graph."
                    ),
                },
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
                EMBEDDING_COORDINATES_CSV
            ),
            "sha256": sha256(
                EMBEDDING_COORDINATES_CSV
            ),
        },
        {
            "role": (
                "Gate3G_embedding_summary"
            ),
            "file": relative(
                EMBEDDING_SUMMARY_CSV
            ),
            "sha256": sha256(
                EMBEDDING_SUMMARY_CSV
            ),
        },
    ]

    write_csv(
        SOURCE_MANIFEST_CSV,
        source_rows,
    )

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed in audit_gates.items()
    )

    end_lines = "\n".join(
        (
            f"### {row['end']}\n\n"
            f"- Fits screened: **{row['fits_screened']}**\n"
            f"- Passing local-preoptimization fits: "
            f"**{row['local_preoptimization_pass_count']}**\n"
            f"- Best rigid RMS/max deviation: "
            f"**{float(row['best_rigid_RMS_deviation_nm']):.6f}/"
            f"{float(row['best_rigid_max_deviation_nm']):.6f} nm**\n"
            f"- Best similarity RMS/max deviation: "
            f"**{float(row['best_similarity_RMS_deviation_nm']):.6f}/"
            f"{float(row['best_similarity_max_deviation_nm']):.6f} nm**\n"
            f"- Best similarity principal strain: "
            f"**{float(row['best_similarity_principal_strain']):.6f}**\n"
            f"- Best affine RMS/max deviation: "
            f"**{float(row['best_affine_RMS_deviation_nm']):.6f}/"
            f"{float(row['best_affine_max_deviation_nm']):.6f} nm**\n"
            f"- Best affine principal strain/anisotropy: "
            f"**{float(row['best_affine_principal_strain']):.6f}/"
            f"{float(row['best_affine_anisotropy']):.6f}**\n"
            f"- Local constrained optimization feasible: "
            f"**{row['local_constrained_optimization_feasible']}**"
        )
        for row in end_summary_rows
    )

    REPORT_MD.write_text(
        f"""# R2 Direct-Junction Geometric Lower-Bound Audit

## Scope

This stage determines whether any cyclic mapping between the
15 alternating seed sites and the 15 complementary annulus sites can
support a direct B-N junction with moderate deformation.

No coordinates were replaced. No molecular topology, force-field
parameters, minimization, MD, or QM calculation was generated.

## Search

- Fits per end:
  **{EXPECTED_FITS_PER_END}**
- Total fits:
  **{len(fit_rows)}**
- Mapping variables:
  seed parity, circumferential orientation, discrete rotation and
  transformation chirality.
- Transformation models:
  rigid, similarity and affine.

{end_lines}

## Acceptance thresholds

- Attachment RMS deviation:
  **≤ {MAX_ATTACHMENT_RMS_DEVIATION_NM:.3f} nm**
- Attachment maximum deviation:
  **≤ {MAX_ATTACHMENT_MAX_DEVIATION_NM:.3f} nm**
- Annulus internal-bond maximum deviation:
  **≤ {MAX_ANNULUS_INTERNAL_BOND_DEVIATION_NM:.3f} nm**
- Maximum principal strain:
  **≤ {MAX_PRINCIPAL_STRAIN:.3f}**
- Maximum affine anisotropy:
  **≤ {MAX_AFFINE_ANISOTROPY:.3f}**
- Annulus-center offset:
  **≤ {MAX_CENTER_OFFSET_NM:.3f} nm**
- Axial gap:
  **positive and ≤ {MAX_GAP_NM:.3f} nm**

## Audit gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed audit-integrity gates:
  **{'NONE' if not failed_audit_gates else ' | '.join(failed_audit_gates)}**
- Local constrained optimization authorized:
  **{'YES' if optimization_authorized else 'NO'}**
- Current direct B-N attachment graph retained:
  **{'YES' if optimization_authorized else 'NO'}**
- Coordinate update authorized:
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

## Interpretation

A negative result means that the current graph cannot be repaired by a
local geometric relaxation without requiring excessive annulus strain,
anisotropy, junction displacement, or nonlocal bond deformation.

In that case, the direct seed-to-annulus B-N edges must be replaced by
an explicit bridging junction or by a different attachment topology.
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 direct-junction geometric "
        "lower-bound audit completed."
    )

    print(
        "Mapping/transform fits screened total: "
        f"{len(fit_rows)}"
    )

    for row in end_summary_rows:
        print(
            f"{row['end']} fits / local-pass count: "
            f"{row['fits_screened']}/"
            f"{row['local_preoptimization_pass_count']}"
        )

        print(
            f"{row['end']} best rigid RMS/max/gap: "
            f"{float(row['best_rigid_RMS_deviation_nm']):.6f}/"
            f"{float(row['best_rigid_max_deviation_nm']):.6f}/"
            f"{float(row['best_rigid_gap_nm']):.6f} nm"
        )

        print(
            f"{row['end']} best similarity RMS/max/"
            "strain/annulus-bond-maxdev: "
            f"{float(row['best_similarity_RMS_deviation_nm']):.6f}/"
            f"{float(row['best_similarity_max_deviation_nm']):.6f}/"
            f"{float(row['best_similarity_principal_strain']):.6f}/"
            f"{float(row['best_similarity_annulus_bond_max_deviation_nm']):.6f}"
        )

        print(
            f"{row['end']} best affine RMS/max/"
            "strain/anisotropy: "
            f"{float(row['best_affine_RMS_deviation_nm']):.6f}/"
            f"{float(row['best_affine_max_deviation_nm']):.6f}/"
            f"{float(row['best_affine_principal_strain']):.6f}/"
            f"{float(row['best_affine_anisotropy']):.6f}"
        )

        print(
            f"{row['end']} local constrained "
            "optimization feasible: "
            f"{row['local_constrained_optimization_feasible']}"
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
        "Local constrained optimization authorized: "
        f"{'YES' if optimization_authorized else 'NO'}"
    )

    print(
        "Current direct B-N attachment graph retained: "
        f"{'YES' if optimization_authorized else 'NO'}"
    )

    print(
        "Coordinate update authorized: NO"
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
        MAPPING_FITS_CSV,
        BEST_FITS_CSV,
        END_SUMMARY_CSV,
        SUMMARY_CSV,
        GATES_CSV,
        AUDIT_JSON,
        SOURCE_MANIFEST_CSV,
        REPORT_MD,
    ):
        print(
            f"Wrote: {relative(path)}"
        )


if __name__ == "__main__":
    main()
