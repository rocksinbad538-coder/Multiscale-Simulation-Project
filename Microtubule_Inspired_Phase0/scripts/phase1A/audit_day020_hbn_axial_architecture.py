#!/usr/bin/env python3

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]

ACCEPTED_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute"
)

GRO_PATH = (
    ACCEPTED_ROOT
    / "nvt_100ps_frozenSolute.gro"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day020_confined_water_axial_radial_density/"
    "hbn_architecture_audit"
)

PLANES_CSV = (
    OUTPUT_ROOT
    / "hbn_axial_planes.csv"
)

SEGMENTS_CSV = (
    OUTPUT_ROOT
    / "hbn_axial_segments.csv"
)

GAPS_CSV = (
    OUTPUT_ROOT
    / "hbn_axial_gaps.csv"
)

PYRENE_CSV = (
    OUTPUT_ROOT
    / "pyrene_positions_vs_hbn_architecture.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "hbn_architecture_summary.csv"
)

FIGURE_STEM = (
    OUTPUT_ROOT
    / "figure_day020_hbn_axial_architecture"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "HBN_AXIAL_ARCHITECTURE_AUDIT_DAY020.md"
)

EXPECTED_HBN_ATOMS = 1680
EXPECTED_PYR_ATOMS = 104
EXPECTED_PYR_RESIDUES = 4

PLANE_CLUSTER_TOLERANCE_NM = 0.020
MINIMUM_SEGMENT_BREAK_NM = 0.35
SEGMENT_BREAK_FACTOR = 4.0
HISTOGRAM_BIN_WIDTH_NM = 0.05

DPI = 400


def log(message: str = "") -> None:
    print(message, flush=True)


def relative(path: Path) -> str:
    return str(path.relative_to(PROJECT_ROOT))


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
) -> None:
    if not rows:
        rows = [{"status": "no_rows"}]

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_gro(
    path: Path,
) -> list[dict[str, object]]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    if len(lines) < 3:
        raise RuntimeError(
            f"Invalid GRO file: {path}"
        )

    atom_count = int(
        lines[1].strip()
    )

    atom_lines = lines[
        2 : 2 + atom_count
    ]

    if len(atom_lines) != atom_count:
        raise RuntimeError(
            "GRO atom-count mismatch"
        )

    atoms: list[
        dict[str, object]
    ] = []

    for atom_index, line in enumerate(
        atom_lines,
        start=1,
    ):
        if len(line) < 44:
            raise RuntimeError(
                f"Malformed GRO line "
                f"for atom {atom_index}"
            )

        atoms.append(
            {
                "atom_index": atom_index,
                "residue_number": int(
                    line[0:5]
                ),
                "residue_name": (
                    line[5:10].strip()
                ),
                "atom_name": (
                    line[10:15].strip()
                ),
                "coordinate_nm": np.asarray(
                    [
                        float(line[20:28]),
                        float(line[28:36]),
                        float(line[36:44]),
                    ],
                    dtype=np.float64,
                ),
            }
        )

    return atoms


def orient_axis(
    axis: np.ndarray,
) -> np.ndarray:
    axis = axis / np.linalg.norm(
        axis
    )

    largest_component = int(
        np.argmax(
            np.abs(axis)
        )
    )

    if axis[
        largest_component
    ] < 0.0:
        axis = -axis

    return axis


def determine_axis(
    coordinates: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    center = np.mean(
        coordinates,
        axis=0,
    )

    centered = (
        coordinates
        - center
    )

    covariance = (
        centered.T
        @ centered
        / coordinates.shape[0]
    )

    eigenvalues, eigenvectors = (
        np.linalg.eigh(
            covariance
        )
    )

    order = np.argsort(
        eigenvalues
    )[::-1]

    eigenvalues = (
        eigenvalues[
            order
        ]
    )

    axis = orient_axis(
        eigenvectors[
            :,
            order[0]
        ]
    )

    return (
        center,
        axis,
        eigenvalues,
    )


def project_coordinates(
    coordinates: np.ndarray,
    center: np.ndarray,
    axis: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:
    relative = (
        coordinates
        - center
    )

    axial = relative @ axis

    radial_vectors = (
        relative
        - np.outer(
            axial,
            axis,
        )
    )

    radial = np.linalg.norm(
        radial_vectors,
        axis=1,
    )

    return axial, radial


def cluster_axial_planes(
    axial_positions_nm: np.ndarray,
) -> tuple[
    list[np.ndarray],
    np.ndarray,
]:
    order = np.argsort(
        axial_positions_nm
    )

    clusters: list[
        list[int]
    ] = []

    current_cluster = [
        int(order[0])
    ]

    previous_value = float(
        axial_positions_nm[
            order[0]
        ]
    )

    for raw_index in order[1:]:
        index = int(raw_index)

        value = float(
            axial_positions_nm[
                index
            ]
        )

        if (
            value
            - previous_value
        ) <= (
            PLANE_CLUSTER_TOLERANCE_NM
        ):
            current_cluster.append(
                index
            )
        else:
            clusters.append(
                current_cluster
            )

            current_cluster = [
                index
            ]

        previous_value = value

    clusters.append(
        current_cluster
    )

    cluster_arrays = [
        np.asarray(
            cluster,
            dtype=np.int64,
        )
        for cluster in clusters
    ]

    plane_centers = np.asarray(
        [
            np.mean(
                axial_positions_nm[
                    cluster
                ]
            )
            for cluster
            in cluster_arrays
        ],
        dtype=np.float64,
    )

    return (
        cluster_arrays,
        plane_centers,
    )


def determine_spacing_and_breaks(
    plane_centers_nm: np.ndarray,
) -> tuple[
    np.ndarray,
    float,
    float,
    np.ndarray,
]:
    if plane_centers_nm.size < 3:
        raise RuntimeError(
            "Too few HBN axial planes"
        )

    spacings = np.diff(
        plane_centers_nm
    )

    if not np.all(
        spacings > 0.0
    ):
        raise RuntimeError(
            "HBN planes are not ordered"
        )

    upper_quartile = float(
        np.quantile(
            spacings,
            0.75,
        )
    )

    local_spacings = spacings[
        spacings
        <= upper_quartile
    ]

    if local_spacings.size == 0:
        raise RuntimeError(
            "Could not determine local "
            "HBN plane spacing"
        )

    typical_spacing_nm = float(
        np.median(
            local_spacings
        )
    )

    break_threshold_nm = max(
        MINIMUM_SEGMENT_BREAK_NM,
        SEGMENT_BREAK_FACTOR
        * typical_spacing_nm,
    )

    break_indices = np.where(
        spacings
        > break_threshold_nm
    )[0]

    return (
        spacings,
        typical_spacing_nm,
        break_threshold_nm,
        break_indices,
    )


def build_segments(
    plane_clusters: list[np.ndarray],
    plane_centers_nm: np.ndarray,
    break_indices: np.ndarray,
    typical_spacing_nm: float,
    axial_positions_nm: np.ndarray,
    radial_positions_nm: np.ndarray,
) -> tuple[
    list[dict[str, object]],
    np.ndarray,
]:
    segment_plane_ranges = []

    start_plane = 0

    for break_index in break_indices:
        end_plane = int(
            break_index
        )

        segment_plane_ranges.append(
            (
                start_plane,
                end_plane,
            )
        )

        start_plane = (
            end_plane + 1
        )

    segment_plane_ranges.append(
        (
            start_plane,
            plane_centers_nm.size - 1,
        )
    )

    atom_segment_ids = np.zeros(
        axial_positions_nm.size,
        dtype=np.int64,
    )

    rows: list[
        dict[str, object]
    ] = []

    for segment_id, (
        first_plane,
        last_plane,
    ) in enumerate(
        segment_plane_ranges,
        start=1,
    ):
        member_indices = np.concatenate(
            plane_clusters[
                first_plane : last_plane + 1
            ]
        )

        atom_segment_ids[
            member_indices
        ] = segment_id

        segment_axial = (
            axial_positions_nm[
                member_indices
            ]
        )

        segment_radial = (
            radial_positions_nm[
                member_indices
            ]
        )

        lower_boundary_nm = (
            plane_centers_nm[
                first_plane
            ]
            - 0.5
            * typical_spacing_nm
        )

        upper_boundary_nm = (
            plane_centers_nm[
                last_plane
            ]
            + 0.5
            * typical_spacing_nm
        )

        rows.append(
            {
                "segment_id": segment_id,
                "first_plane_id": (
                    first_plane + 1
                ),
                "last_plane_id": (
                    last_plane + 1
                ),
                "plane_count": (
                    last_plane
                    - first_plane
                    + 1
                ),
                "atom_count": (
                    member_indices.size
                ),
                "minimum_atom_z_nm": float(
                    np.min(
                        segment_axial
                    )
                ),
                "maximum_atom_z_nm": float(
                    np.max(
                        segment_axial
                    )
                ),
                "lower_boundary_nm": (
                    lower_boundary_nm
                ),
                "upper_boundary_nm": (
                    upper_boundary_nm
                ),
                "boundary_span_nm": (
                    upper_boundary_nm
                    - lower_boundary_nm
                ),
                "mean_radius_nm": float(
                    np.mean(
                        segment_radial
                    )
                ),
                "median_radius_nm": float(
                    np.median(
                        segment_radial
                    )
                ),
                "p05_radius_nm": float(
                    np.percentile(
                        segment_radial,
                        5.0,
                    )
                ),
                "p95_radius_nm": float(
                    np.percentile(
                        segment_radial,
                        95.0,
                    )
                ),
            }
        )

    if np.any(
        atom_segment_ids == 0
    ):
        raise RuntimeError(
            "Some HBN atoms were not assigned "
            "to a segment"
        )

    return rows, atom_segment_ids


def build_gap_rows(
    segment_rows: list[
        dict[str, object]
    ],
    plane_centers_nm: np.ndarray,
    break_indices: np.ndarray,
) -> list[dict[str, object]]:
    rows: list[
        dict[str, object]
    ] = []

    for gap_id, break_index in enumerate(
        break_indices,
        start=1,
    ):
        left_segment = (
            segment_rows[
                gap_id - 1
            ]
        )

        right_segment = (
            segment_rows[
                gap_id
            ]
        )

        left_plane_nm = float(
            plane_centers_nm[
                int(
                    break_index
                )
            ]
        )

        right_plane_nm = float(
            plane_centers_nm[
                int(
                    break_index
                )
                + 1
            ]
        )

        lower_boundary_nm = float(
            left_segment[
                "upper_boundary_nm"
            ]
        )

        upper_boundary_nm = float(
            right_segment[
                "lower_boundary_nm"
            ]
        )

        rows.append(
            {
                "gap_id": gap_id,
                "left_segment_id": gap_id,
                "right_segment_id": (
                    gap_id + 1
                ),
                "left_terminal_plane_nm": (
                    left_plane_nm
                ),
                "right_terminal_plane_nm": (
                    right_plane_nm
                ),
                "plane_to_plane_gap_nm": (
                    right_plane_nm
                    - left_plane_nm
                ),
                "lower_boundary_nm": (
                    lower_boundary_nm
                ),
                "upper_boundary_nm": (
                    upper_boundary_nm
                ),
                "clear_gap_width_nm": (
                    upper_boundary_nm
                    - lower_boundary_nm
                ),
                "gap_center_nm": (
                    0.5
                    * (
                        lower_boundary_nm
                        + upper_boundary_nm
                    )
                ),
            }
        )

    return rows


def classify_axial_position(
    z_nm: float,
    segment_rows: list[
        dict[str, object]
    ],
    gap_rows: list[
        dict[str, object]
    ],
) -> tuple[str, int | str]:
    for segment in segment_rows:
        if (
            float(
                segment[
                    "lower_boundary_nm"
                ]
            )
            <= z_nm
            <= float(
                segment[
                    "upper_boundary_nm"
                ]
            )
        ):
            return (
                "inside_hbn_segment",
                int(
                    segment[
                        "segment_id"
                    ]
                ),
            )

    for gap in gap_rows:
        if (
            float(
                gap[
                    "lower_boundary_nm"
                ]
            )
            < z_nm
            < float(
                gap[
                    "upper_boundary_nm"
                ]
            )
        ):
            return (
                "inside_hbn_gap",
                int(
                    gap[
                        "gap_id"
                    ]
                ),
            )

    return (
        "outside_hbn_axial_extent",
        "",
    )


def build_pyrene_rows(
    pyrene_atoms: list[
        dict[str, object]
    ],
    center: np.ndarray,
    axis: np.ndarray,
    segment_rows: list[
        dict[str, object]
    ],
    gap_rows: list[
        dict[str, object]
    ],
) -> list[dict[str, object]]:
    by_residue: dict[
        int,
        list[dict[str, object]],
    ] = defaultdict(list)

    for atom in pyrene_atoms:
        by_residue[
            int(
                atom[
                    "residue_number"
                ]
            )
        ].append(atom)

    if len(by_residue) != (
        EXPECTED_PYR_RESIDUES
    ):
        raise RuntimeError(
            f"Expected "
            f"{EXPECTED_PYR_RESIDUES} "
            f"PYR residues, found "
            f"{len(by_residue)}"
        )

    rows = []

    for pyrene_id, (
        residue_number,
        atoms,
    ) in enumerate(
        sorted(
            by_residue.items()
        ),
        start=1,
    ):
        coordinates = np.stack(
            [
                atom[
                    "coordinate_nm"
                ]
                for atom in atoms
            ]
        )

        pyrene_center = np.mean(
            coordinates,
            axis=0,
        )

        relative = (
            pyrene_center
            - center
        )

        axial_position = float(
            relative @ axis
        )

        radial_vector = (
            relative
            - axial_position
            * axis
        )

        radial_position = float(
            np.linalg.norm(
                radial_vector
            )
        )

        axial_class, associated_id = (
            classify_axial_position(
                axial_position,
                segment_rows,
                gap_rows,
            )
        )

        rows.append(
            {
                "pyrene_label": (
                    f"PYR_{pyrene_id}"
                ),
                "gro_residue_number": (
                    residue_number
                ),
                "atom_count": len(atoms),
                "center_x_nm": (
                    pyrene_center[0]
                ),
                "center_y_nm": (
                    pyrene_center[1]
                ),
                "center_z_nm": (
                    pyrene_center[2]
                ),
                "axial_position_nm": (
                    axial_position
                ),
                "radial_position_nm": (
                    radial_position
                ),
                "axial_architecture_class": (
                    axial_class
                ),
                "associated_segment_or_gap_id": (
                    associated_id
                ),
            }
        )

    return rows


def build_figure(
    hbn_axial_nm: np.ndarray,
    hbn_radial_nm: np.ndarray,
    pyrene_rows: list[
        dict[str, object]
    ],
    segment_rows: list[
        dict[str, object]
    ],
    gap_rows: list[
        dict[str, object]
    ],
    typical_spacing_nm: float,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.5,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "legend.fontsize": 9,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "axes.linewidth": 0.9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    figure, axes = plt.subplots(
        2,
        1,
        figsize=(10.5, 8.0),
        sharex=True,
    )

    axis_a, axis_b = axes

    histogram_min = (
        float(
            np.min(
                hbn_axial_nm
            )
        )
        - HISTOGRAM_BIN_WIDTH_NM
    )

    histogram_max = (
        float(
            np.max(
                hbn_axial_nm
            )
        )
        + HISTOGRAM_BIN_WIDTH_NM
    )

    histogram_edges = np.arange(
        histogram_min,
        histogram_max
        + HISTOGRAM_BIN_WIDTH_NM,
        HISTOGRAM_BIN_WIDTH_NM,
    )

    axis_a.hist(
        hbn_axial_nm,
        bins=histogram_edges,
    )

    axis_a.set_ylabel(
        "HBN atoms per axial bin"
    )

    axis_a.set_title(
        "Axial occupancy of the HBN scaffold"
    )

    axis_a.grid(
        True,
        alpha=0.22,
    )

    axis_a.text(
        0.02,
        0.95,
        "(a)",
        transform=axis_a.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
    )

    axis_b.scatter(
        hbn_axial_nm,
        hbn_radial_nm,
        s=7,
        alpha=0.45,
        label="HBN atoms",
    )

    for row in pyrene_rows:
        axial_position = float(
            row[
                "axial_position_nm"
            ]
        )

        radial_position = float(
            row[
                "radial_position_nm"
            ]
        )

        axis_b.scatter(
            axial_position,
            radial_position,
            marker="D",
            s=48,
            edgecolor="black",
            linewidth=0.8,
            label=(
                "PYR centers"
                if row[
                    "pyrene_label"
                ] == "PYR_1"
                else None
            ),
            zorder=5,
        )

        axis_b.annotate(
            str(
                row[
                    "pyrene_label"
                ]
            ),
            (
                axial_position,
                radial_position,
            ),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8.5,
        )

    for segment in segment_rows:
        lower = float(
            segment[
                "lower_boundary_nm"
            ]
        )

        upper = float(
            segment[
                "upper_boundary_nm"
            ]
        )

        for axis_plot in axes:
            axis_plot.axvspan(
                lower,
                upper,
                alpha=0.10,
            )

        axis_a.text(
            0.5 * (
                lower + upper
            ),
            0.88
            * axis_a.get_ylim()[1],
            (
                f"HBN segment "
                f"{segment['segment_id']}"
            ),
            ha="center",
            va="top",
            fontsize=8.5,
        )

    for gap in gap_rows:
        lower = float(
            gap[
                "lower_boundary_nm"
            ]
        )

        upper = float(
            gap[
                "upper_boundary_nm"
            ]
        )

        for axis_plot in axes:
            axis_plot.axvspan(
                lower,
                upper,
                alpha=0.06,
                hatch="//",
                edgecolor="0.45",
            )

    axis_b.set_xlabel(
        "Coordinate along HBN principal axis (nm)"
    )

    axis_b.set_ylabel(
        "Radial coordinate (nm)"
    )

    axis_b.set_title(
        "HBN and pyrene geometry in axial–radial space"
    )

    axis_b.grid(
        True,
        alpha=0.22,
    )

    axis_b.legend(
        frameon=False,
    )

    axis_b.text(
        0.02,
        0.95,
        "(b)",
        transform=axis_b.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
    )

    figure.suptitle(
        (
            "Geometry-derived axial architecture "
            "of the frozen HBN–pyrene assembly"
        ),
        fontsize=14,
        y=0.985,
    )

    figure.text(
        0.99,
        0.01,
        (
            "Typical HBN plane spacing: "
            f"{typical_spacing_nm:.4f} nm"
        ),
        ha="right",
        va="bottom",
        fontsize=8.5,
    )

    figure.tight_layout(
        rect=(
            0.0,
            0.03,
            1.0,
            0.96,
        )
    )

    figure.savefig(
        FIGURE_STEM.with_suffix(
            ".png"
        ),
        dpi=DPI,
        bbox_inches="tight",
    )

    figure.savefig(
        FIGURE_STEM.with_suffix(
            ".pdf"
        ),
        bbox_inches="tight",
    )

    plt.close(figure)


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not GRO_PATH.exists():
        raise RuntimeError(
            f"Missing GRO input: "
            f"{GRO_PATH}"
        )

    atoms = parse_gro(
        GRO_PATH
    )

    hbn_atoms = [
        atom
        for atom in atoms
        if atom[
            "residue_name"
        ] == "HBN"
    ]

    pyrene_atoms = [
        atom
        for atom in atoms
        if atom[
            "residue_name"
        ] == "PYR"
    ]

    if len(hbn_atoms) != (
        EXPECTED_HBN_ATOMS
    ):
        raise RuntimeError(
            f"Expected {EXPECTED_HBN_ATOMS} "
            f"HBN atoms, found "
            f"{len(hbn_atoms)}"
        )

    if len(pyrene_atoms) != (
        EXPECTED_PYR_ATOMS
    ):
        raise RuntimeError(
            f"Expected {EXPECTED_PYR_ATOMS} "
            f"PYR atoms, found "
            f"{len(pyrene_atoms)}"
        )

    hbn_coordinates = np.stack(
        [
            atom[
                "coordinate_nm"
            ]
            for atom in hbn_atoms
        ]
    )

    (
        hbn_center,
        hbn_axis,
        pca_eigenvalues,
    ) = determine_axis(
        hbn_coordinates
    )

    (
        hbn_axial_nm,
        hbn_radial_nm,
    ) = project_coordinates(
        hbn_coordinates,
        hbn_center,
        hbn_axis,
    )

    (
        plane_clusters,
        plane_centers_nm,
    ) = cluster_axial_planes(
        hbn_axial_nm
    )

    (
        plane_spacings_nm,
        typical_spacing_nm,
        break_threshold_nm,
        break_indices,
    ) = determine_spacing_and_breaks(
        plane_centers_nm
    )

    (
        segment_rows,
        atom_segment_ids,
    ) = build_segments(
        plane_clusters,
        plane_centers_nm,
        break_indices,
        typical_spacing_nm,
        hbn_axial_nm,
        hbn_radial_nm,
    )

    gap_rows = build_gap_rows(
        segment_rows,
        plane_centers_nm,
        break_indices,
    )

    pyrene_rows = build_pyrene_rows(
        pyrene_atoms,
        hbn_center,
        hbn_axis,
        segment_rows,
        gap_rows,
    )

    plane_rows: list[
        dict[str, object]
    ] = []

    for plane_id, (
        cluster,
        center_nm,
    ) in enumerate(
        zip(
            plane_clusters,
            plane_centers_nm,
        ),
        start=1,
    ):
        segment_ids = np.unique(
            atom_segment_ids[
                cluster
            ]
        )

        if segment_ids.size != 1:
            raise RuntimeError(
                "An axial plane was assigned "
                "to multiple segments"
            )

        spacing_to_next = (
            float(
                plane_centers_nm[
                    plane_id
                ]
                - center_nm
            )
            if plane_id
            < plane_centers_nm.size
            else ""
        )

        plane_rows.append(
            {
                "plane_id": plane_id,
                "segment_id": int(
                    segment_ids[0]
                ),
                "center_z_nm": float(
                    center_nm
                ),
                "minimum_z_nm": float(
                    np.min(
                        hbn_axial_nm[
                            cluster
                        ]
                    )
                ),
                "maximum_z_nm": float(
                    np.max(
                        hbn_axial_nm[
                            cluster
                        ]
                    )
                ),
                "atom_count": (
                    cluster.size
                ),
                "spacing_to_next_plane_nm": (
                    spacing_to_next
                ),
            }
        )

    write_csv(
        PLANES_CSV,
        plane_rows,
    )

    write_csv(
        SEGMENTS_CSV,
        segment_rows,
    )

    write_csv(
        GAPS_CSV,
        gap_rows,
    )

    write_csv(
        PYRENE_CSV,
        pyrene_rows,
    )

    largest_gap_nm = (
        max(
            float(
                row[
                    "clear_gap_width_nm"
                ]
            )
            for row in gap_rows
        )
        if gap_rows
        else 0.0
    )

    summary_rows = [
        {
            "source_GRO": relative(
                GRO_PATH
            ),
            "HBN_atom_count": len(
                hbn_atoms
            ),
            "PYR_atom_count": len(
                pyrene_atoms
            ),
            "HBN_center_x_nm": (
                hbn_center[0]
            ),
            "HBN_center_y_nm": (
                hbn_center[1]
            ),
            "HBN_center_z_nm": (
                hbn_center[2]
            ),
            "axis_x": hbn_axis[0],
            "axis_y": hbn_axis[1],
            "axis_z": hbn_axis[2],
            "PCA_eigenvalue_1_nm2": (
                pca_eigenvalues[0]
            ),
            "PCA_eigenvalue_2_nm2": (
                pca_eigenvalues[1]
            ),
            "PCA_eigenvalue_3_nm2": (
                pca_eigenvalues[2]
            ),
            "axial_plane_count": (
                plane_centers_nm.size
            ),
            "typical_plane_spacing_nm": (
                typical_spacing_nm
            ),
            "segment_break_threshold_nm": (
                break_threshold_nm
            ),
            "detected_segment_count": (
                len(segment_rows)
            ),
            "detected_gap_count": (
                len(gap_rows)
            ),
            "largest_clear_gap_nm": (
                largest_gap_nm
            ),
            "minimum_HBN_radius_nm": float(
                np.min(
                    hbn_radial_nm
                )
            ),
            "mean_HBN_radius_nm": float(
                np.mean(
                    hbn_radial_nm
                )
            ),
            "maximum_HBN_radius_nm": float(
                np.max(
                    hbn_radial_nm
                )
            ),
            "all_HBN_atoms_segmented_pass": bool(
                np.all(
                    atom_segment_ids > 0
                )
            ),
            "overall_validation_pass": True,
        }
    ]

    write_csv(
        SUMMARY_CSV,
        summary_rows,
    )

    build_figure(
        hbn_axial_nm,
        hbn_radial_nm,
        pyrene_rows,
        segment_rows,
        gap_rows,
        typical_spacing_nm,
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 HBN Axial "
            "Architecture Audit\n\n"
        )

        handle.write(
            "## Purpose\n\n"
        )

        handle.write(
            "This audit reconstructs the actual "
            "axial architecture of the frozen HBN "
            "scaffold directly from atomic "
            "coordinates. It does not assume that "
            "the scaffold is one continuous "
            "cylindrical wall.\n\n"
        )

        handle.write(
            "## Principal axis\n\n"
        )

        handle.write(
            f"- HBN axis: "
            f"({hbn_axis[0]:.8f}, "
            f"{hbn_axis[1]:.8f}, "
            f"{hbn_axis[2]:.8f}).\n"
        )

        handle.write(
            f"- HBN axial planes: "
            f"{plane_centers_nm.size}.\n"
        )

        handle.write(
            f"- Typical neighboring-plane spacing: "
            f"{typical_spacing_nm:.6f} nm.\n"
        )

        handle.write(
            f"- Segment-break threshold: "
            f"{break_threshold_nm:.6f} nm.\n\n"
        )

        handle.write(
            "## Detected scaffold segments\n\n"
        )

        for row in segment_rows:
            handle.write(
                f"- Segment "
                f"{row['segment_id']}: "
                f"{row['lower_boundary_nm']:.6f} to "
                f"{row['upper_boundary_nm']:.6f} nm; "
                f"{row['atom_count']} atoms; "
                f"mean radius "
                f"{row['mean_radius_nm']:.6f} nm.\n"
            )

        handle.write(
            "\n## Detected axial gaps\n\n"
        )

        if gap_rows:
            for row in gap_rows:
                handle.write(
                    f"- Gap {row['gap_id']}: "
                    f"{row['lower_boundary_nm']:.6f} to "
                    f"{row['upper_boundary_nm']:.6f} nm; "
                    f"clear width "
                    f"{row['clear_gap_width_nm']:.6f} nm.\n"
                )
        else:
            handle.write(
                "- No large axial gap was detected.\n"
            )

        handle.write(
            "\n## Pyrene positions\n\n"
        )

        for row in pyrene_rows:
            handle.write(
                f"- {row['pyrene_label']}: "
                f"z={row['axial_position_nm']:.6f} nm, "
                f"r={row['radial_position_nm']:.6f} nm, "
                f"class="
                f"{row['axial_architecture_class']}.\n"
            )

        handle.write(
            "\n## Consequence for water-region "
            "classification\n\n"
        )

        handle.write(
            "Water regions must be defined using "
            "the detected HBN segments and gaps. "
            "The previous continuous-cylinder "
            "partition is retained only as a "
            "diagnostic calculation and must not "
            "be used as the final physical "
            "classification unless this audit "
            "detects a single continuous segment.\n"
        )

    expected_outputs = (
        PLANES_CSV,
        SEGMENTS_CSV,
        GAPS_CSV,
        PYRENE_CSV,
        SUMMARY_CSV,
        FIGURE_STEM.with_suffix(
            ".png"
        ),
        FIGURE_STEM.with_suffix(
            ".pdf"
        ),
        REPORT_MD,
    )

    missing = [
        path
        for path in expected_outputs
        if (
            not path.exists()
            or path.stat().st_size == 0
        )
    ]

    if missing:
        raise RuntimeError(
            "Missing or empty outputs:\n"
            + "\n".join(
                str(path)
                for path in missing
            )
        )

    log(
        "Day020 HBN axial-architecture "
        "audit completed."
    )

    log(
        f"HBN axial planes: "
        f"{plane_centers_nm.size}"
    )

    log(
        "Typical plane spacing: "
        f"{typical_spacing_nm:.6f} nm"
    )

    log(
        "Segment-break threshold: "
        f"{break_threshold_nm:.6f} nm"
    )

    log(
        f"Detected HBN segments: "
        f"{len(segment_rows)}"
    )

    for row in segment_rows:
        log(
            f"Segment {row['segment_id']}: "
            f"{row['lower_boundary_nm']:.6f} to "
            f"{row['upper_boundary_nm']:.6f} nm; "
            f"{row['atom_count']} atoms"
        )

    log(
        f"Detected HBN gaps: "
        f"{len(gap_rows)}"
    )

    for row in gap_rows:
        log(
            f"Gap {row['gap_id']}: "
            f"{row['lower_boundary_nm']:.6f} to "
            f"{row['upper_boundary_nm']:.6f} nm; "
            f"clear width "
            f"{row['clear_gap_width_nm']:.6f} nm"
        )

    for row in pyrene_rows:
        log(
            f"{row['pyrene_label']}: "
            f"z={row['axial_position_nm']:.6f} nm, "
            f"r={row['radial_position_nm']:.6f} nm, "
            f"{row['axial_architecture_class']}"
        )

    log(
        "All HBN atoms segmented: PASS"
    )

    if len(segment_rows) == 1:
        log(
            "Architecture interpretation: "
            "continuous HBN segment"
        )
    else:
        log(
            "Architecture interpretation: "
            "discontinuous segmented HBN scaffold"
        )

    log(
        "Overall validation: PASS"
    )

    log(
        f"Wrote: "
        f"{relative(OUTPUT_ROOT)}"
    )


if __name__ == "__main__":
    main()
