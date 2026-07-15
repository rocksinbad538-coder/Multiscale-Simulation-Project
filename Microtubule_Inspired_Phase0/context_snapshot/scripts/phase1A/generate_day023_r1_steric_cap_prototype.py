#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

DAY023_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design"
)

REFERENCE_ROOT = (
    DAY023_ROOT
    / "01_r0_t0_reference"
)

OUTPUT_ROOT = (
    DAY023_ROOT
    / "02_r1_steric_cap_prototype"
)

REFERENCE_GRO = (
    REFERENCE_ROOT
    / "r0_accepted_t0_hydrated_system.gro"
)

REFERENCE_SUMMARY = (
    REFERENCE_ROOT
    / "r0_t0_geometry_and_hydration_summary.csv"
)

CANDIDATE_CSV = (
    OUTPUT_ROOT
    / "r1_steric_cap_candidate_scan.csv"
)

SELECTED_CAP_GRO = (
    OUTPUT_ROOT
    / "r1_selected_steric_caps_only.gro"
)

SELECTED_SYSTEM_GRO = (
    OUTPUT_ROOT
    / "r1_t0_hydrated_with_steric_caps_geometry_only.gro"
)

REMOVED_WATER_CSV = (
    OUTPUT_ROOT
    / "r1_removed_water_molecules_due_to_cap_overlap.csv"
)

SELECTED_CAP_CSV = (
    OUTPUT_ROOT
    / "r1_selected_cap_atom_coordinates.csv"
)

SELECTED_JSON = (
    OUTPUT_ROOT
    / "r1_selected_steric_cap_definition.json"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r1_steric_cap_prototype_summary.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R1_STERIC_CAP_PROTOTYPE_DAY023.md"
)

EXPECTED_ATOMS = 68320
HBN_ATOMS = 1680
PYR_ATOMS = 104
SOLUTE_ATOMS = HBN_ATOMS + PYR_ATOMS
WATER_SITES = 4
EXPECTED_WATERS = 16634
AUTHORITATIVE_LUMEN_WATERS = 428

OFFSETS_NM = (
    0.20,
    0.25,
    0.30,
)

LATTICE_SPACINGS_NM = (
    0.18,
    0.20,
    0.22,
)

CAP_EDGE_MARGIN_NM = 0.16

# Provisional geometric parameters only. They will be replaced by
# explicitly validated nonbonded parameters before any MD execution.
CAP_EFFECTIVE_WATER_EXCLUSION_NM = 0.17
R0_END_EFFECTIVE_WATER_EXCLUSION_NM = 0.17
CAP_WATER_PRUNING_CUTOFF_NM = 0.22

MIN_CAP_SOLUTE_DISTANCE_NM = 0.20
MAX_COVERAGE_HOLE_NM = 0.14
MIN_RETAINED_LUMEN_FRACTION = 0.90
MAX_REMOVED_TOTAL_WATER_FRACTION = 0.02
MIN_AXIAL_STERIC_OVERLAP_NM = 0.02
MIN_DISK_EDGE_MARGIN_NM = 0.10

COVERAGE_GRID_SPACING_NM = 0.025


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
            f"Missing or empty file: {path}"
        )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            block = handle.read(
                1024 * 1024
            )

            if not block:
                break

            digest.update(block)

    return digest.hexdigest()


def read_single_csv_row(
    path: Path,
) -> dict[str, str]:
    require_file(path)

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(
            csv.DictReader(handle)
        )

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected one row in {path}; "
            f"found {len(rows)}"
        )

    return rows[0]


def parse_box(
    values: list[float],
) -> np.ndarray:
    if len(values) == 3:
        return np.diag(
            np.asarray(
                values,
                dtype=float,
            )
        )

    if len(values) == 9:
        return np.asarray(
            [
                [
                    values[0],
                    values[3],
                    values[4],
                ],
                [
                    values[5],
                    values[1],
                    values[6],
                ],
                [
                    values[7],
                    values[8],
                    values[2],
                ],
            ],
            dtype=float,
        )

    raise RuntimeError(
        "Unsupported GRO box representation: "
        f"{len(values)} fields"
    )


def read_gro(
    path: Path,
) -> tuple[
    str,
    list[dict[str, Any]],
    np.ndarray,
]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    if len(lines) < 3:
        raise RuntimeError(
            f"Malformed GRO file: {path}"
        )

    title = lines[0]
    natoms = int(
        lines[1].strip()
    )

    if len(lines) < natoms + 3:
        raise RuntimeError(
            f"Incomplete GRO file: {path}"
        )

    atoms: list[dict[str, Any]] = []

    for zero_index, line in enumerate(
        lines[2 : 2 + natoms]
    ):
        if len(line) < 44:
            raise RuntimeError(
                f"Malformed GRO atom line "
                f"{zero_index + 1}"
            )

        atoms.append(
            {
                "original_index": zero_index,
                "resid": int(
                    line[0:5]
                ),
                "resname": line[
                    5:10
                ].strip(),
                "atomname": line[
                    10:15
                ].strip(),
                "atomnum": int(
                    line[15:20]
                ),
                "position": np.asarray(
                    [
                        float(
                            line[20:28]
                        ),
                        float(
                            line[28:36]
                        ),
                        float(
                            line[36:44]
                        ),
                    ],
                    dtype=float,
                ),
            }
        )

    box = parse_box(
        [
            float(value)
            for value in lines[
                2 + natoms
            ].split()
        ]
    )

    return (
        title,
        atoms,
        box,
    )


def orthorhombic_lengths(
    box: np.ndarray,
) -> np.ndarray:
    diagonal = np.diag(
        np.diag(
            box
        )
    )

    if not np.allclose(
        box,
        diagonal,
        atol=1.0e-8,
        rtol=0.0,
    ):
        raise RuntimeError(
            "R1 prototype generation currently "
            "requires an orthorhombic box."
        )

    lengths = np.diag(
        box
    ).copy()

    if np.any(
        lengths <= 0.0
    ):
        raise RuntimeError(
            f"Invalid box lengths: {lengths}"
        )

    return lengths


def minimum_image(
    displacement: np.ndarray,
    box_lengths: np.ndarray,
) -> np.ndarray:
    return (
        displacement
        - box_lengths
        * np.round(
            displacement
            / box_lengths
        )
    )


def circular_center(
    positions: np.ndarray,
    box_lengths: np.ndarray,
) -> np.ndarray:
    fractional = (
        positions
        / box_lengths
    )

    angles = (
        2.0
        * np.pi
        * fractional
    )

    mean_angle = np.arctan2(
        np.mean(
            np.sin(
                angles
            ),
            axis=0,
        ),
        np.mean(
            np.cos(
                angles
            ),
            axis=0,
        ),
    )

    return (
        (
            mean_angle
            / (
                2.0
                * np.pi
            )
        )
        % 1.0
    ) * box_lengths


def unwrap_about_center(
    positions: np.ndarray,
    center: np.ndarray,
    box_lengths: np.ndarray,
) -> np.ndarray:
    return (
        center
        + minimum_image(
            positions
            - center,
            box_lengths,
        )
    )


def deterministic_axis(
    coordinates: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:
    centered = (
        coordinates
        - np.mean(
            coordinates,
            axis=0,
        )
    )

    covariance = np.cov(
        centered,
        rowvar=False,
        bias=True,
    )

    eigenvalues, eigenvectors = np.linalg.eigh(
        covariance
    )

    axis = eigenvectors[
        :,
        np.argmax(
            eigenvalues
        )
    ]

    dominant = int(
        np.argmax(
            np.abs(
                axis
            )
        )
    )

    if axis[
        dominant
    ] < 0.0:
        axis = -axis

    axis /= np.linalg.norm(
        axis
    )

    return (
        axis,
        eigenvalues,
    )


def perpendicular_basis(
    axis: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:
    axis = np.asarray(
        axis,
        dtype=float,
    )

    reference_candidates = (
        np.asarray(
            [1.0, 0.0, 0.0]
        ),
        np.asarray(
            [0.0, 1.0, 0.0]
        ),
        np.asarray(
            [0.0, 0.0, 1.0]
        ),
    )

    reference = min(
        reference_candidates,
        key=lambda vector: abs(
            float(
                np.dot(
                    axis,
                    vector,
                )
            )
        ),
    )

    first = np.cross(
        axis,
        reference,
    )

    first /= np.linalg.norm(
        first
    )

    second = np.cross(
        axis,
        first,
    )

    second /= np.linalg.norm(
        second
    )

    return (
        first,
        second,
    )


def triangular_disk(
    radius_nm: float,
    spacing_nm: float,
) -> np.ndarray:
    row_spacing = (
        spacing_nm
        * math.sqrt(
            3.0
        )
        / 2.0
    )

    maximum_row = int(
        math.ceil(
            radius_nm
            / row_spacing
        )
    ) + 2

    maximum_column = int(
        math.ceil(
            radius_nm
            / spacing_nm
        )
    ) + 3

    points = []

    for row in range(
        -maximum_row,
        maximum_row + 1,
    ):
        y = (
            row
            * row_spacing
        )

        shift = (
            0.5
            * spacing_nm
            if row % 2
            else 0.0
        )

        for column in range(
            -maximum_column,
            maximum_column + 1,
        ):
            x = (
                column
                * spacing_nm
                + shift
            )

            radius = math.hypot(
                x,
                y,
            )

            if radius <= (
                radius_nm
                + 1.0e-12
            ):
                points.append(
                    (
                        x,
                        y,
                    )
                )

    unique = sorted(
        {
            (
                round(
                    x,
                    12,
                ),
                round(
                    y,
                    12,
                ),
            )
            for x, y in points
        },
        key=lambda point: (
            round(
                math.hypot(
                    point[0],
                    point[1],
                ),
                12,
            ),
            round(
                math.atan2(
                    point[1],
                    point[0],
                ),
                12,
            ),
        ),
    )

    array = np.asarray(
        unique,
        dtype=float,
    )

    if len(
        array
    ) == 0:
        raise RuntimeError(
            "Generated cap disk contains no points."
        )

    return array


def coverage_hole_radius(
    disk_points: np.ndarray,
    accessible_radius_nm: float,
) -> float:
    values = np.arange(
        -accessible_radius_nm,
        accessible_radius_nm
        + 0.5
        * COVERAGE_GRID_SPACING_NM,
        COVERAGE_GRID_SPACING_NM,
    )

    grid_x, grid_y = np.meshgrid(
        values,
        values,
        indexing="xy",
    )

    mask = (
        grid_x
        * grid_x
        + grid_y
        * grid_y
        <= accessible_radius_nm
        * accessible_radius_nm
    )

    samples = np.column_stack(
        (
            grid_x[
                mask
            ],
            grid_y[
                mask
            ],
        )
    )

    minimum_distances = np.full(
        len(
            samples
        ),
        np.inf,
        dtype=float,
    )

    chunk_size = 2048

    for start in range(
        0,
        len(
            samples
        ),
        chunk_size,
    ):
        stop = min(
            start
            + chunk_size,
            len(
                samples
            ),
        )

        displacement = (
            samples[
                start:stop,
                None,
                :,
            ]
            - disk_points[
                None,
                :,
                :,
            ]
        )

        distances_squared = np.sum(
            displacement
            * displacement,
            axis=2,
        )

        minimum_distances[
            start:stop
        ] = np.sqrt(
            np.min(
                distances_squared,
                axis=1,
            )
        )

    return float(
        np.max(
            minimum_distances
        )
    )


def nearest_distances(
    query: np.ndarray,
    reference: np.ndarray,
    box_lengths: np.ndarray,
    query_chunk: int = 512,
) -> np.ndarray:
    if (
        len(
            query
        ) == 0
        or len(
            reference
        ) == 0
    ):
        return np.asarray(
            [],
            dtype=float,
        )

    minima = np.full(
        len(
            query
        ),
        np.inf,
        dtype=float,
    )

    for start in range(
        0,
        len(
            query
        ),
        query_chunk,
    ):
        stop = min(
            start
            + query_chunk,
            len(
                query
            ),
        )

        displacement = (
            query[
                start:stop,
                None,
                :,
            ]
            - reference[
                None,
                :,
                :,
            ]
        )

        displacement = minimum_image(
            displacement,
            box_lengths,
        )

        distances_squared = np.sum(
            displacement
            * displacement,
            axis=2,
        )

        minima[
            start:stop
        ] = np.sqrt(
            np.min(
                distances_squared,
                axis=1,
            )
        )

    return minima


def build_cap_positions(
    center: np.ndarray,
    axis: np.ndarray,
    basis_first: np.ndarray,
    basis_second: np.ndarray,
    axial_low: float,
    axial_high: float,
    offset_nm: float,
    disk_points: np.ndarray,
    box_lengths: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:
    lower_plane_center = (
        center
        + (
            axial_low
            - offset_nm
        )
        * axis
    )

    upper_plane_center = (
        center
        + (
            axial_high
            + offset_nm
        )
        * axis
    )

    in_plane = (
        np.outer(
            disk_points[
                :,
                0
            ],
            basis_first,
        )
        + np.outer(
            disk_points[
                :,
                1
            ],
            basis_second,
        )
    )

    lower = np.mod(
        lower_plane_center
        + in_plane,
        box_lengths,
    )

    upper = np.mod(
        upper_plane_center
        + in_plane,
        box_lengths,
    )

    return (
        lower,
        upper,
    )


def water_chunks(
    atoms: list[dict[str, Any]],
) -> list[list[dict[str, Any]]]:
    if len(
        atoms
    ) % WATER_SITES != 0:
        raise RuntimeError(
            "Water atom count is not divisible by four."
        )

    chunks = []

    for start in range(
        0,
        len(
            atoms
        ),
        WATER_SITES,
    ):
        chunk = atoms[
            start:
            start + WATER_SITES
        ]

        residue_identity = {
            (
                atom[
                    "resid"
                ],
                atom[
                    "resname"
                ],
            )
            for atom in chunk
        }

        if len(
            residue_identity
        ) != 1:
            raise RuntimeError(
                "Water ordering is not residue-consistent "
                f"at water index {start // WATER_SITES}."
            )

        chunks.append(
            chunk
        )

    return chunks


def oxygen_from_water_chunk(
    chunk: list[dict[str, Any]],
) -> dict[str, Any]:
    oxygen_candidates = [
        atom
        for atom in chunk
        if (
            atom[
                "atomname"
            ].upper().startswith(
                "O"
            )
            or atom[
                "atomname"
            ].upper()
            in {
                "OW",
                "OH2",
            }
        )
    ]

    if len(
        oxygen_candidates
    ) == 1:
        return oxygen_candidates[
            0
        ]

    return chunk[
        0
    ]


def gro_atom_line(
    resid: int,
    resname: str,
    atomname: str,
    atomnum: int,
    position: np.ndarray,
) -> str:
    return (
        f"{resid % 100000:5d}"
        f"{resname[:5]:<5s}"
        f"{atomname[:5]:>5s}"
        f"{atomnum % 100000:5d}"
        f"{position[0]:8.3f}"
        f"{position[1]:8.3f}"
        f"{position[2]:8.3f}"
    )


def write_gro(
    path: Path,
    title: str,
    atom_records: list[
        tuple[
            int,
            str,
            str,
            np.ndarray,
        ]
    ],
    box_lengths: np.ndarray,
) -> None:
    lines = [
        title,
        str(
            len(
                atom_records
            )
        ),
    ]

    for atom_number, (
        resid,
        resname,
        atomname,
        position,
    ) in enumerate(
        atom_records,
        start=1,
    ):
        lines.append(
            gro_atom_line(
                resid,
                resname,
                atomname,
                atom_number,
                position,
            )
        )

    lines.append(
        "".join(
            f"{value:10.5f}"
            for value in box_lengths
        )
    )

    path.write_text(
        "\n".join(
            lines
        )
        + "\n",
        encoding="utf-8",
    )


def write_csv(
    path: Path,
    rows: list[dict[str, Any]],
) -> None:
    if not rows:
        raise RuntimeError(
            f"No rows available for {path}"
        )

    fieldnames = list(
        rows[
            0
        ].keys()
    )

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
        writer.writerows(
            rows
        )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    require_file(
        REFERENCE_GRO
    )

    require_file(
        REFERENCE_SUMMARY
    )

    reference_summary = read_single_csv_row(
        REFERENCE_SUMMARY
    )

    if reference_summary.get(
        "authoritative_R1_start_state_accepted",
        "",
    ).strip().lower() not in {
        "true",
        "yes",
        "1",
    }:
        raise RuntimeError(
            "The R0 t=0 state is not authorized "
            "as the R1 starting state."
        )

    (
        reference_title,
        atoms,
        box,
    ) = read_gro(
        REFERENCE_GRO
    )

    if len(
        atoms
    ) != EXPECTED_ATOMS:
        raise RuntimeError(
            f"Expected {EXPECTED_ATOMS} atoms; "
            f"found {len(atoms)}"
        )

    box_lengths = orthorhombic_lengths(
        box
    )

    hbn_atoms = atoms[
        :HBN_ATOMS
    ]

    pyr_atoms = atoms[
        HBN_ATOMS:
        SOLUTE_ATOMS
    ]

    water_atoms = atoms[
        SOLUTE_ATOMS:
    ]

    waters = water_chunks(
        water_atoms
    )

    if len(
        waters
    ) != EXPECTED_WATERS:
        raise RuntimeError(
            f"Expected {EXPECTED_WATERS} waters; "
            f"found {len(waters)}"
        )

    hbn_positions_wrapped = np.asarray(
        [
            atom[
                "position"
            ]
            for atom in hbn_atoms
        ],
        dtype=float,
    )

    pyr_positions = np.asarray(
        [
            atom[
                "position"
            ]
            for atom in pyr_atoms
        ],
        dtype=float,
    )

    water_oxygen_positions = np.asarray(
        [
            oxygen_from_water_chunk(
                chunk
            )[
                "position"
            ]
            for chunk in waters
        ],
        dtype=float,
    )

    tube_center_wrapped = circular_center(
        hbn_positions_wrapped,
        box_lengths,
    )

    hbn_positions = unwrap_about_center(
        hbn_positions_wrapped,
        tube_center_wrapped,
        box_lengths,
    )

    tube_axis, _ = deterministic_axis(
        hbn_positions
    )

    basis_first, basis_second = perpendicular_basis(
        tube_axis
    )

    centered_hbn = (
        hbn_positions
        - tube_center_wrapped
    )

    hbn_axial = (
        centered_hbn
        @ tube_axis
    )

    hbn_perpendicular = (
        centered_hbn
        - np.outer(
            hbn_axial,
            tube_axis,
        )
    )

    hbn_radial = np.linalg.norm(
        hbn_perpendicular,
        axis=1,
    )

    axial_low = float(
        np.quantile(
            hbn_axial,
            0.01,
        )
    )

    axial_high = float(
        np.quantile(
            hbn_axial,
            0.99,
        )
    )

    wall_radius_q99 = float(
        np.quantile(
            hbn_radial,
            0.99,
        )
    )

    wall_radius_median = float(
        np.median(
            hbn_radial
        )
    )

    accessible_radius_nm = float(
        reference_summary[
            "provisional_accessible_radius_nm"
        ]
    )

    disk_radius_nm = (
        wall_radius_q99
        + CAP_EDGE_MARGIN_NM
    )

    disk_edge_margin_nm = (
        disk_radius_nm
        - wall_radius_q99
    )

    water_relative = minimum_image(
        water_oxygen_positions
        - tube_center_wrapped,
        box_lengths,
    )

    water_axial = (
        water_relative
        @ tube_axis
    )

    water_perpendicular = (
        water_relative
        - np.outer(
            water_axial,
            tube_axis,
        )
    )

    water_radial = np.linalg.norm(
        water_perpendicular,
        axis=1,
    )

    lumen_mask = (
        (
            water_axial
            >= axial_low
        )
        & (
            water_axial
            <= axial_high
        )
        & (
            water_radial
            <= accessible_radius_nm
        )
    )

    geometric_lumen_count = int(
        np.count_nonzero(
            lumen_mask
        )
    )

    if geometric_lumen_count != AUTHORITATIVE_LUMEN_WATERS:
        raise RuntimeError(
            "The prototype workflow did not reproduce "
            "the authoritative t=0 lumen occupancy: "
            f"{geometric_lumen_count}/"
            f"{AUTHORITATIVE_LUMEN_WATERS}"
        )

    solute_positions = np.vstack(
        (
            hbn_positions_wrapped,
            pyr_positions,
        )
    )

    candidate_rows: list[
        dict[str, Any]
    ] = []

    candidate_objects: list[
        dict[str, Any]
    ] = []

    for offset_nm in OFFSETS_NM:
        axial_steric_overlap_nm = (
            CAP_EFFECTIVE_WATER_EXCLUSION_NM
            + R0_END_EFFECTIVE_WATER_EXCLUSION_NM
            - offset_nm
        )

        for spacing_nm in LATTICE_SPACINGS_NM:
            disk_points = triangular_disk(
                disk_radius_nm,
                spacing_nm,
            )

            (
                lower_cap,
                upper_cap,
            ) = build_cap_positions(
                tube_center_wrapped,
                tube_axis,
                basis_first,
                basis_second,
                axial_low,
                axial_high,
                offset_nm,
                disk_points,
                box_lengths,
            )

            all_cap_positions = np.vstack(
                (
                    lower_cap,
                    upper_cap,
                )
            )

            coverage_hole_nm = coverage_hole_radius(
                disk_points,
                accessible_radius_nm,
            )

            cap_to_hbn = nearest_distances(
                all_cap_positions,
                hbn_positions_wrapped,
                box_lengths,
            )

            cap_to_pyr = nearest_distances(
                all_cap_positions,
                pyr_positions,
                box_lengths,
            )

            water_to_caps = nearest_distances(
                water_oxygen_positions,
                all_cap_positions,
                box_lengths,
            )

            remove_mask = (
                water_to_caps
                < CAP_WATER_PRUNING_CUTOFF_NM
            )

            removed_total = int(
                np.count_nonzero(
                    remove_mask
                )
            )

            removed_lumen = int(
                np.count_nonzero(
                    remove_mask
                    & lumen_mask
                )
            )

            retained_lumen = (
                geometric_lumen_count
                - removed_lumen
            )

            retained_lumen_fraction = (
                retained_lumen
                / geometric_lumen_count
            )

            removed_total_fraction = (
                removed_total
                / EXPECTED_WATERS
            )

            minimum_cap_hbn_nm = float(
                np.min(
                    cap_to_hbn
                )
            )

            minimum_cap_pyr_nm = float(
                np.min(
                    cap_to_pyr
                )
            )

            minimum_cap_water_nm = float(
                np.min(
                    water_to_caps
                )
            )

            gates = {
                "cap_HBN_distance": (
                    minimum_cap_hbn_nm
                    >= MIN_CAP_SOLUTE_DISTANCE_NM
                ),
                "cap_PYR_distance": (
                    minimum_cap_pyr_nm
                    >= MIN_CAP_SOLUTE_DISTANCE_NM
                ),
                "planar_coverage": (
                    coverage_hole_nm
                    <= MAX_COVERAGE_HOLE_NM
                ),
                "lumen_retention": (
                    retained_lumen_fraction
                    >= MIN_RETAINED_LUMEN_FRACTION
                ),
                "total_water_pruning": (
                    removed_total_fraction
                    <= MAX_REMOVED_TOTAL_WATER_FRACTION
                ),
                "axial_steric_overlap": (
                    axial_steric_overlap_nm
                    >= MIN_AXIAL_STERIC_OVERLAP_NM
                ),
                "disk_edge_margin": (
                    disk_edge_margin_nm
                    >= MIN_DISK_EDGE_MARGIN_NM
                ),
                "lower_upper_symmetry": (
                    len(
                        lower_cap
                    )
                    == len(
                        upper_cap
                    )
                ),
            }

            failed_gates = [
                name
                for name, passed
                in gates.items()
                if not passed
            ]

            passed = (
                len(
                    failed_gates
                )
                == 0
            )

            row = {
                "candidate_id": (
                    f"offset_{offset_nm:.2f}_"
                    f"spacing_{spacing_nm:.2f}"
                ),
                "offset_nm": offset_nm,
                "spacing_nm": spacing_nm,
                "disk_radius_nm": disk_radius_nm,
                "disk_edge_margin_nm": (
                    disk_edge_margin_nm
                ),
                "beads_per_cap": len(
                    disk_points
                ),
                "total_cap_beads": len(
                    all_cap_positions
                ),
                "coverage_hole_nm": (
                    coverage_hole_nm
                ),
                "axial_steric_overlap_nm": (
                    axial_steric_overlap_nm
                ),
                "minimum_cap_HBN_distance_nm": (
                    minimum_cap_hbn_nm
                ),
                "minimum_cap_PYR_distance_nm": (
                    minimum_cap_pyr_nm
                ),
                "minimum_cap_waterO_distance_before_pruning_nm": (
                    minimum_cap_water_nm
                ),
                "removed_total_water_molecules": (
                    removed_total
                ),
                "removed_total_water_fraction": (
                    removed_total_fraction
                ),
                "removed_lumen_water_molecules": (
                    removed_lumen
                ),
                "retained_lumen_water_molecules": (
                    retained_lumen
                ),
                "retained_lumen_fraction": (
                    retained_lumen_fraction
                ),
                "candidate_pass": passed,
                "failed_gates": (
                    " | ".join(
                        failed_gates
                    )
                ),
            }

            candidate_rows.append(
                row
            )

            candidate_objects.append(
                {
                    "row": row,
                    "disk_points": disk_points,
                    "lower_cap": lower_cap,
                    "upper_cap": upper_cap,
                    "all_cap_positions": (
                        all_cap_positions
                    ),
                    "remove_mask": remove_mask,
                    "water_to_caps": (
                        water_to_caps
                    ),
                }
            )

    write_csv(
        CANDIDATE_CSV,
        candidate_rows,
    )

    passing = [
        candidate
        for candidate in candidate_objects
        if bool(
            candidate[
                "row"
            ][
                "candidate_pass"
            ]
        )
    ]

    if not passing:
        REPORT_MD.write_text(
            """# R1 Steric Cap Prototype

No candidate passed all geometric gates.

Review `r1_steric_cap_candidate_scan.csv`.
No topology or MD execution is authorized.
""",
            encoding="utf-8",
        )

        raise RuntimeError(
            "No R1 steric-cap candidate passed."
        )

    selected = sorted(
        passing,
        key=lambda candidate: (
            -int(
                candidate[
                    "row"
                ][
                    "retained_lumen_water_molecules"
                ]
            ),
            int(
                candidate[
                    "row"
                ][
                    "removed_total_water_molecules"
                ]
            ),
            float(
                candidate[
                    "row"
                ][
                    "offset_nm"
                ]
            ),
            int(
                candidate[
                    "row"
                ][
                    "total_cap_beads"
                ]
            ),
            float(
                candidate[
                    "row"
                ][
                    "coverage_hole_nm"
                ]
            ),
        ),
    )[0]

    selected_row = dict(
        selected[
            "row"
        ]
    )

    selected_row[
        "selected"
    ] = True

    lower_cap = selected[
        "lower_cap"
    ]

    upper_cap = selected[
        "upper_cap"
    ]

    all_cap_positions = selected[
        "all_cap_positions"
    ]

    remove_mask = selected[
        "remove_mask"
    ]

    retained_waters = [
        chunk
        for water_index, chunk in enumerate(
            waters
        )
        if not bool(
            remove_mask[
                water_index
            ]
        )
    ]

    removed_water_rows = []

    lower_cap_distances = nearest_distances(
        water_oxygen_positions,
        lower_cap,
        box_lengths,
    )

    upper_cap_distances = nearest_distances(
        water_oxygen_positions,
        upper_cap,
        box_lengths,
    )

    for water_index, removed in enumerate(
        remove_mask
    ):
        if not removed:
            continue

        oxygen = oxygen_from_water_chunk(
            waters[
                water_index
            ]
        )

        lower_distance = float(
            lower_cap_distances[
                water_index
            ]
        )

        upper_distance = float(
            upper_cap_distances[
                water_index
            ]
        )

        if (
            lower_distance
            < CAP_WATER_PRUNING_CUTOFF_NM
            and upper_distance
            < CAP_WATER_PRUNING_CUTOFF_NM
        ):
            cap_side = "BOTH"
        elif (
            lower_distance
            < upper_distance
        ):
            cap_side = "LOWER"
        else:
            cap_side = "UPPER"

        removed_water_rows.append(
            {
                "original_water_index_zero_based": (
                    water_index
                ),
                "original_resid": oxygen[
                    "resid"
                ],
                "resname": oxygen[
                    "resname"
                ],
                "oxygen_atomname": oxygen[
                    "atomname"
                ],
                "oxygen_x_nm": oxygen[
                    "position"
                ][
                    0
                ],
                "oxygen_y_nm": oxygen[
                    "position"
                ][
                    1
                ],
                "oxygen_z_nm": oxygen[
                    "position"
                ][
                    2
                ],
                "minimum_lower_cap_distance_nm": (
                    lower_distance
                ),
                "minimum_upper_cap_distance_nm": (
                    upper_distance
                ),
                "assigned_cap_side": (
                    cap_side
                ),
                "was_authoritative_lumen_water": (
                    bool(
                        lumen_mask[
                            water_index
                        ]
                    )
                ),
            }
        )

    if removed_water_rows:
        write_csv(
            REMOVED_WATER_CSV,
            removed_water_rows,
        )
    else:
        REMOVED_WATER_CSV.write_text(
            (
                "original_water_index_zero_based,"
                "original_resid,resname,"
                "oxygen_atomname,oxygen_x_nm,"
                "oxygen_y_nm,oxygen_z_nm,"
                "minimum_lower_cap_distance_nm,"
                "minimum_upper_cap_distance_nm,"
                "assigned_cap_side,"
                "was_authoritative_lumen_water\n"
            ),
            encoding="utf-8",
        )

    cap_coordinate_rows = []

    for side, positions in (
        (
            "LOWER",
            lower_cap,
        ),
        (
            "UPPER",
            upper_cap,
        ),
    ):
        for local_index, position in enumerate(
            positions,
            start=1,
        ):
            cap_coordinate_rows.append(
                {
                    "side": side,
                    "local_atom_index": (
                        local_index
                    ),
                    "atomname": (
                        f"C{local_index:04d}"
                    ),
                    "x_nm": position[
                        0
                    ],
                    "y_nm": position[
                        1
                    ],
                    "z_nm": position[
                        2
                    ],
                }
            )

    write_csv(
        SELECTED_CAP_CSV,
        cap_coordinate_rows,
    )

    cap_only_records: list[
        tuple[
            int,
            str,
            str,
            np.ndarray,
        ]
    ] = []

    for side_index, (
        resname,
        positions,
    ) in enumerate(
        (
            (
                "CPL",
                lower_cap,
            ),
            (
                "CPU",
                upper_cap,
            ),
        ),
        start=1,
    ):
        for local_index, position in enumerate(
            positions,
            start=1,
        ):
            cap_only_records.append(
                (
                    side_index,
                    resname,
                    f"C{local_index:04d}",
                    position,
                )
            )

    write_gro(
        SELECTED_CAP_GRO,
        (
            "R1 neutral frozen steric cap geometry; "
            f"offset={selected_row['offset_nm']:.2f} nm; "
            f"spacing={selected_row['spacing_nm']:.2f} nm"
        ),
        cap_only_records,
        box_lengths,
    )

    full_records: list[
        tuple[
            int,
            str,
            str,
            np.ndarray,
        ]
    ] = []

    for atom in (
        hbn_atoms
        + pyr_atoms
    ):
        full_records.append(
            (
                atom[
                    "resid"
                ],
                atom[
                    "resname"
                ],
                atom[
                    "atomname"
                ],
                atom[
                    "position"
                ],
            )
        )

    for chunk in retained_waters:
        for atom in chunk:
            full_records.append(
                (
                    atom[
                        "resid"
                    ],
                    atom[
                        "resname"
                    ],
                    atom[
                        "atomname"
                    ],
                    atom[
                        "position"
                    ],
                )
            )

    maximum_existing_resid = max(
        atom[
            "resid"
        ]
        for atom in atoms
    )

    lower_cap_resid = (
        maximum_existing_resid
        + 1
    )

    upper_cap_resid = (
        maximum_existing_resid
        + 2
    )

    for local_index, position in enumerate(
        lower_cap,
        start=1,
    ):
        full_records.append(
            (
                lower_cap_resid,
                "CPL",
                f"C{local_index:04d}",
                position,
            )
        )

    for local_index, position in enumerate(
        upper_cap,
        start=1,
    ):
        full_records.append(
            (
                upper_cap_resid,
                "CPU",
                f"C{local_index:04d}",
                position,
            )
        )

    write_gro(
        SELECTED_SYSTEM_GRO,
        (
            "R1 geometry-only prototype derived from "
            "accepted R0 t=0; no force field assigned"
        ),
        full_records,
        box_lengths,
    )

    retained_water_count = len(
        retained_waters
    )

    final_atom_count = (
        SOLUTE_ATOMS
        + WATER_SITES
        * retained_water_count
        + len(
            all_cap_positions
        )
    )

    if len(
        full_records
    ) != final_atom_count:
        raise RuntimeError(
            "Derived atom-count accounting failed."
        )

    removed_lumen_count = int(
        selected_row[
            "removed_lumen_water_molecules"
        ]
    )

    retained_lumen_count = (
        AUTHORITATIVE_LUMEN_WATERS
        - removed_lumen_count
    )

    selected_definition = {
        "prototype_role": (
            "neutral_frozen_steric_positive_control"
        ),
        "chemical_architecture_claim": False,
        "source_R0_t0_gro": relative(
            REFERENCE_GRO
        ),
        "source_R0_t0_sha256": sha256_file(
            REFERENCE_GRO
        ),
        "selected_candidate_id": (
            selected_row[
                "candidate_id"
            ]
        ),
        "tube_axis": [
            float(value)
            for value in tube_axis
        ],
        "tube_center_wrapped_nm": [
            float(value)
            for value in tube_center_wrapped
        ],
        "axial_low_nm": axial_low,
        "axial_high_nm": axial_high,
        "wall_radius_median_nm": (
            wall_radius_median
        ),
        "wall_radius_q99_nm": (
            wall_radius_q99
        ),
        "accessible_radius_nm": (
            accessible_radius_nm
        ),
        "cap_offset_nm": float(
            selected_row[
                "offset_nm"
            ]
        ),
        "cap_lattice_spacing_nm": float(
            selected_row[
                "spacing_nm"
            ]
        ),
        "cap_disk_radius_nm": (
            disk_radius_nm
        ),
        "cap_beads_per_end": len(
            lower_cap
        ),
        "total_cap_beads": len(
            all_cap_positions
        ),
        "provisional_cap_water_exclusion_nm": (
            CAP_EFFECTIVE_WATER_EXCLUSION_NM
        ),
        "water_pruning_cutoff_nm": (
            CAP_WATER_PRUNING_CUTOFF_NM
        ),
        "initial_water_molecules": (
            EXPECTED_WATERS
        ),
        "removed_water_molecules": int(
            selected_row[
                "removed_total_water_molecules"
            ]
        ),
        "retained_water_molecules": (
            retained_water_count
        ),
        "initial_lumen_water_molecules": (
            AUTHORITATIVE_LUMEN_WATERS
        ),
        "removed_lumen_water_molecules": (
            removed_lumen_count
        ),
        "retained_lumen_water_molecules": (
            retained_lumen_count
        ),
        "derived_atom_count": (
            final_atom_count
        ),
        "topology_assigned": False,
        "MD_execution_authorized": False,
    }

    SELECTED_JSON.write_text(
        json.dumps(
            selected_definition,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    summary = {
        **selected_row,
        "prototype_role": (
            "NEUTRAL_FROZEN_STERIC_POSITIVE_CONTROL"
        ),
        "reference_R0_sha256": sha256_file(
            REFERENCE_GRO
        ),
        "caps_only_sha256": sha256_file(
            SELECTED_CAP_GRO
        ),
        "derived_system_sha256": sha256_file(
            SELECTED_SYSTEM_GRO
        ),
        "initial_system_atoms": (
            EXPECTED_ATOMS
        ),
        "initial_water_molecules": (
            EXPECTED_WATERS
        ),
        "retained_water_molecules": (
            retained_water_count
        ),
        "derived_system_atoms": (
            final_atom_count
        ),
        "selected_cap_beads_per_end": (
            len(
                lower_cap
            )
        ),
        "selected_total_cap_beads": (
            len(
                all_cap_positions
            )
        ),
        "retained_authoritative_lumen_waters": (
            retained_lumen_count
        ),
        "geometry_prototype_accepted": True,
        "topology_generation_authorized": True,
        "MD_execution_authorized": False,
        "required_next_step": (
            "DEFINE_R1_CAP_NONBONDED_MODEL_AND_BUILD_TOPOLOGY"
        ),
    }

    write_csv(
        SUMMARY_CSV,
        [
            summary
        ],
    )

    passing_count = len(
        passing
    )

    report_candidate_lines = "\n".join(
        (
            f"- `{row['candidate_id']}`: "
            f"{'PASS' if row['candidate_pass'] else 'FAIL'}; "
            f"beads/end={row['beads_per_cap']}; "
            f"coverage hole={row['coverage_hole_nm']:.4f} nm; "
            f"removed waters={row['removed_total_water_molecules']}; "
            f"retained lumen="
            f"{row['retained_lumen_water_molecules']}/"
            f"{AUTHORITATIVE_LUMEN_WATERS}"
        )
        for row in candidate_rows
    )

    REPORT_MD.write_text(
        f"""# R1 Steric Cap Prototype

## Scientific role

R1 is a **neutral, frozen steric positive control** designed to test
whether blocking the two axial exits prevents the progressive water
depletion observed in R0.

R1 is not yet a chemically realizable final device architecture.
No physical cap atom type, bonded model, or nonbonded parameter set has
been assigned at this stage.

## R0 reference

- Source:
  `{relative(REFERENCE_GRO)}`
- Source SHA256:
  `{sha256_file(REFERENCE_GRO)}`
- Initial lumen waters:
  **{AUTHORITATIVE_LUMEN_WATERS}**
- Tube axis:
  **({tube_axis[0]:.8f},
  {tube_axis[1]:.8f},
  {tube_axis[2]:.8f})**
- Robust axial planes:
  **{axial_low:.6f}/{axial_high:.6f} nm**
- Wall radius q99:
  **{wall_radius_q99:.6f} nm**
- Provisional accessible radius:
  **{accessible_radius_nm:.6f} nm**

## Candidate scan

- Axial offsets:
  **{', '.join(f'{value:.2f}' for value in OFFSETS_NM)} nm**
- Triangular-lattice spacings:
  **{', '.join(f'{value:.2f}' for value in LATTICE_SPACINGS_NM)} nm**
- Disk radius:
  **{disk_radius_nm:.6f} nm**
- Candidates evaluated:
  **{len(candidate_rows)}**
- Candidates passing:
  **{passing_count}**

{report_candidate_lines}

## Selected candidate

- Candidate:
  **{selected_row['candidate_id']}**
- Offset:
  **{selected_row['offset_nm']:.3f} nm**
- Lattice spacing:
  **{selected_row['spacing_nm']:.3f} nm**
- Beads per cap:
  **{len(lower_cap)}**
- Total cap beads:
  **{len(all_cap_positions)}**
- Maximum planar coverage hole:
  **{selected_row['coverage_hole_nm']:.6f} nm**
- Provisional axial steric overlap:
  **{selected_row['axial_steric_overlap_nm']:.6f} nm**
- Minimum cap-HBN distance:
  **{selected_row['minimum_cap_HBN_distance_nm']:.6f} nm**
- Minimum cap-PYR distance:
  **{selected_row['minimum_cap_PYR_distance_nm']:.6f} nm**

## Water preservation

- Initial waters:
  **{EXPECTED_WATERS}**
- Removed waters:
  **{selected_row['removed_total_water_molecules']}**
- Retained waters:
  **{retained_water_count}**
- Initial lumen waters:
  **{AUTHORITATIVE_LUMEN_WATERS}**
- Removed lumen waters:
  **{removed_lumen_count}**
- Retained lumen waters:
  **{retained_lumen_count}**
- Retained lumen fraction:
  **{selected_row['retained_lumen_fraction']:.6f}**

Only complete TIP4P/2005 water molecules were removed.

## Derived files

- Cap-only geometry:
  `{relative(SELECTED_CAP_GRO)}`
- Geometry-only full system:
  `{relative(SELECTED_SYSTEM_GRO)}`
- Removed-water audit:
  `{relative(REMOVED_WATER_CSV)}`
- Cap coordinate table:
  `{relative(SELECTED_CAP_CSV)}`
- Machine-readable definition:
  `{relative(SELECTED_JSON)}`

## Decision

- Geometry prototype accepted: **YES**
- Topology construction authorized: **YES**
- Energy minimization authorized: **NO**
- MD execution authorized: **NO**
- QM execution authorized: **NO**
- Required next step:
  `DEFINE_R1_CAP_NONBONDED_MODEL_AND_BUILD_TOPOLOGY`

Before simulation, the cap model must explicitly define:

1. zero net charge;
2. cap-cap exclusions or zero cap-cap interaction;
3. a water-oxygen steric interaction;
4. interactions with HBN and PYR;
5. frozen coordinate groups;
6. static energy and overlap validation.
""",
        encoding="utf-8",
    )

    print(
        "Day023 R1 steric-cap prototype generation completed."
    )

    print(
        "Candidates evaluated / passing: "
        f"{len(candidate_rows)}/{passing_count}"
    )

    print(
        "Selected candidate: "
        f"{selected_row['candidate_id']}"
    )

    print(
        "Selected offset / spacing / disk radius: "
        f"{selected_row['offset_nm']:.3f}/"
        f"{selected_row['spacing_nm']:.3f}/"
        f"{disk_radius_nm:.6f} nm"
    )

    print(
        "Cap beads per end / total: "
        f"{len(lower_cap)}/"
        f"{len(all_cap_positions)}"
    )

    print(
        "Coverage hole / axial steric overlap: "
        f"{selected_row['coverage_hole_nm']:.6f}/"
        f"{selected_row['axial_steric_overlap_nm']:.6f} nm"
    )

    print(
        "Minimum cap-HBN / cap-PYR distance: "
        f"{selected_row['minimum_cap_HBN_distance_nm']:.6f}/"
        f"{selected_row['minimum_cap_PYR_distance_nm']:.6f} nm"
    )

    print(
        "Initial / removed / retained waters: "
        f"{EXPECTED_WATERS}/"
        f"{selected_row['removed_total_water_molecules']}/"
        f"{retained_water_count}"
    )

    print(
        "Initial / removed / retained lumen waters: "
        f"{AUTHORITATIVE_LUMEN_WATERS}/"
        f"{removed_lumen_count}/"
        f"{retained_lumen_count}"
    )

    print(
        "Retained lumen fraction: "
        f"{selected_row['retained_lumen_fraction']:.6f}"
    )

    print(
        "Derived system atoms: "
        f"{final_atom_count}"
    )

    print(
        "Geometry prototype accepted: YES"
    )

    print(
        "Topology generation authorized: YES"
    )

    print(
        "MD execution authorized: NO"
    )

    print(
        "Required next step: "
        "DEFINE_R1_CAP_NONBONDED_MODEL_AND_BUILD_TOPOLOGY"
    )

    print(
        f"Wrote: {relative(CANDIDATE_CSV)}"
    )

    print(
        f"Wrote: {relative(SELECTED_CAP_GRO)}"
    )

    print(
        f"Wrote: {relative(SELECTED_SYSTEM_GRO)}"
    )

    print(
        f"Wrote: {relative(REMOVED_WATER_CSV)}"
    )

    print(
        f"Wrote: {relative(SELECTED_CAP_CSV)}"
    )

    print(
        f"Wrote: {relative(SELECTED_JSON)}"
    )

    print(
        f"Wrote: {relative(SUMMARY_CSV)}"
    )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )


if __name__ == "__main__":
    main()
