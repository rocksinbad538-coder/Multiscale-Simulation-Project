#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

DAY023_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design"
)

OUTPUT_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "01_r2_parent_rim_chemical_audit"
)

SELECTION_SUMMARY = (
    DAY023_ROOT
    / "20_r1_r2_architecture_comparison_and_selection"
    / "r1_r2_architecture_selection_summary.csv"
)

GEOMETRY_SUMMARY = (
    DAY023_ROOT
    / "12_r2_partial_cap_geometry_design"
    / "r2_partial_cap_geometry_summary.csv"
)

STATIC_SUMMARY = (
    DAY023_ROOT
    / "13_r2_topology_static_scan"
    / "r2_topology_static_scan_summary.csv"
)

SYSTEM_GRO = (
    DAY023_ROOT
    / "15_r2_frozen_solute_nvt_20ps_preparation"
    / "r2_frozen_solute_nvt_20ps_input.gro"
)

SYSTEM_TPR = (
    DAY023_ROOT
    / "15_r2_frozen_solute_nvt_20ps_preparation"
    / "r2_frozen_solute_nvt_20ps.tpr"
)

TPR_DUMP = (
    OUTPUT_ROOT
    / "r2_parent_system_tpr_dump.txt"
)

TPR_DUMP_STDERR = (
    OUTPUT_ROOT
    / "r2_parent_system_tpr_dump.stderr.log"
)

HBN_ATOMS_CSV = (
    OUTPUT_ROOT
    / "r2_parent_hbn_atoms.csv"
)

HBN_BONDS_CSV = (
    OUTPUT_ROOT
    / "r2_parent_hbn_geometry_derived_bonds.csv"
)

RIM_ATOMS_CSV = (
    OUTPUT_ROOT
    / "r2_parent_terminal_rim_atoms.csv"
)

END_SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_parent_terminal_end_summary.csv"
)

CONSTRAINTS_CSV = (
    OUTPUT_ROOT
    / "r2_chemical_end_rim_design_constraints.csv"
)

CANDIDATES_CSV = (
    OUTPUT_ROOT
    / "r2_preliminary_chemical_candidate_classes.csv"
)

SOURCE_MANIFEST_CSV = (
    OUTPUT_ROOT
    / "r2_parent_rim_source_manifest.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_parent_rim_chemical_audit_summary.csv"
)

GATES_CSV = (
    OUTPUT_ROOT
    / "r2_parent_rim_chemical_audit_gates.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_PARENT_RIM_AND_CHEMICAL_CONSTRAINTS_AUDIT_DAY024.md"
)

EXPECTED_SELECTION_DECISION = (
    "R2_SELECTED_AS_PRIMARY_PARTIAL_CAP_SCREENING_ARCHITECTURE"
)

EXPECTED_ATOMS = 68332
EXPECTED_HBN_ATOMS = 1680
EXPECTED_B_ATOMS = 840
EXPECTED_N_ATOMS = 840

EXPECTED_PYRENE_ATOMS = 104
EXPECTED_WATERS = 16565
WATER_SITES = 4
EXPECTED_WATER_ATOMS = (
    EXPECTED_WATERS
    * WATER_SITES
)

EXPECTED_CAPS = 288
EXPECTED_CAPS_PER_END = 144

EXPECTED_EDGE_ATOMS = 120
EXPECTED_EDGE_ATOMS_PER_END = 60
EXPECTED_EDGE_B_PER_END = 30
EXPECTED_EDGE_N_PER_END = 30

EXPECTED_INTERIOR_ATOMS = (
    EXPECTED_HBN_ATOMS
    - EXPECTED_EDGE_ATOMS
)

EXPECTED_GEOMETRY_BONDS = (
    (
        EXPECTED_EDGE_ATOMS
        * 2
        + EXPECTED_INTERIOR_ATOMS
        * 3
    )
    // 2
)

BOND_SEARCH_MIN_NM = 0.115
BOND_SEARCH_MAX_NM = 0.175

BOND_VALID_MIN_NM = 0.125
BOND_VALID_MAX_NM = 0.165

MAX_END_AXIAL_STD_NM = 0.020
MAX_END_RADIUS_DIFFERENCE_NM = 0.020

MIN_EDGE_ALTERNATION_FRACTION = 0.95

MAX_ZERO_CHARGE_TOLERANCE_E = 1.0e-8

PRIMARY_DECISION = (
    "R2_PARENT_RIM_AND_CHEMICAL_CONSTRAINTS_AUDITED"
)


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


def locate_gmx() -> str:
    executable = shutil.which("gmx")

    if executable:
        return executable

    fallback = Path(
        "/usr/local/gromacs/bin/gmx"
    )

    if fallback.exists():
        return str(fallback)

    raise RuntimeError(
        "Could not locate GROMACS."
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(
                1024 * 1024
            ),
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
                    field: row.get(
                        field,
                        "",
                    )
                    for field in fields
                }
            )


def parse_float(
    row: dict[str, str],
    key: str,
) -> float:
    try:
        value = float(
            row[key]
        )
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
            f"Non-finite field {key!r}"
        )

    return value


def read_gro(
    path: Path,
) -> tuple[
    list[dict[str, Any]],
    np.ndarray,
]:
    require_file(path)

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    if len(lines) < 3:
        raise RuntimeError(
            f"Malformed GRO file: {path}"
        )

    atom_count = int(
        lines[1].strip()
    )

    if len(lines) < atom_count + 3:
        raise RuntimeError(
            f"Incomplete GRO file: {path}"
        )

    atoms: list[dict[str, Any]] = []

    for index, line in enumerate(
        lines[
            2:
            2 + atom_count
        ]
    ):
        if len(line) < 44:
            raise RuntimeError(
                f"Malformed GRO atom line "
                f"{index + 1}"
            )

        atoms.append(
            {
                "global_index_0based": index,
                "global_index_1based": index + 1,
                "resid": int(
                    line[0:5]
                ),
                "resname": line[
                    5:10
                ].strip(),
                "atomname": line[
                    10:15
                ].strip(),
                "gro_atom_number": int(
                    line[15:20]
                ),
                "position_nm": np.asarray(
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

    box_values = [
        float(value)
        for value in lines[
            atom_count + 2
        ].split()
    ]

    if len(box_values) != 3:
        raise RuntimeError(
            "An orthorhombic GRO box is required."
        )

    box = np.asarray(
        box_values,
        dtype=float,
    )

    return atoms, box


def minimum_image(
    displacement: np.ndarray,
    box: np.ndarray,
) -> np.ndarray:
    return (
        displacement
        - box
        * np.round(
            displacement
            / box
        )
    )


def parse_hbn_tpr_atoms(
    dump_text: str,
) -> list[dict[str, Any]]:
    section_starts = list(
        re.finditer(
            r"(?m)^\s*moltype\s*"
            r"\(\s*\d+\s*\):\s*$",
            dump_text,
        )
    )

    hbn_section: str | None = None

    for index, match in enumerate(
        section_starts
    ):
        stop = (
            section_starts[
                index + 1
            ].start()
            if index + 1
            < len(section_starts)
            else len(dump_text)
        )

        section = dump_text[
            match.start():
            stop
        ]

        if re.search(
            r'(?m)^\s*name="HBN"\s*$',
            section,
        ):
            hbn_section = section
            break

    if hbn_section is None:
        raise RuntimeError(
            "Could not locate the HBN moltype "
            "in the TPR dump."
        )

    atom_pattern = re.compile(
        r"atom\[\s*(\d+)\s*\]\s*="
        r"\{([^}]]+)\}"
    )

    rows: list[dict[str, Any]] = []

    for match in atom_pattern.finditer(
        hbn_section
    ):
        local_index = int(
            match.group(1)
        )

        body = match.group(2)

        q_match = re.search(
            r"(?:^|,\s*)q=\s*"
            r"([-+0-9.eE]+)",
            body,
        )

        mass_match = re.search(
            r"(?:^|,\s*)m=\s*"
            r"([-+0-9.eE]+)",
            body,
        )

        atomic_number_match = re.search(
            r"atomnumber=\s*(\d+)",
            body,
        )

        if (
            q_match is None
            or mass_match is None
            or atomic_number_match is None
        ):
            raise RuntimeError(
                "Could not parse an HBN atom record "
                f"from the TPR dump: {match.group(0)}"
            )

        atomic_number = int(
            atomic_number_match.group(1)
        )

        if atomic_number == 5:
            element = "B"
        elif atomic_number == 7:
            element = "N"
        else:
            element = (
                f"Z{atomic_number}"
            )

        rows.append(
            {
                "local_index_0based": (
                    local_index
                ),
                "element": element,
                "atomic_number": (
                    atomic_number
                ),
                "mass_u": float(
                    mass_match.group(1)
                ),
                "charge_e": float(
                    q_match.group(1)
                ),
            }
        )

    rows.sort(
        key=lambda row: int(
            row[
                "local_index_0based"
            ]
        )
    )

    return rows


def determine_tube_axis(
    positions: np.ndarray,
) -> tuple[
    np.ndarray,
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

    covariance = np.cov(
        centered.T
    )

    eigenvalues, eigenvectors = (
        np.linalg.eigh(
            covariance
        )
    )

    order = np.argsort(
        eigenvalues
    )

    axis = eigenvectors[
        :,
        order[-1]
    ]

    axis /= np.linalg.norm(
        axis
    )

    return (
        center,
        axis,
        eigenvalues[
            order
        ],
    )


def perpendicular_basis(
    axis: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:
    reference = (
        np.asarray(
            [1.0, 0.0, 0.0]
        )
        if abs(axis[0]) < 0.85
        else np.asarray(
            [0.0, 1.0, 0.0]
        )
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

    return first, second


def build_geometry_bonds(
    positions: np.ndarray,
    elements: np.ndarray,
    box: np.ndarray,
) -> list[dict[str, Any]]:
    b_indices = np.flatnonzero(
        elements == "B"
    )

    n_indices = np.flatnonzero(
        elements == "N"
    )

    b_positions = positions[
        b_indices
    ]

    n_positions = positions[
        n_indices
    ]

    displacement = (
        b_positions[
            :,
            None,
            :,
        ]
        - n_positions[
            None,
            :,
            :,
        ]
    )

    displacement = minimum_image(
        displacement,
        box,
    )

    distances = np.linalg.norm(
        displacement,
        axis=2,
    )

    pair_rows = []

    matches = np.argwhere(
        (
            distances
            >= BOND_SEARCH_MIN_NM
        )
        & (
            distances
            <= BOND_SEARCH_MAX_NM
        )
    )

    for b_local, n_local in matches:
        first = int(
            b_indices[
                b_local
            ]
        )

        second = int(
            n_indices[
                n_local
            ]
        )

        pair_rows.append(
            {
                "bond_index": (
                    len(pair_rows) + 1
                ),
                "atom_i_local_0based": (
                    first
                ),
                "atom_j_local_0based": (
                    second
                ),
                "atom_i_local_1based": (
                    first + 1
                ),
                "atom_j_local_1based": (
                    second + 1
                ),
                "element_i": "B",
                "element_j": "N",
                "distance_nm": float(
                    distances[
                        b_local,
                        n_local
                    ]
                ),
            }
        )

    return pair_rows


def end_geometry_metrics(
    label: str,
    indices: np.ndarray,
    positions: np.ndarray,
    elements: np.ndarray,
    center: np.ndarray,
    axis: np.ndarray,
    basis_1: np.ndarray,
    basis_2: np.ndarray,
) -> tuple[
    dict[str, Any],
    list[dict[str, Any]],
]:
    if len(indices) == 0:
        raise RuntimeError(
            f"No atoms found for {label}."
        )

    end_positions = positions[
        indices
    ]

    end_center = np.mean(
        end_positions,
        axis=0,
    )

    relative = (
        end_positions
        - end_center
    )

    axial = (
        relative
        @ axis
    )

    x_plane = (
        relative
        @ basis_1
    )

    y_plane = (
        relative
        @ basis_2
    )

    radius = np.sqrt(
        x_plane * x_plane
        + y_plane * y_plane
    )

    angles = np.arctan2(
        y_plane,
        x_plane,
    )

    order = np.argsort(
        angles
    )

    ordered_indices = (
        indices[
            order
        ]
    )

    ordered_angles = (
        angles[
            order
        ]
    )

    ordered_elements = (
        elements[
            ordered_indices
        ]
    )

    wrapped_angles = np.concatenate(
        (
            ordered_angles,
            [
                ordered_angles[0]
                + 2.0 * math.pi
            ],
        )
    )

    angular_gaps = np.diff(
        wrapped_angles
    )

    alternation = [
        ordered_elements[index]
        != ordered_elements[
            (index + 1)
            % len(
                ordered_elements
            )
        ]
        for index in range(
            len(
                ordered_elements
            )
        )
    ]

    alternation_fraction = float(
        np.mean(
            alternation
        )
    )

    atom_rows = []

    for sequence_index, local_index in enumerate(
        ordered_indices
    ):
        atom_rows.append(
            {
                "end": label,
                "circumferential_order": (
                    sequence_index
                ),
                "hbn_local_index_0based": int(
                    local_index
                ),
                "hbn_local_index_1based": int(
                    local_index
                    + 1
                ),
                "element": str(
                    elements[
                        local_index
                    ]
                ),
                "x_nm": float(
                    positions[
                        local_index,
                        0
                    ]
                ),
                "y_nm": float(
                    positions[
                        local_index,
                        1
                    ]
                ),
                "z_nm": float(
                    positions[
                        local_index,
                        2
                    ]
                ),
                "axial_coordinate_nm": float(
                    (
                        positions[
                            local_index
                        ]
                        - center
                    )
                    @ axis
                ),
                "radius_from_end_center_nm": float(
                    radius[
                        order[
                            sequence_index
                        ]
                    ]
                ),
                "circumferential_angle_rad": float(
                    ordered_angles[
                        sequence_index
                    ]
                ),
            }
        )

    metrics = {
        "end": label,
        "atom_count": len(indices),
        "B_count": int(
            np.count_nonzero(
                elements[
                    indices
                ]
                == "B"
            )
        ),
        "N_count": int(
            np.count_nonzero(
                elements[
                    indices
                ]
                == "N"
            )
        ),
        "center_x_nm": float(
            end_center[0]
        ),
        "center_y_nm": float(
            end_center[1]
        ),
        "center_z_nm": float(
            end_center[2]
        ),
        "center_axial_coordinate_nm": float(
            (
                end_center
                - center
            )
            @ axis
        ),
        "axial_standard_deviation_nm": float(
            np.std(
                axial
            )
        ),
        "radius_mean_nm": float(
            np.mean(
                radius
            )
        ),
        "radius_standard_deviation_nm": float(
            np.std(
                radius
            )
        ),
        "radius_minimum_nm": float(
            np.min(
                radius
            )
        ),
        "radius_maximum_nm": float(
            np.max(
                radius
            )
        ),
        "angular_gap_mean_rad": float(
            np.mean(
                angular_gaps
            )
        ),
        "angular_gap_standard_deviation_rad": float(
            np.std(
                angular_gaps
            )
        ),
        "angular_gap_maximum_rad": float(
            np.max(
                angular_gaps
            )
        ),
        "element_alternation_fraction": (
            alternation_fraction
        ),
    }

    return metrics, atom_rows


def split_cap_groups(
    cap_positions: np.ndarray,
    center: np.ndarray,
    axis: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:
    projections = (
        cap_positions
        - center
    ) @ axis

    order = np.argsort(
        projections
    )

    half = len(
        cap_positions
    ) // 2

    return (
        order[:half],
        order[half:],
    )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        SELECTION_SUMMARY,
        GEOMETRY_SUMMARY,
        STATIC_SUMMARY,
        SYSTEM_GRO,
        SYSTEM_TPR,
    ):
        require_file(required)

    selection = read_single_csv_row(
        SELECTION_SUMMARY
    )

    geometry_summary = (
        read_single_csv_row(
            GEOMETRY_SUMMARY
        )
    )

    static_summary = (
        read_single_csv_row(
            STATIC_SUMMARY
        )
    )

    if (
        selection.get(
            "decision",
            "",
        )
        != EXPECTED_SELECTION_DECISION
    ):
        raise RuntimeError(
            "R2 is not in the accepted "
            "architecture-selection state."
        )

    target_aperture_diameter_nm = (
        parse_float(
            selection,
            "R2_effective_aperture_diameter_nm",
        )
    )

    target_aperture_radius_nm = (
        target_aperture_diameter_nm
        / 2.0
    )

    target_open_area_fraction = (
        parse_float(
            selection,
            "R2_open_area_fraction",
        )
    )

    validated_minimum_cap_ow_nm = (
        parse_float(
            selection,
            "R2_minimum_CAP_OW_distance_nm",
        )
    )

    atoms, box = read_gro(
        SYSTEM_GRO
    )

    if len(atoms) != EXPECTED_ATOMS:
        raise RuntimeError(
            "Unexpected system atom count: "
            f"{len(atoms)}/{EXPECTED_ATOMS}"
        )

    hbn_global_indices = np.asarray(
        [
            atom[
                "global_index_0based"
            ]
            for atom in atoms
            if atom[
                "resname"
            ].upper()
            == "HBN"
        ],
        dtype=int,
    )

    pyr_global_indices = np.asarray(
        [
            atom[
                "global_index_0based"
            ]
            for atom in atoms
            if atom[
                "resname"
            ].upper().startswith(
                "PYR"
            )
        ],
        dtype=int,
    )

    water_global_indices = np.asarray(
        [
            atom[
                "global_index_0based"
            ]
            for atom in atoms
            if atom[
                "resname"
            ].upper()
            == "SOL"
        ],
        dtype=int,
    )

    cap_global_indices = np.asarray(
        [
            atom[
                "global_index_0based"
            ]
            for atom in atoms
            if (
                atom[
                    "resname"
                ].upper().startswith(
                    "CAP"
                )
                or atom[
                    "atomname"
                ].upper()
                == "CAP"
            )
        ],
        dtype=int,
    )

    if len(hbn_global_indices) != EXPECTED_HBN_ATOMS:
        raise RuntimeError(
            "Unexpected HBN atom count: "
            f"{len(hbn_global_indices)}/"
            f"{EXPECTED_HBN_ATOMS}"
        )

    if not np.array_equal(
        hbn_global_indices,
        np.arange(
            EXPECTED_HBN_ATOMS,
            dtype=int,
        ),
    ):
        raise RuntimeError(
            "HBN atoms are not the first contiguous "
            "1680 atoms in the current system."
        )

    gmx = locate_gmx()

    dump = subprocess.run(
        [
            gmx,
            "dump",
            "-s",
            str(SYSTEM_TPR),
        ],
        cwd=OUTPUT_ROOT,
        env={
            **os.environ,
            "GMX_MAXBACKUP": "-1",
        },
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )

    TPR_DUMP.write_text(
        dump.stdout,
        encoding="utf-8",
    )

    TPR_DUMP_STDERR.write_text(
        dump.stderr,
        encoding="utf-8",
    )

    if dump.returncode != 0:
        raise RuntimeError(
            "Could not dump the accepted R2 TPR."
        )

    hbn_tpr_atoms = (
        parse_hbn_tpr_atoms(
            dump.stdout
        )
    )

    if len(hbn_tpr_atoms) != EXPECTED_HBN_ATOMS:
        raise RuntimeError(
            "Unexpected HBN atom count in TPR dump: "
            f"{len(hbn_tpr_atoms)}/"
            f"{EXPECTED_HBN_ATOMS}"
        )

    elements = np.asarray(
        [
            row["element"]
            for row in hbn_tpr_atoms
        ],
        dtype=object,
    )

    charges = np.asarray(
        [
            float(
                row[
                    "charge_e"
                ]
            )
            for row in hbn_tpr_atoms
        ],
        dtype=float,
    )

    atomic_numbers = np.asarray(
        [
            int(
                row[
                    "atomic_number"
                ]
            )
            for row in hbn_tpr_atoms
        ],
        dtype=int,
    )

    hbn_positions = np.asarray(
        [
            atoms[index][
                "position_nm"
            ]
            for index in hbn_global_indices
        ],
        dtype=float,
    )

    (
        tube_center,
        tube_axis,
        pca_eigenvalues,
    ) = determine_tube_axis(
        hbn_positions
    )

    basis_1, basis_2 = (
        perpendicular_basis(
            tube_axis
        )
    )

    bond_rows = build_geometry_bonds(
        hbn_positions,
        elements,
        box,
    )

    degrees = np.zeros(
        EXPECTED_HBN_ATOMS,
        dtype=int,
    )

    adjacency: list[
        list[int]
    ] = [
        []
        for _ in range(
            EXPECTED_HBN_ATOMS
        )
    ]

    for row in bond_rows:
        first = int(
            row[
                "atom_i_local_0based"
            ]
        )

        second = int(
            row[
                "atom_j_local_0based"
            ]
        )

        degrees[first] += 1
        degrees[second] += 1

        adjacency[first].append(
            second
        )

        adjacency[second].append(
            first
        )

    bond_distances = np.asarray(
        [
            float(
                row[
                    "distance_nm"
                ]
            )
            for row in bond_rows
        ],
        dtype=float,
    )

    edge_indices = np.flatnonzero(
        degrees == 2
    )

    interior_indices = np.flatnonzero(
        degrees == 3
    )

    anomalous_indices = np.flatnonzero(
        (degrees < 2)
        | (degrees > 3)
    )

    edge_projections = (
        hbn_positions[
            edge_indices
        ]
        - tube_center
    ) @ tube_axis

    edge_order = np.argsort(
        edge_projections
    )

    lower_edge_indices = edge_indices[
        edge_order[
            :
            len(edge_order) // 2
        ]
    ]

    upper_edge_indices = edge_indices[
        edge_order[
            len(edge_order) // 2:
        ]
    ]

    lower_metrics, lower_rows = (
        end_geometry_metrics(
            "LOWER",
            lower_edge_indices,
            hbn_positions,
            elements,
            tube_center,
            tube_axis,
            basis_1,
            basis_2,
        )
    )

    upper_metrics, upper_rows = (
        end_geometry_metrics(
            "UPPER",
            upper_edge_indices,
            hbn_positions,
            elements,
            tube_center,
            tube_axis,
            basis_1,
            basis_2,
        )
    )

    rim_rows = (
        lower_rows
        + upper_rows
    )

    cap_positions = np.asarray(
        [
            atoms[index][
                "position_nm"
            ]
            for index in cap_global_indices
        ],
        dtype=float,
    )

    (
        lower_cap_local,
        upper_cap_local,
    ) = split_cap_groups(
        cap_positions,
        tube_center,
        tube_axis,
    )

    lower_cap_positions = (
        cap_positions[
            lower_cap_local
        ]
    )

    upper_cap_positions = (
        cap_positions[
            upper_cap_local
        ]
    )

    lower_cap_axial = float(
        np.mean(
            (
                lower_cap_positions
                - tube_center
            )
            @ tube_axis
        )
    )

    upper_cap_axial = float(
        np.mean(
            (
                upper_cap_positions
                - tube_center
            )
            @ tube_axis
        )
    )

    lower_cap_offset_nm = abs(
        lower_cap_axial
        - float(
            lower_metrics[
                "center_axial_coordinate_nm"
            ]
        )
    )

    upper_cap_offset_nm = abs(
        upper_cap_axial
        - float(
            upper_metrics[
                "center_axial_coordinate_nm"
            ]
        )
    )

    lower_radius_nm = float(
        lower_metrics[
            "radius_mean_nm"
        ]
    )

    upper_radius_nm = float(
        upper_metrics[
            "radius_mean_nm"
        ]
    )

    mean_parent_rim_radius_nm = (
        0.5
        * (
            lower_radius_nm
            + upper_radius_nm
        )
    )

    required_radial_occlusion_nm = (
        mean_parent_rim_radius_nm
        - target_aperture_radius_nm
    )

    median_bn_bond_nm = float(
        np.median(
            bond_distances
        )
    )

    mean_bn_bond_nm = float(
        np.mean(
            bond_distances
        )
    )

    annulus_area_nm2 = (
        math.pi
        * (
            mean_parent_rim_radius_nm
            * mean_parent_rim_radius_nm
            - target_aperture_radius_nm
            * target_aperture_radius_nm
        )
    )

    hbn_area_per_atom_nm2 = (
        3.0
        * math.sqrt(3.0)
        / 4.0
        * median_bn_bond_nm
        * median_bn_bond_nm
    )

    estimated_monolayer_annulus_atoms = (
        annulus_area_nm2
        / hbn_area_per_atom_nm2
    )

    edge_edge_bonds = 0
    edge_interior_bonds = 0

    edge_set = set(
        int(index)
        for index in edge_indices
    )

    for row in bond_rows:
        first = int(
            row[
                "atom_i_local_0based"
            ]
        )

        second = int(
            row[
                "atom_j_local_0based"
            ]
        )

        first_edge = (
            first in edge_set
        )

        second_edge = (
            second in edge_set
        )

        if first_edge and second_edge:
            edge_edge_bonds += 1
        elif first_edge or second_edge:
            edge_interior_bonds += 1

    atom_rows = []

    end_label_by_index = {
        int(row[
            "hbn_local_index_0based"
        ]): row["end"]
        for row in rim_rows
    }

    for local_index in range(
        EXPECTED_HBN_ATOMS
    ):
        position = hbn_positions[
            local_index
        ]

        atom_rows.append(
            {
                "hbn_local_index_0based": (
                    local_index
                ),
                "hbn_local_index_1based": (
                    local_index + 1
                ),
                "system_global_index_0based": (
                    int(
                        hbn_global_indices[
                            local_index
                        ]
                    )
                ),
                "system_global_index_1based": (
                    int(
                        hbn_global_indices[
                            local_index
                        ]
                    )
                    + 1
                ),
                "element": str(
                    elements[
                        local_index
                    ]
                ),
                "atomic_number": int(
                    atomic_numbers[
                        local_index
                    ]
                ),
                "charge_e": float(
                    charges[
                        local_index
                    ]
                ),
                "x_nm": float(
                    position[0]
                ),
                "y_nm": float(
                    position[1]
                ),
                "z_nm": float(
                    position[2]
                ),
                "coordination_number": int(
                    degrees[
                        local_index
                    ]
                ),
                "is_terminal_rim_atom": (
                    local_index
                    in edge_set
                ),
                "terminal_end": (
                    end_label_by_index.get(
                        local_index,
                        "",
                    )
                ),
                "neighbor_local_indices_1based": (
                    " ".join(
                        str(
                            neighbor + 1
                        )
                        for neighbor
                        in sorted(
                            adjacency[
                                local_index
                            ]
                        )
                    )
                ),
            }
        )

    write_csv(
        HBN_ATOMS_CSV,
        atom_rows,
    )

    write_csv(
        HBN_BONDS_CSV,
        bond_rows,
    )

    write_csv(
        RIM_ATOMS_CSV,
        rim_rows,
    )

    write_csv(
        END_SUMMARY_CSV,
        [
            lower_metrics,
            upper_metrics,
        ],
    )

    design_constraints = [
        {
            "constraint": (
                "parent_HBN_atoms"
            ),
            "value": EXPECTED_HBN_ATOMS,
            "unit": "atoms",
            "basis": (
                "Accepted R2 parent scaffold"
            ),
            "status": "FIXED",
        },
        {
            "constraint": (
                "terminal_anchor_atoms_per_end"
            ),
            "value": (
                EXPECTED_EDGE_ATOMS_PER_END
            ),
            "unit": "atoms/end",
            "basis": (
                "Geometry-derived degree-2 rim"
            ),
            "status": "FIXED",
        },
        {
            "constraint": (
                "terminal_B_anchor_atoms_per_end"
            ),
            "value": (
                EXPECTED_EDGE_B_PER_END
            ),
            "unit": "atoms/end",
            "basis": (
                "Balanced armchair-like BN rim"
            ),
            "status": "FIXED",
        },
        {
            "constraint": (
                "terminal_N_anchor_atoms_per_end"
            ),
            "value": (
                EXPECTED_EDGE_N_PER_END
            ),
            "unit": "atoms/end",
            "basis": (
                "Balanced armchair-like BN rim"
            ),
            "status": "FIXED",
        },
        {
            "constraint": (
                "target_effective_aperture_diameter"
            ),
            "value": (
                target_aperture_diameter_nm
            ),
            "unit": "nm",
            "basis": (
                "Validated R2 screening geometry"
            ),
            "status": "TARGET",
        },
        {
            "constraint": (
                "target_effective_aperture_radius"
            ),
            "value": (
                target_aperture_radius_nm
            ),
            "unit": "nm",
            "basis": (
                "Half of validated aperture diameter"
            ),
            "status": "TARGET",
        },
        {
            "constraint": (
                "validated_open_area_fraction"
            ),
            "value": (
                target_open_area_fraction
            ),
            "unit": "dimensionless",
            "basis": (
                "Validated R2 screening geometry"
            ),
            "status": "TARGET",
        },
        {
            "constraint": (
                "parent_rim_mean_radius"
            ),
            "value": (
                mean_parent_rim_radius_nm
            ),
            "unit": "nm",
            "basis": (
                "Mean of both terminal rims"
            ),
            "status": "MEASURED",
        },
        {
            "constraint": (
                "required_radial_occlusion"
            ),
            "value": (
                required_radial_occlusion_nm
            ),
            "unit": "nm",
            "basis": (
                "Parent rim radius minus target "
                "aperture radius"
            ),
            "status": "REQUIRED",
        },
        {
            "constraint": (
                "annular_area_to_replace_per_end"
            ),
            "value": (
                annulus_area_nm2
            ),
            "unit": "nm^2",
            "basis": (
                "Geometric annulus between parent "
                "rim and target pore"
            ),
            "status": "MEASURED",
        },
        {
            "constraint": (
                "estimated_hBN_monolayer_atoms_per_end"
            ),
            "value": (
                estimated_monolayer_annulus_atoms
            ),
            "unit": "atoms/end",
            "basis": (
                "Annulus area divided by h-BN "
                "area per atom"
            ),
            "status": "SCREENING_ESTIMATE",
        },
        {
            "constraint": (
                "validated_steric_beads_per_end"
            ),
            "value": (
                EXPECTED_CAPS_PER_END
            ),
            "unit": "beads/end",
            "basis": (
                "Accepted R2 screening model"
            ),
            "status": "REFERENCE",
        },
        {
            "constraint": (
                "minimum_validated_CAP_OW_distance"
            ),
            "value": (
                validated_minimum_cap_ow_nm
            ),
            "unit": "nm",
            "basis": (
                "Validated R2 50 ps trajectory"
            ),
            "status": "LOWER_REFERENCE",
        },
        {
            "constraint": (
                "initial_candidate_total_charge"
            ),
            "value": 0,
            "unit": "e",
            "basis": (
                "Minimal-perturbation first candidate"
            ),
            "status": "DESIGN_TARGET",
        },
        {
            "constraint": (
                "unpaired_electrons_in_initial_candidate"
            ),
            "value": 0,
            "unit": "electrons",
            "basis": (
                "No radical or spin-active chemistry "
                "at this gate"
            ),
            "status": "DESIGN_TARGET",
        },
        {
            "constraint": (
                "new_MD_or_QM_execution"
            ),
            "value": 0,
            "unit": "authorized calculations",
            "basis": (
                "Static chemical-design gate only"
            ),
            "status": "BLOCKED",
        },
    ]

    write_csv(
        CONSTRAINTS_CSV,
        design_constraints,
    )

    candidate_rows = [
        {
            "candidate_id": (
                "C0_SIMPLE_EDGE_PASSIVATION_ONLY"
            ),
            "structural_concept": (
                "Passivate the existing BNNT rim "
                "using only terminal H, OH, or NHx groups."
            ),
            "chemical_basis": (
                "BN edge passivation is chemically "
                "documented."
            ),
            "geometric_assessment": (
                "Insufficient as a direct R2-cap "
                "replacement because the required "
                f"inward radial closure is "
                f"{required_radial_occlusion_nm:.6f} nm."
            ),
            "priority": 0,
            "status": (
                "REJECT_AS_STANDALONE_CAP_REPLACEMENT"
            ),
            "principal_risk": (
                "Leaves a pore much wider than the "
                "validated 0.84 nm aperture."
            ),
            "literature_basis": (
                "10.1021/jp805790s; "
                "10.1039/C9NA00530G"
            ),
        },
        {
            "candidate_id": (
                "C1_LINKED_BN_ANNULAR_NANOFLAKE"
            ),
            "structural_concept": (
                "Coaxial atomically thin BN annular "
                "nanoflake placed at the validated cap "
                "plane and attached to the terminal rim "
                "through a chemically defined linker "
                "network."
            ),
            "chemical_basis": (
                "Preserves the BN material family and "
                "permits H or mixed edge passivation."
            ),
            "geometric_assessment": (
                "Closest atomistic analogue to the "
                "validated planar annular steric cap; "
                f"estimated monolayer population "
                f"{estimated_monolayer_annulus_atoms:.1f} "
                "atoms/end."
            ),
            "priority": 1,
            "status": (
                "ADVANCE_TO_EXPLICIT_GEOMETRY_AND_"
                "JUNCTION_DESIGN"
            ),
            "principal_risk": (
                "Outer-rim linker chemistry and junction "
                "strain must be resolved; a direct "
                "90-degree seamless sp2 seam is not "
                "assumed."
            ),
            "literature_basis": (
                "10.1021/jp805790s; "
                "10.1039/C9NA00530G"
            ),
        },
        {
            "candidate_id": (
                "C2_RIGID_ORGANIC_OR_HYBRID_MACROCYCLE"
            ),
            "structural_concept": (
                "Rigid aromatic or heteroaromatic annulus "
                "with an approximately 0.84 nm central "
                "pore, multiply anchored to the BNNT rim."
            ),
            "chemical_basis": (
                "Covalent amine and organic "
                "functionalization of BNNTs is documented."
            ),
            "geometric_assessment": (
                "Potentially synthetically modular, but "
                "requires explicit anchor placement and "
                "force-field parameterization."
            ),
            "priority": 2,
            "status": (
                "ADVANCE_AS_FALLBACK_CANDIDATE_CLASS"
            ),
            "principal_risk": (
                "Greater dielectric, dipolar, and "
                "electronic perturbation near the lumen."
            ),
            "literature_basis": (
                "10.1021/ja063653+"
            ),
        },
        {
            "candidate_id": (
                "C3_INWARD_TETHERED_FUNCTIONAL_CORONA"
            ),
            "structural_concept": (
                "Multiple inward-pointing tethered "
                "aromatic or polar groups attached around "
                "the terminal BN rim."
            ),
            "chemical_basis": (
                "Uses established covalent BNNT "
                "functionalization concepts."
            ),
            "geometric_assessment": (
                "Can reduce the effective aperture but "
                "does not naturally reproduce the rigid "
                "planar exclusion surface."
            ),
            "priority": 3,
            "status": (
                "DEFER_UNLESS_RIGID_ANNULAR_CANDIDATES_FAIL"
            ),
            "principal_risk": (
                "Flexibility, pore intermittency, "
                "collapse, and strong hydration dependence."
            ),
            "literature_basis": (
                "10.1021/ja063653+"
            ),
        },
        {
            "candidate_id": (
                "C4_METAL_OR_SPIN_ACTIVE_CAP"
            ),
            "structural_concept": (
                "Metal-coordinated or radical-containing "
                "terminal annulus."
            ),
            "chemical_basis": (
                "Potentially realizable in a broader "
                "device architecture."
            ),
            "geometric_assessment": (
                "Outside the current minimal-perturbation "
                "confinement gate."
            ),
            "priority": 4,
            "status": (
                "REJECT_AT_CURRENT_GATE"
            ),
            "principal_risk": (
                "Introduces new electronic and spin-active "
                "degrees of freedom before they are "
                "physically justified."
            ),
            "literature_basis": (
                "Not required for current R2 gate"
            ),
        },
    ]

    write_csv(
        CANDIDATES_CSV,
        candidate_rows,
    )

    lower_end_radius = float(
        lower_metrics[
            "radius_mean_nm"
        ]
    )

    upper_end_radius = float(
        upper_metrics[
            "radius_mean_nm"
        ]
    )

    b_count = int(
        np.count_nonzero(
            elements == "B"
        )
    )

    n_count = int(
        np.count_nonzero(
            elements == "N"
        )
    )

    zero_charge_count = int(
        np.count_nonzero(
            np.abs(
                charges
            )
            <= MAX_ZERO_CHARGE_TOLERANCE_E
        )
    )

    gates = {
        "R2_architecture_selection_is_valid": (
            selection.get(
                "decision",
                "",
            )
            == EXPECTED_SELECTION_DECISION
        ),
        "system_has_68332_atoms": (
            len(atoms)
            == EXPECTED_ATOMS
        ),
        "HBN_has_1680_atoms": (
            len(
                hbn_global_indices
            )
            == EXPECTED_HBN_ATOMS
        ),
        "HBN_atoms_are_first_and_contiguous": (
            np.array_equal(
                hbn_global_indices,
                np.arange(
                    EXPECTED_HBN_ATOMS,
                    dtype=int,
                ),
            )
        ),
        "pyrene_atom_count_is_104": (
            len(
                pyr_global_indices
            )
            == EXPECTED_PYRENE_ATOMS
        ),
        "water_atom_count_is_66260": (
            len(
                water_global_indices
            )
            == EXPECTED_WATER_ATOMS
        ),
        "cap_atom_count_is_288": (
            len(
                cap_global_indices
            )
            == EXPECTED_CAPS
        ),
        "cap_split_is_144_per_end": (
            len(
                lower_cap_local
            )
            == EXPECTED_CAPS_PER_END
            and len(
                upper_cap_local
            )
            == EXPECTED_CAPS_PER_END
        ),
        "TPR_dump_return_code_zero": (
            dump.returncode == 0
        ),
        "TPR_HBN_atom_count_is_1680": (
            len(
                hbn_tpr_atoms
            )
            == EXPECTED_HBN_ATOMS
        ),
        "HBN_contains_840_B_and_840_N": (
            b_count
            == EXPECTED_B_ATOMS
            and n_count
            == EXPECTED_N_ATOMS
        ),
        "HBN_contains_only_B_and_N": (
            set(
                elements.tolist()
            )
            == {
                "B",
                "N",
            }
        ),
        "current_HBN_screening_charges_are_zero": (
            zero_charge_count
            == EXPECTED_HBN_ATOMS
        ),
        "geometry_bond_count_is_expected": (
            len(
                bond_rows
            )
            == EXPECTED_GEOMETRY_BONDS
        ),
        "geometry_has_120_degree2_edge_atoms": (
            len(
                edge_indices
            )
            == EXPECTED_EDGE_ATOMS
        ),
        "geometry_has_1560_degree3_interior_atoms": (
            len(
                interior_indices
            )
            == EXPECTED_INTERIOR_ATOMS
        ),
        "geometry_has_no_coordination_anomalies": (
            len(
                anomalous_indices
            )
            == 0
        ),
        "geometry_BN_bond_lengths_are_plausible": (
            len(
                bond_distances
            )
            > 0
            and float(
                np.min(
                    bond_distances
                )
            )
            >= BOND_VALID_MIN_NM
            and float(
                np.max(
                    bond_distances
                )
            )
            <= BOND_VALID_MAX_NM
        ),
        "lower_end_has_60_atoms": (
            int(
                lower_metrics[
                    "atom_count"
                ]
            )
            == EXPECTED_EDGE_ATOMS_PER_END
        ),
        "upper_end_has_60_atoms": (
            int(
                upper_metrics[
                    "atom_count"
                ]
            )
            == EXPECTED_EDGE_ATOMS_PER_END
        ),
        "lower_end_has_30B_and_30N": (
            int(
                lower_metrics[
                    "B_count"
                ]
            )
            == EXPECTED_EDGE_B_PER_END
            and int(
                lower_metrics[
                    "N_count"
                ]
            )
            == EXPECTED_EDGE_N_PER_END
        ),
        "upper_end_has_30B_and_30N": (
            int(
                upper_metrics[
                    "B_count"
                ]
            )
            == EXPECTED_EDGE_B_PER_END
            and int(
                upper_metrics[
                    "N_count"
                ]
            )
            == EXPECTED_EDGE_N_PER_END
        ),
        "lower_end_is_planar": (
            float(
                lower_metrics[
                    "axial_standard_deviation_nm"
                ]
            )
            <= MAX_END_AXIAL_STD_NM
        ),
        "upper_end_is_planar": (
            float(
                upper_metrics[
                    "axial_standard_deviation_nm"
                ]
            )
            <= MAX_END_AXIAL_STD_NM
        ),
        "terminal_radii_are_symmetric": (
            abs(
                lower_end_radius
                - upper_end_radius
            )
            <= MAX_END_RADIUS_DIFFERENCE_NM
        ),
        "lower_end_element_sequence_is_alternating": (
            float(
                lower_metrics[
                    "element_alternation_fraction"
                ]
            )
            >= MIN_EDGE_ALTERNATION_FRACTION
        ),
        "upper_end_element_sequence_is_alternating": (
            float(
                upper_metrics[
                    "element_alternation_fraction"
                ]
            )
            >= MIN_EDGE_ALTERNATION_FRACTION
        ),
        "validated_aperture_is_finite_and_open": (
            math.isfinite(
                target_aperture_diameter_nm
            )
            and target_aperture_diameter_nm
            > 0.0
            and target_aperture_radius_nm
            < mean_parent_rim_radius_nm
        ),
        "validated_open_area_fraction_is_valid": (
            math.isfinite(
                target_open_area_fraction
            )
            and 0.0
            < target_open_area_fraction
            < 1.0
        ),
        "required_radial_occlusion_is_positive": (
            required_radial_occlusion_nm
            > 0.0
        ),
        "monolayer_annulus_population_is_comparable_to_R2_beads": (
            (
                0.70
                * EXPECTED_CAPS_PER_END
            )
            <= estimated_monolayer_annulus_atoms
            <= (
                1.30
                * EXPECTED_CAPS_PER_END
            )
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
        PRIMARY_DECISION
        if accepted
        else
        "R2_PARENT_RIM_CHEMICAL_AUDIT_REQUIRES_REVIEW"
    )

    required_next_step = (
        "DEFINE_AND_RANK_R2_EXPLICIT_END_RIM_CHEMISTRY_CANDIDATES"
        if accepted
        else
        "REVIEW_R2_PARENT_RIM_CHEMICAL_AUDIT_FAILURES"
    )

    source_manifest = [
        {
            "role": (
                "R2_architecture_selection"
            ),
            "file": relative(
                SELECTION_SUMMARY
            ),
            "sha256": sha256(
                SELECTION_SUMMARY
            ),
        },
        {
            "role": (
                "R2_geometry_summary"
            ),
            "file": relative(
                GEOMETRY_SUMMARY
            ),
            "sha256": sha256(
                GEOMETRY_SUMMARY
            ),
        },
        {
            "role": (
                "R2_static_summary"
            ),
            "file": relative(
                STATIC_SUMMARY
            ),
            "sha256": sha256(
                STATIC_SUMMARY
            ),
        },
        {
            "role": (
                "R2_parent_coordinates"
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
                "R2_parent_TPR"
            ),
            "file": relative(
                SYSTEM_TPR
            ),
            "sha256": sha256(
                SYSTEM_TPR
            ),
        },
    ]

    write_csv(
        SOURCE_MANIFEST_CSV,
        source_manifest,
    )

    summary = {
        "decision": decision,
        "system_atoms": len(atoms),
        "HBN_atoms": (
            len(
                hbn_global_indices
            )
        ),
        "HBN_B_atoms": b_count,
        "HBN_N_atoms": n_count,
        "HBN_zero_charge_atoms": (
            zero_charge_count
        ),
        "geometry_derived_BN_bonds": (
            len(
                bond_rows
            )
        ),
        "expected_geometry_BN_bonds": (
            EXPECTED_GEOMETRY_BONDS
        ),
        "BN_bond_mean_nm": (
            mean_bn_bond_nm
        ),
        "BN_bond_median_nm": (
            median_bn_bond_nm
        ),
        "BN_bond_minimum_nm": float(
            np.min(
                bond_distances
            )
        ),
        "BN_bond_maximum_nm": float(
            np.max(
                bond_distances
            )
        ),
        "degree2_terminal_atoms": (
            len(
                edge_indices
            )
        ),
        "degree3_interior_atoms": (
            len(
                interior_indices
            )
        ),
        "coordination_anomaly_atoms": (
            len(
                anomalous_indices
            )
        ),
        "edge_edge_bonds": (
            edge_edge_bonds
        ),
        "edge_interior_bonds": (
            edge_interior_bonds
        ),
        "lower_end_atoms": (
            lower_metrics[
                "atom_count"
            ]
        ),
        "lower_end_B_atoms": (
            lower_metrics[
                "B_count"
            ]
        ),
        "lower_end_N_atoms": (
            lower_metrics[
                "N_count"
            ]
        ),
        "lower_end_radius_mean_nm": (
            lower_metrics[
                "radius_mean_nm"
            ]
        ),
        "lower_end_axial_std_nm": (
            lower_metrics[
                "axial_standard_deviation_nm"
            ]
        ),
        "lower_end_alternation_fraction": (
            lower_metrics[
                "element_alternation_fraction"
            ]
        ),
        "upper_end_atoms": (
            upper_metrics[
                "atom_count"
            ]
        ),
        "upper_end_B_atoms": (
            upper_metrics[
                "B_count"
            ]
        ),
        "upper_end_N_atoms": (
            upper_metrics[
                "N_count"
            ]
        ),
        "upper_end_radius_mean_nm": (
            upper_metrics[
                "radius_mean_nm"
            ]
        ),
        "upper_end_axial_std_nm": (
            upper_metrics[
                "axial_standard_deviation_nm"
            ]
        ),
        "upper_end_alternation_fraction": (
            upper_metrics[
                "element_alternation_fraction"
            ]
        ),
        "lower_cap_plane_offset_nm": (
            lower_cap_offset_nm
        ),
        "upper_cap_plane_offset_nm": (
            upper_cap_offset_nm
        ),
        "target_aperture_diameter_nm": (
            target_aperture_diameter_nm
        ),
        "target_aperture_radius_nm": (
            target_aperture_radius_nm
        ),
        "target_open_area_fraction": (
            target_open_area_fraction
        ),
        "parent_rim_mean_radius_nm": (
            mean_parent_rim_radius_nm
        ),
        "required_radial_occlusion_nm": (
            required_radial_occlusion_nm
        ),
        "annular_area_per_end_nm2": (
            annulus_area_nm2
        ),
        "estimated_monolayer_hBN_atoms_per_end": (
            estimated_monolayer_annulus_atoms
        ),
        "validated_R2_beads_per_end": (
            EXPECTED_CAPS_PER_END
        ),
        "validated_minimum_CAP_OW_distance_nm": (
            validated_minimum_cap_ow_nm
        ),
        "current_parent_rim_is_chemically_passivated": False,
        "current_parent_rim_is_final_chemistry": False,
        "static_candidate_design_authorized": (
            accepted
        ),
        "geometry_generation_authorized": False,
        "energy_minimization_authorized": False,
        "new_MD_authorized": False,
        "mobile_MD_authorized": False,
        "multitemperature_MD_authorized": False,
        "QM_recalculation_authorized": False,
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

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed
        in gates.items()
    )

    candidate_lines = "\n".join(
        (
            f"- `{row['candidate_id']}` — "
            f"priority {row['priority']}: "
            f"**{row['status']}**"
        )
        for row in candidate_rows
    )

    REPORT_MD.write_text(
        f"""# R2 Parent Rim and Chemical-Constraint Audit

## Scope

This audit characterizes the accepted R2 parent BNNT and translates
the validated steric-cap geometry into chemical and atomistic design
constraints.

No coordinates, topology, trajectory, checkpoint, or accepted
simulation result were modified.

No minimization, molecular dynamics, DFT, TDDFT, or other quantum
calculation was executed.

## Parent system

- Total atoms:
  **{len(atoms)}**
- HBN atoms:
  **{len(hbn_global_indices)}**
- B/N atoms:
  **{b_count}/{n_count}**
- Pyrene atoms:
  **{len(pyr_global_indices)}**
- Water atoms:
  **{len(water_global_indices)}**
- Cap beads:
  **{len(cap_global_indices)}**

The accepted HBN screening topology assigns zero charge to all
{EXPECTED_HBN_ATOMS} B/N atoms. This is retained as a record of the
screening model and is not interpreted as final chemical
electrostatics.

## Geometry-derived BN connectivity

- Bonds:
  **{len(bond_rows)}**
- Expected bonds:
  **{EXPECTED_GEOMETRY_BONDS}**
- Mean/median BN distance:
  **{mean_bn_bond_nm:.6f}/{median_bn_bond_nm:.6f} nm**
- Minimum/maximum BN distance:
  **{np.min(bond_distances):.6f}/{np.max(bond_distances):.6f} nm**
- Degree-2 terminal atoms:
  **{len(edge_indices)}**
- Degree-3 interior atoms:
  **{len(interior_indices)}**
- Coordination anomalies:
  **{len(anomalous_indices)}**

The degree-2 atoms represent unsaturated terminal-rim sites in the
current unpassivated parent model. Frozen-coordinate stability does
not itself establish chemical stability of these sites.

## Terminal rims

### Lower end

- Atoms:
  **{lower_metrics['atom_count']}**
- B/N:
  **{lower_metrics['B_count']}/{lower_metrics['N_count']}**
- Mean radius:
  **{lower_metrics['radius_mean_nm']:.6f} nm**
- Axial standard deviation:
  **{lower_metrics['axial_standard_deviation_nm']:.6f} nm**
- Element-alternation fraction:
  **{lower_metrics['element_alternation_fraction']:.6f}**
- Cap-plane offset:
  **{lower_cap_offset_nm:.6f} nm**

### Upper end

- Atoms:
  **{upper_metrics['atom_count']}**
- B/N:
  **{upper_metrics['B_count']}/{upper_metrics['N_count']}**
- Mean radius:
  **{upper_metrics['radius_mean_nm']:.6f} nm**
- Axial standard deviation:
  **{upper_metrics['axial_standard_deviation_nm']:.6f} nm**
- Element-alternation fraction:
  **{upper_metrics['element_alternation_fraction']:.6f}**
- Cap-plane offset:
  **{upper_cap_offset_nm:.6f} nm**

## Validated R2 target

- Effective aperture diameter:
  **{target_aperture_diameter_nm:.6f} nm**
- Effective aperture radius:
  **{target_aperture_radius_nm:.6f} nm**
- Open-area fraction:
  **{target_open_area_fraction:.6f}**
- Mean parent-rim radius:
  **{mean_parent_rim_radius_nm:.6f} nm**
- Required radial occlusion:
  **{required_radial_occlusion_nm:.6f} nm**
- Annular area per end:
  **{annulus_area_nm2:.6f} nm²**
- Estimated one-layer h-BN annular population:
  **{estimated_monolayer_annulus_atoms:.3f} atoms/end**
- Validated steric beads:
  **{EXPECTED_CAPS_PER_END} beads/end**
- Validated minimum CAP–OW separation:
  **{validated_minimum_cap_ow_nm:.6f} nm**

Simple H/OH/NHx passivation of the existing terminal atoms cannot by
itself reproduce the validated aperture because the required inward
radial closure is approximately
**{required_radial_occlusion_nm:.6f} nm**.

## Preliminary candidate classes

{candidate_lines}

These are candidate classes, not selected final chemistries.

The leading geometric analogue is a separately defined annular
nanoflake or rigid macrocycle positioned at the validated cap plane
and attached through an explicit chemically valid junction. A direct
unstrained 90-degree seamless sp2 junction is not assumed.

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Static candidate-design work authorized:
  **{'YES' if accepted else 'NO'}**
- Explicit geometry generation authorized:
  **NO**
- Energy minimization authorized:
  **NO**
- New MD authorized:
  **NO**
- Mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `{required_next_step}`
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 parent-rim and chemical-constraint "
        "audit completed."
    )

    print(
        "System / HBN / PYR / water / cap atoms: "
        f"{len(atoms)}/"
        f"{len(hbn_global_indices)}/"
        f"{len(pyr_global_indices)}/"
        f"{len(water_global_indices)}/"
        f"{len(cap_global_indices)}"
    )

    print(
        "HBN B / N / zero-charge atoms: "
        f"{b_count}/"
        f"{n_count}/"
        f"{zero_charge_count}"
    )

    print(
        "Geometry-derived BN bonds / expected: "
        f"{len(bond_rows)}/"
        f"{EXPECTED_GEOMETRY_BONDS}"
    )

    print(
        "BN bond mean/median/min/max: "
        f"{mean_bn_bond_nm:.6f}/"
        f"{median_bn_bond_nm:.6f}/"
        f"{np.min(bond_distances):.6f}/"
        f"{np.max(bond_distances):.6f} nm"
    )

    print(
        "Degree-2 edge / degree-3 interior / "
        "anomalous atoms: "
        f"{len(edge_indices)}/"
        f"{len(interior_indices)}/"
        f"{len(anomalous_indices)}"
    )

    print(
        "Lower rim atoms / B / N / radius / "
        "axial std / alternation: "
        f"{lower_metrics['atom_count']}/"
        f"{lower_metrics['B_count']}/"
        f"{lower_metrics['N_count']}/"
        f"{lower_metrics['radius_mean_nm']:.6f}/"
        f"{lower_metrics['axial_standard_deviation_nm']:.6f}/"
        f"{lower_metrics['element_alternation_fraction']:.6f}"
    )

    print(
        "Upper rim atoms / B / N / radius / "
        "axial std / alternation: "
        f"{upper_metrics['atom_count']}/"
        f"{upper_metrics['B_count']}/"
        f"{upper_metrics['N_count']}/"
        f"{upper_metrics['radius_mean_nm']:.6f}/"
        f"{upper_metrics['axial_standard_deviation_nm']:.6f}/"
        f"{upper_metrics['element_alternation_fraction']:.6f}"
    )

    print(
        "Lower / upper cap-plane offsets: "
        f"{lower_cap_offset_nm:.6f}/"
        f"{upper_cap_offset_nm:.6f} nm"
    )

    print(
        "Target aperture diameter / radius / "
        "open-area fraction: "
        f"{target_aperture_diameter_nm:.6f}/"
        f"{target_aperture_radius_nm:.6f}/"
        f"{target_open_area_fraction:.6f}"
    )

    print(
        "Parent rim radius / required radial occlusion: "
        f"{mean_parent_rim_radius_nm:.6f}/"
        f"{required_radial_occlusion_nm:.6f} nm"
    )

    print(
        "Annulus area / estimated hBN atoms / "
        "validated beads per end: "
        f"{annulus_area_nm2:.6f} nm^2 / "
        f"{estimated_monolayer_annulus_atoms:.3f} / "
        f"{EXPECTED_CAPS_PER_END}"
    )

    print(
        "Current parent rim chemically passivated: NO"
    )

    print(
        "Current parent rim is final chemistry: NO"
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
        "Static candidate-design work authorized: "
        f"{'YES' if accepted else 'NO'}"
    )

    print(
        "Explicit geometry generation authorized: NO"
    )

    print(
        "Energy minimization authorized: NO"
    )

    print(
        "New MD authorized: NO"
    )

    print(
        "Mobile MD authorized: NO"
    )

    print(
        "Multitemperature MD authorized: NO"
    )

    print(
        "QM recalculation authorized: NO"
    )

    print(
        f"Required next step: {required_next_step}"
    )

    for path in (
        TPR_DUMP,
        HBN_ATOMS_CSV,
        HBN_BONDS_CSV,
        RIM_ATOMS_CSV,
        END_SUMMARY_CSV,
        CONSTRAINTS_CSV,
        CANDIDATES_CSV,
        SOURCE_MANIFEST_CSV,
        SUMMARY_CSV,
        GATES_CSV,
        REPORT_MD,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not accepted:
        raise RuntimeError(
            "R2 parent-rim chemical audit "
            "requires review."
        )


if __name__ == "__main__":
    main()
