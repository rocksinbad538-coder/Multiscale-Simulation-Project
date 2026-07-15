#!/usr/bin/env python3

from __future__ import annotations

import csv
import os
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

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

TPR_PATH = (
    ACCEPTED_ROOT
    / "nvt_100ps_frozenSolute.tpr"
)

XTC_PATH = (
    ACCEPTED_ROOT
    / "nvt_100ps_frozenSolute.xtc"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/"
    "day020_confined_water_axial_radial_density"
)

INDEX_PATH = (
    OUTPUT_ROOT
    / "confined_water_analysis_groups.ndx"
)

GEOMETRY_CSV = (
    OUTPUT_ROOT
    / "nanotube_geometry_summary.csv"
)

PYRENE_CSV = (
    OUTPUT_ROOT
    / "pyrene_geometry_summary.csv"
)

WATER_ROLES_CSV = (
    OUTPUT_ROOT
    / "water_atom_roles.csv"
)

PARAMETERS_CSV = (
    OUTPUT_ROOT
    / "densmap_parameters.csv"
)

DENSITY_DAT = (
    OUTPUT_ROOT
    / "water_oxygen_axial_radial_density.dat"
)

DENSITY_XPM = (
    OUTPUT_ROOT
    / "water_oxygen_axial_radial_density.xpm"
)

DENSITY_LOG = (
    OUTPUT_ROOT
    / "gromacs_densmap.log"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "CONFINED_WATER_AXIAL_RADIAL_DENSITY_DAY020.md"
)

EXPECTED_HBN_ATOMS = 1680
EXPECTED_PYR_ATOMS = 104
EXPECTED_PYR_RESIDUES = 4
EXPECTED_WATER_RESIDUES = 16634
EXPECTED_TOTAL_ATOMS = 68320

END_GROUP_FRACTION = 0.10
BIN_WIDTH_NM = 0.05
AXIAL_PADDING_NM = 1.00
RADIAL_PADDING_NM = 1.50

BEGIN_TIME_PS = 0.0
END_TIME_PS = 100.0


def log(message: str = "") -> None:
    print(message, flush=True)


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
) -> None:
    if not rows:
        raise RuntimeError(
            f"No rows available for {path}"
        )

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


def relative(path: Path) -> str:
    return str(path.relative_to(PROJECT_ROOT))


def find_gromacs() -> str:
    configured = os.environ.get("GMX_BIN")

    candidates = []

    if configured:
        candidates.append(configured)

    candidates.extend(
        [
            "/usr/local/gromacs/bin/gmx",
            "gmx",
            "gmx_mpi",
        ]
    )

    for candidate in candidates:
        resolved = shutil.which(candidate)

        if resolved is not None:
            return resolved

    raise RuntimeError(
        "Could not locate a GROMACS executable. "
        "Set GMX_BIN explicitly."
    )


def parse_box(
    values: list[float],
) -> np.ndarray:
    if len(values) == 3:
        return np.diag(values)

    if len(values) == 9:
        # GRO triclinic ordering:
        # v1x v2y v3z v1y v1z v2x v2z v3x v3y
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
            dtype=np.float64,
        )

    raise RuntimeError(
        f"Unsupported GRO box with "
        f"{len(values)} entries"
    )


def parse_gro(
    path: Path,
) -> tuple[
    list[dict[str, object]],
    np.ndarray,
]:
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

    if atom_count != EXPECTED_TOTAL_ATOMS:
        raise RuntimeError(
            f"Expected {EXPECTED_TOTAL_ATOMS} atoms, "
            f"found {atom_count}"
        )

    atoms: list[
        dict[str, object]
    ] = []

    for line_index, line in enumerate(
        lines[2 : 2 + atom_count],
        start=1,
    ):
        if len(line) < 44:
            raise RuntimeError(
                f"Malformed GRO atom line "
                f"{line_index}"
            )

        residue_number = int(
            line[0:5]
        )

        residue_name = (
            line[5:10].strip()
        )

        atom_name = (
            line[10:15].strip()
        )

        atom_number = int(
            line[15:20]
        )

        coordinate = np.asarray(
            [
                float(line[20:28]),
                float(line[28:36]),
                float(line[36:44]),
            ],
            dtype=np.float64,
        )

        atoms.append(
            {
                "atom_index": line_index,
                "gro_atom_number": atom_number,
                "residue_number": residue_number,
                "residue_name": residue_name,
                "atom_name": atom_name,
                "coordinate_nm": coordinate,
            }
        )

    box_values = [
        float(item)
        for item in lines[
            2 + atom_count
        ].split()
    ]

    box = parse_box(
        box_values
    )

    return atoms, box


def select_atoms(
    atoms: list[dict[str, object]],
    residue_name: str,
) -> list[dict[str, object]]:
    return [
        atom
        for atom in atoms
        if atom["residue_name"]
        == residue_name
    ]


def orient_axis(
    axis: np.ndarray,
) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)

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


def calculate_nanotube_geometry(
    hbn_atoms: list[dict[str, object]],
) -> dict[str, object]:
    coordinates = np.stack(
        [
            atom["coordinate_nm"]
            for atom in hbn_atoms
        ]
    )

    center = np.mean(
        coordinates,
        axis=0,
    )

    centered = coordinates - center

    covariance = (
        centered.T
        @ centered
        / coordinates.shape[0]
    )

    eigenvalues, eigenvectors = (
        np.linalg.eigh(covariance)
    )

    order = np.argsort(
        eigenvalues
    )[::-1]

    eigenvalues = eigenvalues[
        order
    ]

    eigenvectors = eigenvectors[
        :,
        order,
    ]

    axis = orient_axis(
        eigenvectors[:, 0]
    )

    axial_projection = (
        centered @ axis
    )

    radial_vectors = (
        centered
        - np.outer(
            axial_projection,
            axis,
        )
    )

    radial_distance = np.linalg.norm(
        radial_vectors,
        axis=1,
    )

    axial_min = float(
        np.min(
            axial_projection
        )
    )

    axial_max = float(
        np.max(
            axial_projection
        )
    )

    axial_span = (
        axial_max
        - axial_min
    )

    lower_cutoff = np.quantile(
        axial_projection,
        END_GROUP_FRACTION,
    )

    upper_cutoff = np.quantile(
        axial_projection,
        1.0
        - END_GROUP_FRACTION,
    )

    minus_indices = [
        int(
            hbn_atoms[index][
                "atom_index"
            ]
        )
        for index in np.where(
            axial_projection
            <= lower_cutoff
        )[0]
    ]

    plus_indices = [
        int(
            hbn_atoms[index][
                "atom_index"
            ]
        )
        for index in np.where(
            axial_projection
            >= upper_cutoff
        )[0]
    ]

    if (
        len(minus_indices) < 20
        or len(plus_indices) < 20
    ):
        raise RuntimeError(
            "Axis endpoint groups are too small"
        )

    nearest_cartesian_axis = (
        "xyz"[
            int(
                np.argmax(
                    np.abs(axis)
                )
            )
        ]
    )

    return {
        "center_nm": center,
        "axis_unit_vector": axis,
        "pca_eigenvalues_nm2": (
            eigenvalues
        ),
        "axial_projection_nm": (
            axial_projection
        ),
        "radial_distance_nm": (
            radial_distance
        ),
        "axial_min_nm": axial_min,
        "axial_max_nm": axial_max,
        "axial_span_nm": axial_span,
        "mean_wall_radius_nm": float(
            np.mean(
                radial_distance
            )
        ),
        "median_wall_radius_nm": float(
            np.median(
                radial_distance
            )
        ),
        "p05_wall_radius_nm": float(
            np.percentile(
                radial_distance,
                5.0,
            )
        ),
        "p95_wall_radius_nm": float(
            np.percentile(
                radial_distance,
                95.0,
            )
        ),
        "axis_minus_indices": (
            minus_indices
        ),
        "axis_plus_indices": (
            plus_indices
        ),
        "nearest_cartesian_axis": (
            nearest_cartesian_axis
        ),
    }


def identify_water_roles(
    water_atoms: list[dict[str, object]],
) -> tuple[
    list[int],
    list[int],
    list[int],
    list[dict[str, object]],
]:
    by_residue: dict[
        int,
        list[dict[str, object]],
    ] = defaultdict(list)

    for atom in water_atoms:
        by_residue[
            int(
                atom["residue_number"]
            )
        ].append(atom)

    if len(by_residue) != (
        EXPECTED_WATER_RESIDUES
    ):
        raise RuntimeError(
            f"Expected {EXPECTED_WATER_RESIDUES} "
            f"water residues, found "
            f"{len(by_residue)}"
        )

    oxygen_indices: list[int] = []
    hydrogen_indices: list[int] = []
    virtual_site_indices: list[int] = []

    name_counts: dict[
        str,
        int,
    ] = defaultdict(int)

    role_counts: dict[
        tuple[str, str],
        int,
    ] = defaultdict(int)

    for residue_number, residue_atoms in sorted(
        by_residue.items()
    ):
        oxygen_candidates = [
            atom
            for atom in residue_atoms
            if str(
                atom["atom_name"]
            ).upper().startswith("O")
        ]

        hydrogen_candidates = [
            atom
            for atom in residue_atoms
            if str(
                atom["atom_name"]
            ).upper().startswith("H")
        ]

        remaining = [
            atom
            for atom in residue_atoms
            if atom
            not in oxygen_candidates
            and atom
            not in hydrogen_candidates
        ]

        if len(oxygen_candidates) != 1:
            raise RuntimeError(
                f"Water residue {residue_number} "
                f"has {len(oxygen_candidates)} "
                f"oxygen candidates"
            )

        if len(hydrogen_candidates) != 2:
            raise RuntimeError(
                f"Water residue {residue_number} "
                f"has {len(hydrogen_candidates)} "
                f"hydrogen candidates"
            )

        if len(remaining) != 1:
            raise RuntimeError(
                f"Water residue {residue_number} "
                f"has {len(remaining)} virtual-site "
                f"candidates"
            )

        oxygen_indices.append(
            int(
                oxygen_candidates[0][
                    "atom_index"
                ]
            )
        )

        hydrogen_indices.extend(
            int(
                atom["atom_index"]
            )
            for atom in hydrogen_candidates
        )

        virtual_site_indices.append(
            int(
                remaining[0][
                    "atom_index"
                ]
            )
        )

        for atom in residue_atoms:
            atom_name = str(
                atom["atom_name"]
            )

            name_counts[
                atom_name
            ] += 1

        for atom in oxygen_candidates:
            role_counts[
                (
                    str(
                        atom["atom_name"]
                    ),
                    "oxygen",
                )
            ] += 1

        for atom in hydrogen_candidates:
            role_counts[
                (
                    str(
                        atom["atom_name"]
                    ),
                    "hydrogen",
                )
            ] += 1

        for atom in remaining:
            role_counts[
                (
                    str(
                        atom["atom_name"]
                    ),
                    "virtual_site",
                )
            ] += 1

    role_rows = [
        {
            "atom_name": atom_name,
            "assigned_role": role,
            "atom_count": count,
        }
        for (
            atom_name,
            role,
        ), count in sorted(
            role_counts.items()
        )
    ]

    return (
        oxygen_indices,
        hydrogen_indices,
        virtual_site_indices,
        role_rows,
    )


def split_pyrene_residues(
    pyrene_atoms: list[dict[str, object]],
    geometry: dict[str, object],
) -> tuple[
    dict[str, list[int]],
    list[dict[str, object]],
]:
    by_residue: dict[
        int,
        list[dict[str, object]],
    ] = defaultdict(list)

    for atom in pyrene_atoms:
        by_residue[
            int(
                atom["residue_number"]
            )
        ].append(atom)

    if len(by_residue) != (
        EXPECTED_PYR_RESIDUES
    ):
        raise RuntimeError(
            f"Expected {EXPECTED_PYR_RESIDUES} "
            f"PYR residues, found "
            f"{len(by_residue)}"
        )

    center = np.asarray(
        geometry["center_nm"],
        dtype=np.float64,
    )

    axis = np.asarray(
        geometry[
            "axis_unit_vector"
        ],
        dtype=np.float64,
    )

    groups: dict[
        str,
        list[int],
    ] = {}

    rows: list[
        dict[str, object]
    ] = []

    ordered_residues = sorted(
        by_residue
    )

    for pyrene_index, residue_number in enumerate(
        ordered_residues,
        start=1,
    ):
        atoms = by_residue[
            residue_number
        ]

        coordinates = np.stack(
            [
                atom["coordinate_nm"]
                for atom in atoms
            ]
        )

        pyrene_center = np.mean(
            coordinates,
            axis=0,
        )

        relative_vector = (
            pyrene_center
            - center
        )

        axial_position = float(
            relative_vector @ axis
        )

        radial_vector = (
            relative_vector
            - axial_position
            * axis
        )

        radial_position = float(
            np.linalg.norm(
                radial_vector
            )
        )

        group_name = (
            f"PYR_{pyrene_index}"
        )

        groups[group_name] = [
            int(
                atom["atom_index"]
            )
            for atom in atoms
        ]

        rows.append(
            {
                "pyrene_label": group_name,
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
            }
        )

    return groups, rows


def write_ndx_group(
    handle,
    group_name: str,
    indices: list[int],
) -> None:
    handle.write(
        f"[ {group_name} ]\n"
    )

    for start in range(
        0,
        len(indices),
        15,
    ):
        chunk = indices[
            start : start + 15
        ]

        handle.write(
            " ".join(
                f"{index:6d}"
                for index in chunk
            )
        )

        handle.write("\n")

    handle.write("\n")


def write_index(
    geometry: dict[str, object],
    hbn_indices: list[int],
    pyrene_indices: list[int],
    pyrene_groups: dict[
        str,
        list[int],
    ],
    oxygen_indices: list[int],
    hydrogen_indices: list[int],
    virtual_site_indices: list[int],
) -> None:
    with INDEX_PATH.open(
        "w",
        encoding="utf-8",
    ) as handle:
        # These first three groups are used by densmap.
        write_ndx_group(
            handle,
            "AxisMinus",
            list(
                geometry[
                    "axis_minus_indices"
                ]
            ),
        )

        write_ndx_group(
            handle,
            "AxisPlus",
            list(
                geometry[
                    "axis_plus_indices"
                ]
            ),
        )

        write_ndx_group(
            handle,
            "Water_O",
            oxygen_indices,
        )

        write_ndx_group(
            handle,
            "HBN",
            hbn_indices,
        )

        write_ndx_group(
            handle,
            "PYR_all",
            pyrene_indices,
        )

        for group_name in sorted(
            pyrene_groups
        ):
            write_ndx_group(
                handle,
                group_name,
                pyrene_groups[
                    group_name
                ],
            )

        write_ndx_group(
            handle,
            "Water_H",
            hydrogen_indices,
        )

        write_ndx_group(
            handle,
            "Water_M",
            virtual_site_indices,
        )


def run_densmap(
    gromacs: str,
    amax_nm: float,
    rmax_nm: float,
) -> None:
    command = [
        gromacs,
        "densmap",
        "-f",
        str(XTC_PATH),
        "-s",
        str(TPR_PATH),
        "-n",
        str(INDEX_PATH),
        "-od",
        str(DENSITY_DAT),
        "-o",
        str(DENSITY_XPM),
        "-b",
        f"{BEGIN_TIME_PS:.6f}",
        "-e",
        f"{END_TIME_PS:.6f}",
        "-bin",
        f"{BIN_WIDTH_NM:.6f}",
        "-amax",
        f"{amax_nm:.6f}",
        "-rmax",
        f"{rmax_nm:.6f}",
        "-unit",
        "nm-3",
        "-sums",
    ]

    log("")
    log("Running:")
    log(" ".join(command))
    log("")

    with DENSITY_LOG.open(
        "w",
        encoding="utf-8",
    ) as log_handle:
        process = subprocess.Popen(
            command,
            cwd=PROJECT_ROOT,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        assert process.stdin is not None
        assert process.stdout is not None

        process.stdin.write(
            "AxisMinus\n"
            "AxisPlus\n"
            "Water_O\n"
        )

        process.stdin.close()

        for line in process.stdout:
            print(
                line,
                end="",
                flush=True,
            )

            log_handle.write(line)
            log_handle.flush()

        return_code = process.wait()

    if return_code != 0:
        raise RuntimeError(
            f"gmx densmap failed with "
            f"return code {return_code}"
        )

    if (
        not DENSITY_DAT.exists()
        or DENSITY_DAT.stat().st_size
        == 0
    ):
        raise RuntimeError(
            "densmap DAT output is missing "
            "or empty"
        )

    if (
        not DENSITY_XPM.exists()
        or DENSITY_XPM.stat().st_size
        == 0
    ):
        raise RuntimeError(
            "densmap XPM output is missing "
            "or empty"
        )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required_path in (
        GRO_PATH,
        TPR_PATH,
        XTC_PATH,
    ):
        if not required_path.exists():
            raise RuntimeError(
                f"Missing required file: "
                f"{required_path}"
            )

    atoms, box = parse_gro(
        GRO_PATH
    )

    hbn_atoms = select_atoms(
        atoms,
        "HBN",
    )

    pyrene_atoms = select_atoms(
        atoms,
        "PYR",
    )

    water_atoms = select_atoms(
        atoms,
        "SOL",
    )

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

    geometry = (
        calculate_nanotube_geometry(
            hbn_atoms
        )
    )

    (
        oxygen_indices,
        hydrogen_indices,
        virtual_site_indices,
        water_role_rows,
    ) = identify_water_roles(
        water_atoms
    )

    (
        pyrene_groups,
        pyrene_rows,
    ) = split_pyrene_residues(
        pyrene_atoms,
        geometry,
    )

    hbn_indices = [
        int(
            atom["atom_index"]
        )
        for atom in hbn_atoms
    ]

    pyrene_indices = [
        int(
            atom["atom_index"]
        )
        for atom in pyrene_atoms
    ]

    write_index(
        geometry,
        hbn_indices,
        pyrene_indices,
        pyrene_groups,
        oxygen_indices,
        hydrogen_indices,
        virtual_site_indices,
    )

    box_vector_lengths = (
        np.linalg.norm(
            box,
            axis=1,
        )
    )

    half_minimum_box_length = (
        0.5
        * float(
            np.min(
                box_vector_lengths
            )
        )
    )

    amax_nm = (
        0.5
        * float(
            geometry[
                "axial_span_nm"
            ]
        )
        + AXIAL_PADDING_NM
    )

    proposed_rmax_nm = (
        float(
            geometry[
                "mean_wall_radius_nm"
            ]
        )
        + RADIAL_PADDING_NM
    )

    rmax_nm = min(
        proposed_rmax_nm,
        half_minimum_box_length
        - 0.10,
    )

    if rmax_nm <= float(
        geometry[
            "p95_wall_radius_nm"
        ]
    ):
        raise RuntimeError(
            "The available radial range does "
            "not extend beyond the HBN wall"
        )

    axis = np.asarray(
        geometry[
            "axis_unit_vector"
        ]
    )

    eigenvalues = np.asarray(
        geometry[
            "pca_eigenvalues_nm2"
        ]
    )

    center = np.asarray(
        geometry[
            "center_nm"
        ]
    )

    geometry_rows = [
        {
            "source_gro": relative(
                GRO_PATH
            ),
            "HBN_atom_count": (
                len(hbn_atoms)
            ),
            "center_x_nm": center[0],
            "center_y_nm": center[1],
            "center_z_nm": center[2],
            "axis_x": axis[0],
            "axis_y": axis[1],
            "axis_z": axis[2],
            "nearest_cartesian_axis": (
                geometry[
                    "nearest_cartesian_axis"
                ]
            ),
            "pca_eigenvalue_1_nm2": (
                eigenvalues[0]
            ),
            "pca_eigenvalue_2_nm2": (
                eigenvalues[1]
            ),
            "pca_eigenvalue_3_nm2": (
                eigenvalues[2]
            ),
            "axial_min_nm": (
                geometry[
                    "axial_min_nm"
                ]
            ),
            "axial_max_nm": (
                geometry[
                    "axial_max_nm"
                ]
            ),
            "axial_span_nm": (
                geometry[
                    "axial_span_nm"
                ]
            ),
            "mean_wall_radius_nm": (
                geometry[
                    "mean_wall_radius_nm"
                ]
            ),
            "median_wall_radius_nm": (
                geometry[
                    "median_wall_radius_nm"
                ]
            ),
            "p05_wall_radius_nm": (
                geometry[
                    "p05_wall_radius_nm"
                ]
            ),
            "p95_wall_radius_nm": (
                geometry[
                    "p95_wall_radius_nm"
                ]
            ),
            "axis_minus_atom_count": len(
                geometry[
                    "axis_minus_indices"
                ]
            ),
            "axis_plus_atom_count": len(
                geometry[
                    "axis_plus_indices"
                ]
            ),
            "box_vector_1_nm": (
                box_vector_lengths[0]
            ),
            "box_vector_2_nm": (
                box_vector_lengths[1]
            ),
            "box_vector_3_nm": (
                box_vector_lengths[2]
            ),
        }
    ]

    parameter_rows = [
        {
            "parameter": (
                "begin_time_ps"
            ),
            "value": BEGIN_TIME_PS,
        },
        {
            "parameter": (
                "end_time_ps"
            ),
            "value": END_TIME_PS,
        },
        {
            "parameter": (
                "bin_width_nm"
            ),
            "value": BIN_WIDTH_NM,
        },
        {
            "parameter": (
                "axial_half_range_nm"
            ),
            "value": amax_nm,
        },
        {
            "parameter": (
                "radial_maximum_nm"
            ),
            "value": rmax_nm,
        },
        {
            "parameter": (
                "axis_endpoint_fraction"
            ),
            "value": END_GROUP_FRACTION,
        },
        {
            "parameter": (
                "density_unit"
            ),
            "value": "nm^-3",
        },
        {
            "parameter": (
                "water_position"
            ),
            "value": (
                "one oxygen site "
                "per SOL residue"
            ),
        },
    ]

    write_csv(
        GEOMETRY_CSV,
        geometry_rows,
    )

    write_csv(
        PYRENE_CSV,
        pyrene_rows,
    )

    write_csv(
        WATER_ROLES_CSV,
        water_role_rows,
    )

    write_csv(
        PARAMETERS_CSV,
        parameter_rows,
    )

    gromacs = find_gromacs()

    log(
        "Day020 confined-water "
        "axial-radial density analysis"
    )

    log(
        f"GROMACS: {gromacs}"
    )

    log(
        f"HBN atoms: {len(hbn_atoms)}"
    )

    log(
        f"PYR atoms: {len(pyrene_atoms)}"
    )

    log(
        f"Water molecules: "
        f"{len(oxygen_indices)}"
    )

    log(
        "Nanotube center (nm): "
        f"{center[0]:.6f}, "
        f"{center[1]:.6f}, "
        f"{center[2]:.6f}"
    )

    log(
        "Nanotube axis: "
        f"{axis[0]:.8f}, "
        f"{axis[1]:.8f}, "
        f"{axis[2]:.8f}"
    )

    log(
        "Nearest Cartesian axis: "
        f"{geometry['nearest_cartesian_axis']}"
    )

    log(
        "Axial span: "
        f"{float(geometry['axial_span_nm']):.6f} nm"
    )

    log(
        "Mean HBN wall radius: "
        f"{float(geometry['mean_wall_radius_nm']):.6f} nm"
    )

    log(
        "Densmap axial half-range: "
        f"{amax_nm:.6f} nm"
    )

    log(
        "Densmap radial maximum: "
        f"{rmax_nm:.6f} nm"
    )

    run_densmap(
        gromacs,
        amax_nm,
        rmax_nm,
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day020 Confined-Water "
            "Axial–Radial Density\n\n"
        )

        handle.write(
            "## Scope\n\n"
        )

        handle.write(
            "This analysis characterizes the "
            "number-density distribution of one "
            "oxygen site per TIP4P/2005 water "
            "molecule around the fixed HBN "
            "nanotube axis.\n\n"
        )

        handle.write(
            "## Geometry\n\n"
        )

        handle.write(
            f"- HBN atoms: "
            f"{len(hbn_atoms)}.\n"
        )

        handle.write(
            f"- Nanotube axis: "
            f"({axis[0]:.8f}, "
            f"{axis[1]:.8f}, "
            f"{axis[2]:.8f}).\n"
        )

        handle.write(
            f"- Axial span: "
            f"{float(geometry['axial_span_nm']):.6f} nm.\n"
        )

        handle.write(
            f"- Mean geometric HBN wall radius: "
            f"{float(geometry['mean_wall_radius_nm']):.6f} nm.\n"
        )

        handle.write(
            f"- Water molecules represented: "
            f"{len(oxygen_indices)}.\n\n"
        )

        handle.write(
            "## Density map\n\n"
        )

        handle.write(
            f"- Time interval: "
            f"{BEGIN_TIME_PS:.1f}–"
            f"{END_TIME_PS:.1f} ps.\n"
        )

        handle.write(
            f"- Grid spacing: "
            f"{BIN_WIDTH_NM:.3f} nm.\n"
        )

        handle.write(
            f"- Axial half-range: "
            f"{amax_nm:.6f} nm.\n"
        )

        handle.write(
            f"- Radial maximum: "
            f"{rmax_nm:.6f} nm.\n"
        )

        handle.write(
            "- Density unit: nm^-3.\n"
        )

        handle.write(
            "- The first two index groups "
            "define the fixed nanotube axis; "
            "the third group contains one "
            "water oxygen per molecule.\n\n"
        )

        handle.write(
            "## Interpretation boundary\n\n"
        )

        handle.write(
            "The map characterizes solvent "
            "organization around a frozen solute. "
            "It does not measure coupled "
            "water–solute conformational dynamics "
            "or scaffold thermal stability.\n"
        )

    log("")
    log(
        "Day020 confined-water "
        "axial-radial density completed."
    )

    log(
        f"Index: "
        f"{relative(INDEX_PATH)}"
    )

    log(
        f"Density DAT: "
        f"{relative(DENSITY_DAT)}"
    )

    log(
        f"Density XPM: "
        f"{relative(DENSITY_XPM)}"
    )

    log(
        f"Wrote: "
        f"{relative(OUTPUT_ROOT)}"
    )


if __name__ == "__main__":
    main()
