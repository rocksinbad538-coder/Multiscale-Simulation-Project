#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
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

R1_STATIC_ROOT = (
    DAY023_ROOT
    / "05_r1_full_static_scan"
)

R1_SELECTED_ROOT = (
    R1_STATIC_ROOT
    / "selected"
)

R1_PROTOTYPE_ROOT = (
    DAY023_ROOT
    / "02_r1_steric_cap_prototype"
)

R2_GEOMETRY_ROOT = (
    DAY023_ROOT
    / "12_r2_partial_cap_geometry_design"
)

R2_SELECTED_GEOMETRY_ROOT = (
    R2_GEOMETRY_ROOT
    / "selected"
)

OUTPUT_ROOT = (
    DAY023_ROOT
    / "13_r2_topology_static_scan"
)

SELECTED_ROOT = (
    OUTPUT_ROOT
    / "selected"
)

R1_SELECTED_TOPOLOGY = (
    R1_SELECTED_ROOT
    / "r1_selected_cap_model.top"
)

R1_PROTOTYPE_JSON = (
    R1_PROTOTYPE_ROOT
    / "r1_selected_steric_cap_definition.json"
)

R2_GEOMETRY_SUMMARY = (
    R2_GEOMETRY_ROOT
    / "r2_partial_cap_geometry_summary.csv"
)

R2_DEFINITION_JSON = (
    R2_SELECTED_GEOMETRY_ROOT
    / "r2_selected_partial_cap_definition.json"
)

R2_SYSTEM_GRO = (
    R2_SELECTED_GEOMETRY_ROOT
    / "r2_selected_partial_cap_geometry_only.gro"
)

R2_CAPS_GRO = (
    R2_SELECTED_GEOMETRY_ROOT
    / "r2_selected_partial_caps_only.gro"
)

R2_TOPOLOGY = (
    SELECTED_ROOT
    / "r2_selected_cap_model.top"
)

R2_STATIC_MDP = (
    SELECTED_ROOT
    / "r2_static_single_point.mdp"
)

R2_INDEX = (
    SELECTED_ROOT
    / "r2_static_energy_groups.ndx"
)

R2_TPR = (
    SELECTED_ROOT
    / "r2_static_single_point.tpr"
)

GROMPP_LOG = (
    OUTPUT_ROOT
    / "r2_static_grompp.log"
)

TPR_DUMP = (
    OUTPUT_ROOT
    / "r2_static_single_point_tpr_dump.txt"
)

TPR_DUMP_STDERR = (
    OUTPUT_ROOT
    / "r2_static_single_point_tpr_dump.stderr.log"
)

DEFFNM = (
    OUTPUT_ROOT
    / "r2_static_single_point"
)

MDRUN_CONSOLE = (
    OUTPUT_ROOT
    / "r2_static_single_point_mdrun_console.log"
)

ENERGY_MENU_LOG = (
    OUTPUT_ROOT
    / "r2_static_energy_menu.txt"
)

STATIC_ENERGY_CSV = (
    OUTPUT_ROOT
    / "r2_static_gromacs_energy.csv"
)

STRENGTH_SCAN_CSV = (
    OUTPUT_ROOT
    / "r2_cap_water_strength_scan.csv"
)

RADIAL_PROFILE_CSV = (
    OUTPUT_ROOT
    / "r2_selected_5kBT_aperture_radial_profile.csv"
)

AXIAL_PROFILE_CSV = (
    OUTPUT_ROOT
    / "r2_selected_5kBT_aperture_axial_profile.csv"
)

TOPOLOGY_AUDIT_CSV = (
    OUTPUT_ROOT
    / "r2_topology_audit.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_topology_static_scan_summary.csv"
)

GATE_CSV = (
    OUTPUT_ROOT
    / "r2_topology_static_scan_gates.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_TOPOLOGY_STATIC_CAP_WATER_SCAN_DAY023.md"
)

EXPECTED_ATOMS = 68332

HBN_ATOMS = 1680
PYR_ATOMS = 104
SOLUTE_ATOMS = HBN_ATOMS + PYR_ATOMS

WATERS = 16565
WATER_SITES = 4
WATER_ATOMS = WATERS * WATER_SITES

CAPS_PER_END = 144
TOTAL_CAPS = 288

CAPL_START = SOLUTE_ATOMS + WATER_ATOMS
CAPL_STOP = CAPL_START + CAPS_PER_END

CAPU_START = CAPL_STOP
CAPU_STOP = CAPU_START + CAPS_PER_END

TEMPERATURE_K = 300.0
BOLTZMANN_KJ_MOL_K = 0.00831446261815324
KBT_KJ_MOL = (
    BOLTZMANN_KJ_MOL_K
    * TEMPERATURE_K
)

SELECTED_STRENGTH_KBT = 5.0
SIGMA_CAP_OW_NM = -0.17

EXPECTED_SELECTED_EPSILON_KJ_MOL = (
    SELECTED_STRENGTH_KBT
    * KBT_KJ_MOL
    / 4.0
)

REFERENCE_R1_EPSILON_KJ_MOL = (
    3.1179234818
)

STRENGTHS_KBT = (
    5.0,
    10.0,
    20.0,
    40.0,
)

RADIAL_SCAN_MAX_NM = 0.80
RADIAL_SCAN_POINTS = 801

AXIAL_SCAN_HALF_WIDTH_NM = 0.60
AXIAL_SCAN_POINTS = 601

MIN_CAP_OW_DISTANCE_NM = 0.2195

MIN_ACTUAL_APERTURE_RADIUS_NM = 0.20
MAX_ACTUAL_APERTURE_RADIUS_NM = 0.45

MAX_LOWER_UPPER_APERTURE_DIFFERENCE_NM = 0.01

MAX_CENTERLINE_BARRIER_KBT = 5.0

MAX_GROMACS_ANALYTIC_ABSOLUTE_ERROR_KJ_MOL = 0.10
MAX_GROMACS_ANALYTIC_RELATIVE_ERROR = 5.0e-4

ZERO_ENERGY_TOLERANCE_KJ_MOL = 1.0e-4

MAX_SELECTED_CAP_SOL_LJ_KJ_MOL = 100.0
MAX_SELECTED_WATER_FORCE_KJ_MOL_NM = 250.0


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


def run_command(
    command: list[str],
    *,
    cwd: Path,
    input_text: str | None = None,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=cwd,
        env={
            **os.environ,
            "GMX_MAXBACKUP": "-1",
        },
        input=input_text,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def run_live(
    command: list[str],
    *,
    cwd: Path,
    log_path: Path,
) -> int:
    process = subprocess.Popen(
        command,
        cwd=cwd,
        env={
            **os.environ,
            "GMX_MAXBACKUP": "-1",
        },
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,
    )

    if process.stdout is None:
        raise RuntimeError(
            "Could not capture process output."
        )

    with log_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        for line in process.stdout:
            print(
                line,
                end="",
                flush=True,
            )

            handle.write(line)
            handle.flush()

    return process.wait()


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


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
                    key: row.get(
                        key,
                        "",
                    )
                    for key in fields
                }
            )


def parse_box(
    values: list[float],
) -> np.ndarray:
    if len(values) != 3:
        raise RuntimeError(
            "An orthorhombic box is required."
        )

    box = np.asarray(
        values,
        dtype=float,
    )

    if (
        not np.all(
            np.isfinite(box)
        )
        or np.any(box <= 0.0)
    ):
        raise RuntimeError(
            "Invalid box."
        )

    return box


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

    atom_count = int(
        lines[1].strip()
    )

    if len(lines) < atom_count + 3:
        raise RuntimeError(
            f"Incomplete GRO file: {path}"
        )

    atoms = []

    for index, line in enumerate(
        lines[
            2:
            2 + atom_count
        ]
    ):
        if len(line) < 44:
            raise RuntimeError(
                f"Malformed atom line "
                f"{index + 1} in {path}"
            )

        atoms.append(
            {
                "index": index,
                "resid": int(
                    line[0:5]
                ),
                "resname": line[
                    5:10
                ].strip(),
                "atomname": line[
                    10:15
                ].strip(),
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
                2 + atom_count
            ].split()
        ]
    )

    return (
        lines[0],
        atoms,
        box,
    )


def positions(
    atoms: list[dict[str, Any]],
) -> np.ndarray:
    return np.asarray(
        [
            atom["position"]
            for atom in atoms
        ],
        dtype=float,
    )


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


def section_name(line: str) -> str | None:
    match = re.match(
        r"^\s*\[\s*([^\]]+?)\s*\]\s*$",
        line,
    )

    if match is None:
        return None

    return (
        match.group(1)
        .strip()
        .lower()
    )


def first_data_line(
    lines: list[str],
    start: int,
    stop: int,
) -> int:
    for index in range(
        start,
        stop,
    ):
        stripped = lines[
            index
        ].strip()

        if (
            stripped
            and not stripped.startswith(";")
            and not stripped.startswith("#")
        ):
            return index

    raise RuntimeError(
        "No topology data line was found."
    )


def molecule_blocks(
    lines: list[str],
) -> dict[
    str,
    tuple[int, int],
]:
    starts = []

    for index, line in enumerate(lines):
        if section_name(line) == "moleculetype":
            starts.append(index)

    blocks: dict[
        str,
        tuple[int, int],
    ] = {}

    for block_number, start in enumerate(
        starts
    ):
        stop = (
            starts[
                block_number + 1
            ]
            if block_number + 1 < len(starts)
            else len(lines)
        )

        for index in range(
            start + 1,
            stop,
        ):
            current_section = section_name(
                lines[index]
            )

            if current_section in {
                "system",
                "molecules",
            }:
                stop = index
                break

        name_index = first_data_line(
            lines,
            start + 1,
            stop,
        )

        name = (
            lines[
                name_index
            ]
            .split(";", 1)[0]
            .split()[0]
        )

        blocks[name] = (
            start,
            stop,
        )

    return blocks


def replace_cap_atoms(
    lines: list[str],
    molecule_name: str,
    atom_count: int,
) -> list[str]:
    blocks = molecule_blocks(lines)

    if molecule_name not in blocks:
        raise RuntimeError(
            f"Molecule type {molecule_name} was not found."
        )

    block_start, block_stop = blocks[
        molecule_name
    ]

    atoms_header = None

    for index in range(
        block_start,
        block_stop,
    ):
        if section_name(
            lines[index]
        ) == "atoms":
            atoms_header = index
            break

    if atoms_header is None:
        raise RuntimeError(
            f"No [ atoms ] section for {molecule_name}."
        )

    atom_section_stop = block_stop

    for index in range(
        atoms_header + 1,
        block_stop,
    ):
        if section_name(
            lines[index]
        ) is not None:
            atom_section_stop = index
            break

    data_indices = []

    for index in range(
        atoms_header + 1,
        atom_section_stop,
    ):
        stripped = (
            lines[index]
            .split(";", 1)[0]
            .strip()
        )

        if stripped:
            data_indices.append(index)

    if not data_indices:
        raise RuntimeError(
            f"No CAP atom template in {molecule_name}."
        )

    template_tokens = (
        lines[
            data_indices[0]
        ]
        .split(";", 1)[0]
        .split()
    )

    if len(template_tokens) < 7:
        raise RuntimeError(
            "Unexpected CAP atom-table format."
        )

    atom_type = template_tokens[1]
    # Coordinate and topology atom names must match.
    # The interaction type remains CAP; the user-facing atom name
    # is standardized to CAP in both CAPL and CAPU.
    atom_name = "CAP"
    charge = template_tokens[6]

    mass = (
        template_tokens[7]
        if len(template_tokens) >= 8
        else None
    )

    new_atom_lines = []

    for atom_number in range(
        1,
        atom_count + 1,
    ):
        fields = [
            f"{atom_number:6d}",
            f"{atom_type:<10}",
            f"{1:6d}",
            f"{molecule_name:<8}",
            f"{atom_name:<8}",
            f"{atom_number:6d}",
            f"{charge:>14}",
        ]

        if mass is not None:
            fields.append(
                f"{mass:>14}"
            )

        new_atom_lines.append(
            " ".join(fields)
        )

    preserved_prefix = lines[
        atoms_header + 1:
        data_indices[0]
    ]

    preserved_suffix = lines[
        data_indices[-1] + 1:
        atom_section_stop
    ]

    return (
        lines[:atoms_header + 1]
        + preserved_prefix
        + new_atom_lines
        + preserved_suffix
        + lines[atom_section_stop:]
    )


def update_molecule_counts(
    lines: list[str],
) -> list[str]:
    molecules_header = None

    for index, line in enumerate(lines):
        if section_name(line) == "molecules":
            molecules_header = index
            break

    if molecules_header is None:
        raise RuntimeError(
            "No [ molecules ] section was found."
        )

    expected = {
        "SOL": WATERS,
        "CAPL": 1,
        "CAPU": 1,
    }

    found = {
        key: 0
        for key in expected
    }

    result = list(lines)

    for index in range(
        molecules_header + 1,
        len(result),
    ):
        stripped = (
            result[index]
            .split(";", 1)[0]
            .strip()
        )

        if (
            not stripped
            or stripped.startswith("#")
        ):
            continue

        tokens = stripped.split()

        if len(tokens) < 2:
            continue

        name = tokens[0]

        if name not in expected:
            continue

        comment = ""

        if ";" in result[index]:
            comment = (
                " ;"
                + result[index]
                .split(";", 1)[1]
            )

        result[index] = (
            f"{name:<12} "
            f"{expected[name]:8d}"
            f"{comment}"
        )

        found[name] += 1

    if found != {
        "SOL": 1,
        "CAPL": 1,
        "CAPU": 1,
    }:
        raise RuntimeError(
            "Unexpected molecule-count entries: "
            f"{found}"
        )

    return result


def absolutize_existing_includes(
    lines: list[str],
    source_parent: Path,
) -> list[str]:
    output = []

    pattern = re.compile(
        r'^(\s*#include\s+")([^"]+)(".*)$'
    )

    for line in lines:
        match = pattern.match(line)

        if match is None:
            output.append(line)
            continue

        include_path = Path(
            match.group(2)
        )

        if include_path.is_absolute():
            output.append(line)
            continue

        candidate = (
            source_parent
            / include_path
        ).resolve()

        if candidate.exists():
            output.append(
                match.group(1)
                + str(candidate)
                + match.group(3)
            )
        else:
            output.append(line)

    return output


def topology_molecule_counts(
    lines: list[str],
) -> dict[str, int]:
    header = None

    for index, line in enumerate(lines):
        if section_name(line) == "molecules":
            header = index
            break

    if header is None:
        raise RuntimeError(
            "No [ molecules ] section."
        )

    counts: dict[str, int] = {}

    for line in lines[
        header + 1:
    ]:
        if section_name(line) is not None:
            continue

        stripped = (
            line
            .split(";", 1)[0]
            .strip()
        )

        if (
            not stripped
            or stripped.startswith("#")
        ):
            continue

        tokens = stripped.split()

        if len(tokens) >= 2:
            try:
                counts[
                    tokens[0]
                ] = int(
                    tokens[1]
                )
            except ValueError:
                pass

    return counts


def cap_atom_count(
    lines: list[str],
    molecule_name: str,
) -> int:
    blocks = molecule_blocks(lines)

    if molecule_name not in blocks:
        raise RuntimeError(
            f"Missing {molecule_name} block."
        )

    block_start, block_stop = blocks[
        molecule_name
    ]

    atoms_header = None

    for index in range(
        block_start,
        block_stop,
    ):
        if section_name(
            lines[index]
        ) == "atoms":
            atoms_header = index
            break

    if atoms_header is None:
        raise RuntimeError(
            f"No atoms section for {molecule_name}."
        )

    count = 0

    for index in range(
        atoms_header + 1,
        block_stop,
    ):
        if section_name(
            lines[index]
        ) is not None:
            break

        stripped = (
            lines[index]
            .split(";", 1)[0]
            .strip()
        )

        if stripped:
            count += 1

    return count


def parse_cap_ow_override(
    lines: list[str],
) -> dict[str, float | str]:
    active_section = None
    matches = []

    for line in lines:
        current = section_name(line)

        if current is not None:
            active_section = current
            continue

        if active_section != "nonbond_params":
            continue

        stripped = (
            line
            .split(";", 1)[0]
            .strip()
        )

        if not stripped:
            continue

        tokens = stripped.split()

        if len(tokens) < 5:
            continue

        pair = {
            tokens[0].upper(),
            tokens[1].upper(),
        }

        if pair != {
            "CAP",
            "OW",
        }:
            continue

        matches.append(tokens)

    if len(matches) != 1:
        raise RuntimeError(
            "Expected exactly one CAP–OW nonbonded override; "
            f"found {len(matches)}."
        )

    tokens = matches[0]

    return {
        "type_i": tokens[0],
        "type_j": tokens[1],
        "function": int(
            tokens[2]
        ),
        "sigma_nm": float(
            tokens[3]
        ),
        "epsilon_kJ_mol": float(
            tokens[4]
        ),
    }


def write_index(
    path: Path,
) -> None:
    groups = {
        "System": range(
            1,
            EXPECTED_ATOMS + 1,
        ),
        "HBN": range(
            1,
            HBN_ATOMS + 1,
        ),
        "PYR": range(
            HBN_ATOMS + 1,
            SOLUTE_ATOMS + 1,
        ),
        "HBN_PYR": range(
            1,
            SOLUTE_ATOMS + 1,
        ),
        "SOL": range(
            SOLUTE_ATOMS + 1,
            CAPL_START + 1,
        ),
        "CAPL": range(
            CAPL_START + 1,
            CAPL_STOP + 1,
        ),
        "CAPU": range(
            CAPU_START + 1,
            CAPU_STOP + 1,
        ),
        "CAPS": range(
            CAPL_START + 1,
            CAPU_STOP + 1,
        ),
    }

    lines = []

    for name, indices in groups.items():
        lines.append(
            f"[ {name} ]"
        )

        row = []

        for atom_number in indices:
            row.append(
                f"{atom_number:6d}"
            )

            if len(row) == 15:
                lines.append(
                    " ".join(row)
                )

                row = []

        if row:
            lines.append(
                " ".join(row)
            )

        lines.append("")

    path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def water_oxygen_indices(
    atoms: list[dict[str, Any]],
) -> np.ndarray:
    indices = []

    for water_index in range(WATERS):
        start = (
            SOLUTE_ATOMS
            + water_index
            * WATER_SITES
        )

        chunk = atoms[
            start:
            start + WATER_SITES
        ]

        matches = [
            start + local_index
            for local_index, atom
            in enumerate(chunk)
            if atom["atomname"] == "OW"
        ]

        if len(matches) != 1:
            raise RuntimeError(
                "Could not identify exactly one OW "
                f"for water {water_index}."
            )

        indices.append(
            matches[0]
        )

    return np.asarray(
        indices,
        dtype=int,
    )


def analytic_cap_water(
    water_positions: np.ndarray,
    cap_positions: np.ndarray,
    box: np.ndarray,
    c12: float,
) -> dict[str, Any]:
    total_energy = 0.0

    water_forces = np.zeros_like(
        water_positions
    )

    minimum_distance = math.inf

    chunk_size = 256

    for start in range(
        0,
        len(water_positions),
        chunk_size,
    ):
        stop = min(
            start + chunk_size,
            len(water_positions),
        )

        displacement = (
            water_positions[
                start:stop,
                None,
                :,
            ]
            - cap_positions[
                None,
                :,
                :,
            ]
        )

        displacement = minimum_image(
            displacement,
            box,
        )

        distance_squared = np.sum(
            displacement * displacement,
            axis=2,
        )

        distances = np.sqrt(
            distance_squared
        )

        minimum_distance = min(
            minimum_distance,
            float(
                np.min(distances)
            ),
        )

        inverse_r12 = (
            distance_squared ** -6
        )

        total_energy += float(
            np.sum(
                c12
                * inverse_r12
            )
        )

        inverse_r14 = (
            distance_squared ** -7
        )

        pair_force_vectors = (
            12.0
            * c12
            * displacement
            * inverse_r14[
                :,
                :,
                None,
            ]
        )

        water_forces[
            start:stop
        ] = np.sum(
            pair_force_vectors,
            axis=1,
        )

    force_norms = np.linalg.norm(
        water_forces,
        axis=1,
    )

    return {
        "energy_kJ_mol": (
            total_energy
        ),
        "minimum_distance_nm": (
            minimum_distance
        ),
        "maximum_water_force_kJ_mol_nm": float(
            np.max(
                force_norms
            )
        ),
        "mean_water_force_kJ_mol_nm": float(
            np.mean(
                force_norms
            )
        ),
    }


def orthogonal_basis(
    axis: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:
    axis = (
        axis
        / np.linalg.norm(axis)
    )

    trial = np.asarray(
        [1.0, 0.0, 0.0],
        dtype=float,
    )

    if abs(
        np.dot(
            axis,
            trial,
        )
    ) > 0.90:
        trial = np.asarray(
            [0.0, 1.0, 0.0],
            dtype=float,
        )

    first = np.cross(
        axis,
        trial,
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


def probe_energy(
    points: np.ndarray,
    cap_positions: np.ndarray,
    box: np.ndarray,
    c12: float,
) -> np.ndarray:
    displacement = (
        points[
            :,
            None,
            :,
        ]
        - cap_positions[
            None,
            :,
            :,
        ]
    )

    displacement = minimum_image(
        displacement,
        box,
    )

    distance_squared = np.sum(
        displacement
        * displacement,
        axis=2,
    )

    return np.sum(
        c12
        * distance_squared ** -6,
        axis=1,
    )


def first_threshold_crossing(
    coordinate: np.ndarray,
    values: np.ndarray,
    threshold: float,
) -> float:
    if values[0] > threshold:
        return 0.0

    above = np.flatnonzero(
        values >= threshold
    )

    if len(above) == 0:
        return float(
            coordinate[-1]
        )

    index = int(
        above[0]
    )

    if index == 0:
        return float(
            coordinate[0]
        )

    x0 = float(
        coordinate[
            index - 1
        ]
    )

    x1 = float(
        coordinate[index]
    )

    y0 = float(
        values[
            index - 1
        ]
    )

    y1 = float(
        values[index]
    )

    if math.isclose(
        y1,
        y0,
        rel_tol=0.0,
        abs_tol=1.0e-30,
    ):
        return x1

    fraction = (
        threshold - y0
    ) / (
        y1 - y0
    )

    return (
        x0
        + fraction
        * (
            x1 - x0
        )
    )


def aperture_profiles(
    *,
    center: np.ndarray,
    axis: np.ndarray,
    lower_caps: np.ndarray,
    upper_caps: np.ndarray,
    box: np.ndarray,
    c12: float,
    threshold_kJ_mol: float,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, float],
]:
    axis = (
        axis
        / np.linalg.norm(axis)
    )

    radial_basis, _ = orthogonal_basis(
        axis
    )

    center_relative_lower = minimum_image(
        lower_caps - center,
        box,
    )

    center_relative_upper = minimum_image(
        upper_caps - center,
        box,
    )

    lower_axial = float(
        np.mean(
            center_relative_lower
            @ axis
        )
    )

    upper_axial = float(
        np.mean(
            center_relative_upper
            @ axis
        )
    )

    lower_plane_center = (
        center
        + lower_axial
        * axis
    )

    upper_plane_center = (
        center
        + upper_axial
        * axis
    )

    radial_coordinates = np.linspace(
        0.0,
        RADIAL_SCAN_MAX_NM,
        RADIAL_SCAN_POINTS,
    )

    lower_radial_points = (
        lower_plane_center[
            None,
            :,
        ]
        + radial_coordinates[
            :,
            None,
        ]
        * radial_basis[
            None,
            :,
        ]
    )

    upper_radial_points = (
        upper_plane_center[
            None,
            :,
        ]
        + radial_coordinates[
            :,
            None,
        ]
        * radial_basis[
            None,
            :,
        ]
    )

    lower_radial_energy = probe_energy(
        lower_radial_points,
        lower_caps,
        box,
        c12,
    )

    upper_radial_energy = probe_energy(
        upper_radial_points,
        upper_caps,
        box,
        c12,
    )

    lower_aperture_radius = first_threshold_crossing(
        radial_coordinates,
        lower_radial_energy,
        threshold_kJ_mol,
    )

    upper_aperture_radius = first_threshold_crossing(
        radial_coordinates,
        upper_radial_energy,
        threshold_kJ_mol,
    )

    radial_rows = []

    for index, radius in enumerate(
        radial_coordinates
    ):
        radial_rows.append(
            {
                "radius_nm": float(radius),
                "lower_energy_kJ_mol": float(
                    lower_radial_energy[index]
                ),
                "lower_energy_kBT": float(
                    lower_radial_energy[index]
                    / KBT_KJ_MOL
                ),
                "upper_energy_kJ_mol": float(
                    upper_radial_energy[index]
                ),
                "upper_energy_kBT": float(
                    upper_radial_energy[index]
                    / KBT_KJ_MOL
                ),
            }
        )

    axial_offsets = np.linspace(
        -AXIAL_SCAN_HALF_WIDTH_NM,
        AXIAL_SCAN_HALF_WIDTH_NM,
        AXIAL_SCAN_POINTS,
    )

    lower_axial_points = (
        lower_plane_center[
            None,
            :,
        ]
        + axial_offsets[
            :,
            None,
        ]
        * axis[
            None,
            :,
        ]
    )

    upper_axial_points = (
        upper_plane_center[
            None,
            :,
        ]
        + axial_offsets[
            :,
            None,
        ]
        * axis[
            None,
            :,
        ]
    )

    lower_axial_energy = probe_energy(
        lower_axial_points,
        lower_caps,
        box,
        c12,
    )

    upper_axial_energy = probe_energy(
        upper_axial_points,
        upper_caps,
        box,
        c12,
    )

    axial_rows = []

    for index, offset in enumerate(
        axial_offsets
    ):
        axial_rows.append(
            {
                "offset_from_cap_plane_nm": float(
                    offset
                ),
                "lower_centerline_energy_kJ_mol": float(
                    lower_axial_energy[index]
                ),
                "lower_centerline_energy_kBT": float(
                    lower_axial_energy[index]
                    / KBT_KJ_MOL
                ),
                "upper_centerline_energy_kJ_mol": float(
                    upper_axial_energy[index]
                ),
                "upper_centerline_energy_kBT": float(
                    upper_axial_energy[index]
                    / KBT_KJ_MOL
                ),
            }
        )

    metrics = {
        "lower_cap_axial_coordinate_nm": (
            lower_axial
        ),
        "upper_cap_axial_coordinate_nm": (
            upper_axial
        ),
        "lower_actual_aperture_radius_5kBT_nm": (
            lower_aperture_radius
        ),
        "upper_actual_aperture_radius_5kBT_nm": (
            upper_aperture_radius
        ),
        "mean_actual_aperture_radius_5kBT_nm": (
            0.5
            * (
                lower_aperture_radius
                + upper_aperture_radius
            )
        ),
        "lower_centerline_maximum_barrier_kBT": float(
            np.max(
                lower_axial_energy
            )
            / KBT_KJ_MOL
        ),
        "upper_centerline_maximum_barrier_kBT": float(
            np.max(
                upper_axial_energy
            )
            / KBT_KJ_MOL
        ),
    }

    return (
        radial_rows,
        axial_rows,
        metrics,
    )


def parse_tpr_atoms(
    dump_text: str,
) -> int:
    match = re.search(
        r"^\s*natoms\s*=\s*(\d+)",
        dump_text,
        flags=re.MULTILINE,
    )

    if match is None:
        raise RuntimeError(
            "Could not parse TPR atom count."
        )

    return int(
        match.group(1)
    )


def parse_energy_menu(
    output: str,
) -> dict[str, int]:
    terms: dict[str, int] = {}

    pattern = re.compile(
        r"(?<!\S)(\d+)\s+"
        r"([A-Za-z][A-Za-z0-9_.:+\-()]*)"
    )

    for line in output.splitlines():
        for match in pattern.finditer(line):
            terms[
                match.group(2)
            ] = int(
                match.group(1)
            )

    return terms


def resolve_energy_term(
    terms: dict[str, int],
    exact: tuple[str, ...],
    required_tokens: tuple[str, ...],
) -> tuple[
    str,
    int,
]:
    for candidate in exact:
        if candidate in terms:
            return (
                candidate,
                terms[candidate],
            )

    for name, number in terms.items():
        if all(
            token in name
            for token in required_tokens
        ):
            return (
                name,
                number,
            )

    raise RuntimeError(
        "Could not resolve energy term. "
        f"Exact={exact}; tokens={required_tokens}"
    )


def extract_single_energy(
    gmx: str,
    edr: Path,
    term_number: int,
    output_path: Path,
) -> float:
    completed = run_command(
        [
            gmx,
            "energy",
            "-f",
            str(edr),
            "-o",
            str(output_path),
            "-xvg",
            "none",
            "-dp",
        ],
        cwd=OUTPUT_ROOT,
        input_text=(
            f"{term_number}\n"
            "0\n"
        ),
    )

    if (
        completed.returncode != 0
        or not output_path.exists()
        or output_path.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Could not extract energy term {term_number}."
        )

    rows = []

    for line in output_path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        stripped = line.strip()

        if (
            not stripped
            or stripped.startswith("#")
            or stripped.startswith("@")
        ):
            continue

        fields = stripped.split()

        if len(fields) >= 2:
            rows.append(
                (
                    float(fields[0]),
                    float(fields[1]),
                )
            )

    if not rows:
        raise RuntimeError(
            f"No energy values in {output_path}."
        )

    return float(
        rows[-1][1]
    )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    SELECTED_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        R1_SELECTED_TOPOLOGY,
        R1_PROTOTYPE_JSON,
        R2_GEOMETRY_SUMMARY,
        R2_DEFINITION_JSON,
        R2_SYSTEM_GRO,
        R2_CAPS_GRO,
    ):
        require_file(required)

    geometry_summary = read_single_csv_row(
        R2_GEOMETRY_SUMMARY
    )

    if (
        geometry_summary.get(
            "decision"
        )
        !=
        "R2_PARTIAL_CAP_GEOMETRY_STATIC_GATE_PASSED"
    ):
        raise RuntimeError(
            "R2 geometry has not passed its static gate."
        )

    if not parse_bool(
        geometry_summary.get(
            "topology_generation_authorized",
            "false",
        )
    ):
        raise RuntimeError(
            "R2 topology generation is not authorized."
        )

    r2_definition = json.loads(
        R2_DEFINITION_JSON.read_text(
            encoding="utf-8"
        )
    )

    r1_prototype = json.loads(
        R1_PROTOTYPE_JSON.read_text(
            encoding="utf-8"
        )
    )

    _, system_atoms, box = read_gro(
        R2_SYSTEM_GRO
    )

    _, cap_atoms, cap_box = read_gro(
        R2_CAPS_GRO
    )

    if len(system_atoms) != EXPECTED_ATOMS:
        raise RuntimeError(
            "Unexpected R2 atom count: "
            f"{len(system_atoms)}/"
            f"{EXPECTED_ATOMS}"
        )

    if len(cap_atoms) != TOTAL_CAPS:
        raise RuntimeError(
            "Unexpected R2 cap count: "
            f"{len(cap_atoms)}/"
            f"{TOTAL_CAPS}"
        )

    if not np.allclose(
        box,
        cap_box,
        rtol=0.0,
        atol=1.0e-6,
    ):
        raise RuntimeError(
            "R2 system and cap boxes differ."
        )

    source_topology_lines = (
        R1_SELECTED_TOPOLOGY.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines()
    )

    topology_lines = replace_cap_atoms(
        source_topology_lines,
        "CAPL",
        CAPS_PER_END,
    )

    topology_lines = replace_cap_atoms(
        topology_lines,
        "CAPU",
        CAPS_PER_END,
    )

    topology_lines = update_molecule_counts(
        topology_lines
    )

    topology_lines = absolutize_existing_includes(
        topology_lines,
        R1_SELECTED_TOPOLOGY.parent,
    )

    R2_TOPOLOGY.write_text(
        "\n".join(
            topology_lines
        )
        + "\n",
        encoding="utf-8",
    )

    generated_lines = R2_TOPOLOGY.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    molecule_counts = topology_molecule_counts(
        generated_lines
    )

    capl_topology_atoms = cap_atom_count(
        generated_lines,
        "CAPL",
    )

    capu_topology_atoms = cap_atom_count(
        generated_lines,
        "CAPU",
    )

    cap_ow_override = parse_cap_ow_override(
        generated_lines
    )

    sigma_nm = float(
        cap_ow_override[
            "sigma_nm"
        ]
    )

    epsilon_kJ_mol = float(
        cap_ow_override[
            "epsilon_kJ_mol"
        ]
    )

    c12_selected = (
        4.0
        * epsilon_kJ_mol
        * abs(
            sigma_nm
        ) ** 12
    )

    write_index(
        R2_INDEX
    )

    R2_STATIC_MDP.write_text(
        """integrator               = md
dt                       = 0.0005
nsteps                   = 0
continuation             = yes
gen-vel                  = no

cutoff-scheme            = Verlet
nstlist                  = 20
rlist                    = 1.20

coulombtype              = PME
rcoulomb                 = 1.20
fourierspacing           = 0.12
pme-order                = 4
ewald-rtol               = 1e-5

vdwtype                  = Cut-off
rvdw                     = 1.20
vdw-modifier             = Potential-shift
DispCorr                 = EnerPres

constraints              = h-bonds
constraint-algorithm     = lincs
lincs-iter               = 1
lincs-order              = 4

pbc                      = xyz

nstlog                   = 1
nstenergy                = 1
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstxout-compressed       = 0

energygrps               = HBN_PYR CAPS SOL
""",
        encoding="utf-8",
    )

    gmx = locate_gmx()

    grompp = run_command(
        [
            gmx,
            "grompp",
            "-f",
            str(R2_STATIC_MDP),
            "-c",
            str(R2_SYSTEM_GRO),
            "-p",
            str(R2_TOPOLOGY),
            "-n",
            str(R2_INDEX),
            "-o",
            str(R2_TPR),
            "-po",
            str(
                OUTPUT_ROOT
                / "r2_static_single_point_processed.mdp"
            ),
            "-maxwarn",
            "0",
        ],
        cwd=OUTPUT_ROOT,
    )

    GROMPP_LOG.write_text(
        grompp.stdout,
        encoding="utf-8",
    )

    if (
        grompp.returncode != 0
        or not R2_TPR.exists()
        or R2_TPR.stat().st_size == 0
    ):
        raise RuntimeError(
            "R2 grompp failed. "
            f"See {GROMPP_LOG}"
        )

    dump = subprocess.run(
        [
            gmx,
            "dump",
            "-s",
            str(R2_TPR),
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
            "Could not dump the R2 TPR."
        )

    tpr_atoms = parse_tpr_atoms(
        dump.stdout
    )

    mdrun_return_code = run_live(
        [
            gmx,
            "mdrun",
            "-s",
            str(R2_TPR),
            "-deffnm",
            str(DEFFNM),
            "-ntmpi",
            "1",
            "-ntomp",
            "4",
        ],
        cwd=OUTPUT_ROOT,
        log_path=MDRUN_CONSOLE,
    )

    edr = Path(
        str(DEFFNM)
        + ".edr"
    )

    log = Path(
        str(DEFFNM)
        + ".log"
    )

    # A zero-step mdrun is a single-point energy evaluation.
    # GROMACS is not required to write a final-coordinate GRO
    # because no integration step or coordinate update occurs.
    # The evaluated coordinates remain those in R2_SYSTEM_GRO.
    for product in (
        edr,
        log,
    ):
        require_file(product)

    menu_probe = run_command(
        [
            gmx,
            "energy",
            "-f",
            str(edr),
            "-o",
            str(
                OUTPUT_ROOT
                / "r2_energy_menu_probe.xvg"
            ),
        ],
        cwd=OUTPUT_ROOT,
        input_text="0\n",
    )

    ENERGY_MENU_LOG.write_text(
        menu_probe.stdout,
        encoding="utf-8",
    )

    energy_terms = parse_energy_menu(
        menu_probe.stdout
    )

    term_requests = {
        "Potential": (
            (
                "Potential",
            ),
            (
                "Potential",
            ),
        ),
        "LJ_CAPS_SOL": (
            (
                "LJ-SR:CAPS-SOL",
                "LJ-SR:SOL-CAPS",
            ),
            (
                "LJ-SR",
                "CAPS",
                "SOL",
            ),
        ),
        "Coul_CAPS_SOL": (
            (
                "Coul-SR:CAPS-SOL",
                "Coul-SR:SOL-CAPS",
            ),
            (
                "Coul-SR",
                "CAPS",
                "SOL",
            ),
        ),
        "LJ_HBNPYR_CAPS": (
            (
                "LJ-SR:HBN_PYR-CAPS",
                "LJ-SR:CAPS-HBN_PYR",
            ),
            (
                "LJ-SR",
                "HBN_PYR",
                "CAPS",
            ),
        ),
        "Coul_HBNPYR_CAPS": (
            (
                "Coul-SR:HBN_PYR-CAPS",
                "Coul-SR:CAPS-HBN_PYR",
            ),
            (
                "Coul-SR",
                "HBN_PYR",
                "CAPS",
            ),
        ),
    }

    resolved_terms: dict[
        str,
        tuple[str, int],
    ] = {}

    for label, (
        exact,
        tokens,
    ) in term_requests.items():
        resolved_terms[label] = (
            resolve_energy_term(
                energy_terms,
                exact,
                tokens,
            )
        )

    gromacs_energies = {}

    for label, (
        term_name,
        term_number,
    ) in resolved_terms.items():
        output_xvg = (
            OUTPUT_ROOT
            / (
                "r2_static_"
                + label.lower()
                + ".xvg"
            )
        )

        value = extract_single_energy(
            gmx,
            edr,
            term_number,
            output_xvg,
        )

        gromacs_energies[
            label
        ] = {
            "term_name": term_name,
            "term_number": term_number,
            "value_kJ_mol": value,
        }

    write_csv(
        STATIC_ENERGY_CSV,
        [
            {
                "label": label,
                **values,
            }
            for label, values
            in gromacs_energies.items()
        ],
    )

    system_positions = positions(
        system_atoms
    )

    oxygen_indices = water_oxygen_indices(
        system_atoms
    )

    water_positions = (
        system_positions[
            oxygen_indices
        ]
    )

    cap_positions = (
        system_positions[
            CAPL_START:
            CAPU_STOP
        ]
    )

    lower_caps = (
        cap_positions[
            :CAPS_PER_END
        ]
    )

    upper_caps = (
        cap_positions[
            CAPS_PER_END:
        ]
    )

    selected_analytic = analytic_cap_water(
        water_positions,
        cap_positions,
        box,
        c12_selected,
    )

    center = np.asarray(
        r1_prototype[
            "tube_center_wrapped_nm"
        ],
        dtype=float,
    )

    axis = np.asarray(
        r1_prototype[
            "tube_axis"
        ],
        dtype=float,
    )

    axis /= np.linalg.norm(
        axis
    )

    (
        radial_rows,
        axial_rows,
        aperture_metrics,
    ) = aperture_profiles(
        center=center,
        axis=axis,
        lower_caps=lower_caps,
        upper_caps=upper_caps,
        box=box,
        c12=c12_selected,
        threshold_kJ_mol=(
            SELECTED_STRENGTH_KBT
            * KBT_KJ_MOL
        ),
    )

    write_csv(
        RADIAL_PROFILE_CSV,
        radial_rows,
    )

    write_csv(
        AXIAL_PROFILE_CSV,
        axial_rows,
    )

    strength_rows = []

    for strength_kbt in STRENGTHS_KBT:
        epsilon = (
            strength_kbt
            * KBT_KJ_MOL
            / 4.0
        )

        c12 = (
            4.0
            * epsilon
            * abs(
                sigma_nm
            ) ** 12
        )

        analytic = analytic_cap_water(
            water_positions,
            cap_positions,
            box,
            c12,
        )

        (
            _,
            _,
            strength_aperture,
        ) = aperture_profiles(
            center=center,
            axis=axis,
            lower_caps=lower_caps,
            upper_caps=upper_caps,
            box=box,
            c12=c12,
            threshold_kJ_mol=(
                SELECTED_STRENGTH_KBT
                * KBT_KJ_MOL
            ),
        )

        strength_rows.append(
            {
                "strength_kBT_at_sigma": (
                    strength_kbt
                ),
                "epsilon_kJ_mol": (
                    epsilon
                ),
                "C12_kJ_mol_nm12": (
                    c12
                ),
                "initial_CAP_OW_energy_kJ_mol": (
                    analytic[
                        "energy_kJ_mol"
                    ]
                ),
                "minimum_CAP_OW_distance_nm": (
                    analytic[
                        "minimum_distance_nm"
                    ]
                ),
                "maximum_water_force_kJ_mol_nm": (
                    analytic[
                        "maximum_water_force_kJ_mol_nm"
                    ]
                ),
                "lower_actual_5kBT_aperture_radius_nm": (
                    strength_aperture[
                        "lower_actual_aperture_radius_5kBT_nm"
                    ]
                ),
                "upper_actual_5kBT_aperture_radius_nm": (
                    strength_aperture[
                        "upper_actual_aperture_radius_5kBT_nm"
                    ]
                ),
                "lower_centerline_maximum_barrier_kBT": (
                    strength_aperture[
                        "lower_centerline_maximum_barrier_kBT"
                    ]
                ),
                "upper_centerline_maximum_barrier_kBT": (
                    strength_aperture[
                        "upper_centerline_maximum_barrier_kBT"
                    ]
                ),
            }
        )

    write_csv(
        STRENGTH_SCAN_CSV,
        strength_rows,
    )

    gromacs_cap_sol = float(
        gromacs_energies[
            "LJ_CAPS_SOL"
        ][
            "value_kJ_mol"
        ]
    )

    analytic_cap_sol = float(
        selected_analytic[
            "energy_kJ_mol"
        ]
    )

    absolute_error = abs(
        gromacs_cap_sol
        - analytic_cap_sol
    )

    relative_error = (
        absolute_error
        / max(
            abs(
                analytic_cap_sol
            ),
            1.0e-12,
        )
    )

    mdrun_text = (
        MDRUN_CONSOLE.read_text(
            encoding="utf-8",
            errors="replace",
        )
        + "\n"
        + log.read_text(
            encoding="utf-8",
            errors="replace",
        )
    )

    instability_patterns = (
        r"\bnan\b",
        r"fatal error",
        r"segmentation fault",
        r"lincs warning",
        r"constraint warning",
        r"water molecule cannot be settled",
        r"water molecule can not be settled",
    )

    instability_hits = [
        pattern
        for pattern in instability_patterns
        if re.search(
            pattern,
            mdrun_text,
            flags=re.IGNORECASE,
        )
    ]

    conventional_finish_marker = bool(
        re.search(
            r"Finished mdrun",
            mdrun_text,
            flags=re.IGNORECASE,
        )
    )

    zero_step_static_completion = (
        mdrun_return_code == 0
        and edr.is_file()
        and edr.stat().st_size > 0
        and log.is_file()
        and log.stat().st_size > 0
        and bool(
            re.search(
                r"\b0\s+steps,\s+0\.0\s+ps\.",
                mdrun_text,
                flags=re.IGNORECASE,
            )
        )
    )

    finished = (
        conventional_finish_marker
        or zero_step_static_completion
    )

    topology_audit_rows = [
        {
            "item": "HBN_molecules",
            "observed": molecule_counts.get(
                "HBN",
                "",
            ),
            "expected": 1,
        },
        {
            "item": "PYR_molecules",
            "observed": molecule_counts.get(
                "PYR",
                "",
            ),
            "expected": 4,
        },
        {
            "item": "SOL_molecules",
            "observed": molecule_counts.get(
                "SOL",
                "",
            ),
            "expected": WATERS,
        },
        {
            "item": "CAPL_molecules",
            "observed": molecule_counts.get(
                "CAPL",
                "",
            ),
            "expected": 1,
        },
        {
            "item": "CAPU_molecules",
            "observed": molecule_counts.get(
                "CAPU",
                "",
            ),
            "expected": 1,
        },
        {
            "item": "CAPL_atoms",
            "observed": capl_topology_atoms,
            "expected": CAPS_PER_END,
        },
        {
            "item": "CAPU_atoms",
            "observed": capu_topology_atoms,
            "expected": CAPS_PER_END,
        },
        {
            "item": "CAP_OW_sigma_nm",
            "observed": sigma_nm,
            "expected": SIGMA_CAP_OW_NM,
        },
        {
            "item": "CAP_OW_epsilon_kJ_mol",
            "observed": epsilon_kJ_mol,
            "expected": (
                REFERENCE_R1_EPSILON_KJ_MOL
            ),
        },
        {
            "item": "TPR_atoms",
            "observed": tpr_atoms,
            "expected": EXPECTED_ATOMS,
        },
    ]

    write_csv(
        TOPOLOGY_AUDIT_CSV,
        topology_audit_rows,
    )

    lower_aperture = float(
        aperture_metrics[
            "lower_actual_aperture_radius_5kBT_nm"
        ]
    )

    upper_aperture = float(
        aperture_metrics[
            "upper_actual_aperture_radius_5kBT_nm"
        ]
    )

    maximum_centerline_barrier = max(
        float(
            aperture_metrics[
                "lower_centerline_maximum_barrier_kBT"
            ]
        ),
        float(
            aperture_metrics[
                "upper_centerline_maximum_barrier_kBT"
            ]
        ),
    )

    gates = {
        "R2_geometry_gate_passed": (
            geometry_summary.get(
                "decision"
            )
            ==
            "R2_PARTIAL_CAP_GEOMETRY_STATIC_GATE_PASSED"
        ),
        "R2_geometry_authorized_topology": (
            parse_bool(
                geometry_summary.get(
                    "topology_generation_authorized",
                    "false",
                )
            )
        ),
        "R2_GRO_has_68332_atoms": (
            len(system_atoms)
            == EXPECTED_ATOMS
        ),
        "R2_caps_GRO_has_288_beads": (
            len(cap_atoms)
            == TOTAL_CAPS
        ),
        "CAPL_topology_has_144_atoms": (
            capl_topology_atoms
            == CAPS_PER_END
        ),
        "CAPU_topology_has_144_atoms": (
            capu_topology_atoms
            == CAPS_PER_END
        ),
        "SOL_count_is_16565": (
            molecule_counts.get(
                "SOL"
            )
            == WATERS
        ),
        "CAPL_and_CAPU_counts_are_one": (
            molecule_counts.get(
                "CAPL"
            )
            == 1
            and molecule_counts.get(
                "CAPU"
            )
            == 1
        ),
        "CAP_OW_function_is_one": (
            int(
                cap_ow_override[
                    "function"
                ]
            )
            == 1
        ),
        "CAP_OW_sigma_is_negative_0p17nm": (
            math.isclose(
                sigma_nm,
                SIGMA_CAP_OW_NM,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
        ),
        "CAP_OW_epsilon_matches_R1_selected_model": (
            math.isclose(
                epsilon_kJ_mol,
                REFERENCE_R1_EPSILON_KJ_MOL,
                rel_tol=0.0,
                abs_tol=1.0e-9,
            )
        ),
        "grompp_return_code_zero": (
            grompp.returncode == 0
        ),
        "TPR_dump_return_code_zero": (
            dump.returncode == 0
        ),
        "TPR_has_68332_atoms": (
            tpr_atoms
            == EXPECTED_ATOMS
        ),
        "mdrun_return_code_zero": (
            mdrun_return_code == 0
        ),
        "mdrun_finished": finished,
        "no_instability_signatures": (
            len(
                instability_hits
            )
            == 0
        ),
        "minimum_CAP_OW_distance_is_valid": (
            float(
                selected_analytic[
                    "minimum_distance_nm"
                ]
            )
            >= MIN_CAP_OW_DISTANCE_NM
        ),
        "selected_CAP_SOL_energy_is_finite": (
            math.isfinite(
                analytic_cap_sol
            )
            and math.isfinite(
                gromacs_cap_sol
            )
        ),
        "selected_CAP_SOL_energy_is_below_100kJmol": (
            gromacs_cap_sol
            <= MAX_SELECTED_CAP_SOL_LJ_KJ_MOL
        ),
        "selected_maximum_water_force_is_below_250": (
            float(
                selected_analytic[
                    "maximum_water_force_kJ_mol_nm"
                ]
            )
            <= MAX_SELECTED_WATER_FORCE_KJ_MOL_NM
        ),
        "GROMACS_and_analytic_CAP_SOL_agree": (
            absolute_error
            <= MAX_GROMACS_ANALYTIC_ABSOLUTE_ERROR_KJ_MOL
            or relative_error
            <= MAX_GROMACS_ANALYTIC_RELATIVE_ERROR
        ),
        "CAP_SOL_Coulomb_is_zero": (
            abs(
                float(
                    gromacs_energies[
                        "Coul_CAPS_SOL"
                    ][
                        "value_kJ_mol"
                    ]
                )
            )
            <= ZERO_ENERGY_TOLERANCE_KJ_MOL
        ),
        "CAP_HBNPYR_LJ_is_zero": (
            abs(
                float(
                    gromacs_energies[
                        "LJ_HBNPYR_CAPS"
                    ][
                        "value_kJ_mol"
                    ]
                )
            )
            <= ZERO_ENERGY_TOLERANCE_KJ_MOL
        ),
        "CAP_HBNPYR_Coulomb_is_zero": (
            abs(
                float(
                    gromacs_energies[
                        "Coul_HBNPYR_CAPS"
                    ][
                        "value_kJ_mol"
                    ]
                )
            )
            <= ZERO_ENERGY_TOLERANCE_KJ_MOL
        ),
        "lower_actual_aperture_radius_is_0p20_to_0p45nm": (
            MIN_ACTUAL_APERTURE_RADIUS_NM
            <= lower_aperture
            <= MAX_ACTUAL_APERTURE_RADIUS_NM
        ),
        "upper_actual_aperture_radius_is_0p20_to_0p45nm": (
            MIN_ACTUAL_APERTURE_RADIUS_NM
            <= upper_aperture
            <= MAX_ACTUAL_APERTURE_RADIUS_NM
        ),
        "lower_upper_aperture_radii_are_symmetric": (
            abs(
                lower_aperture
                - upper_aperture
            )
            <= MAX_LOWER_UPPER_APERTURE_DIFFERENCE_NM
        ),
        "central_axial_path_barrier_is_at_most_5kBT": (
            maximum_centerline_barrier
            <= MAX_CENTERLINE_BARRIER_KBT
        ),
    }

    failed_gates = [
        name
        for name, passed
        in gates.items()
        if not passed
    ]

    accepted = (
        len(failed_gates) == 0
    )

    decision = (
        "R2_TOPOLOGY_AND_STATIC_CAP_WATER_MODEL_VALIDATED"
        if accepted
        else
        "R2_TOPOLOGY_STATIC_SCAN_REQUIRES_REVIEW"
    )

    required_next_step = (
        "RUN_R2_WATER_ONLY_ENERGY_MINIMIZATION"
        if accepted
        else
        "REVIEW_R2_TOPOLOGY_STATIC_SCAN_GATE_FAILURES"
    )

    summary = {
        "decision": decision,
        "R2_atoms": (
            len(system_atoms)
        ),
        "R2_waters": WATERS,
        "CAPL_atoms": (
            capl_topology_atoms
        ),
        "CAPU_atoms": (
            capu_topology_atoms
        ),
        "CAP_total_atoms": TOTAL_CAPS,
        "CAP_OW_sigma_nm": sigma_nm,
        "CAP_OW_epsilon_kJ_mol": (
            epsilon_kJ_mol
        ),
        "CAP_OW_C12_kJ_mol_nm12": (
            c12_selected
        ),
        "minimum_CAP_OW_distance_nm": (
            selected_analytic[
                "minimum_distance_nm"
            ]
        ),
        "analytic_CAP_SOL_LJ_kJ_mol": (
            analytic_cap_sol
        ),
        "GROMACS_CAP_SOL_LJ_kJ_mol": (
            gromacs_cap_sol
        ),
        "GROMACS_analytic_absolute_error_kJ_mol": (
            absolute_error
        ),
        "GROMACS_analytic_relative_error": (
            relative_error
        ),
        "maximum_water_force_kJ_mol_nm": (
            selected_analytic[
                "maximum_water_force_kJ_mol_nm"
            ]
        ),
        "CAP_SOL_Coulomb_kJ_mol": (
            gromacs_energies[
                "Coul_CAPS_SOL"
            ][
                "value_kJ_mol"
            ]
        ),
        "CAP_HBNPYR_LJ_kJ_mol": (
            gromacs_energies[
                "LJ_HBNPYR_CAPS"
            ][
                "value_kJ_mol"
            ]
        ),
        "CAP_HBNPYR_Coulomb_kJ_mol": (
            gromacs_energies[
                "Coul_HBNPYR_CAPS"
            ][
                "value_kJ_mol"
            ]
        ),
        "lower_actual_aperture_radius_5kBT_nm": (
            lower_aperture
        ),
        "upper_actual_aperture_radius_5kBT_nm": (
            upper_aperture
        ),
        "mean_actual_aperture_radius_5kBT_nm": (
            aperture_metrics[
                "mean_actual_aperture_radius_5kBT_nm"
            ]
        ),
        "lower_centerline_maximum_barrier_kBT": (
            aperture_metrics[
                "lower_centerline_maximum_barrier_kBT"
            ]
        ),
        "upper_centerline_maximum_barrier_kBT": (
            aperture_metrics[
                "upper_centerline_maximum_barrier_kBT"
            ]
        ),
        "failed_gates": (
            " | ".join(
                failed_gates
            )
        ),
        "water_only_energy_minimization_authorized": (
            accepted
        ),
        "short_frozen_solute_NVT_authorized": False,
        "long_mobile_MD_authorized": False,
        "multitemperature_MD_authorized": False,
        "QM_recalculation_authorized": False,
        "required_next_step": (
            required_next_step
        ),
    }

    write_csv(
        SUMMARY_CSV,
        [summary],
    )

    write_csv(
        GATE_CSV,
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
        for name, passed in gates.items()
    )

    REPORT_MD.write_text(
        f"""# R2 Topology and Static CAP–Water Scan

## Scope

The R2 symmetric partial-cap geometry was converted into a complete
GROMACS topology by reusing the validated R1 neutral CAP model.

No molecular dynamics or energy minimization was performed. The
GROMACS run contained zero integration steps and served only as an
independent static energy evaluation.

## System

- Total atoms:
  **{len(system_atoms)}**
- Waters:
  **{WATERS}**
- CAPL/CAPU beads:
  **{capl_topology_atoms}/{capu_topology_atoms}**
- Total cap beads:
  **{TOTAL_CAPS}**

## Reused CAP–OW model

- Sigma:
  **{sigma_nm:.12f} nm**
- Epsilon:
  **{epsilon_kJ_mol:.12f} kJ/mol**
- C12:
  **{c12_selected:.12e} kJ mol^-1 nm^12**

The negative sigma suppresses the C6 term under combination-rule 2,
leaving the validated purely repulsive C12/r12 interaction.

## Static initial-state results

- Minimum CAP–OW distance:
  **{selected_analytic['minimum_distance_nm']:.6f} nm**
- Analytic CAP–SOL LJ:
  **{analytic_cap_sol:.9f} kJ/mol**
- GROMACS CAP–SOL LJ:
  **{gromacs_cap_sol:.9f} kJ/mol**
- Absolute difference:
  **{absolute_error:.9e} kJ/mol**
- Relative difference:
  **{relative_error:.9e}**
- Maximum total CAP force on one water oxygen:
  **{selected_analytic['maximum_water_force_kJ_mol_nm']:.6f}
  kJ mol^-1 nm^-1**
- CAP–SOL Coulomb:
  **{float(gromacs_energies['Coul_CAPS_SOL']['value_kJ_mol']):.9f}
  kJ/mol**
- CAP–HBN/PYR LJ:
  **{float(gromacs_energies['LJ_HBNPYR_CAPS']['value_kJ_mol']):.9f}
  kJ/mol**
- CAP–HBN/PYR Coulomb:
  **{float(gromacs_energies['Coul_HBNPYR_CAPS']['value_kJ_mol']):.9f}
  kJ/mol**

## Actual potential-defined aperture

At the selected 5 kBT CAP–OW model:

- Lower aperture radius:
  **{lower_aperture:.6f} nm**
- Upper aperture radius:
  **{upper_aperture:.6f} nm**
- Mean aperture diameter:
  **{2.0 * float(aperture_metrics['mean_actual_aperture_radius_5kBT_nm']):.6f}
  nm**
- Lower maximum centerline barrier:
  **{float(aperture_metrics['lower_centerline_maximum_barrier_kBT']):.6f}
  kBT**
- Upper maximum centerline barrier:
  **{float(aperture_metrics['upper_centerline_maximum_barrier_kBT']):.6f}
  kBT**

The potential-defined aperture incorporates the summed interaction of
all cap beads and is therefore stricter than the nearest-bead geometric
estimate.

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Water-only minimization authorized:
  **{'YES' if accepted else 'NO'}**
- Short frozen-solute NVT authorized:
  **NO**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `{required_next_step}`

R2 remains a frozen neutral steric screening design. Static acceptance
does not establish water retention, exchange kinetics, chemical
realizability, or device suitability.
""",
        encoding="utf-8",
    )

    print()
    print(
        "Day023 R2 topology construction and "
        "static CAP-water scan completed."
    )

    print(
        "Grompp / TPR dump / static mdrun return codes: "
        f"{grompp.returncode}/"
        f"{dump.returncode}/"
        f"{mdrun_return_code}"
    )

    print(
        "System atoms / waters / CAPL / CAPU: "
        f"{len(system_atoms)}/"
        f"{WATERS}/"
        f"{capl_topology_atoms}/"
        f"{capu_topology_atoms}"
    )

    print(
        "CAP-OW sigma / epsilon / C12: "
        f"{sigma_nm:.12f} nm / "
        f"{epsilon_kJ_mol:.12f} kJ/mol / "
        f"{c12_selected:.12e} kJ mol^-1 nm^12"
    )

    print(
        "Minimum CAP-OW distance: "
        f"{selected_analytic['minimum_distance_nm']:.6f} nm"
    )

    print(
        "Analytic / GROMACS CAP-SOL LJ: "
        f"{analytic_cap_sol:.9f}/"
        f"{gromacs_cap_sol:.9f} kJ/mol"
    )

    print(
        "Absolute / relative energy difference: "
        f"{absolute_error:.9e}/"
        f"{relative_error:.9e}"
    )

    print(
        "Maximum water force from CAPS: "
        f"{selected_analytic['maximum_water_force_kJ_mol_nm']:.6f} "
        "kJ mol^-1 nm^-1"
    )

    print(
        "CAP-SOL Coulomb / CAP-HBNPYR LJ / "
        "CAP-HBNPYR Coulomb: "
        f"{float(gromacs_energies['Coul_CAPS_SOL']['value_kJ_mol']):.9f}/"
        f"{float(gromacs_energies['LJ_HBNPYR_CAPS']['value_kJ_mol']):.9f}/"
        f"{float(gromacs_energies['Coul_HBNPYR_CAPS']['value_kJ_mol']):.9f} "
        "kJ/mol"
    )

    print(
        "Actual lower/upper 5kBT aperture radius: "
        f"{lower_aperture:.6f}/"
        f"{upper_aperture:.6f} nm"
    )

    print(
        "Actual mean 5kBT aperture diameter: "
        f"{2.0 * float(aperture_metrics['mean_actual_aperture_radius_5kBT_nm']):.6f} nm"
    )

    print(
        "Lower/upper centerline maximum barrier: "
        f"{float(aperture_metrics['lower_centerline_maximum_barrier_kBT']):.6f}/"
        f"{float(aperture_metrics['upper_centerline_maximum_barrier_kBT']):.6f} kBT"
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
        "Water-only energy minimization authorized: "
        f"{'YES' if accepted else 'NO'}"
    )

    print(
        "Short frozen-solute NVT authorized: NO"
    )

    print(
        "Long mobile MD authorized: NO"
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

    print(
        f"Wrote: {relative(R2_TOPOLOGY)}"
    )

    print(
        f"Wrote: {relative(R2_STATIC_MDP)}"
    )

    print(
        f"Wrote: {relative(R2_INDEX)}"
    )

    print(
        f"Wrote: {relative(R2_TPR)}"
    )

    print(
        f"Wrote: {relative(STATIC_ENERGY_CSV)}"
    )

    print(
        f"Wrote: {relative(STRENGTH_SCAN_CSV)}"
    )

    print(
        f"Wrote: {relative(RADIAL_PROFILE_CSV)}"
    )

    print(
        f"Wrote: {relative(AXIAL_PROFILE_CSV)}"
    )

    print(
        f"Wrote: {relative(TOPOLOGY_AUDIT_CSV)}"
    )

    print(
        f"Wrote: {relative(SUMMARY_CSV)}"
    )

    print(
        f"Wrote: {relative(GATE_CSV)}"
    )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )

    if not accepted:
        raise RuntimeError(
            "R2 topology and static CAP-water scan "
            "requires review."
        )


if __name__ == "__main__":
    main()
