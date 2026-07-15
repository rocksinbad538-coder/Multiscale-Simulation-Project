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

BASELINE_TOP = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol/"
    "execution/08_nvt_mobile_100ps/"
    "08_nvt_mobile_100ps_processed.top"
)

DAY023_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design"
)

PROTOTYPE_ROOT = (
    DAY023_ROOT
    / "02_r1_steric_cap_prototype"
)

MODEL_ROOT = (
    DAY023_ROOT
    / "04_r1_negative_sigma_validation"
)

OUTPUT_ROOT = (
    DAY023_ROOT
    / "05_r1_full_static_scan"
)

R1_GRO = (
    PROTOTYPE_ROOT
    / "r1_t0_hydrated_with_steric_caps_geometry_only.gro"
)

PROTOTYPE_DEFINITION = (
    PROTOTYPE_ROOT
    / "r1_selected_steric_cap_definition.json"
)

PROTOTYPE_SUMMARY = (
    PROTOTYPE_ROOT
    / "r1_steric_cap_prototype_summary.csv"
)

CORRECTED_CONTRACT = (
    MODEL_ROOT
    / "r1_cap_nonbonded_model_contract_corrected.json"
)

SCAN_CSV = (
    OUTPUT_ROOT
    / "r1_full_static_cap_water_scan.csv"
)

PAIR_VALIDATION_CSV = (
    OUTPUT_ROOT
    / "r1_full_static_pair_validation.csv"
)

SELECTED_MODEL_JSON = (
    OUTPUT_ROOT
    / "r1_selected_cap_model.json"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R1_FULL_STATIC_CAP_WATER_VALIDATION_DAY023.md"
)

SELECTED_ROOT = (
    OUTPUT_ROOT
    / "selected"
)

EXPECTED_R1_ATOMS = 68314
HBN_ATOMS = 1680
PYR_ATOMS = 104
SOLUTE_ATOMS = HBN_ATOMS + PYR_ATOMS
WATER_SITES = 4

TEMPERATURE_K = 300.0
GAS_CONSTANT_KJ_MOL_K = 8.31446261815324e-3
KBT_KJ_MOL = (
    TEMPERATURE_K
    * GAS_CONSTANT_KJ_MOL_K
)

TARGET_LEVELS_KBT = (
    5.0,
    10.0,
    20.0,
    40.0,
)

PAIR_SIGMA_NM = -0.17
STATIC_CUTOFF_NM = 1.0

MIN_INITIAL_CAP_WATER_DISTANCE_NM = 0.218
MIN_HOLE_BARRIER_KBT = 100.0
MAX_ANALYTIC_GMX_RELATIVE_ERROR = 5.0e-4
MAX_ZERO_INTERACTION_ABS_KJ_MOL = 1.0e-5


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


def active_line(raw_line: str) -> str:
    return raw_line.split(
        ";",
        1,
    )[0].strip()


def section_name(raw_line: str) -> str | None:
    match = re.match(
        r"^\[\s*([^\]]+?)\s*\]$",
        active_line(raw_line),
    )

    if match is None:
        return None

    return (
        match.group(1)
        .strip()
        .lower()
    )


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
        for field in row:
            if field not in fields:
                fields.append(field)

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


def parse_baseline_topology(
    path: Path,
) -> dict[str, Any]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    section_indices: dict[
        str,
        list[int],
    ] = {}

    for index, line in enumerate(lines):
        name = section_name(line)

        if name is None:
            continue

        section_indices.setdefault(
            name,
            [],
        ).append(index)

    for required in (
        "atomtypes",
        "moleculetype",
        "system",
        "molecules",
    ):
        if required not in section_indices:
            raise RuntimeError(
                f"Baseline topology lacks [{required}]."
            )

    first_molecule_index = (
        section_indices[
            "moleculetype"
        ][0]
    )

    system_index = (
        section_indices[
            "system"
        ][0]
    )

    molecules_index = (
        section_indices[
            "molecules"
        ][0]
    )

    if not (
        first_molecule_index
        < system_index
        < molecules_index
    ):
        raise RuntimeError(
            "Unexpected topology directive order."
        )

    molecule_definitions: dict[
        str,
        dict[str, Any],
    ] = {}

    system_molecules = []

    current_section = ""
    current_molecule: str | None = None

    for raw_line in lines:
        detected = section_name(raw_line)

        if detected is not None:
            current_section = detected
            continue

        line = active_line(raw_line)

        if (
            not line
            or line.startswith("#")
        ):
            continue

        fields = line.split()

        if current_section == "moleculetype":
            if len(fields) < 2:
                continue

            current_molecule = fields[0]

            molecule_definitions[
                current_molecule
            ] = {
                "nrexcl": int(fields[1]),
                "atom_count": 0,
            }

        elif (
            current_section == "atoms"
            and current_molecule is not None
        ):
            try:
                int(fields[0])
            except (
                IndexError,
                ValueError,
            ):
                continue

            molecule_definitions[
                current_molecule
            ][
                "atom_count"
            ] += 1

        elif current_section == "molecules":
            if len(fields) < 2:
                continue

            try:
                count = int(fields[1])
            except ValueError:
                continue

            system_molecules.append(
                {
                    "name": fields[0],
                    "count": count,
                }
            )

    return {
        "lines": lines,
        "first_molecule_index": (
            first_molecule_index
        ),
        "system_index": system_index,
        "molecules_index": molecules_index,
        "molecule_definitions": (
            molecule_definitions
        ),
        "system_molecules": (
            system_molecules
        ),
    }


def insert_cap_parameters(
    prefix: list[str],
    *,
    water_oxygen_type: str,
    pair_epsilon_kj_mol: float,
) -> list[str]:
    result = list(prefix)

    sections = []

    for index, line in enumerate(result):
        name = section_name(line)

        if name is not None:
            sections.append(
                (
                    index,
                    name,
                )
            )

    atomtype_sections = [
        item
        for item in sections
        if item[1] == "atomtypes"
    ]

    if not atomtype_sections:
        raise RuntimeError(
            "No atomtypes section in topology prefix."
        )

    atomtype_start = atomtype_sections[
        -1
    ][0]

    atomtype_end = len(result)

    for index, _ in sections:
        if index > atomtype_start:
            atomtype_end = index
            break

    cap_atomtype_line = (
        "CAP     0     12.011000     0.000000     "
        "A     0.000000000000e+00     "
        "0.000000000000e+00"
    )

    result.insert(
        atomtype_end,
        cap_atomtype_line,
    )

    result.extend(
        [
            "",
            "[ nonbond_params ]",
            (
                "; i     j     funct     sigma_nm"
                "                 epsilon_kJ_mol"
            ),
            (
                f"CAP     {water_oxygen_type:<8s} "
                f"1     {PAIR_SIGMA_NM:.12e}     "
                f"{pair_epsilon_kj_mol:.12e}"
            ),
            "",
        ]
    )

    return result


def cap_molecule_block(
    molecule_name: str,
    residue_name: str,
    atom_count: int,
) -> list[str]:
    lines = [
        "",
        "[ moleculetype ]",
        "; name     nrexcl",
        f"{molecule_name:<10s} 0",
        "",
        "[ atoms ]",
        (
            "; nr   type   resnr   residue   atom"
            "    cgnr   charge   mass"
        ),
    ]

    for atom_index in range(
        1,
        atom_count + 1,
    ):
        atom_name = (
            f"C{atom_index:04d}"
        )

        lines.append(
            f"{atom_index:6d} "
            f"{'CAP':<7s} "
            f"{1:6d} "
            f"{residue_name:<7s} "
            f"{atom_name:<6s} "
            f"{atom_index:6d} "
            f"{0.0:12.6f} "
            f"{12.011:12.6f}"
        )

    lines.append("")

    return lines


def build_candidate_topology(
    parsed: dict[str, Any],
    destination: Path,
    *,
    water_molecule_name: str,
    water_oxygen_type: str,
    retained_waters: int,
    cap_beads_per_end: int,
    pair_epsilon_kj_mol: float,
) -> tuple[
    int,
    list[dict[str, Any]],
]:
    lines = parsed[
        "lines"
    ]

    first_molecule_index = parsed[
        "first_molecule_index"
    ]

    system_index = parsed[
        "system_index"
    ]

    molecules_index = parsed[
        "molecules_index"
    ]

    prefix = insert_cap_parameters(
        lines[
            :first_molecule_index
        ],
        water_oxygen_type=water_oxygen_type,
        pair_epsilon_kj_mol=(
            pair_epsilon_kj_mol
        ),
    )

    baseline_molecule_blocks = lines[
        first_molecule_index:
        system_index
    ]

    system_block = lines[
        system_index:
        molecules_index
    ]

    proposed_molecules = []
    calculated_atoms = 0

    for entry in parsed[
        "system_molecules"
    ]:
        name = entry[
            "name"
        ]

        count = (
            retained_waters
            if name == water_molecule_name
            else entry[
                "count"
            ]
        )

        if (
            name
            not in parsed[
                "molecule_definitions"
            ]
        ):
            raise RuntimeError(
                f"No definition for molecule {name}."
            )

        atoms_per_molecule = int(
            parsed[
                "molecule_definitions"
            ][
                name
            ][
                "atom_count"
            ]
        )

        calculated_atoms += (
            count
            * atoms_per_molecule
        )

        proposed_molecules.append(
            {
                "name": name,
                "count": count,
                "atoms_per_molecule": (
                    atoms_per_molecule
                ),
            }
        )

    for name in (
        "CAPL",
        "CAPU",
    ):
        proposed_molecules.append(
            {
                "name": name,
                "count": 1,
                "atoms_per_molecule": (
                    cap_beads_per_end
                ),
            }
        )

        calculated_atoms += (
            cap_beads_per_end
        )

    output_lines = []

    output_lines.extend(prefix)
    output_lines.extend(
        baseline_molecule_blocks
    )

    output_lines.extend(
        cap_molecule_block(
            "CAPL",
            "CPL",
            cap_beads_per_end,
        )
    )

    output_lines.extend(
        cap_molecule_block(
            "CAPU",
            "CPU",
            cap_beads_per_end,
        )
    )

    output_lines.extend(
        system_block
    )

    output_lines.extend(
        [
            "[ molecules ]",
            "; molecule     count",
        ]
    )

    for entry in proposed_molecules:
        output_lines.append(
            f"{entry['name']:<16s} "
            f"{entry['count']}"
        )

    output_lines.append("")

    destination.write_text(
        "\n".join(
            output_lines
        )
        + "\n",
        encoding="utf-8",
    )

    return (
        calculated_atoms,
        proposed_molecules,
    )


def write_index(
    path: Path,
    *,
    retained_waters: int,
    cap_beads_per_end: int,
) -> dict[str, tuple[int, int]]:
    water_atoms = (
        retained_waters
        * WATER_SITES
    )

    ranges = {
        "HBN": (
            1,
            HBN_ATOMS,
        ),
        "PYR": (
            HBN_ATOMS + 1,
            SOLUTE_ATOMS,
        ),
        "HBN_PYR": (
            1,
            SOLUTE_ATOMS,
        ),
        "SOL": (
            SOLUTE_ATOMS + 1,
            SOLUTE_ATOMS
            + water_atoms,
        ),
        "CAPL": (
            SOLUTE_ATOMS
            + water_atoms
            + 1,
            SOLUTE_ATOMS
            + water_atoms
            + cap_beads_per_end,
        ),
        "CAPU": (
            SOLUTE_ATOMS
            + water_atoms
            + cap_beads_per_end
            + 1,
            SOLUTE_ATOMS
            + water_atoms
            + 2
            * cap_beads_per_end,
        ),
        "CAPS": (
            SOLUTE_ATOMS
            + water_atoms
            + 1,
            SOLUTE_ATOMS
            + water_atoms
            + 2
            * cap_beads_per_end,
        ),
        "System": (
            1,
            SOLUTE_ATOMS
            + water_atoms
            + 2
            * cap_beads_per_end,
        ),
    }

    lines = []

    for name, (
        first_atom,
        last_atom,
    ) in ranges.items():
        lines.append(
            f"[ {name} ]"
        )

        indices = list(
            range(
                first_atom,
                last_atom + 1,
            )
        )

        for start in range(
            0,
            len(indices),
            15,
        ):
            lines.append(
                " ".join(
                    str(value)
                    for value in indices[
                        start:
                        start + 15
                    ]
                )
            )

        lines.append("")

    path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    return ranges


def write_static_mdp(path: Path) -> None:
    path.write_text(
        f"""integrator               = md
dt                       = 0.001
nsteps                   = 1
continuation             = no

cutoff-scheme            = Verlet
nstlist                  = 1
verlet-buffer-tolerance  = -1
rlist                    = {STATIC_CUTOFF_NM:.3f}

coulombtype              = Cut-off
coulomb-modifier         = None
rcoulomb                 = {STATIC_CUTOFF_NM:.3f}

vdwtype                  = Cut-off
vdw-modifier             = None
rvdw                     = {STATIC_CUTOFF_NM:.3f}
DispCorr                 = no

constraints              = none
pbc                      = xyz
periodic-molecules       = no

tcoupl                   = no
pcoupl                   = no
gen-vel                  = no
comm-mode                = none

nstcalcenergy            = 1
nstenergy                = 1
nstlog                   = 1
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstxout-compressed       = 0

energygrps               = CAPS SOL HBN_PYR
""",
        encoding="utf-8",
    )


def read_gro_positions(
    path: Path,
    *,
    retained_waters: int,
    cap_beads_per_end: int,
    water_oxygen_atom_name: str,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    atom_count = int(
        lines[1].strip()
    )

    if atom_count != EXPECTED_R1_ATOMS:
        raise RuntimeError(
            f"Unexpected R1 GRO atom count: "
            f"{atom_count}/{EXPECTED_R1_ATOMS}"
        )

    positions = []
    atom_names = []

    for line in lines[
        2:
        2 + atom_count
    ]:
        atom_names.append(
            line[
                10:15
            ].strip()
        )

        positions.append(
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
            ]
        )

    coordinates = np.asarray(
        positions,
        dtype=float,
    )

    box_values = [
        float(value)
        for value in lines[
            2 + atom_count
        ].split()
    ]

    if len(box_values) != 3:
        raise RuntimeError(
            "Static scan requires an orthorhombic box."
        )

    box_lengths = np.asarray(
        box_values,
        dtype=float,
    )

    water_start = SOLUTE_ATOMS
    water_stop = (
        water_start
        + retained_waters
        * WATER_SITES
    )

    water_oxygen_indices = []

    for water_index in range(
        retained_waters
    ):
        start = (
            water_start
            + water_index
            * WATER_SITES
        )

        chunk_names = atom_names[
            start:
            start + WATER_SITES
        ]

        matches = [
            local_index
            for local_index, name
            in enumerate(chunk_names)
            if name
            == water_oxygen_atom_name
        ]

        if len(matches) != 1:
            raise RuntimeError(
                "Could not identify exactly one water "
                f"oxygen in molecule {water_index}: "
                f"{chunk_names}"
            )

        water_oxygen_indices.append(
            start
            + matches[0]
        )

    cap_start = water_stop
    cap_stop = (
        cap_start
        + 2
        * cap_beads_per_end
    )

    if cap_stop != atom_count:
        raise RuntimeError(
            "Coordinate segment accounting failed."
        )

    return (
        coordinates[
            water_oxygen_indices
        ],
        coordinates[
            cap_start:
            cap_stop
        ],
        box_lengths,
    )


def analytic_cap_water_metrics(
    water_oxygen_positions: np.ndarray,
    cap_positions: np.ndarray,
    box_lengths: np.ndarray,
    *,
    c12_kj_mol_nm12: float,
) -> dict[str, float]:
    water_forces = np.zeros_like(
        water_oxygen_positions
    )

    cap_forces = np.zeros_like(
        cap_positions
    )

    total_energy = 0.0
    pair_count = 0
    minimum_distance = math.inf
    maximum_pair_energy = 0.0

    chunk_size = 512

    for start in range(
        0,
        len(water_oxygen_positions),
        chunk_size,
    ):
        stop = min(
            start + chunk_size,
            len(water_oxygen_positions),
        )

        displacement = (
            water_oxygen_positions[
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

        displacement -= (
            box_lengths
            * np.round(
                displacement
                / box_lengths
            )
        )

        distance_squared = np.sum(
            displacement
            * displacement,
            axis=2,
        )

        distances = np.sqrt(
            distance_squared
        )

        local_minimum = float(
            np.min(distances)
        )

        minimum_distance = min(
            minimum_distance,
            local_minimum,
        )

        mask = (
            distances
            < STATIC_CUTOFF_NM
        )

        pair_count += int(
            np.count_nonzero(mask)
        )

        safe_distances = np.where(
            mask,
            distances,
            1.0,
        )

        inverse_r12 = np.where(
            mask,
            safe_distances ** -12,
            0.0,
        )

        pair_energies = (
            c12_kj_mol_nm12
            * inverse_r12
        )

        total_energy += float(
            np.sum(
                pair_energies
            )
        )

        maximum_pair_energy = max(
            maximum_pair_energy,
            float(
                np.max(
                    pair_energies
                )
            ),
        )

        force_prefactor = np.where(
            mask,
            12.0
            * c12_kj_mol_nm12
            / safe_distances ** 14,
            0.0,
        )

        pair_forces = (
            force_prefactor[
                :,
                :,
                None,
            ]
            * displacement
        )

        water_forces[
            start:stop
        ] += np.sum(
            pair_forces,
            axis=1,
        )

        cap_forces -= np.sum(
            pair_forces,
            axis=0,
        )

    water_force_norms = np.linalg.norm(
        water_forces,
        axis=1,
    )

    cap_force_norms = np.linalg.norm(
        cap_forces,
        axis=1,
    )

    return {
        "pair_count_within_cutoff": (
            float(pair_count)
        ),
        "minimum_cap_water_distance_nm": (
            minimum_distance
        ),
        "analytic_total_cap_water_energy_kJ_mol": (
            total_energy
        ),
        "maximum_single_pair_energy_kJ_mol": (
            maximum_pair_energy
        ),
        "maximum_water_force_kJ_mol_nm": (
            float(
                np.max(
                    water_force_norms
                )
            )
        ),
        "rms_water_force_kJ_mol_nm": (
            float(
                np.sqrt(
                    np.mean(
                        water_force_norms
                        * water_force_norms
                    )
                )
            )
        ),
        "maximum_cap_force_kJ_mol_nm": (
            float(
                np.max(
                    cap_force_norms
                )
            )
        ),
    }


def parse_energy_menu(
    output: str,
) -> dict[str, int]:
    terms: dict[str, int] = {}

    pattern = re.compile(
        r"(?<!\S)(\d+)\s+"
        r"([A-Za-z][A-Za-z0-9_.:+\-]*)"
    )

    for line in output.splitlines():
        for match in pattern.finditer(line):
            terms[
                match.group(2)
            ] = int(
                match.group(1)
            )

    return terms


def energy_menu(
    gmx: str,
    edr: Path,
    cwd: Path,
) -> dict[str, int]:
    probe = run_command(
        [
            gmx,
            "energy",
            "-f",
            str(edr),
            "-o",
            str(
                cwd
                / "energy_menu_probe.xvg"
            ),
        ],
        cwd=cwd,
        input_text="0\n",
    )

    menu_path = (
        cwd
        / "energy_menu.txt"
    )

    menu_path.write_text(
        probe.stdout,
        encoding="utf-8",
    )

    terms = parse_energy_menu(
        probe.stdout
    )

    if not terms:
        raise RuntimeError(
            f"Could not parse energy menu in {cwd}"
        )

    return terms


def select_cross_term(
    terms: dict[str, int],
    prefix: str,
    group_a: str,
    group_b: str,
) -> tuple[
    str | None,
    int | None,
]:
    exact_candidates = (
        f"{prefix}:{group_a}-{group_b}",
        f"{prefix}:{group_b}-{group_a}",
    )

    for name in exact_candidates:
        if name in terms:
            return (
                name,
                terms[name],
            )

    for name, number in terms.items():
        if not name.startswith(
            f"{prefix}:"
        ):
            continue

        suffix = name.split(
            ":",
            1,
        )[1]

        if (
            group_a in suffix
            and group_b in suffix
        ):
            return (
                name,
                number,
            )

    return (
        None,
        None,
    )


def extract_energy_value(
    gmx: str,
    edr: Path,
    cwd: Path,
    term_number: int,
    output_name: str,
) -> float:
    xvg = (
        cwd
        / output_name
    )

    completed = run_command(
        [
            gmx,
            "energy",
            "-f",
            str(edr),
            "-o",
            str(xvg),
            "-xvg",
            "none",
            "-dp",
        ],
        cwd=cwd,
        input_text=(
            f"{term_number}\n"
            "0\n"
        ),
    )

    if (
        completed.returncode != 0
        or not xvg.exists()
    ):
        raise RuntimeError(
            f"Could not extract energy term "
            f"{term_number} in {cwd}"
        )

    numeric_rows = []

    for raw_line in xvg.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        line = raw_line.strip()

        if (
            not line
            or line.startswith("#")
            or line.startswith("@")
        ):
            continue

        fields = line.split()

        if len(fields) >= 2:
            numeric_rows.append(
                (
                    float(fields[0]),
                    float(fields[1]),
                )
            )

    if not numeric_rows:
        raise RuntimeError(
            f"No numeric energy values in {xvg}"
        )

    return float(
        numeric_rows[
            -1
        ][1]
    )


def optional_cross_energy(
    gmx: str,
    edr: Path,
    cwd: Path,
    terms: dict[str, int],
    prefix: str,
    group_a: str,
    group_b: str,
    output_name: str,
) -> tuple[
    str,
    float,
]:
    name, number = select_cross_term(
        terms,
        prefix,
        group_a,
        group_b,
    )

    if (
        name is None
        or number is None
    ):
        return (
            "ABSENT",
            0.0,
        )

    value = extract_energy_value(
        gmx,
        edr,
        cwd,
        number,
        output_name,
    )

    return (
        name,
        value,
    )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        BASELINE_TOP,
        R1_GRO,
        PROTOTYPE_DEFINITION,
        PROTOTYPE_SUMMARY,
        CORRECTED_CONTRACT,
    ):
        require_file(required)

    prototype = json.loads(
        PROTOTYPE_DEFINITION.read_text(
            encoding="utf-8"
        )
    )

    contract = json.loads(
        CORRECTED_CONTRACT.read_text(
            encoding="utf-8"
        )
    )

    prototype_summary = (
        read_single_csv_row(
            PROTOTYPE_SUMMARY
        )
    )

    if (
        contract.get(
            "decision"
        )
        !=
        "PURE_R12_NEGATIVE_SIGMA_OVERRIDE_VALIDATED"
    ):
        raise RuntimeError(
            "Corrected nonbonded contract is not "
            "authorized for full-system validation."
        )

    if not bool(
        contract.get(
            "all_microvalidation_candidates_pass"
        )
    ):
        raise RuntimeError(
            "The CAP-OW microvalidation did not pass."
        )

    retained_waters = int(
        prototype[
            "retained_water_molecules"
        ]
    )

    cap_beads_per_end = int(
        prototype[
            "cap_beads_per_end"
        ]
    )

    water_molecule_name = str(
        contract[
            "water_molecule_name"
        ]
    )

    water_oxygen_type = str(
        contract[
            "water_oxygen_atom_type"
        ]
    )

    water_oxygen_atom_name = str(
        contract[
            "water_oxygen_atom_name"
        ]
    )

    coverage_hole_nm = float(
        prototype_summary[
            "coverage_hole_nm"
        ]
    )

    pruning_cutoff_nm = float(
        prototype[
            "water_pruning_cutoff_nm"
        ]
    )

    parsed_topology = (
        parse_baseline_topology(
            BASELINE_TOP
        )
    )

    gmx = locate_gmx()

    water_oxygen_positions, cap_positions, box_lengths = (
        read_gro_positions(
            R1_GRO,
            retained_waters=(
                retained_waters
            ),
            cap_beads_per_end=(
                cap_beads_per_end
            ),
            water_oxygen_atom_name=(
                water_oxygen_atom_name
            ),
        )
    )

    scan_rows = []
    pair_rows = []
    candidate_objects = []

    index_path = (
        OUTPUT_ROOT
        / "r1_static_groups.ndx"
    )

    group_ranges = write_index(
        index_path,
        retained_waters=retained_waters,
        cap_beads_per_end=(
            cap_beads_per_end
        ),
    )

    for target_level_kbt in (
        TARGET_LEVELS_KBT
    ):
        candidate_name = (
            f"target_{target_level_kbt:.0f}kBT"
        )

        candidate_root = (
            OUTPUT_ROOT
            / candidate_name
        )

        candidate_root.mkdir(
            parents=True,
            exist_ok=True,
        )

        pair_epsilon_kj_mol = (
            target_level_kbt
            * KBT_KJ_MOL
            / 4.0
        )

        c12_kj_mol_nm12 = (
            4.0
            * pair_epsilon_kj_mol
            * abs(
                PAIR_SIGMA_NM
            ) ** 12
        )

        topology_path = (
            candidate_root
            / "r1_full_system.top"
        )

        calculated_atoms, proposed_molecules = (
            build_candidate_topology(
                parsed_topology,
                topology_path,
                water_molecule_name=(
                    water_molecule_name
                ),
                water_oxygen_type=(
                    water_oxygen_type
                ),
                retained_waters=(
                    retained_waters
                ),
                cap_beads_per_end=(
                    cap_beads_per_end
                ),
                pair_epsilon_kj_mol=(
                    pair_epsilon_kj_mol
                ),
            )
        )

        if (
            calculated_atoms
            != EXPECTED_R1_ATOMS
        ):
            raise RuntimeError(
                "Generated topology atom count is "
                f"{calculated_atoms}, expected "
                f"{EXPECTED_R1_ATOMS}."
            )

        mdp_path = (
            candidate_root
            / "static_rerun.mdp"
        )

        write_static_mdp(
            mdp_path
        )

        tpr_path = (
            candidate_root
            / "r1_static.tpr"
        )

        processed_top_path = (
            candidate_root
            / "r1_static_processed.top"
        )

        mdout_path = (
            candidate_root
            / "mdout.mdp"
        )

        grompp = run_command(
            [
                gmx,
                "grompp",
                "-f",
                str(mdp_path),
                "-c",
                str(R1_GRO),
                "-p",
                str(topology_path),
                "-n",
                str(index_path),
                "-o",
                str(tpr_path),
                "-po",
                str(mdout_path),
                "-pp",
                str(
                    processed_top_path
                ),
                "-maxwarn",
                "0",
            ],
            cwd=candidate_root,
        )

        (
            candidate_root
            / "grompp.log"
        ).write_text(
            grompp.stdout,
            encoding="utf-8",
        )

        if (
            grompp.returncode != 0
            or not tpr_path.exists()
        ):
            raise RuntimeError(
                f"grompp failed for {candidate_name}. "
                f"See {candidate_root / 'grompp.log'}"
            )

        rerun = run_command(
            [
                gmx,
                "mdrun",
                "-s",
                str(tpr_path),
                "-rerun",
                str(R1_GRO),
                "-deffnm",
                "static",
                "-ntmpi",
                "1",
                "-ntomp",
                "1",
            ],
            cwd=candidate_root,
        )

        (
            candidate_root
            / "static_rerun.log"
        ).write_text(
            rerun.stdout,
            encoding="utf-8",
        )

        edr_path = (
            candidate_root
            / "static.edr"
        )

        if (
            rerun.returncode != 0
            or not edr_path.exists()
        ):
            raise RuntimeError(
                f"Static rerun failed for "
                f"{candidate_name}."
            )

        terms = energy_menu(
            gmx,
            edr_path,
            candidate_root,
        )

        cap_sol_lj_name, cap_sol_lj_number = (
            select_cross_term(
                terms,
                "LJ-SR",
                "CAPS",
                "SOL",
            )
        )

        if (
            cap_sol_lj_name is None
            or cap_sol_lj_number is None
        ):
            raise RuntimeError(
                "Could not resolve the CAPS-SOL "
                f"LJ-SR term for {candidate_name}."
            )

        measured_cap_water_energy = (
            extract_energy_value(
                gmx,
                edr_path,
                candidate_root,
                cap_sol_lj_number,
                "cap_sol_lj.xvg",
            )
        )

        (
            cap_sol_coul_name,
            cap_sol_coul_energy,
        ) = optional_cross_energy(
            gmx,
            edr_path,
            candidate_root,
            terms,
            "Coul-SR",
            "CAPS",
            "SOL",
            "cap_sol_coul.xvg",
        )

        (
            cap_solute_lj_name,
            cap_solute_lj_energy,
        ) = optional_cross_energy(
            gmx,
            edr_path,
            candidate_root,
            terms,
            "LJ-SR",
            "CAPS",
            "HBN_PYR",
            "cap_solute_lj.xvg",
        )

        (
            cap_solute_coul_name,
            cap_solute_coul_energy,
        ) = optional_cross_energy(
            gmx,
            edr_path,
            candidate_root,
            terms,
            "Coul-SR",
            "CAPS",
            "HBN_PYR",
            "cap_solute_coul.xvg",
        )

        analytic = (
            analytic_cap_water_metrics(
                water_oxygen_positions,
                cap_positions,
                box_lengths,
                c12_kj_mol_nm12=(
                    c12_kj_mol_nm12
                ),
            )
        )

        analytic_energy = float(
            analytic[
                "analytic_total_cap_water_energy_kJ_mol"
            ]
        )

        absolute_error = abs(
            measured_cap_water_energy
            - analytic_energy
        )

        relative_error = (
            absolute_error
            / abs(
                analytic_energy
            )
            if abs(
                analytic_energy
            )
            > 1.0e-12
            else absolute_error
        )

        energy_at_pruning_cutoff_kbt = (
            target_level_kbt
            * (
                abs(
                    PAIR_SIGMA_NM
                )
                / pruning_cutoff_nm
            ) ** 12
        )

        barrier_at_coverage_hole_kbt = (
            target_level_kbt
            * (
                abs(
                    PAIR_SIGMA_NM
                )
                / coverage_hole_nm
            ) ** 12
        )

        gates = {
            "grompp": (
                grompp.returncode == 0
            ),
            "static_rerun": (
                rerun.returncode == 0
            ),
            "topology_atom_count": (
                calculated_atoms
                == EXPECTED_R1_ATOMS
            ),
            "minimum_initial_cap_water_distance": (
                float(
                    analytic[
                        "minimum_cap_water_distance_nm"
                    ]
                )
                >=
                MIN_INITIAL_CAP_WATER_DISTANCE_NM
            ),
            "analytic_GROMACS_energy_agreement": (
                relative_error
                <=
                MAX_ANALYTIC_GMX_RELATIVE_ERROR
            ),
            "zero_cap_solute_LJ": (
                abs(
                    cap_solute_lj_energy
                )
                <=
                MAX_ZERO_INTERACTION_ABS_KJ_MOL
            ),
            "zero_cap_solute_Coulomb": (
                abs(
                    cap_solute_coul_energy
                )
                <=
                MAX_ZERO_INTERACTION_ABS_KJ_MOL
            ),
            "zero_cap_water_Coulomb": (
                abs(
                    cap_sol_coul_energy
                )
                <=
                MAX_ZERO_INTERACTION_ABS_KJ_MOL
            ),
            "finite_static_energy": (
                math.isfinite(
                    measured_cap_water_energy
                )
            ),
            "finite_static_forces": (
                math.isfinite(
                    float(
                        analytic[
                            "maximum_water_force_kJ_mol_nm"
                        ]
                    )
                )
                and math.isfinite(
                    float(
                        analytic[
                            "maximum_cap_force_kJ_mol_nm"
                        ]
                    )
                )
            ),
            "coverage_hole_barrier": (
                barrier_at_coverage_hole_kbt
                >=
                MIN_HOLE_BARRIER_KBT
            ),
        }

        failed_gates = [
            name
            for name, passed
            in gates.items()
            if not passed
        ]

        candidate_pass = (
            len(
                failed_gates
            )
            == 0
        )

        row = {
            "candidate": candidate_name,
            "target_energy_kBT_at_0p17nm": (
                target_level_kbt
            ),
            "pair_sigma_nm": (
                PAIR_SIGMA_NM
            ),
            "pair_epsilon_kJ_mol": (
                pair_epsilon_kj_mol
            ),
            "C6_kJ_mol_nm6": 0.0,
            "C12_kJ_mol_nm12": (
                c12_kj_mol_nm12
            ),
            "topology_atoms": (
                calculated_atoms
            ),
            "grompp_return_code": (
                grompp.returncode
            ),
            "rerun_return_code": (
                rerun.returncode
            ),
            "cap_SOL_LJ_term": (
                cap_sol_lj_name
            ),
            "GROMACS_cap_water_energy_kJ_mol": (
                measured_cap_water_energy
            ),
            "analytic_cap_water_energy_kJ_mol": (
                analytic_energy
            ),
            "energy_absolute_error_kJ_mol": (
                absolute_error
            ),
            "energy_relative_error": (
                relative_error
            ),
            "cap_SOL_Coulomb_term": (
                cap_sol_coul_name
            ),
            "cap_SOL_Coulomb_kJ_mol": (
                cap_sol_coul_energy
            ),
            "cap_HBN_PYR_LJ_term": (
                cap_solute_lj_name
            ),
            "cap_HBN_PYR_LJ_kJ_mol": (
                cap_solute_lj_energy
            ),
            "cap_HBN_PYR_Coulomb_term": (
                cap_solute_coul_name
            ),
            "cap_HBN_PYR_Coulomb_kJ_mol": (
                cap_solute_coul_energy
            ),
            "minimum_cap_water_distance_nm": (
                analytic[
                    "minimum_cap_water_distance_nm"
                ]
            ),
            "pair_count_within_cutoff": (
                int(
                    analytic[
                        "pair_count_within_cutoff"
                    ]
                )
            ),
            "maximum_single_pair_energy_kJ_mol": (
                analytic[
                    "maximum_single_pair_energy_kJ_mol"
                ]
            ),
            "maximum_water_force_kJ_mol_nm": (
                analytic[
                    "maximum_water_force_kJ_mol_nm"
                ]
            ),
            "rms_water_force_kJ_mol_nm": (
                analytic[
                    "rms_water_force_kJ_mol_nm"
                ]
            ),
            "maximum_cap_force_kJ_mol_nm": (
                analytic[
                    "maximum_cap_force_kJ_mol_nm"
                ]
            ),
            "energy_at_pruning_cutoff_kBT": (
                energy_at_pruning_cutoff_kbt
            ),
            "barrier_at_coverage_hole_kBT": (
                barrier_at_coverage_hole_kbt
            ),
            "candidate_pass": (
                candidate_pass
            ),
            "failed_gates": (
                " | ".join(
                    failed_gates
                )
            ),
        }

        scan_rows.append(row)

        for gate_name, passed in gates.items():
            pair_rows.append(
                {
                    "candidate": (
                        candidate_name
                    ),
                    "gate": gate_name,
                    "pass": passed,
                }
            )

        candidate_objects.append(
            {
                "row": row,
                "root": candidate_root,
                "topology": topology_path,
                "tpr": tpr_path,
                "mdp": mdp_path,
                "processed_topology": (
                    processed_top_path
                ),
                "proposed_molecules": (
                    proposed_molecules
                ),
            }
        )

    write_csv(
        SCAN_CSV,
        scan_rows,
    )

    write_csv(
        PAIR_VALIDATION_CSV,
        pair_rows,
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
        raise RuntimeError(
            "No full-system R1 cap candidate "
            "passed the static gates."
        )

    # Select the weakest validated repulsion. This minimizes
    # perturbation of the initial solvent while retaining a
    # large steric barrier at the largest planar coverage hole.
    selected = sorted(
        passing,
        key=lambda candidate: float(
            candidate[
                "row"
            ][
                "target_energy_kBT_at_0p17nm"
            ]
        ),
    )[0]

    selected_row = selected[
        "row"
    ]

    SELECTED_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    selected_topology = (
        SELECTED_ROOT
        / "r1_selected_cap_model.top"
    )

    selected_tpr = (
        SELECTED_ROOT
        / "r1_selected_static_validation.tpr"
    )

    selected_mdp = (
        SELECTED_ROOT
        / "r1_selected_static_validation.mdp"
    )

    selected_processed_top = (
        SELECTED_ROOT
        / "r1_selected_static_processed.top"
    )

    selected_index = (
        SELECTED_ROOT
        / "r1_selected_groups.ndx"
    )

    shutil.copy2(
        selected[
            "topology"
        ],
        selected_topology,
    )

    shutil.copy2(
        selected[
            "tpr"
        ],
        selected_tpr,
    )

    shutil.copy2(
        selected[
            "mdp"
        ],
        selected_mdp,
    )

    shutil.copy2(
        selected[
            "processed_topology"
        ],
        selected_processed_top,
    )

    shutil.copy2(
        index_path,
        selected_index,
    )

    selected_model = {
        "decision": (
            "FULL_R1_STATIC_CAP_WATER_MODEL_VALIDATED"
        ),
        "selection_principle": (
            "WEAKEST_VALIDATED_REPULSION_WITH_"
            "COVERAGE_HOLE_BARRIER_AT_LEAST_100_KBT"
        ),
        "selected_candidate": (
            selected_row[
                "candidate"
            ]
        ),
        "target_energy_kBT_at_0p17nm": (
            selected_row[
                "target_energy_kBT_at_0p17nm"
            ]
        ),
        "pair_sigma_nm": (
            selected_row[
                "pair_sigma_nm"
            ]
        ),
        "pair_epsilon_kJ_mol": (
            selected_row[
                "pair_epsilon_kJ_mol"
            ]
        ),
        "C6_kJ_mol_nm6": 0.0,
        "C12_kJ_mol_nm12": (
            selected_row[
                "C12_kJ_mol_nm12"
            ]
        ),
        "minimum_cap_water_distance_nm": (
            selected_row[
                "minimum_cap_water_distance_nm"
            ]
        ),
        "initial_cap_water_energy_kJ_mol": (
            selected_row[
                "GROMACS_cap_water_energy_kJ_mol"
            ]
        ),
        "maximum_water_force_kJ_mol_nm": (
            selected_row[
                "maximum_water_force_kJ_mol_nm"
            ]
        ),
        "energy_at_pruning_cutoff_kBT": (
            selected_row[
                "energy_at_pruning_cutoff_kBT"
            ]
        ),
        "barrier_at_coverage_hole_kBT": (
            selected_row[
                "barrier_at_coverage_hole_kBT"
            ]
        ),
        "retained_water_molecules": (
            retained_waters
        ),
        "cap_beads_per_end": (
            cap_beads_per_end
        ),
        "R1_atom_count": (
            EXPECTED_R1_ATOMS
        ),
        "R1_geometry": relative(
            R1_GRO
        ),
        "selected_topology": relative(
            selected_topology
        ),
        "selected_static_TPR": relative(
            selected_tpr
        ),
        "selected_index": relative(
            selected_index
        ),
        "caps_must_remain_frozen": True,
        "screening_solute_must_remain_frozen": True,
        "energy_minimization_authorized": True,
        "MD_execution_authorized": False,
        "required_next_step": (
            "PREPARE_AND_RUN_R1_FROZEN_SOLUTE_"
            "ENERGY_MINIMIZATION"
        ),
    }

    SELECTED_MODEL_JSON.write_text(
        json.dumps(
            selected_model,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    candidate_lines = "\n".join(
        (
            f"- `{row['candidate']}`: "
            f"CAP–water energy="
            f"{row['GROMACS_cap_water_energy_kJ_mol']:.6f} "
            f"kJ/mol; min distance="
            f"{row['minimum_cap_water_distance_nm']:.6f} nm; "
            f"maximum water force="
            f"{row['maximum_water_force_kJ_mol_nm']:.3f} "
            f"kJ mol^-1 nm^-1; hole barrier="
            f"{row['barrier_at_coverage_hole_kBT']:.1f} kBT; "
            f"**{'PASS' if row['candidate_pass'] else 'FAIL'}**"
        )
        for row in scan_rows
    )

    REPORT_MD.write_text(
        f"""# R1 Full-System Static Cap–Water Validation

## Scope

Four complete R1 topologies were constructed from the validated
68,314-atom geometry. Each topology contains:

- the accepted HBN and PYR definitions;
- {retained_waters} TIP4P/2005 waters;
- one 163-bead lower cap;
- one 163-bead upper cap;
- zero cap charge;
- zero CAP–CAP, CAP–HBN, CAP–PYR, CAP–water-H, and CAP–water-M
  interaction;
- one explicit pure-r12 CAP–OW interaction.

No energy minimization or molecular dynamics was performed.

## Static validation

The full system was evaluated using `gmx grompp` and
`gmx mdrun -rerun`. The CAP–water energy was independently calculated
as the explicit sum of all CAP–OW C12/r^12 contributions inside the
1.0 nm static cutoff.

{candidate_lines}

## Selected model

- Candidate:
  **{selected_row['candidate']}**
- Selection rule:
  **weakest validated repulsion satisfying a minimum 100 kBT barrier
  at the largest coverage hole**
- Pair sigma:
  **{selected_row['pair_sigma_nm']:.6f} nm**
- Pair epsilon:
  **{selected_row['pair_epsilon_kJ_mol']:.10f} kJ/mol**
- C12:
  **{selected_row['C12_kJ_mol_nm12']:.12e}
  kJ mol^-1 nm^12**
- Initial minimum CAP–OW distance:
  **{selected_row['minimum_cap_water_distance_nm']:.6f} nm**
- Initial CAP–water energy:
  **{selected_row['GROMACS_cap_water_energy_kJ_mol']:.6f} kJ/mol**
- Maximum initial water force from the cap:
  **{selected_row['maximum_water_force_kJ_mol_nm']:.3f}
  kJ mol^-1 nm^-1**
- Energy at the 0.22 nm pruning boundary:
  **{selected_row['energy_at_pruning_cutoff_kBT']:.6f} kBT**
- Barrier at the largest 0.114878 nm coverage hole:
  **{selected_row['barrier_at_coverage_hole_kBT']:.3f} kBT**
- GROMACS/analytic relative energy error:
  **{selected_row['energy_relative_error']:.6e}**

## Interaction checks

- CAP–HBN/PYR LJ energy:
  **{selected_row['cap_HBN_PYR_LJ_kJ_mol']:.12f} kJ/mol**
- CAP–HBN/PYR Coulomb energy:
  **{selected_row['cap_HBN_PYR_Coulomb_kJ_mol']:.12f} kJ/mol**
- CAP–water Coulomb energy:
  **{selected_row['cap_SOL_Coulomb_kJ_mol']:.12f} kJ/mol**

## Decision

- Full R1 topology validated: **YES**
- Selected cap model:
  **{selected_row['candidate']}**
- Cap atoms must remain frozen: **YES**
- HBN and PYR must remain frozen during initial screening: **YES**
- Energy minimization authorized: **YES**
- Molecular dynamics authorized: **NO**
- Required next step:
  `PREPARE_AND_RUN_R1_FROZEN_SOLUTE_ENERGY_MINIMIZATION`

Selected files:

- `{relative(selected_topology)}`
- `{relative(selected_tpr)}`
- `{relative(selected_index)}`
- `{relative(SELECTED_MODEL_JSON)}`
""",
        encoding="utf-8",
    )

    print(
        "Day023 R1 full-system static cap-water "
        "validation completed."
    )

    print(
        "Candidates evaluated / passing: "
        f"{len(scan_rows)}/"
        f"{len(passing)}"
    )

    for row in scan_rows:
        print(
            "Candidate "
            f"{row['target_energy_kBT_at_0p17nm']:.0f} kBT: "
            "GROMACS/analytic energy="
            f"{row['GROMACS_cap_water_energy_kJ_mol']:.6f}/"
            f"{row['analytic_cap_water_energy_kJ_mol']:.6f} "
            "kJ/mol; min distance="
            f"{row['minimum_cap_water_distance_nm']:.6f} nm; "
            "max water force="
            f"{row['maximum_water_force_kJ_mol_nm']:.3f}; "
            "hole barrier="
            f"{row['barrier_at_coverage_hole_kBT']:.1f} kBT; "
            f"{'PASS' if row['candidate_pass'] else 'FAIL'}"
        )

    print(
        "Selected candidate: "
        f"{selected_row['candidate']}"
    )

    print(
        "Selected sigma / epsilon / C12: "
        f"{selected_row['pair_sigma_nm']:.6f} nm / "
        f"{selected_row['pair_epsilon_kJ_mol']:.10f} kJ/mol / "
        f"{selected_row['C12_kJ_mol_nm12']:.12e}"
    )

    print(
        "Selected initial cap-water energy: "
        f"{selected_row['GROMACS_cap_water_energy_kJ_mol']:.6f} "
        "kJ/mol"
    )

    print(
        "Selected maximum initial water force: "
        f"{selected_row['maximum_water_force_kJ_mol_nm']:.3f} "
        "kJ mol^-1 nm^-1"
    )

    print(
        "Selected coverage-hole barrier: "
        f"{selected_row['barrier_at_coverage_hole_kBT']:.3f} kBT"
    )

    print(
        "CAP-HBN/PYR LJ / Coulomb: "
        f"{selected_row['cap_HBN_PYR_LJ_kJ_mol']:.12f}/"
        f"{selected_row['cap_HBN_PYR_Coulomb_kJ_mol']:.12f} "
        "kJ/mol"
    )

    print(
        "CAP-water Coulomb: "
        f"{selected_row['cap_SOL_Coulomb_kJ_mol']:.12f} "
        "kJ/mol"
    )

    print(
        "Decision: "
        "FULL_R1_STATIC_CAP_WATER_MODEL_VALIDATED"
    )

    print(
        "Energy minimization authorized: YES"
    )

    print(
        "MD execution authorized: NO"
    )

    print(
        "Required next step: "
        "PREPARE_AND_RUN_R1_FROZEN_SOLUTE_ENERGY_MINIMIZATION"
    )

    print(
        f"Wrote: {relative(SCAN_CSV)}"
    )

    print(
        f"Wrote: {relative(PAIR_VALIDATION_CSV)}"
    )

    print(
        f"Wrote: {relative(SELECTED_MODEL_JSON)}"
    )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )

    print(
        f"Wrote: {relative(selected_topology)}"
    )

    print(
        f"Wrote: {relative(selected_tpr)}"
    )

    print(
        f"Wrote: {relative(selected_index)}"
    )


if __name__ == "__main__":
    main()
