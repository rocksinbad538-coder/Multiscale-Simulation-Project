#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

PROTOCOL_ROOT = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

PROCESSED_TOP = (
    PROTOCOL_ROOT
    / "execution/08_nvt_mobile_100ps/"
    "08_nvt_mobile_100ps_processed.top"
)

R1_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design/"
    "02_r1_steric_cap_prototype"
)

R1_GRO = (
    R1_ROOT
    / "r1_t0_hydrated_with_steric_caps_geometry_only.gro"
)

R1_DEFINITION = (
    R1_ROOT
    / "r1_selected_steric_cap_definition.json"
)

OUTPUT_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design/"
    "03_r1_topology_model"
)

TOPOLOGY_INVENTORY_CSV = (
    OUTPUT_ROOT
    / "r1_baseline_topology_inventory.csv"
)

ATOMTYPE_INVENTORY_CSV = (
    OUTPUT_ROOT
    / "r1_active_atomtype_inventory.csv"
)

PROPOSED_MOLECULES_CSV = (
    OUTPUT_ROOT
    / "r1_proposed_molecule_counts.csv"
)

REPULSION_CALIBRATION_CSV = (
    OUTPUT_ROOT
    / "r1_cap_water_repulsion_calibration.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r1_cap_topology_requirement_summary.csv"
)

CONTRACT_JSON = (
    OUTPUT_ROOT
    / "r1_cap_nonbonded_model_contract.json"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R1_CAP_TOPOLOGY_REQUIREMENTS_DAY023.md"
)

HBN_ATOMS = 1680
PYR_ATOMS = 104
INHERITED_SOLUTE_ATOMS = HBN_ATOMS + PYR_ATOMS
WATER_SITES = 4

EXPECTED_R0_ATOMS = 68320
EXPECTED_R1_ATOMS = 68314

TEMPERATURE_K = 300.0
GAS_CONSTANT_KJ_MOL_K = 8.31446261815324e-3

TARGET_REPULSION_RADIUS_NM = 0.17

ENERGY_LEVELS_KBT = (
    5.0,
    10.0,
    20.0,
    40.0,
)

DIAGNOSTIC_RADII_NM = (
    0.114878,
    0.15,
    0.17,
    0.20,
    0.22,
    0.25,
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
            f"Missing or empty file: {path}"
        )


def active_line(raw_line: str) -> str:
    return raw_line.split(
        ";",
        1,
    )[0].strip()


def parse_float(value: str) -> float:
    parsed = float(value)

    if not math.isfinite(parsed):
        raise ValueError(
            f"Non-finite numeric value: {value}"
        )

    return parsed


def parse_processed_topology(
    path: Path,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "defaults": None,
        "atomtypes": {},
        "nonbond_params": [],
        "moleculetypes": {},
        "molecules": [],
    }

    current_section = ""
    current_molecule: str | None = None

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        line = active_line(
            raw_line
        )

        if not line:
            continue

        if line.startswith("#"):
            continue

        section_match = re.match(
            r"^\[\s*([^\]]+?)\s*\]$",
            line,
        )

        if section_match is not None:
            current_section = (
                section_match.group(1)
                .strip()
                .lower()
            )

            continue

        fields = line.split()

        if current_section == "defaults":
            if result["defaults"] is None:
                if len(fields) < 2:
                    raise RuntimeError(
                        "Malformed [ defaults ] entry."
                    )

                result["defaults"] = {
                    "nbfunc": int(fields[0]),
                    "comb_rule": int(fields[1]),
                    "gen_pairs": (
                        fields[2]
                        if len(fields) > 2
                        else ""
                    ),
                    "fudge_lj": (
                        fields[3]
                        if len(fields) > 3
                        else ""
                    ),
                    "fudge_qq": (
                        fields[4]
                        if len(fields) > 4
                        else ""
                    ),
                    "raw": line,
                }

            continue

        if current_section == "atomtypes":
            if len(fields) < 6:
                continue

            try:
                parameter_v = parse_float(
                    fields[-2]
                )

                parameter_w = parse_float(
                    fields[-1]
                )
            except ValueError:
                continue

            atomtype_name = fields[0]

            result["atomtypes"][
                atomtype_name
            ] = {
                "name": atomtype_name,
                "raw": line,
                "tokens": fields,
                "ptype": fields[-3],
                "parameter_v": parameter_v,
                "parameter_w": parameter_w,
            }

            continue

        if current_section == "nonbond_params":
            if len(fields) < 5:
                continue

            try:
                parameter_v = parse_float(
                    fields[3]
                )

                parameter_w = parse_float(
                    fields[4]
                )
            except ValueError:
                continue

            result["nonbond_params"].append(
                {
                    "type_i": fields[0],
                    "type_j": fields[1],
                    "function": int(fields[2]),
                    "parameter_v": parameter_v,
                    "parameter_w": parameter_w,
                    "raw": line,
                }
            )

            continue

        if current_section == "moleculetype":
            if len(fields) < 2:
                continue

            current_molecule = fields[0]

            result["moleculetypes"][
                current_molecule
            ] = {
                "name": current_molecule,
                "nrexcl": int(fields[1]),
                "atoms": [],
            }

            continue

        if current_section == "atoms":
            if current_molecule is None:
                continue

            if len(fields) < 7:
                continue

            try:
                atom_number = int(
                    fields[0]
                )

                residue_number = int(
                    fields[2]
                )

                charge = parse_float(
                    fields[6]
                )

                mass = (
                    parse_float(
                        fields[7]
                    )
                    if len(fields) > 7
                    else math.nan
                )
            except ValueError:
                continue

            result["moleculetypes"][
                current_molecule
            ][
                "atoms"
            ].append(
                {
                    "number": atom_number,
                    "type": fields[1],
                    "residue_number": residue_number,
                    "residue_name": fields[3],
                    "atom_name": fields[4],
                    "charge_group": fields[5],
                    "charge": charge,
                    "mass": mass,
                    "raw": line,
                }
            )

            continue

        if current_section == "molecules":
            if len(fields) < 2:
                continue

            try:
                molecule_count = int(
                    fields[1]
                )
            except ValueError:
                continue

            result["molecules"].append(
                {
                    "name": fields[0],
                    "count": molecule_count,
                }
            )

    if result["defaults"] is None:
        raise RuntimeError(
            "No [ defaults ] entry was parsed."
        )

    if not result["atomtypes"]:
        raise RuntimeError(
            "No [ atomtypes ] entries were parsed."
        )

    if not result["moleculetypes"]:
        raise RuntimeError(
            "No molecule definitions were parsed."
        )

    if not result["molecules"]:
        raise RuntimeError(
            "No [ molecules ] entries were parsed."
        )

    return result


def read_gro_segments(
    path: Path,
    retained_waters: int,
    cap_beads_per_end: int,
) -> dict[str, Any]:
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

    atom_lines = lines[
        2:
        2 + atom_count
    ]

    if len(atom_lines) != atom_count:
        raise RuntimeError(
            "Incomplete R1 GRO atom list."
        )

    records = []

    for index, line in enumerate(
        atom_lines,
        start=1,
    ):
        if len(line) < 44:
            raise RuntimeError(
                f"Malformed R1 GRO line: {index}"
            )

        records.append(
            {
                "global_index": index,
                "resid": int(
                    line[0:5]
                ),
                "resname": line[
                    5:10
                ].strip(),
                "atomname": line[
                    10:15
                ].strip(),
            }
        )

    water_atom_count = (
        retained_waters
        * WATER_SITES
    )

    inherited_end = (
        INHERITED_SOLUTE_ATOMS
    )

    water_end = (
        inherited_end
        + water_atom_count
    )

    lower_cap_end = (
        water_end
        + cap_beads_per_end
    )

    upper_cap_end = (
        lower_cap_end
        + cap_beads_per_end
    )

    if upper_cap_end != atom_count:
        raise RuntimeError(
            "R1 segment accounting does not match "
            f"the GRO atom count: {upper_cap_end}/"
            f"{atom_count}"
        )

    segments = {
        "inherited_HBN": records[
            :HBN_ATOMS
        ],
        "inherited_PYR": records[
            HBN_ATOMS:
            INHERITED_SOLUTE_ATOMS
        ],
        "retained_SOL": records[
            inherited_end:
            water_end
        ],
        "lower_cap": records[
            water_end:
            lower_cap_end
        ],
        "upper_cap": records[
            lower_cap_end:
            upper_cap_end
        ],
    }

    segment_summary = {}

    for name, segment in segments.items():
        segment_summary[name] = {
            "atom_count": len(
                segment
            ),
            "resnames": dict(
                Counter(
                    atom[
                        "resname"
                    ]
                    for atom in segment
                )
            ),
            "atomnames_first": [
                atom[
                    "atomname"
                ]
                for atom in segment[
                    :5
                ]
            ],
            "first_global_index": (
                segment[
                    0
                ][
                    "global_index"
                ]
                if segment
                else None
            ),
            "last_global_index": (
                segment[
                    -1
                ][
                    "global_index"
                ]
                if segment
                else None
            ),
        }

    water_chunk_consistency = True

    water_segment = segments[
        "retained_SOL"
    ]

    for start in range(
        0,
        len(
            water_segment
        ),
        WATER_SITES,
    ):
        chunk = water_segment[
            start:
            start + WATER_SITES
        ]

        identities = {
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

        if (
            len(chunk) != WATER_SITES
            or len(identities) != 1
        ):
            water_chunk_consistency = False
            break

    return {
        "atom_count": atom_count,
        "segments": segment_summary,
        "water_chunk_consistency": (
            water_chunk_consistency
        ),
    }


def identify_water_molecule(
    topology: dict[str, Any],
) -> tuple[
    str,
    dict[str, Any],
]:
    molecule_counts = {
        entry[
            "name"
        ]: entry[
            "count"
        ]
        for entry in topology[
            "molecules"
        ]
    }

    candidates = []

    for name, definition in topology[
        "moleculetypes"
    ].items():
        count = molecule_counts.get(
            name,
            0,
        )

        atom_count = len(
            definition[
                "atoms"
            ]
        )

        score = 0

        if name.upper() == "SOL":
            score += 100

        if count > 1000:
            score += 20

        if atom_count == WATER_SITES:
            score += 20

        oxygen_like = sum(
            (
                atom[
                    "atom_name"
                ].upper().startswith(
                    "O"
                )
                or "OW"
                in atom[
                    "type"
                ].upper()
            )
            for atom in definition[
                "atoms"
            ]
        )

        if oxygen_like == 1:
            score += 10

        if score > 0:
            candidates.append(
                (
                    score,
                    count,
                    name,
                    definition,
                )
            )

    if not candidates:
        raise RuntimeError(
            "Could not identify the water molecule type."
        )

    candidates.sort(
        key=lambda item: (
            item[0],
            item[1],
        ),
        reverse=True,
    )

    best = candidates[
        0
    ]

    if (
        len(
            candidates
        )
        > 1
        and candidates[
            1
        ][
            0
        ]
        == best[
            0
        ]
    ):
        raise RuntimeError(
            "Water molecule identification is ambiguous: "
            + " | ".join(
                item[
                    2
                ]
                for item in candidates[
                    :5
                ]
            )
        )

    return (
        best[
            2
        ],
        best[
            3
        ],
    )


def identify_water_oxygen(
    water_definition: dict[str, Any],
) -> dict[str, Any]:
    atoms = water_definition[
        "atoms"
    ]

    candidates = [
        atom
        for atom in atoms
        if (
            atom[
                "atom_name"
            ].upper().startswith(
                "O"
            )
            or "OW"
            in atom[
                "type"
            ].upper()
        )
    ]

    if len(candidates) == 1:
        return candidates[
            0
        ]

    negative_heavy = [
        atom
        for atom in atoms
        if (
            atom[
                "charge"
            ]
            < -0.1
            and (
                math.isnan(
                    atom[
                        "mass"
                    ]
                )
                or atom[
                    "mass"
                ]
                > 10.0
            )
        )
    ]

    if len(negative_heavy) == 1:
        return negative_heavy[
            0
        ]

    raise RuntimeError(
        "Could not uniquely identify the "
        "water-oxygen atom."
    )


def write_csv(
    path: Path,
    rows: list[dict[str, Any]],
) -> None:
    if not rows:
        raise RuntimeError(
            f"No rows available for {path}"
        )

    fieldnames = []

    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(
                    key
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

        for row in rows:
            writer.writerow(
                {
                    field: row.get(
                        field,
                        ""
                    )
                    for field in fieldnames
                }
            )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        PROCESSED_TOP,
        R1_GRO,
        R1_DEFINITION,
    ):
        require_file(
            required
        )

    definition = json.loads(
        R1_DEFINITION.read_text(
            encoding="utf-8"
        )
    )

    retained_waters = int(
        definition[
            "retained_water_molecules"
        ]
    )

    cap_beads_per_end = int(
        definition[
            "cap_beads_per_end"
        ]
    )

    total_cap_beads = int(
        definition[
            "total_cap_beads"
        ]
    )

    topology = parse_processed_topology(
        PROCESSED_TOP
    )

    defaults = topology[
        "defaults"
    ]

    gro_audit = read_gro_segments(
        R1_GRO,
        retained_waters,
        cap_beads_per_end,
    )

    water_molecule_name, water_definition = (
        identify_water_molecule(
            topology
        )
    )

    water_oxygen = identify_water_oxygen(
        water_definition
    )

    water_oxygen_type = water_oxygen[
        "type"
    ]

    if (
        water_oxygen_type
        not in topology[
            "atomtypes"
        ]
    ):
        raise RuntimeError(
            "The identified water-oxygen type "
            f"{water_oxygen_type!r} is absent from "
            "[ atomtypes ]."
        )

    molecule_inventory = []

    baseline_total_atoms = 0
    active_atomtypes: set[str] = set()

    for entry in topology[
        "molecules"
    ]:
        name = entry[
            "name"
        ]

        count = entry[
            "count"
        ]

        if name not in topology[
            "moleculetypes"
        ]:
            raise RuntimeError(
                "System molecule lacks a definition: "
                f"{name}"
            )

        definition_entry = topology[
            "moleculetypes"
        ][
            name
        ]

        atoms_per_molecule = len(
            definition_entry[
                "atoms"
            ]
        )

        contribution = (
            count
            * atoms_per_molecule
        )

        baseline_total_atoms += (
            contribution
        )

        molecule_charge = float(
            sum(
                atom[
                    "charge"
                ]
                for atom in definition_entry[
                    "atoms"
                ]
            )
        )

        molecule_inventory.append(
            {
                "molecule": name,
                "count": count,
                "atoms_per_molecule": (
                    atoms_per_molecule
                ),
                "atom_contribution": (
                    contribution
                ),
                "molecule_charge_e": (
                    molecule_charge
                ),
                "nrexcl": definition_entry[
                    "nrexcl"
                ],
                "is_identified_water": (
                    name
                    == water_molecule_name
                ),
            }
        )

        for atom in definition_entry[
            "atoms"
        ]:
            active_atomtypes.add(
                atom[
                    "type"
                ]
            )

    write_csv(
        TOPOLOGY_INVENTORY_CSV,
        molecule_inventory,
    )

    atomtype_rows = []

    for atomtype_name in sorted(
        active_atomtypes
    ):
        if (
            atomtype_name
            not in topology[
                "atomtypes"
            ]
        ):
            raise RuntimeError(
                "Active atom type has no definition: "
                f"{atomtype_name}"
            )

        entry = topology[
            "atomtypes"
        ][
            atomtype_name
        ]

        atomtype_rows.append(
            {
                "atomtype": atomtype_name,
                "ptype": entry[
                    "ptype"
                ],
                "parameter_v": entry[
                    "parameter_v"
                ],
                "parameter_w": entry[
                    "parameter_w"
                ],
                "used_by_water_oxygen": (
                    atomtype_name
                    == water_oxygen_type
                ),
                "raw_definition": entry[
                    "raw"
                ],
            }
        )

    write_csv(
        ATOMTYPE_INVENTORY_CSV,
        atomtype_rows,
    )

    proposed_molecules = []

    proposed_total_atoms = 0

    for entry in topology[
        "molecules"
    ]:
        name = entry[
            "name"
        ]

        original_count = entry[
            "count"
        ]

        proposed_count = (
            retained_waters
            if name
            == water_molecule_name
            else original_count
        )

        atoms_per_molecule = len(
            topology[
                "moleculetypes"
            ][
                name
            ][
                "atoms"
            ]
        )

        contribution = (
            proposed_count
            * atoms_per_molecule
        )

        proposed_total_atoms += (
            contribution
        )

        proposed_molecules.append(
            {
                "order": len(
                    proposed_molecules
                )
                + 1,
                "molecule": name,
                "original_count": (
                    original_count
                ),
                "proposed_count": (
                    proposed_count
                ),
                "atoms_per_molecule": (
                    atoms_per_molecule
                ),
                "proposed_atom_contribution": (
                    contribution
                ),
                "action": (
                    "REPLACE_WATER_COUNT"
                    if name
                    == water_molecule_name
                    else "RETAIN"
                ),
            }
        )

    for cap_name in (
        "CAPL",
        "CAPU",
    ):
        proposed_molecules.append(
            {
                "order": len(
                    proposed_molecules
                )
                + 1,
                "molecule": cap_name,
                "original_count": 0,
                "proposed_count": 1,
                "atoms_per_molecule": (
                    cap_beads_per_end
                ),
                "proposed_atom_contribution": (
                    cap_beads_per_end
                ),
                "action": "ADD_FROZEN_CAP",
            }
        )

        proposed_total_atoms += (
            cap_beads_per_end
        )

    write_csv(
        PROPOSED_MOLECULES_CSV,
        proposed_molecules,
    )

    kbt_kj_mol = (
        GAS_CONSTANT_KJ_MOL_K
        * TEMPERATURE_K
    )

    calibration_rows = []

    for energy_multiple in (
        ENERGY_LEVELS_KBT
    ):
        target_energy = (
            energy_multiple
            * kbt_kj_mol
        )

        c12 = (
            target_energy
            * TARGET_REPULSION_RADIUS_NM
            ** 12
        )

        target_force = (
            12.0
            * target_energy
            / TARGET_REPULSION_RADIUS_NM
        )

        for radius_nm in (
            DIAGNOSTIC_RADII_NM
        ):
            energy_kj_mol = (
                c12
                / radius_nm
                ** 12
            )

            force_kj_mol_nm = (
                12.0
                * c12
                / radius_nm
                ** 13
            )

            calibration_rows.append(
                {
                    "target_energy_kBT_at_0p17nm": (
                        energy_multiple
                    ),
                    "target_energy_kJ_mol": (
                        target_energy
                    ),
                    "C6_kJ_mol_nm6": 0.0,
                    "C12_kJ_mol_nm12": (
                        c12
                    ),
                    "diagnostic_radius_nm": (
                        radius_nm
                    ),
                    "potential_kJ_mol": (
                        energy_kj_mol
                    ),
                    "potential_kBT": (
                        energy_kj_mol
                        / kbt_kj_mol
                    ),
                    "force_kJ_mol_nm": (
                        force_kj_mol_nm
                    ),
                    "target_force_at_0p17nm_kJ_mol_nm": (
                        target_force
                    ),
                }
            )

    write_csv(
        REPULSION_CALIBRATION_CSV,
        calibration_rows,
    )

    nbfunc = int(
        defaults[
            "nbfunc"
        ]
    )

    comb_rule = int(
        defaults[
            "comb_rule"
        ]
    )

    pure_r12_standard_override_feasible = bool(
        nbfunc == 1
        and comb_rule == 1
    )

    existing_cap_atomtype = any(
        atomtype_name.upper()
        in {
            "CAP",
            "CAPL",
            "CAPU",
        }
        for atomtype_name in topology[
            "atomtypes"
        ]
    )

    existing_cap_nonbond_override = any(
        (
            entry[
                "type_i"
            ].upper()
            in {
                "CAP",
                "CAPL",
                "CAPU",
            }
            or entry[
                "type_j"
            ].upper()
            in {
                "CAP",
                "CAPL",
                "CAPU",
            }
        )
        for entry in topology[
            "nonbond_params"
        ]
    )

    gates = {
        "baseline_atom_count_is_68320": (
            baseline_total_atoms
            == EXPECTED_R0_ATOMS
        ),
        "R1_GRO_atom_count_is_68314": (
            gro_audit[
                "atom_count"
            ]
            == EXPECTED_R1_ATOMS
        ),
        "proposed_topology_atom_count_matches_R1_GRO": (
            proposed_total_atoms
            == gro_audit[
                "atom_count"
            ]
        ),
        "water_molecule_has_four_sites": (
            len(
                water_definition[
                    "atoms"
                ]
            )
            == WATER_SITES
        ),
        "water_oxygen_type_is_unique": (
            bool(
                water_oxygen_type
            )
        ),
        "R1_water_chunks_are_consistent": (
            gro_audit[
                "water_chunk_consistency"
            ]
        ),
        "lower_cap_has_expected_resname": (
            set(
                gro_audit[
                    "segments"
                ][
                    "lower_cap"
                ][
                    "resnames"
                ]
            )
            == {
                "CPL"
            }
        ),
        "upper_cap_has_expected_resname": (
            set(
                gro_audit[
                    "segments"
                ][
                    "upper_cap"
                ][
                    "resnames"
                ]
            )
            == {
                "CPU"
            }
        ),
        "no_existing_cap_atomtype_collision": (
            not existing_cap_atomtype
        ),
        "no_existing_cap_nonbond_override_collision": (
            not existing_cap_nonbond_override
        ),
    }

    failed_gates = [
        name
        for name, passed
        in gates.items()
        if not passed
    ]

    topology_contract_valid = (
        len(
            failed_gates
        )
        == 0
    )

    if (
        topology_contract_valid
        and pure_r12_standard_override_feasible
    ):
        decision = (
            "PURE_R12_CAP_WATER_MODEL_FEASIBLE"
        )

        required_next_step = (
            "BUILD_R1_CAP_TOPOLOGY_AND_RUN_STATIC_"
            "REPULSION_ENERGY_SCAN"
        )
    elif topology_contract_valid:
        decision = (
            "TOPOLOGY_CONTRACT_VALID_BUT_PURE_R12_"
            "REQUIRES_ALTERNATIVE_IMPLEMENTATION"
        )

        required_next_step = (
            "SELECT_TABULATED_OR_VALIDATED_"
            "ALTERNATIVE_CAP_INTERACTION"
        )
    else:
        decision = (
            "R1_TOPOLOGY_CONTRACT_BLOCKED"
        )

        required_next_step = (
            "RESOLVE_R1_TOPOLOGY_GATE_FAILURES"
        )

    parameter_semantics = (
        "C6_C12"
        if comb_rule == 1
        else (
            "SIGMA_EPSILON"
            if comb_rule
            in {
                2,
                3,
            }
            else "UNSUPPORTED"
        )
    )

    contract = {
        "decision": decision,
        "source_processed_topology": relative(
            PROCESSED_TOP
        ),
        "source_R1_geometry": relative(
            R1_GRO
        ),
        "nbfunc": nbfunc,
        "combination_rule": comb_rule,
        "parameter_semantics": (
            parameter_semantics
        ),
        "water_molecule_name": (
            water_molecule_name
        ),
        "water_sites_per_molecule": len(
            water_definition[
                "atoms"
            ]
        ),
        "water_oxygen_atom_name": (
            water_oxygen[
                "atom_name"
            ]
        ),
        "water_oxygen_atom_type": (
            water_oxygen_type
        ),
        "water_oxygen_charge_e": (
            water_oxygen[
                "charge"
            ]
        ),
        "retained_water_molecules": (
            retained_waters
        ),
        "cap_beads_per_end": (
            cap_beads_per_end
        ),
        "total_cap_beads": (
            total_cap_beads
        ),
        "proposed_cap_atomtype": "CAP",
        "proposed_cap_charge_e": 0.0,
        "proposed_cap_mass_u": 12.011,
        "proposed_cap_particle_type": "A",
        "proposed_cap_base_parameter_v": 0.0,
        "proposed_cap_base_parameter_w": 0.0,
        "proposed_cap_moleculetypes": [
            "CAPL",
            "CAPU",
        ],
        "proposed_cap_atoms_per_molecule": (
            cap_beads_per_end
        ),
        "proposed_cap_cap_interaction": (
            "ZERO"
        ),
        "proposed_cap_HBN_interaction": (
            "ZERO"
        ),
        "proposed_cap_PYR_interaction": (
            "ZERO"
        ),
        "proposed_cap_water_H_or_M_interaction": (
            "ZERO"
        ),
        "proposed_cap_water_O_interaction": (
            "PURE_R12_C6_ZERO_C12_TO_BE_SELECTED"
            if pure_r12_standard_override_feasible
            else
            "NOT_YET_SELECTED"
        ),
        "target_repulsion_radius_nm": (
            TARGET_REPULSION_RADIUS_NM
        ),
        "temperature_K": (
            TEMPERATURE_K
        ),
        "topology_contract_valid": (
            topology_contract_valid
        ),
        "pure_r12_standard_nonbond_override_feasible": (
            pure_r12_standard_override_feasible
        ),
        "static_energy_scan_required": True,
        "energy_minimization_authorized": False,
        "MD_execution_authorized": False,
        "required_next_step": (
            required_next_step
        ),
    }

    CONTRACT_JSON.write_text(
        json.dumps(
            contract,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    summary = {
        "decision": decision,
        "nbfunc": nbfunc,
        "combination_rule": comb_rule,
        "parameter_semantics": (
            parameter_semantics
        ),
        "gen_pairs": defaults[
            "gen_pairs"
        ],
        "baseline_moleculetype_count": len(
            topology[
                "moleculetypes"
            ]
        ),
        "baseline_system_molecule_entries": len(
            topology[
                "molecules"
            ]
        ),
        "baseline_atomtype_count": len(
            topology[
                "atomtypes"
            ]
        ),
        "active_atomtype_count": len(
            active_atomtypes
        ),
        "explicit_nonbond_override_count": len(
            topology[
                "nonbond_params"
            ]
        ),
        "baseline_total_atoms": (
            baseline_total_atoms
        ),
        "R1_GRO_atoms": gro_audit[
            "atom_count"
        ],
        "proposed_topology_atoms": (
            proposed_total_atoms
        ),
        "water_molecule_name": (
            water_molecule_name
        ),
        "water_sites": len(
            water_definition[
                "atoms"
            ]
        ),
        "water_oxygen_atom_name": (
            water_oxygen[
                "atom_name"
            ]
        ),
        "water_oxygen_atom_type": (
            water_oxygen_type
        ),
        "water_oxygen_charge_e": (
            water_oxygen[
                "charge"
            ]
        ),
        "retained_waters": (
            retained_waters
        ),
        "cap_beads_per_end": (
            cap_beads_per_end
        ),
        "total_cap_beads": (
            total_cap_beads
        ),
        "water_chunk_consistency": (
            gro_audit[
                "water_chunk_consistency"
            ]
        ),
        "topology_contract_valid": (
            topology_contract_valid
        ),
        "pure_r12_standard_override_feasible": (
            pure_r12_standard_override_feasible
        ),
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

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed
        in gates.items()
    )

    molecule_lines = "\n".join(
        (
            f"- `{row['molecule']}`: "
            f"{row['count']} × "
            f"{row['atoms_per_molecule']} atoms "
            f"= {row['atom_contribution']}"
            + (
                " — identified water"
                if row[
                    "is_identified_water"
                ]
                else ""
            )
        )
        for row in molecule_inventory
    )

    proposed_lines = "\n".join(
        (
            f"- `{row['molecule']}`: "
            f"{row['proposed_count']} × "
            f"{row['atoms_per_molecule']} "
            f"= {row['proposed_atom_contribution']} "
            f"({row['action']})"
        )
        for row in proposed_molecules
    )

    REPORT_MD.write_text(
        f"""# R1 Cap Topology Requirements

## Purpose

This audit identifies the exact force-field conventions required to
add the two R1 steric caps without introducing uncontrolled attraction,
charge, or cap–solute interactions.

No topology was modified and no GROMACS preprocessing or dynamics was
performed.

## Baseline nonbonded convention

- `nbfunc`: **{nbfunc}**
- `comb-rule`: **{comb_rule}**
- parameter semantics: **{parameter_semantics}**
- `gen-pairs`: **{defaults['gen_pairs']}**
- active atom types: **{len(active_atomtypes)}**
- explicit `[ nonbond_params ]` entries:
  **{len(topology['nonbond_params'])}**

## Baseline molecule inventory

{molecule_lines}

- Baseline atom count:
  **{baseline_total_atoms}/{EXPECTED_R0_ATOMS}**

## Water model identification

- Molecule type:
  **{water_molecule_name}**
- Sites per molecule:
  **{len(water_definition['atoms'])}**
- Oxygen atom:
  **{water_oxygen['atom_name']}**
- Oxygen atom type:
  **{water_oxygen_type}**
- Oxygen charge:
  **{water_oxygen['charge']:.8f} e**

## R1 coordinate ordering

- R1 GRO atoms:
  **{gro_audit['atom_count']}**
- Retained waters:
  **{retained_waters}**
- Lower-cap atoms:
  **{cap_beads_per_end}**
- Upper-cap atoms:
  **{cap_beads_per_end}**
- Water four-site ordering:
  **{'PASS' if gro_audit['water_chunk_consistency'] else 'FAIL'}**

## Proposed R1 molecule inventory

{proposed_lines}

- Proposed topology atom count:
  **{proposed_total_atoms}**
- R1 coordinate atom count:
  **{gro_audit['atom_count']}**

## Proposed cap contract

The provisional cap model will contain:

- atom type: `CAP`;
- charge: **0 e**;
- mass: **12.011 u**;
- base nonbonded parameters: **zero**;
- cap–cap interaction: **zero**;
- cap–HBN interaction: **zero**;
- cap–PYR interaction: **zero**;
- cap–water H/M interaction: **zero**;
- cap–water O interaction:
  **{'pure repulsive r^-12 candidate' if pure_r12_standard_override_feasible else 'requires an alternative implementation'}**;
- two molecule types: `CAPL` and `CAPU`;
- all cap coordinates frozen during the screening.

No C12 value is authorized yet. Candidate values are recorded in:

`{relative(REPULSION_CALIBRATION_CSV)}`

They must be tested with static energy and force scans before energy
minimization or MD.

## Gates

{gate_lines}

## Decision

- Decision: **{decision}**
- Topology contract valid:
  **{'YES' if topology_contract_valid else 'NO'}**
- Standard pure-r12 override feasible:
  **{'YES' if pure_r12_standard_override_feasible else 'NO'}**
- Energy minimization authorized: **NO**
- MD execution authorized: **NO**
- Required next step:
  `{required_next_step}`
""",
        encoding="utf-8",
    )

    print(
        "Day023 R1 cap-topology requirement audit completed."
    )

    print(
        "Baseline nbfunc / comb-rule / semantics: "
        f"{nbfunc} / {comb_rule} / "
        f"{parameter_semantics}"
    )

    print(
        "Baseline moleculetypes / system entries / "
        "active atomtypes / nonbond overrides: "
        f"{len(topology['moleculetypes'])}/"
        f"{len(topology['molecules'])}/"
        f"{len(active_atomtypes)}/"
        f"{len(topology['nonbond_params'])}"
    )

    print(
        "Baseline total atoms: "
        f"{baseline_total_atoms}/"
        f"{EXPECTED_R0_ATOMS}"
    )

    print(
        "Water molecule / sites / oxygen atom / type: "
        f"{water_molecule_name} / "
        f"{len(water_definition['atoms'])} / "
        f"{water_oxygen['atom_name']} / "
        f"{water_oxygen_type}"
    )

    print(
        "Water oxygen charge: "
        f"{water_oxygen['charge']:.8f} e"
    )

    print(
        "R1 GRO / proposed topology atoms: "
        f"{gro_audit['atom_count']}/"
        f"{proposed_total_atoms}"
    )

    print(
        "Retained waters / cap beads per end / total: "
        f"{retained_waters}/"
        f"{cap_beads_per_end}/"
        f"{total_cap_beads}"
    )

    print(
        "R1 water ordering / cap residue ordering: "
        f"{'PASS' if gro_audit['water_chunk_consistency'] else 'FAIL'} / "
        f"{'PASS' if gates['lower_cap_has_expected_resname'] and gates['upper_cap_has_expected_resname'] else 'FAIL'}"
    )

    print(
        "Topology contract valid: "
        f"{'YES' if topology_contract_valid else 'NO'}"
    )

    print(
        "Pure r^-12 standard override feasible: "
        f"{'YES' if pure_r12_standard_override_feasible else 'NO'}"
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
        f"Decision: {decision}"
    )

    print(
        "Energy minimization authorized: NO"
    )

    print(
        "MD execution authorized: NO"
    )

    print(
        f"Required next step: {required_next_step}"
    )

    print(
        f"Wrote: {relative(TOPOLOGY_INVENTORY_CSV)}"
    )

    print(
        f"Wrote: {relative(ATOMTYPE_INVENTORY_CSV)}"
    )

    print(
        f"Wrote: {relative(PROPOSED_MOLECULES_CSV)}"
    )

    print(
        f"Wrote: {relative(REPULSION_CALIBRATION_CSV)}"
    )

    print(
        f"Wrote: {relative(SUMMARY_CSV)}"
    )

    print(
        f"Wrote: {relative(CONTRACT_JSON)}"
    )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )

    if failed_gates:
        raise RuntimeError(
            "R1 topology requirements remain blocked."
        )


if __name__ == "__main__":
    main()
