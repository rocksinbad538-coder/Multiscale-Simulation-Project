#!/usr/bin/env python3
"""
DAY040 / D040-A4

Local hydrogen topology-augmentation design.

This block determines whether the six retained H atoms in the adopted
37-atom working charge model are physical passivants that must exist in
the final classical topology. It identifies their B/N parents using the
authoritative adopted-fragment connectivity whenever available and
cross-checks the assignment geometrically.

No accepted topology is modified.
No coordinates are inserted.
No atom type or bonded parameter is adopted.
No GROMACS calculation is executed.
"""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from resp_common import load_json, require_file, sha256


ROOT = Path(__file__).resolve().parents[2]

A3_REPORT = (
    ROOT
    / "runs/phase1A/day040_phase1A_G_local_mapping_feasibility"
    / "QM_F06_UPPER_V7A_R1_LOCAL_37_ATOM_MAPPING_FEASIBILITY.json"
)

EXPECTED_A3_DECISION = (
    "D040_A3_LOCAL_37_ATOM_MAPPING_FEASIBILITY_PASS_"
    "HYDROGEN_TOPOLOGY_AUGMENTATION_REVIEW_AUTHORIZED"
)

CLOSURE_POINTER = (
    ROOT
    / "runs/phase1A"
    / "LATEST_PHASE1A_F_CHARGE_MODEL_CLOSURE.txt"
)

NOMINAL_EDGES = (
    ROOT
    / "runs/phase1A/day035_qm_f06_upper_v7a_r1_coordinate_adoption"
    / "QM_F06_UPPER_V7A_ADOPTED_nominal_edges.csv"
)

PROVENANCE_MAP = (
    ROOT
    / "runs/phase1A/day035_qm_f06_upper_v7a_r1_coordinate_adoption"
    / "QM_F06_UPPER_V7A_ADOPTED_atom_role_provenance_map.csv"
)

OUTPUT_DIR = (
    ROOT
    / "runs/phase1A/day040_phase1A_G_hydrogen_augmentation_design"
)

PARENT_MAP_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_HYDROGEN_PARENT_MAPPING.csv"
)

PARAMETER_REQUIREMENTS_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_HYDROGEN_PARAMETER_REQUIREMENTS.csv"
)

REPORT_JSON = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_HYDROGEN_AUGMENTATION_DESIGN.json"
)

EXPECTED_H_COUNT = 6
EXPECTED_TOTAL_COUNT = 37

GEOMETRIC_PARENT_CUTOFF_A = 1.35
GEOMETRIC_SECOND_PARENT_GAP_A = 0.15

B_H_REVIEW_RANGE_A = (0.90, 1.35)
N_H_REVIEW_RANGE_A = (0.85, 1.25)


def detect_column(
    columns: list[str],
    candidates: tuple[str, ...],
    description: str,
) -> str:
    exact = {
        column.lower(): column
        for column in columns
    }

    for candidate in candidates:
        if candidate.lower() in exact:
            return exact[candidate.lower()]

    normalized = {
        "".join(
            character
            for character in column.lower()
            if character.isalnum()
        ): column
        for column in columns
    }

    for candidate in candidates:
        key = "".join(
            character
            for character in candidate.lower()
            if character.isalnum()
        )

        if key in normalized:
            return normalized[key]

    raise RuntimeError(
        f"Could not identify {description}.\n"
        f"Available columns: {columns}"
    )


def load_csv(path: Path) -> list[dict]:
    require_file(path)

    with path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(
            csv.DictReader(handle)
        )


def load_adopted_atoms(path: Path) -> list[dict]:
    rows = load_csv(path)

    rows.sort(
        key=lambda row: int(
            row["real_atom_sequence_index"]
        )
    )

    parsed = []

    for row in rows:
        parsed.append(
            {
                "real_atom_sequence_index": int(
                    row["real_atom_sequence_index"]
                ),
                "original_atom_index_0based": int(
                    row["original_atom_index_0based"]
                ),
                "original_atom_index_1based": int(
                    row["original_atom_index_1based"]
                ),
                "atom_id": row["atom_id"],
                "element": row["element"],
                "atom_role": row["atom_role"],
                "node_type": row["node_type"],
                "transfer_status": row["transfer_status"],
                "x_A": float(row["x_A"]),
                "y_A": float(row["y_A"]),
                "z_A": float(row["z_A"]),
                "charge_e": float(
                    row["adopted_working_charge_e"]
                ),
            }
        )

    return parsed


def distance_A(first: dict, second: dict) -> float:
    first_xyz = np.asarray(
        [
            first["x_A"],
            first["y_A"],
            first["z_A"],
        ],
        dtype=float,
    )

    second_xyz = np.asarray(
        [
            second["x_A"],
            second["y_A"],
            second["z_A"],
        ],
        dtype=float,
    )

    return float(
        np.linalg.norm(
            first_xyz - second_xyz
        )
    )


def classify_bond_review(
    parent_element: str,
    bond_length_A: float,
) -> tuple[str, bool]:
    if parent_element == "B":
        lower, upper = B_H_REVIEW_RANGE_A
        bond_class = "B-H"
    elif parent_element == "N":
        lower, upper = N_H_REVIEW_RANGE_A
        bond_class = "N-H"
    else:
        return "UNKNOWN-H", False

    return (
        bond_class,
        lower <= bond_length_A <= upper,
    )


print("=" * 100)
print("DAY040 / D040-A4 — LOCAL HYDROGEN TOPOLOGY-AUGMENTATION DESIGN")
print("=" * 100)


print("\n[1] UPSTREAM AUTHORIZATION")

require_file(A3_REPORT)
require_file(CLOSURE_POINTER)
require_file(NOMINAL_EDGES)
require_file(PROVENANCE_MAP)

a3_report = load_json(
    A3_REPORT
)

if (
    a3_report.get("decision")
    != EXPECTED_A3_DECISION
):
    raise RuntimeError(
        "Unexpected D040-A3 decision.\n"
        f"Observed: {a3_report.get('decision')}"
    )

if (
    a3_report.get(
        "authorizations",
        {},
    ).get(
        "hydrogen_topology_augmentation_review_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Hydrogen topology-augmentation review is not authorized"
    )

closure_dir = (
    ROOT
    / CLOSURE_POINTER.read_text(
        encoding="utf-8"
    ).strip()
)

charges_csv = (
    closure_dir
    / "QM_F06_UPPER_V7A_R1_PHASE1A_F_ADOPTED_WORKING_CHARGES.csv"
)

require_file(charges_csv)

print("D040_A3_decision_gate = PASS")
print("hydrogen_augmentation_review_gate = PASS")
print("topology_modification_blocked_gate = PASS")
print("parameter_adoption_blocked_gate = PASS")


print("\n[2] LOAD ADOPTED ATOMS")

atoms = load_adopted_atoms(
    charges_csv
)

if len(atoms) != EXPECTED_TOTAL_COUNT:
    raise RuntimeError(
        f"Expected {EXPECTED_TOTAL_COUNT} adopted atoms"
    )

atom_by_id = {
    atom["atom_id"]: atom
    for atom in atoms
}

hydrogens = [
    atom
    for atom in atoms
    if atom["element"] == "H"
]

heavy_atoms = [
    atom
    for atom in atoms
    if atom["element"] in {"B", "N"}
]

if len(hydrogens) != EXPECTED_H_COUNT:
    raise RuntimeError(
        f"Expected {EXPECTED_H_COUNT} retained H atoms"
    )

print(f"adopted_atom_count = {len(atoms)}")
print(f"retained_hydrogen_count = {len(hydrogens)}")
print(f"retained_BN_count = {len(heavy_atoms)}")

for hydrogen in hydrogens:
    print(
        f"H atom_id={hydrogen['atom_id']} "
        f"role={hydrogen['atom_role']} "
        f"charge_e={hydrogen['charge_e']:.12g}"
    )


print("\n[3] LOAD AUTHORITATIVE NOMINAL EDGES")

edge_rows = load_csv(
    NOMINAL_EDGES
)

if not edge_rows:
    raise RuntimeError(
        "Nominal-edge table is empty"
    )

edge_columns = list(
    edge_rows[0].keys()
)

first_atom_col = detect_column(
    edge_columns,
    (
        "first_atom",
        "atom1",
        "first_atom_id",
        "atom_i",
    ),
    "first atom column",
)

second_atom_col = detect_column(
    edge_columns,
    (
        "second_atom",
        "atom2",
        "second_atom_id",
        "atom_j",
    ),
    "second atom column",
)

distance_col = detect_column(
    edge_columns,
    (
        "distance_A",
        "distance",
        "bond_length_A",
    ),
    "edge distance column",
)

edge_type_col = detect_column(
    edge_columns,
    (
        "edge_type",
        "bond_type",
        "type",
    ),
    "edge type column",
)

provenance_col = detect_column(
    edge_columns,
    (
        "provenance",
        "source",
        "basis",
    ),
    "edge provenance column",
)

print(f"nominal_edge_count = {len(edge_rows)}")
print(f"first_atom_column = {first_atom_col}")
print(f"second_atom_column = {second_atom_col}")
print(f"distance_column = {distance_col}")
print(f"edge_type_column = {edge_type_col}")
print(f"provenance_column = {provenance_col}")


print("\n[4] HYDROGEN PARENT ASSIGNMENT")

parent_records = []

for hydrogen in hydrogens:
    hydrogen_id = hydrogen["atom_id"]

    connectivity_candidates = []

    for edge in edge_rows:
        first_id = edge[first_atom_col]
        second_id = edge[second_atom_col]

        if first_id == hydrogen_id:
            neighbor_id = second_id
        elif second_id == hydrogen_id:
            neighbor_id = first_id
        else:
            continue

        neighbor = atom_by_id.get(
            neighbor_id
        )

        if neighbor is None:
            continue

        if neighbor["element"] not in {
            "B",
            "N",
        }:
            continue

        connectivity_candidates.append(
            {
                "parent_atom_id": neighbor_id,
                "parent_element": neighbor["element"],
                "edge_distance_A": float(
                    edge[distance_col]
                ),
                "edge_type": edge[edge_type_col],
                "edge_provenance": edge[
                    provenance_col
                ],
                "parent_atom": neighbor,
            }
        )

    geometric_candidates = sorted(
        [
            {
                "parent_atom_id": heavy["atom_id"],
                "parent_element": heavy["element"],
                "distance_A": distance_A(
                    hydrogen,
                    heavy,
                ),
                "parent_atom": heavy,
            }
            for heavy in heavy_atoms
        ],
        key=lambda record: record[
            "distance_A"
        ],
    )

    nearest = geometric_candidates[0]
    second_nearest = geometric_candidates[1]

    if len(connectivity_candidates) == 1:
        selected = connectivity_candidates[0]
        assignment_basis = (
            "AUTHORITATIVE_NOMINAL_EDGE"
        )
        selected_distance_A = distance_A(
            hydrogen,
            selected["parent_atom"],
        )
        nominal_edge_distance_A = selected[
            "edge_distance_A"
        ]
        edge_type = selected["edge_type"]
        edge_provenance = selected[
            "edge_provenance"
        ]
    elif len(connectivity_candidates) > 1:
        connectivity_candidates.sort(
            key=lambda record: distance_A(
                hydrogen,
                record["parent_atom"],
            )
        )
        selected = connectivity_candidates[0]
        assignment_basis = (
            "MULTIPLE_NOMINAL_EDGES_NEAREST_SELECTED"
        )
        selected_distance_A = distance_A(
            hydrogen,
            selected["parent_atom"],
        )
        nominal_edge_distance_A = selected[
            "edge_distance_A"
        ]
        edge_type = selected["edge_type"]
        edge_provenance = selected[
            "edge_provenance"
        ]
    else:
        selected = nearest
        assignment_basis = (
            "GEOMETRIC_NEAREST_BN_FALLBACK"
        )
        selected_distance_A = nearest[
            "distance_A"
        ]
        nominal_edge_distance_A = None
        edge_type = ""
        edge_provenance = ""

    (
        bond_class,
        bond_length_review_pass,
    ) = classify_bond_review(
        selected["parent_element"],
        selected_distance_A,
    )

    nearest_matches_selected = (
        nearest["parent_atom_id"]
        == selected["parent_atom_id"]
    )

    second_parent_gap_A = (
        second_nearest["distance_A"]
        - nearest["distance_A"]
    )

    geometrically_unique = (
        nearest["distance_A"]
        <= GEOMETRIC_PARENT_CUTOFF_A
        and second_parent_gap_A
        >= GEOMETRIC_SECOND_PARENT_GAP_A
    )

    parent_record = {
        "hydrogen_real_atom_sequence_index": (
            hydrogen[
                "real_atom_sequence_index"
            ]
        ),
        "hydrogen_original_atom_index_0based": (
            hydrogen[
                "original_atom_index_0based"
            ]
        ),
        "hydrogen_atom_id": hydrogen_id,
        "hydrogen_role": hydrogen[
            "atom_role"
        ],
        "hydrogen_charge_e": (
            hydrogen["charge_e"]
        ),
        "parent_atom_id": (
            selected["parent_atom_id"]
        ),
        "parent_element": (
            selected["parent_element"]
        ),
        "parent_real_atom_sequence_index": (
            selected["parent_atom"][
                "real_atom_sequence_index"
            ]
        ),
        "parent_original_atom_index_0based": (
            selected["parent_atom"][
                "original_atom_index_0based"
            ]
        ),
        "parent_charge_e": (
            selected["parent_atom"][
                "charge_e"
            ]
        ),
        "assignment_basis": (
            assignment_basis
        ),
        "nominal_connectivity_candidate_count": (
            len(connectivity_candidates)
        ),
        "selected_distance_A": (
            selected_distance_A
        ),
        "nominal_edge_distance_A": (
            nominal_edge_distance_A
        ),
        "edge_type": edge_type,
        "edge_provenance": (
            edge_provenance
        ),
        "nearest_geometric_parent_atom_id": (
            nearest["parent_atom_id"]
        ),
        "nearest_geometric_parent_element": (
            nearest["parent_element"]
        ),
        "nearest_geometric_distance_A": (
            nearest["distance_A"]
        ),
        "second_nearest_parent_atom_id": (
            second_nearest[
                "parent_atom_id"
            ]
        ),
        "second_nearest_distance_A": (
            second_nearest[
                "distance_A"
            ]
        ),
        "second_parent_gap_A": (
            second_parent_gap_A
        ),
        "nearest_matches_selected": (
            nearest_matches_selected
        ),
        "geometrically_unique": (
            geometrically_unique
        ),
        "bond_class": bond_class,
        "bond_length_review_pass": (
            bond_length_review_pass
        ),
    }

    parent_records.append(
        parent_record
    )

    print(
        f"H={hydrogen_id:<28s} "
        f"parent={selected['parent_atom_id']:<28s} "
        f"parent_element={selected['parent_element']} "
        f"distance_A={selected_distance_A:.9f} "
        f"basis={assignment_basis} "
        f"nearest_match={nearest_matches_selected} "
        f"bond_class={bond_class} "
        f"length_gate={'PASS' if bond_length_review_pass else 'FAIL'}"
    )


print("\n[5] PARENT-MAPPING CONSISTENCY")

parent_ids = [
    record["parent_atom_id"]
    for record in parent_records
]

unique_parent_count = len(
    set(parent_ids)
)

duplicate_parent_ids = sorted(
    {
        parent_id
        for parent_id in parent_ids
        if parent_ids.count(
            parent_id
        )
        > 1
    }
)

authoritative_edge_count = sum(
    record[
        "assignment_basis"
    ]
    == "AUTHORITATIVE_NOMINAL_EDGE"
    for record in parent_records
)

fallback_count = sum(
    record[
        "assignment_basis"
    ]
    == "GEOMETRIC_NEAREST_BN_FALLBACK"
    for record in parent_records
)

nearest_agreement_count = sum(
    record[
        "nearest_matches_selected"
    ]
    for record in parent_records
)

geometrically_unique_count = sum(
    record[
        "geometrically_unique"
    ]
    for record in parent_records
)

bond_length_pass_count = sum(
    record[
        "bond_length_review_pass"
    ]
    for record in parent_records
)

B_H_count = sum(
    record["bond_class"]
    == "B-H"
    for record in parent_records
)

N_H_count = sum(
    record["bond_class"]
    == "N-H"
    for record in parent_records
)

print(f"parent_record_count = {len(parent_records)}")
print(f"unique_parent_count = {unique_parent_count}")
print(f"duplicate_parent_ids = {duplicate_parent_ids}")
print(f"authoritative_edge_count = {authoritative_edge_count}")
print(f"geometric_fallback_count = {fallback_count}")
print(f"nearest_agreement_count = {nearest_agreement_count}/6")
print(f"geometrically_unique_count = {geometrically_unique_count}/6")
print(f"bond_length_pass_count = {bond_length_pass_count}/6")
print(f"B_H_count = {B_H_count}")
print(f"N_H_count = {N_H_count}")


print("\n[6] PARAMETER REQUIREMENTS")

parameter_requirements = []

for bond_class, parent_element in (
    ("B-H", "B"),
    ("N-H", "N"),
):
    relevant = [
        record
        for record in parent_records
        if record["bond_class"]
        == bond_class
    ]

    if not relevant:
        continue

    lengths = np.asarray(
        [
            record[
                "selected_distance_A"
            ]
            for record in relevant
        ],
        dtype=float,
    )

    parameter_requirements.append(
        {
            "interaction_class": (
                bond_class
            ),
            "parent_element": (
                parent_element
            ),
            "hydrogen_element": "H",
            "instance_count": (
                len(relevant)
            ),
            "observed_length_mean_A": float(
                np.mean(lengths)
            ),
            "observed_length_std_A": float(
                np.std(lengths)
            ),
            "observed_length_min_A": float(
                np.min(lengths)
            ),
            "observed_length_max_A": float(
                np.max(lengths)
            ),
            "required_new_hydrogen_atom_type": True,
            "required_nonbonded_parameters": (
                "MASS_SIGMA_EPSILON_OR_REFERENCED_BASE_TYPE"
            ),
            "required_bond_parameters": True,
            "required_angle_parameters": True,
            "required_improper_review": True,
            "parameter_adoption_status": (
                "NOT_AUTHORIZED"
            ),
        }
    )

for requirement in parameter_requirements:
    print(
        f"class={requirement['interaction_class']} "
        f"count={requirement['instance_count']} "
        f"mean_A={requirement['observed_length_mean_A']:.9f} "
        f"range_A=[{requirement['observed_length_min_A']:.9f},"
        f"{requirement['observed_length_max_A']:.9f}]"
    )


print("\n[7] PHYSICAL-PASSIVANT INTERPRETATION")

# Physical-hydrogen status must be established from chemical evidence,
# not exclusively from the atom_role label. Two retained hydrogens are
# ORIGINAL_FRAGMENT_ATOM records, but each has an authoritative nominal
# B/N-H edge, agrees with the nearest geometric B/N parent, and has a
# chemically admissible bond length.
physical_hydrogen_evidence_count = sum(
    (
        record[
            "nominal_connectivity_candidate_count"
        ]
        >= 1
        and record[
            "nearest_matches_selected"
        ]
        and record[
            "bond_length_review_pass"
        ]
    )
    for record in parent_records
)

all_H_have_parent_gate = (
    len(parent_records)
    == EXPECTED_H_COUNT
    and all(
        record[
            "parent_atom_id"
        ]
        for record in parent_records
    )
)

all_H_are_physical_passivants_gate = (
    physical_hydrogen_evidence_count
    == EXPECTED_H_COUNT
)

parent_assignment_scientifically_supported_gate = (
    all_H_have_parent_gate
    and authoritative_edge_count
    == EXPECTED_H_COUNT
    and nearest_agreement_count
    == EXPECTED_H_COUNT
    and geometrically_unique_count
    == EXPECTED_H_COUNT
    and bond_length_pass_count
    == EXPECTED_H_COUNT
)

hydrogen_augmentation_design_review_authorized = (
    all_H_are_physical_passivants_gate
    and parent_assignment_scientifically_supported_gate
)

print(
    "physical_hydrogen_evidence_count = "
    f"{physical_hydrogen_evidence_count}/"
    f"{EXPECTED_H_COUNT}"
)
print(
    "all_H_are_physical_passivants = "
    f"{all_H_are_physical_passivants_gate}"
)
print(
    "all_H_have_parent = "
    f"{all_H_have_parent_gate}"
)
print(
    "parent_assignment_scientifically_supported = "
    f"{parent_assignment_scientifically_supported_gate}"
)
print(
    "hydrogen_augmentation_design_review_authorized = "
    f"{hydrogen_augmentation_design_review_authorized}"
)


print("\n[8] WRITE OUTPUTS")

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

parent_fieldnames = [
    "hydrogen_real_atom_sequence_index",
    "hydrogen_original_atom_index_0based",
    "hydrogen_atom_id",
    "hydrogen_role",
    "hydrogen_charge_e",
    "parent_atom_id",
    "parent_element",
    "parent_real_atom_sequence_index",
    "parent_original_atom_index_0based",
    "parent_charge_e",
    "assignment_basis",
    "nominal_connectivity_candidate_count",
    "selected_distance_A",
    "nominal_edge_distance_A",
    "edge_type",
    "edge_provenance",
    "nearest_geometric_parent_atom_id",
    "nearest_geometric_parent_element",
    "nearest_geometric_distance_A",
    "second_nearest_parent_atom_id",
    "second_nearest_distance_A",
    "second_parent_gap_A",
    "nearest_matches_selected",
    "geometrically_unique",
    "bond_class",
    "bond_length_review_pass",
]

with PARENT_MAP_CSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=parent_fieldnames,
    )
    writer.writeheader()
    writer.writerows(
        parent_records
    )

parameter_fieldnames = [
    "interaction_class",
    "parent_element",
    "hydrogen_element",
    "instance_count",
    "observed_length_mean_A",
    "observed_length_std_A",
    "observed_length_min_A",
    "observed_length_max_A",
    "required_new_hydrogen_atom_type",
    "required_nonbonded_parameters",
    "required_bond_parameters",
    "required_angle_parameters",
    "required_improper_review",
    "parameter_adoption_status",
]

with PARAMETER_REQUIREMENTS_CSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=parameter_fieldnames,
    )
    writer.writeheader()
    writer.writerows(
        parameter_requirements
    )


print("\n[9] GATES")

gates = {
    "D040_A3_decision_gate": True,
    "adopted_atom_count_37_gate": (
        len(atoms)
        == EXPECTED_TOTAL_COUNT
    ),
    "retained_hydrogen_count_6_gate": (
        len(hydrogens)
        == EXPECTED_H_COUNT
    ),
    "nominal_edges_loaded_gate": (
        len(edge_rows) > 0
    ),
    "all_H_have_parent_gate": (
        all_H_have_parent_gate
    ),
    "all_H_physical_role_gate": (
        all_H_are_physical_passivants_gate
    ),
    "nearest_parent_agreement_gate": (
        nearest_agreement_count
        == EXPECTED_H_COUNT
    ),
    "bond_length_review_gate": (
        bond_length_pass_count
        == EXPECTED_H_COUNT
    ),
    "parent_map_output_created_gate": (
        PARENT_MAP_CSV.is_file()
        and PARENT_MAP_CSV.stat().st_size > 0
    ),
    "parameter_requirements_output_created_gate": (
        PARAMETER_REQUIREMENTS_CSV.is_file()
        and PARAMETER_REQUIREMENTS_CSV.stat().st_size > 0
    ),
    "no_topology_modified_gate": True,
    "no_coordinates_modified_gate": True,
    "no_atom_type_adopted_gate": True,
    "no_bonded_parameter_adopted_gate": True,
    "force_field_adoption_blocked_gate": True,
    "MD_execution_blocked_gate": True,
}

for gate_name, value in gates.items():
    print(
        f"{gate_name}="
        f"{'PASS' if value else 'FAIL'}"
    )

all_gates_pass = all(
    gates.values()
)


print("\n[10] WRITE REPORT")

decision = (
    "D040_A4_HYDROGEN_AUGMENTATION_DESIGN_PASS_"
    "PARAMETER_SOURCE_AUDIT_AUTHORIZED"
    if all_gates_pass
    else
    "D040_A4_HYDROGEN_AUGMENTATION_DESIGN_"
    "REVIEW_REQUIRED"
)

report = {
    "generated_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "decision": decision,
    "source_identity": {
        "A3_report": {
            "path": str(
                A3_REPORT.relative_to(ROOT)
            ),
            "sha256": sha256(
                A3_REPORT
            ),
        },
        "adopted_charges_csv": {
            "path": str(
                charges_csv.relative_to(ROOT)
            ),
            "sha256": sha256(
                charges_csv
            ),
        },
        "nominal_edges_csv": {
            "path": str(
                NOMINAL_EDGES.relative_to(ROOT)
            ),
            "sha256": sha256(
                NOMINAL_EDGES
            ),
        },
        "provenance_map_csv": {
            "path": str(
                PROVENANCE_MAP.relative_to(ROOT)
            ),
            "sha256": sha256(
                PROVENANCE_MAP
            ),
        },
    },
    "hydrogen_summary": {
        "retained_hydrogen_count": (
            len(hydrogens)
        ),
        "B_H_count": B_H_count,
        "N_H_count": N_H_count,
        "unique_parent_count": (
            unique_parent_count
        ),
        "duplicate_parent_ids": (
            duplicate_parent_ids
        ),
        "authoritative_edge_count": (
            authoritative_edge_count
        ),
        "geometric_fallback_count": (
            fallback_count
        ),
        "nearest_agreement_count": (
            nearest_agreement_count
        ),
        "geometrically_unique_count": (
            geometrically_unique_count
        ),
        "bond_length_pass_count": (
            bond_length_pass_count
        ),
    },
    "parent_records": (
        parent_records
    ),
    "parameter_requirements": (
        parameter_requirements
    ),
    "scientific_interpretation": {
        "all_H_are_physical_passivants": (
            all_H_are_physical_passivants_gate
        ),
        "all_H_have_parent": (
            all_H_have_parent_gate
        ),
        "parent_assignment_scientifically_supported": (
            parent_assignment_scientifically_supported_gate
        ),
        "hydrogen_topology_augmentation_required": True,
        "new_H_atom_type_required": True,
        "B_H_parameters_required": (
            B_H_count > 0
        ),
        "N_H_parameters_required": (
            N_H_count > 0
        ),
    },
    "gates": gates,
    "authorizations": {
        "hydrogen_parent_mapping_review_authorized": (
            all_gates_pass
        ),
        "hydrogen_parameter_source_audit_authorized": (
            all_gates_pass
        ),
        "hydrogen_coordinate_insertion_authorized": False,
        "new_atom_type_definition_authorized": False,
        "bonded_parameter_modification_authorized": False,
        "charge_to_topology_mapping_execution_authorized": False,
        "topology_modification_authorized": False,
        "force_field_adoption_authorized": False,
        "energy_execution_authorized": False,
        "minimization_execution_authorized": False,
        "validation_MD_execution_authorized": False,
        "production_MD_authorized": False,
    },
    "next_required_block": {
        "name": (
            "D040_A5_HYDROGEN_PARAMETER_SOURCE_AUDIT"
        ),
        "required_actions": [
            (
                "Search the accepted parameter libraries and project "
                "history for existing H atom types compatible with "
                "edge-passivated hBN."
            ),
            (
                "Search for B-H and N-H bond, angle and improper "
                "parameters already used or documented in the project."
            ),
            (
                "Compare candidate parameter provenance and units "
                "without adopting any parameter."
            ),
            (
                "Keep coordinate insertion and topology modification "
                "blocked until parameter provenance is resolved."
            ),
        ],
    },
}

REPORT_JSON.write_text(
    json.dumps(
        report,
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)

print(f"parent_map_csv = {PARENT_MAP_CSV}")
print(
    f"parent_map_csv_sha256 = "
    f"{sha256(PARENT_MAP_CSV)}"
)
print(
    f"parameter_requirements_csv = "
    f"{PARAMETER_REQUIREMENTS_CSV}"
)
print(
    "parameter_requirements_csv_sha256 = "
    f"{sha256(PARAMETER_REQUIREMENTS_CSV)}"
)
print(f"report_json = {REPORT_JSON}")
print(
    f"report_json_sha256 = "
    f"{sha256(REPORT_JSON)}"
)


print("\n[11] DECISION")

print(f"decision={decision}")
print(
    "hydrogen_parent_mapping_review_authorized="
    f"{all_gates_pass}"
)
print(
    "hydrogen_parameter_source_audit_authorized="
    f"{all_gates_pass}"
)
print(
    "hydrogen_coordinate_insertion_authorized=False"
)
print(
    "new_atom_type_definition_authorized=False"
)
print(
    "bonded_parameter_modification_authorized=False"
)
print(
    "charge_to_topology_mapping_execution_authorized=False"
)
print(
    "topology_modification_authorized=False"
)
print(
    "force_field_adoption_authorized=False"
)
print(
    "validation_MD_execution_authorized=False"
)
print(
    "production_MD_authorized=False"
)
print("=" * 100)

if not all_gates_pass:
    raise SystemExit(2)
