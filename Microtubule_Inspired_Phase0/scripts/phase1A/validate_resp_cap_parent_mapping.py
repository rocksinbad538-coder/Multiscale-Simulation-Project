#!/usr/bin/env python3
"""
DAY038 / D038-F3

Validate the authoritative artificial-cap to real-parent mapping for
QM_F06_UPPER_V7A_R1 using the adopted nominal-edge table.

No charge redistribution is performed.
RESP Stage 2 remains blocked.
"""

from __future__ import annotations

import csv
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from resp_common import load_json, require_file, sha256


ROOT = Path(__file__).resolve().parents[2]

EXECUTION_PARENT = (
    ROOT
    / "runs/phase1A/day038_resp_stage1_executions"
)

LATEST_POINTER = (
    EXECUTION_PARENT
    / "LATEST_RESP_STAGE1_UPPER_V7A_R1_EXECUTION.txt"
)

NOMINAL_EDGES = (
    ROOT
    / "runs/phase1A/"
      "day035_qm_f06_upper_v7a_r1_coordinate_adoption/"
      "QM_F06_UPPER_V7A_ADOPTED_nominal_edges.csv"
)

ATOM_PROVENANCE = (
    ROOT
    / "runs/phase1A/"
      "day035_qm_f06_upper_v7a_r1_coordinate_adoption/"
      "QM_F06_UPPER_V7A_ADOPTED_atom_role_provenance_map.csv"
)

EXPECTED_UPSTREAM_DECISION = (
    "D038_F2_CAP_PARENT_MAPPING_INVENTORY_PASS_"
    "PARENT_ASSIGNMENTS_REQUIRE_EVIDENCE_REVIEW"
)


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()

    if normalized in {"true", "1", "yes", "y"}:
        return True

    if normalized in {"false", "0", "no", "n", ""}:
        return False

    raise RuntimeError(
        f"Unrecognized Boolean value: {value!r}"
    )


def detect_column(
    fieldnames: list[str],
    candidates: tuple[str, ...],
    label: str,
) -> str:
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate

    raise RuntimeError(
        f"Could not identify {label} column.\n"
        f"Available columns: {fieldnames}"
    )


print("=" * 100)
print("DAY038 / D038-F3 — AUTHORITATIVE CAP-PARENT MAPPING VALIDATION")
print("=" * 100)


print("\n[1] SOURCE FILES")

require_file(LATEST_POINTER)
require_file(NOMINAL_EDGES)
require_file(ATOM_PROVENANCE)

execution_dir = (
    ROOT
    / LATEST_POINTER.read_text(
        encoding="utf-8"
    ).strip()
)

transferability_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1_TRANSFERABILITY.csv"
)

mapping_inventory_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_CAP_PARENT_MAPPING_INVENTORY.json"
)

require_file(transferability_csv)
require_file(mapping_inventory_json)

for path in (
    transferability_csv,
    mapping_inventory_json,
    NOMINAL_EDGES,
    ATOM_PROVENANCE,
):
    print(f"FOUND  {path}")


print("\n[2] UPSTREAM AUTHORIZATION")

inventory = load_json(
    mapping_inventory_json
)

if inventory.get("decision") != EXPECTED_UPSTREAM_DECISION:
    raise RuntimeError(
        "Unexpected cap-mapping inventory decision.\n"
        f"Observed: {inventory.get('decision')}"
    )

authorizations = inventory.get(
    "authorizations",
    {},
)

if (
    authorizations.get(
        "cap_parent_mapping_review_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Cap-parent mapping review is not authorized"
    )

if (
    authorizations.get(
        "cap_charge_redistribution_authorized"
    )
    is not False
):
    raise RuntimeError(
        "Unexpected redistribution authorization"
    )

print("upstream_decision_gate            = PASS")
print("mapping_review_authorized_gate    = PASS")
print("redistribution_blocked_gate       = PASS")


print("\n[3] LOAD RESP ATOM TABLE")

with transferability_csv.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    atom_rows = list(
        csv.DictReader(handle)
    )

if len(atom_rows) != 52:
    raise RuntimeError(
        f"Expected 52 RESP atoms, observed {len(atom_rows)}"
    )

atom_by_id: dict[str, dict] = {}
atom_by_index: dict[int, dict] = {}

for row in atom_rows:
    index = int(row["atom_index_0based"])
    atom_id = row["atom_id"].strip()

    row["atom_index_0based_int"] = index
    row["artificial_cap_bool"] = parse_bool(
        row["artificial_cap"]
    )
    row["RESP_stage1_charge_e_float"] = float(
        row["RESP_stage1_charge_e"]
    )

    if atom_id in atom_by_id:
        raise RuntimeError(
            f"Duplicate atom_id: {atom_id}"
        )

    atom_by_id[atom_id] = row
    atom_by_index[index] = row

caps = [
    row
    for row in atom_rows
    if row["artificial_cap_bool"]
]

real_atoms = [
    row
    for row in atom_rows
    if not row["artificial_cap_bool"]
]

if len(caps) != 15 or len(real_atoms) != 37:
    raise RuntimeError(
        "Unexpected cap/real-atom counts"
    )

print(f"atom_count           = {len(atom_rows)}")
print(f"artificial_cap_count = {len(caps)}")
print(f"real_atom_count       = {len(real_atoms)}")
print("RESP_atom_table_gate = PASS")


print("\n[4] LOAD NOMINAL EDGES")

with NOMINAL_EDGES.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    edge_reader = csv.DictReader(handle)
    edge_rows = list(edge_reader)
    edge_fields = edge_reader.fieldnames or []

print(f"nominal_edge_columns = {edge_fields}")
print(f"nominal_edge_count   = {len(edge_rows)}")

atom1_id_col = detect_column(
    edge_fields,
    (
        "atom_1_id",
        "atom1_id",
        "source_atom_id",
        "node_1_id",
        "atom_i_id",
        "first_atom",
    ),
    "first atom ID",
)

atom1_element_col = detect_column(
    edge_fields,
    (
        "atom_1_element",
        "atom1_element",
        "source_element",
        "node_1_element",
        "atom_i_element",
        "first_element",
    ),
    "first atom element",
)

atom2_id_col = detect_column(
    edge_fields,
    (
        "atom_2_id",
        "atom2_id",
        "target_atom_id",
        "node_2_id",
        "atom_j_id",
        "second_atom",
    ),
    "second atom ID",
)

atom2_element_col = detect_column(
    edge_fields,
    (
        "atom_2_element",
        "atom2_element",
        "target_element",
        "node_2_element",
        "atom_j_element",
        "second_element",
    ),
    "second atom element",
)

distance_col = detect_column(
    edge_fields,
    (
        "distance_A",
        "distance_angstrom",
        "bond_length_A",
        "distance",
    ),
    "edge distance",
)

print(f"atom1_id_column      = {atom1_id_col}")
print(f"atom1_element_column = {atom1_element_col}")
print(f"atom2_id_column      = {atom2_id_col}")
print(f"atom2_element_column = {atom2_element_col}")
print(f"distance_column      = {distance_col}")


print("\n[5] RESOLVE CAP-PARENT EDGES")

cap_ids = {
    row["atom_id"].strip()
    for row in caps
}

cap_edges: dict[str, list[dict]] = {
    cap_id: []
    for cap_id in cap_ids
}

for edge in edge_rows:
    atom1_id = edge[atom1_id_col].strip()
    atom2_id = edge[atom2_id_col].strip()

    if atom1_id in cap_ids:
        cap_edges[atom1_id].append(
            {
                "cap_side": 1,
                "parent_atom_id": atom2_id,
                "parent_element": edge[
                    atom2_element_col
                ].strip(),
                "cap_element": edge[
                    atom1_element_col
                ].strip(),
                "distance_A": float(
                    edge[distance_col]
                ),
                "raw_edge": edge,
            }
        )

    if atom2_id in cap_ids:
        cap_edges[atom2_id].append(
            {
                "cap_side": 2,
                "parent_atom_id": atom1_id,
                "parent_element": edge[
                    atom1_element_col
                ].strip(),
                "cap_element": edge[
                    atom2_element_col
                ].strip(),
                "distance_A": float(
                    edge[distance_col]
                ),
                "raw_edge": edge,
            }
        )


mapping_records = []

for cap in sorted(
    caps,
    key=lambda row: row["atom_index_0based_int"],
):
    cap_id = cap["atom_id"].strip()
    matches = cap_edges[cap_id]

    if len(matches) != 1:
        mapping_status = (
            "UNRESOLVED_NO_EDGE"
            if len(matches) == 0
            else "AMBIGUOUS_MULTIPLE_EDGES"
        )

        parent = None
    else:
        mapping_status = "VERIFIED_SINGLE_NOMINAL_EDGE"
        parent = matches[0]

    if parent is not None:
        parent_id = parent["parent_atom_id"]

        if parent_id not in atom_by_id:
            raise RuntimeError(
                "Nominal-edge parent is not present in the "
                f"52-atom RESP table: {parent_id}"
            )

        parent_row = atom_by_id[parent_id]

        if parent_row["artificial_cap_bool"]:
            raise RuntimeError(
                f"Cap {cap_id} maps to another cap: {parent_id}"
            )

        parent_index = parent_row[
            "atom_index_0based_int"
        ]

        parent_element = parent_row["element"].strip()

        if (
            parent_element
            != parent["parent_element"]
        ):
            raise RuntimeError(
                "Parent-element mismatch.\n"
                f"Cap: {cap_id}\n"
                f"Parent: {parent_id}\n"
                f"Edge element: {parent['parent_element']}\n"
                f"RESP element: {parent_element}"
            )
    else:
        parent_id = None
        parent_index = None
        parent_element = None
        parent_row = None

    record = {
        "cap_atom_index_0based": (
            cap["atom_index_0based_int"]
        ),
        "cap_atom_index_1based": int(
            cap["atom_index_1based"]
        ),
        "cap_atom_id": cap_id,
        "cap_element": cap["element"].strip(),
        "cap_role": cap["atom_role"].strip(),
        "cap_RESP_stage1_charge_e": (
            cap["RESP_stage1_charge_e_float"]
        ),
        "nominal_edge_match_count": len(matches),
        "mapping_status": mapping_status,
        "parent_atom_index_0based": parent_index,
        "parent_atom_index_1based": (
            parent_index + 1
            if parent_index is not None
            else None
        ),
        "parent_atom_id": parent_id,
        "parent_element": parent_element,
        "parent_atom_role": (
            parent_row["atom_role"].strip()
            if parent_row is not None
            else None
        ),
        "parent_RESP_stage1_charge_e": (
            parent_row[
                "RESP_stage1_charge_e_float"
            ]
            if parent_row is not None
            else None
        ),
        "cap_parent_distance_A": (
            parent["distance_A"]
            if parent is not None
            else None
        ),
        "charge_redistribution_authorized": False,
    }

    mapping_records.append(record)

    print(
        f"cap={record['cap_atom_index_0based']:>2} "
        f"id={cap_id:<43} "
        f"q_cap={record['cap_RESP_stage1_charge_e']: .6f} "
        f"parent_index={str(parent_index):>4} "
        f"parent_id={str(parent_id):<30} "
        f"element={str(parent_element):<2} "
        f"distance_A={record['cap_parent_distance_A']} "
        f"status={mapping_status}"
    )


print("\n[6] MAPPING SUMMARY")

verified_records = [
    record
    for record in mapping_records
    if record["mapping_status"]
    == "VERIFIED_SINGLE_NOMINAL_EDGE"
]

unresolved_records = [
    record
    for record in mapping_records
    if record["mapping_status"]
    != "VERIFIED_SINGLE_NOMINAL_EDGE"
]

parent_counts = Counter(
    record["parent_atom_id"]
    for record in verified_records
)

parents_with_multiple_caps = {
    parent_id: count
    for parent_id, count in parent_counts.items()
    if count > 1
}

cap_charge_sum = sum(
    record["cap_RESP_stage1_charge_e"]
    for record in mapping_records
)

print(
    f"verified_mapping_count = "
    f"{len(verified_records)}/{len(mapping_records)}"
)
print(
    f"unresolved_mapping_count = "
    f"{len(unresolved_records)}"
)
print(
    f"unique_parent_atom_count = "
    f"{len(parent_counts)}"
)
print(
    f"parents_with_multiple_caps = "
    f"{parents_with_multiple_caps}"
)
print(
    f"mapped_cap_charge_sum_e = "
    f"{cap_charge_sum:.16g}"
)


print("\n[7] PARENT-LEVEL AGGREGATION")

parent_summary = []

for parent_id in sorted(parent_counts):
    selected = [
        record
        for record in verified_records
        if record["parent_atom_id"] == parent_id
    ]

    cap_sum = sum(
        record["cap_RESP_stage1_charge_e"]
        for record in selected
    )

    parent_row = atom_by_id[parent_id]
    parent_initial_charge = parent_row[
        "RESP_stage1_charge_e_float"
    ]

    summary = {
        "parent_atom_id": parent_id,
        "parent_atom_index_0based": parent_row[
            "atom_index_0based_int"
        ],
        "parent_atom_index_1based": (
            parent_row["atom_index_0based_int"] + 1
        ),
        "parent_element": parent_row[
            "element"
        ].strip(),
        "parent_atom_role": parent_row[
            "atom_role"
        ].strip(),
        "parent_RESP_stage1_charge_e": (
            parent_initial_charge
        ),
        "attached_cap_count": len(selected),
        "attached_cap_atom_indices_0based": [
            record["cap_atom_index_0based"]
            for record in selected
        ],
        "attached_cap_ids": [
            record["cap_atom_id"]
            for record in selected
        ],
        "attached_cap_charge_sum_e": cap_sum,
        "candidate_parent_plus_cap_charge_e": (
            parent_initial_charge + cap_sum
        ),
        "redistribution_authorized": False,
    }

    parent_summary.append(summary)

    print(
        f"parent_index={summary['parent_atom_index_0based']:>2} "
        f"parent_id={parent_id:<30} "
        f"element={summary['parent_element']} "
        f"q_parent={parent_initial_charge: .6f} "
        f"caps={summary['attached_cap_count']} "
        f"q_caps={cap_sum: .6f} "
        f"q_parent_plus_caps="
        f"{summary['candidate_parent_plus_cap_charge_e']: .6f}"
    )


print("\n[8] SCIENTIFIC GATES")

all_distances_finite_gate = all(
    record["cap_parent_distance_A"] is not None
    and math.isfinite(
        record["cap_parent_distance_A"]
    )
    and record["cap_parent_distance_A"] > 0.0
    for record in verified_records
)

all_parents_real_gate = all(
    not atom_by_id[
        record["parent_atom_id"]
    ]["artificial_cap_bool"]
    for record in verified_records
)

mapping_complete_gate = (
    len(verified_records) == 15
    and len(unresolved_records) == 0
)

charge_sum_consistency_gate = (
    abs(cap_charge_sum - 0.54585)
    <= 5.0e-6
)

gates = {
    "upstream_decision_gate": True,
    "mapping_review_authorized_gate": True,
    "atom_count_gate": len(atom_rows) == 52,
    "cap_count_gate": len(caps) == 15,
    "nominal_edge_table_gate": len(edge_rows) > 0,
    "mapping_complete_gate": mapping_complete_gate,
    "all_parents_present_gate": all(
        record["parent_atom_id"] in atom_by_id
        for record in verified_records
    ),
    "all_parents_real_gate": all_parents_real_gate,
    "all_distances_finite_gate": (
        all_distances_finite_gate
    ),
    "charge_sum_consistency_gate": (
        charge_sum_consistency_gate
    ),
    "no_redistribution_performed_gate": all(
        record[
            "charge_redistribution_authorized"
        ]
        is False
        for record in mapping_records
    ),
}

for name, value in gates.items():
    print(
        f"{name} = "
        f"{'PASS' if value else 'FAIL'}"
    )

all_gates_pass = all(gates.values())


print("\n[9] WRITE OUTPUTS")

mapping_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_CAP_PARENT_MAPPING_VALIDATED.csv"
)

mapping_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_CAP_PARENT_MAPPING_VALIDATED.json"
)

parent_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_CAP_PARENT_AGGREGATION.csv"
)

mapping_fieldnames = [
    "cap_atom_index_0based",
    "cap_atom_index_1based",
    "cap_atom_id",
    "cap_element",
    "cap_role",
    "cap_RESP_stage1_charge_e",
    "nominal_edge_match_count",
    "mapping_status",
    "parent_atom_index_0based",
    "parent_atom_index_1based",
    "parent_atom_id",
    "parent_element",
    "parent_atom_role",
    "parent_RESP_stage1_charge_e",
    "cap_parent_distance_A",
    "charge_redistribution_authorized",
]

with mapping_csv.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=mapping_fieldnames,
    )
    writer.writeheader()
    writer.writerows(mapping_records)

parent_fieldnames = [
    "parent_atom_id",
    "parent_atom_index_0based",
    "parent_atom_index_1based",
    "parent_element",
    "parent_atom_role",
    "parent_RESP_stage1_charge_e",
    "attached_cap_count",
    "attached_cap_atom_indices_0based",
    "attached_cap_ids",
    "attached_cap_charge_sum_e",
    "candidate_parent_plus_cap_charge_e",
    "redistribution_authorized",
]

with parent_csv.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=parent_fieldnames,
    )
    writer.writeheader()

    for row in parent_summary:
        csv_row = dict(row)
        csv_row[
            "attached_cap_atom_indices_0based"
        ] = json.dumps(
            csv_row[
                "attached_cap_atom_indices_0based"
            ]
        )
        csv_row["attached_cap_ids"] = json.dumps(
            csv_row["attached_cap_ids"]
        )
        writer.writerow(csv_row)

decision = (
    "D038_F3_CAP_PARENT_MAPPING_VALIDATION_PASS_"
    "REDISTRIBUTION_DESIGN_REVIEW_AUTHORIZED"
    if all_gates_pass
    else
    "D038_F3_CAP_PARENT_MAPPING_VALIDATION_"
    "REVIEW_REQUIRED"
)

report = {
    "generated_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "decision": decision,
    "source_execution_directory": str(
        execution_dir.resolve()
    ),
    "source_identity": {
        "mapping_inventory_json": str(
            mapping_inventory_json.resolve()
        ),
        "mapping_inventory_json_sha256": sha256(
            mapping_inventory_json
        ),
        "nominal_edges": str(
            NOMINAL_EDGES.resolve()
        ),
        "nominal_edges_sha256": sha256(
            NOMINAL_EDGES
        ),
        "atom_provenance": str(
            ATOM_PROVENANCE.resolve()
        ),
        "atom_provenance_sha256": sha256(
            ATOM_PROVENANCE
        ),
        "transferability_csv": str(
            transferability_csv.resolve()
        ),
        "transferability_csv_sha256": sha256(
            transferability_csv
        ),
    },
    "summary": {
        "artificial_cap_count": len(caps),
        "verified_mapping_count": len(
            verified_records
        ),
        "unresolved_mapping_count": len(
            unresolved_records
        ),
        "unique_parent_atom_count": len(
            parent_counts
        ),
        "mapped_cap_charge_sum_e": (
            cap_charge_sum
        ),
        "parents_with_multiple_caps": (
            parents_with_multiple_caps
        ),
    },
    "cap_parent_mappings": mapping_records,
    "parent_aggregation": parent_summary,
    "gates": gates,
    "authorizations": {
        "cap_charge_redistribution_design_review_authorized": (
            all_gates_pass
        ),
        "cap_charge_redistribution_execution_authorized": False,
        "RESP_stage2_execution_authorized": False,
        "charge_adoption_authorized": False,
        "force_field_adoption_authorized": False,
    },
    "outputs": {
        "mapping_csv": str(
            mapping_csv.resolve()
        ),
        "parent_aggregation_csv": str(
            parent_csv.resolve()
        ),
    },
}

mapping_json.write_text(
    json.dumps(
        report,
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)

print(f"mapping_csv = {mapping_csv}")
print(
    f"mapping_csv_sha256 = "
    f"{sha256(mapping_csv)}"
)
print(f"parent_csv = {parent_csv}")
print(
    f"parent_csv_sha256 = "
    f"{sha256(parent_csv)}"
)
print(f"mapping_json = {mapping_json}")
print(
    f"mapping_json_sha256 = "
    f"{sha256(mapping_json)}"
)


print("\n[10] DECISION")

print(f"decision = {decision}")
print(
    "cap_charge_redistribution_design_review_authorized = "
    f"{all_gates_pass}"
)
print(
    "cap_charge_redistribution_execution_authorized = False"
)
print("RESP_stage2_execution_authorized = False")
print("charge_adoption_authorized = False")
print("force_field_adoption_authorized = False")
print("=" * 100)
