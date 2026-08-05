#!/usr/bin/env python3
"""
DAY039 / D039-A10

Focused chemical and topological audit of lambda=4, the only candidate
passing the D039-A9 project-level review filters.

Comparison lambdas:
    3, 4 and 5

No lambda or charge set is adopted.
RESP Stage 2 remains blocked.
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from resp_common import load_json, require_file, sha256


ROOT = Path(__file__).resolve().parents[2]

LATEST_POINTER = (
    ROOT
    / "runs/phase1A/day038_resp_stage1_executions"
    / "LATEST_RESP_STAGE1_UPPER_V7A_R1_EXECUTION.txt"
)

NOMINAL_EDGES = (
    ROOT
    / "runs/phase1A/day035_qm_f06_upper_v7a_r1_coordinate_adoption"
    / "QM_F06_UPPER_V7A_ADOPTED_nominal_edges.csv"
)

EXPECTED_A9_DECISION = (
    "D039_A9_LAMBDA_1_TO_10_REFINEMENT_PASS_"
    "FOCUSED_CANDIDATE_REVIEW_AUTHORIZED"
)

FOCUSED_LAMBDA = 4.0
COMPARISON_LAMBDAS = (3.0, 4.0, 5.0)


def sign_changed(
    initial: float,
    candidate: float,
) -> bool:
    return (
        initial != 0.0
        and candidate != 0.0
        and math.copysign(1.0, initial)
        != math.copysign(1.0, candidate)
    )


def charge_column(value: float) -> str:
    return f"lambda_{value:g}_charge_e"


print("=" * 100)
print("DAY039 / D039-A10 — LAMBDA=4 FOCUSED CHEMICAL AUDIT")
print("=" * 100)


print("\n[1] SOURCE EXECUTION")

require_file(LATEST_POINTER)
require_file(NOMINAL_EDGES)

execution_dir = (
    ROOT
    / LATEST_POINTER.read_text(
        encoding="utf-8"
    ).strip()
)

a9_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA_1_TO_10_REFINEMENT.json"
)

a9_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA_1_TO_10_REFINEMENT.csv"
)

for path in (
    a9_json,
    a9_csv,
    NOMINAL_EDGES,
):
    require_file(path)
    print(f"FOUND  {path}")


print("\n[2] UPSTREAM AUTHORIZATION")

a9_report = load_json(a9_json)

if (
    a9_report.get("decision")
    != EXPECTED_A9_DECISION
):
    raise RuntimeError(
        "Unexpected A9 decision.\n"
        f"Observed: {a9_report.get('decision')}"
    )

admissible_lambdas = [
    float(value)
    for value in a9_report.get(
        "focused_review_admissible_lambdas",
        [],
    )
]

if admissible_lambdas != [FOCUSED_LAMBDA]:
    raise RuntimeError(
        "Expected lambda=4 to be the unique focused-review candidate.\n"
        f"Observed: {admissible_lambdas}"
    )

authorizations = a9_report.get(
    "authorizations",
    {},
)

if (
    authorizations.get(
        "focused_lambda_candidate_review_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Focused lambda review is not authorized"
    )

if (
    authorizations.get(
        "regularization_lambda_adoption_authorized"
    )
    is not False
):
    raise RuntimeError(
        "Unexpected lambda-adoption authorization"
    )

print("A9_decision_gate                    = PASS")
print("unique_lambda4_candidate_gate       = PASS")
print("focused_review_authorized_gate      = PASS")
print("lambda_adoption_blocked_gate        = PASS")


print("\n[3] LOAD REFINED CHARGE TABLE")

with a9_csv.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    rows = list(
        csv.DictReader(handle)
    )

if len(rows) != 37:
    raise RuntimeError(
        f"Expected 37 real atoms, observed {len(rows)}"
    )

required_columns = {
    "real_atom_sequence_index",
    "original_atom_index_0based",
    "atom_id",
    "element",
    "atom_role",
    "RESP_stage1_charge_e",
}

for value in COMPARISON_LAMBDAS:
    required_columns.add(
        charge_column(value)
    )
    required_columns.add(
        f"lambda_{value:g}_delta_e"
    )

missing_columns = (
    required_columns
    - set(rows[0])
)

if missing_columns:
    raise RuntimeError(
        "Missing required columns:\n"
        + "\n".join(
            sorted(missing_columns)
        )
    )

atom_by_id = {}

for row in rows:
    row["real_atom_sequence_index_int"] = int(
        row["real_atom_sequence_index"]
    )

    row["original_atom_index_0based_int"] = int(
        row["original_atom_index_0based"]
    )

    row["RESP_stage1_charge_e_float"] = float(
        row["RESP_stage1_charge_e"]
    )

    for value in COMPARISON_LAMBDAS:
        row[
            f"lambda_{value}_charge_e"
        ] = float(
            row[
                charge_column(value)
            ]
        )

    atom_by_id[
        row["atom_id"].strip()
    ] = row

print(f"real_atom_count = {len(rows)}")
print("candidate_column_contract_gate = PASS")


print("\n[4] LOAD REAL-ATOM TOPOLOGY")

with NOMINAL_EDGES.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    edge_rows = list(
        csv.DictReader(handle)
    )

neighbors: dict[str, list[dict]] = defaultdict(list)

for edge in edge_rows:
    first = edge["first_atom"].strip()
    second = edge["second_atom"].strip()

    if (
        first in atom_by_id
        and second in atom_by_id
    ):
        distance = float(
            edge["distance_A"]
        )

        neighbors[first].append(
            {
                "atom_id": second,
                "distance_A": distance,
                "edge_type": edge["edge_type"],
                "provenance": edge["provenance"],
            }
        )

        neighbors[second].append(
            {
                "atom_id": first,
                "distance_A": distance,
                "edge_type": edge["edge_type"],
                "provenance": edge["provenance"],
            }
        )

real_real_edge_count = sum(
    len(value)
    for value in neighbors.values()
) // 2

print(f"nominal_edge_count   = {len(edge_rows)}")
print(f"real_real_edge_count = {real_real_edge_count}")


print("\n[5] SIGN-CHANGE EVOLUTION")

sign_change_records = {
    value: []
    for value in COMPARISON_LAMBDAS
}

for value in COMPARISON_LAMBDAS:
    print(f"\nlambda={value:g}")

    for row in rows:
        initial = row[
            "RESP_stage1_charge_e_float"
        ]

        candidate = row[
            f"lambda_{value}_charge_e"
        ]

        if sign_changed(
            initial,
            candidate,
        ):
            record = {
                "real_atom_sequence_index": (
                    row[
                        "real_atom_sequence_index_int"
                    ]
                ),
                "original_atom_index_0based": (
                    row[
                        "original_atom_index_0based_int"
                    ]
                ),
                "atom_id": row["atom_id"],
                "element": row["element"],
                "atom_role": row["atom_role"],
                "RESP_stage1_charge_e": initial,
                "candidate_charge_e": candidate,
                "delta_e": (
                    candidate - initial
                ),
            }

            sign_change_records[
                value
            ].append(record)

            print(
                f"  atom_id={record['atom_id']:<30} "
                f"element={record['element']} "
                f"role={record['atom_role']} "
                f"RESP1={initial: .9f} "
                f"candidate={candidate: .9f} "
                f"delta={candidate-initial: .9f}"
            )

    print(
        f"  total_sign_changes="
        f"{len(sign_change_records[value])}"
    )


print("\n[6] LAMBDA=4 PERSISTENT SIGN-CHANGE ATOMS")

lambda4_changes = sign_change_records[
    FOCUSED_LAMBDA
]

if len(lambda4_changes) != 2:
    raise RuntimeError(
        "Expected exactly two lambda=4 sign changes"
    )

persistent_ids = {
    record["atom_id"]
    for record in lambda4_changes
}

for record in lambda4_changes:
    atom_id = record["atom_id"]
    row = atom_by_id[atom_id]

    print(
        f"\natom_id={atom_id} "
        f"element={row['element']} "
        f"role={row['atom_role']} "
        f"original_index="
        f"{row['original_atom_index_0based_int']}"
    )

    print(
        f"  RESP1="
        f"{row['RESP_stage1_charge_e_float']: .9f}"
    )

    for value in COMPARISON_LAMBDAS:
        candidate = row[
            f"lambda_{value}_charge_e"
        ]

        print(
            f"  lambda={value:g} "
            f"charge={candidate: .9f} "
            f"delta="
            f"{candidate-row['RESP_stage1_charge_e_float']: .9f}"
        )

    print(
        f"  real_neighbor_count="
        f"{len(neighbors[atom_id])}"
    )

    for neighbor in sorted(
        neighbors[atom_id],
        key=lambda item: (
            item["distance_A"],
            item["atom_id"],
        ),
    ):
        neighbor_row = atom_by_id[
            neighbor["atom_id"]
        ]

        print(
            f"    neighbor={neighbor['atom_id']} "
            f"element={neighbor_row['element']} "
            f"role={neighbor_row['atom_role']} "
            f"distance_A={neighbor['distance_A']:.6f} "
            f"RESP1="
            f"{neighbor_row['RESP_stage1_charge_e_float']: .6f} "
            f"lambda3="
            f"{neighbor_row['lambda_3.0_charge_e']: .6f} "
            f"lambda4="
            f"{neighbor_row['lambda_4.0_charge_e']: .6f} "
            f"lambda5="
            f"{neighbor_row['lambda_5.0_charge_e']: .6f}"
        )


print("\n[7] ZERO-CROSSING BRACKETS")

zero_crossing_records = []

for row in rows:
    initial = row[
        "RESP_stage1_charge_e_float"
    ]

    q3 = row["lambda_3.0_charge_e"]
    q4 = row["lambda_4.0_charge_e"]
    q5 = row["lambda_5.0_charge_e"]

    crossings = []

    if q3 == 0.0 or q4 == 0.0 or q3 * q4 < 0.0:
        crossings.append(
            {
                "interval": [3.0, 4.0],
                "q_left": q3,
                "q_right": q4,
            }
        )

    if q4 == 0.0 or q5 == 0.0 or q4 * q5 < 0.0:
        crossings.append(
            {
                "interval": [4.0, 5.0],
                "q_left": q4,
                "q_right": q5,
            }
        )

    if crossings:
        record = {
            "atom_id": row["atom_id"],
            "element": row["element"],
            "atom_role": row["atom_role"],
            "RESP_stage1_charge_e": initial,
            "crossings": crossings,
        }

        zero_crossing_records.append(
            record
        )

        print(
            f"atom_id={row['atom_id']:<30} "
            f"element={row['element']} "
            f"role={row['atom_role']} "
            f"crossings={crossings}"
        )

if not zero_crossing_records:
    print("(none between lambda=3 and lambda=5)")


print("\n[8] ELEMENT-LEVEL COMPARISON")

element_records = []

for element in sorted(
    set(
        row["element"]
        for row in rows
    )
):
    selected = [
        row
        for row in rows
        if row["element"] == element
    ]

    initial_sum = sum(
        row[
            "RESP_stage1_charge_e_float"
        ]
        for row in selected
    )

    print(
        f"\nelement={element} "
        f"count={len(selected)} "
        f"RESP1_sum={initial_sum:.9f}"
    )

    element_record = {
        "element": element,
        "count": len(selected),
        "RESP_stage1_sum_e": initial_sum,
    }

    for value in COMPARISON_LAMBDAS:
        candidates = [
            row[
                f"lambda_{value}_charge_e"
            ]
            for row in selected
        ]

        candidate_sum = sum(
            candidates
        )

        deltas = [
            candidate
            - row[
                "RESP_stage1_charge_e_float"
            ]
            for candidate, row in zip(
                candidates,
                selected,
            )
        ]

        changes = sum(
            sign_changed(
                row[
                    "RESP_stage1_charge_e_float"
                ],
                candidate,
            )
            for row, candidate in zip(
                selected,
                candidates,
            )
        )

        print(
            f"  lambda={value:g} "
            f"sum={candidate_sum:.9f} "
            f"delta_sum="
            f"{candidate_sum-initial_sum:.9f} "
            f"delta_RMS="
            f"{math.sqrt(sum(d*d for d in deltas)/len(deltas)):.9f} "
            f"sign_changes={changes}"
        )

        element_record[
            f"lambda_{value}_sum_e"
        ] = candidate_sum

        element_record[
            f"lambda_{value}_delta_sum_e"
        ] = (
            candidate_sum
            - initial_sum
        )

        element_record[
            f"lambda_{value}_delta_RMS_e"
        ] = math.sqrt(
            sum(
                delta * delta
                for delta in deltas
            )
            / len(deltas)
        )

        element_record[
            f"lambda_{value}_sign_change_count"
        ] = changes

    element_records.append(
        element_record
    )


print("\n[9] ROLE-LEVEL COMPARISON")

role_records = []

for role in sorted(
    set(
        row["atom_role"]
        for row in rows
    )
):
    selected = [
        row
        for row in rows
        if row["atom_role"] == role
    ]

    print(
        f"\nrole={role!r} "
        f"count={len(selected)}"
    )

    role_record = {
        "atom_role": role,
        "count": len(selected),
    }

    for value in COMPARISON_LAMBDAS:
        deltas = [
            row[
                f"lambda_{value}_charge_e"
            ]
            - row[
                "RESP_stage1_charge_e_float"
            ]
            for row in selected
        ]

        delta_rms = math.sqrt(
            sum(
                delta * delta
                for delta in deltas
            )
            / len(deltas)
        )

        delta_max = max(
            abs(delta)
            for delta in deltas
        )

        changes = sum(
            sign_changed(
                row[
                    "RESP_stage1_charge_e_float"
                ],
                row[
                    f"lambda_{value}_charge_e"
                ],
            )
            for row in selected
        )

        print(
            f"  lambda={value:g} "
            f"delta_RMS={delta_rms:.9f} "
            f"delta_max={delta_max:.9f} "
            f"sign_changes={changes}"
        )

        role_record[
            f"lambda_{value}_delta_RMS_e"
        ] = delta_rms

        role_record[
            f"lambda_{value}_delta_max_abs_e"
        ] = delta_max

        role_record[
            f"lambda_{value}_sign_change_count"
        ] = changes

    role_records.append(
        role_record
    )


print("\n[10] A9 ELECTROSTATIC CONTEXT")

a9_records = a9_report.get(
    "records",
    []
)

metrics_by_lambda = {
    float(record["regularization_lambda"]): record
    for record in a9_records
    if float(
        record["regularization_lambda"]
    ) in COMPARISON_LAMBDAS
}

if set(metrics_by_lambda) != set(
    COMPARISON_LAMBDAS
):
    raise RuntimeError(
        "Could not recover lambda=3, 4 and 5 metrics"
    )

for value in COMPARISON_LAMBDAS:
    record = metrics_by_lambda[value]

    print(
        f"lambda={value:g} "
        f"RMSE="
        f"{record['electrostatic']['RMSE_au']:.16g} "
        f"R="
        f"{record['electrostatic']['pearson_r']:.16g} "
        f"same_sign="
        f"{record['electrostatic']['same_sign_fraction']:.16g} "
        f"delta_RMS="
        f"{record['charges']['delta_RMS_e']:.16g} "
        f"delta_max="
        f"{record['charges']['delta_max_abs_e']:.16g} "
        f"sign_changes="
        f"{record['charges']['sign_change_count']}"
    )


print("\n[11] WRITE OUTPUTS")

output_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA4_FOCUSED_CHEMICAL_AUDIT.csv"
)

output_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA4_FOCUSED_CHEMICAL_AUDIT.json"
)

fieldnames = [
    "real_atom_sequence_index",
    "original_atom_index_0based",
    "atom_id",
    "element",
    "atom_role",
    "RESP_stage1_charge_e",
    "lambda_3_charge_e",
    "lambda_3_delta_e",
    "lambda_3_sign_changed",
    "lambda_4_charge_e",
    "lambda_4_delta_e",
    "lambda_4_sign_changed",
    "lambda_5_charge_e",
    "lambda_5_delta_e",
    "lambda_5_sign_changed",
    "real_neighbor_count",
    "real_neighbor_ids",
]

with output_csv.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
    )

    writer.writeheader()

    for row in rows:
        initial = row[
            "RESP_stage1_charge_e_float"
        ]

        neighbor_ids = sorted(
            neighbor["atom_id"]
            for neighbor in neighbors[
                row["atom_id"]
            ]
        )

        output = {
            "real_atom_sequence_index": (
                row[
                    "real_atom_sequence_index_int"
                ]
            ),
            "original_atom_index_0based": (
                row[
                    "original_atom_index_0based_int"
                ]
            ),
            "atom_id": row["atom_id"],
            "element": row["element"],
            "atom_role": row["atom_role"],
            "RESP_stage1_charge_e": initial,
            "real_neighbor_count": len(
                neighbor_ids
            ),
            "real_neighbor_ids": json.dumps(
                neighbor_ids
            ),
        }

        for value in COMPARISON_LAMBDAS:
            candidate = row[
                f"lambda_{value}_charge_e"
            ]

            output[
                f"lambda_{value:g}_charge_e"
            ] = candidate

            output[
                f"lambda_{value:g}_delta_e"
            ] = (
                candidate - initial
            )

            output[
                f"lambda_{value:g}_sign_changed"
            ] = sign_changed(
                initial,
                candidate,
            )

        writer.writerow(output)


print("\n[12] SCIENTIFIC GATES")

finite_gate = all(
    math.isfinite(
        row[
            f"lambda_{value}_charge_e"
        ]
    )
    for row in rows
    for value in COMPARISON_LAMBDAS
)

neutrality_gate = all(
    abs(
        sum(
            row[
                f"lambda_{value}_charge_e"
            ]
            for row in rows
        )
    )
    <= 1.0e-10
    for value in COMPARISON_LAMBDAS
)

lambda4_sign_change_gate = (
    len(lambda4_changes) == 2
)

expected_persistent_ids_gate = (
    persistent_ids
    == {
        "A:UPPER:8:4",
        "P:1583",
    }
)

topology_gate = (
    all(
        record["atom_id"]
        in neighbors
        for record in lambda4_changes
    )
)

gates = {
    "upstream_decision_gate": True,
    "unique_lambda4_candidate_gate": True,
    "real_atom_count_gate": (
        len(rows) == 37
    ),
    "candidate_column_contract_gate": True,
    "finite_candidate_gate": finite_gate,
    "neutrality_gate": neutrality_gate,
    "lambda4_two_sign_changes_gate": (
        lambda4_sign_change_gate
    ),
    "expected_persistent_sign_change_ids_gate": (
        expected_persistent_ids_gate
    ),
    "persistent_atom_topology_gate": (
        topology_gate
    ),
    "output_csv_created_gate": (
        output_csv.is_file()
        and output_csv.stat().st_size > 0
    ),
    "no_lambda_adopted_gate": True,
    "no_charge_adopted_gate": True,
    "RESP_stage2_not_executed_gate": True,
}

for name, value in gates.items():
    print(
        f"{name}="
        f"{'PASS' if value else 'FAIL'}"
    )

all_gates_pass = all(
    gates.values()
)


print("\n[13] WRITE JSON REPORT")

decision = (
    "D039_A10_LAMBDA4_FOCUSED_CHEMICAL_AUDIT_PASS_"
    "FINAL_METHOD_REVIEW_AUTHORIZED"
    if all_gates_pass
    else
    "D039_A10_LAMBDA4_FOCUSED_CHEMICAL_AUDIT_"
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
        "A9_json": str(
            a9_json.resolve()
        ),
        "A9_json_sha256": sha256(
            a9_json
        ),
        "A9_csv": str(
            a9_csv.resolve()
        ),
        "A9_csv_sha256": sha256(
            a9_csv
        ),
        "nominal_edges": str(
            NOMINAL_EDGES.resolve()
        ),
        "nominal_edges_sha256": sha256(
            NOMINAL_EDGES
        ),
    },
    "focused_lambda": FOCUSED_LAMBDA,
    "comparison_lambdas": list(
        COMPARISON_LAMBDAS
    ),
    "lambda4_sign_change_records": (
        lambda4_changes
    ),
    "zero_crossing_records": (
        zero_crossing_records
    ),
    "element_records": element_records,
    "role_records": role_records,
    "electrostatic_context": {
        str(value): metrics_by_lambda[value]
        for value in COMPARISON_LAMBDAS
    },
    "gates": gates,
    "authorizations": {
        "lambda4_final_method_scientific_review_authorized": (
            all_gates_pass
        ),
        "regularization_lambda_adoption_authorized": False,
        "constrained_refit_charge_adoption_authorized": False,
        "RESP_stage2_protocol_design_authorized": False,
        "RESP_stage2_execution_authorized": False,
        "charge_adoption_authorized": False,
        "force_field_adoption_authorized": False,
    },
    "outputs": {
        "audit_csv": str(
            output_csv.resolve()
        ),
        "audit_csv_sha256": sha256(
            output_csv
        ),
    },
}

output_json.write_text(
    json.dumps(
        report,
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)

print(f"output_csv={output_csv}")
print(
    f"output_csv_sha256="
    f"{sha256(output_csv)}"
)
print(f"output_json={output_json}")
print(
    f"output_json_sha256="
    f"{sha256(output_json)}"
)


print("\n[14] DECISION")

print(f"decision={decision}")

print(
    "lambda4_final_method_scientific_review_authorized="
    f"{all_gates_pass}"
)

print(
    "regularization_lambda_adoption_authorized=False"
)

print(
    "constrained_refit_charge_adoption_authorized=False"
)

print(
    "RESP_stage2_protocol_design_authorized=False"
)

print(
    "RESP_stage2_execution_authorized=False"
)

print("charge_adoption_authorized=False")
print("force_field_adoption_authorized=False")
print("=" * 100)
