#!/usr/bin/env python3
"""
DAY039 / D039-A7

Chemical and topological plausibility audit of the lambda=1 constrained
37-real-atom charge candidate.

The audit also compares lambda = 0.1, 1 and 10.

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

EXPECTED_A6_DECISION = (
    "D039_A6_REGULARIZATION_SELECTION_REVIEW_PASS_"
    "CANDIDATE_LAMBDA_PROPOSED_NOT_ADOPTED"
)

SELECTED_LAMBDA = 1.0
COMPARISON_LAMBDAS = (0.1, 1.0, 10.0)


def lambda_column(value: float) -> str:
    label = (
        f"{value:.0e}"
        .replace("+", "")
        .replace("-", "m")
    )

    return (
        f"constrained_refit_lambda_{label}_charge_e"
    )


def parse_json_list(value: str) -> list:
    parsed = json.loads(value)

    if not isinstance(parsed, list):
        raise RuntimeError(
            f"Expected JSON list, observed: {value!r}"
        )

    return parsed


print("=" * 100)
print("DAY039 / D039-A7 — LAMBDA=1 CHEMICAL PLAUSIBILITY AUDIT")
print("=" * 100)


print("\n[1] SOURCE EXECUTION")

for path in (
    LATEST_POINTER,
    NOMINAL_EDGES,
):
    require_file(path)

execution_dir = (
    ROOT
    / LATEST_POINTER.read_text(
        encoding="utf-8"
    ).strip()
)

selection_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_REGULARIZATION_SELECTION_REVIEW.json"
)

selection_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_REGULARIZATION_SELECTION_REVIEW.csv"
)

candidate_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_CONSTRAINED_REFIT_REGULARIZATION_PATH.csv"
)

transferability_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1_TRANSFERABILITY.csv"
)

for path in (
    selection_json,
    selection_csv,
    candidate_csv,
    transferability_csv,
):
    require_file(path)
    print(f"FOUND  {path}")


print("\n[2] UPSTREAM AUTHORIZATION")

selection_report = load_json(
    selection_json
)

if (
    selection_report.get("decision")
    != EXPECTED_A6_DECISION
):
    raise RuntimeError(
        "Unexpected A6 decision.\n"
        f"Observed: {selection_report.get('decision')}"
    )

recommended = selection_report.get(
    "recommended_candidate_for_review",
    {},
)

observed_lambda = float(
    recommended.get(
        "regularization_lambda"
    )
)

if not math.isclose(
    observed_lambda,
    SELECTED_LAMBDA,
    rel_tol=0.0,
    abs_tol=1.0e-15,
):
    raise RuntimeError(
        "The A6 recommended candidate is not lambda=1"
    )

authorizations = selection_report.get(
    "authorizations",
    {},
)

if (
    authorizations.get(
        "recommended_lambda_scientific_review_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Lambda scientific review is not authorized"
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

print("A6_decision_gate                    = PASS")
print("recommended_lambda_1_gate           = PASS")
print("scientific_review_authorized_gate    = PASS")
print("lambda_adoption_blocked_gate         = PASS")


print("\n[3] LOAD REAL-ATOM CANDIDATES")

with candidate_csv.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    candidate_rows = list(
        csv.DictReader(handle)
    )

if len(candidate_rows) != 37:
    raise RuntimeError(
        f"Expected 37 retained atoms, observed {len(candidate_rows)}"
    )

required_columns = {
    "real_atom_sequence_index",
    "original_atom_index_0based",
    "original_atom_index_1based",
    "atom_id",
    "element",
    "atom_role",
    "RESP_stage1_charge_e",
}

for value in COMPARISON_LAMBDAS:
    required_columns.add(
        lambda_column(value)
    )

missing = required_columns - set(
    candidate_rows[0]
)

if missing:
    raise RuntimeError(
        f"Missing candidate columns: {sorted(missing)}"
    )

candidate_by_id = {}

for row in candidate_rows:
    atom_id = row["atom_id"].strip()

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
                lambda_column(value)
            ]
        )

    candidate_by_id[atom_id] = row

print(f"real_atom_count = {len(candidate_rows)}")
print("candidate_column_contract_gate = PASS")


print("\n[4] LOAD TOPOLOGICAL GRAPH")

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
        first in candidate_by_id
        and second in candidate_by_id
    ):
        distance = float(
            edge["distance_A"]
        )

        neighbors[first].append(
            {
                "atom_id": second,
                "distance_A": distance,
                "edge_type": edge[
                    "edge_type"
                ],
            }
        )

        neighbors[second].append(
            {
                "atom_id": first,
                "distance_A": distance,
                "edge_type": edge[
                    "edge_type"
                ],
            }
        )

real_real_edge_count = sum(
    len(values)
    for values in neighbors.values()
) // 2

print(f"nominal_edge_count   = {len(edge_rows)}")
print(f"real_real_edge_count = {real_real_edge_count}")


print("\n[5] IDENTIFY LAMBDA=1 SIGN CHANGES")

sign_change_rows = []

for row in candidate_rows:
    initial = row[
        "RESP_stage1_charge_e_float"
    ]

    candidate = row[
        "lambda_1.0_charge_e"
    ]

    sign_changed = (
        initial != 0.0
        and candidate != 0.0
        and math.copysign(
            1.0,
            initial,
        )
        != math.copysign(
            1.0,
            candidate,
        )
    )

    if sign_changed:
        sign_change_rows.append(row)

print(
    f"lambda_1_sign_change_count = "
    f"{len(sign_change_rows)}"
)

if len(sign_change_rows) != 5:
    raise RuntimeError(
        "Expected exactly five lambda=1 sign changes"
    )

for row in sign_change_rows:
    atom_id = row["atom_id"]

    print(
        f"\natom_sequence_index="
        f"{row['real_atom_sequence_index_int']:>2} "
        f"original_index="
        f"{row['original_atom_index_0based_int']:>2} "
        f"atom_id={atom_id} "
        f"element={row['element']} "
        f"role={row['atom_role']}"
    )

    print(
        f"  RESP1="
        f"{row['RESP_stage1_charge_e_float']: .9f}"
    )

    for value in COMPARISON_LAMBDAS:
        charge = row[
            f"lambda_{value}_charge_e"
        ]

        print(
            f"  lambda={value:g} "
            f"charge={charge: .9f} "
            f"delta="
            f"{charge - row['RESP_stage1_charge_e_float']: .9f}"
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
        neighbor_row = candidate_by_id[
            neighbor["atom_id"]
        ]

        print(
            f"    neighbor={neighbor['atom_id']} "
            f"element={neighbor_row['element']} "
            f"role={neighbor_row['atom_role']} "
            f"distance_A={neighbor['distance_A']:.6f} "
            f"RESP1="
            f"{neighbor_row['RESP_stage1_charge_e_float']: .6f} "
            f"lambda1="
            f"{neighbor_row['lambda_1.0_charge_e']: .6f}"
        )


print("\n[6] ELEMENT- AND ROLE-LEVEL SUMMARIES")

summary_records = []

for grouping_name, grouping_key in (
    ("element", "element"),
    ("atom_role", "atom_role"),
):
    groups: dict[str, list[dict]] = defaultdict(list)

    for row in candidate_rows:
        groups[
            row[grouping_key]
        ].append(row)

    print(f"\nGrouping = {grouping_name}")

    for group_value in sorted(groups):
        selected = groups[
            group_value
        ]

        initial_sum = sum(
            row[
                "RESP_stage1_charge_e_float"
            ]
            for row in selected
        )

        print(
            f"\n{grouping_name}={group_value!r} "
            f"count={len(selected)} "
            f"RESP1_sum={initial_sum:.9f}"
        )

        record = {
            "grouping": grouping_name,
            "group_value": group_value,
            "count": len(selected),
            "RESP_stage1_sum_e": initial_sum,
        }

        for value in COMPARISON_LAMBDAS:
            charges = [
                row[
                    f"lambda_{value}_charge_e"
                ]
                for row in selected
            ]

            candidate_sum = sum(
                charges
            )

            deltas = [
                charge
                - row[
                    "RESP_stage1_charge_e_float"
                ]
                for charge, row in zip(
                    charges,
                    selected,
                )
            ]

            sign_change_count = sum(
                row[
                    "RESP_stage1_charge_e_float"
                ]
                != 0.0
                and charge != 0.0
                and math.copysign(
                    1.0,
                    row[
                        "RESP_stage1_charge_e_float"
                    ],
                )
                != math.copysign(
                    1.0,
                    charge,
                )
                for charge, row in zip(
                    charges,
                    selected,
                )
            )

            print(
                f"  lambda={value:g} "
                f"sum={candidate_sum:.9f} "
                f"delta_sum="
                f"{candidate_sum - initial_sum:.9f} "
                f"delta_RMS="
                f"{math.sqrt(sum(d*d for d in deltas)/len(deltas)):.9f} "
                f"sign_changes="
                f"{sign_change_count}"
            )

            record[
                f"lambda_{value}_sum_e"
            ] = candidate_sum

            record[
                f"lambda_{value}_delta_sum_e"
            ] = (
                candidate_sum
                - initial_sum
            )

            record[
                f"lambda_{value}_sign_change_count"
            ] = sign_change_count

        summary_records.append(record)


print("\n[7] LARGEST LAMBDA=1 PERTURBATIONS")

ranked = sorted(
    candidate_rows,
    key=lambda row: abs(
        row["lambda_1.0_charge_e"]
        - row[
            "RESP_stage1_charge_e_float"
        ]
    ),
    reverse=True,
)

largest_records = []

for rank, row in enumerate(
    ranked[:15],
    start=1,
):
    delta = (
        row["lambda_1.0_charge_e"]
        - row[
            "RESP_stage1_charge_e_float"
        ]
    )

    sign_changed = row in sign_change_rows

    record = {
        "rank": rank,
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
        "RESP_stage1_charge_e": (
            row[
                "RESP_stage1_charge_e_float"
            ]
        ),
        "lambda_1_charge_e": (
            row[
                "lambda_1.0_charge_e"
            ]
        ),
        "delta_e": delta,
        "absolute_delta_e": abs(
            delta
        ),
        "sign_changed": (
            sign_changed
        ),
    }

    largest_records.append(
        record
    )

    print(
        f"rank={rank:>2} "
        f"atom_sequence_index="
        f"{record['real_atom_sequence_index']:>2} "
        f"original_index="
        f"{record['original_atom_index_0based']:>2} "
        f"atom_id={record['atom_id']:<30} "
        f"element={record['element']} "
        f"RESP1={record['RESP_stage1_charge_e']: .6f} "
        f"lambda1={record['lambda_1_charge_e']: .6f} "
        f"delta={record['delta_e']: .6f} "
        f"sign_changed={record['sign_changed']}"
    )


print("\n[8] WRITE AUDIT TABLE")

audit_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA1_CHEMICAL_PLAUSIBILITY.csv"
)

fieldnames = [
    "real_atom_sequence_index",
    "original_atom_index_0based",
    "original_atom_index_1based",
    "atom_id",
    "element",
    "atom_role",
    "RESP_stage1_charge_e",
    "lambda_0.1_charge_e",
    "lambda_1_charge_e",
    "lambda_10_charge_e",
    "lambda_1_delta_e",
    "lambda_1_sign_changed",
    "real_neighbor_count",
    "real_neighbor_ids",
]

with audit_csv.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
    )

    writer.writeheader()

    for row in candidate_rows:
        initial = row[
            "RESP_stage1_charge_e_float"
        ]

        candidate = row[
            "lambda_1.0_charge_e"
        ]

        sign_changed = (
            initial != 0.0
            and candidate != 0.0
            and math.copysign(
                1.0,
                initial,
            )
            != math.copysign(
                1.0,
                candidate,
            )
        )

        neighbor_ids = sorted(
            item["atom_id"]
            for item in neighbors[
                row["atom_id"]
            ]
        )

        writer.writerow(
            {
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
                "original_atom_index_1based": int(
                    row[
                        "original_atom_index_1based"
                    ]
                ),
                "atom_id": row["atom_id"],
                "element": row["element"],
                "atom_role": row[
                    "atom_role"
                ],
                "RESP_stage1_charge_e": (
                    initial
                ),
                "lambda_0.1_charge_e": (
                    row[
                        "lambda_0.1_charge_e"
                    ]
                ),
                "lambda_1_charge_e": (
                    candidate
                ),
                "lambda_10_charge_e": (
                    row[
                        "lambda_10.0_charge_e"
                    ]
                ),
                "lambda_1_delta_e": (
                    candidate - initial
                ),
                "lambda_1_sign_changed": (
                    sign_changed
                ),
                "real_neighbor_count": len(
                    neighbor_ids
                ),
                "real_neighbor_ids": json.dumps(
                    neighbor_ids
                ),
            }
        )


print("\n[9] SCIENTIFIC GATES")

finite_candidate_gate = all(
    math.isfinite(
        row[
            f"lambda_{value}_charge_e"
        ]
    )
    for row in candidate_rows
    for value in COMPARISON_LAMBDAS
)

lambda1_neutrality_gate = (
    abs(
        sum(
            row[
                "lambda_1.0_charge_e"
            ]
            for row in candidate_rows
        )
    )
    <= 1.0e-10
)

sign_change_reproduced_gate = (
    len(sign_change_rows) == 5
)

topology_available_gate = (
    real_real_edge_count > 0
)

gates = {
    "upstream_decision_gate": True,
    "recommended_lambda_1_gate": True,
    "real_atom_count_gate": (
        len(candidate_rows) == 37
    ),
    "candidate_column_contract_gate": True,
    "finite_candidate_gate": (
        finite_candidate_gate
    ),
    "lambda1_neutrality_gate": (
        lambda1_neutrality_gate
    ),
    "lambda1_sign_change_reproduced_gate": (
        sign_change_reproduced_gate
    ),
    "topology_available_gate": (
        topology_available_gate
    ),
    "audit_csv_created_gate": (
        audit_csv.is_file()
        and audit_csv.stat().st_size > 0
    ),
    "no_lambda_adopted_gate": True,
    "no_charge_adopted_gate": True,
    "RESP_stage2_not_executed_gate": True,
}

for name, value in gates.items():
    print(
        f"{name} = "
        f"{'PASS' if value else 'FAIL'}"
    )

all_gates_pass = all(
    gates.values()
)


print("\n[10] WRITE JSON REPORT")

audit_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA1_CHEMICAL_PLAUSIBILITY.json"
)

decision = (
    "D039_A7_LAMBDA1_CHEMICAL_PLAUSIBILITY_AUDIT_PASS_"
    "ATOM_LEVEL_SCIENTIFIC_REVIEW_AUTHORIZED"
    if all_gates_pass
    else
    "D039_A7_LAMBDA1_CHEMICAL_PLAUSIBILITY_AUDIT_"
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
        "selection_json": str(
            selection_json.resolve()
        ),
        "selection_json_sha256": sha256(
            selection_json
        ),
        "candidate_csv": str(
            candidate_csv.resolve()
        ),
        "candidate_csv_sha256": sha256(
            candidate_csv
        ),
        "nominal_edges": str(
            NOMINAL_EDGES.resolve()
        ),
        "nominal_edges_sha256": sha256(
            NOMINAL_EDGES
        ),
    },
    "selected_lambda_for_review": (
        SELECTED_LAMBDA
    ),
    "comparison_lambdas": list(
        COMPARISON_LAMBDAS
    ),
    "lambda1_sign_change_count": len(
        sign_change_rows
    ),
    "lambda1_sign_change_atoms": [
        {
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
            "atom_id": row[
                "atom_id"
            ],
            "element": row[
                "element"
            ],
            "atom_role": row[
                "atom_role"
            ],
            "RESP_stage1_charge_e": (
                row[
                    "RESP_stage1_charge_e_float"
                ]
            ),
            "lambda_1_charge_e": (
                row[
                    "lambda_1.0_charge_e"
                ]
            ),
            "neighbor_ids": sorted(
                item["atom_id"]
                for item in neighbors[
                    row["atom_id"]
                ]
            ),
        }
        for row in sign_change_rows
    ],
    "group_summaries": (
        summary_records
    ),
    "largest_lambda1_perturbations": (
        largest_records
    ),
    "gates": gates,
    "authorizations": {
        "lambda1_atom_level_scientific_review_authorized": (
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
            audit_csv.resolve()
        ),
        "audit_csv_sha256": sha256(
            audit_csv
        ),
    },
}

audit_json.write_text(
    json.dumps(
        report,
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)

print(f"audit_csv = {audit_csv}")
print(
    f"audit_csv_sha256 = "
    f"{sha256(audit_csv)}"
)
print(f"audit_json = {audit_json}")
print(
    f"audit_json_sha256 = "
    f"{sha256(audit_json)}"
)


print("\n[11] DECISION")

print(f"decision = {decision}")
print(
    "lambda1_atom_level_scientific_review_authorized = "
    f"{all_gates_pass}"
)
print(
    "regularization_lambda_adoption_authorized = False"
)
print(
    "constrained_refit_charge_adoption_authorized = False"
)
print(
    "RESP_stage2_protocol_design_authorized = False"
)
print(
    "RESP_stage2_execution_authorized = False"
)
print("charge_adoption_authorized = False")
print("force_field_adoption_authorized = False")
print("=" * 100)
