#!/usr/bin/env python3
"""
DAY039 / D039-A8

Scientific comparison of lambda=1 and lambda=10 constrained-refit
candidates for QM_F06_UPPER_V7A_R1.

No candidate is adopted.
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

EXPECTED_A7_DECISION = (
    "D039_A7_LAMBDA1_CHEMICAL_PLAUSIBILITY_AUDIT_PASS_"
    "ATOM_LEVEL_SCIENTIFIC_REVIEW_AUTHORIZED"
)

LAMBDAS = (1.0, 10.0)


def lambda_column(value: float) -> str:
    label = (
        f"{value:.0e}"
        .replace("+", "")
        .replace("-", "m")
    )

    return (
        f"constrained_refit_lambda_{label}_charge_e"
    )


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


print("=" * 100)
print("DAY039 / D039-A8 — LAMBDA 1 VERSUS LAMBDA 10 SCIENTIFIC COMPARISON")
print("=" * 100)


print("\n[1] SOURCE EXECUTION")

require_file(LATEST_POINTER)

execution_dir = (
    ROOT
    / LATEST_POINTER.read_text(
        encoding="utf-8"
    ).strip()
)

a7_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA1_CHEMICAL_PLAUSIBILITY.json"
)

a5_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_CONSTRAINED_REFIT_AUDIT.json"
)

candidate_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_CONSTRAINED_REFIT_REGULARIZATION_PATH.csv"
)

for path in (
    a7_json,
    a5_json,
    candidate_csv,
):
    require_file(path)
    print(f"FOUND  {path}")


print("\n[2] UPSTREAM AUTHORIZATION")

a7_report = load_json(a7_json)

if (
    a7_report.get("decision")
    != EXPECTED_A7_DECISION
):
    raise RuntimeError(
        "Unexpected A7 decision.\n"
        f"Observed: {a7_report.get('decision')}"
    )

authorizations = a7_report.get(
    "authorizations",
    {},
)

if (
    authorizations.get(
        "lambda1_atom_level_scientific_review_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Lambda comparison review is not authorized"
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

print("A7_decision_gate                  = PASS")
print("scientific_comparison_gate        = PASS")
print("lambda_adoption_blocked_gate      = PASS")


print("\n[3] LOAD CANDIDATE CHARGES")

with candidate_csv.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    rows = list(
        csv.DictReader(handle)
    )

if len(rows) != 37:
    raise RuntimeError(
        f"Expected 37 retained atoms, observed {len(rows)}"
    )

for row in rows:
    row["initial_charge_e"] = float(
        row["RESP_stage1_charge_e"]
    )

    for value in LAMBDAS:
        row[
            f"lambda_{value}_charge_e"
        ] = float(
            row[
                lambda_column(value)
            ]
        )

print(f"real_atom_count = {len(rows)}")
print("candidate_contract_gate = PASS")


print("\n[4] LOAD ELECTROSTATIC METRICS")

a5_report = load_json(a5_json)

path_records = a5_report.get(
    "regularization_path",
    [],
)

metrics_by_lambda = {}

for record in path_records:
    value = float(
        record["regularization_lambda"]
    )

    if value in LAMBDAS:
        metrics_by_lambda[value] = record

if set(metrics_by_lambda) != set(LAMBDAS):
    raise RuntimeError(
        "Could not resolve lambda=1 and lambda=10 metrics"
    )

for value in LAMBDAS:
    electrostatic = metrics_by_lambda[
        value
    ]["electrostatic"]

    charges = metrics_by_lambda[
        value
    ]["charges"]

    print(f"\nlambda={value:g}")
    print(
        f"  RMSE_au="
        f"{electrostatic['RMSE_au']}"
    )
    print(
        f"  pearson_r="
        f"{electrostatic['pearson_r']}"
    )
    print(
        f"  same_sign_fraction="
        f"{electrostatic['same_sign_fraction']}"
    )
    print(
        f"  delta_RMS_e="
        f"{charges['delta_RMS_e']}"
    )
    print(
        f"  delta_max_abs_e="
        f"{charges['delta_max_abs_e']}"
    )
    print(
        f"  maximum_absolute_charge_e="
        f"{charges['maximum_absolute_charge_e']}"
    )
    print(
        f"  sign_change_count="
        f"{charges['sign_change_count']}"
    )


print("\n[5] ATOM-LEVEL SIGN CHANGES")

sign_change_records = {
    value: []
    for value in LAMBDAS
}

for value in LAMBDAS:
    print(f"\nlambda={value:g}")

    for row in rows:
        initial = row[
            "initial_charge_e"
        ]

        candidate = row[
            f"lambda_{value}_charge_e"
        ]

        if sign_changed(
            initial,
            candidate,
        ):
            record = {
                "real_atom_sequence_index": int(
                    row[
                        "real_atom_sequence_index"
                    ]
                ),
                "original_atom_index_0based": int(
                    row[
                        "original_atom_index_0based"
                    ]
                ),
                "atom_id": row["atom_id"],
                "element": row["element"],
                "atom_role": row[
                    "atom_role"
                ],
                "RESP_stage1_charge_e": initial,
                "candidate_charge_e": (
                    candidate
                ),
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


print("\n[6] ELEMENT POLARITY AUDIT")

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
        row["initial_charge_e"]
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
        "RESP_stage1_sum_e": (
            initial_sum
        ),
    }

    for value in LAMBDAS:
        candidate_sum = sum(
            row[
                f"lambda_{value}_charge_e"
            ]
            for row in selected
        )

        expected_sign_violations = 0

        if element == "B":
            expected_sign_violations = sum(
                row[
                    f"lambda_{value}_charge_e"
                ] < 0.0
                for row in selected
            )

        elif element == "N":
            expected_sign_violations = sum(
                row[
                    f"lambda_{value}_charge_e"
                ] > 0.0
                for row in selected
            )

        print(
            f"  lambda={value:g} "
            f"sum={candidate_sum:.9f} "
            f"delta_sum="
            f"{candidate_sum-initial_sum:.9f} "
            f"polarity_violations="
            f"{expected_sign_violations}"
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
            f"lambda_{value}_polarity_violations"
        ] = expected_sign_violations

    element_records.append(
        element_record
    )


print("\n[7] ROLE-LEVEL PERTURBATIONS")

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

    record = {
        "atom_role": role,
        "count": len(selected),
    }

    for value in LAMBDAS:
        deltas = [
            row[
                f"lambda_{value}_charge_e"
            ]
            - row["initial_charge_e"]
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

        sign_changes = sum(
            sign_changed(
                row["initial_charge_e"],
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
            f"sign_changes={sign_changes}"
        )

        record[
            f"lambda_{value}_delta_RMS_e"
        ] = delta_rms

        record[
            f"lambda_{value}_delta_max_abs_e"
        ] = delta_max

        record[
            f"lambda_{value}_sign_change_count"
        ] = sign_changes

    role_records.append(record)


print("\n[8] SCIENTIFIC COMPARISON")

lambda1 = metrics_by_lambda[1.0]
lambda10 = metrics_by_lambda[10.0]

lambda1_rmse = float(
    lambda1["electrostatic"]["RMSE_au"]
)

lambda10_rmse = float(
    lambda10["electrostatic"]["RMSE_au"]
)

lambda1_delta_rms = float(
    lambda1["charges"]["delta_RMS_e"]
)

lambda10_delta_rms = float(
    lambda10["charges"]["delta_RMS_e"]
)

rmse_penalty_ratio = (
    lambda10_rmse / lambda1_rmse
)

perturbation_reduction_fraction = (
    1.0
    - lambda10_delta_rms
    / lambda1_delta_rms
)

sign_change_reduction = (
    len(sign_change_records[1.0])
    - len(sign_change_records[10.0])
)

print(
    f"lambda10_to_lambda1_RMSE_ratio = "
    f"{rmse_penalty_ratio:.16g}"
)

print(
    f"lambda10_charge_perturbation_reduction_fraction = "
    f"{perturbation_reduction_fraction:.16g}"
)

print(
    f"lambda10_sign_change_reduction = "
    f"{sign_change_reduction}"
)

lambda10_only_sign_changes = (
    sign_change_records[10.0]
)

lambda10_only_change_is_boundary_B = (
    len(lambda10_only_sign_changes) == 1
    and lambda10_only_sign_changes[0][
        "element"
    ] == "B"
    and "BOUNDARY" in (
        lambda10_only_sign_changes[0][
            "atom_role"
        ]
    )
)

print(
    "lambda10_single_sign_change_is_boundary_B = "
    f"{lambda10_only_change_is_boundary_B}"
)


print("\n[9] WRITE OUTPUTS")

comparison_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA1_VS_LAMBDA10_SCIENTIFIC_COMPARISON.csv"
)

comparison_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA1_VS_LAMBDA10_SCIENTIFIC_COMPARISON.json"
)

fieldnames = [
    "real_atom_sequence_index",
    "original_atom_index_0based",
    "original_atom_index_1based",
    "atom_id",
    "element",
    "atom_role",
    "RESP_stage1_charge_e",
    "lambda_1_charge_e",
    "lambda_1_delta_e",
    "lambda_1_sign_changed",
    "lambda_10_charge_e",
    "lambda_10_delta_e",
    "lambda_10_sign_changed",
]

with comparison_csv.open(
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
            "initial_charge_e"
        ]

        q1 = row[
            "lambda_1.0_charge_e"
        ]

        q10 = row[
            "lambda_10.0_charge_e"
        ]

        writer.writerow(
            {
                "real_atom_sequence_index": int(
                    row[
                        "real_atom_sequence_index"
                    ]
                ),
                "original_atom_index_0based": int(
                    row[
                        "original_atom_index_0based"
                    ]
                ),
                "original_atom_index_1based": int(
                    row[
                        "original_atom_index_1based"
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
                "RESP_stage1_charge_e": initial,
                "lambda_1_charge_e": q1,
                "lambda_1_delta_e": (
                    q1 - initial
                ),
                "lambda_1_sign_changed": (
                    sign_changed(
                        initial,
                        q1,
                    )
                ),
                "lambda_10_charge_e": q10,
                "lambda_10_delta_e": (
                    q10 - initial
                ),
                "lambda_10_sign_changed": (
                    sign_changed(
                        initial,
                        q10,
                    )
                ),
            }
        )


print("\n[10] GATES")

finite_gate = all(
    math.isfinite(
        row[
            f"lambda_{value}_charge_e"
        ]
    )
    for row in rows
    for value in LAMBDAS
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
    for value in LAMBDAS
)

lambda10_fewer_sign_changes_gate = (
    len(sign_change_records[10.0])
    < len(sign_change_records[1.0])
)

lambda10_lower_perturbation_gate = (
    lambda10_delta_rms
    < lambda1_delta_rms
)

lambda1_better_electrostatics_gate = (
    lambda1_rmse
    < lambda10_rmse
)

gates = {
    "upstream_decision_gate": True,
    "candidate_count_gate": (
        len(rows) == 37
    ),
    "finite_candidate_gate": (
        finite_gate
    ),
    "neutrality_gate": (
        neutrality_gate
    ),
    "lambda10_fewer_sign_changes_gate": (
        lambda10_fewer_sign_changes_gate
    ),
    "lambda10_lower_charge_perturbation_gate": (
        lambda10_lower_perturbation_gate
    ),
    "lambda1_better_electrostatics_gate": (
        lambda1_better_electrostatics_gate
    ),
    "comparison_csv_created_gate": (
        comparison_csv.is_file()
        and comparison_csv.stat().st_size > 0
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


print("\n[11] WRITE JSON REPORT")

decision = (
    "D039_A8_LAMBDA1_VS_LAMBDA10_COMPARISON_PASS_"
    "FINAL_LAMBDA_REVIEW_AUTHORIZED"
    if all_gates_pass
    else
    "D039_A8_LAMBDA1_VS_LAMBDA10_COMPARISON_"
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
        "A7_json": str(
            a7_json.resolve()
        ),
        "A7_json_sha256": sha256(
            a7_json
        ),
        "A5_json": str(
            a5_json.resolve()
        ),
        "A5_json_sha256": sha256(
            a5_json
        ),
        "candidate_csv": str(
            candidate_csv.resolve()
        ),
        "candidate_csv_sha256": sha256(
            candidate_csv
        ),
    },
    "electrostatic_metrics": {
        "lambda_1": (
            metrics_by_lambda[1.0][
                "electrostatic"
            ]
        ),
        "lambda_10": (
            metrics_by_lambda[10.0][
                "electrostatic"
            ]
        ),
    },
    "charge_metrics": {
        "lambda_1": (
            metrics_by_lambda[1.0][
                "charges"
            ]
        ),
        "lambda_10": (
            metrics_by_lambda[10.0][
                "charges"
            ]
        ),
    },
    "sign_change_records": {
        "lambda_1": (
            sign_change_records[1.0]
        ),
        "lambda_10": (
            sign_change_records[10.0]
        ),
    },
    "element_records": (
        element_records
    ),
    "role_records": (
        role_records
    ),
    "comparison_summary": {
        "lambda10_to_lambda1_RMSE_ratio": (
            rmse_penalty_ratio
        ),
        "lambda10_charge_perturbation_reduction_fraction": (
            perturbation_reduction_fraction
        ),
        "lambda10_sign_change_reduction": (
            sign_change_reduction
        ),
        "lambda10_single_sign_change_is_boundary_B": (
            lambda10_only_change_is_boundary_B
        ),
    },
    "gates": gates,
    "authorizations": {
        "final_lambda_scientific_review_authorized": (
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
        "comparison_csv": str(
            comparison_csv.resolve()
        ),
        "comparison_csv_sha256": sha256(
            comparison_csv
        ),
    },
}

comparison_json.write_text(
    json.dumps(
        report,
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)

print(f"comparison_csv = {comparison_csv}")
print(
    f"comparison_csv_sha256 = "
    f"{sha256(comparison_csv)}"
)
print(f"comparison_json = {comparison_json}")
print(
    f"comparison_json_sha256 = "
    f"{sha256(comparison_json)}"
)


print("\n[12] DECISION")

print(f"decision = {decision}")
print(
    "final_lambda_scientific_review_authorized = "
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
