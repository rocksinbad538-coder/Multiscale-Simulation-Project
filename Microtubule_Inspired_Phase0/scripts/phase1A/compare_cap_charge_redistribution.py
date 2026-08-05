#!/usr/bin/env python3
"""
DAY039 / D039-A2

Exploratory comparison of artificial-cap charge redistribution
strategies for QM_F06_UPPER_V7A_R1.

Strategies
----------
1. PARENT_ONLY
2. FIRST_SHELL_EQUAL
3. FIRST_SHELL_ABS_Q_WEIGHTED

No strategy is adopted or authorized for RESP Stage 2.
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

EXPECTED_MAPPING_DECISION = (
    "D038_F3_CAP_PARENT_MAPPING_VALIDATION_PASS_"
    "REDISTRIBUTION_DESIGN_REVIEW_AUTHORIZED"
)

EXPECTED_REAL_ATOM_COUNT = 37
EXPECTED_CAP_COUNT = 15
TARGET_TOTAL_CHARGE_E = 0.0
CHARGE_TOLERANCE_E = 5.0e-6
WEIGHT_EPSILON = 1.0e-12


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()

    if normalized in {"true", "1", "yes", "y"}:
        return True

    if normalized in {"false", "0", "no", "n", ""}:
        return False

    raise RuntimeError(
        f"Unrecognized Boolean value: {value!r}"
    )


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def population_std(values: list[float]) -> float:
    center = mean(values)

    return math.sqrt(
        sum((value - center) ** 2 for value in values)
        / len(values)
    )


def summarize_charge_vector(
    original: dict[str, float],
    redistributed: dict[str, float],
) -> dict:
    atom_ids = sorted(original)

    original_values = [
        original[atom_id]
        for atom_id in atom_ids
    ]

    final_values = [
        redistributed[atom_id]
        for atom_id in atom_ids
    ]

    differences = [
        redistributed[atom_id] - original[atom_id]
        for atom_id in atom_ids
    ]

    absolute_differences = [
        abs(value)
        for value in differences
    ]

    sign_changes = [
        atom_id
        for atom_id in atom_ids
        if (
            original[atom_id] != 0.0
            and redistributed[atom_id] != 0.0
            and math.copysign(
                1.0,
                original[atom_id],
            )
            != math.copysign(
                1.0,
                redistributed[atom_id],
            )
        )
    ]

    return {
        "atom_count": len(atom_ids),
        "original_charge_sum_e": sum(
            original_values
        ),
        "final_charge_sum_e": sum(
            final_values
        ),
        "redistributed_charge_sum_e": sum(
            differences
        ),
        "minimum_final_charge_e": min(
            final_values
        ),
        "maximum_final_charge_e": max(
            final_values
        ),
        "maximum_absolute_final_charge_e": max(
            abs(value)
            for value in final_values
        ),
        "difference_mean_e": mean(
            differences
        ),
        "difference_std_e": population_std(
            differences
        ),
        "difference_MAE_e": mean(
            absolute_differences
        ),
        "difference_RMS_e": math.sqrt(
            mean(
                [
                    value * value
                    for value in differences
                ]
            )
        ),
        "difference_max_abs_e": max(
            absolute_differences
        ),
        "changed_atom_count": sum(
            abs(value) > 1.0e-12
            for value in differences
        ),
        "sign_change_count": len(
            sign_changes
        ),
        "sign_change_atom_ids": (
            sign_changes
        ),
    }


print("=" * 100)
print("DAY039 / D039-A2 — CAP-CHARGE REDISTRIBUTION STRATEGY COMPARISON")
print("=" * 100)


print("\n[1] SOURCE FILES")

require_file(LATEST_POINTER)
require_file(NOMINAL_EDGES)

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

mapping_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_CAP_PARENT_MAPPING_VALIDATED.json"
)

require_file(transferability_csv)
require_file(mapping_json)

for path in (
    transferability_csv,
    mapping_json,
    NOMINAL_EDGES,
):
    print(f"FOUND  {path}")


print("\n[2] UPSTREAM AUTHORIZATION")

mapping_report = load_json(
    mapping_json
)

if (
    mapping_report.get("decision")
    != EXPECTED_MAPPING_DECISION
):
    raise RuntimeError(
        "Unexpected cap-parent mapping decision.\n"
        f"Observed: {mapping_report.get('decision')}"
    )

authorizations = mapping_report.get(
    "authorizations",
    {},
)

if (
    authorizations.get(
        "cap_charge_redistribution_design_review_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Redistribution design review is not authorized"
    )

if (
    authorizations.get(
        "cap_charge_redistribution_execution_authorized"
    )
    is not False
):
    raise RuntimeError(
        "Unexpected redistribution-execution authorization"
    )

print("mapping_decision_gate               = PASS")
print("redistribution_design_review_gate   = PASS")
print("redistribution_execution_block_gate = PASS")


print("\n[3] LOAD RESP ATOMS")

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
        f"Expected 52 atoms, observed {len(atom_rows)}"
    )

atom_by_id: dict[str, dict] = {}
real_atom_ids: list[str] = []
cap_atom_ids: list[str] = []

for row in atom_rows:
    atom_id = row["atom_id"].strip()

    row["atom_index_0based_int"] = int(
        row["atom_index_0based"]
    )
    row["charge_e_float"] = float(
        row["RESP_stage1_charge_e"]
    )
    row["artificial_cap_bool"] = parse_bool(
        row["artificial_cap"]
    )

    atom_by_id[atom_id] = row

    if row["artificial_cap_bool"]:
        cap_atom_ids.append(atom_id)
    else:
        real_atom_ids.append(atom_id)

if len(real_atom_ids) != EXPECTED_REAL_ATOM_COUNT:
    raise RuntimeError(
        "Unexpected number of real atoms"
    )

if len(cap_atom_ids) != EXPECTED_CAP_COUNT:
    raise RuntimeError(
        "Unexpected number of caps"
    )

original_real_charges = {
    atom_id: atom_by_id[atom_id][
        "charge_e_float"
    ]
    for atom_id in real_atom_ids
}

cap_charge_sum = sum(
    atom_by_id[atom_id]["charge_e_float"]
    for atom_id in cap_atom_ids
)

print(f"real_atom_count       = {len(real_atom_ids)}")
print(f"artificial_cap_count = {len(cap_atom_ids)}")
print(
    f"original_real_charge_sum_e = "
    f"{sum(original_real_charges.values()):.16g}"
)
print(f"cap_charge_sum_e           = {cap_charge_sum:.16g}")


print("\n[4] LOAD REAL-ATOM GRAPH")

with NOMINAL_EDGES.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    edge_rows = list(
        csv.DictReader(handle)
    )

real_neighbors: dict[str, set[str]] = defaultdict(set)

for edge in edge_rows:
    first = edge["first_atom"].strip()
    second = edge["second_atom"].strip()

    if (
        first in original_real_charges
        and second in original_real_charges
    ):
        real_neighbors[first].add(second)
        real_neighbors[second].add(first)

print(f"nominal_edge_count = {len(edge_rows)}")

real_real_edge_count = sum(
    len(neighbors)
    for neighbors in real_neighbors.values()
) // 2

print(f"real_real_edge_count = {real_real_edge_count}")


print("\n[5] CAP-PARENT MAPPINGS")

mappings = mapping_report.get(
    "cap_parent_mappings",
    [],
)

if len(mappings) != EXPECTED_CAP_COUNT:
    raise RuntimeError(
        "Cap-parent mapping count mismatch"
    )

for record in mappings:
    cap_id = record["cap_atom_id"]
    parent_id = record["parent_atom_id"]

    if cap_id not in atom_by_id:
        raise RuntimeError(
            f"Unknown cap atom: {cap_id}"
        )

    if parent_id not in original_real_charges:
        raise RuntimeError(
            f"Parent is not a real RESP atom: {parent_id}"
        )

    print(
        f"cap={cap_id:<43} "
        f"parent={parent_id:<30} "
        f"q_cap={atom_by_id[cap_id]['charge_e_float']: .6f} "
        f"parent_real_neighbors="
        f"{len(real_neighbors[parent_id])}"
    )


def new_charge_vector() -> dict[str, float]:
    return dict(original_real_charges)


def distribute_charge(
    target: dict[str, float],
    recipients: list[str],
    charge_e: float,
    weights: list[float],
) -> None:
    if not recipients:
        raise RuntimeError(
            "Redistribution recipient list is empty"
        )

    if len(recipients) != len(weights):
        raise RuntimeError(
            "Recipient/weight length mismatch"
        )

    weight_sum = sum(weights)

    if weight_sum <= 0.0:
        raise RuntimeError(
            "Redistribution weights do not have "
            "a positive sum"
        )

    normalized_weights = [
        weight / weight_sum
        for weight in weights
    ]

    distributed = 0.0

    for recipient, weight in zip(
        recipients[:-1],
        normalized_weights[:-1],
    ):
        increment = charge_e * weight
        target[recipient] += increment
        distributed += increment

    # Assign remainder to final recipient so that each
    # cap redistribution is numerically conservative.
    target[recipients[-1]] += (
        charge_e - distributed
    )


print("\n[6] APPLY EXPLORATORY STRATEGIES")

strategy_vectors: dict[str, dict[str, float]] = {
    "PARENT_ONLY": new_charge_vector(),
    "FIRST_SHELL_EQUAL": new_charge_vector(),
    "FIRST_SHELL_ABS_Q_WEIGHTED": (
        new_charge_vector()
    ),
}

strategy_cap_details: dict[
    str,
    list[dict],
] = {
    strategy: []
    for strategy in strategy_vectors
}

for mapping in mappings:
    cap_id = mapping["cap_atom_id"]
    parent_id = mapping["parent_atom_id"]
    cap_charge = atom_by_id[
        cap_id
    ]["charge_e_float"]

    # Strategy 1: direct parent absorption.
    recipients_parent = [parent_id]

    distribute_charge(
        strategy_vectors["PARENT_ONLY"],
        recipients_parent,
        cap_charge,
        [1.0],
    )

    strategy_cap_details[
        "PARENT_ONLY"
    ].append(
        {
            "cap_atom_id": cap_id,
            "cap_charge_e": cap_charge,
            "recipients": recipients_parent,
            "weights": [1.0],
        }
    )

    # Shared recipient shell for strategies 2 and 3.
    recipients_shell = [
        parent_id,
        *sorted(
            real_neighbors[parent_id]
        ),
    ]

    # Strategy 2: equal first-shell distribution.
    equal_weights = [
        1.0
        for _ in recipients_shell
    ]

    distribute_charge(
        strategy_vectors[
            "FIRST_SHELL_EQUAL"
        ],
        recipients_shell,
        cap_charge,
        equal_weights,
    )

    strategy_cap_details[
        "FIRST_SHELL_EQUAL"
    ].append(
        {
            "cap_atom_id": cap_id,
            "cap_charge_e": cap_charge,
            "recipients": recipients_shell,
            "weights": [
                weight / sum(equal_weights)
                for weight in equal_weights
            ],
        }
    )

    # Strategy 3: weights proportional to the
    # magnitude of the Stage 1 charges.
    abs_q_weights = [
        abs(original_real_charges[atom_id])
        + WEIGHT_EPSILON
        for atom_id in recipients_shell
    ]

    distribute_charge(
        strategy_vectors[
            "FIRST_SHELL_ABS_Q_WEIGHTED"
        ],
        recipients_shell,
        cap_charge,
        abs_q_weights,
    )

    total_abs_q_weight = sum(
        abs_q_weights
    )

    strategy_cap_details[
        "FIRST_SHELL_ABS_Q_WEIGHTED"
    ].append(
        {
            "cap_atom_id": cap_id,
            "cap_charge_e": cap_charge,
            "recipients": recipients_shell,
            "weights": [
                weight / total_abs_q_weight
                for weight in abs_q_weights
            ],
        }
    )


print("\n[7] STRATEGY METRICS")

strategy_summaries = {}

for strategy_name, charge_vector in (
    strategy_vectors.items()
):
    summary = summarize_charge_vector(
        original_real_charges,
        charge_vector,
    )

    neutrality_error = (
        summary["final_charge_sum_e"]
        - TARGET_TOTAL_CHARGE_E
    )

    summary["neutrality_error_e"] = (
        neutrality_error
    )
    summary["neutrality_gate"] = (
        abs(neutrality_error)
        <= CHARGE_TOLERANCE_E
    )

    strategy_summaries[
        strategy_name
    ] = summary

    print(f"\nStrategy = {strategy_name}")

    for name, value in summary.items():
        print(f"  {name} = {value}")


print("\n[8] ATOMS WITH LARGEST PERTURBATIONS")

for strategy_name, charge_vector in (
    strategy_vectors.items()
):
    ranked = sorted(
        real_atom_ids,
        key=lambda atom_id: abs(
            charge_vector[atom_id]
            - original_real_charges[atom_id]
        ),
        reverse=True,
    )

    print(f"\nStrategy = {strategy_name}")

    for rank, atom_id in enumerate(
        ranked[:12],
        start=1,
    ):
        original_charge = (
            original_real_charges[atom_id]
        )

        final_charge = charge_vector[atom_id]

        print(
            f"rank={rank:>2} "
            f"atom_index="
            f"{atom_by_id[atom_id]['atom_index_0based_int']:>2} "
            f"atom_id={atom_id:<30} "
            f"element={atom_by_id[atom_id]['element']} "
            f"original={original_charge: .6f} "
            f"final={final_charge: .6f} "
            f"delta={final_charge - original_charge: .6f}"
        )


print("\n[9] WRITE COMPARISON TABLE")

comparison_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_CAP_REDISTRIBUTION_STRATEGY_COMPARISON.csv"
)

fieldnames = [
    "atom_index_0based",
    "atom_index_1based",
    "atom_id",
    "element",
    "atom_role",
    "original_RESP_stage1_charge_e",
]

for strategy_name in strategy_vectors:
    fieldnames.extend(
        [
            f"{strategy_name}_charge_e",
            f"{strategy_name}_delta_e",
        ]
    )

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

    for atom_id in sorted(
        real_atom_ids,
        key=lambda item: atom_by_id[
            item
        ]["atom_index_0based_int"],
    ):
        row = atom_by_id[atom_id]

        output_row = {
            "atom_index_0based": (
                row["atom_index_0based_int"]
            ),
            "atom_index_1based": int(
                row["atom_index_1based"]
            ),
            "atom_id": atom_id,
            "element": row["element"],
            "atom_role": row["atom_role"],
            "original_RESP_stage1_charge_e": (
                original_real_charges[atom_id]
            ),
        }

        for strategy_name, charge_vector in (
            strategy_vectors.items()
        ):
            final_charge = charge_vector[
                atom_id
            ]

            output_row[
                f"{strategy_name}_charge_e"
            ] = final_charge

            output_row[
                f"{strategy_name}_delta_e"
            ] = (
                final_charge
                - original_real_charges[atom_id]
            )

        writer.writerow(output_row)


print("\n[10] SCIENTIFIC GATES")

all_neutral_gate = all(
    summary["neutrality_gate"]
    for summary in strategy_summaries.values()
)

finite_values_gate = all(
    math.isfinite(value)
    for vector in strategy_vectors.values()
    for value in vector.values()
)

real_atom_set_preserved_gate = all(
    set(vector)
    == set(original_real_charges)
    for vector in strategy_vectors.values()
)

gates = {
    "upstream_mapping_gate": True,
    "real_atom_count_gate": (
        len(real_atom_ids)
        == EXPECTED_REAL_ATOM_COUNT
    ),
    "cap_count_gate": (
        len(cap_atom_ids)
        == EXPECTED_CAP_COUNT
    ),
    "mapping_count_gate": (
        len(mappings)
        == EXPECTED_CAP_COUNT
    ),
    "real_atom_set_preserved_gate": (
        real_atom_set_preserved_gate
    ),
    "finite_values_gate": (
        finite_values_gate
    ),
    "all_strategies_neutral_gate": (
        all_neutral_gate
    ),
    "comparison_csv_created_gate": (
        comparison_csv.is_file()
        and comparison_csv.stat().st_size > 0
    ),
    "no_strategy_adopted_gate": True,
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

report_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_CAP_REDISTRIBUTION_STRATEGY_COMPARISON.json"
)

decision = (
    "D039_A2_CAP_REDISTRIBUTION_STRATEGY_COMPARISON_PASS_"
    "SCIENTIFIC_SELECTION_REVIEW_AUTHORIZED"
    if all_gates_pass
    else
    "D039_A2_CAP_REDISTRIBUTION_STRATEGY_COMPARISON_"
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
        "transferability_csv": str(
            transferability_csv.resolve()
        ),
        "transferability_csv_sha256": sha256(
            transferability_csv
        ),
        "cap_parent_mapping_json": str(
            mapping_json.resolve()
        ),
        "cap_parent_mapping_json_sha256": sha256(
            mapping_json
        ),
        "nominal_edges": str(
            NOMINAL_EDGES.resolve()
        ),
        "nominal_edges_sha256": sha256(
            NOMINAL_EDGES
        ),
    },
    "charge_summary": {
        "original_real_charge_sum_e": sum(
            original_real_charges.values()
        ),
        "cap_charge_sum_e": cap_charge_sum,
        "target_total_charge_e": (
            TARGET_TOTAL_CHARGE_E
        ),
    },
    "strategy_summaries": (
        strategy_summaries
    ),
    "strategy_cap_details": (
        strategy_cap_details
    ),
    "gates": gates,
    "authorizations": {
        "redistribution_strategy_scientific_review_authorized": (
            all_gates_pass
        ),
        "redistribution_strategy_adoption_authorized": False,
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

report_json.write_text(
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
print(f"report_json = {report_json}")
print(
    f"report_json_sha256 = "
    f"{sha256(report_json)}"
)


print("\n[12] DECISION")

print(f"decision = {decision}")
print(
    "redistribution_strategy_scientific_review_authorized = "
    f"{all_gates_pass}"
)
print("redistribution_strategy_adoption_authorized = False")
print("RESP_stage2_execution_authorized = False")
print("charge_adoption_authorized = False")
print("force_field_adoption_authorized = False")
print("=" * 100)
