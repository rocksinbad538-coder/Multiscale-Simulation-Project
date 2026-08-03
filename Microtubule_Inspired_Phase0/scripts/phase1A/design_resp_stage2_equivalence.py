#!/usr/bin/env python3
"""
DAY038 / D038-F1B

Evidence-based RESP Stage 2 equivalence-design review for
QM_F06_UPPER_V7A_R1.

This script:
- reads the latest authorized Stage 1 transferability dataset;
- evaluates every candidate equivalence group;
- classifies groups conservatively;
- creates reproducible CSV and JSON design artifacts;
- does not execute RESP Stage 2;
- does not adopt any charges.
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

EXPECTED_TRANSFERABILITY_DECISION = (
    "D038_E2_RESP_STAGE1_TRANSFERABILITY_AND_CAP_"
    "AUDIT_PASS_STAGE2_DESIGN_REVIEW_AUTHORIZED"
)

CHARGE_COLUMN = "RESP_stage1_charge_e"

TIGHT_RANGE_THRESHOLD_E = 0.02
HIGH_RANGE_THRESHOLD_E = 0.10


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


print("=" * 100)
print("DAY038 / D038-F1B — RESP STAGE 2 EQUIVALENCE DESIGN")
print("=" * 100)


print("\n[1] SOURCE EXECUTION")

require_file(LATEST_POINTER)

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

transferability_audit = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1_TRANSFERABILITY_AUDIT.json"
)

require_file(transferability_csv)
require_file(transferability_audit)

print(f"execution_dir          = {execution_dir}")
print(f"transferability_csv    = {transferability_csv}")
print(f"transferability_audit  = {transferability_audit}")


print("\n[2] UPSTREAM AUTHORIZATION")

upstream = load_json(
    transferability_audit
)

if (
    upstream.get("decision")
    != EXPECTED_TRANSFERABILITY_DECISION
):
    raise RuntimeError(
        "Unexpected transferability-audit decision.\n"
        f"Expected: {EXPECTED_TRANSFERABILITY_DECISION}\n"
        f"Observed: {upstream.get('decision')}"
    )

authorizations = upstream.get(
    "authorizations",
    {},
)

if (
    authorizations.get(
        "RESP_stage2_design_review_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Stage 2 design review is not authorized"
    )

if (
    authorizations.get(
        "RESP_stage2_execution_authorized"
    )
    is not False
):
    raise RuntimeError(
        "Unexpected Stage 2 execution authorization"
    )

print("upstream_decision_gate              = PASS")
print("stage2_design_review_authorized_gate = PASS")
print("stage2_execution_blocked_gate        = PASS")


print("\n[3] LOAD ATOM-LEVEL DATA")

with transferability_csv.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    rows = list(
        csv.DictReader(handle)
    )

if len(rows) != 52:
    raise RuntimeError(
        f"Expected 52 atom rows, observed {len(rows)}"
    )

required_columns = {
    "atom_index_0based",
    "atom_index_1based",
    "atom_id",
    "element",
    "atom_role",
    "artificial_cap",
    "transfer_status",
    "RESP_atom_class",
    "candidate_equivalence_key",
    "equivalence_enforced",
    CHARGE_COLUMN,
}

missing_columns = (
    required_columns
    - set(rows[0].keys())
)

if missing_columns:
    raise RuntimeError(
        "Required CSV columns are missing:\n"
        + "\n".join(
            sorted(missing_columns)
        )
    )

print(
    f"available_columns = "
    f"{list(rows[0].keys())}"
)

for row in rows:
    row["atom_index_0based_int"] = int(
        row["atom_index_0based"]
    )
    row["atom_index_1based_int"] = int(
        row["atom_index_1based"]
    )
    row["charge_e_float"] = float(
        row[CHARGE_COLUMN]
    )
    row["artificial_cap_bool"] = parse_bool(
        row["artificial_cap"]
    )
    row["equivalence_enforced_bool"] = parse_bool(
        row["equivalence_enforced"]
    )

    if not math.isfinite(
        row["charge_e_float"]
    ):
        raise RuntimeError(
            "Non-finite Stage 1 charge detected"
        )

indices = [
    row["atom_index_0based_int"]
    for row in rows
]

if indices != list(range(52)):
    raise RuntimeError(
        "Atom order is not exactly 0..51"
    )

print(f"atom_count = {len(rows)}")
print("column_contract_gate = PASS")
print("atom_order_gate      = PASS")


print("\n[4] GROUP CANDIDATE EQUIVALENCES")

groups: dict[str, list[dict]] = defaultdict(list)

for row in rows:
    key = row[
        "candidate_equivalence_key"
    ].strip()

    if not key:
        raise RuntimeError(
            "Blank candidate equivalence key detected"
        )

    groups[key].append(row)

design_records: list[dict] = []

for key in sorted(groups):
    selected = groups[key]

    charges = [
        row["charge_e_float"]
        for row in selected
    ]

    count = len(selected)
    charge_min = min(charges)
    charge_max = max(charges)
    charge_range = charge_max - charge_min
    charge_mean = mean(charges)
    charge_std = population_std(charges)

    all_caps = all(
        row["artificial_cap_bool"]
        for row in selected
    )

    all_transferable = all(
        row["transfer_status"].strip()
        == "TRANSFERABLE_REAL_ATOM"
        for row in selected
    )

    element_set = sorted(
        set(
            row["element"].strip()
            for row in selected
        )
    )

    role_set = sorted(
        set(
            row["atom_role"].strip()
            for row in selected
        )
    )

    enforced_count = sum(
        row["equivalence_enforced_bool"]
        for row in selected
    )

    if count == 1:
        classification = "KEEP_UNIQUE"
        rationale = (
            "Singleton group; no equivalence decision is applicable."
        )

    elif all_caps:
        classification = "FORBIDDEN"
        rationale = (
            "Artificial QM boundary caps are excluded from direct "
            "transfer and must not define transferable Stage 2 "
            "charge equivalences."
        )

    elif not all_transferable:
        classification = "FORBIDDEN"
        rationale = (
            "Group contains one or more non-transferable atoms."
        )

    elif charge_range <= TIGHT_RANGE_THRESHOLD_E:
        classification = "RECOMMENDED_FOR_REVIEW"
        rationale = (
            "Transferable nonsingleton group with Stage 1 charge "
            f"range <= {TIGHT_RANGE_THRESHOLD_E:.3f} e; equality "
            "may be defensible but requires explicit structural "
            "symmetry confirmation."
        )

    elif charge_range >= HIGH_RANGE_THRESHOLD_E:
        classification = "FORBIDDEN"
        rationale = (
            "Large Stage 1 within-group charge dispersion "
            f"(range >= {HIGH_RANGE_THRESHOLD_E:.3f} e) contradicts "
            "electrostatic equivalence."
        )

    else:
        classification = "REQUIRES_SCIENTIFIC_REVIEW"
        rationale = (
            "Intermediate Stage 1 dispersion; equivalence cannot "
            "be justified from charge similarity alone."
        )

    record = {
        "candidate_equivalence_key": key,
        "classification": classification,
        "rationale": rationale,
        "count": count,
        "atom_indices_0based": [
            row["atom_index_0based_int"]
            for row in selected
        ],
        "atom_indices_1based": [
            row["atom_index_1based_int"]
            for row in selected
        ],
        "atom_ids": [
            row["atom_id"]
            for row in selected
        ],
        "elements": element_set,
        "atom_roles": role_set,
        "all_artificial_caps": all_caps,
        "all_transferable_real_atoms": (
            all_transferable
        ),
        "previously_enforced_count": (
            enforced_count
        ),
        "RESP_stage1_charge_mean_e": (
            charge_mean
        ),
        "RESP_stage1_charge_std_e": (
            charge_std
        ),
        "RESP_stage1_charge_min_e": (
            charge_min
        ),
        "RESP_stage1_charge_max_e": (
            charge_max
        ),
        "RESP_stage1_charge_range_e": (
            charge_range
        ),
        "proposed_stage2_equivalence_enforced": False,
    }

    design_records.append(record)


print("\n[5] DESIGN MATRIX")

classification_counts: dict[str, int] = defaultdict(int)

for record in design_records:
    classification_counts[
        record["classification"]
    ] += 1

    print(
        f"{record['classification']:<28} "
        f"count={record['count']:>2} "
        f"range_e={record['RESP_stage1_charge_range_e']:.6f} "
        f"indices={record['atom_indices_0based']} "
        f"key={record['candidate_equivalence_key']}"
    )


print("\n[6] CLASSIFICATION SUMMARY")

for classification in sorted(
    classification_counts
):
    print(
        f"{classification} = "
        f"{classification_counts[classification]}"
    )

nonsingleton_records = [
    record
    for record in design_records
    if record["count"] > 1
]

recommended_records = [
    record
    for record in design_records
    if record["classification"]
    == "RECOMMENDED_FOR_REVIEW"
]

forbidden_records = [
    record
    for record in design_records
    if record["classification"]
    == "FORBIDDEN"
]

review_records = [
    record
    for record in design_records
    if record["classification"]
    == "REQUIRES_SCIENTIFIC_REVIEW"
]

print(
    f"total_group_count = "
    f"{len(design_records)}"
)
print(
    f"nonsingleton_group_count = "
    f"{len(nonsingleton_records)}"
)
print(
    f"recommended_for_review_count = "
    f"{len(recommended_records)}"
)
print(
    f"forbidden_count = "
    f"{len(forbidden_records)}"
)
print(
    f"requires_scientific_review_count = "
    f"{len(review_records)}"
)


print("\n[7] WRITE OUTPUTS")

output_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE2_EQUIVALENCE_DESIGN.json"
)

output_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE2_EQUIVALENCE_DESIGN.csv"
)

csv_fieldnames = [
    "candidate_equivalence_key",
    "classification",
    "rationale",
    "count",
    "atom_indices_0based",
    "atom_indices_1based",
    "atom_ids",
    "elements",
    "atom_roles",
    "all_artificial_caps",
    "all_transferable_real_atoms",
    "previously_enforced_count",
    "RESP_stage1_charge_mean_e",
    "RESP_stage1_charge_std_e",
    "RESP_stage1_charge_min_e",
    "RESP_stage1_charge_max_e",
    "RESP_stage1_charge_range_e",
    "proposed_stage2_equivalence_enforced",
]

with output_csv.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=csv_fieldnames,
    )

    writer.writeheader()

    for record in design_records:
        csv_record = dict(record)

        for key in (
            "atom_indices_0based",
            "atom_indices_1based",
            "atom_ids",
            "elements",
            "atom_roles",
        ):
            csv_record[key] = json.dumps(
                csv_record[key]
            )

        writer.writerow(csv_record)


gates = {
    "upstream_decision_gate": True,
    "stage2_design_review_authorized_gate": True,
    "atom_count_gate": len(rows) == 52,
    "atom_order_gate": indices == list(range(52)),
    "column_contract_gate": True,
    "finite_charge_gate": all(
        math.isfinite(
            row["charge_e_float"]
        )
        for row in rows
    ),
    "candidate_key_complete_gate": all(
        bool(
            row[
                "candidate_equivalence_key"
            ].strip()
        )
        for row in rows
    ),
    "no_automatic_equivalence_enforced_gate": all(
        record[
            "proposed_stage2_equivalence_enforced"
        ]
        is False
        for record in design_records
    ),
}

all_gates_pass = all(
    gates.values()
)

decision = (
    "D038_F1_RESP_STAGE2_EQUIVALENCE_DESIGN_PASS_"
    "NO_STAGE2_EQUIVALENCE_AUTOMATICALLY_AUTHORIZED"
    if all_gates_pass
    else
    "D038_F1_RESP_STAGE2_EQUIVALENCE_DESIGN_"
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
        "transferability_audit": str(
            transferability_audit.resolve()
        ),
        "transferability_audit_sha256": sha256(
            transferability_audit
        ),
    },
    "thresholds": {
        "tight_range_threshold_e": (
            TIGHT_RANGE_THRESHOLD_E
        ),
        "high_range_threshold_e": (
            HIGH_RANGE_THRESHOLD_E
        ),
    },
    "classification_counts": dict(
        classification_counts
    ),
    "summary": {
        "total_group_count": len(
            design_records
        ),
        "nonsingleton_group_count": len(
            nonsingleton_records
        ),
        "recommended_for_review_count": len(
            recommended_records
        ),
        "forbidden_count": len(
            forbidden_records
        ),
        "requires_scientific_review_count": len(
            review_records
        ),
        "automatically_enforced_group_count": 0,
    },
    "groups": design_records,
    "gates": gates,
    "authorizations": {
        "RESP_stage2_protocol_design_authorized": (
            all_gates_pass
        ),
        "RESP_stage2_equivalence_enforcement_authorized": False,
        "RESP_stage2_execution_authorized": False,
        "charge_adoption_authorized": False,
        "force_field_adoption_authorized": False,
    },
    "outputs": {
        "design_csv": str(
            output_csv.resolve()
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

print(f"output_csv = {output_csv}")
print(
    f"output_csv_sha256 = "
    f"{sha256(output_csv)}"
)
print(f"output_json = {output_json}")
print(
    f"output_json_sha256 = "
    f"{sha256(output_json)}"
)


print("\n[8] GATES")

for name, value in gates.items():
    print(
        f"{name} = "
        f"{'PASS' if value else 'FAIL'}"
    )


print("\n[9] DECISION")

print(f"decision = {decision}")
print(
    "RESP_stage2_protocol_design_authorized = "
    f"{all_gates_pass}"
)
print(
    "RESP_stage2_equivalence_enforcement_authorized = False"
)
print("RESP_stage2_execution_authorized = False")
print("charge_adoption_authorized = False")
print("force_field_adoption_authorized = False")
print("=" * 100)
