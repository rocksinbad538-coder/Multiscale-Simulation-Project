#!/usr/bin/env python3
"""
DAY038 / D038-E2

Transferability, artificial-cap, and candidate-equivalence audit for
QM_F06_UPPER_V7A_R1 RESP Stage 1 charges.

This audit does not execute RESP Stage 2 and does not adopt charges.
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

EXPECTED_AUDIT_DECISION = (
    "D038_E1_RESP_STAGE1_CANDIDATE_CHARGES_"
    "SCIENTIFIC_AUDIT_PASS_STAGE2_REMAINS_BLOCKED"
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


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def population_std(values: list[float]) -> float:
    center = mean(values)

    return math.sqrt(
        sum((value - center) ** 2 for value in values)
        / len(values)
    )


def numeric_summary(values: list[float]) -> dict:
    if not values:
        return {
            "count": 0,
            "sum_e": 0.0,
            "mean_e": None,
            "std_e": None,
            "minimum_e": None,
            "maximum_e": None,
            "range_e": None,
        }

    return {
        "count": len(values),
        "sum_e": sum(values),
        "mean_e": mean(values),
        "std_e": population_std(values),
        "minimum_e": min(values),
        "maximum_e": max(values),
        "range_e": max(values) - min(values),
    }


print("=" * 100)
print("DAY038 / D038-E2 — RESP STAGE 1 TRANSFERABILITY AUDIT")
print("=" * 100)


print("\n[1] SOURCE EXECUTION")

require_file(LATEST_POINTER)

execution_dir = (
    ROOT
    / LATEST_POINTER.read_text(
        encoding="utf-8"
    ).strip()
)

stage1_audit_path = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1_AUDIT.json"
)

comparison_path = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1_CHARGE_COMPARISON.csv"
)

require_file(stage1_audit_path)
require_file(comparison_path)

print(f"execution_dir = {execution_dir}")
print(f"stage1_audit  = {stage1_audit_path}")
print(f"comparison    = {comparison_path}")


print("\n[2] UPSTREAM AUTHORIZATION")

stage1_audit = load_json(stage1_audit_path)

if stage1_audit.get("decision") != EXPECTED_AUDIT_DECISION:
    raise RuntimeError(
        "Unexpected Stage 1 audit decision.\n"
        f"Observed: {stage1_audit.get('decision')}"
    )

upstream_authorizations = stage1_audit.get(
    "authorizations",
    {},
)

if (
    upstream_authorizations.get(
        "candidate_charge_interpretation_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Candidate-charge interpretation is not authorized"
    )

if (
    upstream_authorizations.get(
        "RESP_stage2_authorized"
    )
    is not False
):
    raise RuntimeError(
        "Unexpected upstream Stage 2 authorization"
    )

print("upstream_audit_gate             = PASS")
print("interpretation_authorized_gate  = PASS")
print("stage2_remains_blocked_gate     = PASS")


print("\n[3] LOAD ATOM-LEVEL DATA")

with comparison_path.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    rows = list(csv.DictReader(handle))

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
    "node_type",
    "artificial_cap",
    "transfer_status",
    "RESP_atom_class",
    "candidate_equivalence_key",
    "equivalence_enforced",
    "CHELPG_charge_e",
    "RESP_stage1_charge_e",
}

missing_columns = required_columns - set(rows[0])

if missing_columns:
    raise RuntimeError(
        f"Missing required columns: {sorted(missing_columns)}"
    )

for row in rows:
    row["atom_index_0based"] = int(
        row["atom_index_0based"]
    )
    row["atom_index_1based"] = int(
        row["atom_index_1based"]
    )
    row["artificial_cap_bool"] = parse_bool(
        row["artificial_cap"]
    )
    row["equivalence_enforced_bool"] = parse_bool(
        row["equivalence_enforced"]
    )
    row["CHELPG_charge_e_float"] = float(
        row["CHELPG_charge_e"]
    )
    row["RESP_stage1_charge_e_float"] = float(
        row["RESP_stage1_charge_e"]
    )

indices = [
    row["atom_index_0based"]
    for row in rows
]

if indices != list(range(52)):
    raise RuntimeError(
        "Atom order is not exactly 0..51"
    )

print(f"atom_count = {len(rows)}")
print("atom_order_gate = PASS")


print("\n[4] ARTIFICIAL-CAP SUMMARY")

cap_rows = [
    row
    for row in rows
    if row["artificial_cap_bool"]
]

noncap_rows = [
    row
    for row in rows
    if not row["artificial_cap_bool"]
]

cap_resp = [
    row["RESP_stage1_charge_e_float"]
    for row in cap_rows
]

noncap_resp = [
    row["RESP_stage1_charge_e_float"]
    for row in noncap_rows
]

cap_summary = numeric_summary(cap_resp)
noncap_summary = numeric_summary(noncap_resp)

print(f"artificial_cap_count = {len(cap_rows)}")
print(f"noncap_count         = {len(noncap_rows)}")
print(f"artificial_cap_RESP_charge_sum_e = {sum(cap_resp):.16g}")
print(f"noncap_RESP_charge_sum_e         = {sum(noncap_resp):.16g}")

if cap_rows:
    print("\nArtificial-cap atoms:")

    for row in cap_rows:
        print(
            f"atom_0based={row['atom_index_0based']:>2} "
            f"atom_1based={row['atom_index_1based']:>2} "
            f"atom_id={row['atom_id']} "
            f"element={row['element']} "
            f"role={row['atom_role']} "
            f"transfer_status={row['transfer_status']} "
            f"RESP1={row['RESP_stage1_charge_e_float']: .6f}"
        )
else:
    print("(no artificial caps identified)")


print("\n[5] TRANSFER-STATUS GROUPS")

transfer_groups: dict[str, list[dict]] = defaultdict(list)

for row in rows:
    key = row["transfer_status"].strip() or "<blank>"
    transfer_groups[key].append(row)

transfer_summary = {}

for key in sorted(transfer_groups):
    selected = transfer_groups[key]
    charges = [
        row["RESP_stage1_charge_e_float"]
        for row in selected
    ]

    summary = numeric_summary(charges)
    transfer_summary[key] = summary

    print(f"\ntransfer_status = {key!r}")

    for name, value in summary.items():
        print(f"  {name} = {value}")


print("\n[6] ATOM-ROLE GROUPS")

role_groups: dict[str, list[dict]] = defaultdict(list)

for row in rows:
    key = row["atom_role"].strip() or "<blank>"
    role_groups[key].append(row)

role_summary = {}

for key in sorted(role_groups):
    selected = role_groups[key]
    charges = [
        row["RESP_stage1_charge_e_float"]
        for row in selected
    ]

    summary = numeric_summary(charges)
    role_summary[key] = summary

    print(f"\natom_role = {key!r}")

    for name, value in summary.items():
        print(f"  {name} = {value}")


print("\n[7] RESP ATOM CLASSES")

resp_class_groups: dict[str, list[dict]] = defaultdict(list)

for row in rows:
    key = row["RESP_atom_class"].strip() or "<blank>"
    resp_class_groups[key].append(row)

resp_class_summary = {}

for key in sorted(resp_class_groups):
    selected = resp_class_groups[key]
    charges = [
        row["RESP_stage1_charge_e_float"]
        for row in selected
    ]

    summary = numeric_summary(charges)
    resp_class_summary[key] = summary

    print(f"\nRESP_atom_class = {key!r}")

    for name, value in summary.items():
        print(f"  {name} = {value}")


print("\n[8] CANDIDATE EQUIVALENCE GROUPS")

candidate_groups: dict[str, list[dict]] = defaultdict(list)

for row in rows:
    key = (
        row["candidate_equivalence_key"].strip()
        or "<blank>"
    )
    candidate_groups[key].append(row)

candidate_group_summary = {}
nonsingleton_groups = {}
potentially_tight_groups = {}
high_dispersion_groups = {}

TIGHT_RANGE_THRESHOLD_E = 0.02
HIGH_RANGE_THRESHOLD_E = 0.10

for key in sorted(candidate_groups):
    selected = candidate_groups[key]
    charges = [
        row["RESP_stage1_charge_e_float"]
        for row in selected
    ]

    summary = numeric_summary(charges)
    summary["atom_indices_0based"] = [
        row["atom_index_0based"]
        for row in selected
    ]
    summary["elements"] = sorted(
        set(row["element"] for row in selected)
    )
    summary["enforced_count"] = sum(
        row["equivalence_enforced_bool"]
        for row in selected
    )

    candidate_group_summary[key] = summary

    if len(selected) > 1:
        nonsingleton_groups[key] = summary

        if summary["range_e"] <= TIGHT_RANGE_THRESHOLD_E:
            potentially_tight_groups[key] = summary

        if summary["range_e"] >= HIGH_RANGE_THRESHOLD_E:
            high_dispersion_groups[key] = summary

        print(f"\ncandidate_equivalence_key = {key!r}")

        for name, value in summary.items():
            print(f"  {name} = {value}")

        print("  atom_details:")

        for row in selected:
            print(
                f"    atom={row['atom_index_0based']:>2} "
                f"id={row['atom_id']} "
                f"element={row['element']} "
                f"role={row['atom_role']} "
                f"cap={row['artificial_cap_bool']} "
                f"transfer={row['transfer_status']} "
                f"RESP1={row['RESP_stage1_charge_e_float']: .6f}"
            )


print("\n[9] EQUIVALENCE-DISPERSION SUMMARY")

print(
    f"candidate_group_count = "
    f"{len(candidate_groups)}"
)
print(
    f"nonsingleton_candidate_group_count = "
    f"{len(nonsingleton_groups)}"
)
print(
    f"tight_nonsingleton_group_count_range_le_0.02e = "
    f"{len(potentially_tight_groups)}"
)
print(
    f"high_dispersion_group_count_range_ge_0.10e = "
    f"{len(high_dispersion_groups)}"
)

if high_dispersion_groups:
    print("\nHigh-dispersion candidate groups:")

    for key, summary in sorted(
        high_dispersion_groups.items()
    ):
        print(
            f"  {key}: count={summary['count']} "
            f"range_e={summary['range_e']:.6f} "
            f"indices={summary['atom_indices_0based']}"
        )


print("\n[10] SCIENTIFIC GATES")

finite_charge_gate = all(
    math.isfinite(
        row["RESP_stage1_charge_e_float"]
    )
    for row in rows
)

no_enforced_equivalence_gate = all(
    not row["equivalence_enforced_bool"]
    for row in rows
)

cap_classification_complete_gate = all(
    isinstance(
        row["artificial_cap_bool"],
        bool,
    )
    for row in rows
)

transfer_status_complete_gate = all(
    bool(row["transfer_status"].strip())
    for row in rows
)

candidate_key_complete_gate = all(
    bool(row["candidate_equivalence_key"].strip())
    for row in rows
)

gates = {
    "upstream_audit_gate": True,
    "atom_count_gate": len(rows) == 52,
    "atom_order_gate": indices == list(range(52)),
    "finite_charge_gate": finite_charge_gate,
    "no_enforced_equivalence_gate": (
        no_enforced_equivalence_gate
    ),
    "cap_classification_complete_gate": (
        cap_classification_complete_gate
    ),
    "transfer_status_complete_gate": (
        transfer_status_complete_gate
    ),
    "candidate_key_complete_gate": (
        candidate_key_complete_gate
    ),
}

for name, value in gates.items():
    print(
        f"{name} = "
        f"{'PASS' if value else 'FAIL'}"
    )

all_gates_pass = all(gates.values())


print("\n[11] WRITE OUTPUTS")

atom_audit_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1_TRANSFERABILITY.csv"
)

audit_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1_TRANSFERABILITY_AUDIT.json"
)

fieldnames = list(rows[0].keys())

with atom_audit_csv.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        extrasaction="ignore",
    )

    writer.writeheader()
    writer.writerows(rows)

decision = (
    "D038_E2_RESP_STAGE1_TRANSFERABILITY_AND_CAP_"
    "AUDIT_PASS_STAGE2_DESIGN_REVIEW_AUTHORIZED"
    if all_gates_pass
    else
    "D038_E2_RESP_STAGE1_TRANSFERABILITY_AUDIT_"
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
        "stage1_audit": str(
            stage1_audit_path.resolve()
        ),
        "stage1_audit_sha256": sha256(
            stage1_audit_path
        ),
        "stage1_comparison": str(
            comparison_path.resolve()
        ),
        "stage1_comparison_sha256": sha256(
            comparison_path
        ),
    },
    "artificial_caps": {
        "count": len(cap_rows),
        "RESP_charge_summary": cap_summary,
        "noncap_RESP_charge_summary": (
            noncap_summary
        ),
        "atom_indices_0based": [
            row["atom_index_0based"]
            for row in cap_rows
        ],
    },
    "transfer_status_summary": transfer_summary,
    "atom_role_summary": role_summary,
    "RESP_atom_class_summary": resp_class_summary,
    "candidate_equivalence_summary": {
        "all_groups": candidate_group_summary,
        "nonsingleton_groups": nonsingleton_groups,
        "tight_groups_range_le_0.02e": (
            potentially_tight_groups
        ),
        "high_dispersion_groups_range_ge_0.10e": (
            high_dispersion_groups
        ),
    },
    "thresholds": {
        "tight_group_range_e": (
            TIGHT_RANGE_THRESHOLD_E
        ),
        "high_dispersion_range_e": (
            HIGH_RANGE_THRESHOLD_E
        ),
    },
    "gates": gates,
    "authorizations": {
        "RESP_stage2_design_review_authorized": (
            all_gates_pass
        ),
        "RESP_stage2_execution_authorized": False,
        "charge_adoption_authorized": False,
        "force_field_adoption_authorized": False,
    },
    "outputs": {
        "atom_audit_csv": str(
            atom_audit_csv.resolve()
        ),
        "atom_audit_csv_sha256": sha256(
            atom_audit_csv
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

print(f"atom_audit_csv = {atom_audit_csv}")
print(
    f"atom_audit_csv_sha256 = "
    f"{sha256(atom_audit_csv)}"
)
print(f"audit_json = {audit_json}")
print(
    f"audit_json_sha256 = "
    f"{sha256(audit_json)}"
)


print("\n[12] DECISION")

print(f"decision = {decision}")
print(
    "RESP_stage2_design_review_authorized = "
    f"{all_gates_pass}"
)
print("RESP_stage2_execution_authorized = False")
print("charge_adoption_authorized = False")
print("force_field_adoption_authorized = False")
print("=" * 100)
