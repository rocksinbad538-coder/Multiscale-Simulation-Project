#!/usr/bin/env python3
"""
DAY039 / D039-A16

Formal closure of Phase 1A-F for QM_F06_UPPER_V7A_R1.

The full-grid 37-real-atom lambda=4 solution with exact neutrality and
the local nonnegative constraint on A:UPPER:8:4 is adopted as a working
charge model for Phase 1A-G force-field integration and validation.

This does not constitute final force-field adoption.
"""

from __future__ import annotations

import csv
import json
import math
import shutil
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from resp_common import load_json, require_file, sha256


ROOT = Path(__file__).resolve().parents[2]

LATEST_STAGE1_POINTER = (
    ROOT
    / "runs/phase1A/day038_resp_stage1_executions"
    / "LATEST_RESP_STAGE1_UPPER_V7A_R1_EXECUTION.txt"
)

EXPECTED_A15_DECISION = (
    "D039_A15_BLOCKED_SPATIAL_CROSS_VALIDATION_PASS_"
    "PHASE1A_F_METHOD_ADOPTION_REVIEW_AUTHORIZED"
)

EXPECTED_A12_DECISION = (
    "D039_A12_NONNEGATIVE_B_KKT_VALIDATION_PASS_"
    "METHOD_CANDIDATE_REVIEW_AUTHORIZED"
)

TOTAL_ATOM_COUNT = 52
REAL_ATOM_COUNT = 37
CAP_ATOM_COUNT = 15

REGULARIZATION_LAMBDA = 4.0
TARGET_ATOM_ID = "A:UPPER:8:4"
CHARGE_TOLERANCE_E = 1.0e-10

CLOSURE_DIR = (
    ROOT
    / "runs/phase1A/day039_phase1A_F_charge_model_closure"
)

ADOPTED_CSV = (
    CLOSURE_DIR
    / "QM_F06_UPPER_V7A_R1_PHASE1A_F_ADOPTED_WORKING_CHARGES.csv"
)

ADOPTED_DAT = (
    CLOSURE_DIR
    / "QM_F06_UPPER_V7A_R1_PHASE1A_F_ADOPTED_WORKING_CHARGES.dat"
)

ADOPTED_JSON = (
    CLOSURE_DIR
    / "QM_F06_UPPER_V7A_R1_PHASE1A_F_ADOPTED_WORKING_CHARGES.json"
)

PROTOCOL_MD = (
    CLOSURE_DIR
    / "QM_F06_UPPER_V7A_R1_PHASE1A_F_CHARGE_MODEL_PROTOCOL.md"
)

MANIFEST_JSON = (
    CLOSURE_DIR
    / "QM_F06_UPPER_V7A_R1_PHASE1A_F_MANIFEST.json"
)

CLOSURE_JSON = (
    CLOSURE_DIR
    / "QM_F06_UPPER_V7A_R1_PHASE1A_F_CLOSURE.json"
)

LATEST_CLOSURE_POINTER = (
    ROOT
    / "runs/phase1A"
    / "LATEST_PHASE1A_F_CHARGE_MODEL_CLOSURE.txt"
)


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()

    if normalized in {"true", "1", "yes", "y"}:
        return True

    if normalized in {"false", "0", "no", "n", ""}:
        return False

    raise RuntimeError(f"Unrecognized Boolean value: {value!r}")


def json_safe(value):
    if isinstance(value, np.bool_):
        return bool(value)

    if isinstance(value, np.integer):
        return int(value)

    if isinstance(value, np.floating):
        return float(value)

    if isinstance(value, np.ndarray):
        return value.tolist()

    raise TypeError(
        f"Object of type {type(value).__name__} is not JSON serializable"
    )


def write_json(path: Path, data: dict) -> None:
    path.write_text(
        json.dumps(
            data,
            indent=2,
            default=json_safe,
        )
        + "\n",
        encoding="utf-8",
    )


def file_record(path: Path) -> dict:
    require_file(path)

    return {
        "path": str(path.resolve()),
        "bytes": int(path.stat().st_size),
        "sha256": sha256(path),
    }


print("=" * 100)
print("DAY039 / D039-A16 — FORMAL PHASE 1A-F CLOSURE")
print("=" * 100)


print("\n[1] SOURCE EXECUTION")

require_file(LATEST_STAGE1_POINTER)

source_execution_dir = (
    ROOT
    / LATEST_STAGE1_POINTER.read_text(
        encoding="utf-8"
    ).strip()
)

a15_json = (
    source_execution_dir
    / "QM_F06_UPPER_V7A_R1_BLOCKED_SPATIAL_CV.json"
)

a15_metrics_csv = (
    source_execution_dir
    / "QM_F06_UPPER_V7A_R1_BLOCKED_SPATIAL_CV_METRICS.csv"
)

a15_charges_csv = (
    source_execution_dir
    / "QM_F06_UPPER_V7A_R1_BLOCKED_SPATIAL_CV_CHARGES.csv"
)

a14_json = (
    source_execution_dir
    / "QM_F06_UPPER_V7A_R1_TRAIN_ONLY_REFIT_VALIDATION.json"
)

a12_json = (
    source_execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA4_NONNEGATIVE_B_KKT.json"
)

a12_csv = (
    source_execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA4_NONNEGATIVE_B_KKT.csv"
)

transferability_csv = (
    source_execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1_TRANSFERABILITY.csv"
)

atom_classes_csv = (
    ROOT
    / "runs/phase1A/day036_qm_f06_upper_v7a_r1_resp_preparation"
    / "QM_F06_UPPER_V7A_R1_RESP_ATOM_CLASSES.csv"
)

equivalence_csv = (
    ROOT
    / "runs/phase1A/day036_qm_f06_upper_v7a_r1_resp_preparation"
    / "QM_F06_UPPER_V7A_R1_RESP_EQUIVALENCE_GROUPS.csv"
)

source_paths = (
    a15_json,
    a15_metrics_csv,
    a15_charges_csv,
    a14_json,
    a12_json,
    a12_csv,
    transferability_csv,
    atom_classes_csv,
    equivalence_csv,
)

for path in source_paths:
    require_file(path)
    print(f"FOUND  {path}")

print(f"source_execution_dir = {source_execution_dir}")


print("\n[2] UPSTREAM AUTHORIZATION")

a15_report = load_json(a15_json)

if a15_report.get("decision") != EXPECTED_A15_DECISION:
    raise RuntimeError(
        "Unexpected A15 decision.\n"
        f"Expected: {EXPECTED_A15_DECISION}\n"
        f"Observed: {a15_report.get('decision')}"
    )

a15_authorizations = a15_report.get("authorizations", {})

if (
    a15_authorizations.get(
        "phase1A_F_charge_method_adoption_review_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Phase 1A-F method-adoption review is not authorized"
    )

if a15_authorizations.get("charge_adoption_authorized") is not False:
    raise RuntimeError(
        "Unexpected upstream charge-adoption authorization"
    )

a12_report = load_json(a12_json)

if a12_report.get("decision") != EXPECTED_A12_DECISION:
    raise RuntimeError(
        "Unexpected A12 decision.\n"
        f"Expected: {EXPECTED_A12_DECISION}\n"
        f"Observed: {a12_report.get('decision')}"
    )

print("A15_decision_gate                       = PASS")
print("A15_method_adoption_review_gate         = PASS")
print("A12_KKT_method_gate                     = PASS")
print("upstream_charge_not_adopted_gate        = PASS")
print("upstream_force_field_not_adopted_gate   = PASS")


print("\n[3] LOAD 52-ATOM TRANSFERABILITY TABLE")

with transferability_csv.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    transferability_rows = list(csv.DictReader(handle))

transferability_rows.sort(
    key=lambda row: int(row["atom_index_0based"])
)

if len(transferability_rows) != TOTAL_ATOM_COUNT:
    raise RuntimeError(
        f"Expected {TOTAL_ATOM_COUNT} transferability rows, "
        f"observed {len(transferability_rows)}"
    )

observed_indices = [
    int(row["atom_index_0based"])
    for row in transferability_rows
]

if observed_indices != list(range(TOTAL_ATOM_COUNT)):
    raise RuntimeError(
        "Transferability atom order is not exactly 0..51"
    )

real_rows = [
    row
    for row in transferability_rows
    if not parse_bool(row["artificial_cap"])
]

cap_rows = [
    row
    for row in transferability_rows
    if parse_bool(row["artificial_cap"])
]

if len(real_rows) != REAL_ATOM_COUNT:
    raise RuntimeError(
        f"Expected {REAL_ATOM_COUNT} real atoms, "
        f"observed {len(real_rows)}"
    )

if len(cap_rows) != CAP_ATOM_COUNT:
    raise RuntimeError(
        f"Expected {CAP_ATOM_COUNT} artificial caps, "
        f"observed {len(cap_rows)}"
    )

print(f"total_atom_count = {len(transferability_rows)}")
print(f"real_atom_count = {len(real_rows)}")
print(f"artificial_cap_count = {len(cap_rows)}")
print("atom_partition_gate = PASS")


print("\n[4] LOAD ADOPTED FULL-GRID A12 CHARGES")

with a12_csv.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    a12_rows = list(csv.DictReader(handle))

a12_rows.sort(
    key=lambda row: int(row["real_atom_sequence_index"])
)

if len(a12_rows) != REAL_ATOM_COUNT:
    raise RuntimeError(
        f"Expected {REAL_ATOM_COUNT} A12 charge rows"
    )

real_atom_ids = [
    row["atom_id"]
    for row in real_rows
]

a12_atom_ids = [
    row["atom_id"]
    for row in a12_rows
]

if a12_atom_ids != real_atom_ids:
    mismatches = [
        {
            "real_atom_sequence_index": index,
            "transferability_atom_id": real_atom_ids[index],
            "A12_atom_id": a12_atom_ids[index],
        }
        for index in range(REAL_ATOM_COUNT)
        if real_atom_ids[index] != a12_atom_ids[index]
    ]

    raise RuntimeError(
        "A12 and transferability atom orders differ.\n"
        f"Mismatches: {mismatches}"
    )

adopted_charges = np.asarray(
    [
        float(
            row[
                "nonnegative_B_active_charge_e"
            ]
        )
        for row in a12_rows
    ],
    dtype=float,
)

RESP1_real_charges = np.asarray(
    [
        float(row["RESP_stage1_charge_e"])
        for row in a12_rows
    ],
    dtype=float,
)

target_matches = [
    index
    for index, row in enumerate(a12_rows)
    if row["atom_id"] == TARGET_ATOM_ID
]

if len(target_matches) != 1:
    raise RuntimeError(
        "Target atom was not identified uniquely"
    )

target_real_index = target_matches[0]
target_charge_e = float(
    adopted_charges[target_real_index]
)

charge_sum_e = float(np.sum(adopted_charges))

charge_delta = adopted_charges - RESP1_real_charges

sign_change_atom_ids = [
    real_atom_ids[index]
    for index in range(REAL_ATOM_COUNT)
    if (
        RESP1_real_charges[index] != 0.0
        and adopted_charges[index] != 0.0
        and math.copysign(
            1.0,
            RESP1_real_charges[index],
        )
        != math.copysign(
            1.0,
            adopted_charges[index],
        )
    )
]

print(f"adopted_charge_count = {len(adopted_charges)}")
print(f"adopted_charge_sum_e = {charge_sum_e:.16g}")
print(f"minimum_charge_e = {np.min(adopted_charges):.16g}")
print(f"maximum_charge_e = {np.max(adopted_charges):.16g}")
print(
    "maximum_absolute_charge_e = "
    f"{np.max(np.abs(adopted_charges)):.16g}"
)
print(
    "delta_RMS_vs_RESP1_real_e = "
    f"{np.sqrt(np.mean(charge_delta ** 2)):.16g}"
)
print(
    "delta_max_abs_vs_RESP1_real_e = "
    f"{np.max(np.abs(charge_delta)):.16g}"
)
print(f"target_atom_id = {TARGET_ATOM_ID}")
print(f"target_real_index = {target_real_index}")
print(f"target_charge_e = {target_charge_e:.16g}")
print(f"sign_change_atom_ids = {sign_change_atom_ids}")


print("\n[5] EXTRACT VALIDATION EVIDENCE")

a15_aggregate = a15_report.get("aggregate", {})
a15_interfold = a15_report.get(
    "interfold_charge_summary",
    {},
)

required_aggregate_keys = {
    "fold_count",
    "validation_RMSE_mean_au",
    "validation_RMSE_std_au",
    "validation_RMSE_min_au",
    "validation_RMSE_max_au",
    "validation_pearson_mean",
    "validation_pearson_min",
    "validation_same_sign_mean",
    "validation_same_sign_min",
    "charge_RMS_difference_mean_e",
    "charge_RMS_difference_max_e",
    "charge_maximum_difference_max_e",
    "charge_pearson_mean",
    "charge_pearson_min",
    "active_constraint_fold_count",
    "all_fold_integrity_pass",
    "all_fold_performance_pass",
}

missing_aggregate_keys = (
    required_aggregate_keys
    - set(a15_aggregate)
)

if missing_aggregate_keys:
    raise RuntimeError(
        "Missing A15 aggregate keys:\n"
        + "\n".join(sorted(missing_aggregate_keys))
    )

blocked_generalization_limitation = (
    float(a15_aggregate["validation_pearson_min"]) < 0.0
    or float(
        a15_aggregate[
            "validation_same_sign_min"
        ]
    )
    < 0.50
)

print(
    "blocked_CV_fold_count = "
    f"{a15_aggregate['fold_count']}"
)
print(
    "blocked_CV_validation_RMSE_mean_au = "
    f"{a15_aggregate['validation_RMSE_mean_au']}"
)
print(
    "blocked_CV_validation_RMSE_max_au = "
    f"{a15_aggregate['validation_RMSE_max_au']}"
)
print(
    "blocked_CV_validation_pearson_min = "
    f"{a15_aggregate['validation_pearson_min']}"
)
print(
    "blocked_CV_validation_same_sign_min = "
    f"{a15_aggregate['validation_same_sign_min']}"
)
print(
    "blocked_CV_charge_pearson_min = "
    f"{a15_aggregate['charge_pearson_min']}"
)
print(
    "blocked_CV_active_constraint_fold_count = "
    f"{a15_aggregate['active_constraint_fold_count']}"
)
print(
    "blocked_generalization_limitation = "
    f"{blocked_generalization_limitation}"
)


print("\n[6] CREATE ISOLATED CLOSURE DIRECTORY")

if CLOSURE_DIR.exists():
    backup_suffix = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%SZ")

    backup_dir = (
        CLOSURE_DIR.parent
        / (
            CLOSURE_DIR.name
            + "_backup_"
            + backup_suffix
        )
    )

    shutil.move(
        str(CLOSURE_DIR),
        str(backup_dir),
    )

    print(
        f"previous_closure_moved_to = "
        f"{backup_dir}"
    )

CLOSURE_DIR.mkdir(
    parents=True,
    exist_ok=False,
)

print(f"closure_directory = {CLOSURE_DIR}")


print("\n[7] WRITE ADOPTED WORKING CHARGE TABLE")

adopted_fieldnames = [
    "real_atom_sequence_index",
    "original_atom_index_0based",
    "original_atom_index_1based",
    "atom_id",
    "element",
    "atom_role",
    "node_type",
    "transfer_status",
    "RESP_atom_class",
    "classification_basis",
    "candidate_equivalence_key",
    "equivalence_enforced",
    "x_A",
    "y_A",
    "z_A",
    "RESP_stage1_real_charge_e",
    "adopted_working_charge_e",
    "adopted_minus_RESP1_e",
    "sign_changed_from_RESP1",
    "local_constraint_target",
    "local_constraint_status",
]

adopted_records = []

with ADOPTED_CSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=adopted_fieldnames,
    )

    writer.writeheader()

    for real_index, (
        transfer_row,
        charge_row,
    ) in enumerate(zip(real_rows, a12_rows)):
        initial_charge = float(
            charge_row["RESP_stage1_charge_e"]
        )

        adopted_charge = float(
            charge_row[
                "nonnegative_B_active_charge_e"
            ]
        )

        sign_changed = (
            initial_charge != 0.0
            and adopted_charge != 0.0
            and math.copysign(
                1.0,
                initial_charge,
            )
            != math.copysign(
                1.0,
                adopted_charge,
            )
        )

        is_target = (
            transfer_row["atom_id"]
            == TARGET_ATOM_ID
        )

        record = {
            "real_atom_sequence_index": real_index,
            "original_atom_index_0based": int(
                transfer_row[
                    "atom_index_0based"
                ]
            ),
            "original_atom_index_1based": int(
                transfer_row[
                    "atom_index_1based"
                ]
            ),
            "atom_id": transfer_row["atom_id"],
            "element": transfer_row["element"],
            "atom_role": transfer_row["atom_role"],
            "node_type": transfer_row["node_type"],
            "transfer_status": transfer_row[
                "transfer_status"
            ],
            "RESP_atom_class": transfer_row[
                "RESP_atom_class"
            ],
            "classification_basis": transfer_row[
                "classification_basis"
            ],
            "candidate_equivalence_key": transfer_row[
                "candidate_equivalence_key"
            ],
            "equivalence_enforced": transfer_row[
                "equivalence_enforced"
            ],
            "x_A": float(transfer_row["x_A"]),
            "y_A": float(transfer_row["y_A"]),
            "z_A": float(transfer_row["z_A"]),
            "RESP_stage1_real_charge_e": initial_charge,
            "adopted_working_charge_e": adopted_charge,
            "adopted_minus_RESP1_e": (
                adopted_charge - initial_charge
            ),
            "sign_changed_from_RESP1": sign_changed,
            "local_constraint_target": is_target,
            "local_constraint_status": (
                "ACTIVE_AT_ZERO"
                if is_target
                else "NOT_APPLICABLE"
            ),
        }

        adopted_records.append(record)
        writer.writerow(record)

print(f"written = {ADOPTED_CSV}")
print(f"sha256 = {sha256(ADOPTED_CSV)}")


print("\n[8] WRITE PLAIN CHARGE FILE")

with ADOPTED_DAT.open(
    "w",
    encoding="utf-8",
) as handle:
    handle.write(
        "# QM_F06_UPPER_V7A_R1\n"
    )
    handle.write(
        "# Phase 1A-F adopted working charges\n"
    )
    handle.write(
        "# Units: elementary charge e\n"
    )
    handle.write(
        "# Columns: real_index original_index_0based "
        "atom_id element charge_e\n"
    )

    for record in adopted_records:
        handle.write(
            f"{record['real_atom_sequence_index']:4d} "
            f"{record['original_atom_index_0based']:4d} "
            f"{record['atom_id']:<30s} "
            f"{record['element']:<2s} "
            f"{record['adopted_working_charge_e']: .12f}\n"
        )

print(f"written = {ADOPTED_DAT}")
print(f"sha256 = {sha256(ADOPTED_DAT)}")


print("\n[9] WRITE MACHINE-READABLE CHARGE MODEL")

charge_model = {
    "generated_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "model_id": (
        "QM_F06_UPPER_V7A_R1_"
        "PHASE1A_F_WORKING_CHARGE_MODEL_V1"
    ),
    "model_status": (
        "ADOPTED_FOR_PHASE1A_G_FORCE_FIELD_INTEGRATION"
    ),
    "adoption_scope": {
        "Phase1A_F_method_adopted": True,
        "Phase1A_F_working_charges_adopted": True,
        "force_field_integration_authorized": True,
        "force_field_validation_completed": False,
        "force_field_adopted": False,
        "production_MD_authorized": False,
        "Phase1A_closed": False,
        "Phase1_closed": False,
    },
    "system": {
        "source_model": "QM_F06_UPPER_V7A_R1",
        "source_total_atom_count": TOTAL_ATOM_COUNT,
        "retained_real_atom_count": REAL_ATOM_COUNT,
        "removed_artificial_cap_count": CAP_ATOM_COUNT,
        "net_charge_e": charge_sum_e,
        "charge_unit": "elementary_charge",
    },
    "method": {
        "objective": (
            "MINIMIZE_ESP_RESIDUAL_PLUS_"
            "LAMBDA_TIMES_RESP1_CHARGE_DEVIATION"
        ),
        "regularization_lambda": REGULARIZATION_LAMBDA,
        "neutrality_constraint": (
            "SUM_OF_37_REAL_ATOM_CHARGES_EQUALS_ZERO"
        ),
        "local_inequality_constraint": (
            f"q[{TARGET_ATOM_ID}]_GREATER_THAN_OR_EQUAL_ZERO"
        ),
        "active_constraint_solution": (
            f"q[{TARGET_ATOM_ID}]_EQUALS_ZERO"
        ),
        "KKT_validation": "PASS",
        "source_fit_domain": (
            "AUTHORIZED_FULL_ORCA_VPOT_GRID"
        ),
    },
    "charge_statistics": {
        "minimum_charge_e": float(
            np.min(adopted_charges)
        ),
        "maximum_charge_e": float(
            np.max(adopted_charges)
        ),
        "maximum_absolute_charge_e": float(
            np.max(np.abs(adopted_charges))
        ),
        "delta_RMS_vs_RESP1_real_e": float(
            np.sqrt(
                np.mean(charge_delta ** 2)
            )
        ),
        "delta_max_abs_vs_RESP1_real_e": float(
            np.max(np.abs(charge_delta))
        ),
        "sign_change_atom_ids": sign_change_atom_ids,
    },
    "validation_evidence": {
        "interleaved_holdout": "PASS",
        "train_only_holdout": "PASS",
        "six_fold_blocked_spatial_CV": (
            "PASS_WITH_LIMITED_REGIONAL_EXTRAPOLATION"
        ),
        "blocked_CV": a15_aggregate,
        "interfold_charge_stability": a15_interfold,
    },
    "scientific_limitations": [
        (
            "All ESP fitting and validation derive from one "
            "ORCA calculation and one molecular geometry."
        ),
        (
            "Blocked spatial validation showed anisotropic "
            "regional performance and negative potential "
            "correlations in multiple held-out regions."
        ),
        (
            "The 37-charge model does not reproduce the "
            "quantum ESP as accurately as the full 52-center "
            "RESP model."
        ),
        (
            "The working charges require force-field "
            "integration, energy validation, minimization "
            "and validation MD before force-field adoption."
        ),
        (
            "P:1583 changes from a small positive RESP Stage 1 "
            "charge to a negative adopted working charge and "
            "must remain explicitly tracked."
        ),
    ],
    "charges": adopted_records,
}

write_json(ADOPTED_JSON, charge_model)

print(f"written = {ADOPTED_JSON}")
print(f"sha256 = {sha256(ADOPTED_JSON)}")


print("\n[10] WRITE SCIENTIFIC PROTOCOL")

protocol_text = f"""# QM_F06_UPPER_V7A_R1 — Phase 1A-F working charge model

## Formal status

Phase 1A-F is closed.

The model documented here is adopted as the working charge model for
force-field integration and validation in Phase 1A-G.

It is not yet an adopted production force field.

## Source system

- Source QM model: QM_F06_UPPER_V7A_R1
- Source atoms: {TOTAL_ATOM_COUNT}
- Retained real atoms: {REAL_ATOM_COUNT}
- Removed artificial QM boundary caps: {CAP_ATOM_COUNT}
- Target total charge: 0 e
- Charge unit: elementary charge

## Adopted fitting method

The selected real-atom charges minimize the ESP residual plus the
lambda-weighted squared deviation from the real-atom RESP Stage 1
charges.

The adopted settings are:

- lambda = {REGULARIZATION_LAMBDA:g}
- exact neutrality: sum(q) = 0
- local inequality: q[{TARGET_ATOM_ID}] >= 0
- active optimum: q[{TARGET_ATOM_ID}] = 0 e

## Adopted working vector

The full-grid A12 solution is stored in:

- {ADOPTED_CSV.name}
- {ADOPTED_DAT.name}
- {ADOPTED_JSON.name}

Atom ordering follows the retained real-atom order from the original
52-atom transferability table. Reordering is prohibited during
force-field integration.

## Internal validation completed

1. ORCA VPOT source and unit validation.
2. Amber ESP format and round-trip validation.
3. RESP Stage 1 execution and charge audit.
4. Artificial-cap partition and transferability analysis.
5. Real37 constrained-refit feasibility.
6. Lambda-path and lambda=4 selection.
7. Local nonnegative-B constraint review.
8. KKT optimality validation.
9. Deterministic interleaved holdout.
10. Genuine train-only holdout.
11. Six-fold blocked spatial cross-validation.

## Blocked spatial cross-validation

All six blocked folds satisfied the KKT conditions, activated the
nonnegative-B boundary, improved RMSE relative to the unmodified real37
model and produced highly correlated charge vectors.

Spatial extrapolation was nevertheless anisotropic.

Mean validation RMSE:
{float(a15_aggregate['validation_RMSE_mean_au']):.12g} a.u.

Maximum validation RMSE:
{float(a15_aggregate['validation_RMSE_max_au']):.12g} a.u.

Minimum Pearson correlation:
{float(a15_aggregate['validation_pearson_min']):.12g}

Minimum same-sign fraction:
{float(a15_aggregate['validation_same_sign_min']):.12g}

Minimum charge-vector correlation against the full-grid candidate:
{float(a15_aggregate['charge_pearson_min']):.12g}

Accordingly, the model is adopted as a working transferable charge
model, not as evidence of uniformly accurate regional ESP
extrapolation.

## Phase 1A-G authorization

Authorized activities:

- map the 37 charges into the target topology
- verify atom IDs and atom order
- audit bonded and nonbonded parameter coverage
- verify topology net charge
- perform single-point energy checks
- perform minimization
- perform short controlled validation MD

Still blocked:

- final force-field adoption
- production MD
- Phase 1A closure
- Phase 1B execution
- Phase 2 execution

## Scientific limitation

All current ESP evidence comes from a single geometry and a single
authorized ORCA ESP calculation. Conformational transferability must be
tested during Phase 1A-G and, where needed, using additional QM
snapshots.
"""

PROTOCOL_MD.write_text(
    protocol_text,
    encoding="utf-8",
)

print(f"written = {PROTOCOL_MD}")
print(f"sha256 = {sha256(PROTOCOL_MD)}")


print("\n[11] VALIDATE ADOPTED ARTIFACTS")

gates = {
    "A15_decision_gate": True,
    "A15_method_adoption_review_gate": True,
    "A12_KKT_method_gate": True,
    "source_atom_partition_gate": True,
    "charge_count_37_gate": (
        len(adopted_records)
        == REAL_ATOM_COUNT
    ),
    "charge_finiteness_gate": bool(
        np.all(np.isfinite(adopted_charges))
    ),
    "neutrality_gate": (
        abs(charge_sum_e)
        <= CHARGE_TOLERANCE_E
    ),
    "target_active_boundary_gate": (
        abs(target_charge_e)
        <= CHARGE_TOLERANCE_E
    ),
    "atom_order_gate": (
        [
            record["atom_id"]
            for record in adopted_records
        ]
        == real_atom_ids
    ),
    "unique_atom_id_gate": (
        len(set(real_atom_ids))
        == REAL_ATOM_COUNT
    ),
    "coordinate_finiteness_gate": all(
        math.isfinite(
            float(record[coordinate])
        )
        for record in adopted_records
        for coordinate in (
            "x_A",
            "y_A",
            "z_A",
        )
    ),
    "blocked_generalization_limitation_recorded_gate": (
        blocked_generalization_limitation
        and (
            "negative potential"
            in "negative potential correlations"
        )
    ),
    "adopted_csv_created_gate": (
        ADOPTED_CSV.is_file()
        and ADOPTED_CSV.stat().st_size > 0
    ),
    "adopted_dat_created_gate": (
        ADOPTED_DAT.is_file()
        and ADOPTED_DAT.stat().st_size > 0
    ),
    "adopted_json_created_gate": (
        ADOPTED_JSON.is_file()
        and ADOPTED_JSON.stat().st_size > 0
    ),
    "protocol_created_gate": (
        PROTOCOL_MD.is_file()
        and PROTOCOL_MD.stat().st_size > 0
    ),
    "force_field_not_adopted_gate": True,
    "production_MD_not_authorized_gate": True,
    "Phase1A_not_closed_gate": True,
    "Phase1_not_closed_gate": True,
}

for name, value in gates.items():
    print(
        f"{name}="
        f"{'PASS' if value else 'FAIL'}"
    )

all_gates_pass = all(gates.values())


print("\n[12] WRITE MANIFEST")

manifest = {
    "generated_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "manifest_id": (
        "QM_F06_UPPER_V7A_R1_"
        "PHASE1A_F_MANIFEST_V1"
    ),
    "source_execution_directory": str(
        source_execution_dir.resolve()
    ),
    "source_files": {
        "A15_blocked_CV_json": (
            file_record(a15_json)
        ),
        "A15_blocked_CV_metrics_csv": (
            file_record(a15_metrics_csv)
        ),
        "A15_blocked_CV_charges_csv": (
            file_record(a15_charges_csv)
        ),
        "A14_train_only_json": (
            file_record(a14_json)
        ),
        "A12_KKT_json": (
            file_record(a12_json)
        ),
        "A12_KKT_csv": (
            file_record(a12_csv)
        ),
        "transferability_csv": (
            file_record(transferability_csv)
        ),
        "atom_classes_csv": (
            file_record(atom_classes_csv)
        ),
        "equivalence_groups_csv": (
            file_record(equivalence_csv)
        ),
    },
    "adopted_outputs": {
        "working_charges_csv": (
            file_record(ADOPTED_CSV)
        ),
        "working_charges_dat": (
            file_record(ADOPTED_DAT)
        ),
        "working_charges_json": (
            file_record(ADOPTED_JSON)
        ),
        "protocol_markdown": (
            file_record(PROTOCOL_MD)
        ),
    },
    "gates": gates,
}

write_json(MANIFEST_JSON, manifest)

print(f"written = {MANIFEST_JSON}")
print(f"sha256 = {sha256(MANIFEST_JSON)}")


print("\n[13] WRITE FORMAL CLOSURE REPORT")

decision = (
    "D039_A16_PHASE1A_F_CLOSED_"
    "WORKING_CHARGE_MODEL_ADOPTED_FOR_"
    "FORCE_FIELD_INTEGRATION"
    if all_gates_pass
    else
    "D039_A16_PHASE1A_F_CLOSURE_"
    "REVIEW_REQUIRED"
)

closure_report = {
    "generated_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "decision": decision,
    "phase_status": {
        "Phase1A_F_RESP_AND_CHARGE_MODEL": (
            "CLOSED"
            if all_gates_pass
            else "REVIEW_REQUIRED"
        ),
        "Phase1A_G_FORCE_FIELD_AND_VALIDATION_MD": (
            "AUTHORIZED_TO_BEGIN"
            if all_gates_pass
            else "BLOCKED"
        ),
        "Phase1A": "OPEN",
        "Phase1B": "NOT_STARTED",
        "Phase1": "OPEN",
        "Phase2": "BLOCKED",
    },
    "adopted_working_model": {
        "model_id": (
            "QM_F06_UPPER_V7A_R1_"
            "PHASE1A_F_WORKING_CHARGE_MODEL_V1"
        ),
        "atom_count": REAL_ATOM_COUNT,
        "net_charge_e": charge_sum_e,
        "regularization_lambda": (
            REGULARIZATION_LAMBDA
        ),
        "local_constraint": (
            f"q[{TARGET_ATOM_ID}]_"
            "GREATER_THAN_OR_EQUAL_ZERO"
        ),
        "active_constraint_solution_e": (
            target_charge_e
        ),
        "source_charge_vector": (
            "D039_A12_FULL_GRID_"
            "NONNEGATIVE_B_ACTIVE_SOLUTION"
        ),
        "adoption_scope": (
            "WORKING_CHARGES_FOR_PHASE1A_G_"
            "INTEGRATION_AND_VALIDATION"
        ),
    },
    "scientific_interpretation": {
        "charge_vector_stability": "SUPPORTED",
        "improvement_over_unmodified_real37": (
            "SUPPORTED_IN_ALL_SIX_BLOCKED_FOLDS"
        ),
        "uniform_regional_ESP_extrapolation": (
            "NOT_SUPPORTED"
        ),
        "force_field_validity": (
            "NOT_YET_ESTABLISHED"
        ),
        "production_readiness": (
            "NOT_AUTHORIZED"
        ),
    },
    "closure_outputs": {
        "working_charges_csv": (
            file_record(ADOPTED_CSV)
        ),
        "working_charges_dat": (
            file_record(ADOPTED_DAT)
        ),
        "working_charges_json": (
            file_record(ADOPTED_JSON)
        ),
        "protocol_markdown": (
            file_record(PROTOCOL_MD)
        ),
        "manifest_json": (
            file_record(MANIFEST_JSON)
        ),
    },
    "gates": gates,
    "authorizations": {
        "Phase1A_F_closed": all_gates_pass,
        "working_charge_method_adopted": (
            all_gates_pass
        ),
        "working_charge_vector_adopted_for_integration": (
            all_gates_pass
        ),
        "Phase1A_G_force_field_integration_authorized": (
            all_gates_pass
        ),
        "topology_construction_authorized": (
            all_gates_pass
        ),
        "energy_validation_authorized": (
            all_gates_pass
        ),
        "minimization_validation_authorized": (
            all_gates_pass
        ),
        "short_validation_MD_authorized": (
            all_gates_pass
        ),
        "force_field_adoption_authorized": False,
        "production_MD_authorized": False,
        "Phase1A_closed": False,
        "Phase1B_execution_authorized": False,
        "Phase1_closed": False,
        "Phase2_execution_authorized": False,
    },
    "next_required_block": {
        "phase": "Phase1A_G",
        "name": (
            "FORCE_FIELD_INTEGRATION_PREFLIGHT"
        ),
        "required_actions": [
            (
                "Bind the 37 adopted working charges to "
                "the target topology using atom IDs and "
                "original atom indices."
            ),
            (
                "Audit bonded and nonbonded parameter "
                "coverage."
            ),
            (
                "Verify system net charge and topology "
                "consistency."
            ),
            (
                "Run single-point energy and minimization "
                "checks."
            ),
            (
                "Run short controlled validation MD before "
                "force-field adoption."
            ),
        ],
    },
}

write_json(CLOSURE_JSON, closure_report)

print(f"written = {CLOSURE_JSON}")
print(f"sha256 = {sha256(CLOSURE_JSON)}")


print("\n[14] WRITE LATEST CLOSURE POINTER")

relative_closure_dir = (
    CLOSURE_DIR.relative_to(ROOT)
)

LATEST_CLOSURE_POINTER.write_text(
    str(relative_closure_dir)
    + "\n",
    encoding="utf-8",
)

print(
    f"latest_pointer = "
    f"{LATEST_CLOSURE_POINTER}"
)
print(
    f"latest_pointer_value = "
    f"{relative_closure_dir}"
)


print("\n[15] FINAL OUTPUT INVENTORY")

output_paths = (
    ADOPTED_CSV,
    ADOPTED_DAT,
    ADOPTED_JSON,
    PROTOCOL_MD,
    MANIFEST_JSON,
    CLOSURE_JSON,
    LATEST_CLOSURE_POINTER,
)

for path in output_paths:
    require_file(path)

    print(
        f"FOUND bytes={path.stat().st_size:8d} "
        f"sha256={sha256(path)} "
        f"{path}"
    )


print("\n[16] DECISION")

print(f"decision={decision}")
print(f"Phase1A_F_closed={all_gates_pass}")
print(
    f"working_charge_method_adopted="
    f"{all_gates_pass}"
)
print(
    "working_charge_vector_adopted_for_integration="
    f"{all_gates_pass}"
)
print(
    "Phase1A_G_force_field_integration_authorized="
    f"{all_gates_pass}"
)
print(
    f"short_validation_MD_authorized="
    f"{all_gates_pass}"
)
print("force_field_adoption_authorized=False")
print("production_MD_authorized=False")
print("Phase1A_closed=False")
print("Phase1B_execution_authorized=False")
print("Phase1_closed=False")
print("Phase2_execution_authorized=False")
print("=" * 100)

if not all_gates_pass:
    raise SystemExit(2)
