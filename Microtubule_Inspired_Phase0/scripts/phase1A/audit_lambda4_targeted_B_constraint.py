#!/usr/bin/env python3
"""
DAY039 / D039-A11

Targeted chemical-constraint audit for the lambda=4 constrained
37-real-atom electrostatic refit.

The atom A:UPPER:8:4 is the only chemically sensitive B sign inversion
that remains under the focused lambda=4 candidate.

Compared candidates
-------------------
1. Unconstrained lambda=4.
2. lambda=4 with q(A:UPPER:8:4) fixed to 0.0 e.
3. lambda=4 with q(A:UPPER:8:4) fixed to +0.05 e.
4. lambda=4 with q(A:UPPER:8:4) fixed to its RESP Stage 1 value.

All candidates retain exact total neutrality.

No candidate is adopted.
RESP Stage 2 remains blocked.
"""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from resp_common import (
    load_json,
    read_orca_vpot,
    require_file,
    sha256,
)


ROOT = Path(__file__).resolve().parents[2]

LATEST_POINTER = (
    ROOT
    / "runs/phase1A/day038_resp_stage1_executions"
    / "LATEST_RESP_STAGE1_UPPER_V7A_R1_EXECUTION.txt"
)

VPOT = (
    ROOT
    / "runs/phase1A/day036_qm_f06_upper_v7a_r1_esp_executions"
    / "esp_upper_v7a_r1_20260731T174832Z"
    / "esp_upper_v7a_r1.vpot"
)

EXPECTED_A10_DECISION = (
    "D039_A10_LAMBDA4_FOCUSED_CHEMICAL_AUDIT_PASS_"
    "FINAL_METHOD_REVIEW_AUTHORIZED"
)

AUTHORIZED_VPOT_SHA256 = (
    "73df47796c29d0b2a88a03d83efcb723"
    "d4f6e583d2e1789e1ea1a43d82c1064d"
)

REGULARIZATION_LAMBDA = 4.0
TARGET_ATOM_ID = "A:UPPER:8:4"
PERSISTENT_N_ATOM_ID = "P:1583"

CHARGE_TOLERANCE_E = 1.0e-10

TARGET_POLICIES = (
    ("UNCONSTRAINED_LAMBDA4", None),
    ("TARGET_B_FIXED_ZERO", 0.0),
    ("TARGET_B_FIXED_PLUS_0P05", 0.05),
    ("TARGET_B_FIXED_RESP1", "RESP1"),
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


def build_coulomb_matrix(
    grid_xyz_bohr: np.ndarray,
    atom_xyz_bohr: np.ndarray,
    chunk_size: int = 2000,
) -> np.ndarray:
    matrix = np.empty(
        (
            len(grid_xyz_bohr),
            len(atom_xyz_bohr),
        ),
        dtype=float,
    )

    for start in range(
        0,
        len(grid_xyz_bohr),
        chunk_size,
    ):
        stop = min(
            start + chunk_size,
            len(grid_xyz_bohr),
        )

        displacement = (
            grid_xyz_bohr[start:stop, None, :]
            - atom_xyz_bohr[None, :, :]
        )

        distances = np.linalg.norm(
            displacement,
            axis=2,
        )

        if np.any(distances <= 0.0):
            raise RuntimeError(
                "An ESP point coincides with an atomic center"
            )

        matrix[start:stop, :] = 1.0 / distances

    return matrix


def affine_constraint_solution(
    design_matrix: np.ndarray,
    target: np.ndarray,
    reference_charges: np.ndarray,
    regularization_lambda: float,
    constraint_matrix: np.ndarray,
    constraint_values: np.ndarray,
) -> np.ndarray:
    """
    Solve:

        min ||Aq-b||² + lambda ||q-q0||²
        subject to Cq=d

    through an affine null-space parameterization.
    """

    particular_solution, _, _, _ = np.linalg.lstsq(
        constraint_matrix,
        constraint_values,
        rcond=None,
    )

    constraint_residual = (
        constraint_matrix @ particular_solution
        - constraint_values
    )

    if np.max(
        np.abs(
            constraint_residual
        )
    ) > 1.0e-10:
        raise RuntimeError(
            "Could not construct a valid particular solution"
        )

    _, singular_values, vh = np.linalg.svd(
        constraint_matrix,
        full_matrices=True,
    )

    if singular_values.size:
        rank_tolerance = (
            singular_values[0]
            * max(constraint_matrix.shape)
            * np.finfo(float).eps
        )

        constraint_rank = int(
            np.sum(
                singular_values > rank_tolerance
            )
        )
    else:
        constraint_rank = 0

    null_basis = vh[
        constraint_rank:,
        :
    ].T

    reduced_matrix = (
        design_matrix @ null_basis
    )

    shifted_target = (
        target
        - design_matrix @ particular_solution
    )

    sqrt_lambda = math.sqrt(
        regularization_lambda
    )

    augmented_matrix = np.vstack(
        (
            reduced_matrix,
            sqrt_lambda * null_basis,
        )
    )

    augmented_target = np.concatenate(
        (
            shifted_target,
            sqrt_lambda
            * (
                reference_charges
                - particular_solution
            ),
        )
    )

    reduced_solution, _, _, _ = np.linalg.lstsq(
        augmented_matrix,
        augmented_target,
        rcond=None,
    )

    charges = (
        particular_solution
        + null_basis @ reduced_solution
    )

    final_constraint_residual = (
        constraint_matrix @ charges
        - constraint_values
    )

    if np.max(
        np.abs(
            final_constraint_residual
        )
    ) > 1.0e-9:
        raise RuntimeError(
            "Final constrained solution violates its constraints"
        )

    return charges


def electrostatic_metrics(
    reference: np.ndarray,
    candidate: np.ndarray,
) -> dict:
    residual = candidate - reference

    reference_centered = (
        reference - np.mean(reference)
    )

    candidate_centered = (
        candidate - np.mean(candidate)
    )

    denominator = math.sqrt(
        float(
            np.sum(reference_centered ** 2)
            * np.sum(candidate_centered ** 2)
        )
    )

    pearson_r = (
        float(
            np.sum(
                reference_centered
                * candidate_centered
            )
            / denominator
        )
        if denominator > 0.0
        else float("nan")
    )

    return {
        "RMSE_au": float(
            np.sqrt(
                np.mean(residual ** 2)
            )
        ),
        "MAE_au": float(
            np.mean(
                np.abs(residual)
            )
        ),
        "maximum_absolute_error_au": float(
            np.max(
                np.abs(residual)
            )
        ),
        "residual_mean_au": float(
            np.mean(residual)
        ),
        "pearson_r": pearson_r,
        "same_sign_fraction": float(
            np.mean(
                np.sign(reference)
                == np.sign(candidate)
            )
        ),
    }


def charge_metrics(
    reference: np.ndarray,
    candidate: np.ndarray,
    atom_rows: list[dict],
) -> dict:
    delta = candidate - reference

    sign_change_indices = [
        index
        for index in range(
            len(reference)
        )
        if (
            reference[index] != 0.0
            and candidate[index] != 0.0
            and math.copysign(
                1.0,
                reference[index],
            )
            != math.copysign(
                1.0,
                candidate[index],
            )
        )
    ]

    return {
        "charge_sum_e": float(
            np.sum(candidate)
        ),
        "minimum_charge_e": float(
            np.min(candidate)
        ),
        "maximum_charge_e": float(
            np.max(candidate)
        ),
        "maximum_absolute_charge_e": float(
            np.max(
                np.abs(candidate)
            )
        ),
        "delta_MAE_e": float(
            np.mean(
                np.abs(delta)
            )
        ),
        "delta_RMS_e": float(
            np.sqrt(
                np.mean(delta ** 2)
            )
        ),
        "delta_max_abs_e": float(
            np.max(
                np.abs(delta)
            )
        ),
        "sign_change_count": len(
            sign_change_indices
        ),
        "sign_change_real_indices": [
            int(index)
            for index in sign_change_indices
        ],
        "sign_change_atom_ids": [
            atom_rows[index]["atom_id"]
            for index in sign_change_indices
        ],
    }


print("=" * 100)
print("DAY039 / D039-A11 — TARGETED B-CONSTRAINT AUDIT")
print("=" * 100)


print("\n[1] SOURCE EXECUTION")

for path in (
    LATEST_POINTER,
    VPOT,
):
    require_file(path)

if sha256(VPOT) != AUTHORIZED_VPOT_SHA256:
    raise RuntimeError(
        "Authorized VPOT hash mismatch"
    )

execution_dir = (
    ROOT
    / LATEST_POINTER.read_text(
        encoding="utf-8"
    ).strip()
)

a10_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA4_FOCUSED_CHEMICAL_AUDIT.json"
)

transferability_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1_TRANSFERABILITY.csv"
)

for path in (
    a10_json,
    transferability_csv,
):
    require_file(path)
    print(f"FOUND  {path}")


print("\n[2] UPSTREAM AUTHORIZATION")

a10_report = load_json(
    a10_json
)

if (
    a10_report.get("decision")
    != EXPECTED_A10_DECISION
):
    raise RuntimeError(
        "Unexpected A10 decision.\n"
        f"Observed: {a10_report.get('decision')}"
    )

authorizations = a10_report.get(
    "authorizations",
    {},
)

if (
    authorizations.get(
        "lambda4_final_method_scientific_review_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Targeted chemical-constraint review is not authorized"
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

print("A10_decision_gate                     = PASS")
print("targeted_constraint_review_gate       = PASS")
print("lambda_adoption_blocked_gate          = PASS")
print("charge_adoption_blocked_gate          = PASS")


print("\n[3] LOAD VPOT AND RESP DATA")

vpot = read_orca_vpot(
    VPOT
)

atom_xyz_bohr = np.asarray(
    vpot.atom_coordinates_bohr,
    dtype=float,
)

grid_xyz_bohr = np.asarray(
    vpot.grid_coordinates_bohr,
    dtype=float,
)

quantum_esp = np.asarray(
    vpot.grid_potential_au,
    dtype=float,
).reshape(-1)

with transferability_csv.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    all_rows = list(
        csv.DictReader(handle)
    )

all_rows.sort(
    key=lambda row: int(
        row["atom_index_0based"]
    )
)

real_mask = np.asarray(
    [
        not parse_bool(
            row["artificial_cap"]
        )
        for row in all_rows
    ],
    dtype=bool,
)

real_rows = [
    row
    for row, retained in zip(
        all_rows,
        real_mask,
    )
    if retained
]

real_xyz_bohr = atom_xyz_bohr[
    real_mask
]

q0 = np.asarray(
    [
        float(
            row["RESP_stage1_charge_e"]
        )
        for row in real_rows
    ],
    dtype=float,
)

if len(real_rows) != 37:
    raise RuntimeError(
        "Expected 37 retained atoms"
    )

atom_index_by_id = {
    row["atom_id"]: index
    for index, row in enumerate(
        real_rows
    )
}

if TARGET_ATOM_ID not in atom_index_by_id:
    raise RuntimeError(
        f"Target atom not found: {TARGET_ATOM_ID}"
    )

if PERSISTENT_N_ATOM_ID not in atom_index_by_id:
    raise RuntimeError(
        f"Persistent N atom not found: {PERSISTENT_N_ATOM_ID}"
    )

target_index = atom_index_by_id[
    TARGET_ATOM_ID
]

persistent_n_index = atom_index_by_id[
    PERSISTENT_N_ATOM_ID
]

target_RESP1_charge = q0[
    target_index
]

print(f"real_atom_count = {len(real_rows)}")
print(f"grid_point_count = {len(grid_xyz_bohr)}")
print(f"target_atom_id = {TARGET_ATOM_ID}")
print(f"target_real_index = {target_index}")
print(
    f"target_RESP1_charge_e = "
    f"{target_RESP1_charge:.16g}"
)
print(
    f"persistent_N_RESP1_charge_e = "
    f"{q0[persistent_n_index]:.16g}"
)


print("\n[4] BUILD COULOMB MATRIX")

design_matrix = build_coulomb_matrix(
    grid_xyz_bohr,
    real_xyz_bohr,
)

print(
    f"design_matrix_shape = "
    f"{design_matrix.shape}"
)


print("\n[5] SOLVE TARGETED CONSTRAINT POLICIES")

candidate_vectors = {}
candidate_records = []

for policy_name, policy_value in (
    TARGET_POLICIES
):
    if policy_value is None:
        constraint_matrix = np.ones(
            (
                1,
                len(real_rows),
            ),
            dtype=float,
        )

        constraint_values = np.asarray(
            [0.0],
            dtype=float,
        )

        target_fixed_charge = None

    else:
        if policy_value == "RESP1":
            fixed_charge = float(
                target_RESP1_charge
            )
        else:
            fixed_charge = float(
                policy_value
            )

        constraint_matrix = np.zeros(
            (
                2,
                len(real_rows),
            ),
            dtype=float,
        )

        constraint_matrix[0, :] = 1.0
        constraint_matrix[
            1,
            target_index,
        ] = 1.0

        constraint_values = np.asarray(
            [
                0.0,
                fixed_charge,
            ],
            dtype=float,
        )

        target_fixed_charge = fixed_charge

    candidate = affine_constraint_solution(
        design_matrix,
        quantum_esp,
        q0,
        REGULARIZATION_LAMBDA,
        constraint_matrix,
        constraint_values,
    )

    predicted = (
        design_matrix @ candidate
    )

    esp_metrics = electrostatic_metrics(
        quantum_esp,
        predicted,
    )

    q_metrics = charge_metrics(
        q0,
        candidate,
        real_rows,
    )

    record = {
        "policy_name": policy_name,
        "regularization_lambda": (
            REGULARIZATION_LAMBDA
        ),
        "target_atom_id": (
            TARGET_ATOM_ID
        ),
        "target_fixed_charge_e": (
            target_fixed_charge
        ),
        "target_candidate_charge_e": float(
            candidate[target_index]
        ),
        "persistent_N_candidate_charge_e": float(
            candidate[persistent_n_index]
        ),
        "electrostatic": esp_metrics,
        "charges": q_metrics,
    }

    candidate_vectors[
        policy_name
    ] = candidate

    candidate_records.append(
        record
    )

    print(f"\npolicy={policy_name}")

    print(
        f"  target_fixed_charge_e = "
        f"{target_fixed_charge}"
    )
    print(
        f"  target_candidate_charge_e = "
        f"{candidate[target_index]:.16g}"
    )
    print(
        f"  persistent_N_candidate_charge_e = "
        f"{candidate[persistent_n_index]:.16g}"
    )
    print(
        f"  RMSE_au = "
        f"{esp_metrics['RMSE_au']:.16g}"
    )
    print(
        f"  pearson_r = "
        f"{esp_metrics['pearson_r']:.16g}"
    )
    print(
        f"  same_sign_fraction = "
        f"{esp_metrics['same_sign_fraction']:.16g}"
    )
    print(
        f"  delta_RMS_e = "
        f"{q_metrics['delta_RMS_e']:.16g}"
    )
    print(
        f"  delta_max_abs_e = "
        f"{q_metrics['delta_max_abs_e']:.16g}"
    )
    print(
        f"  maximum_absolute_charge_e = "
        f"{q_metrics['maximum_absolute_charge_e']:.16g}"
    )
    print(
        f"  sign_change_count = "
        f"{q_metrics['sign_change_count']}"
    )
    print(
        f"  sign_change_atom_ids = "
        f"{q_metrics['sign_change_atom_ids']}"
    )
    print(
        f"  charge_sum_e = "
        f"{q_metrics['charge_sum_e']:.16g}"
    )


print("\n[6] RELATIVE PENALTIES VERSUS UNCONSTRAINED LAMBDA=4")

reference_record = next(
    record
    for record in candidate_records
    if record["policy_name"]
    == "UNCONSTRAINED_LAMBDA4"
)

reference_rmse = (
    reference_record[
        "electrostatic"
    ]["RMSE_au"]
)

reference_delta_rms = (
    reference_record[
        "charges"
    ]["delta_RMS_e"]
)

for record in candidate_records:
    rmse_ratio = (
        record[
            "electrostatic"
        ]["RMSE_au"]
        / reference_rmse
    )

    delta_rms_ratio = (
        record[
            "charges"
        ]["delta_RMS_e"]
        / reference_delta_rms
    )

    record[
        "relative_to_unconstrained"
    ] = {
        "RMSE_ratio": rmse_ratio,
        "RMSE_increase_fraction": (
            rmse_ratio - 1.0
        ),
        "delta_RMS_ratio": (
            delta_rms_ratio
        ),
        "delta_RMS_change_fraction": (
            delta_rms_ratio - 1.0
        ),
    }

    print(
        f"policy={record['policy_name']} "
        f"RMSE_ratio={rmse_ratio:.16g} "
        f"RMSE_increase_fraction="
        f"{rmse_ratio-1.0:.16g} "
        f"delta_RMS_ratio="
        f"{delta_rms_ratio:.16g} "
        f"sign_changes="
        f"{record['charges']['sign_change_count']}"
    )


print("\n[7] TARGET AND NEIGHBORING ATOM RESPONSE")

target_response_records = []

for atom_id in (
    TARGET_ATOM_ID,
    PERSISTENT_N_ATOM_ID,
):
    index = atom_index_by_id[
        atom_id
    ]

    row = real_rows[index]

    print(
        f"\natom_id={atom_id} "
        f"element={row['element']} "
        f"role={row['atom_role']} "
        f"RESP1={q0[index]: .9f}"
    )

    response = {
        "atom_id": atom_id,
        "element": row["element"],
        "atom_role": row["atom_role"],
        "RESP_stage1_charge_e": float(
            q0[index]
        ),
    }

    for policy_name, candidate in (
        candidate_vectors.items()
    ):
        charge = float(
            candidate[index]
        )

        response[
            f"{policy_name}_charge_e"
        ] = charge

        response[
            f"{policy_name}_delta_e"
        ] = (
            charge
            - q0[index]
        )

        print(
            f"  {policy_name}: "
            f"charge={charge: .9f} "
            f"delta={charge-q0[index]: .9f}"
        )

    target_response_records.append(
        response
    )


print("\n[8] WRITE CANDIDATE CHARGE TABLE")

output_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA4_TARGETED_B_CONSTRAINT.csv"
)

output_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA4_TARGETED_B_CONSTRAINT.json"
)

fieldnames = [
    "real_atom_sequence_index",
    "original_atom_index_0based",
    "atom_id",
    "element",
    "atom_role",
    "RESP_stage1_charge_e",
]

for policy_name, _ in TARGET_POLICIES:
    fieldnames.extend(
        (
            f"{policy_name}_charge_e",
            f"{policy_name}_delta_e",
        )
    )

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

    for index, row in enumerate(
        real_rows
    ):
        output = {
            "real_atom_sequence_index": index,
            "original_atom_index_0based": int(
                row["atom_index_0based"]
            ),
            "atom_id": row["atom_id"],
            "element": row["element"],
            "atom_role": row["atom_role"],
            "RESP_stage1_charge_e": (
                q0[index]
            ),
        }

        for policy_name, candidate in (
            candidate_vectors.items()
        ):
            output[
                f"{policy_name}_charge_e"
            ] = candidate[index]

            output[
                f"{policy_name}_delta_e"
            ] = (
                candidate[index]
                - q0[index]
            )

        writer.writerow(
            output
        )


print("\n[9] SCIENTIFIC GATES")

neutrality_gate = all(
    abs(
        record[
            "charges"
        ]["charge_sum_e"]
    )
    <= CHARGE_TOLERANCE_E
    for record in candidate_records
)

finite_gate = all(
    np.all(
        np.isfinite(
            candidate
        )
    )
    for candidate in (
        candidate_vectors.values()
    )
)

target_constraints_gate = all(
    (
        record[
            "target_fixed_charge_e"
        ]
        is None
    )
    or math.isclose(
        record[
            "target_candidate_charge_e"
        ],
        record[
            "target_fixed_charge_e"
        ],
        rel_tol=0.0,
        abs_tol=1.0e-10,
    )
    for record in candidate_records
)

targeted_nonnegative_gate = all(
    (
        record["policy_name"]
        == "UNCONSTRAINED_LAMBDA4"
    )
    or (
        record[
            "target_candidate_charge_e"
        ]
        >= -1.0e-10
    )
    for record in candidate_records
)

gates = {
    "upstream_decision_gate": True,
    "real_atom_count_gate": (
        len(real_rows) == 37
    ),
    "candidate_policy_count_gate": (
        len(candidate_records) == 4
    ),
    "finite_candidate_gate": (
        finite_gate
    ),
    "neutrality_gate": (
        neutrality_gate
    ),
    "target_constraint_reproduction_gate": (
        target_constraints_gate
    ),
    "targeted_nonnegative_B_gate": (
        targeted_nonnegative_gate
    ),
    "output_csv_created_gate": (
        output_csv.is_file()
        and output_csv.stat().st_size > 0
    ),
    "no_constraint_policy_adopted_gate": True,
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


print("\n[10] WRITE JSON REPORT")

decision = (
    "D039_A11_LAMBDA4_TARGETED_B_CONSTRAINT_AUDIT_PASS_"
    "CONSTRAINT_POLICY_REVIEW_AUTHORIZED"
    if all_gates_pass
    else
    "D039_A11_LAMBDA4_TARGETED_B_CONSTRAINT_AUDIT_"
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
        "A10_json": str(
            a10_json.resolve()
        ),
        "A10_json_sha256": sha256(
            a10_json
        ),
        "transferability_csv": str(
            transferability_csv.resolve()
        ),
        "transferability_csv_sha256": sha256(
            transferability_csv
        ),
        "VPOT": str(
            VPOT.resolve()
        ),
        "VPOT_sha256": sha256(
            VPOT
        ),
    },
    "regularization_lambda": (
        REGULARIZATION_LAMBDA
    ),
    "target_atom": {
        "atom_id": TARGET_ATOM_ID,
        "real_atom_sequence_index": (
            target_index
        ),
        "RESP_stage1_charge_e": float(
            target_RESP1_charge
        ),
    },
    "persistent_N_atom": {
        "atom_id": PERSISTENT_N_ATOM_ID,
        "real_atom_sequence_index": (
            persistent_n_index
        ),
        "RESP_stage1_charge_e": float(
            q0[persistent_n_index]
        ),
    },
    "candidate_records": (
        candidate_records
    ),
    "target_atom_response": (
        target_response_records
    ),
    "gates": gates,
    "authorizations": {
        "targeted_constraint_policy_scientific_review_authorized": (
            all_gates_pass
        ),
        "constraint_policy_adoption_authorized": False,
        "regularization_lambda_adoption_authorized": False,
        "constrained_refit_charge_adoption_authorized": False,
        "RESP_stage2_protocol_design_authorized": False,
        "RESP_stage2_execution_authorized": False,
        "charge_adoption_authorized": False,
        "force_field_adoption_authorized": False,
    },
    "outputs": {
        "candidate_charge_csv": str(
            output_csv.resolve()
        ),
        "candidate_charge_csv_sha256": sha256(
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


print("\n[11] DECISION")

print(f"decision={decision}")

print(
    "targeted_constraint_policy_scientific_review_authorized="
    f"{all_gates_pass}"
)

print(
    "constraint_policy_adoption_authorized=False"
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
