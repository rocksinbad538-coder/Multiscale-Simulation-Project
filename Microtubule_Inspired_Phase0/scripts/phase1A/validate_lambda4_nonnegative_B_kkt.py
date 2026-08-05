#!/usr/bin/env python3
"""
DAY039 / D039-A12

KKT validation of the local inequality-constrained lambda=4 refit:

    minimize ||Aq - V_QM||^2 + lambda ||q - q_RESP1||^2

subject to:

    sum(q) = 0
    q[A:UPPER:8:4] >= 0

The previous targeted audit showed that the equality-active candidate
q[A:UPPER:8:4] = 0 has the best performance among the tested
nonnegative target values.

This script verifies whether that candidate satisfies the
Karush-Kuhn-Tucker conditions for the inequality-constrained convex
problem.

No charge set is adopted.
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

AUTHORIZED_VPOT_SHA256 = (
    "73df47796c29d0b2a88a03d83efcb723"
    "d4f6e583d2e1789e1ea1a43d82c1064d"
)

EXPECTED_A11_DECISION = (
    "D039_A11_LAMBDA4_TARGETED_B_CONSTRAINT_AUDIT_PASS_"
    "CONSTRAINT_POLICY_REVIEW_AUTHORIZED"
)

REGULARIZATION_LAMBDA = 4.0
TARGET_ATOM_ID = "A:UPPER:8:4"

PRIMAL_TOLERANCE = 1.0e-10
DUAL_TOLERANCE = 1.0e-8
STATIONARITY_TOLERANCE = 1.0e-7
COMPLEMENTARITY_TOLERANCE = 1.0e-10


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
                "ESP point coincides with an atomic center"
            )

        matrix[start:stop, :] = 1.0 / distances

    return matrix


def solve_equality_constrained(
    hessian: np.ndarray,
    linear_term: np.ndarray,
    constraint_matrix: np.ndarray,
    constraint_values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Solve the quadratic program:

        min 1/2 q^T H q - f^T q
        subject to Cq = d

    through the KKT system.
    """

    atom_count = hessian.shape[0]
    constraint_count = constraint_matrix.shape[0]

    kkt_matrix = np.block(
        [
            [
                hessian,
                constraint_matrix.T,
            ],
            [
                constraint_matrix,
                np.zeros(
                    (
                        constraint_count,
                        constraint_count,
                    ),
                    dtype=float,
                ),
            ],
        ]
    )

    kkt_rhs = np.concatenate(
        (
            linear_term,
            constraint_values,
        )
    )

    solution = np.linalg.solve(
        kkt_matrix,
        kkt_rhs,
    )

    charges = solution[
        :atom_count
    ]

    equality_multipliers = solution[
        atom_count:
    ]

    return charges, equality_multipliers


def objective_value(
    matrix: np.ndarray,
    target: np.ndarray,
    charges: np.ndarray,
    reference_charges: np.ndarray,
    regularization_lambda: float,
) -> float:
    residual = (
        matrix @ charges
        - target
    )

    charge_delta = (
        charges
        - reference_charges
    )

    return float(
        np.dot(
            residual,
            residual,
        )
        + regularization_lambda
        * np.dot(
            charge_delta,
            charge_delta,
        )
    )


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
        "pearson_r": pearson_r,
        "same_sign_fraction": float(
            np.mean(
                np.sign(reference)
                == np.sign(candidate)
            )
        ),
    }


def json_safe_value(value):
    """
    Convert NumPy scalar and array objects into JSON-safe
    built-in Python objects.
    """

    if isinstance(value, np.bool_):
        return bool(value)

    if isinstance(value, np.integer):
        return int(value)

    if isinstance(value, np.floating):
        return float(value)

    if isinstance(value, np.ndarray):
        return value.tolist()

    raise TypeError(
        f"Object of type {type(value).__name__} "
        "is not JSON serializable"
    )


print("=" * 100)
print("DAY039 / D039-A12 — NONNEGATIVE TARGET-B KKT VALIDATION")
print("=" * 100)


print("\n[1] SOURCE IDENTITY")

for path in (
    LATEST_POINTER,
    VPOT,
):
    require_file(path)

observed_vpot_sha256 = sha256(
    VPOT
)

if (
    observed_vpot_sha256
    != AUTHORIZED_VPOT_SHA256
):
    raise RuntimeError(
        "Authorized VPOT SHA256 mismatch"
    )

execution_dir = (
    ROOT
    / LATEST_POINTER.read_text(
        encoding="utf-8"
    ).strip()
)

a11_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA4_TARGETED_B_CONSTRAINT.json"
)

a11_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA4_TARGETED_B_CONSTRAINT.csv"
)

transferability_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1_TRANSFERABILITY.csv"
)

for path in (
    a11_json,
    a11_csv,
    transferability_csv,
):
    require_file(path)
    print(f"FOUND  {path}")

print(f"VPOT_SHA256 = {observed_vpot_sha256}")


print("\n[2] UPSTREAM AUTHORIZATION")

a11_report = load_json(
    a11_json
)

if (
    a11_report.get("decision")
    != EXPECTED_A11_DECISION
):
    raise RuntimeError(
        "Unexpected A11 decision.\n"
        f"Observed: {a11_report.get('decision')}"
    )

authorizations = a11_report.get(
    "authorizations",
    {},
)

if (
    authorizations.get(
        "targeted_constraint_policy_scientific_review_authorized"
    )
    is not True
):
    raise RuntimeError(
        "KKT validation is not authorized"
    )

if (
    authorizations.get(
        "constraint_policy_adoption_authorized"
    )
    is not False
):
    raise RuntimeError(
        "Unexpected constraint-policy authorization"
    )

print("A11_decision_gate                  = PASS")
print("KKT_validation_review_gate         = PASS")
print("constraint_policy_adoption_blocked = PASS")
print("charge_adoption_blocked_gate       = PASS")


print("\n[3] LOAD VPOT AND REAL-ATOM SYSTEM")

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

target_index = atom_index_by_id[
    TARGET_ATOM_ID
]

print(f"real_atom_count = {len(real_rows)}")
print(f"grid_point_count = {len(grid_xyz_bohr)}")
print(f"target_atom_id = {TARGET_ATOM_ID}")
print(f"target_real_index = {target_index}")
print(
    f"target_RESP_stage1_charge_e = "
    f"{q0[target_index]:.16g}"
)


print("\n[4] BUILD QUADRATIC PROBLEM")

design_matrix = build_coulomb_matrix(
    grid_xyz_bohr,
    real_xyz_bohr,
)

hessian = 2.0 * (
    design_matrix.T @ design_matrix
    + REGULARIZATION_LAMBDA
    * np.eye(
        len(real_rows),
        dtype=float,
    )
)

linear_term = 2.0 * (
    design_matrix.T @ quantum_esp
    + REGULARIZATION_LAMBDA
    * q0
)

hessian_eigenvalues = np.linalg.eigvalsh(
    hessian
)

minimum_hessian_eigenvalue = float(
    np.min(
        hessian_eigenvalues
    )
)

maximum_hessian_eigenvalue = float(
    np.max(
        hessian_eigenvalues
    )
)

print(
    f"design_matrix_shape = "
    f"{design_matrix.shape}"
)
print(
    f"minimum_hessian_eigenvalue = "
    f"{minimum_hessian_eigenvalue:.16g}"
)
print(
    f"maximum_hessian_eigenvalue = "
    f"{maximum_hessian_eigenvalue:.16g}"
)
print(
    f"hessian_condition_number = "
    f"{maximum_hessian_eigenvalue/minimum_hessian_eigenvalue:.16g}"
)


print("\n[5] SOLVE UNCONSTRAINED-IN-SIGN LAMBDA=4")

neutrality_constraint = np.ones(
    (
        1,
        len(real_rows),
    ),
    dtype=float,
)

neutrality_value = np.asarray(
    [0.0],
    dtype=float,
)

unconstrained_charges, unconstrained_multipliers = (
    solve_equality_constrained(
        hessian,
        linear_term,
        neutrality_constraint,
        neutrality_value,
    )
)

unconstrained_target_charge = float(
    unconstrained_charges[
        target_index
    ]
)

unconstrained_objective = objective_value(
    design_matrix,
    quantum_esp,
    unconstrained_charges,
    q0,
    REGULARIZATION_LAMBDA,
)

print(
    f"target_charge_e = "
    f"{unconstrained_target_charge:.16g}"
)
print(
    f"neutrality_residual_e = "
    f"{np.sum(unconstrained_charges):.16g}"
)
print(
    f"objective_value = "
    f"{unconstrained_objective:.16g}"
)
print(
    f"neutrality_multiplier = "
    f"{unconstrained_multipliers[0]:.16g}"
)


print("\n[6] SOLVE ACTIVE-BOUNDARY CANDIDATE q_B=0")

active_constraint_matrix = np.zeros(
    (
        2,
        len(real_rows),
    ),
    dtype=float,
)

active_constraint_matrix[
    0,
    :,
] = 1.0

active_constraint_matrix[
    1,
    target_index,
] = 1.0

active_constraint_values = np.asarray(
    [
        0.0,
        0.0,
    ],
    dtype=float,
)

active_charges, active_multipliers = (
    solve_equality_constrained(
        hessian,
        linear_term,
        active_constraint_matrix,
        active_constraint_values,
    )
)

active_target_charge = float(
    active_charges[
        target_index
    ]
)

active_objective = objective_value(
    design_matrix,
    quantum_esp,
    active_charges,
    q0,
    REGULARIZATION_LAMBDA,
)

print(
    f"target_charge_e = "
    f"{active_target_charge:.16g}"
)
print(
    f"neutrality_residual_e = "
    f"{np.sum(active_charges):.16g}"
)
print(
    f"objective_value = "
    f"{active_objective:.16g}"
)
print(
    f"objective_increase = "
    f"{active_objective-unconstrained_objective:.16g}"
)
print(
    f"objective_increase_fraction = "
    f"{active_objective/unconstrained_objective-1.0:.16g}"
)
print(
    f"neutrality_multiplier = "
    f"{active_multipliers[0]:.16g}"
)
print(
    f"active_equality_multiplier_raw = "
    f"{active_multipliers[1]:.16g}"
)


print("\n[7] KKT CONDITIONS FOR q_TARGET >= 0")

# Write the inequality in standard minimization form:
#
#     g(q) = -q_target <= 0
#
# with multiplier mu >= 0.
#
# Stationarity:
#
#     gradient(f) + nu * 1 - mu * e_target = 0
#
# The equality-constrained KKT solve used:
#
#     gradient(f) + nu * 1 + eta * e_target = 0
#
# Therefore:
#
#     mu = -eta.

gradient = (
    hessian @ active_charges
    - linear_term
)

neutrality_multiplier = float(
    active_multipliers[0]
)

active_equality_multiplier = float(
    active_multipliers[1]
)

inequality_multiplier = (
    -active_equality_multiplier
)

target_basis = np.zeros(
    len(real_rows),
    dtype=float,
)

target_basis[
    target_index
] = 1.0

stationarity_residual = (
    gradient
    + neutrality_multiplier
    * np.ones(
        len(real_rows),
        dtype=float,
    )
    - inequality_multiplier
    * target_basis
)

primal_neutrality_residual = float(
    np.sum(active_charges)
)

primal_inequality_value = float(
    -active_charges[
        target_index
    ]
)

dual_feasibility_value = float(
    inequality_multiplier
)

complementarity_value = float(
    inequality_multiplier
    * active_charges[
        target_index
    ]
)

maximum_stationarity_residual = float(
    np.max(
        np.abs(
            stationarity_residual
        )
    )
)

primal_neutrality_gate = (
    abs(
        primal_neutrality_residual
    )
    <= PRIMAL_TOLERANCE
)

primal_inequality_gate = (
    active_charges[
        target_index
    ]
    >= -PRIMAL_TOLERANCE
)

dual_feasibility_gate = (
    inequality_multiplier
    >= -DUAL_TOLERANCE
)

stationarity_gate = (
    maximum_stationarity_residual
    <= STATIONARITY_TOLERANCE
)

complementarity_gate = (
    abs(
        complementarity_value
    )
    <= COMPLEMENTARITY_TOLERANCE
)

print(
    f"primal_neutrality_residual_e = "
    f"{primal_neutrality_residual:.16g}"
)
print(
    f"primal_inequality_g_value = "
    f"{primal_inequality_value:.16g}"
)
print(
    f"inequality_multiplier_mu = "
    f"{inequality_multiplier:.16g}"
)
print(
    f"maximum_stationarity_residual = "
    f"{maximum_stationarity_residual:.16g}"
)
print(
    f"complementarity_mu_times_q = "
    f"{complementarity_value:.16g}"
)

print(
    f"primal_neutrality_gate = "
    f"{'PASS' if primal_neutrality_gate else 'FAIL'}"
)
print(
    f"primal_inequality_gate = "
    f"{'PASS' if primal_inequality_gate else 'FAIL'}"
)
print(
    f"dual_feasibility_gate = "
    f"{'PASS' if dual_feasibility_gate else 'FAIL'}"
)
print(
    f"stationarity_gate = "
    f"{'PASS' if stationarity_gate else 'FAIL'}"
)
print(
    f"complementarity_gate = "
    f"{'PASS' if complementarity_gate else 'FAIL'}"
)


print("\n[8] ELECTROSTATIC AND CHARGE SUMMARY")

unconstrained_potential = (
    design_matrix
    @ unconstrained_charges
)

active_potential = (
    design_matrix
    @ active_charges
)

unconstrained_metrics = (
    electrostatic_metrics(
        quantum_esp,
        unconstrained_potential,
    )
)

active_metrics = electrostatic_metrics(
    quantum_esp,
    active_potential,
)

unconstrained_delta = (
    unconstrained_charges
    - q0
)

active_delta = (
    active_charges
    - q0
)

unconstrained_sign_changes = [
    real_rows[index]["atom_id"]
    for index in range(
        len(real_rows)
    )
    if (
        q0[index] != 0.0
        and unconstrained_charges[
            index
        ] != 0.0
        and math.copysign(
            1.0,
            q0[index],
        )
        != math.copysign(
            1.0,
            unconstrained_charges[
                index
            ],
        )
    )
]

active_sign_changes = [
    real_rows[index]["atom_id"]
    for index in range(
        len(real_rows)
    )
    if (
        q0[index] != 0.0
        and active_charges[
            index
        ] != 0.0
        and math.copysign(
            1.0,
            q0[index],
        )
        != math.copysign(
            1.0,
            active_charges[
                index
            ],
        )
    )
]

print("\nUNCONSTRAINED_LAMBDA4")

for name, value in (
    unconstrained_metrics.items()
):
    print(f"  {name} = {value}")

print(
    f"  delta_RMS_e = "
    f"{np.sqrt(np.mean(unconstrained_delta**2)):.16g}"
)
print(
    f"  delta_max_abs_e = "
    f"{np.max(np.abs(unconstrained_delta)):.16g}"
)
print(
    f"  sign_change_atom_ids = "
    f"{unconstrained_sign_changes}"
)

print("\nNONNEGATIVE_B_ACTIVE_BOUNDARY")

for name, value in (
    active_metrics.items()
):
    print(f"  {name} = {value}")

print(
    f"  delta_RMS_e = "
    f"{np.sqrt(np.mean(active_delta**2)):.16g}"
)
print(
    f"  delta_max_abs_e = "
    f"{np.max(np.abs(active_delta)):.16g}"
)
print(
    f"  sign_change_atom_ids = "
    f"{active_sign_changes}"
)


print("\n[9] CROSS-CHECK AGAINST A11 q_B=0 CANDIDATE")

with a11_csv.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    a11_rows = list(
        csv.DictReader(handle)
    )

a11_rows.sort(
    key=lambda row: int(
        row["real_atom_sequence_index"]
    )
)

a11_zero_charges = np.asarray(
    [
        float(
            row[
                "TARGET_B_FIXED_ZERO_charge_e"
            ]
        )
        for row in a11_rows
    ],
    dtype=float,
)

maximum_A11_difference_e = float(
    np.max(
        np.abs(
            active_charges
            - a11_zero_charges
        )
    )
)

A11_reproduction_gate = (
    maximum_A11_difference_e
    <= 1.0e-10
)

print(
    f"maximum_A11_zero_candidate_difference_e = "
    f"{maximum_A11_difference_e:.16g}"
)
print(
    f"A11_zero_candidate_reproduction_gate = "
    f"{'PASS' if A11_reproduction_gate else 'FAIL'}"
)


print("\n[10] WRITE CANDIDATE CHARGE TABLE")

output_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA4_NONNEGATIVE_B_KKT.csv"
)

output_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA4_NONNEGATIVE_B_KKT.json"
)

fieldnames = [
    "real_atom_sequence_index",
    "original_atom_index_0based",
    "atom_id",
    "element",
    "atom_role",
    "RESP_stage1_charge_e",
    "unconstrained_lambda4_charge_e",
    "unconstrained_lambda4_delta_e",
    "nonnegative_B_active_charge_e",
    "nonnegative_B_active_delta_e",
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

    for index, row in enumerate(
        real_rows
    ):
        writer.writerow(
            {
                "real_atom_sequence_index": (
                    index
                ),
                "original_atom_index_0based": int(
                    row[
                        "atom_index_0based"
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
                    q0[index]
                ),
                "unconstrained_lambda4_charge_e": (
                    unconstrained_charges[
                        index
                    ]
                ),
                "unconstrained_lambda4_delta_e": (
                    unconstrained_delta[
                        index
                    ]
                ),
                "nonnegative_B_active_charge_e": (
                    active_charges[
                        index
                    ]
                ),
                "nonnegative_B_active_delta_e": (
                    active_delta[
                        index
                    ]
                ),
            }
        )


print("\n[11] FINAL GATES")

strict_convexity_gate = (
    minimum_hessian_eigenvalue
    > 0.0
)

unconstrained_violates_inequality_gate = (
    unconstrained_target_charge
    < 0.0
)

active_boundary_gate = (
    abs(
        active_target_charge
    )
    <= PRIMAL_TOLERANCE
)

all_KKT_gates_pass = all(
    (
        primal_neutrality_gate,
        primal_inequality_gate,
        dual_feasibility_gate,
        stationarity_gate,
        complementarity_gate,
    )
)

gates = {
    "upstream_decision_gate": True,
    "strict_convexity_gate": (
        strict_convexity_gate
    ),
    "unconstrained_solution_violates_nonnegative_B_gate": (
        unconstrained_violates_inequality_gate
    ),
    "active_boundary_qB_zero_gate": (
        active_boundary_gate
    ),
    "primal_neutrality_gate": (
        primal_neutrality_gate
    ),
    "primal_inequality_gate": (
        primal_inequality_gate
    ),
    "dual_feasibility_gate": (
        dual_feasibility_gate
    ),
    "stationarity_gate": (
        stationarity_gate
    ),
    "complementarity_gate": (
        complementarity_gate
    ),
    "all_KKT_conditions_gate": (
        all_KKT_gates_pass
    ),
    "A11_candidate_reproduction_gate": (
        A11_reproduction_gate
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


print("\n[12] WRITE JSON REPORT")

decision = (
    "D039_A12_NONNEGATIVE_B_KKT_VALIDATION_PASS_"
    "METHOD_CANDIDATE_REVIEW_AUTHORIZED"
    if all_gates_pass
    else
    "D039_A12_NONNEGATIVE_B_KKT_VALIDATION_"
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
    "mathematical_problem": {
        "objective": (
            "MINIMIZE_||Aq-V_QM||2_PLUS_"
            "LAMBDA_||q-q_RESP1_REAL||2"
        ),
        "regularization_lambda": (
            REGULARIZATION_LAMBDA
        ),
        "equality_constraint": (
            "SUM_q_EQUALS_ZERO"
        ),
        "inequality_constraint": (
            f"q[{TARGET_ATOM_ID}]_GREATER_THAN_OR_EQUAL_ZERO"
        ),
        "inequality_standard_form": (
            f"-q[{TARGET_ATOM_ID}]_LESS_THAN_OR_EQUAL_ZERO"
        ),
    },
    "source_identity": {
        "A11_json": str(
            a11_json.resolve()
        ),
        "A11_json_sha256": sha256(
            a11_json
        ),
        "A11_csv": str(
            a11_csv.resolve()
        ),
        "A11_csv_sha256": sha256(
            a11_csv
        ),
        "VPOT": str(
            VPOT.resolve()
        ),
        "VPOT_sha256": sha256(
            VPOT
        ),
        "transferability_csv": str(
            transferability_csv.resolve()
        ),
        "transferability_csv_sha256": sha256(
            transferability_csv
        ),
    },
    "convexity": {
        "minimum_hessian_eigenvalue": (
            minimum_hessian_eigenvalue
        ),
        "maximum_hessian_eigenvalue": (
            maximum_hessian_eigenvalue
        ),
        "condition_number": (
            maximum_hessian_eigenvalue
            / minimum_hessian_eigenvalue
        ),
    },
    "unconstrained_solution": {
        "target_charge_e": (
            unconstrained_target_charge
        ),
        "objective_value": (
            unconstrained_objective
        ),
        "electrostatic_metrics": (
            unconstrained_metrics
        ),
        "sign_change_atom_ids": (
            unconstrained_sign_changes
        ),
    },
    "active_boundary_solution": {
        "target_charge_e": (
            active_target_charge
        ),
        "objective_value": (
            active_objective
        ),
        "objective_increase": (
            active_objective
            - unconstrained_objective
        ),
        "objective_increase_fraction": (
            active_objective
            / unconstrained_objective
            - 1.0
        ),
        "electrostatic_metrics": (
            active_metrics
        ),
        "sign_change_atom_ids": (
            active_sign_changes
        ),
    },
    "KKT": {
        "primal_neutrality_residual_e": (
            primal_neutrality_residual
        ),
        "primal_inequality_g_value": (
            primal_inequality_value
        ),
        "inequality_multiplier_mu": (
            inequality_multiplier
        ),
        "maximum_stationarity_residual": (
            maximum_stationarity_residual
        ),
        "complementarity_mu_times_q": (
            complementarity_value
        ),
    },
    "gates": gates,
    "authorizations": {
        "lambda4_nonnegative_B_method_candidate_review_authorized": (
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
        default=json_safe_value,
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


print("\n[13] DECISION")

print(f"decision={decision}")

print(
    "lambda4_nonnegative_B_method_candidate_review_authorized="
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
