#!/usr/bin/env python3
"""
DAY039 / D039-A14

Genuine deterministic train-only hold-out validation for the
37-real-atom electrostatic refit.

Unlike D039-A13, this block refits the candidate using only the
deterministic 80% training subset. The 20% validation subset is not
used in the optimization.

Optimization
------------
minimize:

    ||A_train q - V_train||^2
    + lambda ||q - q_RESP1_real||^2

subject to:

    sum(q) = 0
    q[A:UPPER:8:4] >= 0

The inequality is solved through an active-set decision:

1. Solve with exact neutrality only.
2. If q_target >= 0, the inequality is inactive.
3. If q_target < 0, solve with q_target = 0 as an active constraint.
4. Verify the corresponding KKT conditions.

Models evaluated
----------------
- RESP52 fixed baseline
- REAL37_UNMODIFIED fixed baseline
- TRAIN_ONLY_LAMBDA4_UNCONSTRAINED
- TRAIN_ONLY_LAMBDA4_NONNEGATIVE_B
- FULL_GRID_LAMBDA4_NONNEGATIVE_B from D039-A12

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

EXPECTED_A13_DECISION = (
    "D039_A13_SPATIAL_GENERALIZATION_VALIDATION_PASS_"
    "METHOD_ADOPTION_REVIEW_AUTHORIZED"
)

EXPECTED_ATOM_COUNT = 52
EXPECTED_REAL_COUNT = 37
EXPECTED_GRID_COUNT = 24835

REGULARIZATION_LAMBDA = 4.0
TARGET_ATOM_ID = "A:UPPER:8:4"

VALIDATION_MODULUS = 5
VALIDATION_REMAINDER = 0

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


def json_safe_value(value):
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

        matrix[start:stop, :] = (
            1.0 / distances
        )

    return matrix


def solve_equality_qp(
    hessian: np.ndarray,
    linear_term: np.ndarray,
    constraint_matrix: np.ndarray,
    constraint_values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    atom_count = hessian.shape[0]
    constraint_count = (
        constraint_matrix.shape[0]
    )

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

    return (
        solution[:atom_count],
        solution[atom_count:],
    )


def calculate_metrics(
    reference: np.ndarray,
    candidate: np.ndarray,
) -> dict:
    residual = candidate - reference
    absolute_residual = np.abs(residual)

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
        "point_count": int(
            len(reference)
        ),
        "RMSE_au": float(
            np.sqrt(
                np.mean(residual ** 2)
            )
        ),
        "MAE_au": float(
            np.mean(absolute_residual)
        ),
        "maximum_absolute_error_au": float(
            np.max(absolute_residual)
        ),
        "residual_mean_au": float(
            np.mean(residual)
        ),
        "pearson_r": pearson_r,
        "r_squared": pearson_r ** 2,
        "same_sign_fraction": float(
            np.mean(
                np.sign(reference)
                == np.sign(candidate)
            )
        ),
    }


def charge_comparison(
    reference: np.ndarray,
    candidate: np.ndarray,
    atom_rows: list[dict],
) -> dict:
    delta = candidate - reference

    sign_change_indices = [
        index
        for index in range(len(reference))
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
            np.max(np.abs(candidate))
        ),
        "delta_MAE_e": float(
            np.mean(np.abs(delta))
        ),
        "delta_RMS_e": float(
            np.sqrt(
                np.mean(delta ** 2)
            )
        ),
        "delta_max_abs_e": float(
            np.max(np.abs(delta))
        ),
        "sign_change_count": len(
            sign_change_indices
        ),
        "sign_change_atom_ids": [
            atom_rows[index]["atom_id"]
            for index in sign_change_indices
        ],
    }


def evaluate_model(
    model_name: str,
    full_matrix: np.ndarray,
    charges: np.ndarray,
    quantum_esp: np.ndarray,
    training_indices: np.ndarray,
    validation_indices: np.ndarray,
) -> dict:
    potential = full_matrix @ charges

    training_metrics = calculate_metrics(
        quantum_esp[training_indices],
        potential[training_indices],
    )

    validation_metrics = calculate_metrics(
        quantum_esp[validation_indices],
        potential[validation_indices],
    )

    validation_to_training_ratio = (
        validation_metrics["RMSE_au"]
        / training_metrics["RMSE_au"]
    )

    return {
        "model_name": model_name,
        "training": training_metrics,
        "validation": validation_metrics,
        "generalization_gap_au": (
            validation_metrics["RMSE_au"]
            - training_metrics["RMSE_au"]
        ),
        "validation_to_training_RMSE_ratio": (
            validation_to_training_ratio
        ),
        "relative_generalization_gap": (
            validation_to_training_ratio
            - 1.0
        ),
    }


print("=" * 100)
print("DAY039 / D039-A14 — GENUINE TRAIN-ONLY CONSTRAINED REFIT VALIDATION")
print("=" * 100)


print("\n[1] SOURCE IDENTITY")

for path in (
    LATEST_POINTER,
    VPOT,
):
    require_file(path)

if sha256(VPOT) != AUTHORIZED_VPOT_SHA256:
    raise RuntimeError(
        "Authorized VPOT SHA256 mismatch"
    )

execution_dir = (
    ROOT
    / LATEST_POINTER.read_text(
        encoding="utf-8"
    ).strip()
)

a13_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_SPATIAL_HOLDOUT_GENERALIZATION.json"
)

a12_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA4_NONNEGATIVE_B_KKT.csv"
)

transferability_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1_TRANSFERABILITY.csv"
)

for path in (
    a13_json,
    a12_csv,
    transferability_csv,
):
    require_file(path)
    print(f"FOUND  {path}")

print("source_identity_gate = PASS")


print("\n[2] UPSTREAM AUTHORIZATION")

a13_report = load_json(
    a13_json
)

if (
    a13_report.get("decision")
    != EXPECTED_A13_DECISION
):
    raise RuntimeError(
        "Unexpected A13 decision.\n"
        f"Observed: {a13_report.get('decision')}"
    )

authorizations = a13_report.get(
    "authorizations",
    {},
)

if (
    authorizations.get(
        "lambda4_nonnegative_B_method_adoption_review_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Train-only refit validation is not authorized"
    )

if (
    authorizations.get(
        "charge_adoption_authorized"
    )
    is not False
):
    raise RuntimeError(
        "Unexpected charge-adoption authorization"
    )

print("A13_decision_gate                    = PASS")
print("train_only_refit_review_gate         = PASS")
print("charge_adoption_blocked_gate         = PASS")
print("RESP_stage2_execution_blocked_gate   = PASS")


print("\n[3] LOAD VPOT")

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

if atom_xyz_bohr.shape != (
    EXPECTED_ATOM_COUNT,
    3,
):
    raise RuntimeError(
        f"Unexpected atom shape: {atom_xyz_bohr.shape}"
    )

if grid_xyz_bohr.shape != (
    EXPECTED_GRID_COUNT,
    3,
):
    raise RuntimeError(
        f"Unexpected grid shape: {grid_xyz_bohr.shape}"
    )

print(f"atom_coordinate_shape = {atom_xyz_bohr.shape}")
print(f"grid_coordinate_shape = {grid_xyz_bohr.shape}")
print(f"quantum_ESP_shape = {quantum_esp.shape}")


print("\n[4] LOAD ATOMS AND CHARGES")

with transferability_csv.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    atom_rows = list(
        csv.DictReader(handle)
    )

atom_rows.sort(
    key=lambda row: int(
        row["atom_index_0based"]
    )
)

RESP52_charges = np.asarray(
    [
        float(
            row["RESP_stage1_charge_e"]
        )
        for row in atom_rows
    ],
    dtype=float,
)

real_mask = np.asarray(
    [
        not parse_bool(
            row["artificial_cap"]
        )
        for row in atom_rows
    ],
    dtype=bool,
)

real_rows = [
    row
    for row, retained in zip(
        atom_rows,
        real_mask,
    )
    if retained
]

real_xyz_bohr = atom_xyz_bohr[
    real_mask
]

q0_real = RESP52_charges[
    real_mask
]

if len(real_rows) != EXPECTED_REAL_COUNT:
    raise RuntimeError(
        f"Expected 37 real atoms, observed {len(real_rows)}"
    )

atom_index_by_id = {
    row["atom_id"]: index
    for index, row in enumerate(
        real_rows
    )
}

if TARGET_ATOM_ID not in atom_index_by_id:
    raise RuntimeError(
        f"Missing target atom: {TARGET_ATOM_ID}"
    )

target_index = atom_index_by_id[
    TARGET_ATOM_ID
]

with a12_csv.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    a12_rows = list(
        csv.DictReader(handle)
    )

a12_rows.sort(
    key=lambda row: int(
        row["real_atom_sequence_index"]
    )
)

if [
    row["atom_id"]
    for row in a12_rows
] != [
    row["atom_id"]
    for row in real_rows
]:
    raise RuntimeError(
        "A12 and retained real-atom orders differ"
    )

full_grid_nonnegative_B_charges = np.asarray(
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

print(f"real_atom_count = {len(real_rows)}")
print(f"target_atom_id = {TARGET_ATOM_ID}")
print(f"target_real_index = {target_index}")
print(
    f"target_RESP1_charge_e = "
    f"{q0_real[target_index]:.16g}"
)


print("\n[5] RECONSTRUCT DETERMINISTIC PARTITION")

lexicographic_order = np.lexsort(
    (
        grid_xyz_bohr[:, 2],
        grid_xyz_bohr[:, 1],
        grid_xyz_bohr[:, 0],
    )
)

sorted_positions = np.arange(
    EXPECTED_GRID_COUNT,
    dtype=int,
)

validation_mask = (
    sorted_positions
    % VALIDATION_MODULUS
    == VALIDATION_REMAINDER
)

training_indices = lexicographic_order[
    ~validation_mask
]

validation_indices = lexicographic_order[
    validation_mask
]

training_points = len(
    training_indices
)

validation_points = len(
    validation_indices
)

print(f"training_points = {training_points}")
print(f"validation_points = {validation_points}")
print(
    f"training_fraction = "
    f"{training_points/EXPECTED_GRID_COUNT:.16g}"
)
print(
    f"validation_fraction = "
    f"{validation_points/EXPECTED_GRID_COUNT:.16g}"
)


print("\n[6] BUILD TRAINING AND FULL COULOMB MATRICES")

full_real_matrix = build_coulomb_matrix(
    grid_xyz_bohr,
    real_xyz_bohr,
)

training_matrix = full_real_matrix[
    training_indices,
    :
]

training_target = quantum_esp[
    training_indices
]

full_RESP52_matrix = build_coulomb_matrix(
    grid_xyz_bohr,
    atom_xyz_bohr,
)

print(
    f"full_real_matrix_shape = "
    f"{full_real_matrix.shape}"
)
print(
    f"training_matrix_shape = "
    f"{training_matrix.shape}"
)
print(
    f"full_RESP52_matrix_shape = "
    f"{full_RESP52_matrix.shape}"
)


print("\n[7] BUILD TRAIN-ONLY QUADRATIC PROBLEM")

hessian = 2.0 * (
    training_matrix.T
    @ training_matrix
    + REGULARIZATION_LAMBDA
    * np.eye(
        EXPECTED_REAL_COUNT,
        dtype=float,
    )
)

linear_term = 2.0 * (
    training_matrix.T
    @ training_target
    + REGULARIZATION_LAMBDA
    * q0_real
)

eigenvalues = np.linalg.eigvalsh(
    hessian
)

minimum_eigenvalue = float(
    np.min(eigenvalues)
)

maximum_eigenvalue = float(
    np.max(eigenvalues)
)

print(
    f"minimum_hessian_eigenvalue = "
    f"{minimum_eigenvalue:.16g}"
)
print(
    f"maximum_hessian_eigenvalue = "
    f"{maximum_eigenvalue:.16g}"
)
print(
    f"hessian_condition_number = "
    f"{maximum_eigenvalue/minimum_eigenvalue:.16g}"
)


print("\n[8] TRAIN-ONLY UNCONSTRAINED-IN-SIGN SOLUTION")

neutrality_matrix = np.ones(
    (
        1,
        EXPECTED_REAL_COUNT,
    ),
    dtype=float,
)

neutrality_value = np.asarray(
    [0.0],
    dtype=float,
)

train_unconstrained_charges, train_unconstrained_multipliers = (
    solve_equality_qp(
        hessian,
        linear_term,
        neutrality_matrix,
        neutrality_value,
    )
)

unconstrained_target_charge = float(
    train_unconstrained_charges[
        target_index
    ]
)

print(
    f"target_charge_e = "
    f"{unconstrained_target_charge:.16g}"
)
print(
    f"charge_sum_e = "
    f"{np.sum(train_unconstrained_charges):.16g}"
)
print(
    f"inequality_violated = "
    f"{unconstrained_target_charge < 0.0}"
)


print("\n[9] ACTIVE-SET NONNEGATIVE-B SOLUTION")

if unconstrained_target_charge >= 0.0:
    inequality_status = "INACTIVE"

    train_nonnegative_B_charges = (
        train_unconstrained_charges.copy()
    )

    inequality_multiplier = 0.0

    active_multipliers = (
        train_unconstrained_multipliers
    )

else:
    inequality_status = "ACTIVE_AT_ZERO"

    active_constraint_matrix = np.zeros(
        (
            2,
            EXPECTED_REAL_COUNT,
        ),
        dtype=float,
    )

    active_constraint_matrix[0, :] = 1.0

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

    train_nonnegative_B_charges, active_multipliers = (
        solve_equality_qp(
            hessian,
            linear_term,
            active_constraint_matrix,
            active_constraint_values,
        )
    )

    inequality_multiplier = float(
        -active_multipliers[1]
    )

print(
    f"inequality_status = "
    f"{inequality_status}"
)
print(
    f"target_charge_e = "
    f"{train_nonnegative_B_charges[target_index]:.16g}"
)
print(
    f"charge_sum_e = "
    f"{np.sum(train_nonnegative_B_charges):.16g}"
)
print(
    f"inequality_multiplier_mu = "
    f"{inequality_multiplier:.16g}"
)


print("\n[10] TRAIN-ONLY KKT VALIDATION")

gradient = (
    hessian
    @ train_nonnegative_B_charges
    - linear_term
)

neutrality_multiplier = float(
    active_multipliers[0]
)

target_basis = np.zeros(
    EXPECTED_REAL_COUNT,
    dtype=float,
)

target_basis[
    target_index
] = 1.0

stationarity_residual = (
    gradient
    + neutrality_multiplier
    * np.ones(
        EXPECTED_REAL_COUNT,
        dtype=float,
    )
    - inequality_multiplier
    * target_basis
)

primal_neutrality_residual = float(
    np.sum(
        train_nonnegative_B_charges
    )
)

target_charge = float(
    train_nonnegative_B_charges[
        target_index
    ]
)

primal_inequality_value = (
    -target_charge
)

complementarity_value = (
    inequality_multiplier
    * target_charge
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
    target_charge
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


print("\n[11] OUT-OF-SAMPLE MODEL EVALUATION")

RESP52_potential = (
    full_RESP52_matrix
    @ RESP52_charges
)

REAL37_UNMODIFIED_potential = (
    full_real_matrix
    @ q0_real
)

models = {
    "RESP52": (
        RESP52_potential
    ),
    "REAL37_UNMODIFIED": (
        REAL37_UNMODIFIED_potential
    ),
    "TRAIN_ONLY_LAMBDA4_UNCONSTRAINED": (
        full_real_matrix
        @ train_unconstrained_charges
    ),
    "TRAIN_ONLY_LAMBDA4_NONNEGATIVE_B": (
        full_real_matrix
        @ train_nonnegative_B_charges
    ),
    "FULL_GRID_LAMBDA4_NONNEGATIVE_B": (
        full_real_matrix
        @ full_grid_nonnegative_B_charges
    ),
}

model_results = []

for model_name, potential in (
    models.items()
):
    training_metrics = calculate_metrics(
        quantum_esp[
            training_indices
        ],
        potential[
            training_indices
        ],
    )

    validation_metrics = calculate_metrics(
        quantum_esp[
            validation_indices
        ],
        potential[
            validation_indices
        ],
    )

    ratio = (
        validation_metrics["RMSE_au"]
        / training_metrics["RMSE_au"]
    )

    record = {
        "model_name": model_name,
        "training": training_metrics,
        "validation": validation_metrics,
        "generalization_gap_au": (
            validation_metrics["RMSE_au"]
            - training_metrics["RMSE_au"]
        ),
        "validation_to_training_RMSE_ratio": (
            ratio
        ),
        "relative_generalization_gap": (
            ratio - 1.0
        ),
    }

    model_results.append(record)

    print(f"\nmodel={model_name}")

    print(
        f"  training_RMSE_au = "
        f"{training_metrics['RMSE_au']:.16g}"
    )
    print(
        f"  validation_RMSE_au = "
        f"{validation_metrics['RMSE_au']:.16g}"
    )
    print(
        f"  generalization_gap_au = "
        f"{record['generalization_gap_au']:.16g}"
    )
    print(
        f"  validation_to_training_RMSE_ratio = "
        f"{ratio:.16g}"
    )
    print(
        f"  validation_pearson_r = "
        f"{validation_metrics['pearson_r']:.16g}"
    )
    print(
        f"  validation_same_sign_fraction = "
        f"{validation_metrics['same_sign_fraction']:.16g}"
    )


print("\n[12] VALIDATION RANKING")

validation_ranking = sorted(
    model_results,
    key=lambda record: (
        record[
            "validation"
        ]["RMSE_au"]
    ),
)

for rank, record in enumerate(
    validation_ranking,
    start=1,
):
    record[
        "validation_rank"
    ] = rank

    print(
        f"rank={rank} "
        f"model={record['model_name']} "
        f"validation_RMSE_au="
        f"{record['validation']['RMSE_au']:.16g} "
        f"training_RMSE_au="
        f"{record['training']['RMSE_au']:.16g} "
        f"ratio="
        f"{record['validation_to_training_RMSE_ratio']:.16g}"
    )


print("\n[13] CHARGE STABILITY AGAINST FULL-GRID A12")

charge_delta = (
    train_nonnegative_B_charges
    - full_grid_nonnegative_B_charges
)

charge_stability = {
    "RMS_difference_e": float(
        np.sqrt(
            np.mean(
                charge_delta ** 2
            )
        )
    ),
    "MAE_difference_e": float(
        np.mean(
            np.abs(
                charge_delta
            )
        )
    ),
    "maximum_absolute_difference_e": float(
        np.max(
            np.abs(
                charge_delta
            )
        )
    ),
    "pearson_r": float(
        np.corrcoef(
            train_nonnegative_B_charges,
            full_grid_nonnegative_B_charges,
        )[0, 1]
    ),
}

print(
    f"charge_RMS_difference_e = "
    f"{charge_stability['RMS_difference_e']:.16g}"
)
print(
    f"charge_MAE_difference_e = "
    f"{charge_stability['MAE_difference_e']:.16g}"
)
print(
    f"charge_maximum_absolute_difference_e = "
    f"{charge_stability['maximum_absolute_difference_e']:.16g}"
)
print(
    f"charge_pearson_r = "
    f"{charge_stability['pearson_r']:.16g}"
)


print("\n[14] CHARGE SCIENTIFIC SUMMARY")

train_unconstrained_charge_metrics = charge_comparison(
    q0_real,
    train_unconstrained_charges,
    real_rows,
)

train_nonnegative_charge_metrics = charge_comparison(
    q0_real,
    train_nonnegative_B_charges,
    real_rows,
)

print("\nTRAIN_ONLY_LAMBDA4_UNCONSTRAINED")

for name, value in (
    train_unconstrained_charge_metrics.items()
):
    print(f"  {name} = {value}")

print("\nTRAIN_ONLY_LAMBDA4_NONNEGATIVE_B")

for name, value in (
    train_nonnegative_charge_metrics.items()
):
    print(f"  {name} = {value}")


print("\n[15] WRITE OUTPUTS")

output_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_TRAIN_ONLY_REFIT_VALIDATION.csv"
)

charge_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_TRAIN_ONLY_REFIT_CHARGES.csv"
)

output_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_TRAIN_ONLY_REFIT_VALIDATION.json"
)

metric_fieldnames = [
    "validation_rank",
    "model_name",
    "training_point_count",
    "training_RMSE_au",
    "training_MAE_au",
    "training_pearson_r",
    "training_r_squared",
    "training_same_sign_fraction",
    "validation_point_count",
    "validation_RMSE_au",
    "validation_MAE_au",
    "validation_pearson_r",
    "validation_r_squared",
    "validation_same_sign_fraction",
    "generalization_gap_au",
    "validation_to_training_RMSE_ratio",
    "relative_generalization_gap",
]

with output_csv.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=metric_fieldnames,
    )

    writer.writeheader()

    for record in sorted(
        model_results,
        key=lambda item: (
            item[
                "validation_rank"
            ]
        ),
    ):
        writer.writerow(
            {
                "validation_rank": (
                    record[
                        "validation_rank"
                    ]
                ),
                "model_name": (
                    record["model_name"]
                ),
                "training_point_count": (
                    record[
                        "training"
                    ]["point_count"]
                ),
                "training_RMSE_au": (
                    record[
                        "training"
                    ]["RMSE_au"]
                ),
                "training_MAE_au": (
                    record[
                        "training"
                    ]["MAE_au"]
                ),
                "training_pearson_r": (
                    record[
                        "training"
                    ]["pearson_r"]
                ),
                "training_r_squared": (
                    record[
                        "training"
                    ]["r_squared"]
                ),
                "training_same_sign_fraction": (
                    record[
                        "training"
                    ]["same_sign_fraction"]
                ),
                "validation_point_count": (
                    record[
                        "validation"
                    ]["point_count"]
                ),
                "validation_RMSE_au": (
                    record[
                        "validation"
                    ]["RMSE_au"]
                ),
                "validation_MAE_au": (
                    record[
                        "validation"
                    ]["MAE_au"]
                ),
                "validation_pearson_r": (
                    record[
                        "validation"
                    ]["pearson_r"]
                ),
                "validation_r_squared": (
                    record[
                        "validation"
                    ]["r_squared"]
                ),
                "validation_same_sign_fraction": (
                    record[
                        "validation"
                    ]["same_sign_fraction"]
                ),
                "generalization_gap_au": (
                    record[
                        "generalization_gap_au"
                    ]
                ),
                "validation_to_training_RMSE_ratio": (
                    record[
                        "validation_to_training_RMSE_ratio"
                    ]
                ),
                "relative_generalization_gap": (
                    record[
                        "relative_generalization_gap"
                    ]
                ),
            }
        )

charge_fieldnames = [
    "real_atom_sequence_index",
    "original_atom_index_0based",
    "atom_id",
    "element",
    "atom_role",
    "RESP_stage1_charge_e",
    "train_only_unconstrained_charge_e",
    "train_only_nonnegative_B_charge_e",
    "full_grid_nonnegative_B_charge_e",
    "train_minus_full_nonnegative_B_e",
]

with charge_csv.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=charge_fieldnames,
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
                "atom_id": row["atom_id"],
                "element": row["element"],
                "atom_role": row["atom_role"],
                "RESP_stage1_charge_e": (
                    q0_real[index]
                ),
                "train_only_unconstrained_charge_e": (
                    train_unconstrained_charges[
                        index
                    ]
                ),
                "train_only_nonnegative_B_charge_e": (
                    train_nonnegative_B_charges[
                        index
                    ]
                ),
                "full_grid_nonnegative_B_charge_e": (
                    full_grid_nonnegative_B_charges[
                        index
                    ]
                ),
                "train_minus_full_nonnegative_B_e": (
                    charge_delta[index]
                ),
            }
        )


print("\n[16] SCIENTIFIC GATES")

all_kkt_gate = all(
    (
        primal_neutrality_gate,
        primal_inequality_gate,
        dual_feasibility_gate,
        stationarity_gate,
        complementarity_gate,
    )
)

results_by_name = {
    record["model_name"]: record
    for record in model_results
}

train_candidate = results_by_name[
    "TRAIN_ONLY_LAMBDA4_NONNEGATIVE_B"
]

full_candidate = results_by_name[
    "FULL_GRID_LAMBDA4_NONNEGATIVE_B"
]

real37_baseline = results_by_name[
    "REAL37_UNMODIFIED"
]

true_holdout_improvement_gate = (
    train_candidate[
        "validation"
    ]["RMSE_au"]
    < real37_baseline[
        "validation"
    ]["RMSE_au"]
)

train_only_generalization_gap_gate = (
    abs(
        train_candidate[
            "relative_generalization_gap"
        ]
    )
    <= 0.10
)

full_vs_train_validation_stability_gate = (
    abs(
        train_candidate[
            "validation"
        ]["RMSE_au"]
        / full_candidate[
            "validation"
        ]["RMSE_au"]
        - 1.0
    )
    <= 0.05
)

charge_stability_gate = (
    charge_stability[
        "RMS_difference_e"
    ]
    <= 0.05
    and charge_stability[
        "maximum_absolute_difference_e"
    ]
    <= 0.15
    and charge_stability[
        "pearson_r"
    ]
    >= 0.95
)

gates = {
    "source_identity_gate": True,
    "upstream_decision_gate": True,
    "partition_count_gate": (
        training_points == 19868
        and validation_points == 4967
    ),
    "strict_convexity_gate": (
        minimum_eigenvalue > 0.0
    ),
    "train_only_KKT_gate": (
        all_kkt_gate
    ),
    "true_holdout_improvement_gate": (
        true_holdout_improvement_gate
    ),
    "train_only_generalization_gap_gate": (
        train_only_generalization_gap_gate
    ),
    "full_vs_train_validation_stability_gate": (
        full_vs_train_validation_stability_gate
    ),
    "charge_stability_gate": (
        charge_stability_gate
    ),
    "output_metrics_csv_created_gate": (
        output_csv.is_file()
        and output_csv.stat().st_size > 0
    ),
    "output_charge_csv_created_gate": (
        charge_csv.is_file()
        and charge_csv.stat().st_size > 0
    ),
    "validation_subset_not_used_in_fit_gate": True,
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


print("\n[17] WRITE JSON REPORT")

decision = (
    "D039_A14_TRAIN_ONLY_HOLDOUT_REFIT_VALIDATION_PASS_"
    "METHOD_FINAL_REVIEW_AUTHORIZED"
    if all_gates_pass
    else
    "D039_A14_TRAIN_ONLY_HOLDOUT_REFIT_VALIDATION_"
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
    "validation_design": {
        "fit_points": (
            "DETERMINISTIC_80_PERCENT_TRAINING_ONLY"
        ),
        "validation_points_used_in_fit": False,
        "validation_rule": (
            "LEXICOGRAPHIC_SORT_POSITION_MOD_5_EQUALS_ZERO"
        ),
        "regularization_lambda": (
            REGULARIZATION_LAMBDA
        ),
        "neutrality_constraint": True,
        "target_nonnegative_constraint": (
            TARGET_ATOM_ID
        ),
    },
    "source_identity": {
        "A13_json": str(
            a13_json.resolve()
        ),
        "A13_json_sha256": sha256(
            a13_json
        ),
        "A12_csv": str(
            a12_csv.resolve()
        ),
        "A12_csv_sha256": sha256(
            a12_csv
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
    "partition": {
        "training_point_count": (
            training_points
        ),
        "validation_point_count": (
            validation_points
        ),
    },
    "train_only_optimization": {
        "minimum_hessian_eigenvalue": (
            minimum_eigenvalue
        ),
        "maximum_hessian_eigenvalue": (
            maximum_eigenvalue
        ),
        "condition_number": (
            maximum_eigenvalue
            / minimum_eigenvalue
        ),
        "unconstrained_target_charge_e": (
            unconstrained_target_charge
        ),
        "inequality_status": (
            inequality_status
        ),
        "active_target_charge_e": (
            target_charge
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
    "model_results": (
        model_results
    ),
    "charge_stability_vs_full_grid": (
        charge_stability
    ),
    "train_only_unconstrained_charge_metrics": (
        train_unconstrained_charge_metrics
    ),
    "train_only_nonnegative_charge_metrics": (
        train_nonnegative_charge_metrics
    ),
    "gates": gates,
    "authorizations": {
        "train_only_holdout_validation_authorized": (
            all_gates_pass
        ),
        "lambda4_nonnegative_B_final_method_review_authorized": (
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
        "metrics_csv": str(
            output_csv.resolve()
        ),
        "metrics_csv_sha256": sha256(
            output_csv
        ),
        "charges_csv": str(
            charge_csv.resolve()
        ),
        "charges_csv_sha256": sha256(
            charge_csv
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
print(f"charge_csv={charge_csv}")
print(
    f"charge_csv_sha256="
    f"{sha256(charge_csv)}"
)
print(f"output_json={output_json}")
print(
    f"output_json_sha256="
    f"{sha256(output_json)}"
)


print("\n[18] DECISION")

print(f"decision={decision}")

print(
    "train_only_holdout_validation_authorized="
    f"{all_gates_pass}"
)

print(
    "lambda4_nonnegative_B_final_method_review_authorized="
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
    "RESP_stage2_execution_authorized=False"
)

print("charge_adoption_authorized=False")
print("force_field_adoption_authorized=False")
print("=" * 100)

if not all_gates_pass:
    raise SystemExit(2)
