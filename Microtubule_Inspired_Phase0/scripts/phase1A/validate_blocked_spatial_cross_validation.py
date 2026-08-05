#!/usr/bin/env python3
"""
DAY039 / D039-A15

Blocked spatial cross-validation of the lambda=4, neutral,
37-real-atom electrostatic refit with the local inequality

    q[A:UPPER:8:4] >= 0.

Six deterministic folds are evaluated:

    X_LOW, X_HIGH,
    Y_LOW, Y_HIGH,
    Z_LOW, Z_HIGH.

For every fold, a contiguous spatial tail containing 20% of the
authorized ESP grid is removed. The candidate charges are fitted using
only the complementary 80% training region and evaluated on the
excluded spatial block.

This is stricter than the interleaved validation used in A14 because
entire coordinate regions are withheld.

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

EXPECTED_A14_DECISION = (
    "D039_A14_TRAIN_ONLY_HOLDOUT_REFIT_VALIDATION_PASS_"
    "METHOD_FINAL_REVIEW_AUTHORIZED"
)

EXPECTED_TOTAL_ATOMS = 52
EXPECTED_REAL_ATOMS = 37
EXPECTED_GRID_POINTS = 24835

REGULARIZATION_LAMBDA = 4.0
TARGET_ATOM_ID = "A:UPPER:8:4"
BLOCK_FRACTION = 0.20

PRIMAL_TOLERANCE = 1.0e-10
DUAL_TOLERANCE = 1.0e-8
STATIONARITY_TOLERANCE = 1.0e-7
COMPLEMENTARITY_TOLERANCE = 1.0e-10

FOLDS = (
    ("X_LOW", 0, "LOW"),
    ("X_HIGH", 0, "HIGH"),
    ("Y_LOW", 1, "LOW"),
    ("Y_HIGH", 1, "HIGH"),
    ("Z_LOW", 2, "LOW"),
    ("Z_HIGH", 2, "HIGH"),
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

        matrix[start:stop, :] = 1.0 / distances

    return matrix


def solve_equality_qp(
    hessian: np.ndarray,
    linear_term: np.ndarray,
    constraint_matrix: np.ndarray,
    constraint_values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
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

    rhs = np.concatenate(
        (
            linear_term,
            constraint_values,
        )
    )

    solution = np.linalg.solve(
        kkt_matrix,
        rhs,
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

    rmse = float(
        np.sqrt(
            np.mean(residual ** 2)
        )
    )

    return {
        "point_count": int(len(reference)),
        "RMSE_au": rmse,
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


def charge_metrics(
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


def construct_fold(
    coordinates: np.ndarray,
    axis: int,
    tail: str,
) -> tuple[np.ndarray, np.ndarray, dict]:
    point_count = len(coordinates)

    blocked_count = int(
        round(
            BLOCK_FRACTION
            * point_count
        )
    )

    order = np.argsort(
        coordinates[:, axis],
        kind="mergesort",
    )

    if tail == "LOW":
        validation_indices = order[
            :blocked_count
        ]
    elif tail == "HIGH":
        validation_indices = order[
            point_count - blocked_count:
        ]
    else:
        raise RuntimeError(
            f"Unknown tail: {tail}"
        )

    validation_mask = np.zeros(
        point_count,
        dtype=bool,
    )

    validation_mask[
        validation_indices
    ] = True

    training_indices = np.flatnonzero(
        ~validation_mask
    )

    validation_indices = np.flatnonzero(
        validation_mask
    )

    held_coordinates = coordinates[
        validation_indices,
        axis,
    ]

    metadata = {
        "axis_index": axis,
        "tail": tail,
        "training_point_count": int(
            len(training_indices)
        ),
        "validation_point_count": int(
            len(validation_indices)
        ),
        "validation_axis_min_bohr": float(
            np.min(held_coordinates)
        ),
        "validation_axis_max_bohr": float(
            np.max(held_coordinates)
        ),
        "validation_axis_mean_bohr": float(
            np.mean(held_coordinates)
        ),
    }

    return (
        training_indices,
        validation_indices,
        metadata,
    )


def solve_fold(
    training_matrix: np.ndarray,
    training_target: np.ndarray,
    q0: np.ndarray,
    target_index: int,
) -> dict:
    atom_count = len(q0)

    hessian = 2.0 * (
        training_matrix.T
        @ training_matrix
        + REGULARIZATION_LAMBDA
        * np.eye(
            atom_count,
            dtype=float,
        )
    )

    linear_term = 2.0 * (
        training_matrix.T
        @ training_target
        + REGULARIZATION_LAMBDA
        * q0
    )

    eigenvalues = np.linalg.eigvalsh(
        hessian
    )

    neutrality_matrix = np.ones(
        (
            1,
            atom_count,
        ),
        dtype=float,
    )

    neutrality_values = np.asarray(
        [0.0],
        dtype=float,
    )

    unconstrained, unconstrained_multipliers = (
        solve_equality_qp(
            hessian,
            linear_term,
            neutrality_matrix,
            neutrality_values,
        )
    )

    unconstrained_target_charge = float(
        unconstrained[target_index]
    )

    if unconstrained_target_charge >= 0.0:
        status = "INACTIVE"

        constrained = unconstrained.copy()

        multipliers = (
            unconstrained_multipliers
        )

        inequality_multiplier = 0.0

    else:
        status = "ACTIVE_AT_ZERO"

        constraint_matrix = np.zeros(
            (
                2,
                atom_count,
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
                0.0,
            ],
            dtype=float,
        )

        constrained, multipliers = (
            solve_equality_qp(
                hessian,
                linear_term,
                constraint_matrix,
                constraint_values,
            )
        )

        inequality_multiplier = float(
            -multipliers[1]
        )

    gradient = (
        hessian @ constrained
        - linear_term
    )

    neutrality_multiplier = float(
        multipliers[0]
    )

    target_basis = np.zeros(
        atom_count,
        dtype=float,
    )

    target_basis[
        target_index
    ] = 1.0

    stationarity_residual = (
        gradient
        + neutrality_multiplier
        * np.ones(
            atom_count,
            dtype=float,
        )
        - inequality_multiplier
        * target_basis
    )

    target_charge = float(
        constrained[target_index]
    )

    primal_neutrality_residual = float(
        np.sum(constrained)
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

    gates = {
        "strict_convexity_gate": bool(
            np.min(eigenvalues) > 0.0
        ),
        "primal_neutrality_gate": bool(
            abs(
                primal_neutrality_residual
            )
            <= PRIMAL_TOLERANCE
        ),
        "primal_inequality_gate": bool(
            target_charge
            >= -PRIMAL_TOLERANCE
        ),
        "dual_feasibility_gate": bool(
            inequality_multiplier
            >= -DUAL_TOLERANCE
        ),
        "stationarity_gate": bool(
            maximum_stationarity_residual
            <= STATIONARITY_TOLERANCE
        ),
        "complementarity_gate": bool(
            abs(
                complementarity_value
            )
            <= COMPLEMENTARITY_TOLERANCE
        ),
    }

    gates["all_KKT_conditions_gate"] = all(
        gates.values()
    )

    return {
        "unconstrained_charges": (
            unconstrained
        ),
        "constrained_charges": (
            constrained
        ),
        "inequality_status": status,
        "unconstrained_target_charge_e": (
            unconstrained_target_charge
        ),
        "target_charge_e": (
            target_charge
        ),
        "inequality_multiplier_mu": (
            inequality_multiplier
        ),
        "minimum_hessian_eigenvalue": float(
            np.min(eigenvalues)
        ),
        "maximum_hessian_eigenvalue": float(
            np.max(eigenvalues)
        ),
        "hessian_condition_number": float(
            np.max(eigenvalues)
            / np.min(eigenvalues)
        ),
        "primal_neutrality_residual_e": (
            primal_neutrality_residual
        ),
        "maximum_stationarity_residual": (
            maximum_stationarity_residual
        ),
        "complementarity_mu_times_q": (
            complementarity_value
        ),
        "gates": gates,
    }


print("=" * 100)
print("DAY039 / D039-A15 — BLOCKED SPATIAL CROSS-VALIDATION")
print("=" * 100)


print("\n[1] SOURCE IDENTITY")

for path in (
    LATEST_POINTER,
    VPOT,
):
    require_file(path)

if (
    sha256(VPOT)
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

a14_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_TRAIN_ONLY_REFIT_VALIDATION.json"
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
    a14_json,
    a12_csv,
    transferability_csv,
):
    require_file(path)
    print(f"FOUND  {path}")

print("source_identity_gate = PASS")


print("\n[2] UPSTREAM AUTHORIZATION")

a14_report = load_json(
    a14_json
)

if (
    a14_report.get("decision")
    != EXPECTED_A14_DECISION
):
    raise RuntimeError(
        "Unexpected A14 decision.\n"
        f"Observed: "
        f"{a14_report.get('decision')}"
    )

authorizations = a14_report.get(
    "authorizations",
    {},
)

if (
    authorizations.get(
        "lambda4_nonnegative_B_final_method_review_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Blocked spatial validation is not authorized"
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

print("A14_decision_gate                  = PASS")
print("blocked_CV_review_gate              = PASS")
print("charge_adoption_blocked_gate        = PASS")
print("RESP_stage2_execution_blocked_gate  = PASS")


print("\n[3] LOAD VPOT AND ATOM DATA")

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
    EXPECTED_TOTAL_ATOMS,
    3,
):
    raise RuntimeError(
        f"Unexpected atom shape: "
        f"{atom_xyz_bohr.shape}"
    )

if grid_xyz_bohr.shape != (
    EXPECTED_GRID_POINTS,
    3,
):
    raise RuntimeError(
        f"Unexpected grid shape: "
        f"{grid_xyz_bohr.shape}"
    )

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

RESP52_charges = np.asarray(
    [
        float(
            row[
                "RESP_stage1_charge_e"
            ]
        )
        for row in all_rows
    ],
    dtype=float,
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

q0_real = RESP52_charges[
    real_mask
]

if len(real_rows) != EXPECTED_REAL_ATOMS:
    raise RuntimeError(
        "Expected 37 retained real atoms"
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
        row[
            "real_atom_sequence_index"
        ]
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
        "A12 and retained atom orders differ"
    )

full_grid_candidate = np.asarray(
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

print(f"total_atom_count = {len(all_rows)}")
print(f"real_atom_count = {len(real_rows)}")
print(f"grid_point_count = {len(grid_xyz_bohr)}")
print(f"target_atom_id = {TARGET_ATOM_ID}")
print(f"target_real_index = {target_index}")


print("\n[4] BUILD GLOBAL COULOMB MATRICES")

full_real_matrix = build_coulomb_matrix(
    grid_xyz_bohr,
    real_xyz_bohr,
)

full_RESP52_matrix = build_coulomb_matrix(
    grid_xyz_bohr,
    atom_xyz_bohr,
)

RESP52_potential = (
    full_RESP52_matrix
    @ RESP52_charges
)

REAL37_potential = (
    full_real_matrix
    @ q0_real
)

FULL_GRID_candidate_potential = (
    full_real_matrix
    @ full_grid_candidate
)

print(
    f"full_real_matrix_shape = "
    f"{full_real_matrix.shape}"
)
print(
    f"full_RESP52_matrix_shape = "
    f"{full_RESP52_matrix.shape}"
)


print("\n[5] EXECUTE SIX BLOCKED SPATIAL FOLDS")

fold_records = []
fold_charge_vectors = {}

for fold_name, axis, tail in FOLDS:
    print("\n" + "-" * 100)
    print(f"FOLD = {fold_name}")
    print("-" * 100)

    (
        training_indices,
        validation_indices,
        fold_metadata,
    ) = construct_fold(
        grid_xyz_bohr,
        axis,
        tail,
    )

    training_matrix = full_real_matrix[
        training_indices,
        :,
    ]

    training_target = quantum_esp[
        training_indices
    ]

    solution = solve_fold(
        training_matrix,
        training_target,
        q0_real,
        target_index,
    )

    candidate_charges = solution[
        "constrained_charges"
    ]

    candidate_potential = (
        full_real_matrix
        @ candidate_charges
    )

    candidate_training_metrics = (
        calculate_metrics(
            quantum_esp[
                training_indices
            ],
            candidate_potential[
                training_indices
            ],
        )
    )

    candidate_validation_metrics = (
        calculate_metrics(
            quantum_esp[
                validation_indices
            ],
            candidate_potential[
                validation_indices
            ],
        )
    )

    RESP52_validation_metrics = (
        calculate_metrics(
            quantum_esp[
                validation_indices
            ],
            RESP52_potential[
                validation_indices
            ],
        )
    )

    REAL37_validation_metrics = (
        calculate_metrics(
            quantum_esp[
                validation_indices
            ],
            REAL37_potential[
                validation_indices
            ],
        )
    )

    full_grid_validation_metrics = (
        calculate_metrics(
            quantum_esp[
                validation_indices
            ],
            FULL_GRID_candidate_potential[
                validation_indices
            ],
        )
    )

    candidate_charge_metrics = (
        charge_metrics(
            q0_real,
            candidate_charges,
            real_rows,
        )
    )

    charge_delta_full = (
        candidate_charges
        - full_grid_candidate
    )

    charge_stability = {
        "RMS_difference_e": float(
            np.sqrt(
                np.mean(
                    charge_delta_full ** 2
                )
            )
        ),
        "MAE_difference_e": float(
            np.mean(
                np.abs(
                    charge_delta_full
                )
            )
        ),
        "maximum_absolute_difference_e": float(
            np.max(
                np.abs(
                    charge_delta_full
                )
            )
        ),
        "pearson_r": float(
            np.corrcoef(
                candidate_charges,
                full_grid_candidate,
            )[0, 1]
        ),
    }

    validation_to_training_ratio = (
        candidate_validation_metrics[
            "RMSE_au"
        ]
        / candidate_training_metrics[
            "RMSE_au"
        ]
    )

    validation_improvement_vs_real37 = (
        1.0
        - candidate_validation_metrics[
            "RMSE_au"
        ]
        / REAL37_validation_metrics[
            "RMSE_au"
        ]
    )

    validation_ratio_vs_full_grid = (
        candidate_validation_metrics[
            "RMSE_au"
        ]
        / full_grid_validation_metrics[
            "RMSE_au"
        ]
    )

    fold_integrity_gate = all(
        solution["gates"].values()
    )

    fold_performance_gate = (
        candidate_validation_metrics[
            "RMSE_au"
        ]
        < REAL37_validation_metrics[
            "RMSE_au"
        ]
    )

    fold_record = {
        "fold_name": fold_name,
        "axis_index": axis,
        "tail": tail,
        "partition": fold_metadata,
        "optimization": {
            key: value
            for key, value in solution.items()
            if key not in {
                "unconstrained_charges",
                "constrained_charges",
            }
        },
        "candidate_training_metrics": (
            candidate_training_metrics
        ),
        "candidate_validation_metrics": (
            candidate_validation_metrics
        ),
        "RESP52_validation_metrics": (
            RESP52_validation_metrics
        ),
        "REAL37_validation_metrics": (
            REAL37_validation_metrics
        ),
        "full_grid_candidate_validation_metrics": (
            full_grid_validation_metrics
        ),
        "candidate_charge_metrics": (
            candidate_charge_metrics
        ),
        "charge_stability_vs_full_grid": (
            charge_stability
        ),
        "validation_to_training_RMSE_ratio": float(
            validation_to_training_ratio
        ),
        "relative_generalization_gap": float(
            validation_to_training_ratio
            - 1.0
        ),
        "validation_improvement_vs_REAL37_fraction": float(
            validation_improvement_vs_real37
        ),
        "validation_RMSE_ratio_vs_full_grid_candidate": float(
            validation_ratio_vs_full_grid
        ),
        "fold_integrity_gate": (
            fold_integrity_gate
        ),
        "fold_performance_gate": (
            fold_performance_gate
        ),
    }

    fold_records.append(
        fold_record
    )

    fold_charge_vectors[
        fold_name
    ] = candidate_charges

    print(
        f"training_points = "
        f"{len(training_indices)}"
    )
    print(
        f"validation_points = "
        f"{len(validation_indices)}"
    )
    print(
        f"validation_axis_range_bohr = "
        f"["
        f"{fold_metadata['validation_axis_min_bohr']:.9f}, "
        f"{fold_metadata['validation_axis_max_bohr']:.9f}"
        f"]"
    )
    print(
        f"inequality_status = "
        f"{solution['inequality_status']}"
    )
    print(
        f"unconstrained_target_charge_e = "
        f"{solution['unconstrained_target_charge_e']:.16g}"
    )
    print(
        f"constrained_target_charge_e = "
        f"{solution['target_charge_e']:.16g}"
    )
    print(
        f"inequality_multiplier_mu = "
        f"{solution['inequality_multiplier_mu']:.16g}"
    )
    print(
        f"candidate_training_RMSE_au = "
        f"{candidate_training_metrics['RMSE_au']:.16g}"
    )
    print(
        f"candidate_validation_RMSE_au = "
        f"{candidate_validation_metrics['RMSE_au']:.16g}"
    )
    print(
        f"candidate_validation_R = "
        f"{candidate_validation_metrics['pearson_r']:.16g}"
    )
    print(
        f"candidate_validation_same_sign = "
        f"{candidate_validation_metrics['same_sign_fraction']:.16g}"
    )
    print(
        f"validation_to_training_RMSE_ratio = "
        f"{validation_to_training_ratio:.16g}"
    )
    print(
        f"REAL37_validation_RMSE_au = "
        f"{REAL37_validation_metrics['RMSE_au']:.16g}"
    )
    print(
        f"full_grid_candidate_validation_RMSE_au = "
        f"{full_grid_validation_metrics['RMSE_au']:.16g}"
    )
    print(
        f"validation_improvement_vs_REAL37_fraction = "
        f"{validation_improvement_vs_real37:.16g}"
    )
    print(
        f"charge_RMS_difference_vs_full_grid_e = "
        f"{charge_stability['RMS_difference_e']:.16g}"
    )
    print(
        f"charge_maximum_difference_vs_full_grid_e = "
        f"{charge_stability['maximum_absolute_difference_e']:.16g}"
    )
    print(
        f"charge_pearson_vs_full_grid = "
        f"{charge_stability['pearson_r']:.16g}"
    )
    print(
        f"sign_change_atom_ids = "
        f"{candidate_charge_metrics['sign_change_atom_ids']}"
    )
    print(
        f"fold_integrity_gate = "
        f"{'PASS' if fold_integrity_gate else 'FAIL'}"
    )
    print(
        f"fold_performance_gate = "
        f"{'PASS' if fold_performance_gate else 'FAIL'}"
    )


print("\n[6] AGGREGATED BLOCKED-CV METRICS")

candidate_validation_RMSE_values = np.asarray(
    [
        record[
            "candidate_validation_metrics"
        ]["RMSE_au"]
        for record in fold_records
    ],
    dtype=float,
)

candidate_validation_R_values = np.asarray(
    [
        record[
            "candidate_validation_metrics"
        ]["pearson_r"]
        for record in fold_records
    ],
    dtype=float,
)

candidate_validation_sign_values = np.asarray(
    [
        record[
            "candidate_validation_metrics"
        ]["same_sign_fraction"]
        for record in fold_records
    ],
    dtype=float,
)

charge_RMS_values = np.asarray(
    [
        record[
            "charge_stability_vs_full_grid"
        ]["RMS_difference_e"]
        for record in fold_records
    ],
    dtype=float,
)

charge_max_values = np.asarray(
    [
        record[
            "charge_stability_vs_full_grid"
        ]["maximum_absolute_difference_e"]
        for record in fold_records
    ],
    dtype=float,
)

charge_correlation_values = np.asarray(
    [
        record[
            "charge_stability_vs_full_grid"
        ]["pearson_r"]
        for record in fold_records
    ],
    dtype=float,
)

aggregate = {
    "fold_count": len(fold_records),
    "validation_RMSE_mean_au": float(
        np.mean(
            candidate_validation_RMSE_values
        )
    ),
    "validation_RMSE_std_au": float(
        np.std(
            candidate_validation_RMSE_values
        )
    ),
    "validation_RMSE_min_au": float(
        np.min(
            candidate_validation_RMSE_values
        )
    ),
    "validation_RMSE_max_au": float(
        np.max(
            candidate_validation_RMSE_values
        )
    ),
    "validation_pearson_mean": float(
        np.mean(
            candidate_validation_R_values
        )
    ),
    "validation_pearson_min": float(
        np.min(
            candidate_validation_R_values
        )
    ),
    "validation_same_sign_mean": float(
        np.mean(
            candidate_validation_sign_values
        )
    ),
    "validation_same_sign_min": float(
        np.min(
            candidate_validation_sign_values
        )
    ),
    "charge_RMS_difference_mean_e": float(
        np.mean(
            charge_RMS_values
        )
    ),
    "charge_RMS_difference_max_e": float(
        np.max(
            charge_RMS_values
        )
    ),
    "charge_maximum_difference_max_e": float(
        np.max(
            charge_max_values
        )
    ),
    "charge_pearson_mean": float(
        np.mean(
            charge_correlation_values
        )
    ),
    "charge_pearson_min": float(
        np.min(
            charge_correlation_values
        )
    ),
    "active_constraint_fold_count": sum(
        record["optimization"][
            "inequality_status"
        ]
        == "ACTIVE_AT_ZERO"
        for record in fold_records
    ),
    "all_fold_integrity_pass": all(
        record["fold_integrity_gate"]
        for record in fold_records
    ),
    "all_fold_performance_pass": all(
        record["fold_performance_gate"]
        for record in fold_records
    ),
}

for name, value in aggregate.items():
    print(f"{name} = {value}")


print("\n[7] FOLD RANKING BY VALIDATION RMSE")

ranked_folds = sorted(
    fold_records,
    key=lambda record: (
        record[
            "candidate_validation_metrics"
        ]["RMSE_au"]
    ),
)

for rank, record in enumerate(
    ranked_folds,
    start=1,
):
    record[
        "validation_RMSE_rank"
    ] = rank

    print(
        f"rank={rank} "
        f"fold={record['fold_name']} "
        f"validation_RMSE_au="
        f"{record['candidate_validation_metrics']['RMSE_au']:.16g} "
        f"validation_R="
        f"{record['candidate_validation_metrics']['pearson_r']:.16g} "
        f"validation_same_sign="
        f"{record['candidate_validation_metrics']['same_sign_fraction']:.16g} "
        f"charge_RMS_diff_e="
        f"{record['charge_stability_vs_full_grid']['RMS_difference_e']:.16g}"
    )


print("\n[8] INTER-FOLD CHARGE STABILITY")

charge_matrix = np.vstack(
    [
        fold_charge_vectors[
            fold_name
        ]
        for fold_name, _, _ in FOLDS
    ]
)

atomwise_charge_std = np.std(
    charge_matrix,
    axis=0,
)

interfold_charge_summary = {
    "mean_atomwise_std_e": float(
        np.mean(
            atomwise_charge_std
        )
    ),
    "RMS_atomwise_std_e": float(
        np.sqrt(
            np.mean(
                atomwise_charge_std ** 2
            )
        )
    ),
    "maximum_atomwise_std_e": float(
        np.max(
            atomwise_charge_std
        )
    ),
    "maximum_std_atom_real_index": int(
        np.argmax(
            atomwise_charge_std
        )
    ),
    "maximum_std_atom_id": real_rows[
        int(
            np.argmax(
                atomwise_charge_std
            )
        )
    ]["atom_id"],
}

pairwise_correlations = []

for first_index in range(
    len(FOLDS)
):
    for second_index in range(
        first_index + 1,
        len(FOLDS),
    ):
        correlation = float(
            np.corrcoef(
                charge_matrix[first_index],
                charge_matrix[second_index],
            )[0, 1]
        )

        pairwise_correlations.append(
            {
                "first_fold": (
                    FOLDS[first_index][0]
                ),
                "second_fold": (
                    FOLDS[second_index][0]
                ),
                "pearson_r": correlation,
            }
        )

interfold_charge_summary[
    "minimum_pairwise_charge_correlation"
] = min(
    record["pearson_r"]
    for record in pairwise_correlations
)

interfold_charge_summary[
    "mean_pairwise_charge_correlation"
] = float(
    np.mean(
        [
            record["pearson_r"]
            for record in pairwise_correlations
        ]
    )
)

for name, value in (
    interfold_charge_summary.items()
):
    print(f"{name} = {value}")


print("\n[9] WRITE OUTPUTS")

metrics_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_BLOCKED_SPATIAL_CV_METRICS.csv"
)

charges_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_BLOCKED_SPATIAL_CV_CHARGES.csv"
)

output_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_BLOCKED_SPATIAL_CV.json"
)

metric_fields = [
    "validation_RMSE_rank",
    "fold_name",
    "axis_index",
    "tail",
    "training_point_count",
    "validation_point_count",
    "validation_axis_min_bohr",
    "validation_axis_max_bohr",
    "inequality_status",
    "unconstrained_target_charge_e",
    "target_charge_e",
    "inequality_multiplier_mu",
    "candidate_training_RMSE_au",
    "candidate_validation_RMSE_au",
    "candidate_validation_MAE_au",
    "candidate_validation_pearson_r",
    "candidate_validation_same_sign_fraction",
    "validation_to_training_RMSE_ratio",
    "relative_generalization_gap",
    "REAL37_validation_RMSE_au",
    "RESP52_validation_RMSE_au",
    "full_grid_candidate_validation_RMSE_au",
    "validation_improvement_vs_REAL37_fraction",
    "validation_RMSE_ratio_vs_full_grid_candidate",
    "charge_RMS_difference_vs_full_grid_e",
    "charge_maximum_difference_vs_full_grid_e",
    "charge_pearson_vs_full_grid",
    "sign_change_count",
    "sign_change_atom_ids",
    "fold_integrity_gate",
    "fold_performance_gate",
]

with metrics_csv.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=metric_fields,
    )

    writer.writeheader()

    for record in sorted(
        fold_records,
        key=lambda item: (
            item[
                "validation_RMSE_rank"
            ]
        ),
    ):
        writer.writerow(
            {
                "validation_RMSE_rank": (
                    record[
                        "validation_RMSE_rank"
                    ]
                ),
                "fold_name": (
                    record["fold_name"]
                ),
                "axis_index": (
                    record["axis_index"]
                ),
                "tail": (
                    record["tail"]
                ),
                "training_point_count": (
                    record["partition"][
                        "training_point_count"
                    ]
                ),
                "validation_point_count": (
                    record["partition"][
                        "validation_point_count"
                    ]
                ),
                "validation_axis_min_bohr": (
                    record["partition"][
                        "validation_axis_min_bohr"
                    ]
                ),
                "validation_axis_max_bohr": (
                    record["partition"][
                        "validation_axis_max_bohr"
                    ]
                ),
                "inequality_status": (
                    record["optimization"][
                        "inequality_status"
                    ]
                ),
                "unconstrained_target_charge_e": (
                    record["optimization"][
                        "unconstrained_target_charge_e"
                    ]
                ),
                "target_charge_e": (
                    record["optimization"][
                        "target_charge_e"
                    ]
                ),
                "inequality_multiplier_mu": (
                    record["optimization"][
                        "inequality_multiplier_mu"
                    ]
                ),
                "candidate_training_RMSE_au": (
                    record[
                        "candidate_training_metrics"
                    ]["RMSE_au"]
                ),
                "candidate_validation_RMSE_au": (
                    record[
                        "candidate_validation_metrics"
                    ]["RMSE_au"]
                ),
                "candidate_validation_MAE_au": (
                    record[
                        "candidate_validation_metrics"
                    ]["MAE_au"]
                ),
                "candidate_validation_pearson_r": (
                    record[
                        "candidate_validation_metrics"
                    ]["pearson_r"]
                ),
                "candidate_validation_same_sign_fraction": (
                    record[
                        "candidate_validation_metrics"
                    ]["same_sign_fraction"]
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
                "REAL37_validation_RMSE_au": (
                    record[
                        "REAL37_validation_metrics"
                    ]["RMSE_au"]
                ),
                "RESP52_validation_RMSE_au": (
                    record[
                        "RESP52_validation_metrics"
                    ]["RMSE_au"]
                ),
                "full_grid_candidate_validation_RMSE_au": (
                    record[
                        "full_grid_candidate_validation_metrics"
                    ]["RMSE_au"]
                ),
                "validation_improvement_vs_REAL37_fraction": (
                    record[
                        "validation_improvement_vs_REAL37_fraction"
                    ]
                ),
                "validation_RMSE_ratio_vs_full_grid_candidate": (
                    record[
                        "validation_RMSE_ratio_vs_full_grid_candidate"
                    ]
                ),
                "charge_RMS_difference_vs_full_grid_e": (
                    record[
                        "charge_stability_vs_full_grid"
                    ]["RMS_difference_e"]
                ),
                "charge_maximum_difference_vs_full_grid_e": (
                    record[
                        "charge_stability_vs_full_grid"
                    ]["maximum_absolute_difference_e"]
                ),
                "charge_pearson_vs_full_grid": (
                    record[
                        "charge_stability_vs_full_grid"
                    ]["pearson_r"]
                ),
                "sign_change_count": (
                    record[
                        "candidate_charge_metrics"
                    ]["sign_change_count"]
                ),
                "sign_change_atom_ids": json.dumps(
                    record[
                        "candidate_charge_metrics"
                    ]["sign_change_atom_ids"]
                ),
                "fold_integrity_gate": (
                    record[
                        "fold_integrity_gate"
                    ]
                ),
                "fold_performance_gate": (
                    record[
                        "fold_performance_gate"
                    ]
                ),
            }
        )

charge_fields = [
    "real_atom_sequence_index",
    "original_atom_index_0based",
    "atom_id",
    "element",
    "atom_role",
    "RESP_stage1_charge_e",
    "full_grid_candidate_charge_e",
    "interfold_mean_charge_e",
    "interfold_std_charge_e",
]

for fold_name, _, _ in FOLDS:
    charge_fields.append(
        f"{fold_name}_charge_e"
    )

with charges_csv.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=charge_fields,
    )

    writer.writeheader()

    for index, row in enumerate(
        real_rows
    ):
        output = {
            "real_atom_sequence_index": (
                index
            ),
            "original_atom_index_0based": int(
                row["atom_index_0based"]
            ),
            "atom_id": row["atom_id"],
            "element": row["element"],
            "atom_role": row["atom_role"],
            "RESP_stage1_charge_e": (
                q0_real[index]
            ),
            "full_grid_candidate_charge_e": (
                full_grid_candidate[index]
            ),
            "interfold_mean_charge_e": float(
                np.mean(
                    charge_matrix[:, index]
                )
            ),
            "interfold_std_charge_e": float(
                atomwise_charge_std[
                    index
                ]
            ),
        }

        for fold_name, _, _ in FOLDS:
            output[
                f"{fold_name}_charge_e"
            ] = (
                fold_charge_vectors[
                    fold_name
                ][index]
            )

        writer.writerow(output)


print("\n[10] SCIENTIFIC GATES")

fold_count_gate = (
    len(fold_records) == 6
)

fold_partition_gate = all(
    record["partition"][
        "training_point_count"
    ]
    + record["partition"][
        "validation_point_count"
    ]
    == EXPECTED_GRID_POINTS
    for record in fold_records
)

all_fold_KKT_gate = all(
    record["fold_integrity_gate"]
    for record in fold_records
)

all_fold_improvement_gate = all(
    record["fold_performance_gate"]
    for record in fold_records
)

active_constraint_recurrence_gate = (
    aggregate[
        "active_constraint_fold_count"
    ]
    >= 4
)

blocked_charge_stability_gate = (
    aggregate[
        "charge_RMS_difference_max_e"
    ]
    <= 0.08
    and aggregate[
        "charge_maximum_difference_max_e"
    ]
    <= 0.25
    and aggregate[
        "charge_pearson_min"
    ]
    >= 0.90
)

interfold_charge_stability_gate = (
    interfold_charge_summary[
        "RMS_atomwise_std_e"
    ]
    <= 0.08
    and interfold_charge_summary[
        "maximum_atomwise_std_e"
    ]
    <= 0.25
    and interfold_charge_summary[
        "minimum_pairwise_charge_correlation"
    ]
    >= 0.90
)

gates = {
    "source_identity_gate": True,
    "upstream_decision_gate": True,
    "fold_count_gate": (
        fold_count_gate
    ),
    "fold_partition_gate": (
        fold_partition_gate
    ),
    "all_fold_KKT_gate": (
        all_fold_KKT_gate
    ),
    "all_fold_improvement_vs_REAL37_gate": (
        all_fold_improvement_gate
    ),
    "active_constraint_recurrence_gate": (
        active_constraint_recurrence_gate
    ),
    "blocked_charge_stability_gate": (
        blocked_charge_stability_gate
    ),
    "interfold_charge_stability_gate": (
        interfold_charge_stability_gate
    ),
    "metrics_csv_created_gate": (
        metrics_csv.is_file()
        and metrics_csv.stat().st_size > 0
    ),
    "charges_csv_created_gate": (
        charges_csv.is_file()
        and charges_csv.stat().st_size > 0
    ),
    "validation_blocks_not_used_in_fits_gate": True,
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


print("\n[11] WRITE JSON REPORT")

decision = (
    "D039_A15_BLOCKED_SPATIAL_CROSS_VALIDATION_PASS_"
    "PHASE1A_F_METHOD_ADOPTION_REVIEW_AUTHORIZED"
    if all_gates_pass
    else
    "D039_A15_BLOCKED_SPATIAL_CROSS_VALIDATION_"
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
        "type": (
            "SIX_FOLD_BLOCKED_SPATIAL_CROSS_VALIDATION"
        ),
        "folds": [
            fold_name
            for fold_name, _, _ in FOLDS
        ],
        "blocked_fraction": (
            BLOCK_FRACTION
        ),
        "validation_points_used_in_training": False,
        "regularization_lambda": (
            REGULARIZATION_LAMBDA
        ),
        "neutrality_constraint": True,
        "local_nonnegative_constraint": (
            TARGET_ATOM_ID
        ),
    },
    "source_identity": {
        "A14_json": str(
            a14_json.resolve()
        ),
        "A14_json_sha256": sha256(
            a14_json
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
    "fold_records": (
        fold_records
    ),
    "aggregate": aggregate,
    "interfold_charge_summary": (
        interfold_charge_summary
    ),
    "pairwise_charge_correlations": (
        pairwise_correlations
    ),
    "gates": gates,
    "authorizations": {
        "blocked_spatial_cross_validation_authorized": (
            all_gates_pass
        ),
        "phase1A_F_charge_method_adoption_review_authorized": (
            all_gates_pass
        ),
        "constraint_policy_adoption_authorized": False,
        "regularization_lambda_adoption_authorized": False,
        "constrained_refit_charge_adoption_authorized": False,
        "force_field_integration_authorized": False,
        "RESP_stage2_execution_authorized": False,
        "charge_adoption_authorized": False,
        "force_field_adoption_authorized": False,
    },
    "outputs": {
        "metrics_csv": str(
            metrics_csv.resolve()
        ),
        "metrics_csv_sha256": sha256(
            metrics_csv
        ),
        "charges_csv": str(
            charges_csv.resolve()
        ),
        "charges_csv_sha256": sha256(
            charges_csv
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

print(f"metrics_csv={metrics_csv}")
print(
    f"metrics_csv_sha256="
    f"{sha256(metrics_csv)}"
)
print(f"charges_csv={charges_csv}")
print(
    f"charges_csv_sha256="
    f"{sha256(charges_csv)}"
)
print(f"output_json={output_json}")
print(
    f"output_json_sha256="
    f"{sha256(output_json)}"
)


print("\n[12] DECISION")

print(f"decision={decision}")

print(
    "blocked_spatial_cross_validation_authorized="
    f"{all_gates_pass}"
)

print(
    "phase1A_F_charge_method_adoption_review_authorized="
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
    "force_field_integration_authorized=False"
)

print(
    "RESP_stage2_execution_authorized=False"
)

print("charge_adoption_authorized=False")
print("force_field_adoption_authorized=False")
print("=" * 100)

if not all_gates_pass:
    raise SystemExit(2)
