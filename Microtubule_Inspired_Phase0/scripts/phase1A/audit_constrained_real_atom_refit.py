#!/usr/bin/env python3
"""
DAY039 / D039-A5

Numerical conditioning and regularization-path audit for a constrained
37-real-atom electrostatic refit of QM_F06_UPPER_V7A_R1.

Constraint
----------
sum(q_real) = 0 exactly, through a null-space parameterization.

Objective family
----------------
minimize ||A q - V_QM||^2 + lambda ||q - q_RESP1_real||^2
subject to sum(q) = 0.

The lambda=0 case is the unregularized constrained least-squares
solution.

Scientific policy
-----------------
- Candidate solutions are diagnostic only.
- No candidate is adopted.
- Amber RESP Stage 2 is not executed.
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

FEASIBILITY_SCRIPT = (
    ROOT
    / "scripts/phase1A/evaluate_real_atom_refit_feasibility.py"
)

AUTHORIZED_VPOT_SHA256 = (
    "73df47796c29d0b2a88a03d83efcb723"
    "d4f6e583d2e1789e1ea1a43d82c1064d"
)

EXPECTED_ATOM_COUNT = 52
EXPECTED_REAL_ATOM_COUNT = 37
EXPECTED_CAP_COUNT = 15
TARGET_CHARGE_E = 0.0
CHARGE_TOLERANCE_E = 1.0e-10

REGULARIZATION_LAMBDAS = (
    0.0,
    1.0e-8,
    1.0e-7,
    1.0e-6,
    1.0e-5,
    1.0e-4,
    1.0e-3,
    1.0e-2,
    1.0e-1,
    1.0,
    10.0,
    100.0,
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


def comparison_metrics(
    reference: np.ndarray,
    candidate: np.ndarray,
) -> dict[str, float]:
    residual = candidate - reference
    absolute = np.abs(residual)

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

    reference_std = float(
        np.std(reference)
    )

    rmse = float(
        np.sqrt(
            np.mean(residual ** 2)
        )
    )

    return {
        "RMSE_au": rmse,
        "MAE_au": float(
            np.mean(absolute)
        ),
        "maximum_absolute_error_au": float(
            np.max(absolute)
        ),
        "residual_mean_au": float(
            np.mean(residual)
        ),
        "RMSE_over_reference_std": (
            rmse / reference_std
            if reference_std > 0.0
            else float("nan")
        ),
        "pearson_r": pearson_r,
        "r_squared": pearson_r * pearson_r,
        "same_sign_fraction": float(
            np.mean(
                np.sign(reference)
                == np.sign(candidate)
            )
        ),
    }


def charge_metrics(
    original: np.ndarray,
    candidate: np.ndarray,
) -> dict:
    delta = candidate - original

    sign_change_indices = [
        int(index)
        for index in range(len(original))
        if (
            original[index] != 0.0
            and candidate[index] != 0.0
            and math.copysign(
                1.0,
                original[index],
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
        "delta_mean_e": float(
            np.mean(delta)
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
        "sign_change_real_indices": (
            sign_change_indices
        ),
    }


def constrained_solution(
    design_matrix: np.ndarray,
    target: np.ndarray,
    q_reference: np.ndarray,
    null_basis: np.ndarray,
    regularization_lambda: float,
) -> np.ndarray:
    """
    Solve the exactly neutral constrained problem using q = Z y.

    For lambda > 0:
        min ||A Z y - b||² + lambda ||Z y - q_reference||²
    """

    reduced_matrix = design_matrix @ null_basis

    if regularization_lambda == 0.0:
        y, _, _, _ = np.linalg.lstsq(
            reduced_matrix,
            target,
            rcond=None,
        )
    else:
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
                target,
                sqrt_lambda * q_reference,
            )
        )

        y, _, _, _ = np.linalg.lstsq(
            augmented_matrix,
            augmented_target,
            rcond=None,
        )

    q = null_basis @ y

    # Numerical cleanup of the exact linear constraint.
    q -= (
        np.sum(q) / len(q)
    )

    return q


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
print("DAY039 / D039-A5 — CONSTRAINED REAL-ATOM REFIT AUDIT")
print("=" * 100)


print("\n[1] SOURCE IDENTITY")

for path in (
    LATEST_POINTER,
    VPOT,
    FEASIBILITY_SCRIPT,
):
    require_file(path)
    print(f"FOUND  {path}")

observed_vpot_sha256 = sha256(VPOT)

print(f"VPOT_SHA256 = {observed_vpot_sha256}")

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

transferability_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1_TRANSFERABILITY.csv"
)

electrostatic_audit_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_CAP_REDISTRIBUTION_ELECTROSTATIC_AUDIT.json"
)

for path in (
    transferability_csv,
    electrostatic_audit_json,
):
    require_file(path)
    print(f"FOUND  {path}")


print("\n[2] UPSTREAM AUTHORIZATION")

electrostatic_report = load_json(
    electrostatic_audit_json
)

expected_upstream_decision = (
    "D039_A3_CAP_REDISTRIBUTION_ELECTROSTATIC_AUDIT_PASS_"
    "MULTICRITERIA_SELECTION_REVIEW_AUTHORIZED"
)

if (
    electrostatic_report.get("decision")
    != expected_upstream_decision
):
    raise RuntimeError(
        "Unexpected electrostatic-audit decision.\n"
        f"Observed: {electrostatic_report.get('decision')}"
    )

if (
    electrostatic_report.get(
        "authorizations",
        {},
    ).get(
        "multicriteria_strategy_selection_review_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Constrained-refit review is not authorized"
    )

print("upstream_decision_gate       = PASS")
print("constrained_refit_review_gate = PASS")
print("refit_adoption_block_gate    = PASS")


print("\n[3] LOAD VPOT AND ATOM DATA")

vpot = read_orca_vpot(VPOT)

atom_xyz_bohr = np.asarray(
    vpot.atom_coordinates_bohr,
    dtype=float,
)

grid_xyz_bohr = np.asarray(
    vpot.grid_coordinates_bohr,
    dtype=float,
)

quantum_esp_au = np.asarray(
    vpot.grid_potential_au,
    dtype=float,
).reshape(-1)

with transferability_csv.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    atom_rows = list(
        csv.DictReader(handle)
    )

if len(atom_rows) != EXPECTED_ATOM_COUNT:
    raise RuntimeError(
        f"Expected 52 atoms, observed {len(atom_rows)}"
    )

atom_rows = sorted(
    atom_rows,
    key=lambda row: int(
        row["atom_index_0based"]
    ),
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

cap_mask = ~real_mask

real_rows = [
    row
    for row, keep in zip(
        atom_rows,
        real_mask,
    )
    if keep
]

real_xyz_bohr = atom_xyz_bohr[
    real_mask
]

q0_real = np.asarray(
    [
        float(
            row["RESP_stage1_charge_e"]
        )
        for row in real_rows
    ],
    dtype=float,
)

if len(real_rows) != EXPECTED_REAL_ATOM_COUNT:
    raise RuntimeError(
        "Expected 37 retained real atoms"
    )

if int(np.sum(cap_mask)) != EXPECTED_CAP_COUNT:
    raise RuntimeError(
        "Expected 15 artificial caps"
    )

print(f"real_atom_count = {len(real_rows)}")
print(f"grid_point_count = {len(grid_xyz_bohr)}")
print(
    f"initial_real_charge_sum_e = "
    f"{np.sum(q0_real):.16g}"
)
print(
    f"target_charge_e = "
    f"{TARGET_CHARGE_E:.16g}"
)


print("\n[4] BUILD COULOMB DESIGN MATRIX")

design_matrix = build_coulomb_matrix(
    grid_xyz_bohr,
    real_xyz_bohr,
)

if design_matrix.shape != (
    24835,
    37,
):
    raise RuntimeError(
        "Unexpected Coulomb matrix shape"
    )

print(
    f"design_matrix_shape = "
    f"{design_matrix.shape}"
)
print(
    f"design_matrix_bytes = "
    f"{design_matrix.nbytes}"
)


print("\n[5] EXACT-NEUTRALITY NULL SPACE")

constraint_row = np.ones(
    (
        1,
        EXPECTED_REAL_ATOM_COUNT,
    ),
    dtype=float,
)

_, constraint_singular_values, constraint_vh = (
    np.linalg.svd(
        constraint_row,
        full_matrices=True,
    )
)

null_basis = constraint_vh[1:, :].T

null_residual = (
    constraint_row @ null_basis
)

orthogonality_residual = (
    null_basis.T @ null_basis
    - np.eye(
        EXPECTED_REAL_ATOM_COUNT - 1
    )
)

print(
    f"null_basis_shape = "
    f"{null_basis.shape}"
)
print(
    f"maximum_constraint_null_residual = "
    f"{np.max(np.abs(null_residual)):.16g}"
)
print(
    f"maximum_null_orthogonality_residual = "
    f"{np.max(np.abs(orthogonality_residual)):.16g}"
)


print("\n[6] NUMERICAL CONDITIONING")

reduced_matrix = (
    design_matrix @ null_basis
)

singular_values = np.linalg.svd(
    reduced_matrix,
    compute_uv=False,
)

largest_singular_value = float(
    singular_values[0]
)

smallest_singular_value = float(
    singular_values[-1]
)

machine_epsilon = np.finfo(float).eps

rank_tolerance = (
    largest_singular_value
    * max(reduced_matrix.shape)
    * machine_epsilon
)

numerical_rank = int(
    np.sum(
        singular_values > rank_tolerance
    )
)

condition_number = (
    largest_singular_value
    / smallest_singular_value
    if smallest_singular_value > 0.0
    else float("inf")
)

print(
    f"reduced_matrix_shape = "
    f"{reduced_matrix.shape}"
)
print(
    f"largest_singular_value = "
    f"{largest_singular_value:.16g}"
)
print(
    f"smallest_singular_value = "
    f"{smallest_singular_value:.16g}"
)
print(
    f"rank_tolerance = "
    f"{rank_tolerance:.16g}"
)
print(
    f"numerical_rank = "
    f"{numerical_rank}/36"
)
print(
    f"condition_number = "
    f"{condition_number:.16g}"
)

print("\nSmallest ten singular values:")

for index, value in enumerate(
    singular_values[-10:],
    start=len(singular_values) - 9,
):
    print(
        f"  singular_index={index:>2} "
        f"value={value:.16g}"
    )


print("\n[7] REGULARIZATION PATH")

candidate_records = []
candidate_charge_vectors = {}

for regularization_lambda in (
    REGULARIZATION_LAMBDAS
):
    candidate = constrained_solution(
        design_matrix,
        quantum_esp_au,
        q0_real,
        null_basis,
        regularization_lambda,
    )

    predicted_potential = (
        design_matrix @ candidate
    )

    electrostatic_metrics = (
        comparison_metrics(
            quantum_esp_au,
            predicted_potential,
        )
    )

    candidate_charge_metrics = (
        charge_metrics(
            q0_real,
            candidate,
        )
    )

    record = {
        "regularization_lambda": (
            regularization_lambda
        ),
        "electrostatic": (
            electrostatic_metrics
        ),
        "charges": (
            candidate_charge_metrics
        ),
    }

    candidate_records.append(record)

    candidate_charge_vectors[
        regularization_lambda
    ] = candidate

    print(
        f"\nlambda={regularization_lambda:.1e}"
    )

    print(
        f"  RMSE_au = "
        f"{electrostatic_metrics['RMSE_au']:.16g}"
    )
    print(
        f"  MAE_au = "
        f"{electrostatic_metrics['MAE_au']:.16g}"
    )
    print(
        f"  max_abs_error_au = "
        f"{electrostatic_metrics['maximum_absolute_error_au']:.16g}"
    )
    print(
        f"  pearson_r = "
        f"{electrostatic_metrics['pearson_r']:.16g}"
    )
    print(
        f"  same_sign_fraction = "
        f"{electrostatic_metrics['same_sign_fraction']:.16g}"
    )
    print(
        f"  charge_sum_e = "
        f"{candidate_charge_metrics['charge_sum_e']:.16g}"
    )
    print(
        f"  max_abs_charge_e = "
        f"{candidate_charge_metrics['maximum_absolute_charge_e']:.16g}"
    )
    print(
        f"  delta_RMS_e = "
        f"{candidate_charge_metrics['delta_RMS_e']:.16g}"
    )
    print(
        f"  delta_max_abs_e = "
        f"{candidate_charge_metrics['delta_max_abs_e']:.16g}"
    )
    print(
        f"  sign_change_count = "
        f"{candidate_charge_metrics['sign_change_count']}"
    )


print("\n[8] RELATIVE RANKINGS")

rmse_ranking = sorted(
    candidate_records,
    key=lambda record: (
        record["electrostatic"]["RMSE_au"]
    ),
)

perturbation_ranking = sorted(
    candidate_records,
    key=lambda record: (
        record["charges"]["delta_RMS_e"]
    ),
)

print("\nElectrostatic RMSE ranking:")

for rank, record in enumerate(
    rmse_ranking,
    start=1,
):
    print(
        f"rank={rank:>2} "
        f"lambda={record['regularization_lambda']:.1e} "
        f"RMSE={record['electrostatic']['RMSE_au']:.10g} "
        f"delta_RMS={record['charges']['delta_RMS_e']:.10g} "
        f"max_abs_q="
        f"{record['charges']['maximum_absolute_charge_e']:.10g}"
    )

print("\nCharge-perturbation ranking:")

for rank, record in enumerate(
    perturbation_ranking,
    start=1,
):
    print(
        f"rank={rank:>2} "
        f"lambda={record['regularization_lambda']:.1e} "
        f"delta_RMS={record['charges']['delta_RMS_e']:.10g} "
        f"RMSE={record['electrostatic']['RMSE_au']:.10g} "
        f"max_abs_q="
        f"{record['charges']['maximum_absolute_charge_e']:.10g}"
    )


print("\n[9] WRITE CANDIDATE CHARGE TABLE")

candidate_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_CONSTRAINED_REFIT_REGULARIZATION_PATH.csv"
)

fieldnames = [
    "real_atom_sequence_index",
    "original_atom_index_0based",
    "original_atom_index_1based",
    "atom_id",
    "element",
    "atom_role",
    "RESP_stage1_charge_e",
]

lambda_column_names = {}

for regularization_lambda in (
    REGULARIZATION_LAMBDAS
):
    label = (
        f"lambda_{regularization_lambda:.0e}"
        .replace("+", "")
        .replace("-", "m")
    )

    column_name = (
        f"constrained_refit_{label}_charge_e"
    )

    lambda_column_names[
        regularization_lambda
    ] = column_name

    fieldnames.append(column_name)

with candidate_csv.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
    )

    writer.writeheader()

    for real_index, row in enumerate(
        real_rows
    ):
        output_row = {
            "real_atom_sequence_index": (
                real_index
            ),
            "original_atom_index_0based": int(
                row["atom_index_0based"]
            ),
            "original_atom_index_1based": int(
                row["atom_index_1based"]
            ),
            "atom_id": row["atom_id"],
            "element": row["element"],
            "atom_role": row["atom_role"],
            "RESP_stage1_charge_e": (
                q0_real[real_index]
            ),
        }

        for regularization_lambda in (
            REGULARIZATION_LAMBDAS
        ):
            output_row[
                lambda_column_names[
                    regularization_lambda
                ]
            ] = candidate_charge_vectors[
                regularization_lambda
            ][real_index]

        writer.writerow(output_row)


print("\n[10] SCIENTIFIC GATES")

full_rank_gate = (
    numerical_rank
    == EXPECTED_REAL_ATOM_COUNT - 1
)

finite_condition_gate = (
    math.isfinite(
        condition_number
    )
)

null_space_gate = (
    np.max(
        np.abs(
            null_residual
        )
    )
    <= 1.0e-12
    and np.max(
        np.abs(
            orthogonality_residual
        )
    )
    <= 1.0e-12
)

candidate_finiteness_gate = all(
    np.all(
        np.isfinite(candidate)
    )
    for candidate in (
        candidate_charge_vectors.values()
    )
)

candidate_neutrality_gate = all(
    abs(
        float(
            np.sum(candidate)
        )
        - TARGET_CHARGE_E
    )
    <= CHARGE_TOLERANCE_E
    for candidate in (
        candidate_charge_vectors.values()
    )
)

unregularized_improves_real37_gate = (
    candidate_records[0][
        "electrostatic"
    ]["RMSE_au"]
    < 0.05621033540022653
)

gates = {
    "source_identity_gate": True,
    "upstream_decision_gate": True,
    "design_matrix_shape_gate": (
        design_matrix.shape
        == (24835, 37)
    ),
    "null_space_gate": null_space_gate,
    "full_reduced_rank_gate": (
        full_rank_gate
    ),
    "finite_condition_gate": (
        finite_condition_gate
    ),
    "candidate_finiteness_gate": (
        candidate_finiteness_gate
    ),
    "candidate_neutrality_gate": (
        candidate_neutrality_gate
    ),
    "unregularized_improves_real37_gate": (
        unregularized_improves_real37_gate
    ),
    "candidate_csv_created_gate": (
        candidate_csv.is_file()
        and candidate_csv.stat().st_size > 0
    ),
    "no_candidate_adopted_gate": True,
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

report_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_CONSTRAINED_REFIT_AUDIT.json"
)

decision = (
    "D039_A5_CONSTRAINED_REAL_ATOM_REFIT_AUDIT_PASS_"
    "REGULARIZATION_SELECTION_REVIEW_AUTHORIZED"
    if all_gates_pass
    else
    "D039_A5_CONSTRAINED_REAL_ATOM_REFIT_AUDIT_"
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
        "constraint": "SUM_q_EQUALS_ZERO",
        "real_atom_count": (
            EXPECTED_REAL_ATOM_COUNT
        ),
        "grid_point_count": 24835,
    },
    "source_identity": {
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
        "electrostatic_audit_json": str(
            electrostatic_audit_json.resolve()
        ),
        "electrostatic_audit_json_sha256": sha256(
            electrostatic_audit_json
        ),
        "feasibility_script": str(
            FEASIBILITY_SCRIPT.resolve()
        ),
        "feasibility_script_sha256": sha256(
            FEASIBILITY_SCRIPT
        ),
    },
    "conditioning": {
        "reduced_matrix_shape": list(
            reduced_matrix.shape
        ),
        "numerical_rank": numerical_rank,
        "expected_rank": (
            EXPECTED_REAL_ATOM_COUNT - 1
        ),
        "largest_singular_value": (
            largest_singular_value
        ),
        "smallest_singular_value": (
            smallest_singular_value
        ),
        "condition_number": (
            condition_number
        ),
        "rank_tolerance": (
            rank_tolerance
        ),
        "singular_values": [
            float(value)
            for value in singular_values
        ],
    },
    "regularization_path": (
        candidate_records
    ),
    "electrostatic_RMSE_ranking_lambdas": [
        record["regularization_lambda"]
        for record in rmse_ranking
    ],
    "charge_perturbation_ranking_lambdas": [
        record["regularization_lambda"]
        for record in perturbation_ranking
    ],
    "gates": gates,
    "authorizations": {
        "regularization_selection_review_authorized": (
            all_gates_pass
        ),
        "constrained_refit_candidate_adoption_authorized": False,
        "RESP_stage2_execution_authorized": False,
        "charge_adoption_authorized": False,
        "force_field_adoption_authorized": False,
    },
    "outputs": {
        "candidate_charge_csv": str(
            candidate_csv.resolve()
        ),
        "candidate_charge_csv_sha256": sha256(
            candidate_csv
        ),
    },
}

report_json.write_text(
    json.dumps(
        report,
        indent=2,
        default=json_safe_value,
    )
    + "\n",
    encoding="utf-8",
)

print(f"candidate_csv = {candidate_csv}")
print(
    f"candidate_csv_sha256 = "
    f"{sha256(candidate_csv)}"
)
print(f"report_json = {report_json}")
print(
    f"report_json_sha256 = "
    f"{sha256(report_json)}"
)


print("\n[12] DECISION")

print(f"decision = {decision}")
print(
    "regularization_selection_review_authorized = "
    f"{all_gates_pass}"
)
print(
    "constrained_refit_candidate_adoption_authorized = False"
)
print("RESP_stage2_execution_authorized = False")
print("charge_adoption_authorized = False")
print("force_field_adoption_authorized = False")
print("=" * 100)
