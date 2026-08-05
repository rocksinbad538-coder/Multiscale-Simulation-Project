#!/usr/bin/env python3
"""
DAY039 / D039-A9

Focused regularization refinement between lambda=1 and lambda=10
for the constrained 37-real-atom electrostatic refit.

No lambda or charge set is adopted.
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

EXPECTED_A8_DECISION = (
    "D039_A8_LAMBDA1_VS_LAMBDA10_COMPARISON_PASS_"
    "FINAL_LAMBDA_REVIEW_AUTHORIZED"
)

AUTHORIZED_VPOT_SHA256 = (
    "73df47796c29d0b2a88a03d83efcb723"
    "d4f6e583d2e1789e1ea1a43d82c1064d"
)

LAMBDA_GRID = (
    1.0,
    1.25,
    1.5,
    2.0,
    2.5,
    3.0,
    4.0,
    5.0,
    6.0,
    7.5,
    10.0,
)

CHARGE_TOLERANCE_E = 1.0e-10


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
        (len(grid_xyz_bohr), len(atom_xyz_bohr)),
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
                "ESP point coincides with an atom"
            )

        matrix[start:stop, :] = 1.0 / distances

    return matrix


def solve_candidate(
    matrix: np.ndarray,
    target: np.ndarray,
    q_reference: np.ndarray,
    null_basis: np.ndarray,
    regularization_lambda: float,
) -> np.ndarray:
    reduced = matrix @ null_basis
    sqrt_lambda = math.sqrt(
        regularization_lambda
    )

    augmented_matrix = np.vstack(
        (
            reduced,
            sqrt_lambda * null_basis,
        )
    )

    augmented_target = np.concatenate(
        (
            target,
            sqrt_lambda * q_reference,
        )
    )

    reduced_solution, _, _, _ = np.linalg.lstsq(
        augmented_matrix,
        augmented_target,
        rcond=None,
    )

    charges = (
        null_basis @ reduced_solution
    )

    charges -= (
        np.sum(charges) / len(charges)
    )

    return charges


def potential_metrics(
    reference: np.ndarray,
    candidate: np.ndarray,
) -> dict:
    residual = candidate - reference

    centered_reference = (
        reference - np.mean(reference)
    )

    centered_candidate = (
        candidate - np.mean(candidate)
    )

    denominator = math.sqrt(
        float(
            np.sum(centered_reference ** 2)
            * np.sum(centered_candidate ** 2)
        )
    )

    pearson = (
        float(
            np.sum(
                centered_reference
                * centered_candidate
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
        "pearson_r": pearson,
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


print("=" * 100)
print("DAY039 / D039-A9 — REGULARIZATION REFINEMENT FROM LAMBDA 1 TO 10")
print("=" * 100)


print("\n[1] SOURCES")

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

a8_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA1_VS_LAMBDA10_SCIENTIFIC_COMPARISON.json"
)

transferability_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1_TRANSFERABILITY.csv"
)

for path in (
    a8_json,
    transferability_csv,
):
    require_file(path)
    print(f"FOUND  {path}")


print("\n[2] UPSTREAM AUTHORIZATION")

a8_report = load_json(
    a8_json
)

if (
    a8_report.get("decision")
    != EXPECTED_A8_DECISION
):
    raise RuntimeError(
        "Unexpected A8 decision"
    )

if (
    a8_report.get(
        "authorizations",
        {},
    ).get(
        "final_lambda_scientific_review_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Focused lambda refinement is not authorized"
    )

print("A8_decision_gate             = PASS")
print("lambda_refinement_gate       = PASS")
print("lambda_adoption_block_gate   = PASS")


print("\n[3] LOAD DATA")

vpot = read_orca_vpot(VPOT)

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
)

with transferability_csv.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    rows = list(
        csv.DictReader(handle)
    )

rows.sort(
    key=lambda row: int(
        row["atom_index_0based"]
    )
)

real_mask = np.asarray(
    [
        not parse_bool(
            row["artificial_cap"]
        )
        for row in rows
    ],
    dtype=bool,
)

real_rows = [
    row
    for row, retained in zip(
        rows,
        real_mask,
    )
    if retained
]

real_xyz = atom_xyz_bohr[
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
        "Expected 37 real atoms"
    )

print(f"real_atom_count = {len(real_rows)}")
print(f"grid_point_count = {len(grid_xyz_bohr)}")


print("\n[4] BUILD CONSTRAINED SYSTEM")

matrix = build_coulomb_matrix(
    grid_xyz_bohr,
    real_xyz,
)

constraint = np.ones(
    (1, 37),
    dtype=float,
)

_, _, constraint_vh = np.linalg.svd(
    constraint,
    full_matrices=True,
)

null_basis = constraint_vh[
    1:, :
].T

print(f"matrix_shape = {matrix.shape}")
print(f"null_basis_shape = {null_basis.shape}")


print("\n[5] REFINED LAMBDA PATH")

candidate_vectors = {}
records = []

for value in LAMBDA_GRID:
    charges = solve_candidate(
        matrix,
        quantum_esp,
        q0,
        null_basis,
        value,
    )

    predicted = matrix @ charges

    electrostatic = potential_metrics(
        quantum_esp,
        predicted,
    )

    charge = charge_metrics(
        q0,
        charges,
    )

    candidate_vectors[value] = charges

    record = {
        "regularization_lambda": value,
        "electrostatic": electrostatic,
        "charges": charge,
    }

    records.append(record)

    sign_change_ids = [
        real_rows[index]["atom_id"]
        for index in charge[
            "sign_change_real_indices"
        ]
    ]

    print(f"\nlambda={value:g}")
    print(
        f"  RMSE_au="
        f"{electrostatic['RMSE_au']:.16g}"
    )
    print(
        f"  pearson_r="
        f"{electrostatic['pearson_r']:.16g}"
    )
    print(
        f"  same_sign_fraction="
        f"{electrostatic['same_sign_fraction']:.16g}"
    )
    print(
        f"  delta_RMS_e="
        f"{charge['delta_RMS_e']:.16g}"
    )
    print(
        f"  delta_max_abs_e="
        f"{charge['delta_max_abs_e']:.16g}"
    )
    print(
        f"  maximum_absolute_charge_e="
        f"{charge['maximum_absolute_charge_e']:.16g}"
    )
    print(
        f"  sign_change_count="
        f"{charge['sign_change_count']}"
    )
    print(
        f"  sign_change_atom_ids="
        f"{sign_change_ids}"
    )


print("\n[6] TRANSITION ANALYSIS")

for previous, current in zip(
    records[:-1],
    records[1:],
):
    previous_signs = set(
        previous["charges"][
            "sign_change_real_indices"
        ]
    )

    current_signs = set(
        current["charges"][
            "sign_change_real_indices"
        ]
    )

    removed = sorted(
        previous_signs - current_signs
    )

    added = sorted(
        current_signs - previous_signs
    )

    print(
        f"lambda "
        f"{previous['regularization_lambda']:g}"
        f" -> "
        f"{current['regularization_lambda']:g}: "
        f"sign_changes "
        f"{len(previous_signs)}"
        f" -> "
        f"{len(current_signs)}, "
        f"removed_indices={removed}, "
        f"added_indices={added}"
    )


print("\n[7] REVIEW FILTER")

review_candidates = []

for record in records:
    electrostatic = record[
        "electrostatic"
    ]

    charge = record[
        "charges"
    ]

    # Project-level review criteria.
    passes = (
        charge[
            "sign_change_count"
        ]
        <= 2
        and charge[
            "maximum_absolute_charge_e"
        ]
        <= 1.0
        and charge[
            "delta_max_abs_e"
        ]
        <= 0.30
        and electrostatic[
            "pearson_r"
        ]
        >= 0.25
        and electrostatic[
            "same_sign_fraction"
        ]
        >= 0.55
    )

    record[
        "focused_review_admissible"
    ] = passes

    if passes:
        review_candidates.append(
            record
        )

    print(
        f"lambda="
        f"{record['regularization_lambda']:g} "
        f"focused_review_admissible="
        f"{passes}"
    )


print("\n[8] WRITE OUTPUTS")

output_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA_1_TO_10_REFINEMENT.csv"
)

output_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA_1_TO_10_REFINEMENT.json"
)

fieldnames = [
    "real_atom_sequence_index",
    "original_atom_index_0based",
    "atom_id",
    "element",
    "atom_role",
    "RESP_stage1_charge_e",
]

for value in LAMBDA_GRID:
    fieldnames.extend(
        (
            f"lambda_{value:g}_charge_e",
            f"lambda_{value:g}_delta_e",
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
            "RESP_stage1_charge_e": q0[index],
        }

        for value in LAMBDA_GRID:
            candidate = (
                candidate_vectors[value][index]
            )

            output[
                f"lambda_{value:g}_charge_e"
            ] = candidate

            output[
                f"lambda_{value:g}_delta_e"
            ] = (
                candidate - q0[index]
            )

        writer.writerow(output)


print("\n[9] GATES")

neutrality_gate = all(
    abs(
        record["charges"]["charge_sum_e"]
    )
    <= CHARGE_TOLERANCE_E
    for record in records
)

finite_gate = all(
    all(
        math.isfinite(value)
        for value in (
            record[
                "electrostatic"
            ].values()
        )
    )
    and all(
        math.isfinite(value)
        for key, value in (
            record["charges"].items()
        )
        if key
        != "sign_change_real_indices"
    )
    for record in records
)

gates = {
    "upstream_decision_gate": True,
    "lambda_grid_count_gate": (
        len(records) == len(LAMBDA_GRID)
    ),
    "neutrality_gate": neutrality_gate,
    "finite_candidate_gate": finite_gate,
    "output_csv_created_gate": (
        output_csv.is_file()
        and output_csv.stat().st_size > 0
    ),
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
    "D039_A9_LAMBDA_1_TO_10_REFINEMENT_PASS_"
    "FOCUSED_CANDIDATE_REVIEW_AUTHORIZED"
    if all_gates_pass
    else
    "D039_A9_LAMBDA_1_TO_10_REFINEMENT_"
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
    "lambda_grid": list(
        LAMBDA_GRID
    ),
    "records": records,
    "focused_review_admissible_lambdas": [
        record[
            "regularization_lambda"
        ]
        for record in review_candidates
    ],
    "review_filter": {
        "maximum_sign_change_count": 2,
        "maximum_absolute_charge_e": 1.0,
        "maximum_delta_absolute_e": 0.30,
        "minimum_pearson_r": 0.25,
        "minimum_same_sign_fraction": 0.55,
        "policy": (
            "PROJECT_LEVEL_REVIEW_FILTERS_"
            "NOT_UNIVERSAL_CONSTANTS"
        ),
    },
    "gates": gates,
    "authorizations": {
        "focused_lambda_candidate_review_authorized": (
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
        "refinement_csv": str(
            output_csv.resolve()
        ),
        "refinement_csv_sha256": sha256(
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
    "focused_review_admissible_lambdas="
    f"{[record['regularization_lambda'] for record in review_candidates]}"
)

print(
    "focused_lambda_candidate_review_authorized="
    f"{all_gates_pass}"
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
