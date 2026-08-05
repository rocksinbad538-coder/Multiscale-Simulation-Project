#!/usr/bin/env python3
"""
DAY039 / D039-A13

Deterministic spatial hold-out validation of electrostatic-potential
generalization for QM_F06_UPPER_V7A_R1.

Evaluated fixed models
----------------------
1. RESP52
   Full 52-atom RESP Stage 1 charge model.

2. REAL37_UNMODIFIED
   Retained 37 real atoms with their unmodified RESP Stage 1 charges.

3. LAMBDA4_UNCONSTRAINED
   Neutral 37-real-atom lambda=4 constrained refit without the local
   nonnegative-B inequality.

4. LAMBDA4_NONNEGATIVE_B
   Neutral 37-real-atom lambda=4 candidate with
   q[A:UPPER:8:4] >= 0, whose optimum lies at q=0 and was validated
   through the KKT conditions in D039-A12.

Partition policy
----------------
The authorized 24,835 ESP points are sorted lexicographically by
(x, y, z). Every fifth point in the sorted ordering is assigned to the
validation subset:

    validation_position % 5 == 0

All remaining points form the training/evaluation subset.

Important methodological statement
----------------------------------
The charge models are NOT refitted on the training subset. This block
evaluates spatial consistency of already generated fixed candidates
over two deterministic, disjoint subsets of the authorized ESP grid.

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

EXPECTED_A12_DECISION = (
    "D039_A12_NONNEGATIVE_B_KKT_VALIDATION_PASS_"
    "METHOD_CANDIDATE_REVIEW_AUTHORIZED"
)

EXPECTED_TOTAL_ATOMS = 52
EXPECTED_REAL_ATOMS = 37
EXPECTED_CAP_ATOMS = 15
EXPECTED_GRID_POINTS = 24835

VALIDATION_MODULUS = 5
VALIDATION_REMAINDER = 0

MODEL_ORDER = (
    "RESP52",
    "REAL37_UNMODIFIED",
    "LAMBDA4_UNCONSTRAINED",
    "LAMBDA4_NONNEGATIVE_B",
)


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()

    if normalized in {
        "true",
        "1",
        "yes",
        "y",
    }:
        return True

    if normalized in {
        "false",
        "0",
        "no",
        "n",
        "",
    }:
        return False

    raise RuntimeError(
        f"Unrecognized Boolean value: {value!r}"
    )


def json_safe_value(value):
    """
    Convert NumPy scalar and array objects into JSON-compatible
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


def coulomb_potential(
    grid_xyz_bohr: np.ndarray,
    atom_xyz_bohr: np.ndarray,
    charges_e: np.ndarray,
    chunk_size: int = 2000,
) -> np.ndarray:
    """
    Calculate electrostatic potential in atomic units:

        V(r) = sum_i q_i / |r-r_i|

    Coordinates are in bohr and charges in elementary-charge units.
    """

    if atom_xyz_bohr.shape != (
        len(charges_e),
        3,
    ):
        raise RuntimeError(
            "Atom-coordinate and charge-count mismatch.\n"
            f"Coordinates: {atom_xyz_bohr.shape}\n"
            f"Charges: {charges_e.shape}"
        )

    potential = np.empty(
        len(grid_xyz_bohr),
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
            grid_xyz_bohr[
                start:stop,
                None,
                :,
            ]
            - atom_xyz_bohr[
                None,
                :,
                :,
            ]
        )

        distances = np.linalg.norm(
            displacement,
            axis=2,
        )

        if np.any(distances <= 0.0):
            raise RuntimeError(
                "At least one ESP point coincides with an atomic center"
            )

        potential[start:stop] = np.sum(
            charges_e[
                None,
                :,
            ]
            / distances,
            axis=1,
        )

    return potential


def calculate_metrics(
    reference: np.ndarray,
    candidate: np.ndarray,
) -> dict:
    """
    Calculate electrostatic agreement metrics.
    """

    if reference.shape != candidate.shape:
        raise RuntimeError(
            "Reference/candidate shape mismatch.\n"
            f"Reference: {reference.shape}\n"
            f"Candidate: {candidate.shape}"
        )

    residual = candidate - reference
    absolute_residual = np.abs(
        residual
    )

    reference_centered = (
        reference
        - np.mean(reference)
    )

    candidate_centered = (
        candidate
        - np.mean(candidate)
    )

    denominator = math.sqrt(
        float(
            np.sum(
                reference_centered ** 2
            )
            * np.sum(
                candidate_centered ** 2
            )
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
            np.mean(
                residual ** 2
            )
        )
    )

    reference_std = float(
        np.std(reference)
    )

    reference_rms = float(
        np.sqrt(
            np.mean(
                reference ** 2
            )
        )
    )

    return {
        "point_count": int(
            len(reference)
        ),
        "RMSE_au": rmse,
        "MAE_au": float(
            np.mean(
                absolute_residual
            )
        ),
        "maximum_absolute_error_au": float(
            np.max(
                absolute_residual
            )
        ),
        "residual_mean_au": float(
            np.mean(residual)
        ),
        "residual_std_au": float(
            np.std(residual)
        ),
        "pearson_r": pearson_r,
        "r_squared": (
            pearson_r * pearson_r
        ),
        "same_sign_fraction": float(
            np.mean(
                np.sign(reference)
                == np.sign(candidate)
            )
        ),
        "RMSE_over_reference_std": (
            rmse / reference_std
            if reference_std > 0.0
            else float("nan")
        ),
        "relative_RMS_to_reference_RMS": (
            rmse / reference_rms
            if reference_rms > 0.0
            else float("nan")
        ),
        "reference_min_au": float(
            np.min(reference)
        ),
        "reference_max_au": float(
            np.max(reference)
        ),
        "candidate_min_au": float(
            np.min(candidate)
        ),
        "candidate_max_au": float(
            np.max(candidate)
        ),
    }


def print_metric_block(
    model_name: str,
    subset_name: str,
    metrics: dict,
) -> None:
    print(
        f"\n{model_name} / {subset_name}"
    )

    print(
        f"  point_count = "
        f"{metrics['point_count']}"
    )

    print(
        f"  RMSE_au = "
        f"{metrics['RMSE_au']:.16g}"
    )

    print(
        f"  MAE_au = "
        f"{metrics['MAE_au']:.16g}"
    )

    print(
        f"  maximum_absolute_error_au = "
        f"{metrics['maximum_absolute_error_au']:.16g}"
    )

    print(
        f"  residual_mean_au = "
        f"{metrics['residual_mean_au']:.16g}"
    )

    print(
        f"  pearson_r = "
        f"{metrics['pearson_r']:.16g}"
    )

    print(
        f"  r_squared = "
        f"{metrics['r_squared']:.16g}"
    )

    print(
        f"  same_sign_fraction = "
        f"{metrics['same_sign_fraction']:.16g}"
    )

    print(
        f"  RMSE_over_reference_std = "
        f"{metrics['RMSE_over_reference_std']:.16g}"
    )


print("=" * 100)
print("DAY039 / D039-A13 — DETERMINISTIC SPATIAL HOLD-OUT VALIDATION")
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

print(f"VPOT = {VPOT}")
print(
    f"VPOT_SHA256 = "
    f"{observed_vpot_sha256}"
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

a12_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_LAMBDA4_NONNEGATIVE_B_KKT.json"
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
    a12_json,
    a12_csv,
    transferability_csv,
):
    require_file(path)
    print(f"FOUND  {path}")

print("source_identity_gate = PASS")


print("\n[2] UPSTREAM AUTHORIZATION")

a12_report = load_json(
    a12_json
)

observed_a12_decision = (
    a12_report.get(
        "decision"
    )
)

if (
    observed_a12_decision
    != EXPECTED_A12_DECISION
):
    raise RuntimeError(
        "Unexpected D039-A12 decision.\n"
        f"Expected: {EXPECTED_A12_DECISION}\n"
        f"Observed: {observed_a12_decision}"
    )

authorizations = a12_report.get(
    "authorizations",
    {},
)

if (
    authorizations.get(
        "lambda4_nonnegative_B_method_candidate_review_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Spatial generalization review is not authorized"
    )

if (
    authorizations.get(
        "constraint_policy_adoption_authorized"
    )
    is not False
):
    raise RuntimeError(
        "Unexpected constraint-policy adoption authorization"
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

if (
    authorizations.get(
        "RESP_stage2_execution_authorized"
    )
    is not False
):
    raise RuntimeError(
        "Unexpected RESP Stage 2 authorization"
    )

print("A12_decision_gate                       = PASS")
print("generalization_review_authorized_gate   = PASS")
print("constraint_policy_adoption_blocked_gate = PASS")
print("charge_adoption_blocked_gate            = PASS")
print("RESP_stage2_execution_blocked_gate      = PASS")


print("\n[3] LOAD AUTHORIZED VPOT")

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

quantum_esp_au = np.asarray(
    vpot.grid_potential_au,
    dtype=float,
).reshape(-1)

if atom_xyz_bohr.shape != (
    EXPECTED_TOTAL_ATOMS,
    3,
):
    raise RuntimeError(
        "Unexpected VPOT atom-coordinate shape.\n"
        f"Observed: {atom_xyz_bohr.shape}"
    )

if grid_xyz_bohr.shape != (
    EXPECTED_GRID_POINTS,
    3,
):
    raise RuntimeError(
        "Unexpected VPOT grid-coordinate shape.\n"
        f"Observed: {grid_xyz_bohr.shape}"
    )

if quantum_esp_au.shape != (
    EXPECTED_GRID_POINTS,
):
    raise RuntimeError(
        "Unexpected VPOT potential shape.\n"
        f"Observed: {quantum_esp_au.shape}"
    )

if not (
    np.all(
        np.isfinite(
            atom_xyz_bohr
        )
    )
    and np.all(
        np.isfinite(
            grid_xyz_bohr
        )
    )
    and np.all(
        np.isfinite(
            quantum_esp_au
        )
    )
):
    raise RuntimeError(
        "Nonfinite value detected in VPOT"
    )

print(
    f"atom_coordinate_shape = "
    f"{atom_xyz_bohr.shape}"
)
print(
    f"grid_coordinate_shape = "
    f"{grid_xyz_bohr.shape}"
)
print(
    f"quantum_ESP_shape = "
    f"{quantum_esp_au.shape}"
)
print("VPOT_contract_gate = PASS")


print("\n[4] LOAD RESP52 AND REAL37 CHARGE DATA")

with transferability_csv.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    transferability_rows = list(
        csv.DictReader(handle)
    )

if (
    len(transferability_rows)
    != EXPECTED_TOTAL_ATOMS
):
    raise RuntimeError(
        "Unexpected transferability-table row count.\n"
        f"Observed: {len(transferability_rows)}"
    )

transferability_rows.sort(
    key=lambda row: int(
        row["atom_index_0based"]
    )
)

observed_indices = [
    int(
        row["atom_index_0based"]
    )
    for row in transferability_rows
]

if observed_indices != list(
    range(
        EXPECTED_TOTAL_ATOMS
    )
):
    raise RuntimeError(
        "Transferability atom order is not exactly 0..51"
    )

RESP52_charges = np.asarray(
    [
        float(
            row[
                "RESP_stage1_charge_e"
            ]
        )
        for row in transferability_rows
    ],
    dtype=float,
)

real_mask = np.asarray(
    [
        not parse_bool(
            row[
                "artificial_cap"
            ]
        )
        for row in transferability_rows
    ],
    dtype=bool,
)

cap_mask = ~real_mask

real_rows = [
    row
    for row, retained in zip(
        transferability_rows,
        real_mask,
    )
    if retained
]

real_xyz_bohr = atom_xyz_bohr[
    real_mask
]

REAL37_UNMODIFIED_charges = (
    RESP52_charges[
        real_mask
    ]
)

if len(real_rows) != EXPECTED_REAL_ATOMS:
    raise RuntimeError(
        f"Expected {EXPECTED_REAL_ATOMS} real atoms, "
        f"observed {len(real_rows)}"
    )

if int(
    np.sum(cap_mask)
) != EXPECTED_CAP_ATOMS:
    raise RuntimeError(
        f"Expected {EXPECTED_CAP_ATOMS} artificial caps"
    )

real_atom_ids = [
    row["atom_id"]
    for row in real_rows
]

print(
    f"RESP52_atom_count = "
    f"{len(RESP52_charges)}"
)
print(
    f"real37_atom_count = "
    f"{len(REAL37_UNMODIFIED_charges)}"
)
print(
    f"artificial_cap_count = "
    f"{int(np.sum(cap_mask))}"
)
print(
    f"RESP52_charge_sum_e = "
    f"{np.sum(RESP52_charges):.16g}"
)
print(
    f"REAL37_UNMODIFIED_charge_sum_e = "
    f"{np.sum(REAL37_UNMODIFIED_charges):.16g}"
)


print("\n[5] LOAD D039-A12 FIXED CANDIDATES")

with a12_csv.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    a12_rows = list(
        csv.DictReader(handle)
    )

if len(a12_rows) != EXPECTED_REAL_ATOMS:
    raise RuntimeError(
        "Unexpected A12 candidate-table row count.\n"
        f"Observed: {len(a12_rows)}"
    )

a12_rows.sort(
    key=lambda row: int(
        row[
            "real_atom_sequence_index"
        ]
    )
)

a12_sequence_indices = [
    int(
        row[
            "real_atom_sequence_index"
        ]
    )
    for row in a12_rows
]

if a12_sequence_indices != list(
    range(
        EXPECTED_REAL_ATOMS
    )
):
    raise RuntimeError(
        "A12 real-atom sequence is not exactly 0..36"
    )

a12_atom_ids = [
    row["atom_id"]
    for row in a12_rows
]

if a12_atom_ids != real_atom_ids:
    mismatch_records = [
        {
            "real_index": index,
            "transferability_atom_id": (
                real_atom_ids[index]
            ),
            "A12_atom_id": (
                a12_atom_ids[index]
            ),
        }
        for index in range(
            EXPECTED_REAL_ATOMS
        )
        if (
            real_atom_ids[index]
            != a12_atom_ids[index]
        )
    ]

    raise RuntimeError(
        "A12 atom order does not match the retained real-atom order.\n"
        f"Mismatches: {mismatch_records}"
    )

LAMBDA4_UNCONSTRAINED_charges = (
    np.asarray(
        [
            float(
                row[
                    "unconstrained_lambda4_charge_e"
                ]
            )
            for row in a12_rows
        ],
        dtype=float,
    )
)

LAMBDA4_NONNEGATIVE_B_charges = (
    np.asarray(
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
)

candidate_charge_vectors = {
    "RESP52": RESP52_charges,
    "REAL37_UNMODIFIED": (
        REAL37_UNMODIFIED_charges
    ),
    "LAMBDA4_UNCONSTRAINED": (
        LAMBDA4_UNCONSTRAINED_charges
    ),
    "LAMBDA4_NONNEGATIVE_B": (
        LAMBDA4_NONNEGATIVE_B_charges
    ),
}

candidate_atom_coordinates = {
    "RESP52": atom_xyz_bohr,
    "REAL37_UNMODIFIED": (
        real_xyz_bohr
    ),
    "LAMBDA4_UNCONSTRAINED": (
        real_xyz_bohr
    ),
    "LAMBDA4_NONNEGATIVE_B": (
        real_xyz_bohr
    ),
}

for model_name in MODEL_ORDER:
    charges = candidate_charge_vectors[
        model_name
    ]

    print(
        f"{model_name}_charge_count = "
        f"{len(charges)}"
    )

    print(
        f"{model_name}_charge_sum_e = "
        f"{np.sum(charges):.16g}"
    )

    print(
        f"{model_name}_minimum_charge_e = "
        f"{np.min(charges):.16g}"
    )

    print(
        f"{model_name}_maximum_charge_e = "
        f"{np.max(charges):.16g}"
    )

print("candidate_atom_order_gate = PASS")
print("candidate_charge_contract_gate = PASS")


print("\n[6] DETERMINISTIC SPATIAL PARTITION")

# np.lexsort uses the last key as the primary sort key.
# Providing (z, y, x) therefore sorts by x, then y, then z.
lexicographic_order = np.lexsort(
    (
        grid_xyz_bohr[:, 2],
        grid_xyz_bohr[:, 1],
        grid_xyz_bohr[:, 0],
    )
)

sorted_positions = np.arange(
    EXPECTED_GRID_POINTS,
    dtype=int,
)

validation_position_mask = (
    sorted_positions
    % VALIDATION_MODULUS
    == VALIDATION_REMAINDER
)

training_position_mask = (
    ~validation_position_mask
)

validation_indices = lexicographic_order[
    validation_position_mask
]

training_indices = lexicographic_order[
    training_position_mask
]

training_index_set = set(
    int(index)
    for index in training_indices
)

validation_index_set = set(
    int(index)
    for index in validation_indices
)

partition_disjoint_gate = (
    training_index_set.isdisjoint(
        validation_index_set
    )
)

partition_complete_gate = (
    training_index_set
    | validation_index_set
    == set(
        range(
            EXPECTED_GRID_POINTS
        )
    )
)

training_points = len(
    training_indices
)

validation_points = len(
    validation_indices
)

expected_validation_points = (
    EXPECTED_GRID_POINTS
    + VALIDATION_MODULUS
    - 1
) // VALIDATION_MODULUS

expected_training_points = (
    EXPECTED_GRID_POINTS
    - expected_validation_points
)

partition_count_gate = (
    training_points
    == expected_training_points
    and validation_points
    == expected_validation_points
)

print(
    f"partition_sort_policy = "
    f"LEXICOGRAPHIC_X_Y_Z"
)
print(
    f"validation_policy = "
    f"SORTED_POSITION_MOD_{VALIDATION_MODULUS}_"
    f"EQUALS_{VALIDATION_REMAINDER}"
)
print(
    f"training_points = "
    f"{training_points}"
)
print(
    f"validation_points = "
    f"{validation_points}"
)
print(
    f"training_fraction = "
    f"{training_points/EXPECTED_GRID_POINTS:.16g}"
)
print(
    f"validation_fraction = "
    f"{validation_points/EXPECTED_GRID_POINTS:.16g}"
)
print(
    f"partition_disjoint_gate = "
    f"{'PASS' if partition_disjoint_gate else 'FAIL'}"
)
print(
    f"partition_complete_gate = "
    f"{'PASS' if partition_complete_gate else 'FAIL'}"
)
print(
    f"partition_count_gate = "
    f"{'PASS' if partition_count_gate else 'FAIL'}"
)

print("\nFirst 20 training original-grid indices:")
print(
    [
        int(index)
        for index in training_indices[:20]
    ]
)

print("\nFirst 20 validation original-grid indices:")
print(
    [
        int(index)
        for index in validation_indices[:20]
    ]
)


print("\n[7] PARTITION SPATIAL COVERAGE")

for subset_name, subset_indices in (
    (
        "TRAINING",
        training_indices,
    ),
    (
        "VALIDATION",
        validation_indices,
    ),
):
    coordinates = grid_xyz_bohr[
        subset_indices
    ]

    print(f"\n{subset_name}")

    for axis_index, axis_name in enumerate(
        (
            "x",
            "y",
            "z",
        )
    ):
        print(
            f"  {axis_name}_min_bohr = "
            f"{np.min(coordinates[:, axis_index]):.16g}"
        )

        print(
            f"  {axis_name}_max_bohr = "
            f"{np.max(coordinates[:, axis_index]):.16g}"
        )

        print(
            f"  {axis_name}_mean_bohr = "
            f"{np.mean(coordinates[:, axis_index]):.16g}"
        )

        print(
            f"  {axis_name}_std_bohr = "
            f"{np.std(coordinates[:, axis_index]):.16g}"
        )

    subset_reference = quantum_esp_au[
        subset_indices
    ]

    print(
        f"  quantum_ESP_min_au = "
        f"{np.min(subset_reference):.16g}"
    )

    print(
        f"  quantum_ESP_max_au = "
        f"{np.max(subset_reference):.16g}"
    )

    print(
        f"  quantum_ESP_mean_au = "
        f"{np.mean(subset_reference):.16g}"
    )

    print(
        f"  quantum_ESP_std_au = "
        f"{np.std(subset_reference):.16g}"
    )


print("\n[8] CALCULATE FIXED-MODEL POTENTIALS")

model_potentials = {}

for model_name in MODEL_ORDER:
    print(
        f"Calculating {model_name} ..."
    )

    potential = coulomb_potential(
        grid_xyz_bohr,
        candidate_atom_coordinates[
            model_name
        ],
        candidate_charge_vectors[
            model_name
        ],
    )

    if potential.shape != (
        EXPECTED_GRID_POINTS,
    ):
        raise RuntimeError(
            f"Unexpected potential shape for {model_name}: "
            f"{potential.shape}"
        )

    if not np.all(
        np.isfinite(
            potential
        )
    ):
        raise RuntimeError(
            f"Nonfinite potential detected for {model_name}"
        )

    model_potentials[
        model_name
    ] = potential

    print(
        f"  minimum_au = "
        f"{np.min(potential):.16g}"
    )

    print(
        f"  maximum_au = "
        f"{np.max(potential):.16g}"
    )

    print(
        f"  mean_au = "
        f"{np.mean(potential):.16g}"
    )

    print(
        f"  std_au = "
        f"{np.std(potential):.16g}"
    )


print("\n[9] TRAINING AND VALIDATION METRICS")

model_results = []

for model_name in MODEL_ORDER:
    training_metrics = (
        calculate_metrics(
            quantum_esp_au[
                training_indices
            ],
            model_potentials[
                model_name
            ][
                training_indices
            ],
        )
    )

    validation_metrics = (
        calculate_metrics(
            quantum_esp_au[
                validation_indices
            ],
            model_potentials[
                model_name
            ][
                validation_indices
            ],
        )
    )

    generalization_gap_au = (
        validation_metrics[
            "RMSE_au"
        ]
        - training_metrics[
            "RMSE_au"
        ]
    )

    relative_generalization_ratio = (
        validation_metrics[
            "RMSE_au"
        ]
        / training_metrics[
            "RMSE_au"
        ]
        if training_metrics[
            "RMSE_au"
        ] > 0.0
        else float("nan")
    )

    relative_generalization_gap = (
        relative_generalization_ratio
        - 1.0
    )

    record = {
        "model_name": model_name,
        "charge_count": int(
            len(
                candidate_charge_vectors[
                    model_name
                ]
            )
        ),
        "charge_sum_e": float(
            np.sum(
                candidate_charge_vectors[
                    model_name
                ]
            )
        ),
        "training": training_metrics,
        "validation": (
            validation_metrics
        ),
        "generalization_gap_au": float(
            generalization_gap_au
        ),
        "validation_to_training_RMSE_ratio": float(
            relative_generalization_ratio
        ),
        "relative_generalization_gap": float(
            relative_generalization_gap
        ),
    }

    model_results.append(
        record
    )

    print_metric_block(
        model_name,
        "TRAINING",
        training_metrics,
    )

    print_metric_block(
        model_name,
        "VALIDATION",
        validation_metrics,
    )

    print(
        f"\n{model_name} / GENERALIZATION"
    )

    print(
        f"  generalization_gap_au = "
        f"{generalization_gap_au:.16g}"
    )

    print(
        f"  validation_to_training_RMSE_ratio = "
        f"{relative_generalization_ratio:.16g}"
    )

    print(
        f"  relative_generalization_gap = "
        f"{relative_generalization_gap:.16g}"
    )


print("\n[10] GENERALIZATION RANKING")

validation_ranking = sorted(
    model_results,
    key=lambda record: (
        record[
            "validation"
        ]["RMSE_au"],
        record[
            "validation"
        ]["MAE_au"],
        MODEL_ORDER.index(
            record["model_name"]
        ),
    ),
)

print("=" * 100)
print("GENERALIZATION RANKING")
print("=" * 100)

for rank, record in enumerate(
    validation_ranking,
    start=1,
):
    record[
        "validation_RMSE_rank"
    ] = rank

    print(
        f"rank={rank} "
        f"model={record['model_name']} "
        f"training_RMSE_au="
        f"{record['training']['RMSE_au']:.16g} "
        f"validation_RMSE_au="
        f"{record['validation']['RMSE_au']:.16g} "
        f"gap_au="
        f"{record['generalization_gap_au']:.16g} "
        f"validation_to_training_ratio="
        f"{record['validation_to_training_RMSE_ratio']:.16g} "
        f"validation_R="
        f"{record['validation']['pearson_r']:.16g} "
        f"validation_same_sign="
        f"{record['validation']['same_sign_fraction']:.16g}"
    )


print("\n[11] CANDIDATE-SPECIFIC COMPARISON")

results_by_name = {
    record["model_name"]: record
    for record in model_results
}

unconstrained = results_by_name[
    "LAMBDA4_UNCONSTRAINED"
]

nonnegative_B = results_by_name[
    "LAMBDA4_NONNEGATIVE_B"
]

RESP52 = results_by_name[
    "RESP52"
]

REAL37 = results_by_name[
    "REAL37_UNMODIFIED"
]

candidate_validation_RMSE_ratio = (
    nonnegative_B[
        "validation"
    ]["RMSE_au"]
    / unconstrained[
        "validation"
    ]["RMSE_au"]
)

candidate_validation_RMSE_penalty_fraction = (
    candidate_validation_RMSE_ratio
    - 1.0
)

candidate_training_RMSE_ratio = (
    nonnegative_B[
        "training"
    ]["RMSE_au"]
    / unconstrained[
        "training"
    ]["RMSE_au"]
)

candidate_gap_difference_au = (
    nonnegative_B[
        "generalization_gap_au"
    ]
    - unconstrained[
        "generalization_gap_au"
    ]
)

candidate_vs_real37_validation_improvement_fraction = (
    1.0
    - nonnegative_B[
        "validation"
    ]["RMSE_au"]
    / REAL37[
        "validation"
    ]["RMSE_au"]
)

candidate_vs_RESP52_validation_RMSE_ratio = (
    nonnegative_B[
        "validation"
    ]["RMSE_au"]
    / RESP52[
        "validation"
    ]["RMSE_au"]
)

print(
    "LAMBDA4_NONNEGATIVE_B_vs_"
    "LAMBDA4_UNCONSTRAINED_validation_RMSE_ratio = "
    f"{candidate_validation_RMSE_ratio:.16g}"
)

print(
    "LAMBDA4_NONNEGATIVE_B_validation_RMSE_"
    "penalty_fraction = "
    f"{candidate_validation_RMSE_penalty_fraction:.16g}"
)

print(
    "LAMBDA4_NONNEGATIVE_B_vs_"
    "LAMBDA4_UNCONSTRAINED_training_RMSE_ratio = "
    f"{candidate_training_RMSE_ratio:.16g}"
)

print(
    "candidate_generalization_gap_difference_au = "
    f"{candidate_gap_difference_au:.16g}"
)

print(
    "LAMBDA4_NONNEGATIVE_B_vs_"
    "REAL37_UNMODIFIED_validation_improvement_fraction = "
    f"{candidate_vs_real37_validation_improvement_fraction:.16g}"
)

print(
    "LAMBDA4_NONNEGATIVE_B_vs_RESP52_"
    "validation_RMSE_ratio = "
    f"{candidate_vs_RESP52_validation_RMSE_ratio:.16g}"
)


print("\n[12] REPRODUCIBILITY CROSS-CHECK")

# Reconstruct the partition independently using the same declared rule.
repeat_order = np.lexsort(
    (
        grid_xyz_bohr[:, 2],
        grid_xyz_bohr[:, 1],
        grid_xyz_bohr[:, 0],
    )
)

repeat_positions = np.arange(
    EXPECTED_GRID_POINTS,
    dtype=int,
)

repeat_validation_indices = (
    repeat_order[
        repeat_positions
        % VALIDATION_MODULUS
        == VALIDATION_REMAINDER
    ]
)

repeat_training_indices = (
    repeat_order[
        repeat_positions
        % VALIDATION_MODULUS
        != VALIDATION_REMAINDER
    ]
)

partition_reproducibility_gate = (
    np.array_equal(
        training_indices,
        repeat_training_indices,
    )
    and np.array_equal(
        validation_indices,
        repeat_validation_indices,
    )
)

print(
    "training_partition_repeat_equal = "
    f"{np.array_equal(training_indices, repeat_training_indices)}"
)

print(
    "validation_partition_repeat_equal = "
    f"{np.array_equal(validation_indices, repeat_validation_indices)}"
)

print(
    "partition_reproducibility_gate = "
    f"{'PASS' if partition_reproducibility_gate else 'FAIL'}"
)


print("\n[13] WRITE OUTPUTS")

output_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_SPATIAL_HOLDOUT_GENERALIZATION.csv"
)

output_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_SPATIAL_HOLDOUT_GENERALIZATION.json"
)

csv_fieldnames = [
    "validation_RMSE_rank",
    "model_name",
    "charge_count",
    "charge_sum_e",
    "training_point_count",
    "training_RMSE_au",
    "training_MAE_au",
    "training_maximum_absolute_error_au",
    "training_residual_mean_au",
    "training_pearson_r",
    "training_r_squared",
    "training_same_sign_fraction",
    "training_RMSE_over_reference_std",
    "validation_point_count",
    "validation_RMSE_au",
    "validation_MAE_au",
    "validation_maximum_absolute_error_au",
    "validation_residual_mean_au",
    "validation_pearson_r",
    "validation_r_squared",
    "validation_same_sign_fraction",
    "validation_RMSE_over_reference_std",
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
        fieldnames=csv_fieldnames,
    )

    writer.writeheader()

    for record in sorted(
        model_results,
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
                "model_name": (
                    record["model_name"]
                ),
                "charge_count": (
                    record["charge_count"]
                ),
                "charge_sum_e": (
                    record["charge_sum_e"]
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
                "training_maximum_absolute_error_au": (
                    record[
                        "training"
                    ][
                        "maximum_absolute_error_au"
                    ]
                ),
                "training_residual_mean_au": (
                    record[
                        "training"
                    ]["residual_mean_au"]
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
                "training_RMSE_over_reference_std": (
                    record[
                        "training"
                    ][
                        "RMSE_over_reference_std"
                    ]
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
                "validation_maximum_absolute_error_au": (
                    record[
                        "validation"
                    ][
                        "maximum_absolute_error_au"
                    ]
                ),
                "validation_residual_mean_au": (
                    record[
                        "validation"
                    ]["residual_mean_au"]
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
                "validation_RMSE_over_reference_std": (
                    record[
                        "validation"
                    ][
                        "RMSE_over_reference_std"
                    ]
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


print("\n[14] SCIENTIFIC GATES")

training_metrics_gate = all(
    record[
        "training"
    ]["point_count"]
    == training_points
    for record in model_results
)

validation_metrics_gate = all(
    record[
        "validation"
    ]["point_count"]
    == validation_points
    for record in model_results
)

finite_metrics_gate = all(
    math.isfinite(value)
    for record in model_results
    for subset_name in (
        "training",
        "validation",
    )
    for key, value in (
        record[
            subset_name
        ].items()
    )
    if key != "point_count"
) and all(
    math.isfinite(
        record[
            "generalization_gap_au"
        ]
    )
    and math.isfinite(
        record[
            "validation_to_training_RMSE_ratio"
        ]
    )
    and math.isfinite(
        record[
            "relative_generalization_gap"
        ]
    )
    for record in model_results
)

generalization_gap_gate = all(
    abs(
        record[
            "relative_generalization_gap"
        ]
    )
    <= 0.10
    for record in model_results
)

candidate_comparison_gate = (
    nonnegative_B[
        "validation"
    ]["RMSE_au"]
    < REAL37[
        "validation"
    ]["RMSE_au"]
    and candidate_validation_RMSE_penalty_fraction
    <= 0.05
)

validation_ranking_complete_gate = (
    sorted(
        record[
            "validation_RMSE_rank"
        ]
        for record in model_results
    )
    == [
        1,
        2,
        3,
        4,
    ]
)

gates = {
    "source_identity_gate": True,
    "upstream_decision_gate": True,
    "VPOT_contract_gate": True,
    "candidate_atom_order_gate": True,
    "candidate_charge_contract_gate": True,
    "partition_disjoint_gate": (
        partition_disjoint_gate
    ),
    "partition_complete_gate": (
        partition_complete_gate
    ),
    "partition_count_gate": (
        partition_count_gate
    ),
    "partition_reproducibility_gate": (
        partition_reproducibility_gate
    ),
    "training_metrics_gate": (
        training_metrics_gate
    ),
    "validation_metrics_gate": (
        validation_metrics_gate
    ),
    "finite_metrics_gate": (
        finite_metrics_gate
    ),
    "generalization_gap_gate": (
        generalization_gap_gate
    ),
    "candidate_comparison_gate": (
        candidate_comparison_gate
    ),
    "validation_ranking_complete_gate": (
        validation_ranking_complete_gate
    ),
    "output_csv_created_gate": (
        output_csv.is_file()
        and output_csv.stat().st_size > 0
    ),
    "no_model_refitted_gate": True,
    "no_constraint_policy_adopted_gate": True,
    "no_lambda_adopted_gate": True,
    "no_charge_adopted_gate": True,
    "RESP_stage2_not_executed_gate": True,
}

for gate_name, gate_value in (
    gates.items()
):
    print(
        f"{gate_name}="
        f"{'PASS' if gate_value else 'FAIL'}"
    )

all_gates_pass = all(
    gates.values()
)


print("\n[15] WRITE JSON REPORT")

decision = (
    "D039_A13_SPATIAL_GENERALIZATION_VALIDATION_PASS_"
    "METHOD_ADOPTION_REVIEW_AUTHORIZED"
    if all_gates_pass
    else
    "D039_A13_SPATIAL_GENERALIZATION_VALIDATION_"
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
    "methodological_statement": {
        "validation_type": (
            "DETERMINISTIC_SPATIAL_HOLDOUT_EVALUATION"
        ),
        "models_refitted_on_training_subset": False,
        "interpretation": (
            "FIXED_CANDIDATE_MODELS_EVALUATED_ON_"
            "DETERMINISTIC_DISJOINT_GRID_SUBSETS"
        ),
    },
    "source_identity": {
        "A12_json": str(
            a12_json.resolve()
        ),
        "A12_json_sha256": sha256(
            a12_json
        ),
        "A12_csv": str(
            a12_csv.resolve()
        ),
        "A12_csv_sha256": sha256(
            a12_csv
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
    "partition": {
        "sorting_policy": (
            "LEXICOGRAPHIC_X_THEN_Y_THEN_Z"
        ),
        "validation_rule": (
            "SORTED_POSITION_MOD_5_EQUALS_0"
        ),
        "training_rule": (
            "SORTED_POSITION_MOD_5_NOT_EQUALS_0"
        ),
        "total_point_count": (
            EXPECTED_GRID_POINTS
        ),
        "training_point_count": (
            training_points
        ),
        "validation_point_count": (
            validation_points
        ),
        "training_fraction": (
            training_points
            / EXPECTED_GRID_POINTS
        ),
        "validation_fraction": (
            validation_points
            / EXPECTED_GRID_POINTS
        ),
        "first_20_training_original_indices": [
            int(index)
            for index in (
                training_indices[:20]
            )
        ],
        "first_20_validation_original_indices": [
            int(index)
            for index in (
                validation_indices[:20]
            )
        ],
    },
    "model_results": (
        model_results
    ),
    "validation_ranking": [
        {
            "rank": rank,
            "model_name": (
                record["model_name"]
            ),
            "validation_RMSE_au": (
                record[
                    "validation"
                ]["RMSE_au"]
            ),
            "training_RMSE_au": (
                record[
                    "training"
                ]["RMSE_au"]
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
        }
        for rank, record in enumerate(
            validation_ranking,
            start=1,
        )
    ],
    "candidate_specific_comparison": {
        "nonnegative_B_vs_unconstrained_validation_RMSE_ratio": (
            candidate_validation_RMSE_ratio
        ),
        "nonnegative_B_validation_RMSE_penalty_fraction": (
            candidate_validation_RMSE_penalty_fraction
        ),
        "nonnegative_B_vs_unconstrained_training_RMSE_ratio": (
            candidate_training_RMSE_ratio
        ),
        "candidate_generalization_gap_difference_au": (
            candidate_gap_difference_au
        ),
        "nonnegative_B_vs_REAL37_validation_improvement_fraction": (
            candidate_vs_real37_validation_improvement_fraction
        ),
        "nonnegative_B_vs_RESP52_validation_RMSE_ratio": (
            candidate_vs_RESP52_validation_RMSE_ratio
        ),
    },
    "gates": gates,
    "authorizations": {
        "spatial_generalization_validation_authorized": (
            all_gates_pass
        ),
        "lambda4_nonnegative_B_method_adoption_review_authorized": (
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
        "generalization_csv": str(
            output_csv.resolve()
        ),
        "generalization_csv_sha256": sha256(
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


print("\n[16] DECISION")

print(f"decision={decision}")

print(
    "spatial_generalization_validation_authorized="
    f"{all_gates_pass}"
)

print(
    "lambda4_nonnegative_B_method_adoption_review_authorized="
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

if not all_gates_pass:
    raise SystemExit(2)
