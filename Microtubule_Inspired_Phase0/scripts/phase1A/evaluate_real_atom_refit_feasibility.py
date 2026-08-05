#!/usr/bin/env python3
"""
DAY039 / D039-A4

Evaluate the electrostatic feasibility of refitting charges on the
37 retained real atoms after removing the 15 artificial boundary caps.

Comparisons
-----------
1. Full 52-atom RESP Stage 1 potential versus quantum VPOT.
2. Unmodified 37-real-atom potential versus quantum VPOT.
3. Unmodified 37-real-atom potential versus full 52-atom RESP potential.

No refit is executed.
No charges are modified or adopted.
RESP Stage 2 remains blocked.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from resp_common import (
    read_orca_vpot,
    require_file,
    sha256,
)


ROOT = Path(__file__).resolve().parents[2]

VPOT = (
    ROOT
    / "runs/phase1A/day036_qm_f06_upper_v7a_r1_esp_executions"
    / "esp_upper_v7a_r1_20260731T174832Z"
    / "esp_upper_v7a_r1.vpot"
)

LATEST_POINTER = (
    ROOT
    / "runs/phase1A/day038_resp_stage1_executions"
    / "LATEST_RESP_STAGE1_UPPER_V7A_R1_EXECUTION.txt"
)

AUTHORIZED_VPOT_SHA256 = (
    "73df47796c29d0b2a88a03d83efcb723"
    "d4f6e583d2e1789e1ea1a43d82c1064d"
)

BOHR_TO_ANGSTROM = 0.529177210903


def get_object_attribute(
    obj,
    candidates: tuple[str, ...],
    label: str,
):
    """
    Return the first available non-None object attribute.

    This keeps the diagnostic compatible with the established
    OrcaVpotData API without treating it as a dictionary.
    """

    for name in candidates:
        if hasattr(obj, name):
            value = getattr(obj, name)

            if value is not None:
                print(
                    f"{label}_attribute = {name}"
                )
                return value

    available = [
        name
        for name in dir(obj)
        if not name.startswith("_")
    ]

    raise RuntimeError(
        f"Could not resolve {label}.\n"
        f"Candidate attributes: {candidates}\n"
        f"Available public attributes: {available}"
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


def coulomb_potential(
    grid_xyz_bohr: np.ndarray,
    atom_xyz_bohr: np.ndarray,
    charges_e: np.ndarray,
    chunk_size: int = 2000,
) -> np.ndarray:
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
            grid_xyz_bohr[start:stop, None, :]
            - atom_xyz_bohr[None, :, :]
        )

        distances = np.linalg.norm(
            displacement,
            axis=2,
        )

        if np.any(distances <= 0.0):
            raise RuntimeError(
                "Grid point coincides with an atomic center"
            )

        potential[start:stop] = np.sum(
            charges_e[None, :] / distances,
            axis=1,
        )

    return potential


def comparison_metrics(
    reference: np.ndarray,
    candidate: np.ndarray,
) -> dict[str, float]:
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

    reference_rms = float(
        np.sqrt(
            np.mean(reference ** 2)
        )
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
            np.mean(
                np.abs(residual)
            )
        ),
        "maximum_absolute_error_au": float(
            np.max(
                np.abs(residual)
            )
        ),
        "relative_RMS_to_reference_RMS": (
            rmse / reference_rms
            if reference_rms > 0.0
            else float("nan")
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
        "residual_mean_au": float(
            np.mean(residual)
        ),
    }


def print_metrics(
    title: str,
    metrics: dict[str, float],
) -> None:
    print(f"\n{title}")

    for name, value in metrics.items():
        print(f"  {name} = {value:.16g}")


print("=" * 100)
print("DAY039 / D039-A4 — REAL-ATOM REFIT FEASIBILITY")
print("=" * 100)


print("\n[1] SOURCE IDENTITY")

require_file(VPOT)
require_file(LATEST_POINTER)

observed_vpot_sha256 = sha256(VPOT)

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

transferability_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1_TRANSFERABILITY.csv"
)

require_file(transferability_csv)

print(
    f"transferability_csv = "
    f"{transferability_csv}"
)
print("source_identity_gate = PASS")


print("\n[2] READ ORCA VPOT OBJECT")

vpot = read_orca_vpot(VPOT)

atom_xyz_bohr = np.asarray(
    get_object_attribute(
        vpot,
        (
            "atom_coordinates_bohr",
            "atom_xyz_bohr",
            "atoms_bohr",
        ),
        "atomic_coordinates",
    ),
    dtype=float,
)

grid_xyz_bohr = np.asarray(
    get_object_attribute(
        vpot,
        (
            "grid_coordinates_bohr",
            "grid_xyz_bohr",
            "esp_coordinates_bohr",
        ),
        "grid_coordinates",
    ),
    dtype=float,
)

quantum_esp_au = np.asarray(
    get_object_attribute(
        vpot,
        (
            "potential_au",
            "grid_potential_au",
            "grid_potentials_au",
            "esp_values_au",
            "potentials_au",
        ),
        "quantum_ESP",
    ),
    dtype=float,
)

quantum_esp_au = quantum_esp_au.reshape(-1)

if atom_xyz_bohr.shape != (52, 3):
    raise RuntimeError(
        "Unexpected VPOT atom-coordinate shape.\n"
        f"Observed: {atom_xyz_bohr.shape}"
    )

if grid_xyz_bohr.shape != (24835, 3):
    raise RuntimeError(
        "Unexpected VPOT grid-coordinate shape.\n"
        f"Observed: {grid_xyz_bohr.shape}"
    )

if quantum_esp_au.shape != (24835,):
    raise RuntimeError(
        "Unexpected VPOT potential shape.\n"
        f"Observed: {quantum_esp_au.shape}"
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
print("VPOT_object_contract_gate = PASS")


print("\n[3] LOAD RESP STAGE 1 ATOMS")

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
        f"Expected 52 rows, observed {len(rows)}"
    )

rows = sorted(
    rows,
    key=lambda row: int(
        row["atom_index_0based"]
    ),
)

indices = [
    int(row["atom_index_0based"])
    for row in rows
]

if indices != list(range(52)):
    raise RuntimeError(
        "Atom order is not exactly 0..51"
    )

full_charges = np.asarray(
    [
        float(
            row["RESP_stage1_charge_e"]
        )
        for row in rows
    ],
    dtype=float,
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

cap_mask = ~real_mask

real_xyz_bohr = atom_xyz_bohr[
    real_mask
]

real_charges_unmodified = full_charges[
    real_mask
]

cap_charges = full_charges[
    cap_mask
]

print(f"full_atom_count = {len(rows)}")
print(
    f"real_atom_count = "
    f"{int(np.sum(real_mask))}"
)
print(
    f"artificial_cap_count = "
    f"{int(np.sum(cap_mask))}"
)
print(
    f"full_RESP_charge_sum_e = "
    f"{np.sum(full_charges):.16g}"
)
print(
    f"real_unmodified_charge_sum_e = "
    f"{np.sum(real_charges_unmodified):.16g}"
)
print(
    f"cap_charge_sum_e = "
    f"{np.sum(cap_charges):.16g}"
)

if int(np.sum(real_mask)) != 37:
    raise RuntimeError(
        "Expected 37 real atoms"
    )

if int(np.sum(cap_mask)) != 15:
    raise RuntimeError(
        "Expected 15 artificial caps"
    )

print("atom_partition_gate = PASS")


print("\n[4] CALCULATE POTENTIALS")

full_RESP52_potential = coulomb_potential(
    grid_xyz_bohr,
    atom_xyz_bohr,
    full_charges,
)

real37_unmodified_potential = coulomb_potential(
    grid_xyz_bohr,
    real_xyz_bohr,
    real_charges_unmodified,
)

cap_only_potential = (
    full_RESP52_potential
    - real37_unmodified_potential
)

print(
    f"quantum_ESP_min_au = "
    f"{np.min(quantum_esp_au):.16g}"
)
print(
    f"quantum_ESP_max_au = "
    f"{np.max(quantum_esp_au):.16g}"
)
print(
    f"quantum_ESP_std_au = "
    f"{np.std(quantum_esp_au):.16g}"
)
print(
    f"full_RESP52_potential_std_au = "
    f"{np.std(full_RESP52_potential):.16g}"
)
print(
    f"real37_unmodified_potential_std_au = "
    f"{np.std(real37_unmodified_potential):.16g}"
)
print(
    f"cap_only_potential_std_au = "
    f"{np.std(cap_only_potential):.16g}"
)


print("\n[5] ELECTROSTATIC COMPARISONS")

metrics_full_vs_quantum = comparison_metrics(
    quantum_esp_au,
    full_RESP52_potential,
)

metrics_real_vs_quantum = comparison_metrics(
    quantum_esp_au,
    real37_unmodified_potential,
)

metrics_real_vs_full = comparison_metrics(
    full_RESP52_potential,
    real37_unmodified_potential,
)

print_metrics(
    "A. Full RESP52 versus quantum VPOT",
    metrics_full_vs_quantum,
)

print_metrics(
    "B. Unmodified real37 versus quantum VPOT",
    metrics_real_vs_quantum,
)

print_metrics(
    "C. Unmodified real37 versus full RESP52",
    metrics_real_vs_full,
)


print("\n[6] ERROR ATTRIBUTION")

full_fit_rmse = (
    metrics_full_vs_quantum["RMSE_au"]
)

real_fit_rmse = (
    metrics_real_vs_quantum["RMSE_au"]
)

cap_removal_rmse = (
    metrics_real_vs_full["RMSE_au"]
)

incremental_rmse = (
    real_fit_rmse - full_fit_rmse
)

rmse_ratio_real_to_full_fit = (
    real_fit_rmse / full_fit_rmse
    if full_fit_rmse > 0.0
    else float("nan")
)

print(
    f"full_RESP52_vs_quantum_RMSE_au = "
    f"{full_fit_rmse:.16g}"
)
print(
    f"real37_vs_quantum_RMSE_au = "
    f"{real_fit_rmse:.16g}"
)
print(
    f"real37_vs_full_RESP52_RMSE_au = "
    f"{cap_removal_rmse:.16g}"
)
print(
    f"incremental_RMSE_after_cap_removal_au = "
    f"{incremental_rmse:.16g}"
)
print(
    f"real37_to_full_fit_RMSE_ratio = "
    f"{rmse_ratio_real_to_full_fit:.16g}"
)


print("\n[7] FEASIBILITY INTERPRETATION")

finite_gate = bool(
    np.all(
        np.isfinite(
            quantum_esp_au
        )
    )
    and np.all(
        np.isfinite(
            full_RESP52_potential
        )
    )
    and np.all(
        np.isfinite(
            real37_unmodified_potential
        )
    )
)

reference_fit_reproduced_gate = (
    abs(
        full_fit_rmse
        - 0.002561369228452248
    )
    <= 5.0e-6
)

real_only_degradation_detected_gate = (
    real_fit_rmse > full_fit_rmse
)

refit_has_improvement_target_gate = (
    cap_removal_rmse > 1.0e-6
)

gates = {
    "source_identity_gate": True,
    "VPOT_object_contract_gate": True,
    "atom_partition_gate": True,
    "finite_value_gate": finite_gate,
    "reference_fit_reproduced_gate": (
        reference_fit_reproduced_gate
    ),
    "real_only_degradation_detected_gate": (
        real_only_degradation_detected_gate
    ),
    "refit_has_improvement_target_gate": (
        refit_has_improvement_target_gate
    ),
    "no_refit_executed_gate": True,
    "no_charge_modified_gate": True,
}

for name, value in gates.items():
    print(
        f"{name} = "
        f"{'PASS' if value else 'FAIL'}"
    )

all_gates_pass = all(
    gates.values()
)


print("\n[8] DECISION")

decision = (
    "D039_A4_REAL_ATOM_REFIT_FEASIBILITY_PASS_"
    "CONSTRAINED_REFIT_DESIGN_REVIEW_AUTHORIZED"
    if all_gates_pass
    else
    "D039_A4_REAL_ATOM_REFIT_FEASIBILITY_"
    "REVIEW_REQUIRED"
)

print(f"decision = {decision}")
print(
    "constrained_refit_design_review_authorized = "
    f"{all_gates_pass}"
)
print("constrained_refit_execution_authorized = False")
print("RESP_stage2_execution_authorized = False")
print("charge_adoption_authorized = False")
print("force_field_adoption_authorized = False")
print("=" * 100)
