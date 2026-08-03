#!/usr/bin/env python3
"""
Day038 / D038-A4c

Robust global diagnostic of the ORCA .vpot electrostatic potential
against the potential reconstructed from ORCA CHELPG point charges.

Purpose
-------
1. Distinguish unit/scale errors from the intrinsic CHELPG fitting error.
2. Avoid unstable pointwise ratios near zero potential.
3. Quantify residuals, sign agreement, robust ratios, and outliers.
4. Test candidate interpretations of the stored potential unit.

READ ONLY.

Creates no files.
Modifies no files.
Does not execute ORCA or RESP.
"""

from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Iterable, Sequence


BOHR_TO_ANGSTROM = 0.529177210903
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM

HARTREE_TO_EV = 27.211386245988
HARTREE_TO_KCAL_PER_MOL = 627.5094740631

AUTHORIZED_VPOT_SHA256 = (
    "73df47796c29d0b2a88a03d83efcb723"
    "d4f6e583d2e1789e1ea1a43d82c1064d"
)

EXPECTED_ATOMS = 52
EXPECTED_POINTS = 24_835


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def parse_float(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def mean(values: Sequence[float]) -> float:
    if not values:
        raise ValueError("Cannot calculate mean of an empty sequence")

    return sum(values) / len(values)


def population_std(values: Sequence[float]) -> float:
    center = mean(values)

    return math.sqrt(
        sum((value - center) ** 2 for value in values)
        / len(values)
    )


def rmse(reference: Sequence[float], candidate: Sequence[float]) -> float:
    if len(reference) != len(candidate):
        raise ValueError("RMSE input lengths differ")

    return math.sqrt(
        sum(
            (candidate_value - reference_value) ** 2
            for reference_value, candidate_value
            in zip(reference, candidate)
        )
        / len(reference)
    )


def mae(reference: Sequence[float], candidate: Sequence[float]) -> float:
    if len(reference) != len(candidate):
        raise ValueError("MAE input lengths differ")

    return (
        sum(
            abs(candidate_value - reference_value)
            for reference_value, candidate_value
            in zip(reference, candidate)
        )
        / len(reference)
    )


def quantile(sorted_values: Sequence[float], probability: float) -> float:
    if not sorted_values:
        raise ValueError("Cannot calculate quantile of empty data")

    if not 0.0 <= probability <= 1.0:
        raise ValueError("Quantile probability must lie in [0, 1]")

    position = probability * (len(sorted_values) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)

    if lower == upper:
        return sorted_values[lower]

    weight = position - lower

    return (
        sorted_values[lower] * (1.0 - weight)
        + sorted_values[upper] * weight
    )


def regression(
    independent: Sequence[float],
    dependent: Sequence[float],
) -> dict[str, float]:
    if len(independent) != len(dependent):
        raise ValueError("Regression input lengths differ")

    x_mean = mean(independent)
    y_mean = mean(dependent)

    sxx = sum((x - x_mean) ** 2 for x in independent)
    syy = sum((y - y_mean) ** 2 for y in dependent)

    sxy = sum(
        (x - x_mean) * (y - y_mean)
        for x, y in zip(independent, dependent)
    )

    if sxx == 0.0 or syy == 0.0:
        raise ValueError("Regression variance is zero")

    slope = sxy / sxx
    intercept = y_mean - slope * x_mean
    pearson = sxy / math.sqrt(sxx * syy)

    predicted = [
        intercept + slope * value
        for value in independent
    ]

    return {
        "slope": slope,
        "intercept": intercept,
        "pearson_r": pearson,
        "r_squared": pearson * pearson,
        "affine_rmse": rmse(dependent, predicted),
        "affine_mae": mae(dependent, predicted),
    }


def load_vpot(
    path: Path,
) -> tuple[
    list[tuple[float, float, float]],
    list[tuple[float, float, float, float]],
]:
    lines = path.read_text(encoding="utf-8").splitlines()

    if not lines:
        raise ValueError("VPOT file is empty")

    header = lines[0].split()

    if len(header) != 2:
        raise ValueError(f"Invalid VPOT header: {lines[0]!r}")

    atom_count = int(header[0])
    point_count = int(header[1])

    expected_lines = 1 + atom_count + point_count

    if len(lines) != expected_lines:
        raise ValueError(
            f"VPOT line-count mismatch: "
            f"observed={len(lines)}, expected={expected_lines}"
        )

    atom_rows = [
        tuple(parse_float(token) for token in line.split())
        for line in lines[1:1 + atom_count]
    ]

    grid_rows = [
        tuple(parse_float(token) for token in line.split())
        for line in lines[1 + atom_count:]
    ]

    if any(len(row) != 3 for row in atom_rows):
        raise ValueError("One or more VPOT atom rows are malformed")

    if any(len(row) != 4 for row in grid_rows):
        raise ValueError("One or more VPOT grid rows are malformed")

    if atom_count != EXPECTED_ATOMS:
        raise ValueError(
            f"Expected {EXPECTED_ATOMS} atoms, observed {atom_count}"
        )

    if point_count != EXPECTED_POINTS:
        raise ValueError(
            f"Expected {EXPECTED_POINTS} points, observed {point_count}"
        )

    return atom_rows, grid_rows


def load_pc_chelpg(
    path: Path,
) -> list[tuple[float, float, float, float]]:
    lines = path.read_text(encoding="utf-8").splitlines()

    if len(lines) < 2:
        raise ValueError(".pc_chelpg file is incomplete")

    declared_count = int(lines[0].strip())

    rows = []

    for line_number, line in enumerate(lines[2:], start=3):
        if not line.strip():
            continue

        tokens = line.split()

        if len(tokens) != 5 or tokens[0].upper() != "Q":
            raise ValueError(
                f"Malformed .pc_chelpg row at line {line_number}: "
                f"{line!r}"
            )

        charge = parse_float(tokens[1])
        x_angstrom = parse_float(tokens[2])
        y_angstrom = parse_float(tokens[3])
        z_angstrom = parse_float(tokens[4])

        rows.append(
            (
                charge,
                x_angstrom * ANGSTROM_TO_BOHR,
                y_angstrom * ANGSTROM_TO_BOHR,
                z_angstrom * ANGSTROM_TO_BOHR,
            )
        )

    if declared_count != len(rows):
        raise ValueError(
            f".pc_chelpg count mismatch: "
            f"declared={declared_count}, parsed={len(rows)}"
        )

    if len(rows) != EXPECTED_ATOMS:
        raise ValueError(
            f"Expected {EXPECTED_ATOMS} CHELPG rows, "
            f"observed {len(rows)}"
        )

    return rows


def coulomb_potential(
    point_bohr: tuple[float, float, float],
    charges: Iterable[tuple[float, float, float, float]],
) -> tuple[float, float, int]:
    px, py, pz = point_bohr

    potential = 0.0
    minimum_distance = math.inf
    nearest_atom = -1

    for atom_index, (charge, x, y, z) in enumerate(charges):
        dx = px - x
        dy = py - y
        dz = pz - z

        distance = math.sqrt(dx * dx + dy * dy + dz * dz)

        if distance < minimum_distance:
            minimum_distance = distance
            nearest_atom = atom_index

        if distance <= 1.0e-12:
            raise ValueError(
                f"Grid point coincides with atom {atom_index}"
            )

        potential += charge / distance

    return potential, minimum_distance, nearest_atom


phase_root = Path.cwd()

pointer_path = (
    phase_root
    / "runs/phase1A/day036_qm_f06_upper_v7a_r1_esp_executions"
    / "LATEST_ESP_UPPER_V7A_R1_EXECUTION.txt"
)

if not pointer_path.is_file():
    raise SystemExit(
        f"BLOCKED: missing execution pointer: {pointer_path}"
    )

execution_dir = Path(
    pointer_path.read_text(encoding="utf-8").strip()
)

if not execution_dir.is_absolute():
    execution_dir = phase_root / execution_dir

vpot_path = execution_dir / "esp_upper_v7a_r1.vpot"
pc_path = execution_dir / "esp_upper_v7a_r1.pc_chelpg"

for required_path in (vpot_path, pc_path):
    if not required_path.is_file():
        raise SystemExit(
            f"BLOCKED: required file not found: {required_path}"
        )


observed_vpot_sha256 = sha256_file(vpot_path)

if observed_vpot_sha256 != AUTHORIZED_VPOT_SHA256:
    raise SystemExit(
        "BLOCKED: VPOT SHA256 does not match the authorized source"
    )


atom_rows, grid_rows = load_vpot(vpot_path)
charges = load_pc_chelpg(pc_path)

stored_potential = []
calculated_potential = []
minimum_distances = []
nearest_atoms = []

for stored_value, x, y, z in grid_rows:
    calculated_value, minimum_distance, nearest_atom = (
        coulomb_potential((x, y, z), charges)
    )

    stored_potential.append(stored_value)
    calculated_potential.append(calculated_value)
    minimum_distances.append(minimum_distance)
    nearest_atoms.append(nearest_atom)


residuals = [
    calculated - stored
    for stored, calculated
    in zip(stored_potential, calculated_potential)
]

absolute_residuals = sorted(abs(value) for value in residuals)

stored_mean = mean(stored_potential)
calculated_mean = mean(calculated_potential)

stored_std = population_std(stored_potential)
calculated_std = population_std(calculated_potential)

raw_rmse = rmse(stored_potential, calculated_potential)
raw_mae = mae(stored_potential, calculated_potential)
maximum_absolute_error = max(absolute_residuals)

forward_regression = regression(
    stored_potential,
    calculated_potential,
)

inverse_regression = regression(
    calculated_potential,
    stored_potential,
)


print("=" * 100)
print("DAY038 / D038-A4c — ROBUST GLOBAL VPOT DIAGNOSTICS")
print("=" * 100)

print("\n[1] SOURCE AND DATASET")
print(f"vpot_path              = {vpot_path}")
print(f"vpot_sha256            = {observed_vpot_sha256}")
print(f"atom_count             = {len(atom_rows)}")
print(f"grid_point_count       = {len(grid_rows)}")
print(f"CHELPG_charge_count    = {len(charges)}")
print(
    f"CHELPG_charge_sum_e    = "
    f"{sum(row[0] for row in charges):.16g}"
)


print("\n[2] POTENTIAL DISTRIBUTIONS")
print(
    f"stored_min_raw         = {min(stored_potential):.16g}"
)
print(
    f"stored_max_raw         = {max(stored_potential):.16g}"
)
print(f"stored_mean_raw        = {stored_mean:.16g}")
print(f"stored_std_raw         = {stored_std:.16g}")

print(
    f"calculated_min_au      = {min(calculated_potential):.16g}"
)
print(
    f"calculated_max_au      = {max(calculated_potential):.16g}"
)
print(f"calculated_mean_au     = {calculated_mean:.16g}")
print(f"calculated_std_au      = {calculated_std:.16g}")


print("\n[3] RAW AGREEMENT — SAME-NUMERICAL-UNIT HYPOTHESIS")
print(f"RMSE                   = {raw_rmse:.16g}")
print(f"MAE                    = {raw_mae:.16g}")
print(f"maximum_abs_error      = {maximum_absolute_error:.16g}")
print(
    f"RMSE_over_stored_std   = "
    f"{raw_rmse / stored_std:.16g}"
)
print(
    f"MAE_over_stored_std    = "
    f"{raw_mae / stored_std:.16g}"
)
print(
    f"RMSE_over_stored_range = "
    f"{raw_rmse / (max(stored_potential) - min(stored_potential)):.16g}"
)


print("\n[4] LINEAR ASSOCIATION")
print("Model: calculated_CHELPG_au = intercept + slope * stored_VPOT")

for key, value in forward_regression.items():
    print(f"{key:>20} = {value:.16g}")

print("\nInverse model: stored_VPOT = intercept + slope * calculated_CHELPG_au")

for key, value in inverse_regression.items():
    print(f"{key:>20} = {value:.16g}")


print("\n[5] ABSOLUTE-RESIDUAL QUANTILES")

for probability in (
    0.00,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.99,
    0.999,
    1.00,
):
    print(
        f"q{100 * probability:06.2f} = "
        f"{quantile(absolute_residuals, probability):.16g}"
    )


print("\n[6] SIGN AGREEMENT CONDITIONED ON |STORED POTENTIAL|")

for threshold in (
    0.0,
    1.0e-5,
    1.0e-4,
    5.0e-4,
    1.0e-3,
    2.5e-3,
    5.0e-3,
    1.0e-2,
):
    selected = [
        (stored, calculated)
        for stored, calculated
        in zip(stored_potential, calculated_potential)
        if abs(stored) >= threshold
    ]

    same_sign = sum(
        1
        for stored, calculated in selected
        if (
            (stored > 0.0 and calculated > 0.0)
            or (stored < 0.0 and calculated < 0.0)
            or (stored == 0.0 and calculated == 0.0)
        )
    )

    percentage = (
        100.0 * same_sign / len(selected)
        if selected
        else float("nan")
    )

    print(
        f"threshold={threshold:.1e} "
        f"selected={len(selected):>6} "
        f"same_sign={same_sign:>6} "
        f"percentage={percentage:9.4f}"
    )


print("\n[7] ROBUST CALCULATED/STORED RATIOS")

for threshold in (
    1.0e-4,
    5.0e-4,
    1.0e-3,
    2.5e-3,
    5.0e-3,
    1.0e-2,
):
    ratios = sorted(
        calculated / stored
        for stored, calculated
        in zip(stored_potential, calculated_potential)
        if abs(stored) >= threshold
    )

    if not ratios:
        continue

    print(
        f"threshold={threshold:.1e} "
        f"n={len(ratios):>6} "
        f"median={quantile(ratios, 0.50): .9f} "
        f"q05={quantile(ratios, 0.05): .9f} "
        f"q95={quantile(ratios, 0.95): .9f}"
    )


print("\n[8] CANDIDATE STORED-POTENTIAL UNIT TESTS")
print(
    "Each hypothesis converts the stored column to atomic units "
    "before comparison with q/r."
)

unit_hypotheses = (
    (
        "atomic_unit_Eh_per_e",
        1.0,
    ),
    (
        "electronvolt_per_e",
        1.0 / HARTREE_TO_EV,
    ),
    (
        "volt",
        1.0 / HARTREE_TO_EV,
    ),
    (
        "kcal_per_mol_per_e",
        1.0 / HARTREE_TO_KCAL_PER_MOL,
    ),
)

unit_results = []

for name, stored_to_au_factor in unit_hypotheses:
    converted_stored = [
        value * stored_to_au_factor
        for value in stored_potential
    ]

    hypothesis_rmse = rmse(
        converted_stored,
        calculated_potential,
    )

    hypothesis_mae = mae(
        converted_stored,
        calculated_potential,
    )

    unit_results.append(
        (
            hypothesis_rmse,
            name,
            stored_to_au_factor,
            hypothesis_mae,
        )
    )

    print(
        f"{name:<26} "
        f"factor_to_au={stored_to_au_factor:.12g} "
        f"RMSE_au={hypothesis_rmse:.12g} "
        f"MAE_au={hypothesis_mae:.12g}"
    )

best_unit_result = min(unit_results, key=lambda row: row[0])

print(
    "\nbest_candidate_by_raw_RMSE = "
    f"{best_unit_result[1]}"
)


print("\n[9] NEAREST-ATOM DISTANCE")
sorted_distances = sorted(minimum_distances)

for probability in (
    0.00,
    0.01,
    0.05,
    0.50,
    0.95,
    0.99,
    1.00,
):
    print(
        f"q{100 * probability:06.2f}_bohr = "
        f"{quantile(sorted_distances, probability):.16g}"
    )


print("\n[10] TWENTY LARGEST ABSOLUTE RESIDUALS")

ranked_indices = sorted(
    range(len(residuals)),
    key=lambda index: abs(residuals[index]),
    reverse=True,
)

for rank, index in enumerate(ranked_indices[:20], start=1):
    stored_value = stored_potential[index]
    calculated_value = calculated_potential[index]
    residual = residuals[index]
    _, x, y, z = grid_rows[index]

    print(
        f"rank={rank:>2} "
        f"grid_index={index:>5} "
        f"stored={stored_value: .10e} "
        f"calculated={calculated_value: .10e} "
        f"residual={residual: .10e} "
        f"abs_residual={abs(residual):.10e} "
        f"rmin_bohr={minimum_distances[index]:.8f} "
        f"nearest_atom={nearest_atoms[index]:>2} "
        f"xyz_bohr=({x:.7f},{y:.7f},{z:.7f})"
    )


source_gate = (
    observed_vpot_sha256 == AUTHORIZED_VPOT_SHA256
    and len(atom_rows) == EXPECTED_ATOMS
    and len(grid_rows) == EXPECTED_POINTS
    and len(charges) == EXPECTED_ATOMS
)

atomic_unit_is_best = (
    best_unit_result[1] == "atomic_unit_Eh_per_e"
)

same_scale_gate = (
    0.8 <= forward_regression["slope"] <= 1.2
    and abs(forward_regression["intercept"]) <= 5.0e-3
    and forward_regression["pearson_r"] >= 0.95
)

fit_quality_characterized = (
    math.isfinite(raw_rmse)
    and math.isfinite(raw_mae)
    and math.isfinite(forward_regression["pearson_r"])
)


print("\n[11] GATES")
print(
    f"source_integrity_gate       = "
    f"{'PASS' if source_gate else 'FAIL'}"
)
print(
    f"atomic_unit_best_fit_gate   = "
    f"{'PASS' if atomic_unit_is_best else 'FAIL'}"
)
print(
    f"same_scale_association_gate = "
    f"{'PASS' if same_scale_gate else 'FAIL'}"
)
print(
    f"fit_quality_audit_gate      = "
    f"{'PASS' if fit_quality_characterized else 'FAIL'}"
)


all_pass = (
    source_gate
    and atomic_unit_is_best
    and same_scale_gate
    and fit_quality_characterized
)


print("\n[12] D038-A4c DECISION")
print(
    "decision = "
    + (
        "D038_A4C_VPOT_POTENTIAL_ATOMIC_UNIT_"
        "HYPOTHESIS_SUPPORTED_AND_CHELPG_FIT_CHARACTERIZED"
        if all_pass
        else "D038_A4C_REVIEW_REQUIRED"
    )
)

print("\nInterpretation:")
print(
    "- This audit tests the unit and physical correspondence of the "
    "stored quantum ESP."
)
print(
    "- It does not require the fitted CHELPG point charges to reproduce "
    "the quantum ESP exactly."
)
print(
    "- RESP execution remains blocked."
)
print("\nNo files were created or modified.")
print("=" * 100)
