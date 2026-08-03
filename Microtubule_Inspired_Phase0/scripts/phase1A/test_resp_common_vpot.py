#!/usr/bin/env python3
"""
Day038 / D038-B1

Deterministic integration test for the reusable ORCA VPOT and
CHELPG point-charge functions in resp_common.py.

The test uses the adopted Day035 geometry and the authorized Day036
ORCA ESP artifacts.

READ ONLY with respect to scientific source data.

Creates no scientific output files.
Does not execute ORCA or RESP.
"""

from __future__ import annotations

import math
from pathlib import Path

from resp_common import (
    BOHR_TO_ANGSTROM,
    coordinate_difference_metrics,
    convert_bohr_to_angstrom,
    read_orca_pc_chelpg,
    read_orca_vpot,
    reconstruct_coulomb_potential,
    sha256,
    validate_orca_pc_chelpg,
    validate_orca_vpot,
)


AUTHORIZED_VPOT_SHA256 = (
    "73df47796c29d0b2a88a03d83efcb723"
    "d4f6e583d2e1789e1ea1a43d82c1064d"
)

AUTHORIZED_PC_CHELPG_SHA256 = (
    "a5b634abf509a79d7b377223aae321d68"
    "b25b39b0eff3e69bbba8aca87836981"
)

AUTHORIZED_XYZ_SHA256 = (
    "59cfd417753fbf6e5e4adf78a91761c2"
    "927824c50e24c7410010917d387574b2"
)

EXPECTED_ATOMS = 52
EXPECTED_GRID_POINTS = 24_835


def parse_adopted_xyz(
    path: Path,
) -> tuple[
    tuple[str, ...],
    tuple[tuple[float, float, float], ...],
]:
    lines = path.read_text(
        encoding="utf-8"
    ).splitlines()

    atom_count = int(lines[0].strip())

    body = [
        line.split()
        for line in lines[2:]
        if line.strip()
    ]

    if len(body) != atom_count:
        raise RuntimeError(
            "Adopted XYZ count mismatch"
        )

    elements = tuple(row[0] for row in body)

    coordinates = tuple(
        tuple(float(value) for value in row[1:4])
        for row in body
    )

    return elements, coordinates


phase_root = Path.cwd()

pointer_path = (
    phase_root
    / "runs/phase1A/day036_qm_f06_upper_v7a_r1_esp_executions"
    / "LATEST_ESP_UPPER_V7A_R1_EXECUTION.txt"
)

if not pointer_path.is_file():
    raise SystemExit(
        f"BLOCKED: missing pointer: {pointer_path}"
    )

execution_dir = Path(
    pointer_path.read_text(
        encoding="utf-8"
    ).strip()
)

if not execution_dir.is_absolute():
    execution_dir = phase_root / execution_dir

vpot_path = (
    execution_dir
    / "esp_upper_v7a_r1.vpot"
)

pc_path = (
    execution_dir
    / "esp_upper_v7a_r1.pc_chelpg"
)

xyz_path = (
    phase_root
    / "runs/phase1A/"
      "day035_qm_f06_upper_v7a_r1_coordinate_adoption/"
      "QM_F06_UPPER_V7A_ADOPTED_FINAL.xyz"
)

for required_path in (
    vpot_path,
    pc_path,
    xyz_path,
):
    if not required_path.is_file():
        raise SystemExit(
            f"BLOCKED: missing required file: "
            f"{required_path}"
        )


print("=" * 100)
print("DAY038 / D038-B1 — RESP_COMMON ORCA VPOT INTEGRATION TEST")
print("=" * 100)


print("\n[1] SOURCE IDENTITY")
print(f"vpot_sha256 = {sha256(vpot_path)}")
print(f"pc_sha256   = {sha256(pc_path)}")
print(f"xyz_sha256  = {sha256(xyz_path)}")

source_hash_gate = (
    sha256(vpot_path) == AUTHORIZED_VPOT_SHA256
    and sha256(pc_path)
    == AUTHORIZED_PC_CHELPG_SHA256
    and sha256(xyz_path) == AUTHORIZED_XYZ_SHA256
)

print(
    f"source_hash_gate = "
    f"{'PASS' if source_hash_gate else 'FAIL'}"
)


vpot = read_orca_vpot(
    vpot_path,
    expected_sha256=AUTHORIZED_VPOT_SHA256,
)

pc = read_orca_pc_chelpg(
    pc_path,
    expected_sha256=AUTHORIZED_PC_CHELPG_SHA256,
)


print("\n[2] PARSER OUTPUT")
print(f"vpot_atom_count       = {vpot.atom_count}")
print(
    f"vpot_grid_point_count = "
    f"{vpot.grid_point_count}"
)
print(f"pc_atom_count         = {pc.atom_count}")
print(f"pc_comment            = {pc.comment!r}")


vpot_validation = validate_orca_vpot(
    vpot,
    expected_atom_count=EXPECTED_ATOMS,
    expected_grid_point_count=EXPECTED_GRID_POINTS,
)

pc_validation = validate_orca_pc_chelpg(
    pc,
    expected_atom_count=EXPECTED_ATOMS,
    expected_total_charge_e=0.0,
    total_charge_tolerance_e=1.0e-5,
)


print("\n[3] VPOT VALIDATION SUMMARY")
for key, value in vpot_validation.items():
    print(f"{key} = {value}")


print("\n[4] PC_CHELPG VALIDATION SUMMARY")
for key, value in pc_validation.items():
    print(f"{key} = {value}")


elements, adopted_coordinates_A = (
    parse_adopted_xyz(xyz_path)
)

vpot_atomic_coordinates_A = tuple(
    convert_bohr_to_angstrom(row)
    for row in vpot.atom_coordinates_bohr
)


vpot_vs_xyz = coordinate_difference_metrics(
    adopted_coordinates_A,
    vpot_atomic_coordinates_A,
)

pc_vs_xyz = coordinate_difference_metrics(
    adopted_coordinates_A,
    pc.coordinates_angstrom,
)

pc_vs_vpot = coordinate_difference_metrics(
    vpot_atomic_coordinates_A,
    pc.coordinates_angstrom,
)


print("\n[5] GEOMETRY CORRESPONDENCE")
print(
    f"composition = "
    f"B:{elements.count('B')} "
    f"N:{elements.count('N')} "
    f"H:{elements.count('H')}"
)

print("\nVPOT converted to angstrom versus adopted XYZ")
for key, value in vpot_vs_xyz.items():
    print(f"{key} = {value:.16g}")

print("\nPC_CHELPG versus adopted XYZ")
for key, value in pc_vs_xyz.items():
    print(f"{key} = {value:.16g}")

print("\nPC_CHELPG versus VPOT converted to angstrom")
for key, value in pc_vs_vpot.items():
    print(f"{key} = {value:.16g}")


geometry_gate = (
    vpot_vs_xyz["maximum_absolute_component"]
    <= 1.0e-6
    and pc_vs_xyz["maximum_absolute_component"]
    <= 1.0e-6
    and pc_vs_vpot["maximum_absolute_component"]
    <= 1.0e-6
)

print(
    f"\ngeometry_correspondence_gate = "
    f"{'PASS' if geometry_gate else 'FAIL'}"
)


print("\n[6] FIRST TEN COULOMB RECONSTRUCTIONS")

pilot_coordinates = (
    vpot.grid_coordinates_bohr[:10]
)

(
    calculated_potential,
    nearest_distances,
    nearest_indices,
) = reconstruct_coulomb_potential(
    pilot_coordinates,
    pc.charges_e,
    pc.coordinates_bohr,
)

pilot_abs_errors = []

for index, (
    stored,
    calculated,
    distance,
    nearest_atom,
) in enumerate(
    zip(
        vpot.grid_potential_au[:10],
        calculated_potential,
        nearest_distances,
        nearest_indices,
    )
):
    difference = calculated - stored
    pilot_abs_errors.append(abs(difference))

    print(
        f"index={index:>2} "
        f"stored_au={stored: .10e} "
        f"calculated_au={calculated: .10e} "
        f"difference_au={difference: .10e} "
        f"rmin_bohr={distance:.8f} "
        f"nearest_atom={nearest_atom}"
    )


pilot_gate = (
    max(pilot_abs_errors) <= 2.0e-3
    and all(
        math.isfinite(value)
        for value in calculated_potential
    )
    and min(nearest_distances) > 2.0
)

print(
    f"\npilot_coulomb_gate = "
    f"{'PASS' if pilot_gate else 'FAIL'}"
)


unit_constant_gate = (
    abs(
        BOHR_TO_ANGSTROM
        - 0.529177210903
    )
    <= 1.0e-15
)

count_gate = (
    vpot.atom_count == EXPECTED_ATOMS
    and vpot.grid_point_count
    == EXPECTED_GRID_POINTS
    and pc.atom_count == EXPECTED_ATOMS
)

finite_gate = (
    vpot_validation["finite_value_gate"]
    and pc_validation["finite_value_gate"]
)


print("\n[7] FINAL GATES")
print(
    f"unit_constant_gate          = "
    f"{'PASS' if unit_constant_gate else 'FAIL'}"
)
print(
    f"count_gate                  = "
    f"{'PASS' if count_gate else 'FAIL'}"
)
print(
    f"finite_value_gate           = "
    f"{'PASS' if finite_gate else 'FAIL'}"
)
print(
    f"geometry_correspondence_gate = "
    f"{'PASS' if geometry_gate else 'FAIL'}"
)
print(
    f"pilot_coulomb_gate          = "
    f"{'PASS' if pilot_gate else 'FAIL'}"
)


all_pass = (
    source_hash_gate
    and unit_constant_gate
    and count_gate
    and finite_gate
    and geometry_gate
    and pilot_gate
)


print("\n[8] D038-B1 DECISION")
print(
    "decision = "
    + (
        "D038_B1_RESP_COMMON_ORCA_VPOT_AND_"
        "CHELPG_PARSERS_IMPLEMENTED_AND_VALIDATED"
        if all_pass
        else "D038_B1_BLOCKED_REVIEW_REQUIRED"
    )
)

print("\nRESP execution remains blocked.")
print("=" * 100)
