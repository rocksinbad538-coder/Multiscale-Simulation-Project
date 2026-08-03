#!/usr/bin/env python3

from pathlib import Path
import numpy as np

from resp_common import (
    read_orca_vpot,
    read_amber_esp,
)

ROOT = Path(__file__).resolve().parents[2]

VPOT = (
    ROOT
    / "runs/phase1A/day036_qm_f06_upper_v7a_r1_esp_executions"
    / "esp_upper_v7a_r1_20260731T174832Z"
    / "esp_upper_v7a_r1.vpot"
)

ESP = (
    ROOT
    / "runs/phase1A/day038_resp_generation"
    / "candidate_from_orca_vpot.esp"
)

print("=" * 100)
print("DAY038 / D038-C2")
print("AMBER ESP ROUNDTRIP VALIDATION")
print("=" * 100)
print()

vpot = read_orca_vpot(VPOT)
esp = read_amber_esp(ESP)

print("[1] COUNTS")
print(f"VPOT atoms  = {vpot.atom_count}")
print(f"ESP atoms   = {esp['natoms']}")
print(f"VPOT points = {vpot.grid_point_count}")
print(f"ESP points  = {esp['npoints']}")
print()

atoms_diff = np.abs(
    np.asarray(vpot.atom_coordinates_bohr)
    - esp["atom_xyz_bohr"]
)

grid_diff = np.abs(
    np.asarray(vpot.grid_coordinates_bohr)
    - esp["esp_xyz_bohr"]
)

pot_diff = np.abs(
    np.asarray(vpot.grid_potential_au)
    - esp["esp_values_au"]
)

print("[2] MAXIMUM DIFFERENCES")

print(
    "atom coordinates :",
    atoms_diff.max(),
)

print(
    "grid coordinates :",
    grid_diff.max(),
)

print(
    "ESP values       :",
    pot_diff.max(),
)

print()

print("[3] RMS DIFFERENCES")

print(
    "atom coordinates :",
    np.sqrt(np.mean(atoms_diff ** 2)),
)

print(
    "grid coordinates :",
    np.sqrt(np.mean(grid_diff ** 2)),
)

print(
    "ESP values       :",
    np.sqrt(np.mean(pot_diff ** 2)),
)

print()

if (
    atoms_diff.max() < 1e-12
    and grid_diff.max() < 1e-12
    and pot_diff.max() < 1e-12
):
    decision = (
        "D038_C2_ROUNDTRIP_PASS"
    )
else:
    decision = (
        "D038_C2_ROUNDTRIP_NUMERICAL_DIFFERENCE"
    )

print("[4] DECISION")
print(decision)
