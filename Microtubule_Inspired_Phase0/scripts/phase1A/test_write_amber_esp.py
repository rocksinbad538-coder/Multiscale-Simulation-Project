#!/usr/bin/env python3

from pathlib import Path
import hashlib

from resp_common import (
    read_orca_vpot,
    write_amber_esp,
)

ROOT = Path(__file__).resolve().parents[2]

VPOT = (
    ROOT
    / "runs/phase1A/day036_qm_f06_upper_v7a_r1_esp_executions"
    / "esp_upper_v7a_r1_20260731T174832Z"
    / "esp_upper_v7a_r1.vpot"
)

OUTPUT = (
    ROOT
    / "runs/phase1A/day038_resp_generation"
    / "candidate_from_orca_vpot.esp"
)

AUTHORIZED_SHA = (
    "73df47796c29d0b2a88a03d83efcb723d4f6e583d2e1789e1ea1a43d82c1064d"
)


def sha256(path: Path):

    h = hashlib.sha256()

    with path.open("rb") as f:

        while True:

            block = f.read(1024 * 1024)

            if not block:
                break

            h.update(block)

    return h.hexdigest()


print("=" * 100)
print("DAY038 / D038-C1A")
print("RESP writer integration test")
print("=" * 100)

print()

print("Loading VPOT...")

vpot = read_orca_vpot(VPOT)

print(f"atoms       = {vpot.atom_count}")
print(f"grid points = {vpot.grid_point_count}")

print()

print("Checking SHA256...")

observed = sha256(VPOT)

print(observed)

if observed != AUTHORIZED_SHA:
    raise RuntimeError("Unauthorized VPOT")

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

print()

print("Writing Amber ESP...")

write_amber_esp(
    OUTPUT,
    vpot.atom_coordinates_bohr,
    vpot.grid_coordinates_bohr,
    vpot.grid_potential_au,
)

print()

print("DONE")

print("output =", OUTPUT)

print("sha256 =", sha256(OUTPUT))
