#!/usr/bin/env python3
"""
DAY038 / D038-C1

Generate the first Amber ESP file from the audited ORCA VPOT dataset.

This script does NOT execute RESP.
It only converts the audited VPOT dataset into an Amber-compatible ESP file.
"""

from pathlib import Path
import hashlib

from resp_common import (
    parse_orca_vpot,
    write_amber_esp,
)

AUTHORIZED_SHA = (
    "73df47796c29d0b2a88a03d83efcb723d4f6e583d2e1789e1ea1a43d82c1064d"
)

ROOT = Path(__file__).resolve().parents[2]

VPOT = (
    ROOT /
    "runs/phase1A/day036_qm_f06_upper_v7a_r1_esp_executions/"
    "esp_upper_v7a_r1_20260731T174832Z/"
    "esp_upper_v7a_r1.vpot"
)

OUTPUT = (
    ROOT /
    "runs/phase1A/day038_resp_generation/"
    "candidate_from_orca_vpot.esp"
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            block = f.read(1024 * 1024)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


print("=" * 90)
print("DAY038 / D038-C1")
print("Generate Amber ESP from audited ORCA VPOT")
print("=" * 90)
print()

print("Loading VPOT ...")

dataset = parse_orca_vpot(VPOT)

print(f"atoms        = {dataset['atom_count']}")
print(f"grid_points  = {dataset['grid_point_count']}")
print()

print("Source SHA256")
print(sha256(VPOT))

if sha256(VPOT) != AUTHORIZED_SHA:
    raise RuntimeError("BLOCKED: unauthorized VPOT")

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

print()
print("Writing Amber ESP ...")

write_amber_esp(
    OUTPUT,
    dataset["atom_xyz_bohr"],
    dataset["grid_xyz_bohr"],
    dataset["potential_au"],
)

print()
print("Amber ESP written.")
print()

print("Output file")
print(OUTPUT)

print()
print("Output SHA256")
print(sha256(OUTPUT))

print()
print("Decision")
print("D038_C1_AMBER_ESP_GENERATED")
print("RESP remains blocked.")
