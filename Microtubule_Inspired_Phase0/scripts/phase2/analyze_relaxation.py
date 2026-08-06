#!/usr/bin/env python3

from __future__ import annotations

import math
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]

INITIAL = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_lammps_indexed_model"
    / "PHASE1B_LAMMPS_INDEXED_MODEL.json"
)

FINAL = (
    ROOT
    / "runs"
    / "phase2"
    / "day042_import_test"
    / "minimized.xyz"
)


def read_initial():

    import json

    obj = json.loads(INITIAL.read_text())

    atoms = {}

    for atom in obj["system"]["atoms"]:

        idx = str(atom["lammps_index"])

        atoms[idx] = (

            float(atom["x_A"]),

            float(atom["y_A"]),

            float(atom["z_A"])

        )

    return atoms


def read_xyz():

    atoms = {}

    with open(FINAL) as f:

        lines = f.readlines()

    start = None

    for i, line in enumerate(lines):

        if line.startswith("ITEM: ATOMS"):

            start = i + 1

            break

    if start is None:

        raise RuntimeError("ATOMS section not found.")

    for line in lines[start:]:

        s = line.split()

        idx = s[0]

        atoms[idx] = (

            float(s[2]),

            float(s[3]),

            float(s[4])

        )

    return atoms


initial = read_initial()
final = read_xyz()

displacements = []

for idx in initial:

    xi, yi, zi = initial[idx]
    xf, yf, zf = final[idx]

    d = math.sqrt(
        (xf - xi) ** 2 +
        (yf - yi) ** 2 +
        (zf - zi) ** 2
    )

    displacements.append(d)

rmsd = math.sqrt(
    sum(d * d for d in displacements)
    / len(displacements)
)

print("=" * 90)
print("DAY042 / PHASE2-A8")
print("STRUCTURAL RELAXATION ANALYSIS")
print("=" * 90)
print()

print("Atoms:", len(displacements))
print(f"RMSD               {rmsd:.6f} Å")
print(f"Mean displacement  {sum(displacements)/len(displacements):.6f} Å")
print(f"Maximum shift      {max(displacements):.6f} Å")
