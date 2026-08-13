#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import pathlib
import sys

from md_analysis.paths import get_paths



ROOT = pathlib.Path(__file__).resolve().parents[2]

PATHS = get_paths()

TRAJECTORY = PATHS["XYZ"]

OUT = PATHS["OUT"]

sys.path.append(str(ROOT / "scripts" / "phase2"))

from md_analysis.inertia import (
    principal_moments,
    asphericity,
    acylindricity,
    relative_shape_anisotropy,
)

from md_analysis.io import read_lammps_dump


OUT = PATHS["OUT"]

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

print("="*90)
print("DAY046 / PHASE2-A21")
print("SHAPE ANALYSIS")
print("="*90)
print()

frames = read_lammps_dump(TRAJECTORY)

rows = []

for frame in frames:

    eigvals, eigvecs = principal_moments(frame["atoms"])

    rows.append({

        "timestep": frame["timestep"],

        "I1": float(eigvals[0]),
        "I2": float(eigvals[1]),
        "I3": float(eigvals[2]),

        "Asphericity": float(asphericity(eigvals)),
        "Acylindricity": float(acylindricity(eigvals)),
        "RelativeShapeAnisotropy":
            float(relative_shape_anisotropy(eigvals)),
    })

csvfile = OUT / "shape_analysis.csv"

with open(csvfile,"w",newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=rows[0].keys()
    )

    writer.writeheader()
    writer.writerows(rows)

report = {

    "frames": len(rows),

    "final_asphericity":
        rows[-1]["Asphericity"],

    "final_acylindricity":
        rows[-1]["Acylindricity"],

    "final_relative_shape_anisotropy":
        rows[-1]["RelativeShapeAnisotropy"],

    "maximum_relative_shape_anisotropy":
        max(
            x["RelativeShapeAnisotropy"]
            for x in rows
        ),

    "mean_relative_shape_anisotropy":
        sum(
            x["RelativeShapeAnisotropy"]
            for x in rows
        )/len(rows)

}

outfile = OUT / "SHAPE_REPORT.json"

outfile.write_text(
    json.dumps(report,indent=2)
)

print(csvfile)
print(outfile)
