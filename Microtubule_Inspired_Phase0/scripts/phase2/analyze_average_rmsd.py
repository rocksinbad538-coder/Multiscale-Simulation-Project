#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]

sys.path.append(str(ROOT / "scripts" / "phase2"))

from md_analysis.paths import get_paths
from md_analysis.io import read_lammps_dump
from md_analysis.alignment import aligned_rmsd

PATHS = get_paths()

XYZ = PATHS["XYZ"]
OUT = PATHS["OUT"]

OUT.mkdir(parents=True, exist_ok=True)

frames = read_lammps_dump(XYZ)

natoms = frames[0]["natoms"]

# ------------------------------------------------------------
# Average structure
# ------------------------------------------------------------

average_atoms = []

for i in range(natoms):

    atom = dict(frames[0]["atoms"][i])

    atom["x"] = np.mean(
        [f["atoms"][i]["x"] for f in frames]
    )

    atom["y"] = np.mean(
        [f["atoms"][i]["y"] for f in frames]
    )

    atom["z"] = np.mean(
        [f["atoms"][i]["z"] for f in frames]
    )

    average_atoms.append(atom)

# ------------------------------------------------------------
# RMSD
# ------------------------------------------------------------

rows = []

for frame in frames:

    rows.append({

        "timestep": frame["timestep"],

        "AverageRMSD": aligned_rmsd(
            average_atoms,
            frame["atoms"]
        )

    })

csvfile = OUT/"average_rmsd.csv"

with open(csvfile,"w",newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=rows[0].keys()
    )

    writer.writeheader()

    writer.writerows(rows)

report = {

    "frames": len(rows),

    "mean_average_rmsd_A":
        sum(r["AverageRMSD"] for r in rows)/len(rows),

    "maximum_average_rmsd_A":
        max(r["AverageRMSD"] for r in rows),

    "final_average_rmsd_A":
        rows[-1]["AverageRMSD"]

}

outfile = OUT/"AVERAGE_RMSD_REPORT.json"

outfile.write_text(
    json.dumps(report,indent=2)
)

print("="*90)
print("DAY047 / PHASE2-A37")
print("AVERAGE STRUCTURE RMSD")
print("="*90)
print(csvfile)
print(outfile)
