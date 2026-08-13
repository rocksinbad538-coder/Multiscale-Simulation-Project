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

from md_analysis.io import read_lammps_dump
from md_analysis.alignment import aligned_rmsd

print("="*90)
print("DAY046 / PHASE2-A22")
print("ALIGNED RMSD")
print("="*90)
print()

frames = read_lammps_dump(TRAJECTORY)

reference = frames[0]

rows = []

for frame in frames:

    rows.append({

        "timestep": frame["timestep"],

        "AlignedRMSD": aligned_rmsd(
            reference["atoms"],
            frame["atoms"]
        )

    })

csvfile = OUT / "aligned_rmsd.csv"

with open(csvfile,"w",newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=rows[0].keys()
    )

    writer.writeheader()

    writer.writerows(rows)

report = {

    "frames": len(rows),

    "final_aligned_rmsd_A":
        rows[-1]["AlignedRMSD"],

    "maximum_aligned_rmsd_A":
        max(
            x["AlignedRMSD"]
            for x in rows
        ),

    "mean_aligned_rmsd_A":
        sum(
            x["AlignedRMSD"]
            for x in rows
        )/len(rows)

}

outfile = OUT / "ALIGNED_RMSD_REPORT.json"

outfile.write_text(
    json.dumps(report,indent=2)
)

print(csvfile)
print(outfile)
