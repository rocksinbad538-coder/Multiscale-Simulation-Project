#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import pathlib
import sys

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

rows = []

rows.append({
    "timestep": frames[0]["timestep"],
    "IncrementalRMSD": 0.0,
})

for previous, current in zip(frames[:-1], frames[1:]):

    rows.append({

        "timestep": current["timestep"],

        "IncrementalRMSD": aligned_rmsd(
            previous["atoms"],
            current["atoms"],
        )

    })

csvfile = OUT / "incremental_rmsd.csv"

with open(csvfile, "w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=rows[0].keys()
    )

    writer.writeheader()
    writer.writerows(rows)

report = {

    "frames": len(rows),

    "mean_incremental_rmsd_A":
        sum(r["IncrementalRMSD"] for r in rows)/len(rows),

    "maximum_incremental_rmsd_A":
        max(r["IncrementalRMSD"] for r in rows),

    "final_incremental_rmsd_A":
        rows[-1]["IncrementalRMSD"]

}

outfile = OUT/"INCREMENTAL_RMSD_REPORT.json"

outfile.write_text(
    json.dumps(report, indent=2)
)

print("="*90)
print("DAY047 / PHASE2-A36")
print("INCREMENTAL RMSD")
print("="*90)
print(csvfile)
print(outfile)
