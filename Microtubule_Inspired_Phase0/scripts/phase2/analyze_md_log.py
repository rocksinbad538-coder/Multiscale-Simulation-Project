#!/usr/bin/env python3

from __future__ import annotations

import csv
from md_analysis.paths import get_paths
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]

PATHS = get_paths()

OUT = PATHS["OUT"]

LOG = PATHS["LOG"]

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

rows = []

reading = False

with LOG.open() as f:

    for line in f:

        if "Step" in line and "Temp" in line:
            reading = True
            continue

        if not reading:
            continue

        if line.startswith("Loop time"):
            break

        s = line.split()

        if len(s) != 11:
            continue

        rows.append({

            "Step": int(s[0]),
            "Temp": float(s[1]),
            "PotEng": float(s[2]),
            "KinEng": float(s[3]),
            "TotEng": float(s[4]),
            "Ebond": float(s[5]),
            "Eangle": float(s[6]),
            "Edihed": float(s[7]),
            "Eimproper": float(s[8]),
            "Evdwl": float(s[9]),
            "Press": float(s[10])

        })

if len(rows) == 0:
    raise RuntimeError("No thermodynamic records were parsed.")

csvfile = OUT / "thermodynamics.csv"

with csvfile.open("w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=rows[0].keys()
    )

    writer.writeheader()
    writer.writerows(rows)

report = {

    "frames": len(rows),

    "initial_step": rows[0]["Step"],
    "final_step": rows[-1]["Step"],

    "initial_PE": rows[0]["PotEng"],
    "final_PE": rows[-1]["PotEng"],

    "minimum_PE": min(r["PotEng"] for r in rows),
    "maximum_PE": max(r["PotEng"] for r in rows),

    "mean_temperature":

        sum(r["Temp"] for r in rows) / len(rows),

    "minimum_temperature":

        min(r["Temp"] for r in rows),

    "maximum_temperature":

        max(r["Temp"] for r in rows),

    "mean_pressure":

        sum(r["Press"] for r in rows) / len(rows),

    "minimum_pressure":

        min(r["Press"] for r in rows),

    "maximum_pressure":

        max(r["Press"] for r in rows)

}

outfile = OUT / "THERMODYNAMIC_REPORT.json"

outfile.write_text(
    json.dumps(
        report,
        indent=2
    )
)

print("="*90)
print("DAY045 / PHASE2-A18")
print("THERMODYNAMIC ANALYSIS")
print("="*90)
print()

print("Frames :", report["frames"])
print("Initial step :", report["initial_step"])
print("Final step   :", report["final_step"])
print()

print(csvfile)
print(outfile)

