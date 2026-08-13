#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]

CAMPAIGN = (
    ROOT
    / "runs"
    / "phase2"
    / "campaign"
)

rows = []

for folder in sorted(CAMPAIGN.iterdir()):

    if not folder.is_dir():
        continue

    T = int(folder.name.replace("K",""))

    report = json.loads(
        (folder/"campaign_report.json").read_text()
    )

    rows.append({

        "Temperature_K": T,

        "Minimization": report["in.minimize"]["pass"],

        "Heating": report["in.heating"]["pass"],

        "Equilibration": report["in.nvt"]["pass"],

        "Production": report["in.production"]["pass"],

        "ProductionTime_s":
            report["in.production"]["elapsed_seconds"]

    })

outfile = CAMPAIGN / "campaign_summary.csv"

with open(outfile,"w",newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=rows[0].keys()
    )

    writer.writeheader()

    writer.writerows(rows)

jsonfile = CAMPAIGN / "campaign_summary.json"

jsonfile.write_text(
    json.dumps(rows,indent=2)
)

print("="*90)
print("DAY046 / PHASE2-A26")
print("CAMPAIGN SUMMARY")
print("="*90)

print(outfile)
print(jsonfile)
