#!/usr/bin/env python3

from pathlib import Path
import json
import csv

ROOT = Path(__file__).resolve().parents[2]

CAMPAIGN = ROOT / "runs" / "phase2" / "campaign"

rows = []

for folder in sorted(CAMPAIGN.iterdir()):

    if not folder.is_dir():
        continue

    report = folder / "analysis" / "MD_SCIENTIFIC_REPORT.json"

    if not report.exists():
        continue

    r = json.loads(report.read_text())

    rows.append({

        "Temperature_K": int(folder.name[:-1]),

        "Final_Rg_A":
            r["trajectory"]["final_Rg_A"],

        "Mean_RMSD_A":
            r["trajectory"]["mean_RMSD_A"],

        "Final_PE":
            r["thermodynamics"]["final_PE"],

        "Mean_Temperature_K":
            r["thermodynamics"]["mean_temperature"],

        "Mean_RMSF_A":
            r["rmsf"]["mean_RMSF_A"],

        "Shape_kappa2":
            r["shape"]["mean_relative_shape_anisotropy"],

        "Aligned_RMSD_A":
            r["aligned_rmsd"]["mean_aligned_rmsd_A"]

    })

csvfile = CAMPAIGN / "campaign_scientific_summary.csv"

with open(csvfile,"w",newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=rows[0].keys()
    )

    writer.writeheader()

    writer.writerows(rows)

jsonfile = CAMPAIGN / "campaign_scientific_summary.json"

jsonfile.write_text(
    json.dumps(rows,indent=2)
)

print("="*90)
print("DAY046 / PHASE2-A27")
print("CAMPAIGN SCIENTIFIC SUMMARY")
print("="*90)

print(csvfile)
print(jsonfile)
