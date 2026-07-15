#!/usr/bin/env python3

"""
Day026
Primary Source Reclassification

Initializes the Day026 evidence matrix using the Day025
environment classification.

No scientific decisions are taken here.
The script only prepares the workspace for the literature audit.
"""

from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[2]

INPUT = ROOT / (
    "runs/phase1A/day025_force_field_route_comparison/"
    "39_r2_force_field_route_coverage_and_risk_matrix/"
    "r2_environment_route_coverage_and_risk_matrix.csv"
)

OUTPUT = ROOT / (
    "runs/phase1A/day026_primary_source_reclassification/"
    "r2_primary_source_reclassification_matrix.csv"
)

rows = []

with open(INPUT, newline="") as f:

    reader = csv.DictReader(f)

    for r in reader:

        rows.append({

            "Environment_ID":
                r.get("environment_id",""),

            "Environment_Name":
                r.get("environment_name",""),

            "Chemical_Category":
                r.get("chemical_category",""),

            "Day025_Status":
                r.get("coverage_status",""),

            "Rajan2018":
                "",

            "Bamane2023":
                "",

            "Ghorai2025":
                "",

            "Literature_Coverage":
                "",

            "Recommended_Action":
                "",

            "Scientific_Justification":
                ""

        })

with open(OUTPUT,"w",newline="") as f:

    writer = csv.DictWriter(

        f,

        fieldnames=[

            "Environment_ID",
            "Environment_Name",
            "Chemical_Category",
            "Day025_Status",
            "Rajan2018",
            "Bamane2023",
            "Ghorai2025",
            "Literature_Coverage",
            "Recommended_Action",
            "Scientific_Justification"

        ]

    )

    writer.writeheader()

    writer.writerows(rows)

print()

print("Day026 matrix initialized.")

print("Environments:",len(rows))

print("Output:")

print(OUTPUT)
