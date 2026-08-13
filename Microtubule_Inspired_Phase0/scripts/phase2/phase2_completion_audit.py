#!/usr/bin/env python3

from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[2]

CAMPAIGN = ROOT / "runs" / "phase2" / "campaign"

checks = {}

checks["Campaign summary CSV"] = (
    CAMPAIGN / "campaign_scientific_summary.csv"
).exists()

checks["Campaign summary JSON"] = (
    CAMPAIGN / "campaign_scientific_summary.json"
).exists()

checks["Markdown report"] = (
    CAMPAIGN / "CAMPAIGN_SCIENTIFIC_REPORT.md"
).exists()

figures = CAMPAIGN / "figures"

checks["Figures directory"] = figures.exists()

expected = [
    "Rg_vs_Temperature.png",
    "Mean_RMSD_vs_Temperature.png",
    "Mean_RMSF_vs_Temperature.png",
    "PotentialEnergy_vs_Temperature.png",
    "ShapeAnisotropy_vs_Temperature.png",
    "AlignedRMSD_vs_Temperature.png",
]

checks["All figures"] = all(
    (figures / f).exists()
    for f in expected
)

temps = ["150K","200K","250K","300K","350K"]

checks["All temperature folders"] = all(
    (CAMPAIGN / t).exists()
    for t in temps
)

checks["All campaign reports"] = all(
    (CAMPAIGN / t / "campaign_report.json").exists()
    for t in temps
)

audit = {
    "completed": checks,
    "completed_items": sum(checks.values()),
    "total_items": len(checks),
    "completion_percent":
        100.0 * sum(checks.values()) / len(checks)
}

outfile = CAMPAIGN / "PHASE2_COMPLETION_AUDIT.json"

outfile.write_text(
    json.dumps(audit, indent=2)
)

print("="*90)
print("DAY046 / PHASE2-A30")
print("PHASE2 COMPLETION AUDIT")
print("="*90)

print(json.dumps(audit, indent=2))
