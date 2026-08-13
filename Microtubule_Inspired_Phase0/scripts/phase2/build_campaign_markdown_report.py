#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

CAMPAIGN = ROOT / "runs" / "phase2" / "campaign"

df = pd.read_csv(CAMPAIGN / "campaign_scientific_summary.csv")

report = CAMPAIGN / "CAMPAIGN_SCIENTIFIC_REPORT.md"

lines = []

lines.append("# Phase 2 Multi-temperature Molecular Dynamics Campaign\n")
lines.append("")
lines.append("## Simulation summary\n")
lines.append("")
lines.append(f"- Number of simulations: {len(df)}")
lines.append("- Temperature range: 150–350 K")
lines.append("- Production length: 1 ns")
lines.append("- Force field: Phase1B")
lines.append("- Ensemble: NVT")
lines.append("")

lines.append("## Comparative results\n")
lines.append("")

# -------------------------------------------------------
# Markdown table without external dependencies
# -------------------------------------------------------

cols = list(df.columns)

lines.append("| " + " | ".join(cols) + " |")
lines.append("|" + "|".join(["---"] * len(cols)) + "|")

for _, row in df.iterrows():

    values = []

    for c in cols:

        v = row[c]

        if isinstance(v, float):
            values.append(f"{v:.6f}")
        else:
            values.append(str(v))

    lines.append("| " + " | ".join(values) + " |")

lines.append("")
lines.append("## Preliminary observations\n")
lines.append("")
lines.append("- Temperature control remained close to the target value for all simulations.")
lines.append("- Radius of gyration shows moderate temperature dependence.")
lines.append("- RMSD varies smoothly across temperatures.")
lines.append("- RMSF remains comparatively stable.")
lines.append("- Relative shape anisotropy remains nearly constant.")
lines.append("")
lines.append("## Figures\n")
lines.append("")
lines.append("- figures/Rg_vs_Temperature.png")
lines.append("- figures/Mean_RMSD_vs_Temperature.png")
lines.append("- figures/Mean_RMSF_vs_Temperature.png")
lines.append("- figures/PotentialEnergy_vs_Temperature.png")
lines.append("- figures/ShapeAnisotropy_vs_Temperature.png")
lines.append("- figures/AlignedRMSD_vs_Temperature.png")

report.write_text("\n".join(lines))

print("="*90)
print("DAY046 / PHASE2-A29")
print("CAMPAIGN MARKDOWN REPORT")
print("="*90)
print(report)
