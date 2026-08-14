#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]

CAMPAIGN = ROOT / "runs" / "phase2" / "campaign"

OUT = CAMPAIGN / "convergence_figures"
OUT.mkdir(exist_ok=True)


VARIABLES = [

    ("thermodynamics.csv","Step","PotEng","PotentialEnergy"),

    ("trajectory_summary.csv","timestep","Rg","RadiusGyration"),

    ("aligned_rmsd.csv","timestep","AlignedRMSD","AlignedRMSD"),

]


for filename,xcol,ycol,label in VARIABLES:

    plt.figure(figsize=(8,5))

    for folder in sorted(CAMPAIGN.glob("*K")):

        f = folder/"analysis"/filename

        if not f.exists():
            continue

        df = pd.read_csv(f)

        plt.plot(
            df[xcol],
            df[ycol],
            linewidth=1,
            label=folder.name,
        )

    plt.xlabel("Simulation step")

    plt.ylabel(label)

    plt.title(label+" vs Simulation Time")

    plt.legend()

    plt.tight_layout()

    png = OUT/f"{label}_vs_Time.png"

    pdf = OUT/f"{label}_vs_Time.pdf"

    plt.savefig(png,dpi=300)

    plt.savefig(pdf)

    plt.close()

print("="*90)
print("DAY047 / PHASE2-A35")
print("CONVERGENCE FIGURES")
print("="*90)
print(OUT)
