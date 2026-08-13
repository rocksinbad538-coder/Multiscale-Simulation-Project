#!/usr/bin/env python3

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]

CAMPAIGN = ROOT / "runs" / "phase2" / "campaign"

csvfile = CAMPAIGN / "campaign_scientific_summary.csv"

df = pd.read_csv(csvfile)

OUT = CAMPAIGN / "figures"
OUT.mkdir(exist_ok=True)

plots = [

    ("Final_Rg_A",
     "Final Radius of Gyration (Å)",
     "Rg_vs_Temperature"),

    ("Mean_RMSD_A",
     "Mean RMSD (Å)",
     "Mean_RMSD_vs_Temperature"),

    ("Mean_RMSF_A",
     "Mean RMSF (Å)",
     "Mean_RMSF_vs_Temperature"),

    ("Final_PE",
     "Final Potential Energy",
     "PotentialEnergy_vs_Temperature"),

    ("Shape_kappa2",
     "Mean Relative Shape Anisotropy",
     "ShapeAnisotropy_vs_Temperature"),

    ("Aligned_RMSD_A",
     "Mean Aligned RMSD (Å)",
     "AlignedRMSD_vs_Temperature"),

]

for column, ylabel, name in plots:

    plt.figure(figsize=(6,4))

    plt.plot(
        df["Temperature_K"],
        df[column],
        marker="o",
        linewidth=2
    )

    plt.xlabel("Temperature (K)")
    plt.ylabel(ylabel)

    plt.grid(True)

    plt.tight_layout()

    plt.savefig(
        OUT / f"{name}.png",
        dpi=300
    )

    plt.savefig(
        OUT / f"{name}.pdf"
    )

    plt.close()

print("="*90)
print("DAY046 / PHASE2-A28")
print("CAMPAIGN FIGURES")
print("="*90)

for f in sorted(OUT.iterdir()):
    print(f)
