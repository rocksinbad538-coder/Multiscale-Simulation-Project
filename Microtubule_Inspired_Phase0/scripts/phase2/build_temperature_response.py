#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]

META = ROOT/"runs"/"phase2"/"campaign"/"campaign_meta_summary.csv"

OUT = ROOT/"runs"/"phase2"/"campaign"/"temperature_response"

OUT.mkdir(exist_ok=True)

df = pd.read_csv(META)

plots = [

    ("final_Rg_A",
     "Final Radius of Gyration (Å)",
     "Figure1_FinalRg_vs_Temperature.png"),

    ("mean_RMSF_A",
     "Mean RMSF (Å)",
     "Figure2_MeanRMSF_vs_Temperature.png"),

    ("maximum_RMSF_A",
     "Maximum RMSF (Å)",
     "Figure3_MaxRMSF_vs_Temperature.png"),

    ("mean_aligned_rmsd_A",
     "Mean Aligned RMSD (Å)",
     "Figure4_MeanAlignedRMSD_vs_Temperature.png"),

    ("mean_relative_shape_anisotropy",
     "Mean Relative Shape Anisotropy",
     "Figure5_Anisotropy_vs_Temperature.png"),

    ("maximum_atomic_displacement_A",
     "Maximum Atomic Displacement (Å)",
     "Figure6_MaxDisplacement_vs_Temperature.png"),
]

for column,title,filename in plots:

    plt.figure(figsize=(6,4))

    plt.plot(
        df["Temperature_K"],
        df[column],
        "-o",
        linewidth=2,
        markersize=6
    )

    plt.grid(alpha=0.3)

    plt.xlabel("Temperature (K)")
    plt.ylabel(title)
    plt.title(title)

    plt.tight_layout()

    plt.savefig(
        OUT/filename,
        dpi=300
    )

    plt.close()

print("="*90)
print("DAY049 / PHASE2-B09")
print("TEMPERATURE RESPONSE")
print("="*90)
print(OUT)
