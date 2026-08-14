#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from md_analysis.statistics import running_average, running_std

ROOT = Path(__file__).resolve().parents[2]

CAMPAIGN = ROOT / "runs" / "phase2" / "campaign"

OUT = CAMPAIGN / "running_statistics"

OUT.mkdir(exist_ok=True)

summary = []

for folder in sorted(CAMPAIGN.glob("*K")):

    T = folder.name

    traj = pd.read_csv(
        folder/"analysis"/"trajectory_summary.csv"
    )

    thermo = pd.read_csv(
        folder/"analysis"/"thermodynamics.csv"
    )

    aligned = pd.read_csv(
        folder/"analysis"/"aligned_rmsd.csv"
    )

    n = min(
        len(traj),
        len(thermo),
        len(aligned),
    )

    df = pd.DataFrame({

        "Rg": traj["Rg"][:n],

        "RMSD": traj["RMSD"][:n],

        "AlignedRMSD": aligned["AlignedRMSD"][:n],

        "PotEng": thermo["PotEng"][:n],

        "Temp": thermo["Temp"][:n]

    })

    for col in df.columns:

        df[f"{col}_RunningMean"] = running_average(
            df[col]
        )

        df[f"{col}_RunningStd"] = running_std(
            df[col]
        )

    outfile = OUT / f"{T}_running_statistics.csv"

    df.to_csv(
        outfile,
        index=False
    )

    summary.append({

        "Temperature": T,

        "Final_Rg_RunningMean":
            df["Rg_RunningMean"].iloc[-1],

        "Final_PE_RunningMean":
            df["PotEng_RunningMean"].iloc[-1],

        "Final_RMSD_RunningMean":
            df["RMSD_RunningMean"].iloc[-1],

        "Final_Rg_RunningStd":
            df["Rg_RunningStd"].iloc[-1],

        "Final_PE_RunningStd":
            df["PotEng_RunningStd"].iloc[-1]

    })

    fig, ax = plt.subplots(
        figsize=(8,5)
    )

    ax.plot(
        df["Rg_RunningMean"],
        label="Rg"
    )

    ax.plot(
        df["PotEng_RunningMean"],
        label="Potential Energy"
    )

    ax.set_title(
        f"{T} Running Means"
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        OUT/f"{T}_RunningMeans.png",
        dpi=300
    )

    plt.close(fig)

summary = pd.DataFrame(summary)

summary.to_csv(
    OUT/"running_summary.csv",
    index=False
)

(OUT/"running_summary.json").write_text(
    json.dumps(
        summary.to_dict(
            orient="records"
        ),
        indent=2
    )
)

print("="*90)
print("DAY047 / PHASE2-A44")
print("RUNNING STATISTICS")
print("="*90)
print(OUT)
