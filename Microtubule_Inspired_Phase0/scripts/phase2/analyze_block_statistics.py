#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from md_analysis.statistics import block_statistics

ROOT = Path(__file__).resolve().parents[2]

CAMPAIGN = ROOT / "runs" / "phase2" / "campaign_phase5_corrected"

OUT = CAMPAIGN / "block_statistics"

OUT.mkdir(exist_ok=True)

summary = []

for folder in sorted(CAMPAIGN.glob("*K")):

    temperature = folder.name

    traj = pd.read_csv(
        folder / "analysis" / "trajectory_summary.csv"
    )

    thermo = pd.read_csv(
        folder / "analysis" / "thermodynamics.csv"
    )

    aligned = pd.read_csv(
        folder / "analysis" / "aligned_rmsd.csv"
    )

    df = pd.DataFrame({

        "Rg": traj["Rg"],

        "RMSD": traj["RMSD"],

        "AlignedRMSD": aligned["AlignedRMSD"],

        "PotentialEnergy": thermo["PotEng"],

        "Temperature": thermo["Temp"]

    })

    stats = block_statistics(

        df,

        columns=[
            "Rg",
            "RMSD",
            "AlignedRMSD",
            "PotentialEnergy",
            "Temperature"
        ],

        nblocks=5

    )

    outfile = OUT / f"{temperature}_block_statistics.csv"

    stats.to_csv(
        outfile,
        index=False
    )

    summary.append({

        "Temperature": temperature,

        "Rg_last_block":

            stats.iloc[-1]["Rg_mean"],

        "Rg_std_last_block":

            stats.iloc[-1]["Rg_std"],

        "PE_last_block":

            stats.iloc[-1]["PotentialEnergy_mean"],

        "PE_std_last_block":

            stats.iloc[-1]["PotentialEnergy_std"]

    })

summary = pd.DataFrame(summary)

summary.to_csv(

    OUT / "campaign_block_summary.csv",

    index=False

)

(OUT / "campaign_block_summary.json").write_text(

    json.dumps(

        summary.to_dict(orient="records"),

        indent=2

    )

)

print("="*90)
print("DAY047 / PHASE2-A43")
print("BLOCK STATISTICS")
print("="*90)
print(OUT)
