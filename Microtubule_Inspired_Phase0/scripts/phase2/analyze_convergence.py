#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = ROOT / "runs" / "phase2" / "campaign"

TAIL_WINDOWS = [0.10, 0.20, 0.30]

PE_SLOPE_THRESHOLD = 1e-4
RG_CV_THRESHOLD = 0.02
RMSD_CV_THRESHOLD = 0.05




def find_column(df, candidates):

    for c in candidates:
        if c in df.columns:
            return c

    raise RuntimeError(
        f"None of the candidate columns {candidates} "
        f"were found. Available columns: {list(df.columns)}"
    )


def tail_statistics(series, fraction):

    n = len(series)
    start = int(n * (1.0 - fraction))
    tail = np.asarray(series[start:], dtype=float)

    x = np.arange(len(tail))
    slope = np.polyfit(x, tail, 1)[0]

    mean = float(np.mean(tail))
    std = float(np.std(tail))
    cv = std / abs(mean) if abs(mean) > 1e-12 else 0.0

    return {
        "mean": mean,
        "std": std,
        "cv": cv,
        "slope": float(slope),
    }


summary = []

for folder in sorted(CAMPAIGN.glob("*K")):

    analysis = folder / "analysis"

    thermo = analysis / "thermodynamics.csv"
    traj = analysis / "trajectory_summary.csv"
    rmsd = analysis / "aligned_rmsd.csv"

    if not (thermo.exists() and traj.exists() and rmsd.exists()):
        print(f"Skipping {folder.name}")
        continue

    thermo_df = pd.read_csv(thermo)
    traj_df = pd.read_csv(traj)
    rmsd_df = pd.read_csv(rmsd)

    pe_col = find_column(
        thermo_df,
        ["PotEng","PE","pe"]
    )

    rg_col = find_column(
        traj_df,
        ["Rg","Rg_A","RadiusGyration"]
    )

    rmsd_col = find_column(
        rmsd_df,
        ["AlignedRMSD","Aligned_RMSD","Aligned_RMSD_A","RMSD"]
    )

    pe = thermo_df[pe_col]
    rg = traj_df[rg_col]
    rmsd_series = rmsd_df[rmsd_col]

    result = {
        "Temperature": folder.name
    }

    pe_ok = True
    rg_ok = True
    rmsd_ok = True

    for frac in TAIL_WINDOWS:

        tag = f"{int(frac*100)}"

        pe_stats = tail_statistics(pe, frac)
        rg_stats = tail_statistics(rg, frac)
        rmsd_stats = tail_statistics(rmsd_series, frac)

        result[f"PE_mean_last{tag}"] = pe_stats["mean"]
        result[f"PE_slope_last{tag}"] = pe_stats["slope"]

        result[f"Rg_mean_last{tag}"] = rg_stats["mean"]
        result[f"Rg_cv_last{tag}"] = rg_stats["cv"]

        result[f"RMSD_mean_last{tag}"] = rmsd_stats["mean"]
        result[f"RMSD_cv_last{tag}"] = rmsd_stats["cv"]

        pe_ok &= abs(pe_stats["slope"]) < PE_SLOPE_THRESHOLD
        rg_ok &= rg_stats["cv"] < RG_CV_THRESHOLD
        rmsd_ok &= rmsd_stats["cv"] < RMSD_CV_THRESHOLD

    result["PotentialEnergy_Converged"] = pe_ok
    result["Rg_Converged"] = rg_ok
    result["RMSD_Converged"] = rmsd_ok

    result["Simulation_Stable"] = pe_ok and rg_ok and rmsd_ok

    summary.append(result)


summary_df = pd.DataFrame(summary)

csv_file = CAMPAIGN / "convergence_summary.csv"
json_file = CAMPAIGN / "convergence_summary.json"
md_file = CAMPAIGN / "CONVERGENCE_REPORT.md"

summary_df.to_csv(csv_file, index=False)

json_file.write_text(
    json.dumps(summary, indent=2)
)

with open(md_file, "w") as f:

    f.write("# Phase 2 Convergence Analysis\n\n")

    f.write(summary_df.to_csv(index=False))

    f.write("\n")

    stable = summary_df["Simulation_Stable"].sum()

    f.write(f"Stable simulations: {stable}/{len(summary_df)}\n")


print("=" * 90)
print("DAY047 / PHASE2-A34")
print("CONVERGENCE ANALYSIS")
print("=" * 90)
print(csv_file)
print(json_file)
print(md_file)
