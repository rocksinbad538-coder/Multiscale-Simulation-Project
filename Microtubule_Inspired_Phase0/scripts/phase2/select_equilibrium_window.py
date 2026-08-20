#!/usr/bin/env python3

from pathlib import Path

import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

CAMPAIGN = ROOT/"runs"/"phase2"/"campaign_phase5_corrected"

OUT = CAMPAIGN/"equilibrium_detection"

OUT.mkdir(exist_ok=True)

WINDOW = 500


def running_mean(x, window):

    return (
        pd.Series(x)
        .rolling(
            window,
            center=True,
            min_periods=1
        )
        .mean()
        .to_numpy()
    )


def running_slope(x, y, window):

    x = np.asarray(x)

    y = np.asarray(y)

    slopes = np.zeros(len(y))

    half = window//2

    for i in range(len(y)):

        a=max(0,i-half)

        b=min(len(y),i+half)

        xx=x[a:b]

        yy=y[a:b]

        if len(xx)<5:

            slopes[i]=0.0

            continue

        m=np.polyfit(xx,yy,1)[0]

        slopes[i]=m

    return slopes


results=[]

temperatures=[150,200,250,300,350]

for T in temperatures:

    analysis=CAMPAIGN/f"{T}K"/"analysis"

    traj=pd.read_csv(
        analysis/"trajectory_summary.csv"
    )

    thermo=pd.read_csv(
        analysis/"thermodynamics.csv"
    )

    rg_mean=running_mean(
        traj["Rg"],
        WINDOW
    )

    rmsd_mean=running_mean(
        traj["RMSD"],
        WINDOW
    )

    pe_mean=running_mean(
        thermo["PotEng"],
        WINDOW
    )

    rg_slope=running_slope(
        traj["timestep"],
        rg_mean,
        WINDOW
    )

    rmsd_slope=running_slope(
        traj["timestep"],
        rmsd_mean,
        WINDOW
    )

    pe_slope=running_slope(
        thermo["Step"],
        pe_mean,
        WINDOW
    )

    pd.DataFrame({

        "Step":traj["timestep"],

        "Rg":traj["Rg"],

        "RgRunning":rg_mean,

        "RgSlope":rg_slope,

        "RMSD":traj["RMSD"],

        "RMSDRunning":rmsd_mean,

        "RMSDSlope":rmsd_slope

    }).to_csv(

        OUT/f"{T}K_running.csv",

        index=False

    )

    pd.DataFrame({

        "Step":thermo["Step"],

        "PotEng":thermo["PotEng"],

        "PotEngRunning":pe_mean,

        "PotEngSlope":pe_slope

    }).to_csv(

        OUT/f"{T}K_energy_running.csv",

        index=False

    )

    results.append({

        "Temperature":T,

        "RgSlopeMean":float(np.mean(np.abs(rg_slope))),

        "RMSDSlopeMean":float(np.mean(np.abs(rmsd_slope))),

        "PotEngSlopeMean":float(np.mean(np.abs(pe_slope)))

    })

summary=pd.DataFrame(results)

summary.to_csv(

    OUT/"running_slope_summary.csv",

    index=False

)

summary.to_json(

    OUT/"running_slope_summary.json",

    orient="records",

    indent=2

)

print("="*90)
print("DAY049 / PHASE2-B06")
print("RUNNING-SLOPE ANALYSIS")
print("="*90)
print(OUT)
print()
print(summary)
