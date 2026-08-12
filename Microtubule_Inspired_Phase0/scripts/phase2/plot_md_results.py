#!/usr/bin/env python3

from __future__ import annotations

import pathlib

import pandas as pd
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[2]

TRAJ = ROOT/"runs/phase2/day045_md_analysis/trajectory_summary.csv"
THERMO = ROOT/"runs/phase2/day045_md_analysis/thermodynamics.csv"

OUT = ROOT/"runs/phase2/day045_md_analysis"

traj = pd.read_csv(TRAJ)
thermo = pd.read_csv(THERMO)


def make_plot(x,y,xlabel,ylabel,title,filename):

    plt.figure(figsize=(7,4))

    plt.plot(x,y,linewidth=1.5)

    plt.xlabel(xlabel)

    plt.ylabel(ylabel)

    plt.title(title)

    plt.tight_layout()

    plt.savefig(OUT/filename,dpi=300)

    plt.close()


make_plot(
    traj["timestep"],
    traj["RMSD"],
    "Timestep",
    "RMSD (Å)",
    "RMSD vs Time",
    "RMSD_vs_Time.png"
)

make_plot(
    traj["timestep"],
    traj["Rg"],
    "Timestep",
    "Radius of Gyration (Å)",
    "Radius of Gyration",
    "Rg_vs_Time.png"
)

make_plot(
    thermo["Step"],
    thermo["Temp"],
    "Timestep",
    "Temperature (K)",
    "Temperature",
    "Temperature_vs_Time.png"
)

make_plot(
    thermo["Step"],
    thermo["PotEng"],
    "Timestep",
    "Potential Energy",
    "Potential Energy",
    "PotentialEnergy_vs_Time.png"
)

make_plot(
    thermo["Step"],
    thermo["Press"],
    "Timestep",
    "Pressure",
    "Pressure",
    "Pressure_vs_Time.png"
)

print("="*90)
print("DAY045 / PHASE2-A19")
print("MD PLOTS GENERATED")
print("="*90)

for f in sorted(OUT.glob("*.png")):
    print(f)
