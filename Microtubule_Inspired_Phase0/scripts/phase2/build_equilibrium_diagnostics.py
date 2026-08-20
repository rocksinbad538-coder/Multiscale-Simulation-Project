#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

EQ = ROOT/"runs"/"phase2"/"campaign_phase5_corrected"/"equilibrium_detection"

temperatures=[150,200,250,300,350]

fig,axes=plt.subplots(
    3,
    1,
    figsize=(12,12),
    sharex=False
)

colors={
    150:"tab:blue",
    200:"tab:orange",
    250:"tab:green",
    300:"tab:red",
    350:"tab:purple"
}

for T in temperatures:

    run=pd.read_csv(EQ/f"{T}K_running.csv")

    ene=pd.read_csv(EQ/f"{T}K_energy_running.csv")

    axes[0].plot(
        run["Step"],
        run["RgRunning"],
        label=f"{T} K",
        color=colors[T],
        linewidth=1.5
    )

    axes[1].plot(
        run["Step"],
        run["RMSDRunning"],
        label=f"{T} K",
        color=colors[T],
        linewidth=1.5
    )

    axes[2].plot(
        ene["Step"],
        ene["PotEngRunning"],
        label=f"{T} K",
        color=colors[T],
        linewidth=1.5
    )

axes[0].set_title("Running Radius of Gyration")
axes[0].set_ylabel("Rg (Å)")

axes[1].set_title("Running RMSD")
axes[1].set_ylabel("RMSD (Å)")

axes[2].set_title("Running Potential Energy")
axes[2].set_ylabel("Potential Energy")
axes[2].set_xlabel("Simulation Step")

for ax in axes:
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right")

plt.tight_layout()

outfile=EQ/"EquilibriumDiagnostics.png"

plt.savefig(
    outfile,
    dpi=300
)

plt.close()

print("="*90)
print("DAY049 / PHASE2-B07")
print("EQUILIBRIUM DIAGNOSTICS")
print("="*90)
print(outfile)
