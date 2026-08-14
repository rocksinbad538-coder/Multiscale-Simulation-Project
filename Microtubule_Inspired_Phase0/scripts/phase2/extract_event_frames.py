#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import shutil

ROOT = Path(__file__).resolve().parents[2]

analysis = ROOT/"runs"/"phase2"/"day045_md_analysis"

events = pd.read_csv(
    analysis/"incremental_rmsd_largest_events.csv"
)

traj = ROOT/"runs"/"phase2"/"campaign"/"300K"/"production.xyz"

print("="*90)
print("DAY047 / PHASE2-A39")
print("EVENT FRAME LIST")
print("="*90)

print()

print("Trajectory:")
print(traj)

print()

print("Largest events:")

print(events[["timestep","IncrementalRMSD"]])
