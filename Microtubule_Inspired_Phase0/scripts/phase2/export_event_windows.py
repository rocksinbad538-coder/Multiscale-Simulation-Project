#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import json

ROOT = Path(__file__).resolve().parents[2]

analysis = ROOT/"runs"/"phase2"/"day045_md_analysis"

events = pd.read_csv(
    analysis/"incremental_rmsd_largest_events.csv"
)

WINDOW = 200

selected = []

for ts in events["timestep"]:

    selected.extend([
        ts-WINDOW,
        ts-100,
        ts,
        ts+100,
        ts+WINDOW,
    ])

selected = sorted(set(selected))

outfile = analysis/"event_windows.json"

outfile.write_text(
    json.dumps(selected, indent=2)
)

print("="*90)
print("DAY047 / PHASE2-A40")
print("EVENT WINDOWS")
print("="*90)
print(outfile)
print()
print(selected)
