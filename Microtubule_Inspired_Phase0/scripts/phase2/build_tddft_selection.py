#!/usr/bin/env python3

from pathlib import Path
import json
import shutil
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

REP = ROOT/"runs"/"phase2"/"campaign_phase5_corrected"/"snapshot_clustering"/"cluster_representatives.csv"

SOURCE = ROOT/"runs"/"phase2"/"campaign_phase5_corrected"/"representative_ensemble"

OUT = ROOT/"runs"/"phase2"/"campaign_phase5_corrected"/"tddft_selected_structures"

OUT.mkdir(exist_ok=True)

rep = pd.read_csv(REP)

manifest = []

for _, row in rep.iterrows():

    T = int(row.temperature_K)

    snap = row.snapshot

    src = SOURCE/f"{T}K"/f"{snap}.dump"

    dst = OUT/f"{T}K_{snap}.dump"

    shutil.copy2(src, dst)

    manifest.append({
        "cluster": int(row.cluster),
        "temperature_K": T,
        "snapshot": snap,
        "file": dst.name
    })

with open(
    OUT/"TDDFT_SELECTION.json",
    "w"
) as f:
    json.dump(manifest, f, indent=2)

print("="*90)
print("PHASE3-A04")
print("TDDFT SELECTION")
print("="*90)
print("Selected structures:", len(manifest))
print(OUT)
