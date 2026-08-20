#!/usr/bin/env python3

from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[2]

ENSEMBLE = ROOT/"runs"/"phase2"/"campaign_phase5_corrected"/"representative_ensemble"

OUT = ROOT/"runs"/"phase2"/"campaign_phase5_corrected"/"tddft_manifest"

OUT.mkdir(exist_ok=True)

manifest=[]

for folder in sorted(ENSEMBLE.glob("*K")):

    T=int(folder.name[:-1])

    for snap in sorted(folder.glob("snapshot_*.dump")):

        manifest.append({

            "temperature_K":T,

            "snapshot":snap.stem,

            "path":str(snap.relative_to(ROOT)),

            "electronic_status":"pending",

            "tddft_status":"pending",

            "exciton_status":"pending"

        })

with open(
    OUT/"TDDFT_INPUT_MANIFEST.json",
    "w"
) as f:

    json.dump(
        manifest,
        f,
        indent=2
    )

print("="*90)
print("DAY049 / PHASE2-B12")
print("TDDFT MANIFEST")
print("="*90)
print("Configurations :",len(manifest))
print(OUT)
