#!/usr/bin/env python3

from pathlib import Path
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

ENSEMBLE = ROOT/"runs"/"phase2"/"campaign"/"representative_ensemble"

OUT = ROOT/"runs"/"phase2"/"campaign"/"snapshot_descriptors"

OUT.mkdir(exist_ok=True)

rows=[]

for folder in sorted(ENSEMBLE.glob("*K")):

    T=int(folder.name[:-1])

    for snap in sorted(folder.glob("snapshot_*.dump")):

        with open(snap) as f:

            lines=f.readlines()

        atom_header=None

        for i,l in enumerate(lines):

            if l.startswith("ITEM: ATOMS"):

                atom_header=i
                break

        if atom_header is None:
            raise RuntimeError(snap)

        xyz=[]

        for line in lines[atom_header+1:]:

            s=line.split()

            x=float(s[4])
            y=float(s[5])
            z=float(s[6])

            xyz.append([x,y,z])

        xyz=np.asarray(xyz)

        xyz-=xyz.mean(axis=0)

        vector=xyz.reshape(-1)

        row={
            "temperature_K":T,
            "snapshot":snap.stem
        }

        for i,v in enumerate(vector):

            row[f"f{i:03d}"]=float(v)

        rows.append(row)

df=pd.DataFrame(rows)

df.to_csv(
    OUT/"snapshot_descriptor_matrix.csv",
    index=False
)

print("="*90)
print("DAY049 / PHASE3-A01")
print("SNAPSHOT DESCRIPTOR MATRIX")
print("="*90)

print(df.shape)

print()

print(OUT)
