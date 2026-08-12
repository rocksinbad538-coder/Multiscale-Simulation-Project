#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import pathlib

import matplotlib.pyplot as plt
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]

TRAJECTORY = (
    ROOT
    / "runs"
    / "phase2"
    / "day044_md_protocol"
    / "production.xyz"
)

OUT = (
    ROOT
    / "runs"
    / "phase2"
    / "day045_md_analysis"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

frames = []

with open(TRAJECTORY) as f:

    lines = f.readlines()

i = 0

while i < len(lines):

    if not lines[i].startswith("ITEM: TIMESTEP"):
        i += 1
        continue

    timestep = int(lines[i+1])

    natoms = int(lines[i+3])

    atom_start = i + 9

    atoms = []

    for line in lines[atom_start:atom_start+natoms]:

        s = line.split()

        atoms.append((
            int(s[0]),
            float(s[2]),
            float(s[3]),
            float(s[4]),
        ))

    frames.append(atoms)

    i = atom_start + natoms

natoms = len(frames[0])

means = []

for atom in range(natoms):

    xs = [f[atom][1] for f in frames]
    ys = [f[atom][2] for f in frames]
    zs = [f[atom][3] for f in frames]

    means.append((
        sum(xs)/len(xs),
        sum(ys)/len(ys),
        sum(zs)/len(zs)
    ))

rows=[]

for atom in range(natoms):

    mx,my,mz = means[atom]

    s=0.0

    for frame in frames:

        dx=frame[atom][1]-mx
        dy=frame[atom][2]-my
        dz=frame[atom][3]-mz

        s+=dx*dx+dy*dy+dz*dz

    rmsf=(s/len(frames))**0.5

    rows.append({

        "atom":atom+1,
        "RMSF_A":rmsf

    })

csvfile=OUT/"atom_rmsf.csv"

pd.DataFrame(rows).to_csv(
    csvfile,
    index=False
)

plt.figure(figsize=(8,4))

plt.plot(
    [r["atom"] for r in rows],
    [r["RMSF_A"] for r in rows],
    linewidth=1.5
)

plt.xlabel("Atom")

plt.ylabel("RMSF (Å)")

plt.tight_layout()

plt.savefig(
    OUT/"RMSF_per_atom.png",
    dpi=300
)

plt.close()

report={

    "atom_count":natoms,

    "maximum_RMSF_A":max(r["RMSF_A"] for r in rows),

    "mean_RMSF_A":sum(r["RMSF_A"] for r in rows)/len(rows)

}

json.dump(
    report,
    open(OUT/"RMSF_REPORT.json","w"),
    indent=2
)

print("="*90)
print("DAY045 / PHASE2-A20")
print("RMSF ANALYSIS")
print("="*90)
print()

print(csvfile)
print(OUT/"RMSF_REPORT.json")
print(OUT/"RMSF_per_atom.png")
