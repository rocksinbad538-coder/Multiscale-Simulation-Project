#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

XYZROOT = ROOT/"runs"/"phase2"/"campaign_phase5_corrected"/"tddft_xyz"

rows=[]

for xyz in sorted(XYZROOT.glob("*.xyz")):

    with open(xyz) as f:

        natoms = int(f.readline())
        comment = f.readline()

        elements=[]

        for line in f:

            if not line.strip():
                continue

            elements.append(line.split()[0])

    c = Counter(elements)

    row={

        "structure":xyz.stem,

        "natoms":natoms

    }

    for e,n in sorted(c.items()):
        row[e]=n

    rows.append(row)

df=pd.DataFrame(rows).fillna(0)

OUT=ROOT/"runs"/"phase2"/"campaign_phase5_corrected"/"audit"

OUT.mkdir(exist_ok=True)

csv=OUT/"xyz_composition_audit.csv"

df.to_csv(csv,index=False)

print("="*90)
print("PHASE4-C01")
print("XYZ COMPOSITION AUDIT")
print("="*90)

print(df)

print()

print(csv)
