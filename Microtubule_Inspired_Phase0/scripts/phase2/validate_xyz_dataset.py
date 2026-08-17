#!/usr/bin/env python3

from pathlib import Path
import json
import math

ROOT = Path(__file__).resolve().parents[2]

XYZDIR = ROOT/"runs"/"phase2"/"campaign"/"tddft_xyz"

summary=[]

for xyz in sorted(XYZDIR.glob("*.xyz")):

    with open(xyz) as f:

        lines=f.readlines()

    natoms=int(lines[0])

    coords=[]

    composition={}

    valid=True

    for line in lines[2:]:

        s=line.split()

        if len(s)!=4:

            valid=False
            continue

        elem=s[0]

        composition[elem]=composition.get(elem,0)+1

        x=float(s[1])
        y=float(s[2])
        z=float(s[3])

        if not(
            math.isfinite(x)
            and
            math.isfinite(y)
            and
            math.isfinite(z)
        ):
            valid=False

        coords.append((x,y,z))

    summary.append({

        "file":xyz.name,

        "natoms":natoms,

        "composition":composition,

        "valid":valid,

        "parsed_atoms":len(coords)

    })

with open(
    XYZDIR/"XYZ_VALIDATION.json",
    "w"
) as f:

    json.dump(
        summary,
        f,
        indent=2
    )

print("="*90)
print("PHASE3-A06")
print("XYZ VALIDATION")
print("="*90)

for s in summary:

    print(
        f"{s['file']:30s}",
        f"atoms={s['natoms']:2d}",
        f"parsed={s['parsed_atoms']:2d}",
        f"valid={s['valid']}"
    )

print()

print(XYZDIR/"XYZ_VALIDATION.json")
