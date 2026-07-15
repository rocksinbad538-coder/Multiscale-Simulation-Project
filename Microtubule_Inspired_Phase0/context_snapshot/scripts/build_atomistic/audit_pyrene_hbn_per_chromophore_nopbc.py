#!/usr/bin/env python3
from pathlib import Path
import csv, math, sys

gro = Path(sys.argv[1])
out = Path(sys.argv[2])

lines = gro.read_text().splitlines()
natoms = int(lines[1].strip())

atoms = []
for line in lines[2:2+natoms]:
    atoms.append({
        "idx": int(line[15:20]),
        "resnr": int(line[0:5]),
        "resname": line[5:10].strip(),
        "name": line[10:15].strip(),
        "x": float(line[20:28]),
        "y": float(line[28:36]),
        "z": float(line[36:44]),
    })

def dist(a,b):
    return math.sqrt(
        (a["x"]-b["x"])**2 +
        (a["y"]-b["y"])**2 +
        (a["z"]-b["z"])**2
    )

hbn = [a for a in atoms if a["resname"] == "HBN"]
pyrs = {}
for a in atoms:
    if a["resname"] == "PYR":
        pyrs.setdefault(a["resnr"], []).append(a)

with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow([
        "source_gro","pyrene_resnr","atom_count",
        "min_pyr_hbn_nm","pyr_atom","hbn_atom",
        "contacts_lt_0p30","contacts_lt_0p35","contacts_lt_0p40","contacts_lt_0p45"
    ])

    for resnr, pyr in sorted(pyrs.items()):
        best = (999,None,None)
        c030 = c035 = c040 = c045 = 0
        for p in pyr:
            for h in hbn:
                d = dist(p,h)
                if d < best[0]:
                    best = (d,p,h)
                if d < 0.30: c030 += 1
                if d < 0.35: c035 += 1
                if d < 0.40: c040 += 1
                if d < 0.45: c045 += 1

        d,p,h = best
        w.writerow([
            str(gro),resnr,len(pyr),
            f"{d:.5f}",f"{p['idx']}:{p['name']}",f"{h['idx']}:{h['name']}",
            c030,c035,c040,c045
        ])

print(out)
