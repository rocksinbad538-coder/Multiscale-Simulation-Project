#!/usr/bin/env python3
from pathlib import Path
import csv, math

gro = Path("runs/phase1A/hybrid_hydrated_gromacs_min/min032.gro")
out = Path("parameters/phase1A/hybrid_hydrated_gromacs/audit_pyrene_hbn_distances.csv")

lines = gro.read_text().splitlines()
natoms = int(lines[1].strip())
box = [float(x) for x in lines[2+natoms].split()[:3]]

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

def pbc_delta(dx,L):
    return dx - L*round(dx/L)

def dist(a,b):
    dx = pbc_delta(a["x"]-b["x"], box[0])
    dy = pbc_delta(a["y"]-b["y"], box[1])
    dz = pbc_delta(a["z"]-b["z"], box[2])
    return math.sqrt(dx*dx + dy*dy + dz*dz)

hbn = [a for a in atoms if a["resname"] == "HBN"]
pyr = [a for a in atoms if a["resname"] == "PYR"]

with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["pyrene_resnr","min_pyr_hbn_nm","pyr_atom","hbn_atom","contacts_lt_0p35","contacts_lt_0p40","contacts_lt_0p45"])
    for resnr in sorted(set(a["resnr"] for a in pyr)):
        mol = [a for a in pyr if a["resnr"] == resnr]
        best = (999,None,None)
        c35 = c40 = c45 = 0
        for p in mol:
            for h in hbn:
                d = dist(p,h)
                if d < best[0]:
                    best = (d,p,h)
                if d < 0.35: c35 += 1
                if d < 0.40: c40 += 1
                if d < 0.45: c45 += 1
        d,p,h = best
        w.writerow([
            resnr,
            f"{d:.5f}",
            f"{p['idx']}:{p['name']}",
            f"{h['idx']}:{h['name']}",
            c35,c40,c45
        ])

print(out.read_text())
