#!/usr/bin/env python3
from pathlib import Path
import csv, math

gro = Path("runs/phase1A/hybrid_hydrated_gromacs_min/min032.gro")
out = Path("parameters/phase1A/hybrid_hydrated_gromacs/audit_pyrene_force_region.csv")

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
    return math.sqrt(sum([
        pbc_delta(a["x"]-b["x"], box[0])**2,
        pbc_delta(a["y"]-b["y"], box[1])**2,
        pbc_delta(a["z"]-b["z"], box[2])**2,
    ]))

targets = [1773, 1774]
with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["target_atom","neighbor_atom","target_res","target_name","neighbor_res","neighbor_name","distance_nm"])

    for tid in targets:
        t = next(a for a in atoms if a["idx"] == tid)
        neigh = []
        for a in atoms:
            if a["idx"] == tid:
                continue
            d = dist(t,a)
            if d < 0.50:
                neigh.append((d,a))
        for d,a in sorted(neigh):
            w.writerow([tid,a["idx"],t["resname"],t["name"],a["resname"],a["name"],f"{d:.5f}"])

print(f"Wrote {out}")
