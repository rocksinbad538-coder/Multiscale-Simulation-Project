#!/usr/bin/env python3
from pathlib import Path
import csv, math

gro = Path("runs/phase1A/hybrid_hydrated_gromacs_min/min032.gro")
out = Path("parameters/phase1A/hybrid_hydrated_gromacs/audit_pyrene_molecule_geometries.csv")

lines = gro.read_text().splitlines()
natoms = int(lines[1].strip())
atoms = []

for line in lines[2:2+natoms]:
    resname = line[5:10].strip()
    if resname == "PYR":
        atoms.append({
            "idx": int(line[15:20]),
            "resnr": int(line[0:5]),
            "name": line[10:15].strip(),
            "x": float(line[20:28]),
            "y": float(line[28:36]),
            "z": float(line[36:44]),
        })

def dist(a,b):
    return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)

with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["chrom_resnr","atom_count","min_pair_nm","min_pair","max_pair_nm","mean_CC_bond_like_nm","n_CC_bond_like"])

    for resnr in sorted(set(a["resnr"] for a in atoms)):
        mol = [a for a in atoms if a["resnr"] == resnr]
        pairs = []
        cc_bond_like = []
        for i,a in enumerate(mol):
            for b in mol[i+1:]:
                d = dist(a,b)
                pairs.append((d,a,b))
                if a["name"].startswith("C") and b["name"].startswith("C") and 0.13 <= d <= 0.16:
                    cc_bond_like.append(d)
        minp = min(pairs, key=lambda x: x[0])
        maxp = max(pairs, key=lambda x: x[0])
        mean_cc = sum(cc_bond_like)/len(cc_bond_like) if cc_bond_like else float("nan")
        w.writerow([
            resnr,
            len(mol),
            f"{minp[0]:.5f}",
            f"{minp[1]['idx']}:{minp[1]['name']}-{minp[2]['idx']}:{minp[2]['name']}",
            f"{maxp[0]:.5f}",
            f"{mean_cc:.5f}",
            len(cc_bond_like),
        ])

print(out.read_text())
