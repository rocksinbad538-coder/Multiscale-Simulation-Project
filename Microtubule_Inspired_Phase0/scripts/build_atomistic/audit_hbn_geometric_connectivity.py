#!/usr/bin/env python3
from pathlib import Path
import csv, math

gro = Path("runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_min/min_gap45_pyr5shift_clean032_hbnSoftLJ_posresPYR50k.gro")
out = Path("runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_min/audit_hbn_geometric_connectivity.csv")

lines = gro.read_text().splitlines()
natoms = int(lines[1].strip())

atoms = []
for line in lines[2:2+natoms]:
    resname = line[5:10].strip()
    if resname != "HBN":
        continue
    atoms.append({
        "idx": int(line[15:20]),
        "name": line[10:15].strip(),
        "elem": line[10:15].strip()[0],
        "x": float(line[20:28]),
        "y": float(line[28:36]),
        "z": float(line[36:44]),
    })

def dist(a,b):
    return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)

pairs = []
degree = {a["idx"]: 0 for a in atoms}

for i,a in enumerate(atoms):
    for b in atoms[i+1:]:
        d = dist(a,b)
        if d < 0.20:  # nm = 2.0 Å, candidate covalent BN/BH/NH range
            pairs.append((a["idx"], b["idx"], a["name"], b["name"], a["elem"], b["elem"], d))
            degree[a["idx"]] += 1
            degree[b["idx"]] += 1

with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["hbn_atoms", len(atoms)])
    w.writerow(["candidate_bonds_lt_0p20nm", len(pairs)])
    w.writerow(["degree0_atoms", sum(1 for v in degree.values() if v == 0)])
    w.writerow(["degree1_atoms", sum(1 for v in degree.values() if v == 1)])
    w.writerow(["degree2_atoms", sum(1 for v in degree.values() if v == 2)])
    w.writerow(["degree3_atoms", sum(1 for v in degree.values() if v == 3)])
    w.writerow(["degree_gt3_atoms", sum(1 for v in degree.values() if v > 3)])
    w.writerow([])
    w.writerow(["atom_i","atom_j","name_i","name_j","elem_i","elem_j","distance_nm"])
    for p in sorted(pairs, key=lambda x: x[-1]):
        w.writerow([p[0],p[1],p[2],p[3],p[4],p[5],f"{p[6]:.5f}"])

print(out)
print(out.read_text().splitlines()[:15])
