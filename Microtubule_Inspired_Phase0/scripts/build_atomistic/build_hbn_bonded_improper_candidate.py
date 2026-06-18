#!/usr/bin/env python3
from pathlib import Path
import math

base = Path("parameters/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032")
src = base / "hbn_bonded_candidate_kang2000.itp"
gro = Path("runs/phase1A/hybrid_hbnBonded_kang2000_min/min_hbnBonded_kang2000.gro")
out = base / "hbn_bonded_candidate_kang2000_improper100.itp"

k_improper = 100.0  # kJ/mol/rad^2, first conservative out-of-plane stiffness

lines = gro.read_text().splitlines()
natoms = int(lines[1].strip())

atoms = []
for line in lines[2:2+natoms]:
    if line[5:10].strip() != "HBN":
        continue
    atoms.append({
        "local": len(atoms) + 1,
        "name": line[10:15].strip(),
        "x": float(line[20:28]),
        "y": float(line[28:36]),
        "z": float(line[36:44]),
    })

def dist(a,b):
    return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)

neighbors = {a["local"]: [] for a in atoms}

for i,a in enumerate(atoms):
    for b in atoms[i+1:]:
        d = dist(a,b)
        if d < 0.20:
            neighbors[a["local"]].append(b["local"])
            neighbors[b["local"]].append(a["local"])

impropers = []
for center, ns in neighbors.items():
    if len(ns) == 3:
        i, j, k = sorted(ns)
        # GROMACS proper/improper function type 2: harmonic improper
        # ordering: neighbor1 center neighbor2 neighbor3
        impropers.append((i, center, j, k))

text = src.read_text().rstrip()

# remove old POSRES block temporarily
pre, posres = text.split("#ifdef POSRES_HBN", 1)

with out.open("w") as f:
    f.write(pre.rstrip() + "\n\n")
    f.write("[ dihedrals ]\n")
    f.write("; ai  aj  ak  al  funct  phi0(deg)  k\n")
    for i,j,k,l in impropers:
        f.write(f"{i:6d} {j:6d} {k:6d} {l:6d}  2  {0.0:10.4f}  {k_improper:10.4f}\n")
    f.write("\n#ifdef POSRES_HBN")
    f.write(posres)

print("Wrote", out)
print("Improper terms:", len(impropers))
