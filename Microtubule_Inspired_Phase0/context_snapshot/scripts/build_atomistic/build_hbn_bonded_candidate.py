#!/usr/bin/env python3
from pathlib import Path
import math

base = Path("parameters/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032")
src_itp = base / "hbn_fixed_dummy.itp"
gro = Path("runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_min/min_gap45_pyr5shift_clean032_hbnSoftLJ_posresPYR50k.gro")
out_itp = base / "hbn_bonded_candidate.itp"

# Force constants: intentionally conservative first candidate
# Units: GROMACS bonds: kJ mol-1 nm-2; angles: kJ mol-1 rad-2
k_bond = 200000.0
k_angle = 500.0

lines = gro.read_text().splitlines()
natoms = int(lines[1].strip())

atoms = []
for line in lines[2:2+natoms]:
    if line[5:10].strip() != "HBN":
        continue
    idx_global = int(line[15:20])
    local_idx = len(atoms) + 1
    name = line[10:15].strip()
    atoms.append({
        "global": idx_global,
        "local": local_idx,
        "name": name,
        "elem": name[0],
        "x": float(line[20:28]),
        "y": float(line[28:36]),
        "z": float(line[36:44]),
    })

def dist(a,b):
    return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)

# Candidate covalent bonds from geometry
bonds = []
neighbors = {a["local"]: [] for a in atoms}

for i,a in enumerate(atoms):
    for b in atoms[i+1:]:
        d = dist(a,b)
        if d < 0.20:
            bonds.append((a["local"], b["local"], d))
            neighbors[a["local"]].append((b["local"], d))
            neighbors[b["local"]].append((a["local"], d))

# Angles i-j-k where j is central atom
angles = []
seen = set()
for j, neigh in neighbors.items():
    ns = [x[0] for x in neigh]
    for a in range(len(ns)):
        for b in range(a+1, len(ns)):
            i, k = ns[a], ns[b]
            key = tuple(sorted((i,j,k)))
            if key in seen:
                continue
            seen.add(key)

            ai = atoms[i-1]
            aj = atoms[j-1]
            ak = atoms[k-1]

            v1 = (ai["x"]-aj["x"], ai["y"]-aj["y"], ai["z"]-aj["z"])
            v2 = (ak["x"]-aj["x"], ak["y"]-aj["y"], ak["z"]-aj["z"])
            n1 = math.sqrt(sum(x*x for x in v1))
            n2 = math.sqrt(sum(x*x for x in v2))
            c = sum(v1[m]*v2[m] for m in range(3))/(n1*n2)
            c = max(-1.0, min(1.0, c))
            theta = math.degrees(math.acos(c))
            angles.append((i,j,k,theta))

src = src_itp.read_text().splitlines()

# Keep everything until before POSRES block, remove existing position restraints from source
kept = []
for line in src:
    if line.strip().startswith("#ifdef POSRES_HBN"):
        break
    kept.append(line)

with out_itp.open("w") as f:
    f.write("\n".join(kept).rstrip() + "\n\n")

    f.write("[ bonds ]\n")
    f.write("; ai  aj  funct  r0(nm)      k\n")
    for i,j,d in bonds:
        f.write(f"{i:6d} {j:6d}  1  {d:10.5f}  {k_bond:10.1f}\n")

    f.write("\n[ angles ]\n")
    f.write("; ai  aj  ak  funct  theta(deg)  k\n")
    for i,j,k,theta in angles:
        f.write(f"{i:6d} {j:6d} {k:6d}  1  {theta:10.4f}  {k_angle:10.1f}\n")

    f.write("\n#ifdef POSRES_HBN\n")
    f.write("[ position_restraints ]\n")
    f.write("; atom  type      fx       fy       fz\n")
    for i in range(1, len(atoms)+1):
        f.write(f"{i:6d}     1  1000000  1000000  1000000\n")
    f.write("#endif\n")

print("Wrote", out_itp)
print("HBN atoms:", len(atoms))
print("Bonds:", len(bonds))
print("Angles:", len(angles))
print("Degree 0:", sum(1 for v in neighbors.values() if len(v)==0))
print("Degree 1:", sum(1 for v in neighbors.values() if len(v)==1))
print("Degree 2:", sum(1 for v in neighbors.values() if len(v)==2))
print("Degree 3:", sum(1 for v in neighbors.values() if len(v)==3))
print("Degree >3:", sum(1 for v in neighbors.values() if len(v)>3))
