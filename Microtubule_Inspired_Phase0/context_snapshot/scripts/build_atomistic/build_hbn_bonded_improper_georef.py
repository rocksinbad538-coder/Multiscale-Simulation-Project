#!/usr/bin/env python3
from pathlib import Path
import math

base = Path("parameters/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032")
src = base / "hbn_bonded_candidate_kang2000.itp"
gro = Path("runs/phase1A/hybrid_hbnBonded_kang2000_min/min_hbnBonded_kang2000.gro")
out = base / "hbn_bonded_candidate_kang2000_improperGeo100.itp"

k_improper = 100.0

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

def sub(a,b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def dot(a,b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def cross(a,b):
    return (
        a[1]*b[2]-a[2]*b[1],
        a[2]*b[0]-a[0]*b[2],
        a[0]*b[1]-a[1]*b[0],
    )

def norm(a):
    return math.sqrt(dot(a,a))

def dihedral(p1,p2,p3,p4):
    b1 = sub(p2,p1)
    b2 = sub(p3,p2)
    b3 = sub(p4,p3)

    n1 = cross(b1,b2)
    n2 = cross(b2,b3)

    if norm(n1) < 1e-12 or norm(n2) < 1e-12:
        return 0.0

    n1 = tuple(x/norm(n1) for x in n1)
    n2 = tuple(x/norm(n2) for x in n2)
    b2u = tuple(x/norm(b2) for x in b2)

    m1 = cross(n1,b2u)
    x = dot(n1,n2)
    y = dot(m1,n2)
    return math.degrees(math.atan2(y,x))

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
        i, k, l = sorted(ns)
        # ordering: neighbor1, center, neighbor2, neighbor3
        ai = atoms[i-1]
        aj = atoms[center-1]
        ak = atoms[k-1]
        al = atoms[l-1]
        p1 = (ai["x"],ai["y"],ai["z"])
        p2 = (aj["x"],aj["y"],aj["z"])
        p3 = (ak["x"],ak["y"],ak["z"])
        p4 = (al["x"],al["y"],al["z"])
        phi0 = dihedral(p1,p2,p3,p4)
        impropers.append((i, center, k, l, phi0))

text = src.read_text().rstrip()
pre, posres = text.split("#ifdef POSRES_HBN", 1)

with out.open("w") as f:
    f.write(pre.rstrip() + "\n\n")
    f.write("[ dihedrals ]\n")
    f.write("; ai  aj  ak  al  funct  phi0(deg)  k\n")
    for i,j,k,l,phi0 in impropers:
        f.write(f"{i:6d} {j:6d} {k:6d} {l:6d}  2  {phi0:10.4f}  {k_improper:10.4f}\n")
    f.write("\n#ifdef POSRES_HBN")
    f.write(posres)

phis = [x[4] for x in impropers]
print("Wrote", out)
print("Improper terms:", len(impropers))
print("phi0 min/max/mean:", min(phis), max(phis), sum(phis)/len(phis))
