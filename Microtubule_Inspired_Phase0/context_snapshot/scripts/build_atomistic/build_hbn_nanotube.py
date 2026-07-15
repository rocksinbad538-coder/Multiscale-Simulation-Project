#!/usr/bin/env python3
import math, csv
from pathlib import Path
from collections import Counter, defaultdict

OUT = Path("structures/phase1A")
OUT.mkdir(parents=True, exist_ok=True)

# -----------------------------
# Parameters
# -----------------------------
bn = 1.45             # target B-N bond length, Angstrom
target_diam = 24.0    # ~2.4 nm
target_len  = 60.0    # ~6.0 nm
passivate = False

bh = 1.20
nh = 1.02

# 2D honeycomb BN lattice:
# B at R = i*a1 + j*a2
# N at R + d1
sqrt3 = math.sqrt(3.0)
a1 = (sqrt3*bn, 0.0)
a2 = (0.5*sqrt3*bn, 1.5*bn)

dvecs = [
    (0.0, bn),
    (0.5*sqrt3*bn, -0.5*bn),
    (-0.5*sqrt3*bn, -0.5*bn),
]

# choose circumference commensurate with a1
m = max(6, round(math.pi * target_diam / (sqrt3 * bn)))
C = m * sqrt3 * bn
R = C / (2.0 * math.pi)

nrows = max(4, round(target_len / (1.5 * bn)))

def wrap2tube(x, y):
    x = x % C
    theta = 2.0 * math.pi * x / C
    return (R*math.cos(theta), R*math.sin(theta), y)

atoms = []
sites2d = []

# generate rectangular BN sheet and roll x direction
for j in range(nrows):
    for i in range(m):
        bx = i*a1[0] + j*a2[0]
        by = i*a1[1] + j*a2[1]
        nx = bx + dvecs[0][0]
        ny = by + dvecs[0][1]

        X,Y,Z = wrap2tube(bx, by)
        atoms.append({"elem":"B","x":X,"y":Y,"z":Z})
        sites2d.append(("B", bx % C, by))

        X,Y,Z = wrap2tube(nx, ny)
        atoms.append({"elem":"N","x":X,"y":Y,"z":Z})
        sites2d.append(("N", nx % C, ny))

def dist(a,b):
    return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)

def bond_list(atoms):
    bonds = []
    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            ei, ej = atoms[i]["elem"], atoms[j]["elem"]
            d = dist(atoms[i], atoms[j])
            if {ei,ej} == {"B","N"} and 1.20 <= d <= 1.70:
                bonds.append((i,j,d,"B-N"))
            elif {ei,ej} == {"B","H"} and 0.90 <= d <= 1.35:
                bonds.append((i,j,d,"B-H"))
            elif {ei,ej} == {"N","H"} and 0.80 <= d <= 1.20:
                bonds.append((i,j,d,"N-H"))
    return bonds

bonds = bond_list(atoms)
deg = Counter()
for i,j,d,t in bonds:
    deg[i]+=1
    deg[j]+=1

# -----------------------------
# Passivation
# -----------------------------
if passivate:
    h_atoms = []

    # use actual local radial/tangential/axial geometry from missing 2D neighbor directions
    existing_positions = [(round(x,6), round(y,6), elem) for elem,x,y in sites2d]

    def has_N_at(x,y):
        x = x % C
        for elem, xx, yy in sites2d:
            if elem == "N":
                dx = min(abs(xx-x), C-abs(xx-x))
                if dx < 1e-4 and abs(yy-y) < 1e-4:
                    return True
        return False

    def has_B_at(x,y):
        x = x % C
        for elem, xx, yy in sites2d:
            if elem == "B":
                dx = min(abs(xx-x), C-abs(xx-x))
                if dx < 1e-4 and abs(yy-y) < 1e-4:
                    return True
        return False

    nat0 = len(atoms)

    for idx in range(nat0):
        elem, x2, y2 = sites2d[idx]
        if elem not in ("B","N"):
            continue

        missing_dirs = []

        if elem == "B":
            for dx,dy in dvecs:
                if not has_N_at(x2+dx, y2+dy):
                    missing_dirs.append((dx,dy))
            hlen = bh
        else:
            for dx,dy in dvecs:
                if not has_B_at(x2-dx, y2-dy):
                    missing_dirs.append((-dx,-dy))
            hlen = nh

        for dx,dy in missing_dirs:
            norm = math.sqrt(dx*dx + dy*dy)
            ux, uy = dx/norm, dy/norm
            hx2 = x2 + ux*hlen
            hy2 = y2 + uy*hlen
            X,Y,Z = wrap2tube(hx2, hy2)
            h_atoms.append({"elem":"H","x":X,"y":Y,"z":Z})

    atoms.extend(h_atoms)

bonds = bond_list(atoms)
deg = Counter()
for i,j,d,t in bonds:
    deg[i]+=1
    deg[j]+=1

# -----------------------------
# Audits
# -----------------------------
overlaps = []
for i in range(len(atoms)):
    for j in range(i+1, len(atoms)):
        d = dist(atoms[i], atoms[j])
        ei, ej = atoms[i]["elem"], atoms[j]["elem"]

        # nonbonded hard threshold
        bonded = False
        for bi,bj,bd,bt in bonds:
            if (bi == i and bj == j) or (bi == j and bj == i):
                bonded = True
                break

        if not bonded and d < 0.90:
            overlaps.append((i+1,j+1,ei,ej,d))

bad_bn = [
    (i+1,a["elem"],deg[i]) for i,a in enumerate(atoms)
    if a["elem"] in ("B","N") and deg[i] != 3
]

vals = defaultdict(list)
for _,_,d,t in bonds:
    vals[t].append(d)

# XYZ
with (OUT/"hbn_nt_initial.xyz").open("w") as f:
    f.write(f"{len(atoms)}\n")
    f.write(f"h-BN nanotube armchair-like; m={m}; R={R:.6f} A; BN={bn:.3f} A\n")
    for a in atoms:
        f.write(f'{a["elem"]:2s} {a["x"]: .8f} {a["y"]: .8f} {a["z"]: .8f}\n')

# PDB
with (OUT/"hbn_nt_initial.pdb").open("w") as f:
    f.write("REMARK h-BN nanotube atomistic model, Phase 1A\n")
    for k,a in enumerate(atoms,1):
        f.write(
            f"HETATM{k:5d} {a['elem']:>2s}   HBN A   1    "
            f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
            f"  1.00  0.00          {a['elem']:>2s}\n"
        )
    f.write("END\n")

with (OUT/"audit_bonds.txt").open("w") as f:
    for i,j,d,t in bonds:
        f.write(f"{i+1:6d} {j+1:6d} {t:4s} {d:10.5f}\n")

with (OUT/"audit_connectivity.txt").open("w") as f:
    for i,a in enumerate(atoms):
        f.write(f"{i+1:6d} {a['elem']:2s} degree={deg[i]}\n")

with (OUT/"audit_distances.txt").open("w") as f:
    for t in ["B-N","B-H","N-H"]:
        ds = vals.get(t,[])
        if ds:
            f.write(f"{t}: count={len(ds)} min={min(ds):.5f} mean={sum(ds)/len(ds):.5f} max={max(ds):.5f}\n")
        else:
            f.write(f"{t}: count=0\n")

with (OUT/"audit_overlaps.txt").open("w") as f:
    if overlaps:
        for row in overlaps:
            f.write("%6d %6d %2s %2s %10.5f\n" % row)
    else:
        f.write("No nonbonded overlaps below 0.90 A\n")

counts = Counter(a["elem"] for a in atoms)

with (OUT/"audit_summary.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["m_circumference_cells",m])
    w.writerow(["rows",nrows])
    w.writerow(["radius_A",f"{R:.6f}"])
    w.writerow(["diameter_A",f"{2*R:.6f}"])
    w.writerow(["atoms_total",len(atoms)])
    w.writerow(["B_atoms",counts["B"]])
    w.writerow(["N_atoms",counts["N"]])
    w.writerow(["H_atoms",counts["H"]])
    w.writerow(["bonds_total",len(bonds)])
    w.writerow(["BN_bonds",len(vals.get("B-N",[]))])
    w.writerow(["BH_bonds",len(vals.get("B-H",[]))])
    w.writerow(["NH_bonds",len(vals.get("N-H",[]))])
    w.writerow(["overlap_count",len(overlaps)])
    w.writerow(["non_degree3_BN_after_passivation",len(bad_bn)])

print("Wrote corrected Phase 1A h-BN nanotube files to", OUT)
print("m =", m, "rows =", nrows, "diameter_A =", 2*R, "atoms =", len(atoms))
print("overlaps =", len(overlaps))
print("non_degree3_BN_after_passivation =", len(bad_bn))
