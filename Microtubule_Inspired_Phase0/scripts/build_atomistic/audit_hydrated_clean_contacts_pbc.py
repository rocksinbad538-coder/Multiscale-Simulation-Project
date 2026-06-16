#!/usr/bin/env python3
from pathlib import Path
import csv, math

gro = Path("runs/phase1A/hybrid_hydrated_gromacs_min/min_clean.gro")
if not gro.exists():
    gro = Path("parameters/phase1A/hybrid_hydrated_gromacs/hbn_pyrene_4_tip4p2005_solvated_clean.gro")

out = Path("parameters/phase1A/hybrid_hydrated_gromacs/audit_hydrated_clean_contacts_pbc.csv")

lines = gro.read_text().splitlines()
natoms = int(lines[1].strip())
box_vals = [float(x) for x in lines[2+natoms].split()[:3]]

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

def pbc_delta(dx, L):
    return dx - L * round(dx / L)

def dist(a,b):
    dx = pbc_delta(a["x"]-b["x"], box_vals[0])
    dy = pbc_delta(a["y"]-b["y"], box_vals[1])
    dz = pbc_delta(a["z"]-b["z"], box_vals[2])
    return math.sqrt(dx*dx + dy*dy + dz*dz)

solute = [a for a in atoms if a["resname"] in ("HBN","PYR")]
pyr = [a for a in atoms if a["resname"] == "PYR"]
hbn = [a for a in atoms if a["resname"] == "HBN"]
water = [a for a in atoms if a["resname"] == "SOL"]
ow = [a for a in water if a["name"] == "OW"]

def min_pair(g1,g2):
    best = (999,None,None)
    for a in g1:
        for b in g2:
            d = dist(a,b)
            if d < best[0]:
                best = (d,a,b)
    return best

target = next((a for a in atoms if a["idx"] == 1773), None)
neighbors = []
if target:
    for a in atoms:
        if a["idx"] != target["idx"]:
            d = dist(target,a)
            if d < 0.45:
                neighbors.append((d,a))

pairs = [
    ("OW_PYR", *min_pair(ow,pyr)),
    ("WATERATOM_PYR", *min_pair(water,pyr)),
    ("OW_HBN", *min_pair(ow,hbn)),
    ("WATERATOM_HBN", *min_pair(water,hbn)),
]

with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["source_gro", str(gro)])
    w.writerow(["natoms", natoms])
    w.writerow(["box_nm", ";".join(map(str,box_vals))])
    for label,d,a,b in pairs:
        w.writerow([f"min_{label}_nm", f"{d:.5f}"])
        w.writerow([f"min_{label}_pair", f"{a['idx']}:{a['resname']}:{a['name']} - {b['idx']}:{b['resname']}:{b['name']}"])
    if target:
        w.writerow(["target_atom_1773", f"{target['idx']}:{target['resname']}:{target['name']}"])
        w.writerow(["neighbors_1773_within_0p45nm_count", len(neighbors)])
        for d,a in sorted(neighbors)[:40]:
            w.writerow(["neighbor_1773", f"{d:.5f} nm -> {a['idx']}:{a['resname']}:{a['name']}"])

print(out.read_text())
