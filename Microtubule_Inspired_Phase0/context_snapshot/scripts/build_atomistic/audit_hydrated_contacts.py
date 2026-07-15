#!/usr/bin/env python3
from pathlib import Path
import csv, math

gro = Path("parameters/phase1A/hybrid_hydrated_gromacs/hbn_pyrene_4_tip4p2005_solvated.gro")
out = Path("parameters/phase1A/hybrid_hydrated_gromacs/audit_hydrated_contacts.csv")

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

solute = [a for a in atoms if a["resname"] in ("HBN", "PYR")]
pyr = [a for a in atoms if a["resname"] == "PYR"]
hbn = [a for a in atoms if a["resname"] == "HBN"]
waters = [a for a in atoms if a["resname"] == "SOL"]

def dist(a,b):
    return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)

def min_pair(group1, group2):
    best = (999.0, None, None)
    for a in group1:
        for b in group2:
            d = dist(a,b)
            if d < best[0]:
                best = (d,a,b)
    return best

water_ow = [a for a in waters if a["name"] == "OW"]
water_all = waters

pairs = []
for label, g1, g2 in [
    ("OW_PYR", water_ow, pyr),
    ("WATERATOM_PYR", water_all, pyr),
    ("OW_HBN", water_ow, hbn),
    ("WATERATOM_HBN", water_all, hbn),
]:
    d,a,b = min_pair(g1,g2)
    pairs.append((label,d,a,b))

# neighbors near max-force atom from previous run: 1773
target = next((a for a in atoms if a["idx"] == 1773), None)
neighbors = []
if target:
    for a in atoms:
        if a["idx"] == target["idx"]:
            continue
        d = dist(target,a)
        if d < 0.35:
            neighbors.append((d,a))

with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["natoms", natoms])
    w.writerow(["hbn_atoms", len(hbn)])
    w.writerow(["pyrene_atoms", len(pyr)])
    w.writerow(["water_atoms", len(waters)])
    w.writerow(["water_molecules", len(waters)//4])
    for label,d,a,b in pairs:
        w.writerow([f"min_{label}_nm", f"{d:.5f}"])
        w.writerow([f"min_{label}_atom_A", f"{a['idx']}:{a['resname']}:{a['name']}"])
        w.writerow([f"min_{label}_atom_B", f"{b['idx']}:{b['resname']}:{b['name']}"])
    if target:
        w.writerow(["target_atom_1773", f"{target['idx']}:{target['resname']}:{target['name']}"])
        w.writerow(["neighbors_within_0p35nm_count", len(neighbors)])
        for d,a in sorted(neighbors)[:20]:
            w.writerow(["neighbor_1773", f"{d:.5f} nm -> {a['idx']}:{a['resname']}:{a['name']}"])

print(out.read_text())
