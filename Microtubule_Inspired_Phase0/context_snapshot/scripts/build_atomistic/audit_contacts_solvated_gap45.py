#!/usr/bin/env python3
from pathlib import Path
import csv, math

gro = Path("parameters/phase1A/hybrid_hydrated_gromacs/hbn_pyrene_4_tip4p2005_solvated_gap45.gro")
out = Path("parameters/phase1A/hybrid_hydrated_gromacs/audit_contacts_solvated_gap45.csv")

lines = gro.read_text().splitlines()
natoms = int(lines[1].strip())
box = [float(x) for x in lines[2 + natoms].split()[:3]]

atoms = []
for line in lines[2:2 + natoms]:
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

def dist(a, b):
    return math.sqrt(
        pbc_delta(a["x"] - b["x"], box[0])**2 +
        pbc_delta(a["y"] - b["y"], box[1])**2 +
        pbc_delta(a["z"] - b["z"], box[2])**2
    )

hbn = [a for a in atoms if a["resname"] == "HBN"]
pyr = [a for a in atoms if a["resname"] == "PYR"]
water = [a for a in atoms if a["resname"] == "SOL"]
ow = [a for a in water if a["name"] == "OW"]
solute = hbn + pyr

def minpair(A, B):
    best = (999.0, None, None)
    for a in A:
        for b in B:
            d = dist(a, b)
            if d < best[0]:
                best = (d, a, b)
    return best

items = {
    "min_OW_solute": minpair(ow, solute),
    "min_WATERATOM_solute": minpair(water, solute),
    "min_PYR_HBN": minpair(pyr, hbn),
}

with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric", "value"])
    w.writerow(["source_gro", str(gro)])
    w.writerow(["natoms", natoms])
    w.writerow(["water_molecules", len(water)//4])
    w.writerow(["box_nm", ";".join(str(x) for x in box)])
    for label, (d, a, b) in items.items():
        w.writerow([label + "_nm", f"{d:.5f}"])
        w.writerow([label + "_pair", f"{a['idx']}:{a['resname']}:{a['name']} - {b['idx']}:{b['resname']}:{b['name']}"])

print(out.read_text())
