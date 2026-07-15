#!/usr/bin/env python3
from pathlib import Path
import math, csv, re

inp = Path("parameters/phase1A/hybrid_hydrated_gromacs/hbn_pyrene_4_tip4p2005_solvated_gap45.gro")
top_in = Path("parameters/phase1A/hybrid_hydrated_gromacs/hbn_pyrene_4_hydratable.top")

outgro = Path("parameters/phase1A/hybrid_hydrated_gromacs/hbn_pyrene_4_tip4p2005_solvated_gap45_clean032.gro")
outtop = Path("parameters/phase1A/hybrid_hydrated_gromacs/hbn_pyrene_4_hydratable_gap45_clean032.top")
audit = Path("parameters/phase1A/hybrid_hydrated_gromacs/audit_hydration_gap45_clean032.csv")

cutoff_nm = 0.32

lines = inp.read_text().splitlines()
title = lines[0]
natoms = int(lines[1].strip())
atom_lines = lines[2:2+natoms]
box = lines[2+natoms]

atoms = []
for line in atom_lines:
    atoms.append({
        "line": line,
        "idx": int(line[15:20]),
        "resnr": int(line[0:5]),
        "resname": line[5:10].strip(),
        "name": line[10:15].strip(),
        "x": float(line[20:28]),
        "y": float(line[28:36]),
        "z": float(line[36:44]),
    })

solute = [a for a in atoms if a["resname"] in ("HBN", "PYR")]
waters = [a for a in atoms if a["resname"] == "SOL"]

waters_by_resnr = {}
for a in waters:
    waters_by_resnr.setdefault(a["resnr"], []).append(a)

def dist(a,b):
    return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)

remove_resnrs = set()
min_contact = 999.0
min_pair = None

for w in waters:
    for s in solute:
        d = dist(w,s)
        if d < min_contact:
            min_contact = d
            min_pair = (w,s)
        if d < cutoff_nm:
            remove_resnrs.add(w["resnr"])
            break

kept_atoms = []
for a in atoms:
    if a["resname"] == "SOL" and a["resnr"] in remove_resnrs:
        continue
    kept_atoms.append(a)

kept_sol_mols = len([a for a in kept_atoms if a["resname"] == "SOL"]) // 4

# Renumber atom serials only; keep residue IDs readable
new_lines = []
for new_idx, a in enumerate(kept_atoms, start=1):
    line = a["line"]
    new_line = line[:15] + f"{new_idx:5d}" + line[20:]
    new_lines.append(new_line)

outgro.write_text("\n".join([title + " CLEANED", f"{len(kept_atoms):5d}"] + new_lines + [box]) + "\n")

# Update topology SOL count
top_text = top_in.read_text()
new_top_lines = []
for line in top_text.splitlines():
    parts = line.split()
    if len(parts) == 2 and parts[0] == "SOL":
        new_top_lines.append(f"SOL {kept_sol_mols}")
    else:
        new_top_lines.append(line)
outtop.write_text("\n".join(new_top_lines) + "\n")

# Recompute min contacts after cleaning
kept_solute = [a for a in kept_atoms if a["resname"] in ("HBN", "PYR")]
kept_water = [a for a in kept_atoms if a["resname"] == "SOL"]

post_min = 999.0
post_pair = None
for w in kept_water:
    for s in kept_solute:
        d = dist(w,s)
        if d < post_min:
            post_min = d
            post_pair = (w,s)

with audit.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["input_atoms", natoms])
    w.writerow(["output_atoms", len(kept_atoms)])
    w.writerow(["input_water_molecules", len(waters)//4])
    w.writerow(["removed_water_molecules", len(remove_resnrs)])
    w.writerow(["kept_water_molecules", kept_sol_mols])
    w.writerow(["cutoff_nm", cutoff_nm])
    w.writerow(["pre_min_water_solute_nm", f"{min_contact:.5f}"])
    w.writerow(["post_min_water_solute_nm", f"{post_min:.5f}"])
    if min_pair:
        wa, sa = min_pair
        w.writerow(["pre_min_pair", f"{wa['idx']}:{wa['resname']}:{wa['name']} - {sa['idx']}:{sa['resname']}:{sa['name']}"])
    if post_pair:
        wa, sa = post_pair
        w.writerow(["post_min_pair", f"{wa['idx']}:{wa['resname']}:{wa['name']} - {sa['idx']}:{sa['resname']}:{sa['name']}"])

print(audit.read_text())
