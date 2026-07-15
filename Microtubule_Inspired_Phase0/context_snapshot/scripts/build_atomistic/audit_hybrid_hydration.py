#!/usr/bin/env python3
from pathlib import Path
import csv, re, math

base = Path("parameters/phase1A/hybrid_hydrated_gromacs")
gro = base/"hbn_pyrene_4_tip4p2005_solvated.gro"
top = base/"hbn_pyrene_4_hydratable.top"
log = Path("logs/phase1A/solvate_hybrid_tip4p2005.log")
out = base/"audit_hybrid_hydration.csv"

text = log.read_text(errors="ignore") if log.exists() else ""
fatal = bool(re.search(r"Fatal error|ERROR|Segmentation fault", text, re.I))

natoms = "NA"
n_hbn = 0
n_pyr = 0
n_sol_atoms = 0
box = ("NA","NA","NA")

if gro.exists():
    lines = gro.read_text().splitlines()
    natoms = int(lines[1].strip())
    atom_lines = lines[2:2+natoms]
    box = tuple(lines[2+natoms].split()[:3])
    for line in atom_lines:
        resname = line[5:10].strip()
        if resname == "HBN":
            n_hbn += 1
        elif resname == "PYR":
            n_pyr += 1
        elif resname == "SOL":
            n_sol_atoms += 1

n_sol = n_sol_atoms // 4 if isinstance(n_sol_atoms, int) else "NA"
expected_atoms = 1680 + 104 + (4*n_sol if isinstance(n_sol, int) else 0)

# Parse SOL count from topology
top_text = top.read_text(errors="ignore") if top.exists() else ""
sol_count_top = "NA"
for line in top_text.splitlines():
    parts = line.split()
    if len(parts) == 2 and parts[0] == "SOL":
        sol_count_top = int(parts[1])

with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["solvated_gro_exists", gro.exists()])
    w.writerow(["top_exists", top.exists()])
    w.writerow(["solvate_fatal_or_error_detected", fatal])
    w.writerow(["natoms_total", natoms])
    w.writerow(["hbn_atoms", n_hbn])
    w.writerow(["pyrene_atoms_total", n_pyr])
    w.writerow(["sol_atoms", n_sol_atoms])
    w.writerow(["sol_molecules_from_gro", n_sol])
    w.writerow(["sol_molecules_from_top", sol_count_top])
    w.writerow(["expected_atoms_from_counts", expected_atoms])
    w.writerow(["count_consistent", expected_atoms == natoms])
    w.writerow(["box_x_nm", box[0]])
    w.writerow(["box_y_nm", box[1]])
    w.writerow(["box_z_nm", box[2]])

print(out.read_text())
