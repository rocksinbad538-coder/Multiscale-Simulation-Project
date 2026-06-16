#!/usr/bin/env python3
from pathlib import Path
import csv, re

base = Path("parameters/phase1A/hybrid_dry_gromacs")
run = Path("runs/phase1A/hybrid_dry_gromacs_min")
log = Path("logs/phase1A/grompp_hybrid_dry_gromacs_min.log")

text = log.read_text(errors="ignore") if log.exists() else ""

fatal = bool(re.search(r"Fatal error|ERROR|Segmentation fault", text, re.I))
note_count = len(re.findall(r"NOTE \d+", text))
warning_count = len(re.findall(r"WARNING \d+", text))

gro = base/"hbn_pyrene_4_dry.gro"
natoms_gro = "NA"
if gro.exists():
    lines = gro.read_text().splitlines()
    natoms_gro = int(lines[1].strip())

out = base/"audit_hybrid_dry_gromacs_topology.csv"
with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["gro_exists", gro.exists()])
    w.writerow(["top_exists", (base/"hbn_pyrene_4_dry.top").exists()])
    w.writerow(["tpr_exists", (run/"min.tpr").exists()])
    w.writerow(["natoms_gro", natoms_gro])
    w.writerow(["fatal_or_error_detected", fatal])
    w.writerow(["note_count", note_count])
    w.writerow(["warning_count", warning_count])

print(out.read_text())
