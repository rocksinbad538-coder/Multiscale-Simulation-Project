#!/usr/bin/env python3
from pathlib import Path
import re, csv

outf = Path("runs/phase1A/pyrene_gaff2_min/min_pyrene.out")
rst = Path("runs/phase1A/pyrene_gaff2_min/min_pyrene.rst")
pdb = Path("runs/phase1A/pyrene_gaff2_min/pyrene_minimized.pdb")
audit = Path("runs/phase1A/pyrene_gaff2_min/audit_pyrene_minimization.csv")

text = outf.read_text(errors="ignore") if outf.exists() else ""

completed = "Total CPU time" in text or "FINAL RESULTS" in text
has_error = bool(re.search(r"ERROR|LINMIN FAILURE|NaN|Inf", text, re.I))

energies = re.findall(r"ENERGY\s+=\s+([-\d.Ee+]+)", text)
final_energy = energies[-1] if energies else "NA"

with audit.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["min_output_exists",outf.exists()])
    w.writerow(["restart_exists",rst.exists()])
    w.writerow(["minimized_pdb_exists",pdb.exists()])
    w.writerow(["completed_detected",completed])
    w.writerow(["error_detected",has_error])
    w.writerow(["final_energy_raw",final_energy])

print(audit.read_text())
