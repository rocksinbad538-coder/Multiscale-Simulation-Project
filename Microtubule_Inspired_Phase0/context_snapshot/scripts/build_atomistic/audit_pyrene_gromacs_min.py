#!/usr/bin/env python3
from pathlib import Path
import csv, re

run = Path("runs/phase1A/pyrene_gromacs_min")
logs = [
    Path("logs/phase1A/grompp_pyrene_gromacs_min_pbc.log"),
    Path("logs/phase1A/mdrun_pyrene_gromacs_min_pbc.log"),
    run/"min.log",
]

files = {
    "tpr_exists": run/"min.tpr",
    "gro_exists": run/"min.gro",
    "edr_exists": run/"min.edr",
    "log_exists": run/"min.log",
}

text = "\n".join(p.read_text(errors="ignore") for p in logs if p.exists())

fatal_error = bool(re.search(r"Fatal error|Segmentation fault", text, re.I))
nan_detected = bool(re.search(r"(^|[\s=,:])[-+]?nan($|[\s,;])", text, re.I))

converged = bool(re.search(r"Steepest Descents converged to Fmax", text, re.I))
pot_match = re.findall(r"Potential Energy\s+=\s+([-\d.Ee+]+)", text)
fmax_match = re.findall(r"Maximum force\s+=\s+([-\d.Ee+]+)", text)

fatal = fatal_error or nan_detected

out = run/"audit_pyrene_gromacs_min.csv"
with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    for k,p in files.items():
        w.writerow([k, p.exists()])
    w.writerow(["fatal_error_detected", fatal_error])
    w.writerow(["nan_detected", nan_detected])
    w.writerow(["fatal_or_nan_detected", fatal])
    w.writerow(["minimization_converged_detected", converged])
    w.writerow(["final_potential_energy_kJmol", pot_match[-1] if pot_match else "NA"])
    w.writerow(["final_max_force_kJmol_nm", fmax_match[-1] if fmax_match else "NA"])

print(out.read_text())
