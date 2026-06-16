#!/usr/bin/env python3
from pathlib import Path
import csv, re

run = Path("runs/phase1A/hybrid_hydrated_gromacs_min")
logs = [
    Path("logs/phase1A/grompp_hybrid_hydrated_min.log"),
    Path("logs/phase1A/mdrun_hybrid_hydrated_min.log"),
    run/"min.log",
]

text = "\n".join(p.read_text(errors="ignore") for p in logs if p.exists())

fatal_error = bool(re.search(r"Fatal error|Segmentation fault", text, re.I))
nan_detected = bool(re.search(r"(^|[\s=,:])[-+]?nan($|[\s,;])", text, re.I))
converged = bool(re.search(r"Steepest Descents converged to Fmax", text, re.I))

pot_match = re.findall(r"Potential Energy\s+=\s+([-\d.Ee+]+)", text)
fmax_match = re.findall(r"Maximum force\s+=\s+([-\d.Ee+]+)", text)

out = run/"audit_hybrid_hydrated_min.csv"
with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["tpr_exists", (run/"min.tpr").exists()])
    w.writerow(["gro_exists", (run/"min.gro").exists()])
    w.writerow(["edr_exists", (run/"min.edr").exists()])
    w.writerow(["log_exists", (run/"min.log").exists()])
    w.writerow(["fatal_error_detected", fatal_error])
    w.writerow(["nan_detected", nan_detected])
    w.writerow(["fatal_or_nan_detected", fatal_error or nan_detected])
    w.writerow(["minimization_converged_detected", converged])
    w.writerow(["final_potential_energy_kJmol", pot_match[-1] if pot_match else "NA"])
    w.writerow(["final_max_force_kJmol_nm", fmax_match[-1] if fmax_match else "NA"])

print(out.read_text())
