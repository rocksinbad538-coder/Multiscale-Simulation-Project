#!/usr/bin/env python3
from pathlib import Path
import csv, re

run = Path("runs/phase1A/hybrid_hydrated_gromacs_min")
logs = [
    Path("logs/phase1A/grompp_hybrid_hydrated_min_gap45_clean032_ext.log"),
    Path("logs/phase1A/mdrun_hybrid_hydrated_min_gap45_clean032_ext.log"),
    run/"min_gap45_clean032_ext.log",
]

text = "\n".join(p.read_text(errors="ignore") for p in logs if p.exists())

fatal = bool(re.search(r"Fatal error|Segmentation fault", text, re.I))
nan = bool(re.search(r"(^|[\s=,:])[-+]?nan($|[\s,;])", text, re.I))
converged = bool(re.search(r"Steepest Descents converged to Fmax", text, re.I))

pot = re.findall(r"Potential Energy\s+=\s+([-\d.Ee+]+)", text)
fmax = re.findall(r"Maximum force\s+=\s+([-\d.Ee+]+)", text)
atom = re.findall(r"Maximum force\s+=\s+[-\d.Ee+]+\s+on atom\s+(\d+)", text)

out = run/"audit_hybrid_hydrated_min_gap45_clean032_ext.csv"
with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["tpr_exists",(run/"min_gap45_clean032_ext.tpr").exists()])
    w.writerow(["gro_exists",(run/"min_gap45_clean032_ext.gro").exists()])
    w.writerow(["edr_exists",(run/"min_gap45_clean032_ext.edr").exists()])
    w.writerow(["log_exists",(run/"min_gap45_clean032_ext.log").exists()])
    w.writerow(["fatal_error_detected",fatal])
    w.writerow(["nan_detected",nan])
    w.writerow(["fatal_or_nan_detected",fatal or nan])
    w.writerow(["minimization_converged_detected",converged])
    w.writerow(["final_potential_energy_kJmol",pot[-1] if pot else "NA"])
    w.writerow(["final_max_force_kJmol_nm",fmax[-1] if fmax else "NA"])
    w.writerow(["final_max_force_atom",atom[-1] if atom else "NA"])

print(out.read_text())
