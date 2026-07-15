#!/usr/bin/env python3
from pathlib import Path
import csv, re

run = Path("runs/phase1A/hybrid_hydrated_gap45_clean032_nvt")
logs = [
    Path("logs/phase1A/grompp_nvt_gap45_clean032_100ps.log"),
    Path("logs/phase1A/mdrun_nvt_gap45_clean032_100ps.log"),
    run/"nvt_100ps.log",
]

text = "\n".join(p.read_text(errors="ignore") for p in logs if p.exists())

fatal = bool(re.search(r"Fatal error|Segmentation fault|LINCS WARNING", text, re.I))
nan = bool(re.search(r"(^|[\s=,:])[-+]?nan($|[\s,;])", text, re.I))
finished = "Finished mdrun" in text
temp = re.findall(r"Temperature\s*\n\s*[-\d.Ee+]+\s+[-\d.Ee+]+\s+[-\d.Ee+]+\s+[-\d.Ee+]+\s+([-\d.Ee+]+)", text)
drift = re.findall(r"Conserved energy drift:\s*([-\d.Ee+]+)", text)
perf = re.findall(r"Performance:\s*([-\d.Ee+]+)", text)

out = run/"audit_nvt_100ps.csv"
with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["tpr_exists",(run/"nvt_100ps.tpr").exists()])
    w.writerow(["gro_exists",(run/"nvt_100ps.gro").exists()])
    w.writerow(["xtc_exists",(run/"nvt_100ps.xtc").exists()])
    w.writerow(["edr_exists",(run/"nvt_100ps.edr").exists()])
    w.writerow(["log_exists",(run/"nvt_100ps.log").exists()])
    w.writerow(["finished_mdrun_detected",finished])
    w.writerow(["fatal_or_lincs_detected",fatal])
    w.writerow(["nan_detected",nan])
    w.writerow(["average_temperature_K",temp[-1] if temp else "NA"])
    w.writerow(["conserved_energy_drift_kJmol_ps_atom",drift[-1] if drift else "NA"])
    w.writerow(["performance_ns_per_day",perf[-1] if perf else "NA"])

print(out.read_text())
