#!/usr/bin/env python3
from pathlib import Path
import csv, re

run = Path("runs/phase1A/hybrid_hydrated_gromacs_min")
log = Path("logs/phase1A/grompp_hybrid_hydrated_min.log")
text = log.read_text(errors="ignore") if log.exists() else ""

fatal = bool(re.search(r"Fatal error|ERROR|Segmentation fault", text, re.I))
warning_count = len(re.findall(r"WARNING \\d+", text))
note_count = len(re.findall(r"NOTE \\d+", text))

out = run/"audit_hybrid_hydrated_grompp.csv"
with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["tpr_exists", (run/"min.tpr").exists()])
    w.writerow(["fatal_or_error_detected", fatal])
    w.writerow(["warning_count", warning_count])
    w.writerow(["note_count", note_count])

print(out.read_text())
