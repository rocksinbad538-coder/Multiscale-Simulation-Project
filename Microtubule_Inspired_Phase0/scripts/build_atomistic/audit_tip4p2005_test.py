#!/usr/bin/env python3
from pathlib import Path
import csv, re

run = Path("runs/phase1A/tip4p2005_test")
log = Path("logs/phase1A/grompp_tip4p2005_test.log")
text = log.read_text(errors="ignore") if log.exists() else ""

fatal = bool(re.search(r"Fatal error|ERROR|Segmentation fault", text, re.I))
warning_count = len(re.findall(r"WARNING \d+", text))
note_count = len(re.findall(r"NOTE \d+", text))

out = run/"audit_tip4p2005_test.csv"
with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["tip4p2005_itp_exists", Path("parameters/phase1A/water_tip4p2005/tip4p2005.itp").exists()])
    w.writerow(["tip4p2005_atomtypes_exists", Path("parameters/phase1A/water_tip4p2005/tip4p2005_atomtypes.itp").exists()])
    w.writerow(["tip4p2005_gro_exists", Path("parameters/phase1A/water_tip4p2005/tip4p2005.gro").exists()])
    w.writerow(["tpr_exists", (run/"min.tpr").exists()])
    w.writerow(["fatal_or_error_detected", fatal])
    w.writerow(["warning_count", warning_count])
    w.writerow(["note_count", note_count])

print(out.read_text())
