#!/usr/bin/env python3
from pathlib import Path
import csv, re

mol2 = Path("parameters/phase1A/pyrene_gaff2/pyrene_gaff2.mol2")
frcmod = Path("parameters/phase1A/pyrene_gaff2/pyrene_gaff2.frcmod")
tleap_log = Path("logs/phase1A/tleap_pyrene_gaff2.log")
out = Path("parameters/phase1A/pyrene_gaff2/audit_pyrene_gaff2.csv")

if not mol2.exists():
    raise SystemExit("ERROR: pyrene_gaff2.mol2 not found")

atoms = []
charges = []
in_atom = False

for line in mol2.read_text().splitlines():
    if line.startswith("@<TRIPOS>ATOM"):
        in_atom = True
        continue
    if line.startswith("@<TRIPOS>BOND"):
        in_atom = False
        continue
    if in_atom and line.strip():
        parts = line.split()
        if len(parts) >= 9:
            idx = int(parts[0])
            name = parts[1]
            atom_type = parts[5]
            charge = float(parts[8])
            atoms.append((idx, name, atom_type, charge))
            charges.append(charge)

total_charge = sum(charges)
n_atoms = len(atoms)
types = sorted(set(a[2] for a in atoms))

tleap_text = tleap_log.read_text(errors="ignore") if tleap_log.exists() else ""
frcmod_text = frcmod.read_text(errors="ignore") if frcmod.exists() else ""

has_tleap_errors = bool(re.search(r"\bERROR\b|Fatal|Could not", tleap_text, re.I))
has_missing_params = bool(re.search(r"ATTN|MISSING|need revision", frcmod_text, re.I))

with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["mol2_exists",mol2.exists()])
    w.writerow(["frcmod_exists",frcmod.exists()])
    w.writerow(["tleap_log_exists",tleap_log.exists()])
    w.writerow(["atom_count",n_atoms])
    w.writerow(["unique_atom_types",";".join(types)])
    w.writerow(["total_charge_e",f"{total_charge:.8f}"])
    w.writerow(["abs_total_charge_e",f"{abs(total_charge):.8e}"])
    w.writerow(["tleap_errors_detected",has_tleap_errors])
    w.writerow(["frcmod_missing_or_attention_detected",has_missing_params])

print(out.read_text())
