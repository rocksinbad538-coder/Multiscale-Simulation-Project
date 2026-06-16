#!/usr/bin/env python3
from pathlib import Path
import csv, re

mol2 = Path("parameters/phase1A/pyrene_gaff2/pyrene_gaff2_neutral.mol2")
frcmod = Path("parameters/phase1A/pyrene_gaff2/pyrene_gaff2.frcmod")
tleap_log = Path("logs/phase1A/tleap_pyrene_gaff2_neutral.log")
prmtop = Path("parameters/phase1A/pyrene_gaff2/pyrene_neutral.prmtop")
inpcrd = Path("parameters/phase1A/pyrene_gaff2/pyrene_neutral.inpcrd")
out = Path("parameters/phase1A/pyrene_gaff2/audit_pyrene_gaff2_neutral.csv")

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
            atoms.append((int(parts[0]), parts[1], parts[5], float(parts[8])))
            charges.append(float(parts[8]))

tleap_text = tleap_log.read_text(errors="ignore") if tleap_log.exists() else ""
frcmod_text = frcmod.read_text(errors="ignore") if frcmod.exists() else ""

has_tleap_errors = bool(re.search(r"\bERROR\b|Fatal|Could not", tleap_text, re.I))
has_missing_params = bool(re.search(r"ATTN|MISSING|need revision", frcmod_text, re.I))

with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric","value"])
    w.writerow(["neutral_mol2_exists",mol2.exists()])
    w.writerow(["frcmod_exists",frcmod.exists()])
    w.writerow(["prmtop_exists",prmtop.exists()])
    w.writerow(["inpcrd_exists",inpcrd.exists()])
    w.writerow(["atom_count",len(atoms)])
    w.writerow(["unique_atom_types",";".join(sorted(set(a[2] for a in atoms)))])
    w.writerow(["total_charge_e",f"{sum(charges):.8f}"])
    w.writerow(["abs_total_charge_e",f"{abs(sum(charges)):.8e}"])
    w.writerow(["tleap_errors_detected",has_tleap_errors])
    w.writerow(["frcmod_missing_or_attention_detected",has_missing_params])

print(out.read_text())
