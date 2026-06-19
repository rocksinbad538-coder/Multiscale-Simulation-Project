from pathlib import Path
import csv

root = Path("runs/phase1A/day012_electronic_prep/orca_inputs")
out = root / "orca_input_audit.csv"

valid_elements = {"H", "C", "N", "O", "B"}
bad_elements = {"M", "MW", "X"}

rows = []

for inp in sorted(root.glob("**/*.inp")):
    if inp.name.endswith(".inp"):
        natoms = 0
        bad = 0
        has_charge_mult = False
        has_tddft = False
        has_pbe0 = False

        for line in inp.read_text().splitlines():
            s = line.strip()
            parts = s.split()

            if s.lower().startswith("* xyz"):
                has_charge_mult = True

            if "%tddft" in s.lower():
                has_tddft = True

            if "PBE0" in s:
                has_pbe0 = True

            if len(parts) == 4:
                elem = parts[0]
                if elem in valid_elements:
                    natoms += 1
                if elem in bad_elements:
                    bad += 1

        rows.append({
            "input_file": str(inp),
            "qm_atoms": natoms,
            "bad_virtual_sites": bad,
            "has_charge_multiplicity": has_charge_mult,
            "has_tddft_block": has_tddft,
            "has_pbe0": has_pbe0,
        })

with out.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=[
        "input_file",
        "qm_atoms",
        "bad_virtual_sites",
        "has_charge_multiplicity",
        "has_tddft_block",
        "has_pbe0",
    ])
    w.writeheader()
    w.writerows(rows)

print(out.read_text())
