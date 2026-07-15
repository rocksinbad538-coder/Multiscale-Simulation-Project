#!/usr/bin/env python3

from pathlib import Path

inp = Path("systems/inorganic/bn_like_scaffold_water/outputs/bn_like_scaffold_water_30000w_nvt300_hold_extended_contained.data")
out = Path("systems/inorganic/bn_neutralized_scaffold_water/outputs/bn_neutralized_scaffold_water_30000w_initial_from_bn_reference.data")

lines = inp.read_text().splitlines()
new = []
section = None

for line in lines:
    s = line.strip()

    if s in {"Masses", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs", "Velocities", "Bonds", "Angles"}:
        section = s
        new.append(line)
        continue

    if s == "Atoms # full" or s == "Atoms":
        section = "Atoms"
        new.append(line)
        continue

    if section == "Atoms":
        parts = s.split()
        if len(parts) >= 7 and parts[0].isdigit():
            atom_type = int(parts[2])
            if atom_type in (1, 2):
                parts[3] = "0.000000"
                new.append(" ".join(parts))
                continue

    new.append(line)

out.write_text("\n".join(new) + "\n")
print(f"Wrote {out}")
