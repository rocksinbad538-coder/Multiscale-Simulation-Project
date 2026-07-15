#!/usr/bin/env python3

from pathlib import Path
import re

inp = Path("systems/inorganic/bn_like_scaffold_water/outputs/bn_like_scaffold_water_30000w_nvt300_hold_extended_contained.data")
out = Path("systems/inorganic/carbon_like_scaffold_water/outputs/carbon_like_scaffold_water_30000w_initial_from_bn_reference.data")

text = inp.read_text().splitlines()

new = []
section = None

for line in text:
    s = line.strip()

    if s == "Masses":
        section = "Masses"
        new.append(line)
        continue
    if s == "Atoms # full" or s == "Atoms":
        section = "Atoms"
        new.append(line)
        continue
    if s in {"Velocities", "Bonds", "Angles", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs"}:
        section = s
        new.append(line)
        continue

    if section == "Masses":
        parts = s.split()
        if len(parts) >= 2 and parts[0].isdigit():
            typ = int(parts[0])
            if typ == 1:
                new.append("1 12.011 # carbon-like neutral scaffold")
                continue
            elif typ == 2:
                new.append("2 12.011 # unused scaffold type placeholder")
                continue
        new.append(line)
        continue

    if section == "Atoms":
        parts = s.split()
        if len(parts) >= 7 and parts[0].isdigit():
            atom_id = parts[0]
            mol_id = parts[1]
            atom_type = int(parts[2])
            charge = parts[3]

            # Convert scaffold types 1 and 2 to neutral carbon-like type 1.
            if atom_type in (1, 2):
                parts[2] = "1"
                parts[3] = "0.000000"
                new.append(" ".join(parts))
                continue

        new.append(line)
        continue

    new.append(line)

out.write_text("\n".join(new) + "\n")
print(f"Wrote {out}")
