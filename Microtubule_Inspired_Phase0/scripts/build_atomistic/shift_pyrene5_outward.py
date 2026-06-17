#!/usr/bin/env python3
from pathlib import Path
import math

inp = Path("parameters/phase1A/hybrid_dry_gromacs/hbn_pyrene_4_dry.gro")
out = Path("parameters/phase1A/hybrid_dry_gromacs/hbn_pyrene_4_dry_pyr5_shifted.gro")

shift_nm = 0.12  # 1.2 Å outward only for pyrene residue 5

lines = inp.read_text().splitlines()
natoms = int(lines[1].strip())
atom_lines = lines[2:2+natoms]
box = lines[2+natoms]

new_lines = []
for line in atom_lines:
    resnr = int(line[0:5])
    resname = line[5:10].strip()
    x = float(line[20:28])
    y = float(line[28:36])
    z = float(line[36:44])

    if resname == "PYR" and resnr == 5:
        cx = x - 3.697265  # half of 7.39453
        cy = y - 3.697265
        r = math.sqrt(cx*cx + cy*cy)
        ux, uy = cx/r, cy/r
        x += shift_nm * ux
        y += shift_nm * uy

    new_lines.append(line[:20] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[44:])

out.write_text("\n".join([lines[0], f"{natoms:5d}"] + new_lines + [box]) + "\n")
print("Wrote", out)
