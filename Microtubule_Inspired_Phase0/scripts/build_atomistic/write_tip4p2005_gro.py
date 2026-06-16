#!/usr/bin/env python3
from pathlib import Path

out = Path("parameters/phase1A/water_tip4p2005/tip4p2005.gro")

# Coordinates in nm. Geometry is rounded to .gro precision.
atoms = [
    (1, "SOL", "OW",  1,  1.500, 1.500, 1.500),
    (1, "SOL", "HW1", 2,  1.596, 1.500, 1.500),
    (1, "SOL", "HW2", 3,  1.476, 1.593, 1.500),
    (1, "SOL", "MW",  4,  1.515, 1.512, 1.500),
]

with out.open("w") as f:
    f.write("TIP4P/2005 single water in 3 nm box\n")
    f.write(f"{len(atoms):5d}\n")
    for resnr, resname, atomname, atomnr, x, y, z in atoms:
        f.write(f"{resnr:5d}{resname:<5s}{atomname:>5s}{atomnr:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
    f.write(f"{3.00000:10.5f}{3.00000:10.5f}{3.00000:10.5f}\n")

print(out.read_text())
