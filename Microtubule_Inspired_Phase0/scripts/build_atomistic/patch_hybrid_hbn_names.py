#!/usr/bin/env python3
from pathlib import Path

gro = Path("parameters/phase1A/hybrid_dry_gromacs/hbn_pyrene_4_dry.gro")
itp = Path("parameters/phase1A/hybrid_dry_gromacs/hbn_fixed_dummy.itp")

# Patch GRO atom names for first 1680 HBN atoms from B/N to B0001/N0002...
lines = gro.read_text().splitlines()
title = lines[0]
natoms = int(lines[1].strip())
atom_lines = lines[2:2+natoms]
box = lines[2+natoms]

new_atom_lines = []
for i, line in enumerate(atom_lines, start=1):
    if i <= 1680:
        elem = "B" if i % 2 == 1 else "N"
        name = f"{elem}{i:04d}"[-5:]
        # GROMACS .gro fixed columns:
        # resnr 1-5, resname 6-10, atomname 11-15, atomnr 16-20, coords
        new_line = line[:10] + f"{name:>5s}" + line[15:]
        new_atom_lines.append(new_line)
    else:
        new_atom_lines.append(line)

gro.write_text("\n".join([title, f"{natoms:5d}"] + new_atom_lines + [box]) + "\n")

print("Patched HBN atom names in", gro)
