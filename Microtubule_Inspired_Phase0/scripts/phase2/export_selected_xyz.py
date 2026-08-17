#!/usr/bin/env python3

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

SOURCE = ROOT/"runs"/"phase2"/"campaign"/"tddft_selected_structures"

OUT = ROOT/"runs"/"phase2"/"campaign"/"tddft_xyz"

OUT.mkdir(exist_ok=True)

TYPE_MAP = {
    1: "B",
    2: "N",
    3: "H",
    4: "O",
}

for dump in sorted(SOURCE.glob("*.dump")):

    with open(dump) as f:
        lines = f.readlines()

    atom_start = None

    for i, line in enumerate(lines):
        if line.startswith("ITEM: ATOMS"):
            atom_start = i + 1
            break

    atoms = []

    for line in lines[atom_start:]:
        s = line.split()

        if len(s) < 7:
            continue

        typ = int(s[2])

        x = float(s[4])
        y = float(s[5])
        z = float(s[6])

        atoms.append(
            (
                TYPE_MAP.get(typ, "X"),
                x,
                y,
                z
            )
        )

    xyz = OUT/(dump.stem + ".xyz")

    with open(xyz, "w") as f:

        f.write(f"{len(atoms)}\n")
        f.write(f"{dump.stem}\n")

        for e, x, y, z in atoms:

            f.write(
                f"{e:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n"
            )

print("="*90)
print("PHASE3-A05")
print("XYZ EXPORT COMPLETED")
print("="*90)
print(OUT)
