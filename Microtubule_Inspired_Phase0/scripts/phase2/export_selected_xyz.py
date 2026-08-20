#!/usr/bin/env python3

from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]

SOURCE = (
    ROOT
    / "runs"
    / "phase2"
    / "campaign_phase5_corrected"
    / "tddft_selected_structures"
)

OUT = (
    ROOT
    / "runs"
    / "phase2"
    / "campaign_phase5_corrected"
    / "tddft_xyz"
)

TYPE_MAP = {
    1: "B",
    2: "N",
    3: "H",   # HB
    4: "H",   # HN
}

EXPECTED_NATOMS = 37


def read_single_dump(path):

    lines = path.read_text().splitlines()

    atom_header_index = None

    for i, line in enumerate(lines):
        if line.startswith("ITEM: ATOMS"):
            atom_header_index = i
            break

    if atom_header_index is None:
        raise RuntimeError(
            f"{path}: ITEM: ATOMS header not found."
        )

    columns = lines[atom_header_index].split()[2:]
    col = {
        name: i
        for i, name in enumerate(columns)
    }

    required = ["id", "type"]

    for name in required:
        if name not in col:
            raise RuntimeError(
                f"{path}: required column '{name}' missing."
            )

    if all(name in col for name in ("xu", "yu", "zu")):
        xyz_names = ("xu", "yu", "zu")
        coordinate_source = "unwrapped"
    elif all(name in col for name in ("x", "y", "z")):
        xyz_names = ("x", "y", "z")
        coordinate_source = "wrapped"
    else:
        raise RuntimeError(
            f"{path}: neither xu/yu/zu nor x/y/z available."
        )

    atoms = []

    for line in lines[atom_header_index + 1:]:

        if not line.strip():
            continue

        if line.startswith("ITEM:"):
            break

        fields = line.split()

        atom_id = int(fields[col["id"]])
        atom_type = int(fields[col["type"]])

        if atom_type not in TYPE_MAP:
            raise RuntimeError(
                f"{path}: unknown LAMMPS atom type {atom_type}."
            )

        x = float(fields[col[xyz_names[0]]])
        y = float(fields[col[xyz_names[1]]])
        z = float(fields[col[xyz_names[2]]])

        atoms.append(
            (
                atom_id,
                TYPE_MAP[atom_type],
                x,
                y,
                z,
            )
        )

    atoms.sort(key=lambda row: row[0])

    if len(atoms) != EXPECTED_NATOMS:
        raise RuntimeError(
            f"{path}: expected {EXPECTED_NATOMS} atoms, "
            f"parsed {len(atoms)}."
        )

    ids = [row[0] for row in atoms]

    if ids != list(range(1, EXPECTED_NATOMS + 1)):
        raise RuntimeError(
            f"{path}: unexpected atom-ID sequence."
        )

    return atoms, coordinate_source


if OUT.exists():
    shutil.rmtree(OUT)

OUT.mkdir(parents=True, exist_ok=True)

dump_files = sorted(
    SOURCE.glob("*.dump")
)

if not dump_files:
    raise RuntimeError(
        f"No selected TDDFT dumps found in {SOURCE}"
    )

for dump in dump_files:

    atoms, coordinate_source = read_single_dump(
        dump
    )

    xyz = OUT / f"{dump.stem}.xyz"

    with xyz.open("w") as f:

        f.write(f"{len(atoms)}\n")

        f.write(
            f"{dump.stem} "
            f"coordinates={coordinate_source}\n"
        )

        for _, element, x, y, z in atoms:

            f.write(
                f"{element:2s} "
                f"{x:15.8f} "
                f"{y:15.8f} "
                f"{z:15.8f}\n"
            )

    print(
        f"{dump.name} -> {xyz.name} "
        f"atoms={len(atoms)} "
        f"coordinates={coordinate_source}"
    )

print()
print("="*90)
print("PHASE5-D50")
print("TDDFT XYZ EXPORT")
print("="*90)
print(f"Structures: {len(dump_files)}")
print(OUT)
print("TDDFT_XYZ_EXPORTER_READY")
