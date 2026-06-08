#!/usr/bin/env python3

from pathlib import Path
import argparse


SECTION_NAMES = {
    "Masses",
    "Pair Coeffs",
    "Bond Coeffs",
    "Angle Coeffs",
    "Atoms",
    "Velocities",
    "Bonds",
    "Angles",
}


def read_sections(path: Path):
    lines = path.read_text(errors="ignore").splitlines()

    header = []
    sections = {}
    current = None

    for line in lines:
        stripped = line.strip()

        if stripped in SECTION_NAMES or stripped.startswith("Atoms #") or stripped.startswith("Pair Coeffs #"):
            name = stripped.split("#")[0].strip()
            current = name
            sections[current] = []
            continue

        if current is None:
            header.append(line)
        else:
            sections[current].append(line)

    return header, sections


def clean_section_lines(lines):
    clean = []
    for ln in lines:
        s = ln.strip()
        if not s:
            continue
        if s.startswith("#"):
            continue
        clean.append(ln)
    return clean


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-data", required=True)
    ap.add_argument("--output-data", required=True)
    ap.add_argument("--oxygen-type", type=int, default=3)
    ap.add_argument("--hydrogen-type", type=int, default=4)
    args = ap.parse_args()

    inp = Path(args.input_data)
    out = Path(args.output_data)
    out.parent.mkdir(parents=True, exist_ok=True)

    if not inp.exists():
        raise FileNotFoundError(f"Input data file not found: {inp}")

    _header, sections = read_sections(inp)

    required = ["Atoms", "Bonds", "Angles"]
    missing = [s for s in required if s not in sections]
    if missing:
        raise RuntimeError(f"Missing required LAMMPS data sections: {missing}")

    atom_lines = clean_section_lines(sections["Atoms"])
    vel_lines = clean_section_lines(sections.get("Velocities", []))
    bond_lines = clean_section_lines(sections["Bonds"])
    angle_lines = clean_section_lines(sections["Angles"])

    keep_types = {args.oxygen_type, args.hydrogen_type}

    old_to_new = {}
    atom_records = []

    for ln in atom_lines:
        p = ln.split()
        if len(p) < 7:
            continue

        old_id = int(p[0])
        atom_type = int(p[2])

        if atom_type not in keep_types:
            continue

        old_to_new[old_id] = len(old_to_new) + 1
        atom_records.append((old_id, p))

    new_atom_lines = []
    for old_id, p in atom_records:
        p = p.copy()
        p[0] = str(old_to_new[old_id])
        new_atom_lines.append(" ".join(p))

    new_vel_lines = []
    for ln in vel_lines:
        p = ln.split()
        if len(p) < 4:
            continue

        old_id = int(p[0])
        if old_id in old_to_new:
            p[0] = str(old_to_new[old_id])
            new_vel_lines.append(" ".join(p))

    new_bond_lines = []
    bond_id = 1
    for ln in bond_lines:
        p = ln.split()
        if len(p) < 4:
            continue

        a = int(p[2])
        b = int(p[3])

        if a in old_to_new and b in old_to_new:
            new_bond_lines.append(f"{bond_id} {p[1]} {old_to_new[a]} {old_to_new[b]}")
            bond_id += 1

    new_angle_lines = []
    angle_id = 1
    for ln in angle_lines:
        p = ln.split()
        if len(p) < 5:
            continue

        a = int(p[2])
        b = int(p[3])
        c = int(p[4])

        if a in old_to_new and b in old_to_new and c in old_to_new:
            new_angle_lines.append(
                f"{angle_id} {p[1]} {old_to_new[a]} {old_to_new[b]} {old_to_new[c]}"
            )
            angle_id += 1

    n_atoms = len(new_atom_lines)
    n_bonds = len(new_bond_lines)
    n_angles = len(new_angle_lines)

    if n_atoms == 0:
        raise RuntimeError("No water atoms were selected. Check oxygen/hydrogen type IDs.")

    if n_bonds == 0 or n_angles == 0:
        raise RuntimeError("No water bonds/angles were selected. Check topology and atom type IDs.")

    text = []
    text.append("LAMMPS water-only control data generated from hydrated BN-like scaffold-water state")
    text.append("")
    text.append(f"{n_atoms} atoms")
    text.append(f"{n_bonds} bonds")
    text.append(f"{n_angles} angles")
    text.append("")
    text.append("4 atom types")
    text.append("1 bond types")
    text.append("1 angle types")
    text.append("")
    text.append("-140 140 xlo xhi")
    text.append("-140 140 ylo yhi")
    text.append("-120 120 zlo zhi")
    text.append("")
    text.append("Masses")
    text.append("")
    text.append("1 10.811")
    text.append("2 14.007")
    text.append("3 15.9994")
    text.append("4 1.008")
    text.append("")
    text.append("Atoms # full")
    text.append("")
    text.extend(new_atom_lines)
    text.append("")
    text.append("Velocities")
    text.append("")
    text.extend(new_vel_lines)
    text.append("")
    text.append("Bonds")
    text.append("")
    text.extend(new_bond_lines)
    text.append("")
    text.append("Angles")
    text.append("")
    text.extend(new_angle_lines)
    text.append("")

    out.write_text("\n".join(text))

    print(f"Wrote {out}")
    print(f"Water atoms: {n_atoms}")
    print(f"Water molecules inferred from angles: {n_angles}")
    print(f"Bonds: {n_bonds}")
    print(f"Angles: {n_angles}")


if __name__ == "__main__":
    main()
