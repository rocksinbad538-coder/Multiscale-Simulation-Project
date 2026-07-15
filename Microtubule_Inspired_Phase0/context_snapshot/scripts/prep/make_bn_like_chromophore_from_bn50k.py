#!/usr/bin/env python3

from pathlib import Path
import math
import re


INPUT = Path(
    "systems/inorganic/bn_like_scaffold_water/outputs/"
    "bn_like_scaffold_water_30000w_nvt300_hold_contained_50k.data"
)

OUTPUT = Path(
    "systems/hybrid/bn_like_chromophore_scaffold_water/outputs/"
    "bn_like_chromophore_12dipoles_30000w_initial_from_bn50k.data"
)

N_SITES = 12
R_SITE = 68.0
DIPOLE_SEPARATION_A = 2.0
Q_POS = 0.25
Q_NEG = -0.25
TYPE_POS = 5
TYPE_NEG = 6
MASS_POS = 120.0
MASS_NEG = 120.0


SECTION_NAMES = {
    "Masses",
    "Pair Coeffs",
    "Bond Coeffs",
    "Angle Coeffs",
    "Dihedral Coeffs",
    "Improper Coeffs",
    "Atoms",
    "Velocities",
    "Bonds",
    "Angles",
    "Dihedrals",
    "Impropers",
}


def is_section_header(line):
    stripped = line.strip()
    if not stripped:
        return False
    base = stripped.split("#")[0].strip()
    return base in SECTION_NAMES


def find_section(lines, name):
    start = None
    for i, line in enumerate(lines):
        if line.strip().split("#")[0].strip() == name:
            start = i
            break
    if start is None:
        raise RuntimeError(f"Could not find section: {name}")

    end = len(lines)
    for j in range(start + 1, len(lines)):
        if is_section_header(lines[j]):
            end = j
            break
    return start, end


def parse_atom_line(line):
    parts = line.split()
    if len(parts) < 7:
        return None
    try:
        atom_id = int(parts[0])
        mol_id = int(parts[1])
        atom_type = int(parts[2])
        q = float(parts[3])
        x = float(parts[4])
        y = float(parts[5])
        z = float(parts[6])
        ix = int(parts[7]) if len(parts) > 7 else 0
        iy = int(parts[8]) if len(parts) > 8 else 0
        iz = int(parts[9]) if len(parts) > 9 else 0
        return atom_id, mol_id, atom_type, q, x, y, z, ix, iy, iz
    except ValueError:
        return None


def collect_atoms(lines):
    start, end = find_section(lines, "Atoms")
    atom_records = []
    for line in lines[start + 1:end]:
        rec = parse_atom_line(line)
        if rec is not None:
            atom_records.append(rec)
    if not atom_records:
        raise RuntimeError("No atoms parsed from Atoms section.")
    return atom_records


def update_header_counts(lines, atom_count_add, new_atom_types):
    out = []
    for line in lines:
        if re.match(r"^\s*\d+\s+atoms\s*$", line):
            n = int(line.split()[0]) + atom_count_add
            out.append(f"{n} atoms\n")
        elif re.match(r"^\s*\d+\s+atom types\s*$", line):
            out.append(f"{new_atom_types} atom types\n")
        else:
            out.append(line)
    return out


def add_masses(lines):
    start, end = find_section(lines, "Masses")
    before = lines[:end]
    after = lines[end:]

    existing = set()
    for line in lines[start + 1:end]:
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            existing.add(int(parts[0]))

    additions = []
    if TYPE_POS not in existing:
        additions.append(f"{TYPE_POS} {MASS_POS:.6f} # chromophore_pos\n")
    if TYPE_NEG not in existing:
        additions.append(f"{TYPE_NEG} {MASS_NEG:.6f} # chromophore_neg\n")

    if additions:
        if before and before[-1].strip() != "":
            before.append("\n")
        before.extend(additions)
        before.append("\n")

    return before + after


def generate_chromophore_atoms(max_atom_id, max_mol_id):
    atoms = []
    velocities = []

    z_levels = [-60.0, -20.0, 20.0, 60.0]
    n_theta = 3
    site_index = 0

    for iz, zc in enumerate(z_levels):
        for it in range(n_theta):
            if site_index >= N_SITES:
                break

            theta = 2.0 * math.pi * (it / n_theta) + iz * math.pi / n_theta
            x = R_SITE * math.cos(theta)
            y = R_SITE * math.sin(theta)

            z_pos = zc - 0.5 * DIPOLE_SEPARATION_A
            z_neg = zc + 0.5 * DIPOLE_SEPARATION_A

            mol = max_mol_id + site_index + 1

            atom_id_pos = max_atom_id + 2 * site_index + 1
            atom_id_neg = max_atom_id + 2 * site_index + 2

            atoms.append(
                f"{atom_id_pos} {mol} {TYPE_POS} {Q_POS:.6f} "
                f"{x:.6f} {y:.6f} {z_pos:.6f} 0 0 0\n"
            )
            atoms.append(
                f"{atom_id_neg} {mol} {TYPE_NEG} {Q_NEG:.6f} "
                f"{x:.6f} {y:.6f} {z_neg:.6f} 0 0 0\n"
            )

            velocities.append(f"{atom_id_pos} 0.000000 0.000000 0.000000\n")
            velocities.append(f"{atom_id_neg} 0.000000 0.000000 0.000000\n")

            site_index += 1

    return atoms, velocities


def append_to_section(lines, section_name, new_lines):
    start, end = find_section(lines, section_name)
    before = lines[:end]
    after = lines[end:]

    if before and before[-1].strip() != "":
        before.append("\n")
    before.extend(new_lines)
    before.append("\n")

    return before + after


def main():
    if not INPUT.exists():
        raise FileNotFoundError(INPUT)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    lines = INPUT.read_text().splitlines(keepends=True)

    atom_records = collect_atoms(lines)
    max_atom_id = max(a[0] for a in atom_records)
    max_mol_id = max(a[1] for a in atom_records)

    new_atoms, new_velocities = generate_chromophore_atoms(max_atom_id, max_mol_id)

    lines = update_header_counts(lines, atom_count_add=len(new_atoms), new_atom_types=6)
    lines = add_masses(lines)
    lines = append_to_section(lines, "Atoms", new_atoms)
    lines = append_to_section(lines, "Velocities", new_velocities)

    OUTPUT.write_text("".join(lines))

    print("Created hybrid/chromophore-bearing model")
    print(f"Input:  {INPUT}")
    print(f"Output: {OUTPUT}")
    print(f"Added pseudo-atoms: {len(new_atoms)}")
    print(f"Added chromophore dipoles: {len(new_atoms)//2}")
    print(f"Atom types added: {TYPE_POS}, {TYPE_NEG}")
    print(f"Net added charge: {(len(new_atoms)//2)*(Q_POS+Q_NEG):.6f} e")


if __name__ == "__main__":
    main()
