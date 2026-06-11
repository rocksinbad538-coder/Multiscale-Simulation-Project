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
    "bn_like_chromophore_12dipoles_carved_initial_from_bn50k.data"
)

N_SITES = 12

# Keep sites slightly inside the nominal inner wall, but not too deep into dense water.
R_SITE = 66.0

DIPOLE_SEPARATION_A = 2.0

# Start with smaller charges for initial stability screening.
# We can ramp later after geometry is stable.
Q_POS = 0.10
Q_NEG = -0.10

TYPE_POS = 5
TYPE_NEG = 6
MASS_POS = 120.0
MASS_NEG = 120.0

WATER_O_TYPE = 3
WATER_H_TYPE = 4

# Remove any entire water molecule whose oxygen is closer than this
# to any chromophore pseudo-site.
EXCLUSION_RADIUS_A = 4.5

SECTION_HEADERS = {
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


def clean_header_name(line):
    return line.strip().split("#")[0].strip()


def find_section(lines, name):
    start = None
    for i, line in enumerate(lines):
        if clean_header_name(line) == name:
            start = i
            break
    if start is None:
        raise RuntimeError(f"Section not found: {name}")

    end = len(lines)
    for j in range(start + 1, len(lines)):
        if clean_header_name(lines[j]) in SECTION_HEADERS:
            end = j
            break

    return start, end


def extract_records(lines, section_name, min_fields):
    start, end = find_section(lines, section_name)
    records = []
    for line in lines[start + 1:end]:
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) >= min_fields:
            try:
                int(parts[0])
                records.append(s)
            except ValueError:
                pass
    return records


def parse_atom(rec):
    p = rec.split()
    return {
        "id": int(p[0]),
        "mol": int(p[1]),
        "type": int(p[2]),
        "q": float(p[3]),
        "x": float(p[4]),
        "y": float(p[5]),
        "z": float(p[6]),
        "ix": int(p[7]) if len(p) > 7 else 0,
        "iy": int(p[8]) if len(p) > 8 else 0,
        "iz": int(p[9]) if len(p) > 9 else 0,
        "raw": rec,
    }


def parse_bond(rec):
    p = rec.split()
    return int(p[0]), int(p[1]), int(p[2]), int(p[3]), rec


def parse_angle(rec):
    p = rec.split()
    return int(p[0]), int(p[1]), int(p[2]), int(p[3]), int(p[4]), rec


def chromophore_sites(max_atom_id, max_mol_id):
    atoms = []
    velocities = []
    coords = []

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
                f"{x:.6f} {y:.6f} {z_pos:.6f} 0 0 0"
            )
            atoms.append(
                f"{atom_id_neg} {mol} {TYPE_NEG} {Q_NEG:.6f} "
                f"{x:.6f} {y:.6f} {z_neg:.6f} 0 0 0"
            )

            velocities.append(f"{atom_id_pos} 0.000000 0.000000 0.000000")
            velocities.append(f"{atom_id_neg} 0.000000 0.000000 0.000000")

            coords.append((x, y, z_pos))
            coords.append((x, y, z_neg))

            site_index += 1

    return atoms, velocities, coords


def dist2(a, b):
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return dx * dx + dy * dy + dz * dz


def main():
    if not INPUT.exists():
        raise FileNotFoundError(INPUT)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    lines = INPUT.read_text().splitlines()

    atom_records_raw = extract_records(lines, "Atoms", 7)
    vel_records_raw = extract_records(lines, "Velocities", 4)
    bond_records_raw = extract_records(lines, "Bonds", 4)
    angle_records_raw = extract_records(lines, "Angles", 5)

    atoms = [parse_atom(r) for r in atom_records_raw]
    max_atom_id = max(a["id"] for a in atoms)
    max_mol_id = max(a["mol"] for a in atoms)

    new_atoms, new_vels, chrom_coords = chromophore_sites(max_atom_id, max_mol_id)

    # Identify water molecules to remove based on oxygen distance to any chromophore pseudo-site.
    cutoff2 = EXCLUSION_RADIUS_A ** 2
    mols_to_remove = set()

    for a in atoms:
        if a["type"] != WATER_O_TYPE:
            continue
        pos = (a["x"], a["y"], a["z"])
        if any(dist2(pos, c) < cutoff2 for c in chrom_coords):
            mols_to_remove.add(a["mol"])

    keep_atom_ids = set()
    kept_atoms = []
    removed_atoms = []

    for a in atoms:
        if a["mol"] in mols_to_remove and a["type"] in (WATER_O_TYPE, WATER_H_TYPE):
            removed_atoms.append(a)
        else:
            kept_atoms.append(a)
            keep_atom_ids.add(a["id"])

    # Add chromophore pseudo-atoms.
    for rec in new_atoms:
        p = rec.split()
        keep_atom_ids.add(int(p[0]))

    # Filter velocities.
    vel_by_id = {}
    for rec in vel_records_raw:
        p = rec.split()
        vel_by_id[int(p[0])] = rec

    kept_vels = []
    for a in kept_atoms:
        if a["id"] in vel_by_id:
            kept_vels.append(vel_by_id[a["id"]])
        else:
            kept_vels.append(f'{a["id"]} 0.000000 0.000000 0.000000')

    kept_vels.extend(new_vels)

    # Filter bonds and angles that reference removed atoms.
    kept_bonds = []
    for rec in bond_records_raw:
        bond_id, bond_type, a1, a2, raw = parse_bond(rec)
        if a1 in keep_atom_ids and a2 in keep_atom_ids:
            kept_bonds.append(raw)

    kept_angles = []
    for rec in angle_records_raw:
        angle_id, angle_type, a1, a2, a3, raw = parse_angle(rec)
        if a1 in keep_atom_ids and a2 in keep_atom_ids and a3 in keep_atom_ids:
            kept_angles.append(raw)

    kept_atom_lines = [a["raw"] for a in kept_atoms] + new_atoms

    n_atoms = len(kept_atom_lines)
    n_bonds = len(kept_bonds)
    n_angles = len(kept_angles)

    # Preserve box and header lines up to Masses, but rewrite counts.
    header = []
    for line in lines:
        if clean_header_name(line) == "Masses":
            break

        if re.match(r"^\s*\d+\s+atoms\s*$", line):
            header.append(f"{n_atoms} atoms")
        elif re.match(r"^\s*\d+\s+bonds\s*$", line):
            header.append(f"{n_bonds} bonds")
        elif re.match(r"^\s*\d+\s+angles\s*$", line):
            header.append(f"{n_angles} angles")
        elif re.match(r"^\s*\d+\s+atom types\s*$", line):
            header.append("6 atom types")
        else:
            header.append(line)

    # Ensure required header counts are correct even if not captured.
    text = "\n".join(header)

    sections = []

    sections += [
        "Masses",
        "",
        "1 10.811000 # B_like_scaffold",
        "2 14.007000 # N_like_scaffold",
        "3 15.999400 # water_O",
        "4 1.008000 # water_H",
        "5 120.000000 # chromophore_pos",
        "6 120.000000 # chromophore_neg",
        "",
        "Pair Coeffs",
        "",
        "1 0.050000 3.400000 # B_like_scaffold",
        "2 0.050000 3.400000 # N_like_scaffold",
        "3 0.152100 3.150700 # water_O",
        "4 0.000000 0.000000 # water_H",
        "5 0.120000 3.500000 # chromophore_pos",
        "6 0.120000 3.500000 # chromophore_neg",
        "",
        "Bond Coeffs # harmonic",
        "",
        "1 450.0 0.9572",
        "",
        "Angle Coeffs # harmonic",
        "",
        "1 55.0 104.52",
        "",
        "Atoms # full",
        "",
    ]

    sections += kept_atom_lines
    sections += ["", "Velocities", ""]
    sections += kept_vels
    sections += ["", "Bonds", ""]
    sections += kept_bonds
    sections += ["", "Angles", ""]
    sections += kept_angles
    sections += [""]

    OUTPUT.write_text(text + "\n" + "\n".join(sections))

    print("Created carved hybrid/chromophore-bearing model")
    print(f"Input:  {INPUT}")
    print(f"Output: {OUTPUT}")
    print(f"Chromophore dipoles: {N_SITES}")
    print(f"Chromophore pseudo-atoms added: {len(new_atoms)}")
    print(f"Chromophore charge magnitude: +/- {Q_POS:.3f} e")
    print(f"Exclusion radius around chromophore pseudo-sites: {EXCLUSION_RADIUS_A:.2f} A")
    print(f"Water molecules removed: {len(mols_to_remove)}")
    print(f"Final atoms: {n_atoms}")
    print(f"Final bonds: {n_bonds}")
    print(f"Final angles: {n_angles}")
    print(f"Net added chromophore charge: {(N_SITES * (Q_POS + Q_NEG)):.6f} e")


if __name__ == "__main__":
    main()
