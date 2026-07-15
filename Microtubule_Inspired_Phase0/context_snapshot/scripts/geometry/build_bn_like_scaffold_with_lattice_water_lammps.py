#!/usr/bin/env python3
"""
Build a BN-like tubular scaffold with confined TIP3P-like water placed on a lattice.

Purpose:
- Avoid slow random sequential insertion at high water density.
- Generate larger confined-water systems efficiently.
- Write XYZ and LAMMPS data files.

This is still an initial-configuration generator. The lattice-like water placement
must later be relaxed by minimization/equilibration in LAMMPS.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np


MASS_B = 10.811
MASS_N = 14.007
MASS_O = 15.9994
MASS_H = 1.008

Q_B = +0.40
Q_N = -0.40
Q_O = -0.834
Q_H = +0.417

OH_BOND_A = 0.9572
HOH_ANGLE_DEG = 104.52


def random_unit_vector(rng: np.random.Generator) -> np.ndarray:
    v = rng.normal(size=3)
    n = np.linalg.norm(v)
    if n == 0:
        return random_unit_vector(rng)
    return v / n


def orthogonal_unit_vector(v: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    trial = random_unit_vector(rng)
    u = np.cross(v, trial)
    n = np.linalg.norm(u)
    if n < 1e-10:
        return orthogonal_unit_vector(v, rng)
    return u / n


def build_water_geometry(O: np.ndarray, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    bisector = random_unit_vector(rng)
    perp = orthogonal_unit_vector(bisector, rng)

    half_angle = np.deg2rad(HOH_ANGLE_DEG / 2.0)
    h1_dir = np.cos(half_angle) * bisector + np.sin(half_angle) * perp
    h2_dir = np.cos(half_angle) * bisector - np.sin(half_angle) * perp

    H1 = O + OH_BOND_A * h1_dir
    H2 = O + OH_BOND_A * h2_dir
    return O, H1, H2


def build_scaffold(
    outer_diameter_nm: float,
    lumen_diameter_nm: float,
    length_nm: float,
    n_theta: int,
    n_z: int,
    n_radial: int,
):
    r_outer_nm = outer_diameter_nm / 2.0
    r_inner_nm = lumen_diameter_nm / 2.0

    radii_nm = np.linspace(r_inner_nm, r_outer_nm, n_radial)
    theta_vals = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    z_vals_nm = np.linspace(-length_nm / 2.0, length_nm / 2.0, n_z)

    labels, coords, types, charges, mol_ids = [], [], [], [], []

    for iz, z_nm in enumerate(z_vals_nm):
        for ir, r_nm in enumerate(radii_nm):
            for it, theta in enumerate(theta_vals):
                x_nm = r_nm * np.cos(theta)
                y_nm = r_nm * np.sin(theta)

                if (iz + ir + it) % 2 == 0:
                    labels.append("B")
                    types.append(1)
                    charges.append(Q_B)
                else:
                    labels.append("N")
                    types.append(2)
                    charges.append(Q_N)

                coords.append([10.0 * x_nm, 10.0 * y_nm, 10.0 * z_nm])
                mol_ids.append(0)

    return labels, np.asarray(coords), np.asarray(types), np.asarray(charges), np.asarray(mol_ids)


def lattice_oxygen_positions(
    target_n_water: int,
    lumen_diameter_nm: float,
    length_nm: float,
    wall_clearance_A: float,
    spacing_A: float,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)

    lumen_radius_A = lumen_diameter_nm * 10.0 / 2.0
    usable_radius_A = lumen_radius_A - wall_clearance_A

    z_min = -length_nm * 10.0 / 2.0 + wall_clearance_A
    z_max =  length_nm * 10.0 / 2.0 - wall_clearance_A

    x_vals = np.arange(-usable_radius_A, usable_radius_A + 1e-9, spacing_A)
    y_vals = np.arange(-usable_radius_A, usable_radius_A + 1e-9, spacing_A)
    z_vals = np.arange(z_min, z_max + 1e-9, spacing_A)

    positions = []
    for z in z_vals:
        for x in x_vals:
            for y in y_vals:
                if x*x + y*y <= usable_radius_A * usable_radius_A:
                    # Small jitter to avoid a perfectly crystalline water oxygen lattice.
                    jitter = rng.uniform(-0.15, 0.15, size=3) * spacing_A
                    pos = np.array([x, y, z]) + jitter

                    # Ensure jitter did not move outside usable cylinder.
                    if pos[0]**2 + pos[1]**2 <= usable_radius_A**2 and z_min <= pos[2] <= z_max:
                        positions.append(pos)

    positions = np.asarray(positions, dtype=float)

    if len(positions) < target_n_water:
        raise RuntimeError(
            f"Only {len(positions)} lattice sites available, but {target_n_water} waters requested. "
            f"Decrease spacing_A or target_n_water."
        )

    selected = rng.choice(len(positions), size=target_n_water, replace=False)
    return positions[selected]


def build_waters(
    n_water: int,
    lumen_diameter_nm: float,
    length_nm: float,
    wall_clearance_A: float,
    spacing_A: float,
    seed: int,
):
    rng = np.random.default_rng(seed)
    O_positions = lattice_oxygen_positions(
        target_n_water=n_water,
        lumen_diameter_nm=lumen_diameter_nm,
        length_nm=length_nm,
        wall_clearance_A=wall_clearance_A,
        spacing_A=spacing_A,
        seed=seed,
    )

    labels, coords, types, charges, mol_ids = [], [], [], [], []
    bonds, angles = [], []

    atom_counter = 0
    for iw, O_pos in enumerate(O_positions, start=1):
        O, H1, H2 = build_water_geometry(O_pos, rng)

        o_id = atom_counter + 1
        h1_id = atom_counter + 2
        h2_id = atom_counter + 3

        labels.extend(["O", "H", "H"])
        coords.extend([O, H1, H2])
        types.extend([3, 4, 4])
        charges.extend([Q_O, Q_H, Q_H])
        mol_ids.extend([iw, iw, iw])

        bonds.append((1, o_id, h1_id))
        bonds.append((1, o_id, h2_id))
        angles.append((1, h1_id, o_id, h2_id))

        atom_counter += 3

    return labels, np.asarray(coords), np.asarray(types), np.asarray(charges), np.asarray(mol_ids), bonds, angles


def write_xyz(labels: list[str], coords_A: np.ndarray, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as f:
        f.write(f"{len(labels)}\n")
        f.write("BN-like scaffold with lattice-confined TIP3P-like water; coordinates in Angstrom\n")
        for label, (x, y, z) in zip(labels, coords_A):
            f.write(f"{label:2s} {x:15.6f} {y:15.6f} {z:15.6f}\n")


def write_lammps_data(
    labels,
    coords_A,
    atom_types,
    charges,
    mol_ids,
    bonds,
    angles,
    output_path: Path,
    padding_A: float = 20.0,
):
    output_path.parent.mkdir(parents=True, exist_ok=True)

    xlo = float(coords_A[:, 0].min() - padding_A)
    xhi = float(coords_A[:, 0].max() + padding_A)
    ylo = float(coords_A[:, 1].min() - padding_A)
    yhi = float(coords_A[:, 1].max() + padding_A)
    zlo = float(coords_A[:, 2].min() - padding_A)
    zhi = float(coords_A[:, 2].max() + padding_A)

    with output_path.open("w") as f:
        f.write("LAMMPS data file: BN-like scaffold with lattice-confined water\n\n")
        f.write(f"{len(labels)} atoms\n")
        f.write(f"{len(bonds)} bonds\n")
        f.write(f"{len(angles)} angles\n\n")
        f.write("4 atom types\n")
        f.write("1 bond types\n")
        f.write("1 angle types\n\n")
        f.write(f"{xlo:15.6f} {xhi:15.6f} xlo xhi\n")
        f.write(f"{ylo:15.6f} {yhi:15.6f} ylo yhi\n")
        f.write(f"{zlo:15.6f} {zhi:15.6f} zlo zhi\n\n")

        f.write("Masses\n\n")
        f.write(f"1 {MASS_B:.6f} # B-like\n")
        f.write(f"2 {MASS_N:.6f} # N-like\n")
        f.write(f"3 {MASS_O:.6f} # O water\n")
        f.write(f"4 {MASS_H:.6f} # H water\n\n")

        f.write("Atoms # full\n\n")
        for i, (mol_id, atom_type, q, xyz) in enumerate(zip(mol_ids, atom_types, charges, coords_A), start=1):
            x, y, z = xyz
            f.write(f"{i:8d} {int(mol_id):8d} {int(atom_type):3d} {q:12.6f} {x:15.6f} {y:15.6f} {z:15.6f}\n")

        f.write("\nBonds\n\n")
        for ib, (bt, a1, a2) in enumerate(bonds, start=1):
            f.write(f"{ib:8d} {bt:3d} {a1:8d} {a2:8d}\n")

        f.write("\nAngles\n\n")
        for ia, (at, a1, a2, a3) in enumerate(angles, start=1):
            f.write(f"{ia:8d} {at:3d} {a1:8d} {a2:8d} {a3:8d}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outer-diameter-nm", type=float, default=24.0)
    parser.add_argument("--lumen-diameter-nm", type=float, default=14.0)
    parser.add_argument("--length-nm", type=float, default=20.0)
    parser.add_argument("--n-theta", type=int, default=72)
    parser.add_argument("--n-z", type=int, default=60)
    parser.add_argument("--n-radial", type=int, default=4)
    parser.add_argument("--n-water", type=int, required=True)
    parser.add_argument("--wall-clearance-A", type=float, default=3.0)
    parser.add_argument("--spacing-A", type=float, default=3.1)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--xyz-output", type=Path, required=True)
    parser.add_argument("--data-output", type=Path, required=True)
    args = parser.parse_args()

    sc_labels, sc_coords, sc_types, sc_charges, sc_mol_ids = build_scaffold(
        args.outer_diameter_nm, args.lumen_diameter_nm, args.length_nm,
        args.n_theta, args.n_z, args.n_radial
    )

    water_labels, water_coords, water_types, water_charges, water_mol_ids, bonds, angles = build_waters(
        n_water=args.n_water,
        lumen_diameter_nm=args.lumen_diameter_nm,
        length_nm=args.length_nm,
        wall_clearance_A=args.wall_clearance_A,
        spacing_A=args.spacing_A,
        seed=args.seed,
    )

    water_mol_ids = water_mol_ids + 1000

    shift = len(sc_labels)
    shifted_bonds = [(bt, a1 + shift, a2 + shift) for bt, a1, a2 in bonds]
    shifted_angles = [(at, a1 + shift, a2 + shift, a3 + shift) for at, a1, a2, a3 in angles]

    labels = sc_labels + water_labels
    coords = np.vstack([sc_coords, water_coords])
    atom_types = np.concatenate([sc_types, water_types])
    charges = np.concatenate([sc_charges, water_charges])
    mol_ids = np.concatenate([sc_mol_ids, water_mol_ids])

    write_xyz(labels, coords, args.xyz_output)
    write_lammps_data(labels, coords, atom_types, charges, mol_ids, shifted_bonds, shifted_angles, args.data_output)

    print("BN-like scaffold + lattice-confined water generation complete.")
    print(f"Scaffold atoms: {len(sc_labels)}")
    print(f"Water molecules: {args.n_water}")
    print(f"Water atoms: {len(water_labels)}")
    print(f"Total atoms: {len(labels)}")
    print(f"Bonds: {len(shifted_bonds)}")
    print(f"Angles: {len(shifted_angles)}")
    print(f"Total charge: {np.sum(charges):.6f} e")
    print(f"Lattice spacing: {args.spacing_A:.3f} A")
    print(f"XYZ output: {args.xyz_output}")
    print(f"LAMMPS data output: {args.data_output}")


if __name__ == "__main__":
    main()
