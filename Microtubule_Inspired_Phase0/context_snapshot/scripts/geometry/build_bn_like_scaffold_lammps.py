#!/usr/bin/env python3
"""
Build a simplified BN-like tubular scaffold for Phase 0.

Purpose:
- Generate a geometry-equivalent inorganic tubular scaffold.
- Assign alternating B/N-like atom types.
- Write both XYZ and LAMMPS data formats.
- This is a Phase 0 screening scaffold, not yet a chemically exact BNNT.

Units:
- Input dimensions: nm
- Output coordinates: Angstrom
- LAMMPS units expected later: real
- LAMMPS atom_style: charge
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np


MASS_B = 10.811
MASS_N = 14.007


def build_shell_points(
    outer_diameter_nm: float,
    lumen_diameter_nm: float,
    length_nm: float,
    n_theta: int,
    n_z: int,
    n_radial: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r_outer = outer_diameter_nm / 2.0
    r_inner = lumen_diameter_nm / 2.0

    if r_outer <= r_inner:
        raise ValueError("Outer diameter must be larger than lumen diameter.")

    radii = np.linspace(r_inner, r_outer, n_radial)
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    zvals = np.linspace(-length_nm / 2.0, length_nm / 2.0, n_z)

    coords = []
    type_ids = []
    charges = []
    labels = []

    atom_id = 0
    for iz, z in enumerate(zvals):
        for ir, r in enumerate(radii):
            for it, t in enumerate(theta):
                x = r * np.cos(t)
                y = r * np.sin(t)

                # Alternating B/N-like assignment.
                # This is only a screening scaffold pattern.
                parity = (iz + ir + it) % 2

                if parity == 0:
                    atom_type = 1
                    charge = +0.40
                    label = "B"
                else:
                    atom_type = 2
                    charge = -0.40
                    label = "N"

                coords.append([x, y, z])
                type_ids.append(atom_type)
                charges.append(charge)
                labels.append(label)
                atom_id += 1

    coords_nm = np.asarray(coords, dtype=float)
    type_ids_arr = np.asarray(type_ids, dtype=int)
    charges_arr = np.asarray(charges, dtype=float)
    labels_arr = np.asarray(labels, dtype=str)

    return coords_nm, type_ids_arr, charges_arr, labels_arr


def write_xyz(coords_nm: np.ndarray, labels: np.ndarray, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    coords_A = coords_nm * 10.0

    with output_path.open("w") as f:
        f.write(f"{len(coords_A)}\n")
        f.write("BN-like Phase 0 scaffold; coordinates in Angstrom\n")
        for label, (x, y, z) in zip(labels, coords_A):
            f.write(f"{label:2s} {x:15.6f} {y:15.6f} {z:15.6f}\n")


def write_lammps_data(
    coords_nm: np.ndarray,
    type_ids: np.ndarray,
    charges: np.ndarray,
    output_path: Path,
    padding_A: float = 20.0,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    coords_A = coords_nm * 10.0

    xlo = float(coords_A[:, 0].min() - padding_A)
    xhi = float(coords_A[:, 0].max() + padding_A)
    ylo = float(coords_A[:, 1].min() - padding_A)
    yhi = float(coords_A[:, 1].max() + padding_A)
    zlo = float(coords_A[:, 2].min() - padding_A)
    zhi = float(coords_A[:, 2].max() + padding_A)

    with output_path.open("w") as f:
        f.write("LAMMPS data file: BN-like Phase 0 tubular scaffold\n\n")
        f.write(f"{len(coords_A)} atoms\n")
        f.write("2 atom types\n\n")

        f.write(f"{xlo:15.6f} {xhi:15.6f} xlo xhi\n")
        f.write(f"{ylo:15.6f} {yhi:15.6f} ylo yhi\n")
        f.write(f"{zlo:15.6f} {zhi:15.6f} zlo zhi\n\n")

        f.write("Masses\n\n")
        f.write(f"1 {MASS_B:.6f} # B-like\n")
        f.write(f"2 {MASS_N:.6f} # N-like\n\n")

        f.write("Atoms # charge\n\n")
        for i, ((x, y, z), atom_type, q) in enumerate(zip(coords_A, type_ids, charges), start=1):
            f.write(f"{i:8d} {atom_type:3d} {q:12.6f} {x:15.6f} {y:15.6f} {z:15.6f}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outer-diameter-nm", type=float, default=24.0)
    parser.add_argument("--lumen-diameter-nm", type=float, default=14.0)
    parser.add_argument("--length-nm", type=float, default=20.0)
    parser.add_argument("--n-theta", type=int, default=72)
    parser.add_argument("--n-z", type=int, default=60)
    parser.add_argument("--n-radial", type=int, default=4)
    parser.add_argument("--xyz-output", type=Path, required=True)
    parser.add_argument("--data-output", type=Path, required=True)
    args = parser.parse_args()

    coords, type_ids, charges, labels = build_shell_points(
        outer_diameter_nm=args.outer_diameter_nm,
        lumen_diameter_nm=args.lumen_diameter_nm,
        length_nm=args.length_nm,
        n_theta=args.n_theta,
        n_z=args.n_z,
        n_radial=args.n_radial,
    )

    write_xyz(coords, labels, args.xyz_output)
    write_lammps_data(coords, type_ids, charges, args.data_output)

    n_b = int(np.sum(type_ids == 1))
    n_n = int(np.sum(type_ids == 2))
    q_total = float(np.sum(charges))

    print("BN-like scaffold generation complete.")
    print(f"Atoms total: {len(coords)}")
    print(f"B-like atoms: {n_b}")
    print(f"N-like atoms: {n_n}")
    print(f"Total charge: {q_total:.6f} e")
    print(f"XYZ output: {args.xyz_output}")
    print(f"LAMMPS data output: {args.data_output}")


if __name__ == "__main__":
    main()
