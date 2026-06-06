#!/usr/bin/env python3
"""
Build a simple cylindrical shell point-cloud geometry for Phase 0.

This is not yet a force-field-ready molecular model. It is a geometric scaffold
used to verify dimensions, coordinate conventions, visualization, and downstream
candidate construction.

Units:
    nm for input dimensions.
    Angstrom for output XYZ coordinates, because many visualization tools expect Å.

Output:
    systems/<candidate>/geometry/<name>.xyz
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np


def build_cylindrical_shell(
    outer_diameter_nm: float,
    lumen_diameter_nm: float,
    length_nm: float,
    n_theta: int,
    n_z: int,
    n_radial: int,
) -> np.ndarray:
    """Generate points in a cylindrical shell."""
    r_outer = outer_diameter_nm / 2.0
    r_inner = lumen_diameter_nm / 2.0

    if r_inner <= 0:
        raise ValueError("lumen_diameter_nm must be positive.")
    if r_outer <= r_inner:
        raise ValueError("outer_diameter_nm must be larger than lumen_diameter_nm.")
    if length_nm <= 0:
        raise ValueError("length_nm must be positive.")

    radii = np.linspace(r_inner, r_outer, n_radial)
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    z = np.linspace(-length_nm / 2.0, length_nm / 2.0, n_z)

    coords = []
    for zz in z:
        for rr in radii:
            for tt in theta:
                x = rr * np.cos(tt)
                y = rr * np.sin(tt)
                coords.append([x, y, zz])

    return np.asarray(coords, dtype=float)


def write_xyz(coords_nm: np.ndarray, output_path: Path, atom_label: str = "X") -> None:
    """Write coordinates to XYZ in Angstrom."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    coords_angstrom = coords_nm * 10.0

    with output_path.open("w") as f:
        f.write(f"{len(coords_angstrom)}\n")
        f.write("Generic cylindrical shell geometry for Phase 0; coordinates in Angstrom\n")
        for x, y, z in coords_angstrom:
            f.write(f"{atom_label:2s} {x:15.6f} {y:15.6f} {z:15.6f}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outer-diameter-nm", type=float, default=24.0)
    parser.add_argument("--lumen-diameter-nm", type=float, default=14.0)
    parser.add_argument("--length-nm", type=float, default=20.0)
    parser.add_argument("--n-theta", type=int, default=72)
    parser.add_argument("--n-z", type=int, default=60)
    parser.add_argument("--n-radial", type=int, default=4)
    parser.add_argument("--atom-label", type=str, default="X")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("systems/inorganic/geometry/generic_tube.xyz"),
    )
    args = parser.parse_args()

    coords = build_cylindrical_shell(
        outer_diameter_nm=args.outer_diameter_nm,
        lumen_diameter_nm=args.lumen_diameter_nm,
        length_nm=args.length_nm,
        n_theta=args.n_theta,
        n_z=args.n_z,
        n_radial=args.n_radial,
    )
    write_xyz(coords, args.output, args.atom_label)

    print(f"Wrote {len(coords)} shell points to: {args.output}")
    print("Geometry summary:")
    print(f"  outer diameter: {args.outer_diameter_nm:.3f} nm")
    print(f"  lumen diameter: {args.lumen_diameter_nm:.3f} nm")
    print(f"  length:         {args.length_nm:.3f} nm")


if __name__ == "__main__":
    main()
