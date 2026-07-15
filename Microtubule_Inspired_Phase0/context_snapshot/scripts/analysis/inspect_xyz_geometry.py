#!/usr/bin/env python3
"""
Inspect an XYZ geometry file generated for Phase 0.

Purpose:
- Read XYZ coordinates.
- Report coordinate bounds.
- Estimate radial range.
- Estimate outer diameter, lumen diameter proxy, and length.
- Save a simple radial cross-section plot.

Assumption:
- Coordinates are stored in Angstrom.
- Tube axis is z.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    with path.open("r") as f:
        lines = f.readlines()

    try:
        n_atoms = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError(f"First line of {path} is not a valid XYZ atom count.") from exc

    labels = []
    coords = []

    for line in lines[2:2 + n_atoms]:
        parts = line.split()
        if len(parts) < 4:
            continue
        labels.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    coords_arr = np.asarray(coords, dtype=float)

    if coords_arr.shape[0] != n_atoms:
        raise ValueError(
            f"XYZ atom count mismatch: header={n_atoms}, parsed={coords_arr.shape[0]}"
        )

    return labels, coords_arr


def inspect_geometry(coords_angstrom: np.ndarray) -> dict[str, float]:
    x = coords_angstrom[:, 0]
    y = coords_angstrom[:, 1]
    z = coords_angstrom[:, 2]
    r = np.sqrt(x**2 + y**2)

    summary = {
        "n_points": float(coords_angstrom.shape[0]),
        "x_min_A": float(np.min(x)),
        "x_max_A": float(np.max(x)),
        "y_min_A": float(np.min(y)),
        "y_max_A": float(np.max(y)),
        "z_min_A": float(np.min(z)),
        "z_max_A": float(np.max(z)),
        "r_min_A": float(np.min(r)),
        "r_max_A": float(np.max(r)),
        "estimated_lumen_diameter_nm": float(2.0 * np.min(r) / 10.0),
        "estimated_outer_diameter_nm": float(2.0 * np.max(r) / 10.0),
        "estimated_length_nm": float((np.max(z) - np.min(z)) / 10.0),
    }

    return summary


def save_cross_section_plot(coords_angstrom: np.ndarray, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    x_nm = coords_angstrom[:, 0] / 10.0
    y_nm = coords_angstrom[:, 1] / 10.0

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x_nm, y_nm, s=2)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_title("Phase 0 generic tubular cross-section")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xyz", type=Path, required=True)
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=Path("results/geometry/generic_tube_summary.csv"),
    )
    parser.add_argument(
        "--figure",
        type=Path,
        default=Path("figures/geometry/generic_tube_cross_section.png"),
    )
    args = parser.parse_args()

    labels, coords = read_xyz(args.xyz)
    summary = inspect_geometry(coords)

    args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([summary]).to_csv(args.summary_csv, index=False)

    save_cross_section_plot(coords, args.figure)

    print("Geometry inspection complete.")
    print(f"Input XYZ:        {args.xyz}")
    print(f"Summary CSV:      {args.summary_csv}")
    print(f"Cross-section:    {args.figure}")
    print()
    for key, value in summary.items():
        if key == "n_points":
            print(f"{key}: {int(value)}")
        else:
            print(f"{key}: {value:.6f}")


if __name__ == "__main__":
    main()
