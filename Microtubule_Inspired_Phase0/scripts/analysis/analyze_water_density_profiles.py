#!/usr/bin/env python3
"""
Analyze radial and axial water oxygen distributions from LAMMPS dump trajectories.

Phase 0 utility:
- Reads selected LAMMPS custom dump frames.
- Uses water oxygen atoms, assumed type 3.
- Computes radial histogram r = sqrt(x^2 + y^2).
- Computes axial histogram z.
- Writes CSV summaries and figures.

This is intended for workflow/stability diagnostics, not final production physics.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_lammpstrj(path: Path, oxygen_type: int = 3):
    with path.open("r") as f:
        while True:
            line = f.readline()
            if not line:
                break

            if not line.startswith("ITEM: TIMESTEP"):
                continue

            timestep = int(f.readline().strip())

            line = f.readline()
            if not line.startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError(f"Unexpected dump format near timestep {timestep}")

            n_atoms = int(f.readline().strip())

            line = f.readline()
            if not line.startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"Unexpected dump format near timestep {timestep}")

            # skip box bounds
            for _ in range(3):
                f.readline()

            header = f.readline().strip()
            if not header.startswith("ITEM: ATOMS"):
                raise ValueError(f"Unexpected atom header near timestep {timestep}: {header}")

            cols = header.split()[2:]
            col_index = {name: i for i, name in enumerate(cols)}

            required = ["type", "x", "y", "z"]
            for r in required:
                if r not in col_index:
                    raise ValueError(f"Column {r} not found in dump atom columns: {cols}")

            xyz = []
            for _ in range(n_atoms):
                parts = f.readline().split()
                atype = int(parts[col_index["type"]])
                if atype == oxygen_type:
                    x = float(parts[col_index["x"]])
                    y = float(parts[col_index["y"]])
                    z = float(parts[col_index["z"]])
                    xyz.append((x, y, z))

            yield timestep, np.asarray(xyz, dtype=float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trajectory", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--oxygen-type", type=int, default=3)
    ap.add_argument("--r-max-A", type=float, default=80.0)
    ap.add_argument("--z-min-A", type=float, default=-120.0)
    ap.add_argument("--z-max-A", type=float, default=120.0)
    ap.add_argument("--n-r-bins", type=int, default=80)
    ap.add_argument("--n-z-bins", type=int, default=120)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    traj = Path(args.trajectory)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    r_bins = np.linspace(0.0, args.r_max_A, args.n_r_bins + 1)
    z_bins = np.linspace(args.z_min_A, args.z_max_A, args.n_z_bins + 1)

    radial_records = []
    axial_records = []

    frame_count = 0
    for timestep, oxy in parse_lammpstrj(traj, oxygen_type=args.oxygen_type):
        if oxy.size == 0:
            continue

        r = np.sqrt(oxy[:, 0] ** 2 + oxy[:, 1] ** 2)
        z = oxy[:, 2]

        r_hist, _ = np.histogram(r, bins=r_bins)
        z_hist, _ = np.histogram(z, bins=z_bins)

        r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
        z_centers = 0.5 * (z_bins[:-1] + z_bins[1:])

        # Normalize to fractions per bin. This is not volumetric density yet.
        r_frac = r_hist / r_hist.sum()
        z_frac = z_hist / z_hist.sum()

        for c, val in zip(r_centers, r_frac):
            radial_records.append(
                {"run": args.label, "timestep": timestep, "r_A": c, "fraction": val}
            )

        for c, val in zip(z_centers, z_frac):
            axial_records.append(
                {"run": args.label, "timestep": timestep, "z_A": c, "fraction": val}
            )

        frame_count += 1

    radial = pd.DataFrame(radial_records)
    axial = pd.DataFrame(axial_records)

    radial_csv = outdir / f"{args.label}_radial_profile.csv"
    axial_csv = outdir / f"{args.label}_axial_profile.csv"

    radial.to_csv(radial_csv, index=False)
    axial.to_csv(axial_csv, index=False)

    # Plot initial, middle, final frames
    def selected_timesteps(df):
        ts = sorted(df["timestep"].unique())
        if len(ts) <= 3:
            return ts
        return [ts[0], ts[len(ts)//2], ts[-1]]

    plt.figure(figsize=(8, 5))
    for ts in selected_timesteps(radial):
        sub = radial[radial["timestep"] == ts]
        plt.plot(sub["r_A"], sub["fraction"], label=f"step {ts}")
    plt.xlabel("Water oxygen radial position r, Å")
    plt.ylabel("Fraction per radial bin")
    plt.title(f"Radial water oxygen distribution: {args.label}")
    plt.legend()
    plt.tight_layout()
    radial_png = outdir / f"{args.label}_radial_profile.png"
    plt.savefig(radial_png, dpi=200)
    plt.close()

    plt.figure(figsize=(8, 5))
    for ts in selected_timesteps(axial):
        sub = axial[axial["timestep"] == ts]
        plt.plot(sub["z_A"], sub["fraction"], label=f"step {ts}")
    plt.xlabel("Water oxygen axial position z, Å")
    plt.ylabel("Fraction per axial bin")
    plt.title(f"Axial water oxygen distribution: {args.label}")
    plt.legend()
    plt.tight_layout()
    axial_png = outdir / f"{args.label}_axial_profile.png"
    plt.savefig(axial_png, dpi=200)
    plt.close()

    print("Water density/profile analysis complete.")
    print(f"Frames processed: {frame_count}")
    print(f"Radial CSV: {radial_csv}")
    print(f"Axial CSV:   {axial_csv}")
    print(f"Radial PNG: {radial_png}")
    print(f"Axial PNG:  {axial_png}")


if __name__ == "__main__":
    main()
