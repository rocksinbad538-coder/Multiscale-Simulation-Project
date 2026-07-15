#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_frames(path: Path, oxygen_type: int = 3):
    with path.open("r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue

            timestep = int(f.readline().strip())

            if not f.readline().startswith("ITEM: NUMBER OF ATOMS"):
                raise RuntimeError("Unexpected dump format")
            n_atoms = int(f.readline().strip())

            if not f.readline().startswith("ITEM: BOX BOUNDS"):
                raise RuntimeError("Unexpected dump format")

            bounds = []
            for _ in range(3):
                lo, hi = map(float, f.readline().split()[:2])
                bounds.append((lo, hi))

            header = f.readline().strip()
            cols = header.split()[2:]
            idx = {c: i for i, c in enumerate(cols)}

            needed = ["type", "x", "y", "z"]
            for c in needed:
                if c not in idx:
                    raise RuntimeError(f"Missing column {c}; columns are {cols}")

            xyz = []
            for _ in range(n_atoms):
                parts = f.readline().split()
                if int(parts[idx["type"]]) == oxygen_type:
                    xyz.append([
                        float(parts[idx["x"]]),
                        float(parts[idx["y"]]),
                        float(parts[idx["z"]]),
                    ])

            yield timestep, np.asarray(xyz), bounds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trajectory", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--oxygen-type", type=int, default=3)
    ap.add_argument("--r-max-A", type=float, default=80.0)
    ap.add_argument("--z-half-length-A", type=float, default=100.0)
    ap.add_argument("--n-r-bins", type=int, default=80)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    r_edges = np.linspace(0.0, args.r_max_A, args.n_r_bins + 1)
    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
    z_length = 2.0 * args.z_half_length_A

    shell_volumes = np.pi * (r_edges[1:]**2 - r_edges[:-1]**2) * z_length

    records = []

    for timestep, oxy, _bounds in read_frames(Path(args.trajectory), args.oxygen_type):
        r = np.sqrt(oxy[:, 0]**2 + oxy[:, 1]**2)
        z = oxy[:, 2]

        inside_z = np.abs(z) <= args.z_half_length_A
        r_inside = r[inside_z]

        counts, _ = np.histogram(r_inside, bins=r_edges)
        number_density = counts / shell_volumes

        for rc, count, rho in zip(r_centers, counts, number_density):
            records.append({
                "run": args.label,
                "timestep": timestep,
                "r_A": rc,
                "count": int(count),
                "shell_volume_A3": float(shell_volumes[np.where(r_centers == rc)][0]),
                "number_density_O_per_A3": float(rho),
            })

    df = pd.DataFrame(records)
    csv = out / f"{args.label}_radial_number_density.csv"
    df.to_csv(csv, index=False)

    ts = sorted(df["timestep"].unique())
    selected = [ts[0], ts[len(ts)//2], ts[-1]] if len(ts) > 2 else ts

    plt.figure(figsize=(8, 5))
    for t in selected:
        sub = df[df["timestep"] == t]
        plt.plot(sub["r_A"], sub["number_density_O_per_A3"], label=f"step {t}")
    plt.xlabel("Radial position r, Å")
    plt.ylabel("Water oxygen number density, O/Å³")
    plt.title(f"Volume-corrected radial water oxygen density: {args.label}")
    plt.legend()
    plt.tight_layout()

    png = out / f"{args.label}_radial_number_density.png"
    plt.savefig(png, dpi=200)
    plt.close()

    print("Radial number-density analysis complete.")
    print(f"CSV: {csv}")
    print(f"PNG: {png}")


if __name__ == "__main__":
    main()
