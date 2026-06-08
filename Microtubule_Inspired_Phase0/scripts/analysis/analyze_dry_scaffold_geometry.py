#!/usr/bin/env python3

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_lammps_dump(path):
    with Path(path).open() as f:
        while True:
            line = f.readline()
            if not line:
                break

            if not line.startswith("ITEM: TIMESTEP"):
                continue

            timestep = int(f.readline().strip())

            assert f.readline().startswith("ITEM: NUMBER OF ATOMS")
            n_atoms = int(f.readline().strip())

            assert f.readline().startswith("ITEM: BOX BOUNDS")
            bounds = [list(map(float, f.readline().split()[:2])) for _ in range(3)]

            header = f.readline().split()
            assert header[0] == "ITEM:" and header[1] == "ATOMS"
            fields = header[2:]
            idx = {name: i for i, name in enumerate(fields)}

            xyz = np.zeros((n_atoms, 3), dtype=float)
            types = np.zeros(n_atoms, dtype=int)

            for i in range(n_atoms):
                parts = f.readline().split()
                types[i] = int(parts[idx["type"]])
                xyz[i, 0] = float(parts[idx["x"]])
                xyz[i, 1] = float(parts[idx["y"]])
                xyz[i, 2] = float(parts[idx["z"]])

            yield timestep, types, xyz, bounds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trajectory", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    records = []

    for timestep, types, xyz, bounds in read_lammps_dump(args.trajectory):
        r = np.sqrt(xyz[:, 0]**2 + xyz[:, 1]**2)
        z = xyz[:, 2]

        records.append({
            "run": args.label,
            "timestep": timestep,
            "n_atoms": len(xyz),
            "r_mean_A": float(r.mean()),
            "r_std_A": float(r.std()),
            "r_min_A": float(r.min()),
            "r_max_A": float(r.max()),
            "z_mean_A": float(z.mean()),
            "z_min_A": float(z.min()),
            "z_max_A": float(z.max()),
            "z_span_A": float(z.max() - z.min()),
            "fraction_outside_outer_radius_120A": float(np.mean(r > 120.0)),
            "fraction_inside_lumen_radius_70A": float(np.mean(r < 70.0)),
            "fraction_outside_z_box_120A": float(np.mean(np.abs(z) > 120.0)),
        })

    df = pd.DataFrame(records)

    csv_path = outdir / f"{args.label}_geometry_summary.csv"
    df.to_csv(csv_path, index=False)

    # Post-transient summary using step >= 3500 for consistency with the contained water analysis.
    stable = df[df["timestep"] >= 3500].copy()
    post = {
        "run": args.label,
        "n_records_total": len(df),
        "n_records_post_transient": len(stable),
        "post_transient_step_min": stable["timestep"].min(),
        "post_transient_step_max": stable["timestep"].max(),
    }

    for c in [
        "r_mean_A",
        "r_std_A",
        "r_min_A",
        "r_max_A",
        "z_min_A",
        "z_max_A",
        "z_span_A",
        "fraction_outside_outer_radius_120A",
        "fraction_inside_lumen_radius_70A",
        "fraction_outside_z_box_120A",
    ]:
        post[f"{c}_mean"] = stable[c].mean()
        post[f"{c}_std"] = stable[c].std()
        post[f"{c}_min"] = stable[c].min()
        post[f"{c}_max"] = stable[c].max()
        post[f"{c}_final"] = stable[c].iloc[-1]

    post_path = outdir / f"{args.label}_post_transient_geometry_summary.csv"
    pd.DataFrame([post]).to_csv(post_path, index=False)

    plt.figure(figsize=(8, 5))
    plt.plot(df["timestep"], df["r_mean_A"], label="mean radius")
    plt.plot(df["timestep"], df["r_min_A"], label="min radius")
    plt.plot(df["timestep"], df["r_max_A"], label="max radius")
    plt.xlabel("Timestep")
    plt.ylabel("Radius, Å")
    plt.title(f"Dry scaffold radial geometry: {args.label}")
    plt.legend()
    plt.tight_layout()
    radial_png = outdir / f"{args.label}_radial_geometry.png"
    plt.savefig(radial_png, dpi=200)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(df["timestep"], df["z_span_A"], label="z span")
    plt.plot(df["timestep"], df["z_min_A"], label="z min")
    plt.plot(df["timestep"], df["z_max_A"], label="z max")
    plt.xlabel("Timestep")
    plt.ylabel("Axial geometry, Å")
    plt.title(f"Dry scaffold axial geometry: {args.label}")
    plt.legend()
    plt.tight_layout()
    axial_png = outdir / f"{args.label}_axial_geometry.png"
    plt.savefig(axial_png, dpi=200)
    plt.close()

    print("Dry scaffold geometry analysis complete.")
    print(f"CSV: {csv_path}")
    print(f"Post-transient CSV: {post_path}")
    print(f"Radial PNG: {radial_png}")
    print(f"Axial PNG: {axial_png}")


if __name__ == "__main__":
    main()
