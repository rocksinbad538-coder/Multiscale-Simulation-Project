#!/usr/bin/env python3

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_lammps_dump(path):
    path = Path(path)

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
                raise RuntimeError("Expected ITEM: NUMBER OF ATOMS")
            n_atoms = int(f.readline().strip())

            line = f.readline()
            if not line.startswith("ITEM: BOX BOUNDS"):
                raise RuntimeError("Expected ITEM: BOX BOUNDS")

            bounds = []
            for _ in range(3):
                lo, hi = map(float, f.readline().split()[:2])
                bounds.append((lo, hi))

            atom_header = f.readline().strip()
            if not atom_header.startswith("ITEM: ATOMS"):
                raise RuntimeError(f"Unexpected atom header: {atom_header}")

            cols = atom_header.split()[2:]
            idx = {c: i for i, c in enumerate(cols)}

            required = ["id", "mol", "type", "x", "y", "z"]
            missing = [c for c in required if c not in idx]
            if missing:
                raise RuntimeError(f"Missing columns in dump: {missing}. Available: {cols}")

            data = np.empty((n_atoms, 6), dtype=float)

            for i in range(n_atoms):
                p = f.readline().split()
                data[i, 0] = int(p[idx["id"]])
                data[i, 1] = int(p[idx["mol"]])
                data[i, 2] = int(p[idx["type"]])
                data[i, 3] = float(p[idx["x"]])
                data[i, 4] = float(p[idx["y"]])
                data[i, 5] = float(p[idx["z"]])

            frame = pd.DataFrame(data, columns=["id", "mol", "type", "x", "y", "z"])
            frame[["id", "mol", "type"]] = frame[["id", "mol", "type"]].astype(int)

            yield timestep, frame, bounds


def compute_water_dipoles_fast(frame, oxygen_type, hydrogen_type):
    oxy = frame.loc[
        frame["type"] == oxygen_type,
        ["mol", "x", "y", "z"]
    ].copy()

    hyd = frame.loc[
        frame["type"] == hydrogen_type,
        ["mol", "x", "y", "z"]
    ].copy()

    if oxy.empty or hyd.empty:
        raise RuntimeError("No oxygen or hydrogen atoms found. Check atom type IDs.")

    # Keep only valid water molecules: one oxygen and two hydrogens per molecule.
    o_counts = oxy.groupby("mol").size()
    h_counts = hyd.groupby("mol").size()

    valid_mols = o_counts.index[
        (o_counts == 1) & (h_counts.reindex(o_counts.index).fillna(0).astype(int) == 2)
    ]

    oxy = oxy[oxy["mol"].isin(valid_mols)].copy()
    hyd = hyd[hyd["mol"].isin(valid_mols)].copy()

    # Hydrogen midpoint per molecule.
    h_mid = (
        hyd.groupby("mol", sort=False)[["x", "y", "z"]]
        .mean()
        .rename(columns={"x": "hx", "y": "hy", "z": "hz"})
        .reset_index()
    )

    oxy = oxy.rename(columns={"x": "ox", "y": "oy", "z": "oz"})
    merged = oxy.merge(h_mid, on="mol", how="inner")

    mu = merged[["hx", "hy", "hz"]].to_numpy() - merged[["ox", "oy", "oz"]].to_numpy()
    mu_norm = np.linalg.norm(mu, axis=1)

    good = mu_norm > 0.0
    merged = merged.loc[good].copy()
    mu = mu[good]
    mu_norm = mu_norm[good]

    mu_hat = mu / mu_norm[:, None]

    ox = merged["ox"].to_numpy()
    oy = merged["oy"].to_numpy()
    oz = merged["oz"].to_numpy()

    r = np.sqrt(ox**2 + oy**2)
    r_safe = np.where(r > 1.0e-12, r, np.nan)

    cos_z = mu_hat[:, 2]
    s_z = 0.5 * (3.0 * cos_z**2 - 1.0)

    cos_radial = (mu_hat[:, 0] * ox + mu_hat[:, 1] * oy) / r_safe

    out = pd.DataFrame({
        "mol": merged["mol"].to_numpy(dtype=int),
        "oxygen_x_A": ox,
        "oxygen_y_A": oy,
        "oxygen_z_A": oz,
        "oxygen_r_A": r,
        "mu_x": mu_hat[:, 0],
        "mu_y": mu_hat[:, 1],
        "mu_z": mu_hat[:, 2],
        "cos_theta_z": cos_z,
        "abs_cos_theta_z": np.abs(cos_z),
        "S_z": s_z,
        "cos_theta_radial": cos_radial,
        "abs_cos_theta_radial": np.abs(cos_radial),
    })

    return out


def selected_timesteps(values):
    values = sorted(values)
    if len(values) <= 3:
        return values
    return [values[0], values[len(values)//2], values[-1]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trajectory", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--oxygen-type", type=int, default=3)
    ap.add_argument("--hydrogen-type", type=int, default=4)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_records = []
    summary_records = []

    for iframe, (timestep, frame, _bounds) in enumerate(read_lammps_dump(args.trajectory), start=1):
        print(f"[{args.label}] processing frame {iframe}, step {timestep}", flush=True)

        dip = compute_water_dipoles_fast(frame, args.oxygen_type, args.hydrogen_type)
        dip.insert(0, "run", args.label)
        dip.insert(1, "timestep", timestep)

        all_records.append(dip)

        summary_records.append({
            "run": args.label,
            "timestep": timestep,
            "n_water": len(dip),
            "cos_theta_z_mean": dip["cos_theta_z"].mean(),
            "cos_theta_z_std": dip["cos_theta_z"].std(),
            "abs_cos_theta_z_mean": dip["abs_cos_theta_z"].mean(),
            "S_z_mean": dip["S_z"].mean(),
            "S_z_std": dip["S_z"].std(),
            "cos_theta_radial_mean": dip["cos_theta_radial"].mean(),
            "cos_theta_radial_std": dip["cos_theta_radial"].std(),
            "abs_cos_theta_radial_mean": dip["abs_cos_theta_radial"].mean(),
        })

    records = pd.concat(all_records, ignore_index=True)
    summary = pd.DataFrame(summary_records)

    records_csv = outdir / f"{args.label}_water_dipole_orientation_records.csv"
    summary_csv = outdir / f"{args.label}_water_dipole_orientation_summary.csv"

    records.to_csv(records_csv, index=False)
    summary.to_csv(summary_csv, index=False)

    chosen = selected_timesteps(summary["timestep"].unique())

    plt.figure(figsize=(8, 5))
    for ts in chosen:
        sub = records[records["timestep"] == ts]
        plt.hist(
            sub["cos_theta_z"],
            bins=60,
            density=True,
            histtype="step",
            label=f"step {ts}",
        )
    plt.xlabel("cos(theta_z)")
    plt.ylabel("Probability density")
    plt.title(f"Water dipole orientation relative to z axis: {args.label}")
    plt.legend()
    plt.tight_layout()
    cosz_png = outdir / f"{args.label}_cos_theta_z_distribution.png"
    plt.savefig(cosz_png, dpi=200)
    plt.close()

    plt.figure(figsize=(8, 5))
    for ts in chosen:
        sub = records[records["timestep"] == ts]
        plt.hist(
            sub["cos_theta_radial"].dropna(),
            bins=60,
            density=True,
            histtype="step",
            label=f"step {ts}",
        )
    plt.xlabel("cos(theta_radial)")
    plt.ylabel("Probability density")
    plt.title(f"Water dipole orientation relative to radial direction: {args.label}")
    plt.legend()
    plt.tight_layout()
    cosr_png = outdir / f"{args.label}_cos_theta_radial_distribution.png"
    plt.savefig(cosr_png, dpi=200)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(summary["timestep"], summary["S_z_mean"], label="S_z")
    plt.axhline(0.0, linestyle="--", linewidth=1)
    plt.xlabel("LAMMPS step")
    plt.ylabel("Axial dipolar order parameter S_z")
    plt.title(f"Axial water dipolar order: {args.label}")
    plt.legend()
    plt.tight_layout()
    sz_png = outdir / f"{args.label}_S_z_vs_time.png"
    plt.savefig(sz_png, dpi=200)
    plt.close()

    print("Water dipole orientation analysis complete.")
    print(f"Records CSV: {records_csv}")
    print(f"Summary CSV: {summary_csv}")
    print(f"cos(theta_z) PNG: {cosz_png}")
    print(f"cos(theta_radial) PNG: {cosr_png}")
    print(f"S_z PNG: {sz_png}")


if __name__ == "__main__":
    main()
