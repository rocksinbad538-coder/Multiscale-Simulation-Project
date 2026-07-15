#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def final_timestep(df):
    return sorted(df["timestep"].unique())[-1]


def main():
    outdir = Path("figures/day004")
    outdir.mkdir(parents=True, exist_ok=True)

    conf_scaf = pd.read_csv("results/confinement/nvt300_fieldZ_scaffold_water_confinement_summary.csv")
    conf_water = pd.read_csv("results/confinement/nvt300_fieldZ_water_only_confinement_summary.csv")

    dip_scaf = pd.read_csv("results/dipoles/nvt300_fieldZ_scaffold_water_water_dipole_orientation_summary.csv")
    dip_water = pd.read_csv("results/dipoles/nvt300_fieldZ_water_only_water_dipole_orientation_summary.csv")

    rho_scaf = pd.read_csv("results/profiles/nvt300_fieldZ_scaffold_water_radial_number_density.csv")
    rho_water = pd.read_csv("results/profiles/nvt300_fieldZ_water_only_radial_number_density.csv")

    ts_scaf = final_timestep(rho_scaf)
    ts_water = final_timestep(rho_water)

    rho_scaf_f = rho_scaf[rho_scaf["timestep"] == ts_scaf]
    rho_water_f = rho_water[rho_water["timestep"] == ts_water]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    axes[0, 0].plot(
        dip_scaf["timestep"],
        dip_scaf["S_z_mean"],
        label="scaffold-water + fieldZ",
    )
    axes[0, 0].plot(
        dip_water["timestep"],
        dip_water["S_z_mean"],
        label="water-only + fieldZ",
    )
    axes[0, 0].axhline(0.0, linestyle="--", linewidth=1)
    axes[0, 0].set_title("Axial dipolar order under fieldZ")
    axes[0, 0].set_xlabel("LAMMPS step")
    axes[0, 0].set_ylabel("S_z")
    axes[0, 0].legend()

    axes[0, 1].plot(
        dip_scaf["timestep"],
        dip_scaf["cos_theta_z_mean"],
        label="scaffold-water + fieldZ",
    )
    axes[0, 1].plot(
        dip_water["timestep"],
        dip_water["cos_theta_z_mean"],
        label="water-only + fieldZ",
    )
    axes[0, 1].axhline(0.0, linestyle="--", linewidth=1)
    axes[0, 1].set_title("Mean axial dipole projection under fieldZ")
    axes[0, 1].set_xlabel("LAMMPS step")
    axes[0, 1].set_ylabel("mean cos(theta_z)")
    axes[0, 1].legend()

    axes[1, 0].plot(
        rho_scaf_f["r_A"],
        rho_scaf_f["number_density_O_per_A3"],
        label=f"scaffold-water, step {ts_scaf}",
    )
    axes[1, 0].plot(
        rho_water_f["r_A"],
        rho_water_f["number_density_O_per_A3"],
        label=f"water-only, step {ts_water}",
    )
    axes[1, 0].axvline(70.0, linestyle="--", linewidth=1)
    axes[1, 0].set_title("Final radial O-water density under fieldZ")
    axes[1, 0].set_xlabel("Radial position r, Å")
    axes[1, 0].set_ylabel("O-water number density, O/Å³")
    axes[1, 0].legend()

    axes[1, 1].plot(
        conf_scaf["timestep"],
        conf_scaf["water_r_max_A"],
        label="scaffold-water + fieldZ",
    )
    axes[1, 1].plot(
        conf_water["timestep"],
        conf_water["water_r_max_A"],
        label="water-only + fieldZ",
    )
    axes[1, 1].axhline(70.0, linestyle="--", linewidth=1)
    axes[1, 1].set_title("Maximum radial oxygen position under fieldZ")
    axes[1, 1].set_xlabel("LAMMPS step")
    axes[1, 1].set_ylabel("Max r, Å")
    axes[1, 1].legend()

    fig.suptitle("Day 004 fieldZ response: scaffold-water vs water-only", fontsize=15)
    fig.tight_layout()

    png = outdir / "day004_fieldZ_response_scaffold_water_vs_water_only.png"
    fig.savefig(png, dpi=200)
    plt.close(fig)

    print(f"Wrote {png}")


if __name__ == "__main__":
    main()
