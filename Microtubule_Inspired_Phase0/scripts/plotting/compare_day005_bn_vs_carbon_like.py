#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def final_frame(df):
    return df[df["timestep"] == df["timestep"].max()]


def main():
    outdir = Path("figures/day005")
    outdir.mkdir(parents=True, exist_ok=True)

    bn_rho = pd.read_csv("results/profiles/nvt300_hold_extended_contained_radial_number_density.csv")
    carbon_rho = pd.read_csv("results/profiles/nvt300_contained_carbon_like_20k_radial_number_density.csv")
    water_rho = pd.read_csv("results/profiles/nvt300_water_only_contained_radial_number_density.csv")

    bn_conf = pd.read_csv("results/confinement/bn_like_scaffold_water_30000w_nvt300_hold_extended_contained_confinement_summary.csv")
    carbon_conf = pd.read_csv("results/confinement/carbon_like_scaffold_water_30000w_nvt300_contained_20k_confinement_summary.csv")
    water_conf = pd.read_csv("results/confinement/water_only_30000w_nvt300_contained_confinement_summary.csv")

    bn_dip = pd.read_csv("results/dipoles/nvt300_hold_extended_contained_water_dipole_orientation_summary.csv")
    carbon_dip = pd.read_csv("results/dipoles/nvt300_contained_carbon_like_20k_water_dipole_orientation_summary.csv")
    water_dip = pd.read_csv("results/dipoles/nvt300_water_only_contained_water_dipole_orientation_summary.csv")

    bn_rho_f = final_frame(bn_rho)
    carbon_rho_f = final_frame(carbon_rho)
    water_rho_f = final_frame(water_rho)

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    axes[0, 0].plot(bn_rho_f["r_A"], bn_rho_f["number_density_O_per_A3"], label="BN-like polar scaffold")
    axes[0, 0].plot(carbon_rho_f["r_A"], carbon_rho_f["number_density_O_per_A3"], label="carbon-like neutral scaffold")
    axes[0, 0].plot(water_rho_f["r_A"], water_rho_f["number_density_O_per_A3"], label="water-only")
    axes[0, 0].axvline(70.0, linestyle="--", linewidth=1)
    axes[0, 0].set_title("Final radial water oxygen density")
    axes[0, 0].set_xlabel("Radial position r, Å")
    axes[0, 0].set_ylabel("O-water number density, O/Å³")
    axes[0, 0].legend()

    axes[0, 1].plot(bn_conf["timestep"], bn_conf["water_r_max_A"], label="BN-like polar scaffold")
    axes[0, 1].plot(carbon_conf["timestep"], carbon_conf["water_r_max_A"], label="carbon-like neutral scaffold")
    axes[0, 1].plot(water_conf["timestep"], water_conf["water_r_max_A"], label="water-only")
    axes[0, 1].axhline(70.0, linestyle="--", linewidth=1)
    axes[0, 1].set_title("Maximum radial water oxygen position")
    axes[0, 1].set_xlabel("LAMMPS step")
    axes[0, 1].set_ylabel("Max r, Å")
    axes[0, 1].legend()

    axes[1, 0].plot(bn_dip["timestep"], bn_dip["S_z_mean"], label="BN-like polar scaffold")
    axes[1, 0].plot(carbon_dip["timestep"], carbon_dip["S_z_mean"], label="carbon-like neutral scaffold")
    axes[1, 0].plot(water_dip["timestep"], water_dip["S_z_mean"], label="water-only")
    axes[1, 0].axhline(0.0, linestyle="--", linewidth=1)
    axes[1, 0].set_title("Axial water dipolar order without applied field")
    axes[1, 0].set_xlabel("LAMMPS step")
    axes[1, 0].set_ylabel("S_z")
    axes[1, 0].legend()

    axes[1, 1].plot(bn_dip["timestep"], bn_dip["abs_cos_theta_z_mean"], label="BN-like polar scaffold")
    axes[1, 1].plot(carbon_dip["timestep"], carbon_dip["abs_cos_theta_z_mean"], label="carbon-like neutral scaffold")
    axes[1, 1].plot(water_dip["timestep"], water_dip["abs_cos_theta_z_mean"], label="water-only")
    axes[1, 1].axhline(0.5, linestyle="--", linewidth=1)
    axes[1, 1].set_title("Mean absolute axial dipole projection")
    axes[1, 1].set_xlabel("LAMMPS step")
    axes[1, 1].set_ylabel("mean |cos(theta_z)|")
    axes[1, 1].legend()

    fig.suptitle("Day 005 Phase 0 matched-time comparison: BN-like vs carbon-like vs water-only", fontsize=15)
    fig.tight_layout()

    out = outdir / "day005_matched_time_bn_like_vs_carbon_like_vs_water_only_comparison.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
