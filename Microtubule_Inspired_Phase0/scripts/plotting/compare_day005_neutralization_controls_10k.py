#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def truncate(df, max_step=10000):
    return df[df["timestep"] <= max_step].copy()


def final_frame(df):
    return df[df["timestep"] == df["timestep"].max()]


def main():
    outdir = Path("figures/day005")
    outdir.mkdir(parents=True, exist_ok=True)

    bn_rho = truncate(pd.read_csv("results/profiles/nvt300_hold_extended_contained_radial_number_density.csv"))
    bn_neut_rho = pd.read_csv("results/profiles/nvt300_contained_bn_neutralized_radial_number_density.csv")
    carbon_rho = truncate(pd.read_csv("results/profiles/nvt300_contained_carbon_like_20k_radial_number_density.csv"))
    water_rho = truncate(pd.read_csv("results/profiles/nvt300_water_only_contained_radial_number_density.csv"))

    bn_conf = truncate(pd.read_csv("results/confinement/bn_like_scaffold_water_30000w_nvt300_hold_extended_contained_confinement_summary.csv"))
    bn_neut_conf = pd.read_csv("results/confinement/bn_neutralized_scaffold_water_30000w_nvt300_contained_confinement_summary.csv")
    carbon_conf = truncate(pd.read_csv("results/confinement/carbon_like_scaffold_water_30000w_nvt300_contained_20k_confinement_summary.csv"))
    water_conf = truncate(pd.read_csv("results/confinement/water_only_30000w_nvt300_contained_confinement_summary.csv"))

    bn_dip = truncate(pd.read_csv("results/dipoles/nvt300_hold_extended_contained_water_dipole_orientation_summary.csv"))
    bn_neut_dip = pd.read_csv("results/dipoles/nvt300_contained_bn_neutralized_water_dipole_orientation_summary.csv")
    carbon_dip = truncate(pd.read_csv("results/dipoles/nvt300_contained_carbon_like_20k_water_dipole_orientation_summary.csv"))
    water_dip = truncate(pd.read_csv("results/dipoles/nvt300_water_only_contained_water_dipole_orientation_summary.csv"))

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    for df, label in [
        (final_frame(bn_rho), "BN-like polar"),
        (final_frame(bn_neut_rho), "BN-neutralized"),
        (final_frame(carbon_rho), "carbon-like neutral"),
        (final_frame(water_rho), "water-only"),
    ]:
        axes[0, 0].plot(df["r_A"], df["number_density_O_per_A3"], label=label)

    axes[0, 0].axvline(70.0, linestyle="--", linewidth=1)
    axes[0, 0].set_title("Final radial water oxygen density at 10k")
    axes[0, 0].set_xlabel("Radial position r, Å")
    axes[0, 0].set_ylabel("O-water number density, O/Å³")
    axes[0, 0].legend()

    for df, label in [
        (bn_conf, "BN-like polar"),
        (bn_neut_conf, "BN-neutralized"),
        (carbon_conf, "carbon-like neutral"),
        (water_conf, "water-only"),
    ]:
        axes[0, 1].plot(df["timestep"], df["water_r_max_A"], label=label)

    axes[0, 1].axhline(70.0, linestyle="--", linewidth=1)
    axes[0, 1].set_title("Maximum radial water oxygen position")
    axes[0, 1].set_xlabel("LAMMPS step")
    axes[0, 1].set_ylabel("Max r, Å")
    axes[0, 1].legend()

    for df, label in [
        (bn_dip, "BN-like polar"),
        (bn_neut_dip, "BN-neutralized"),
        (carbon_dip, "carbon-like neutral"),
        (water_dip, "water-only"),
    ]:
        axes[1, 0].plot(df["timestep"], df["S_z_mean"], label=label)

    axes[1, 0].axhline(0.0, linestyle="--", linewidth=1)
    axes[1, 0].set_title("Axial water dipolar order without applied field")
    axes[1, 0].set_xlabel("LAMMPS step")
    axes[1, 0].set_ylabel("S_z")
    axes[1, 0].legend()

    for df, label in [
        (bn_dip, "BN-like polar"),
        (bn_neut_dip, "BN-neutralized"),
        (carbon_dip, "carbon-like neutral"),
        (water_dip, "water-only"),
    ]:
        axes[1, 1].plot(df["timestep"], df["abs_cos_theta_z_mean"], label=label)

    axes[1, 1].axhline(0.5, linestyle="--", linewidth=1)
    axes[1, 1].set_title("Mean absolute axial dipole projection")
    axes[1, 1].set_xlabel("LAMMPS step")
    axes[1, 1].set_ylabel("mean |cos(theta_z)|")
    axes[1, 1].legend()

    fig.suptitle("Day 005 charge/polarity control at 10k: BN-like vs BN-neutralized vs carbon-like vs water-only", fontsize=14)
    fig.tight_layout()

    out = outdir / "day005_charge_polarity_controls_10k_comparison.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
