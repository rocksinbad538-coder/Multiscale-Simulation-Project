#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def final_frame(df):
    return df[df["timestep"] == df["timestep"].max()]


def main():
    outdir = Path("figures/day005")
    outdir.mkdir(parents=True, exist_ok=True)

    # FieldZ radial density
    bn_rho = pd.read_csv("results/profiles/nvt300_fieldZ_scaffold_water_radial_number_density.csv")
    bn_neut_rho = pd.read_csv("results/profiles/nvt300_fieldZ_bn_neutralized_radial_number_density.csv")
    carbon_rho = pd.read_csv("results/profiles/nvt300_fieldZ_carbon_like_radial_number_density.csv")
    water_rho = pd.read_csv("results/profiles/nvt300_fieldZ_water_only_radial_number_density.csv")

    # FieldZ confinement
    bn_conf = pd.read_csv("results/confinement/nvt300_fieldZ_scaffold_water_confinement_summary.csv")
    bn_neut_conf = pd.read_csv("results/confinement/bn_neutralized_scaffold_water_30000w_fieldZ_confinement_summary.csv")
    carbon_conf = pd.read_csv("results/confinement/carbon_like_scaffold_water_30000w_fieldZ_confinement_summary.csv")
    water_conf = pd.read_csv("results/confinement/nvt300_fieldZ_water_only_confinement_summary.csv")

    # FieldZ dipoles
    bn_dip = pd.read_csv("results/dipoles/nvt300_fieldZ_scaffold_water_water_dipole_orientation_summary.csv")
    bn_neut_dip = pd.read_csv("results/dipoles/nvt300_fieldZ_bn_neutralized_water_dipole_orientation_summary.csv")
    carbon_dip = pd.read_csv("results/dipoles/nvt300_fieldZ_carbon_like_water_dipole_orientation_summary.csv")
    water_dip = pd.read_csv("results/dipoles/nvt300_fieldZ_water_only_water_dipole_orientation_summary.csv")

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    for df, label in [
        (final_frame(bn_rho), "BN-like polar + fieldZ"),
        (final_frame(bn_neut_rho), "BN-neutralized + fieldZ"),
        (final_frame(carbon_rho), "carbon-like neutral + fieldZ"),
        (final_frame(water_rho), "water-only + fieldZ"),
    ]:
        axes[0, 0].plot(df["r_A"], df["number_density_O_per_A3"], label=label)

    axes[0, 0].axvline(70.0, linestyle="--", linewidth=1)
    axes[0, 0].set_title("Final radial water oxygen density under fieldZ")
    axes[0, 0].set_xlabel("Radial position r, Å")
    axes[0, 0].set_ylabel("O-water number density, O/Å³")
    axes[0, 0].legend(fontsize=8)

    for df, label in [
        (bn_conf, "BN-like polar + fieldZ"),
        (bn_neut_conf, "BN-neutralized + fieldZ"),
        (carbon_conf, "carbon-like neutral + fieldZ"),
        (water_conf, "water-only + fieldZ"),
    ]:
        axes[0, 1].plot(df["timestep"], df["water_r_max_A"], label=label)

    axes[0, 1].axhline(70.0, linestyle="--", linewidth=1)
    axes[0, 1].set_title("Maximum radial water oxygen position under fieldZ")
    axes[0, 1].set_xlabel("LAMMPS step")
    axes[0, 1].set_ylabel("Max r, Å")
    axes[0, 1].legend(fontsize=8)

    for df, label in [
        (bn_dip, "BN-like polar + fieldZ"),
        (bn_neut_dip, "BN-neutralized + fieldZ"),
        (carbon_dip, "carbon-like neutral + fieldZ"),
        (water_dip, "water-only + fieldZ"),
    ]:
        axes[1, 0].plot(df["timestep"], df["S_z_mean"], label=label)

    axes[1, 0].axhline(0.0, linestyle="--", linewidth=1)
    axes[1, 0].set_title("Axial water dipolar order under fieldZ")
    axes[1, 0].set_xlabel("LAMMPS step")
    axes[1, 0].set_ylabel("S_z")
    axes[1, 0].legend(fontsize=8)

    for df, label in [
        (bn_dip, "BN-like polar + fieldZ"),
        (bn_neut_dip, "BN-neutralized + fieldZ"),
        (carbon_dip, "carbon-like neutral + fieldZ"),
        (water_dip, "water-only + fieldZ"),
    ]:
        axes[1, 1].plot(df["timestep"], df["cos_theta_z_mean"], label=label)

    axes[1, 1].axhline(0.0, linestyle="--", linewidth=1)
    axes[1, 1].set_title("Mean axial dipole projection under fieldZ")
    axes[1, 1].set_xlabel("LAMMPS step")
    axes[1, 1].set_ylabel("mean cos(theta_z)")
    axes[1, 1].legend(fontsize=8)

    fig.suptitle("Day 005 fieldZ charge/polarity controls", fontsize=15)
    fig.tight_layout()

    out = outdir / "day005_fieldZ_charge_polarity_controls_comparison.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)

    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
