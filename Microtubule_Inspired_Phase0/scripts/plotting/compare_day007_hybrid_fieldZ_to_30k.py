#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def read_csv(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    print(f"Reading {path}")
    return pd.read_csv(path)


def main():
    outdir = Path("figures/day007")
    outdir.mkdir(parents=True, exist_ok=True)

    conf_ff_5k = read_csv(
        "results/confinement/"
        "bn_like_chromophore_12dipoles_carved_nvt300_contained_5k_confinement_summary.csv"
    )
    conf_ff_20k = read_csv(
        "results/confinement/"
        "bn_like_chromophore_12dipoles_carved_nvt300_contained_extend_5k_to_20k_confinement_summary.csv"
    )
    conf_fz_10k = read_csv(
        "results/confinement/"
        "bn_like_chromophore_12dipoles_carved_nvt300_fieldZ_contained_10k_confinement_summary.csv"
    )
    conf_fz_20k = read_csv(
        "results/confinement/"
        "bn_like_chromophore_12dipoles_carved_nvt300_fieldZ_contained_extend_10k_to_20k_confinement_summary.csv"
    )
    conf_fz_30k = read_csv(
        "results/confinement/"
        "bn_like_chromophore_12dipoles_carved_nvt300_fieldZ_contained_extend_20k_to_30k_confinement_summary.csv"
    )

    dip_ff_5k = read_csv(
        "results/dipoles/"
        "nvt300_contained_bn_like_chromophore_12dipoles_carved_5k_water_dipole_orientation_summary.csv"
    )
    dip_ff_20k = read_csv(
        "results/dipoles/"
        "nvt300_contained_bn_like_chromophore_12dipoles_carved_extend_5k_to_20k_water_dipole_orientation_summary.csv"
    )
    dip_fz_10k = read_csv(
        "results/dipoles/"
        "nvt300_fieldZ_contained_bn_like_chromophore_12dipoles_carved_10k_water_dipole_orientation_summary.csv"
    )
    dip_fz_20k = read_csv(
        "results/dipoles/"
        "nvt300_fieldZ_contained_bn_like_chromophore_12dipoles_carved_extend_10k_to_20k_water_dipole_orientation_summary.csv"
    )
    dip_fz_30k = read_csv(
        "results/dipoles/"
        "nvt300_fieldZ_contained_bn_like_chromophore_12dipoles_carved_extend_20k_to_30k_water_dipole_orientation_summary.csv"
    )

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    # A: max radial confinement
    axes[0, 0].plot(conf_ff_5k["timestep"], conf_ff_5k["water_r_max_A"], label="field-free 0-5k")
    axes[0, 0].plot(conf_ff_20k["timestep"], conf_ff_20k["water_r_max_A"], label="field-free 5k-20k")
    axes[0, 0].plot(conf_fz_10k["timestep"], conf_fz_10k["water_r_max_A"], label="fieldZ 0-10k")
    axes[0, 0].plot(conf_fz_20k["timestep"], conf_fz_20k["water_r_max_A"], label="fieldZ 10k-20k")
    axes[0, 0].plot(conf_fz_30k["timestep"], conf_fz_30k["water_r_max_A"], label="fieldZ 20k-30k")
    axes[0, 0].axhline(70.0, linestyle="--", linewidth=1, label="nominal lumen radius")
    axes[0, 0].set_title("Maximum radial water oxygen position")
    axes[0, 0].set_xlabel("step")
    axes[0, 0].set_ylabel("Max r, Å")
    axes[0, 0].legend(fontsize=7)

    # B: fraction inside
    axes[0, 1].plot(conf_ff_5k["timestep"], conf_ff_5k["fraction_inside_lumen_segment"], label="field-free 0-5k")
    axes[0, 1].plot(conf_ff_20k["timestep"], conf_ff_20k["fraction_inside_lumen_segment"], label="field-free 5k-20k")
    axes[0, 1].plot(conf_fz_10k["timestep"], conf_fz_10k["fraction_inside_lumen_segment"], label="fieldZ 0-10k")
    axes[0, 1].plot(conf_fz_20k["timestep"], conf_fz_20k["fraction_inside_lumen_segment"], label="fieldZ 10k-20k")
    axes[0, 1].plot(conf_fz_30k["timestep"], conf_fz_30k["fraction_inside_lumen_segment"], label="fieldZ 20k-30k")
    axes[0, 1].set_title("Fraction inside nominal lumen segment")
    axes[0, 1].set_xlabel("step")
    axes[0, 1].set_ylabel("Fraction inside")
    axes[0, 1].legend(fontsize=7)

    # C: mean cos(theta_z)
    axes[1, 0].plot(dip_ff_5k["timestep"], dip_ff_5k["cos_theta_z_mean"], label="field-free 0-5k")
    axes[1, 0].plot(dip_ff_20k["timestep"], dip_ff_20k["cos_theta_z_mean"], label="field-free 5k-20k")
    axes[1, 0].plot(dip_fz_10k["timestep"], dip_fz_10k["cos_theta_z_mean"], label="fieldZ 0-10k")
    axes[1, 0].plot(dip_fz_20k["timestep"], dip_fz_20k["cos_theta_z_mean"], label="fieldZ 10k-20k")
    axes[1, 0].plot(dip_fz_30k["timestep"], dip_fz_30k["cos_theta_z_mean"], label="fieldZ 20k-30k")
    axes[1, 0].axhline(0.0, linestyle="--", linewidth=1)
    axes[1, 0].set_title("Axial water-dipole projection")
    axes[1, 0].set_xlabel("step")
    axes[1, 0].set_ylabel("mean cos(theta_z)")
    axes[1, 0].legend(fontsize=7)

    # D: S_z
    axes[1, 1].plot(dip_ff_5k["timestep"], dip_ff_5k["S_z_mean"], label="field-free 0-5k")
    axes[1, 1].plot(dip_ff_20k["timestep"], dip_ff_20k["S_z_mean"], label="field-free 5k-20k")
    axes[1, 1].plot(dip_fz_10k["timestep"], dip_fz_10k["S_z_mean"], label="fieldZ 0-10k")
    axes[1, 1].plot(dip_fz_20k["timestep"], dip_fz_20k["S_z_mean"], label="fieldZ 10k-20k")
    axes[1, 1].plot(dip_fz_30k["timestep"], dip_fz_30k["S_z_mean"], label="fieldZ 20k-30k")
    axes[1, 1].axhline(0.0, linestyle="--", linewidth=1)
    axes[1, 1].set_title("Axial dipolar order parameter")
    axes[1, 1].set_xlabel("step")
    axes[1, 1].set_ylabel("S_z")
    axes[1, 1].legend(fontsize=7)

    fig.suptitle("Day 007 lead hybrid candidate: field-free and fieldZ validation to 30k", fontsize=15)
    fig.tight_layout()

    out = outdir / "day007_hybrid_carved_fieldfree_and_fieldZ_to_30k.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)

    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
