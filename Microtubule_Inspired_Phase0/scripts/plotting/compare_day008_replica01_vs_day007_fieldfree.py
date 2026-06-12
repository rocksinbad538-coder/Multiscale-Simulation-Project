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
    outdir = Path("figures/day008")
    outdir.mkdir(parents=True, exist_ok=True)

    conf_day007 = read_csv(
        "results/confinement/"
        "bn_like_chromophore_12dipoles_carved_nvt300_contained_extend_5k_to_20k_confinement_summary.csv"
    )
    conf_rep01 = read_csv(
        "results/confinement/"
        "bn_like_chromophore_12dipoles_carved_nvt300_contained_replica01_20k_confinement_summary.csv"
    )

    dip_day007 = read_csv(
        "results/dipoles/"
        "nvt300_contained_bn_like_chromophore_12dipoles_carved_extend_5k_to_20k_water_dipole_orientation_summary.csv"
    )
    dip_rep01 = read_csv(
        "results/dipoles/"
        "nvt300_contained_bn_like_chromophore_12dipoles_carved_replica01_20k_water_dipole_orientation_summary.csv"
    )

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    axes[0, 0].plot(conf_day007["timestep"], conf_day007["water_r_max_A"], label="Day 007 field-free 5k-20k")
    axes[0, 0].plot(conf_rep01["timestep"], conf_rep01["water_r_max_A"], label="Day 008 replica01 field-free 0-20k")
    axes[0, 0].axhline(70.0, linestyle="--", linewidth=1, label="nominal lumen radius")
    axes[0, 0].set_title("Maximum radial water oxygen position")
    axes[0, 0].set_xlabel("step")
    axes[0, 0].set_ylabel("Max r, Å")
    axes[0, 0].legend(fontsize=8)

    axes[0, 1].plot(conf_day007["timestep"], conf_day007["fraction_inside_lumen_segment"], label="Day 007 field-free")
    axes[0, 1].plot(conf_rep01["timestep"], conf_rep01["fraction_inside_lumen_segment"], label="Day 008 replica01")
    axes[0, 1].set_title("Fraction inside nominal lumen segment")
    axes[0, 1].set_xlabel("step")
    axes[0, 1].set_ylabel("Fraction inside")
    axes[0, 1].legend(fontsize=8)

    axes[1, 0].plot(dip_day007["timestep"], dip_day007["cos_theta_z_mean"], label="Day 007 field-free")
    axes[1, 0].plot(dip_rep01["timestep"], dip_rep01["cos_theta_z_mean"], label="Day 008 replica01")
    axes[1, 0].axhline(0.0, linestyle="--", linewidth=1)
    axes[1, 0].set_title("Field-free axial water-dipole projection")
    axes[1, 0].set_xlabel("step")
    axes[1, 0].set_ylabel("mean cos(theta_z)")
    axes[1, 0].legend(fontsize=8)

    axes[1, 1].plot(dip_day007["timestep"], dip_day007["S_z_mean"], label="Day 007 field-free")
    axes[1, 1].plot(dip_rep01["timestep"], dip_rep01["S_z_mean"], label="Day 008 replica01")
    axes[1, 1].axhline(0.0, linestyle="--", linewidth=1)
    axes[1, 1].set_title("Field-free axial dipolar order parameter")
    axes[1, 1].set_xlabel("step")
    axes[1, 1].set_ylabel("S_z")
    axes[1, 1].legend(fontsize=8)

    fig.suptitle("Day 008 field-free reproducibility: lead hybrid candidate", fontsize=15)
    fig.tight_layout()

    out = outdir / "day008_hybrid_fieldfree_replica01_vs_day007.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)

    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
