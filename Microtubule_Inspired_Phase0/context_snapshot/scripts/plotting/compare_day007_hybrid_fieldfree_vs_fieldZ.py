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

    hyb_ff = read_csv(
        "results/confinement/"
        "bn_like_chromophore_12dipoles_carved_nvt300_contained_extend_5k_to_20k_confinement_summary.csv"
    )
    hyb_fz = read_csv(
        "results/confinement/"
        "bn_like_chromophore_12dipoles_carved_nvt300_fieldZ_contained_10k_confinement_summary.csv"
    )

    dip_ff = read_csv(
        "results/dipoles/"
        "nvt300_contained_bn_like_chromophore_12dipoles_carved_extend_5k_to_20k_water_dipole_orientation_summary.csv"
    )
    dip_fz = read_csv(
        "results/dipoles/"
        "nvt300_fieldZ_contained_bn_like_chromophore_12dipoles_carved_10k_water_dipole_orientation_summary.csv"
    )

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    # A: max radial confinement
    axes[0, 0].plot(
        hyb_ff["timestep"],
        hyb_ff["water_r_max_A"],
        label="Hybrid field-free 5k-20k",
    )
    axes[0, 0].plot(
        hyb_fz["timestep"],
        hyb_fz["water_r_max_A"],
        label="Hybrid fieldZ 0-10k",
    )
    axes[0, 0].axhline(70.0, linestyle="--", linewidth=1, label="nominal lumen radius")
    axes[0, 0].set_title("Maximum radial water oxygen position")
    axes[0, 0].set_xlabel("step")
    axes[0, 0].set_ylabel("Max r, Å")
    axes[0, 0].legend(fontsize=8)

    # B: fraction inside
    axes[0, 1].plot(
        hyb_ff["timestep"],
        hyb_ff["fraction_inside_lumen_segment"],
        label="Hybrid field-free 5k-20k",
    )
    axes[0, 1].plot(
        hyb_fz["timestep"],
        hyb_fz["fraction_inside_lumen_segment"],
        label="Hybrid fieldZ 0-10k",
    )
    axes[0, 1].set_title("Fraction inside nominal lumen segment")
    axes[0, 1].set_xlabel("step")
    axes[0, 1].set_ylabel("Fraction inside")
    axes[0, 1].legend(fontsize=8)

    # C: mean cos(theta_z)
    axes[1, 0].plot(
        dip_ff["timestep"],
        dip_ff["cos_theta_z_mean"],
        label="Hybrid field-free 5k-20k",
    )
    axes[1, 0].plot(
        dip_fz["timestep"],
        dip_fz["cos_theta_z_mean"],
        label="Hybrid fieldZ 0-10k",
    )
    axes[1, 0].axhline(0.0, linestyle="--", linewidth=1)
    axes[1, 0].set_title("Axial water-dipole projection")
    axes[1, 0].set_xlabel("step")
    axes[1, 0].set_ylabel("mean cos(theta_z)")
    axes[1, 0].legend(fontsize=8)

    # D: S_z
    axes[1, 1].plot(
        dip_ff["timestep"],
        dip_ff["S_z_mean"],
        label="Hybrid field-free 5k-20k",
    )
    axes[1, 1].plot(
        dip_fz["timestep"],
        dip_fz["S_z_mean"],
        label="Hybrid fieldZ 0-10k",
    )
    axes[1, 1].axhline(0.0, linestyle="--", linewidth=1)
    axes[1, 1].set_title("Axial dipolar order parameter")
    axes[1, 1].set_xlabel("step")
    axes[1, 1].set_ylabel("S_z")
    axes[1, 1].legend(fontsize=8)

    fig.suptitle("Day 007 carved hybrid candidate: field-free vs fieldZ response", fontsize=15)
    fig.tight_layout()

    out = outdir / "day007_hybrid_carved_fieldfree_vs_fieldZ_response.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)

    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
