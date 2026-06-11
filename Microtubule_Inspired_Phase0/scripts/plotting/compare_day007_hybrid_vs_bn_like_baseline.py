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

    # BN-like field-free baseline, 20k to 50k extension.
    bn_50k = read_csv(
        "results/confinement/"
        "bn_like_scaffold_water_30000w_nvt300_hold_contained_extend_20k_to_50k_confinement_summary.csv"
    )

    # Hybrid carved initial 0-5k.
    hyb_5k = read_csv(
        "results/confinement/"
        "bn_like_chromophore_12dipoles_carved_nvt300_contained_5k_confinement_summary.csv"
    )

    # Hybrid carved extension 5k-20k.
    hyb_20k = read_csv(
        "results/confinement/"
        "bn_like_chromophore_12dipoles_carved_nvt300_contained_extend_5k_to_20k_confinement_summary.csv"
    )

    # Dipole summaries.
    hyb_dip_5k = read_csv(
        "results/dipoles/"
        "nvt300_contained_bn_like_chromophore_12dipoles_carved_5k_water_dipole_orientation_summary.csv"
    )
    hyb_dip_20k = read_csv(
        "results/dipoles/"
        "nvt300_contained_bn_like_chromophore_12dipoles_carved_extend_5k_to_20k_water_dipole_orientation_summary.csv"
    )

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    # A: r_max
    axes[0, 0].plot(
        bn_50k["timestep"],
        bn_50k["water_r_max_A"],
        label="BN-like baseline 20k-50k",
    )
    axes[0, 0].plot(
        hyb_5k["timestep"],
        hyb_5k["water_r_max_A"],
        label="Hybrid carved 0-5k",
    )
    axes[0, 0].plot(
        hyb_20k["timestep"],
        hyb_20k["water_r_max_A"],
        label="Hybrid carved 5k-20k",
    )
    axes[0, 0].axhline(70.0, linestyle="--", linewidth=1, label="nominal lumen radius")
    axes[0, 0].set_title("Maximum radial water oxygen position")
    axes[0, 0].set_xlabel("step")
    axes[0, 0].set_ylabel("Max r, Å")
    axes[0, 0].legend(fontsize=8)

    # B: fraction inside
    axes[0, 1].plot(
        bn_50k["timestep"],
        bn_50k["fraction_inside_lumen_segment"],
        label="BN-like baseline 20k-50k",
    )
    axes[0, 1].plot(
        hyb_5k["timestep"],
        hyb_5k["fraction_inside_lumen_segment"],
        label="Hybrid carved 0-5k",
    )
    axes[0, 1].plot(
        hyb_20k["timestep"],
        hyb_20k["fraction_inside_lumen_segment"],
        label="Hybrid carved 5k-20k",
    )
    axes[0, 1].set_title("Fraction inside nominal lumen segment")
    axes[0, 1].set_xlabel("step")
    axes[0, 1].set_ylabel("Fraction inside")
    axes[0, 1].legend(fontsize=8)

    # C: mean radial position
    axes[1, 0].plot(
        bn_50k["timestep"],
        bn_50k["water_r_mean_A"],
        label="BN-like baseline 20k-50k",
    )
    axes[1, 0].plot(
        hyb_5k["timestep"],
        hyb_5k["water_r_mean_A"],
        label="Hybrid carved 0-5k",
    )
    axes[1, 0].plot(
        hyb_20k["timestep"],
        hyb_20k["water_r_mean_A"],
        label="Hybrid carved 5k-20k",
    )
    axes[1, 0].set_title("Mean radial water oxygen position")
    axes[1, 0].set_xlabel("step")
    axes[1, 0].set_ylabel("Mean r, Å")
    axes[1, 0].legend(fontsize=8)

    # D: water dipole order S_z for hybrid field-free
    axes[1, 1].plot(
        hyb_dip_5k["timestep"],
        hyb_dip_5k["S_z_mean"],
        label="Hybrid carved 0-5k",
    )
    axes[1, 1].plot(
        hyb_dip_20k["timestep"],
        hyb_dip_20k["S_z_mean"],
        label="Hybrid carved 5k-20k",
    )
    axes[1, 1].axhline(0.0, linestyle="--", linewidth=1)
    axes[1, 1].set_title("Field-free axial water-dipole order")
    axes[1, 1].set_xlabel("step")
    axes[1, 1].set_ylabel("S_z")
    axes[1, 1].legend(fontsize=8)

    fig.suptitle("Day 007 hybrid chromophore-bearing BN-like candidate vs BN-like baseline", fontsize=15)
    fig.tight_layout()

    out = outdir / "day007_hybrid_carved_20k_vs_bn_like_baseline.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)

    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
