#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def main():
    outdir = Path("figures/day006")
    outdir.mkdir(parents=True, exist_ok=True)

    # Day 004 / Day 005 BN-like 0-20k
    bn_20k = pd.read_csv(
        "results/confinement/bn_like_scaffold_water_30000w_nvt300_hold_extended_contained_confinement_summary.csv"
    )

    # Day 006 BN-like extension 20k-50k
    bn_ext = pd.read_csv(
        "results/confinement/bn_like_scaffold_water_30000w_nvt300_hold_contained_extend_20k_to_50k_confinement_summary.csv"
    )

    # Day 005 controls, 0-20k
    bn_neut = pd.read_csv(
        "results/confinement/bn_neutralized_scaffold_water_30000w_nvt300_contained_20k_confinement_summary.csv"
    )
    carbon = pd.read_csv(
        "results/confinement/carbon_like_scaffold_water_30000w_nvt300_contained_20k_confinement_summary.csv"
    )
    water = pd.read_csv(
        "results/confinement/water_only_30000w_nvt300_contained_confinement_summary.csv"
    )

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    # Panel A: max radial water oxygen position
    axes[0, 0].plot(bn_20k["timestep"], bn_20k["water_r_max_A"], label="BN-like polar 0-20k")
    axes[0, 0].plot(bn_ext["timestep"], bn_ext["water_r_max_A"], label="BN-like polar 20k-50k")
    axes[0, 0].plot(bn_neut["timestep"], bn_neut["water_r_max_A"], label="BN-neutralized 0-20k")
    axes[0, 0].plot(carbon["timestep"], carbon["water_r_max_A"], label="carbon-like neutral 0-20k")
    axes[0, 0].plot(water["timestep"], water["water_r_max_A"], label="water-only 0-20k")
    axes[0, 0].axhline(70.0, linestyle="--", linewidth=1)
    axes[0, 0].set_title("Maximum radial water oxygen position")
    axes[0, 0].set_xlabel("LAMMPS step")
    axes[0, 0].set_ylabel("Max r, Å")
    axes[0, 0].legend(fontsize=8)

    # Panel B: mean radial position
    axes[0, 1].plot(bn_20k["timestep"], bn_20k["water_r_mean_A"], label="BN-like polar 0-20k")
    axes[0, 1].plot(bn_ext["timestep"], bn_ext["water_r_mean_A"], label="BN-like polar 20k-50k")
    axes[0, 1].plot(bn_neut["timestep"], bn_neut["water_r_mean_A"], label="BN-neutralized 0-20k")
    axes[0, 1].plot(carbon["timestep"], carbon["water_r_mean_A"], label="carbon-like neutral 0-20k")
    axes[0, 1].plot(water["timestep"], water["water_r_mean_A"], label="water-only 0-20k")
    axes[0, 1].set_title("Mean radial water oxygen position")
    axes[0, 1].set_xlabel("LAMMPS step")
    axes[0, 1].set_ylabel("Mean r, Å")
    axes[0, 1].legend(fontsize=8)

    # Panel C: fraction inside lumen segment
    axes[1, 0].plot(bn_20k["timestep"], bn_20k["fraction_inside_lumen_segment"], label="BN-like polar 0-20k")
    axes[1, 0].plot(bn_ext["timestep"], bn_ext["fraction_inside_lumen_segment"], label="BN-like polar 20k-50k")
    axes[1, 0].plot(bn_neut["timestep"], bn_neut["fraction_inside_lumen_segment"], label="BN-neutralized 0-20k")
    axes[1, 0].plot(carbon["timestep"], carbon["fraction_inside_lumen_segment"], label="carbon-like neutral 0-20k")
    axes[1, 0].plot(water["timestep"], water["fraction_inside_lumen_segment"], label="water-only 0-20k")
    axes[1, 0].set_title("Fraction inside nominal lumen segment")
    axes[1, 0].set_xlabel("LAMMPS step")
    axes[1, 0].set_ylabel("Fraction inside")
    axes[1, 0].legend(fontsize=8)

    # Panel D: radial outside-lumen fraction
    axes[1, 1].plot(bn_20k["timestep"], bn_20k["fraction_radial_outside_lumen"], label="BN-like polar 0-20k")
    axes[1, 1].plot(bn_ext["timestep"], bn_ext["fraction_radial_outside_lumen"], label="BN-like polar 20k-50k")
    axes[1, 1].plot(bn_neut["timestep"], bn_neut["fraction_radial_outside_lumen"], label="BN-neutralized 0-20k")
    axes[1, 1].plot(carbon["timestep"], carbon["fraction_radial_outside_lumen"], label="carbon-like neutral 0-20k")
    axes[1, 1].plot(water["timestep"], water["fraction_radial_outside_lumen"], label="water-only 0-20k")
    axes[1, 1].set_title("Fraction radially outside nominal lumen")
    axes[1, 1].set_xlabel("LAMMPS step")
    axes[1, 1].set_ylabel("Radial outside fraction")
    axes[1, 1].legend(fontsize=8)

    fig.suptitle("Day 006 BN-like polar scaffold-water stability extension to 50k", fontsize=15)
    fig.tight_layout()

    out = outdir / "day006_bn_like_50k_stability_vs_day005_controls.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)

    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
