#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def main():
    outdir = Path("figures/day004")
    outdir.mkdir(parents=True, exist_ok=True)

    hyd_sum = pd.read_csv(
        "results/dipoles/nvt300_hold_extended_contained_water_dipole_orientation_summary.csv"
    )
    wat_sum = pd.read_csv(
        "results/dipoles/nvt300_water_only_contained_water_dipole_orientation_summary.csv"
    )

    hyd_rec = pd.read_csv(
        "results/dipoles/nvt300_hold_extended_contained_water_dipole_orientation_records.csv"
    )
    wat_rec = pd.read_csv(
        "results/dipoles/nvt300_water_only_contained_water_dipole_orientation_records.csv"
    )

    hyd_t = hyd_sum["timestep"].max()
    wat_t = wat_sum["timestep"].max()

    hyd_final = hyd_rec[hyd_rec["timestep"] == hyd_t]
    wat_final = wat_rec[wat_rec["timestep"] == wat_t]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    # S_z over time
    axes[0, 0].plot(
        hyd_sum["timestep"],
        hyd_sum["S_z_mean"],
        label="hydrated scaffold-water",
    )
    axes[0, 0].plot(
        wat_sum["timestep"],
        wat_sum["S_z_mean"],
        label="water-only",
    )
    axes[0, 0].axhline(0.0, linestyle="--", linewidth=1)
    axes[0, 0].set_title("Axial dipolar order parameter")
    axes[0, 0].set_xlabel("LAMMPS step")
    axes[0, 0].set_ylabel("S_z")
    axes[0, 0].legend()

    # Final cos(theta_z)
    axes[0, 1].hist(
        hyd_final["cos_theta_z"],
        bins=60,
        density=True,
        histtype="step",
        label=f"hydrated scaffold-water, step {hyd_t}",
    )
    axes[0, 1].hist(
        wat_final["cos_theta_z"],
        bins=60,
        density=True,
        histtype="step",
        label=f"water-only, step {wat_t}",
    )
    axes[0, 1].set_title("Final dipole orientation relative to z")
    axes[0, 1].set_xlabel("cos(theta_z)")
    axes[0, 1].set_ylabel("Probability density")
    axes[0, 1].legend()

    # Mean absolute axial projection
    axes[1, 0].plot(
        hyd_sum["timestep"],
        hyd_sum["abs_cos_theta_z_mean"],
        label="hydrated scaffold-water",
    )
    axes[1, 0].plot(
        wat_sum["timestep"],
        wat_sum["abs_cos_theta_z_mean"],
        label="water-only",
    )
    axes[1, 0].axhline(0.5, linestyle="--", linewidth=1)
    axes[1, 0].set_title("Mean absolute axial dipole projection")
    axes[1, 0].set_xlabel("LAMMPS step")
    axes[1, 0].set_ylabel("mean |cos(theta_z)|")
    axes[1, 0].legend()

    # Mean absolute radial projection
    axes[1, 1].plot(
        hyd_sum["timestep"],
        hyd_sum["abs_cos_theta_radial_mean"],
        label="hydrated scaffold-water",
    )
    axes[1, 1].plot(
        wat_sum["timestep"],
        wat_sum["abs_cos_theta_radial_mean"],
        label="water-only",
    )
    axes[1, 1].axhline(0.5, linestyle="--", linewidth=1)
    axes[1, 1].set_title("Mean absolute radial dipole projection")
    axes[1, 1].set_xlabel("LAMMPS step")
    axes[1, 1].set_ylabel("mean |cos(theta_radial)|")
    axes[1, 1].legend()

    fig.suptitle("Day 004 Phase 0 control comparison: water dipole orientation", fontsize=15)
    fig.tight_layout()

    png = outdir / "day004_scaffold_water_vs_water_only_dipole_orientation.png"
    fig.savefig(png, dpi=200)
    plt.close(fig)

    print(f"Wrote {png}")


if __name__ == "__main__":
    main()
