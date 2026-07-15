#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def pick_final(df, time_col):
    t = df[time_col].max()
    return df[df[time_col] == t].copy(), t


def main():
    outdir = Path("figures/day004")
    outdir.mkdir(parents=True, exist_ok=True)

    hydrated_density_path = Path("results/profiles/nvt300_hold_extended_contained_radial_number_density.csv")
    water_density_path = Path("results/profiles/nvt300_water_only_contained_radial_number_density.csv")

    hydrated_conf_path = Path("results/confinement/bn_like_scaffold_water_30000w_nvt300_hold_extended_contained_confinement_summary.csv")
    water_conf_path = Path("results/confinement/water_only_30000w_nvt300_contained_confinement_summary.csv")

    hydrated_density = pd.read_csv(hydrated_density_path)
    water_density = pd.read_csv(water_density_path)

    hydrated_conf = pd.read_csv(hydrated_conf_path)
    water_conf = pd.read_csv(water_conf_path)

    hyd_final, hyd_t = pick_final(hydrated_density, "timestep")
    wat_final, wat_t = pick_final(water_density, "timestep")

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    # Final radial number density
    axes[0, 0].plot(
        hyd_final["r_A"],
        hyd_final["number_density_O_per_A3"],
        label=f"hydrated scaffold-water, step {hyd_t}",
    )
    axes[0, 0].plot(
        wat_final["r_A"],
        wat_final["number_density_O_per_A3"],
        label=f"water-only, step {wat_t}",
    )
    axes[0, 0].axvline(70.0, linestyle="--", linewidth=1)
    axes[0, 0].set_title("Final volume-corrected radial O-water density")
    axes[0, 0].set_xlabel("Radial position r, Å")
    axes[0, 0].set_ylabel("O-water number density, O/Å³")
    axes[0, 0].legend()

    # Mean radial position
    axes[0, 1].plot(
        hydrated_conf["timestep"],
        hydrated_conf["water_r_mean_A"],
        label="hydrated scaffold-water",
    )
    axes[0, 1].plot(
        water_conf["timestep"],
        water_conf["water_r_mean_A"],
        label="water-only",
    )
    axes[0, 1].set_title("Mean radial position of water oxygen")
    axes[0, 1].set_xlabel("LAMMPS step")
    axes[0, 1].set_ylabel("Mean r, Å")
    axes[0, 1].legend()

    # Max radial position
    axes[1, 0].plot(
        hydrated_conf["timestep"],
        hydrated_conf["water_r_max_A"],
        label="hydrated scaffold-water",
    )
    axes[1, 0].plot(
        water_conf["timestep"],
        water_conf["water_r_max_A"],
        label="water-only",
    )
    axes[1, 0].axhline(70.0, linestyle="--", linewidth=1)
    axes[1, 0].set_title("Maximum radial position of water oxygen")
    axes[1, 0].set_xlabel("LAMMPS step")
    axes[1, 0].set_ylabel("Max r, Å")
    axes[1, 0].legend()

    # Fraction radially outside nominal lumen
    axes[1, 1].plot(
        hydrated_conf["timestep"],
        hydrated_conf["fraction_radial_outside_lumen"],
        label="hydrated scaffold-water",
    )
    axes[1, 1].plot(
        water_conf["timestep"],
        water_conf["fraction_radial_outside_lumen"],
        label="water-only",
    )
    axes[1, 1].set_title("Fraction radially outside nominal 70 Å lumen")
    axes[1, 1].set_xlabel("LAMMPS step")
    axes[1, 1].set_ylabel("Fraction")
    axes[1, 1].legend()

    fig.suptitle("Day 004 Phase 0 control comparison: scaffold-water vs water-only", fontsize=15)
    fig.tight_layout()

    png = outdir / "day004_scaffold_water_vs_water_only_control_comparison.png"
    fig.savefig(png, dpi=200)
    plt.close(fig)

    print(f"Wrote {png}")


if __name__ == "__main__":
    main()
