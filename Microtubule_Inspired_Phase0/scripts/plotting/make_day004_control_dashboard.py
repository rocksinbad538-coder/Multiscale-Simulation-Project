#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def main():
    outdir = Path("figures/day004")
    outdir.mkdir(parents=True, exist_ok=True)

    thermo_path = Path("results/thermo/bn_like_scaffold_water_30000w_thermo_records_with_extended_contained.csv")
    conf_path = Path("results/confinement/bn_like_scaffold_water_30000w_nvt300_hold_extended_contained_confinement_summary.csv")
    dry_geom_path = Path("results/dry_scaffold/nvt300_dry_tethered_bn_like_scaffold_geometry_summary.csv")
    dry_thermo_path = Path("results/thermo/bn_like_scaffold_dry_tethered_300K_thermo_records.csv")

    thermo = pd.read_csv(thermo_path)
    conf = pd.read_csv(conf_path)
    dry_geom = pd.read_csv(dry_geom_path)
    dry_thermo = pd.read_csv(dry_thermo_path)

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    # Hydrated extended contained water temperature
    sub = thermo[thermo["run"] == "nvt300_hold_extended_contained"].copy()
    axes[0, 0].plot(sub["Step"], sub["c_twater"])
    axes[0, 0].axhline(300.0, linestyle="--", linewidth=1)
    axes[0, 0].set_title("Hydrated extended contained hold: water temperature")
    axes[0, 0].set_xlabel("LAMMPS step")
    axes[0, 0].set_ylabel("Water temperature, K")

    # Hydrated extended contained confinement
    axes[0, 1].plot(conf["timestep"], conf["fraction_inside_lumen_segment"], label="inside lumen segment")
    axes[0, 1].plot(conf["timestep"], conf["fraction_radial_outside_lumen"], label="radially outside lumen")
    axes[0, 1].plot(conf["timestep"], conf["fraction_axial_outside_segment"], label="axially outside segment")
    axes[0, 1].set_title("Hydrated extended contained hold: water confinement")
    axes[0, 1].set_xlabel("LAMMPS step")
    axes[0, 1].set_ylabel("Fraction")
    axes[0, 1].legend()

    # Dry scaffold radial geometry
    axes[1, 0].plot(dry_geom["timestep"], dry_geom["r_mean_A"], label="mean radius")
    axes[1, 0].plot(dry_geom["timestep"], dry_geom["r_min_A"], label="min radius")
    axes[1, 0].plot(dry_geom["timestep"], dry_geom["r_max_A"], label="max radius")
    axes[1, 0].axhline(70.0, linestyle="--", linewidth=1)
    axes[1, 0].axhline(120.0, linestyle="--", linewidth=1)
    axes[1, 0].set_title("Dry scaffold control: radial geometry")
    axes[1, 0].set_xlabel("LAMMPS step")
    axes[1, 0].set_ylabel("Radius, Å")
    axes[1, 0].legend()

    # Dry scaffold temperature
    axes[1, 1].plot(dry_thermo["Step"], dry_thermo["Temp"])
    axes[1, 1].axhline(300.0, linestyle="--", linewidth=1)
    axes[1, 1].set_title("Dry scaffold control: tethered temperature")
    axes[1, 1].set_xlabel("LAMMPS step")
    axes[1, 1].set_ylabel("Temperature, K")

    fig.suptitle("Day 004 Phase 0 controls: hydrated baseline and dry scaffold reference", fontsize=15)
    fig.tight_layout()

    png = outdir / "day004_bn_like_hydrated_and_dry_controls_dashboard.png"
    fig.savefig(png, dpi=200)
    plt.close(fig)

    print(f"Wrote {png}")


if __name__ == "__main__":
    main()
