#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def final_frame(df, time_col):
    t = df[time_col].max()
    return df[df[time_col] == t].copy(), t


def main():
    outdir = Path("figures/day004")
    outdir.mkdir(parents=True, exist_ok=True)

    hydrated_radial = pd.read_csv(
        "results/profiles/nvt300_hold_extended_contained_radial_number_density.csv"
    )
    water_radial = pd.read_csv(
        "results/profiles/nvt300_water_only_contained_radial_number_density.csv"
    )

    hydrated_conf = pd.read_csv(
        "results/confinement/bn_like_scaffold_water_30000w_nvt300_hold_extended_contained_confinement_summary.csv"
    )
    water_conf = pd.read_csv(
        "results/confinement/water_only_30000w_nvt300_contained_confinement_summary.csv"
    )

    dry_geom = pd.read_csv(
        "results/dry_scaffold/nvt300_dry_tethered_bn_like_scaffold_geometry_summary.csv"
    )

    hydrated_dip = pd.read_csv(
        "results/dipoles/nvt300_hold_extended_contained_water_dipole_orientation_summary.csv"
    )
    water_dip = pd.read_csv(
        "results/dipoles/nvt300_water_only_contained_water_dipole_orientation_summary.csv"
    )

    autocorr = pd.read_csv(
        "results/dipoles/day004_water_dipole_autocorrelation_comparison.csv"
    )

    h_final, h_t = final_frame(hydrated_radial, "timestep")
    w_final, w_t = final_frame(water_radial, "timestep")

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))

    # Panel 1: final radial number density
    axes[0, 0].plot(
        h_final["r_A"],
        h_final["number_density_O_per_A3"],
        label=f"hydrated scaffold-water, step {h_t}",
    )
    axes[0, 0].plot(
        w_final["r_A"],
        w_final["number_density_O_per_A3"],
        label=f"water-only, step {w_t}",
    )
    axes[0, 0].axvline(70.0, linestyle="--", linewidth=1)
    axes[0, 0].set_title("Final volume-corrected radial O-water density")
    axes[0, 0].set_xlabel("Radial position r, Å")
    axes[0, 0].set_ylabel("O-water number density, O/Å³")
    axes[0, 0].legend()

    # Panel 2: maximum radial water oxygen position
    axes[0, 1].plot(
        hydrated_conf["timestep"],
        hydrated_conf["water_r_max_A"],
        label="hydrated scaffold-water",
    )
    axes[0, 1].plot(
        water_conf["timestep"],
        water_conf["water_r_max_A"],
        label="water-only",
    )
    axes[0, 1].axhline(70.0, linestyle="--", linewidth=1)
    axes[0, 1].set_title("Maximum radial position of water oxygen")
    axes[0, 1].set_xlabel("LAMMPS step")
    axes[0, 1].set_ylabel("Max r, Å")
    axes[0, 1].legend()

    # Panel 3: dry scaffold radial geometry
    axes[1, 0].plot(dry_geom["timestep"], dry_geom["r_mean_A"], label="mean radius")
    axes[1, 0].plot(dry_geom["timestep"], dry_geom["r_min_A"], label="min radius")
    axes[1, 0].plot(dry_geom["timestep"], dry_geom["r_max_A"], label="max radius")
    axes[1, 0].axhline(70.0, linestyle="--", linewidth=1)
    axes[1, 0].axhline(120.0, linestyle="--", linewidth=1)
    axes[1, 0].set_title("Dry scaffold control: radial geometry")
    axes[1, 0].set_xlabel("LAMMPS step")
    axes[1, 0].set_ylabel("Radius, Å")
    axes[1, 0].legend()

    # Panel 4: dipole autocorrelation
    axes[1, 1].plot(
        autocorr["lag_time_ps"],
        autocorr["C_mu_hydrated_scaffold_water"],
        marker="o",
        label="hydrated scaffold-water",
    )
    axes[1, 1].plot(
        autocorr["lag_time_ps"],
        autocorr["C_mu_water_only_contained"],
        marker="o",
        label="water-only contained",
    )
    axes[1, 1].axhline(0.0, linestyle="--", linewidth=1)
    axes[1, 1].set_title("Water dipole orientational autocorrelation")
    axes[1, 1].set_xlabel("Lag time, ps")
    axes[1, 1].set_ylabel("C_mu(t)")
    axes[1, 1].legend()

    fig.suptitle(
        "Day 004 Phase 0 controls: BN-like scaffold-water, dry scaffold, and water-only references",
        fontsize=15,
    )
    fig.tight_layout()

    png = outdir / "day004_phase0_controls_consolidated_dashboard.png"
    fig.savefig(png, dpi=200)
    plt.close(fig)

    print(f"Wrote {png}")


if __name__ == "__main__":
    main()
