#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def main():
    outdir = Path("figures/day004")
    outdir.mkdir(parents=True, exist_ok=True)

    hydrated_path = Path("results/dipoles/nvt300_hold_extended_contained_water_dipole_autocorrelation.csv")
    water_path = Path("results/dipoles/nvt300_water_only_contained_water_dipole_autocorrelation.csv")

    hydrated = pd.read_csv(hydrated_path)
    water = pd.read_csv(water_path)

    merged = hydrated[["lag_time_ps", "C1_mu"]].rename(
        columns={"C1_mu": "C_mu_hydrated_scaffold_water"}
    ).merge(
        water[["lag_time_ps", "C1_mu"]].rename(
            columns={"C1_mu": "C_mu_water_only_contained"}
        ),
        on="lag_time_ps",
        how="inner",
    )

    merged["delta_C_mu_hydrated_minus_water_only"] = (
        merged["C_mu_hydrated_scaffold_water"]
        - merged["C_mu_water_only_contained"]
    )

    out_csv = Path("results/dipoles/day004_water_dipole_autocorrelation_comparison.csv")
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_csv, index=False)

    plt.figure(figsize=(8, 5))
    plt.plot(
        merged["lag_time_ps"],
        merged["C_mu_hydrated_scaffold_water"],
        marker="o",
        label="hydrated scaffold-water",
    )
    plt.plot(
        merged["lag_time_ps"],
        merged["C_mu_water_only_contained"],
        marker="o",
        label="water-only contained",
    )
    plt.axhline(0.0, linestyle="--", linewidth=1)
    plt.xlabel("Lag time, ps")
    plt.ylabel("Dipole orientational autocorrelation")
    plt.title("Day 004 water dipole autocorrelation comparison")
    plt.legend()
    plt.tight_layout()

    out_png = outdir / "day004_water_dipole_autocorrelation_comparison.png"
    plt.savefig(out_png, dpi=200)
    plt.close()

    print(f"Wrote {out_csv}")
    print(f"Wrote {out_png}")
    print(merged.tail().to_string(index=False))


if __name__ == "__main__":
    main()
