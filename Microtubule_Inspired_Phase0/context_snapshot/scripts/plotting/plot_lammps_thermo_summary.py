#!/usr/bin/env python3
"""
Plot Phase 0 LAMMPS thermo records.

Inputs:
- results/thermo/bn_like_scaffold_water_30000w_thermo_records.csv

Outputs:
- water temperature vs cumulative step
- water MSD vs cumulative step
- potential energy vs cumulative step
- pressure vs cumulative step
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


RUN_ORDER = ["nvt50_corrected", "nvt150", "nvt200", "nvt250", "nvt300", "nvt300_hold", "nvt300_hold_extended_contained"]


def add_cumulative_step(df: pd.DataFrame) -> pd.DataFrame:
    out = []
    offset = 0.0

    for run in RUN_ORDER:
        sub = df[df["run"] == run].copy()
        if sub.empty:
            continue
        sub["cumulative_step"] = sub["Step"] + offset
        offset = float(sub["cumulative_step"].max())
        out.append(sub)

    return pd.concat(out, ignore_index=True)


def plot_variable(df: pd.DataFrame, y: str, ylabel: str, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.8))

    for run in RUN_ORDER:
        sub = df[df["run"] == run]
        if sub.empty:
            continue
        ax.plot(sub["cumulative_step"], sub[y], label=run)

    ax.set_xlabel("Cumulative LAMMPS step")
    ax.set_ylabel(ylabel)
    ax.legend(frameon=False, fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--records-csv", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    records_csv = Path(args.records_csv)
    output_dir = Path(args.output_dir)

    df = pd.read_csv(records_csv)
    df = add_cumulative_step(df)

    plot_variable(
        df,
        y="c_twater",
        ylabel="Water temperature, K",
        output=output_dir / "bn_like_30000w_with_hold_water_temperature.png",
    )

    plot_variable(
        df,
        y="c_msdwater[4]",
        ylabel="Water MSD, Å²",
        output=output_dir / "bn_like_30000w_with_hold_water_msd.png",
    )

    plot_variable(
        df,
        y="PotEng",
        ylabel="Potential energy, kcal/mol",
        output=output_dir / "bn_like_30000w_with_hold_potential_energy.png",
    )

    plot_variable(
        df,
        y="Press",
        ylabel="Pressure, atm",
        output=output_dir / "bn_like_30000w_with_hold_pressure.png",
    )

    print("Thermo plots generated:")
    for path in sorted(output_dir.glob("bn_like_30000w_*.png")):
        print(f"  {path}")


if __name__ == "__main__":
    main()
