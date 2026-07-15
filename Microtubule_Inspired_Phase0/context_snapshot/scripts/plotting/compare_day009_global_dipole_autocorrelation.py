#!/usr/bin/env python3

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


FILES = {
    "Day 007 field-free":
        "results/dipole_autocorrelation/"
        "day007_fieldfree_dipole_autocorrelation.csv",
    "Day 008 field-free replica01":
        "results/dipole_autocorrelation/"
        "day008_fieldfree_replica01_dipole_autocorrelation.csv",
    "Day 007 fieldZ":
        "results/dipole_autocorrelation/"
        "day007_fieldZ_0_to_20k_dipole_autocorrelation.csv",
    "Day 008 fieldZ replica01":
        "results/dipole_autocorrelation/"
        "day008_fieldZ_replica01_0_to_20k_dipole_autocorrelation.csv",
}


def main() -> None:
    output_dir = Path("figures/day009")
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, axis = plt.subplots(figsize=(10, 6.5))

    for label, filename in FILES.items():
        path = Path(filename)

        if not path.exists():
            raise FileNotFoundError(path)

        df = pd.read_csv(path)

        axis.plot(
            df["lag_ps"],
            df["C1_connected"],
            marker="o",
            label=label,
        )

    axis.axhline(0.0, linewidth=1)
    axis.set_xlabel("Lag time, ps")
    axis.set_ylabel("Connected dipole autocorrelation")
    axis.set_title(
        "Day 009 — Global water-dipole relaxation\n"
        "field-free and fieldZ reproducibility"
    )
    axis.legend(fontsize=8)
    fig.tight_layout()

    output = (
        output_dir
        / "day009_global_connected_dipole_autocorrelation.png"
    )
    fig.savefig(output, dpi=200)
    plt.close(fig)

    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
