#!/usr/bin/env python3
"""
Analyze radial and axial confinement from LAMMPS custom dump trajectories.

Assumed dump format:
ITEM: ATOMS id type mol q x y z

Phase 0 geometry:
- Tube axis: z
- Inner/lumen radius: 70 Å
- Outer radius: 120 Å
- Nominal segment bounds: z = [-100, 100] Å

The analysis focuses on water oxygen atoms, type 3.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterator

import matplotlib.pyplot as plt
import pandas as pd


def read_lammpstrj_frames(path: Path) -> Iterator[tuple[int, list[tuple[int, int, float, float, float]]]]:
    """
    Yield frames as:
    timestep, [(id, type, x, y, z), ...]
    """
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        while True:
            line = fh.readline()
            if not line:
                break

            if not line.startswith("ITEM: TIMESTEP"):
                continue

            timestep = int(fh.readline().strip())

            line = fh.readline()
            if not line.startswith("ITEM: NUMBER OF ATOMS"):
                raise RuntimeError(f"Unexpected dump format after timestep {timestep}")

            n_atoms = int(fh.readline().strip())

            line = fh.readline()
            if not line.startswith("ITEM: BOX BOUNDS"):
                raise RuntimeError(f"Unexpected dump format before box bounds at timestep {timestep}")

            # Skip box bounds.
            fh.readline()
            fh.readline()
            fh.readline()

            atom_header = fh.readline().strip()
            if not atom_header.startswith("ITEM: ATOMS"):
                raise RuntimeError(f"Unexpected atom header at timestep {timestep}: {atom_header}")

            columns = atom_header.split()[2:]
            col_index = {name: idx for idx, name in enumerate(columns)}

            required = ["id", "type", "x", "y", "z"]
            missing = [name for name in required if name not in col_index]
            if missing:
                raise RuntimeError(f"Missing columns in dump: {missing}")

            atoms: list[tuple[int, int, float, float, float]] = []

            for _ in range(n_atoms):
                parts = fh.readline().split()
                atom_id = int(parts[col_index["id"]])
                atom_type = int(parts[col_index["type"]])
                x = float(parts[col_index["x"]])
                y = float(parts[col_index["y"]])
                z = float(parts[col_index["z"]])
                atoms.append((atom_id, atom_type, x, y, z))

            yield timestep, atoms


def summarize_frame(
    run: str,
    timestep: int,
    atoms: list[tuple[int, int, float, float, float]],
    lumen_radius_a: float,
    outer_radius_a: float,
    z_half_length_a: float,
) -> dict[str, float | str]:
    water_o = [(x, y, z) for _, typ, x, y, z in atoms if typ == 3]
    scaffold = [(x, y, z) for _, typ, x, y, z in atoms if typ in (1, 2)]

    if not water_o:
        raise RuntimeError(f"No water oxygen atoms found in run {run}, timestep {timestep}")

    water_r = [(x * x + y * y) ** 0.5 for x, y, _ in water_o]
    water_abs_z = [abs(z) for _, _, z in water_o]

    scaffold_r = [(x * x + y * y) ** 0.5 for x, y, _ in scaffold] if scaffold else []

    n_water_o = len(water_o)

    inside_lumen = [
        1
        for r, abs_z in zip(water_r, water_abs_z)
        if r <= lumen_radius_a and abs_z <= z_half_length_a
    ]

    radial_outside_lumen = [1 for r in water_r if r > lumen_radius_a]
    axial_outside_segment = [1 for abs_z in water_abs_z if abs_z > z_half_length_a]
    outside_outer = [1 for r in water_r if r > outer_radius_a]

    return {
        "run": run,
        "timestep": timestep,
        "n_water_oxygen": n_water_o,
        "water_r_mean_A": sum(water_r) / n_water_o,
        "water_r_min_A": min(water_r),
        "water_r_max_A": max(water_r),
        "water_abs_z_mean_A": sum(water_abs_z) / n_water_o,
        "water_abs_z_min_A": min(water_abs_z),
        "water_abs_z_max_A": max(water_abs_z),
        "fraction_inside_lumen_segment": len(inside_lumen) / n_water_o,
        "fraction_radial_outside_lumen": len(radial_outside_lumen) / n_water_o,
        "fraction_axial_outside_segment": len(axial_outside_segment) / n_water_o,
        "fraction_radial_outside_outer_radius": len(outside_outer) / n_water_o,
        "scaffold_r_min_A": min(scaffold_r) if scaffold_r else float("nan"),
        "scaffold_r_max_A": max(scaffold_r) if scaffold_r else float("nan"),
        "scaffold_r_mean_A": sum(scaffold_r) / len(scaffold_r) if scaffold_r else float("nan"),
    }


def write_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    if not rows:
        path.write_text("", encoding="utf-8")
        return

    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)

    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def add_cumulative_step(df: pd.DataFrame, run_order: list[str]) -> pd.DataFrame:
    parts = []
    offset = 0.0

    for run in run_order:
        sub = df[df["run"] == run].copy()
        if sub.empty:
            continue
        sub["cumulative_step"] = sub["timestep"] + offset
        offset = float(sub["cumulative_step"].max())
        parts.append(sub)

    return pd.concat(parts, ignore_index=True)


def plot_metric(df: pd.DataFrame, y: str, ylabel: str, output: Path, run_order: list[str]) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.8))

    for run in run_order:
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
    parser.add_argument("--trajectory", action="append", nargs=2, metavar=("RUN", "PATH"), required=True)
    parser.add_argument("--lumen-radius-A", type=float, default=70.0)
    parser.add_argument("--outer-radius-A", type=float, default=120.0)
    parser.add_argument("--z-half-length-A", type=float, default=100.0)
    parser.add_argument("--summary-csv", required=True)
    parser.add_argument("--figure-dir", required=True)
    args = parser.parse_args()

    rows: list[dict[str, float | str]] = []
    run_order: list[str] = []

    for run, path_str in args.trajectory:
        run_order.append(run)
        path = Path(path_str)
        print(f"Reading {run}: {path}")

        for timestep, atoms in read_lammpstrj_frames(path):
            rows.append(
                summarize_frame(
                    run=run,
                    timestep=timestep,
                    atoms=atoms,
                    lumen_radius_a=args.lumen_radius_A,
                    outer_radius_a=args.outer_radius_A,
                    z_half_length_a=args.z_half_length_A,
                )
            )

    summary_csv = Path(args.summary_csv)
    write_csv(summary_csv, rows)

    df = pd.read_csv(summary_csv)
    df = add_cumulative_step(df, run_order)

    figure_dir = Path(args.figure_dir)

    plot_metric(
        df,
        "fraction_inside_lumen_segment",
        "Fraction of water O inside lumen segment",
        figure_dir / "bn_like_30000w_fraction_inside_lumen_segment.png",
        run_order,
    )

    plot_metric(
        df,
        "fraction_radial_outside_lumen",
        "Fraction of water O with r > lumen radius",
        figure_dir / "bn_like_30000w_fraction_radial_outside_lumen.png",
        run_order,
    )

    plot_metric(
        df,
        "fraction_axial_outside_segment",
        "Fraction of water O with |z| > 100 Å",
        figure_dir / "bn_like_30000w_fraction_axial_outside_segment.png",
        run_order,
    )

    plot_metric(
        df,
        "water_r_mean_A",
        "Mean water O radial position, Å",
        figure_dir / "bn_like_30000w_water_mean_radius.png",
        run_order,
    )

    plot_metric(
        df,
        "water_abs_z_mean_A",
        "Mean |z| of water O, Å",
        figure_dir / "bn_like_30000w_water_mean_abs_z.png",
        run_order,
    )

    print("Confinement analysis complete.")
    print(f"Summary CSV: {summary_csv}")
    print(f"Figure directory: {figure_dir}")


if __name__ == "__main__":
    main()
