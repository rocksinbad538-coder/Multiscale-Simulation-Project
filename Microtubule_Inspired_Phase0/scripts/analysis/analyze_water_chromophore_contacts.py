#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


@dataclass
class Frame:
    timestep: int
    ids: np.ndarray
    mols: np.ndarray
    types: np.ndarray
    xyz: np.ndarray


def iter_lammpstrj(path: Path) -> Iterator[Frame]:
    with path.open("r", encoding="utf-8") as handle:
        while True:
            line = handle.readline()
            if not line:
                return

            if not line.startswith("ITEM: TIMESTEP"):
                continue

            timestep = int(handle.readline().strip())

            if not handle.readline().startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError("Expected NUMBER OF ATOMS section")
            natoms = int(handle.readline().strip())

            box_header = handle.readline()
            if not box_header.startswith("ITEM: BOX BOUNDS"):
                raise ValueError("Expected BOX BOUNDS section")

            for _ in range(3):
                handle.readline()

            atom_header = handle.readline().strip().split()[2:]
            required = {"id", "mol", "type", "x", "y", "z"}
            missing = required.difference(atom_header)
            if missing:
                raise ValueError(
                    f"Trajectory is missing required columns: {sorted(missing)}"
                )

            index = {name: atom_header.index(name) for name in required}

            ids = np.empty(natoms, dtype=np.int64)
            mols = np.empty(natoms, dtype=np.int64)
            types = np.empty(natoms, dtype=np.int32)
            xyz = np.empty((natoms, 3), dtype=np.float64)

            for i in range(natoms):
                fields = handle.readline().split()
                ids[i] = int(fields[index["id"]])
                mols[i] = int(fields[index["mol"]])
                types[i] = int(fields[index["type"]])
                xyz[i, 0] = float(fields[index["x"]])
                xyz[i, 1] = float(fields[index["y"]])
                xyz[i, 2] = float(fields[index["z"]])

            yield Frame(timestep, ids, mols, types, xyz)


def minimum_distances_to_sites(
    water_xyz: np.ndarray,
    site_xyz: np.ndarray,
    chunk_size: int = 5000,
) -> np.ndarray:
    if site_xyz.size == 0:
        raise ValueError("No chromophore pseudo-sites were found")

    result = np.empty(len(water_xyz), dtype=np.float64)

    for start in range(0, len(water_xyz), chunk_size):
        stop = min(start + chunk_size, len(water_xyz))
        delta = water_xyz[start:stop, None, :] - site_xyz[None, :, :]
        distance_squared = np.einsum("ijk,ijk->ij", delta, delta)
        result[start:stop] = np.sqrt(distance_squared.min(axis=1))

    return result


def continuous_run_lengths(contact_matrix: np.ndarray) -> np.ndarray:
    runs: list[int] = []

    for column in range(contact_matrix.shape[1]):
        current = 0
        for value in contact_matrix[:, column]:
            if value:
                current += 1
            elif current:
                runs.append(current)
                current = 0
        if current:
            runs.append(current)

    return np.asarray(runs, dtype=np.int64)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trajectory", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("results/contacts"))
    parser.add_argument("--oxygen-type", type=int, default=3)
    parser.add_argument("--site-types", type=int, nargs="+", default=[5, 6])
    parser.add_argument(
        "--cutoffs-A",
        type=float,
        nargs="+",
        default=[3.5, 4.5, 6.0],
    )
    parser.add_argument(
        "--frame-spacing-fs",
        type=float,
        default=50.0,
        help="Physical time between saved frames",
    )
    args = parser.parse_args()

    if not args.trajectory.exists():
        raise FileNotFoundError(args.trajectory)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    cutoffs = sorted(args.cutoffs_A)
    frame_rows: list[dict[str, float | int | str]] = []
    timesteps: list[int] = []
    molecule_ids: np.ndarray | None = None
    contact_history: dict[float, list[np.ndarray]] = {
        cutoff: [] for cutoff in cutoffs
    }

    for frame_number, frame in enumerate(iter_lammpstrj(args.trajectory), start=1):
        water_mask = frame.types == args.oxygen_type
        site_mask = np.isin(frame.types, args.site_types)

        water_xyz = frame.xyz[water_mask]
        water_mols = frame.mols[water_mask]
        site_xyz = frame.xyz[site_mask]

        if molecule_ids is None:
            order = np.argsort(water_mols)
            molecule_ids = water_mols[order]
        else:
            order = np.argsort(water_mols)
            current_ids = water_mols[order]
            if not np.array_equal(current_ids, molecule_ids):
                raise ValueError(
                    "Water molecule IDs changed between trajectory frames"
                )

        water_xyz = water_xyz[order]
        distances = minimum_distances_to_sites(water_xyz, site_xyz)

        row: dict[str, float | int | str] = {
            "label": args.label,
            "frame": frame_number,
            "timestep": frame.timestep,
            "time_from_first_frame_fs": (
                (frame_number - 1) * args.frame_spacing_fs
            ),
            "n_water": len(water_xyz),
            "n_sites": len(site_xyz),
            "minimum_water_site_distance_A": float(distances.min()),
            "mean_water_site_distance_A": float(distances.mean()),
        }

        for cutoff in cutoffs:
            contacts = distances <= cutoff
            contact_history[cutoff].append(contacts)
            row[f"n_contact_le_{cutoff:g}A"] = int(contacts.sum())
            row[f"fraction_contact_le_{cutoff:g}A"] = float(contacts.mean())

        frame_rows.append(row)
        timesteps.append(frame.timestep)

        print(
            f"[{args.label}] frame {frame_number}, step {frame.timestep}, "
            f"minimum distance = {distances.min():.3f} Å"
        )

    if molecule_ids is None or not frame_rows:
        raise RuntimeError("No trajectory frames were read")

    frame_df = pd.DataFrame(frame_rows)
    frame_csv = args.output_dir / f"{args.label}_contact_by_frame.csv"
    frame_df.to_csv(frame_csv, index=False)

    molecule_rows: list[dict[str, float | int | str]] = []
    summary_rows: list[dict[str, float | int | str]] = []

    for cutoff in cutoffs:
        matrix = np.vstack(contact_history[cutoff])
        occupancy_fraction = matrix.mean(axis=0)
        run_lengths = continuous_run_lengths(matrix)

        for index, molecule_id in enumerate(molecule_ids):
            molecule_rows.append(
                {
                    "label": args.label,
                    "water_molecule_id": int(molecule_id),
                    "cutoff_A": cutoff,
                    "contact_frame_fraction": float(occupancy_fraction[index]),
                    "contact_frames": int(matrix[:, index].sum()),
                    "total_frames": matrix.shape[0],
                }
            )

        summary_rows.append(
            {
                "label": args.label,
                "cutoff_A": cutoff,
                "n_frames": matrix.shape[0],
                "n_water": matrix.shape[1],
                "frame_spacing_fs": args.frame_spacing_fs,
                "mean_contact_count": float(matrix.sum(axis=1).mean()),
                "std_contact_count": float(matrix.sum(axis=1).std(ddof=0)),
                "mean_contact_fraction": float(matrix.mean()),
                "waters_ever_contacting": int(matrix.any(axis=0).sum()),
                "fraction_waters_ever_contacting": float(
                    matrix.any(axis=0).mean()
                ),
                "mean_occupancy_fraction_among_all_waters": float(
                    occupancy_fraction.mean()
                ),
                "mean_occupancy_fraction_among_contacting_waters": (
                    float(occupancy_fraction[occupancy_fraction > 0].mean())
                    if np.any(occupancy_fraction > 0)
                    else 0.0
                ),
                "n_continuous_contact_events": int(len(run_lengths)),
                "mean_continuous_residence_fs": (
                    float(run_lengths.mean() * args.frame_spacing_fs)
                    if len(run_lengths)
                    else 0.0
                ),
                "maximum_continuous_residence_fs": (
                    float(run_lengths.max() * args.frame_spacing_fs)
                    if len(run_lengths)
                    else 0.0
                ),
            }
        )

    molecule_df = pd.DataFrame(molecule_rows)
    summary_df = pd.DataFrame(summary_rows)

    molecule_csv = (
        args.output_dir / f"{args.label}_contact_by_water_molecule.csv"
    )
    summary_csv = args.output_dir / f"{args.label}_contact_summary.csv"

    molecule_df.to_csv(molecule_csv, index=False)
    summary_df.to_csv(summary_csv, index=False)

    fig, axis = plt.subplots(figsize=(9, 5.5))
    for cutoff in cutoffs:
        axis.plot(
            frame_df["timestep"],
            frame_df[f"n_contact_le_{cutoff:g}A"],
            marker="o",
            label=f"d ≤ {cutoff:g} Å",
        )

    axis.set_xlabel("step")
    axis.set_ylabel("Water oxygens near chromophore sites")
    axis.set_title(f"Water–chromophore contact occupancy: {args.label}")
    axis.legend()
    fig.tight_layout()

    figure = args.output_dir / f"{args.label}_contact_occupancy.png"
    fig.savefig(figure, dpi=200)
    plt.close(fig)

    print("\nContact analysis complete.")
    print(f"Frame CSV: {frame_csv}")
    print(f"Molecule CSV: {molecule_csv}")
    print(f"Summary CSV: {summary_csv}")
    print(f"Figure: {figure}")
    print("\nSummary:")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
