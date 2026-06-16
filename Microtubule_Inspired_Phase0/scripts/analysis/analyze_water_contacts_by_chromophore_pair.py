#!/usr/bin/env python3

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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
                raise ValueError("Expected NUMBER OF ATOMS")
            natoms = int(handle.readline().strip())

            box_header = handle.readline()
            if not box_header.startswith("ITEM: BOX BOUNDS"):
                raise ValueError("Expected BOX BOUNDS")

            for _ in range(3):
                handle.readline()

            columns = handle.readline().strip().split()[2:]
            required = {"id", "mol", "type", "x", "y", "z"}
            missing = required.difference(columns)
            if missing:
                raise ValueError(
                    f"Missing trajectory columns: {sorted(missing)}"
                )

            col = {name: columns.index(name) for name in required}

            ids = np.empty(natoms, dtype=np.int64)
            mols = np.empty(natoms, dtype=np.int64)
            types = np.empty(natoms, dtype=np.int32)
            xyz = np.empty((natoms, 3), dtype=np.float64)

            for i in range(natoms):
                fields = handle.readline().split()
                ids[i] = int(fields[col["id"]])
                mols[i] = int(fields[col["mol"]])
                types[i] = int(fields[col["type"]])
                xyz[i, 0] = float(fields[col["x"]])
                xyz[i, 1] = float(fields[col["y"]])
                xyz[i, 2] = float(fields[col["z"]])

            yield Frame(
                timestep=timestep,
                ids=ids,
                mols=mols,
                types=types,
                xyz=xyz,
            )


def greedy_minimum_pairing(
    ids_a: np.ndarray,
    xyz_a: np.ndarray,
    ids_b: np.ndarray,
    xyz_b: np.ndarray,
) -> list[tuple[int, int, float]]:
    if len(ids_a) != len(ids_b):
        raise ValueError(
            f"Unequal site counts: type A={len(ids_a)}, type B={len(ids_b)}"
        )

    distances = np.linalg.norm(
        xyz_a[:, None, :] - xyz_b[None, :, :],
        axis=2,
    )

    unused_a = set(range(len(ids_a)))
    unused_b = set(range(len(ids_b)))
    pairs: list[tuple[int, int, float]] = []

    while unused_a:
        best = None

        for ia in unused_a:
            for ib in unused_b:
                candidate = (float(distances[ia, ib]), ia, ib)
                if best is None or candidate < best:
                    best = candidate

        if best is None:
            raise RuntimeError("Pairing failed")

        distance, ia, ib = best
        pairs.append((ia, ib, distance))
        unused_a.remove(ia)
        unused_b.remove(ib)

    pairs.sort(
        key=lambda item: (
            0.5 * (xyz_a[item[0], 2] + xyz_b[item[1], 2]),
            np.arctan2(
                0.5 * (xyz_a[item[0], 1] + xyz_b[item[1], 1]),
                0.5 * (xyz_a[item[0], 0] + xyz_b[item[1], 0]),
            ),
        )
    )

    return pairs


def minimum_distance_to_two_sites(
    water_xyz: np.ndarray,
    site_a: np.ndarray,
    site_b: np.ndarray,
) -> np.ndarray:
    da = np.linalg.norm(water_xyz - site_a[None, :], axis=1)
    db = np.linalg.norm(water_xyz - site_b[None, :], axis=1)
    return np.minimum(da, db)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--trajectory",
        required=True,
        nargs="+",
        type=Path,
        help="One or more consecutive LAMMPS trajectories",
    )
    parser.add_argument("--label", required=True)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/contacts_by_pair"),
    )
    parser.add_argument("--oxygen-type", type=int, default=3)
    parser.add_argument("--positive-site-type", type=int, default=5)
    parser.add_argument("--negative-site-type", type=int, default=6)
    parser.add_argument("--cutoff-A", type=float, default=6.0)
    parser.add_argument("--frame-spacing-fs", type=float, required=True)
    args = parser.parse_args()

    for path in args.trajectory:
        if not path.exists():
            raise FileNotFoundError(path)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    unique_frames: dict[int, Frame] = {}

    for path in args.trajectory:
        print(f"Reading {path}")
        for frame in iter_lammpstrj(path):
            if frame.timestep in unique_frames:
                print(
                    f"Skipping duplicate timestep {frame.timestep} "
                    f"from {path.name}"
                )
                continue
            unique_frames[frame.timestep] = frame

    if not unique_frames:
        raise RuntimeError("No frames were read")

    frames = [unique_frames[t] for t in sorted(unique_frames)]
    reference = frames[0]

    pos_mask = reference.types == args.positive_site_type
    neg_mask = reference.types == args.negative_site_type

    pos_ids = reference.ids[pos_mask]
    neg_ids = reference.ids[neg_mask]
    pos_xyz = reference.xyz[pos_mask]
    neg_xyz = reference.xyz[neg_mask]

    if len(pos_ids) != 12 or len(neg_ids) != 12:
        raise ValueError(
            f"Expected 12+12 chromophore sites, found "
            f"{len(pos_ids)}+{len(neg_ids)}"
        )

    pairing = greedy_minimum_pairing(
        pos_ids,
        pos_xyz,
        neg_ids,
        neg_xyz,
    )

    pair_rows = []

    for pair_index, (ipos, ineg, pair_distance) in enumerate(
        pairing,
        start=1,
    ):
        center = 0.5 * (pos_xyz[ipos] + neg_xyz[ineg])
        pair_rows.append(
            {
                "pair_index": pair_index,
                "positive_atom_id": int(pos_ids[ipos]),
                "negative_atom_id": int(neg_ids[ineg]),
                "pair_distance_A": pair_distance,
                "center_x_A": center[0],
                "center_y_A": center[1],
                "center_z_A": center[2],
                "center_r_A": float(np.hypot(center[0], center[1])),
                "center_phi_rad": float(np.arctan2(center[1], center[0])),
            }
        )

    pair_definition_df = pd.DataFrame(pair_rows)
    pair_definition_csv = (
        args.output_dir / f"{args.label}_pair_definitions.csv"
    )
    pair_definition_df.to_csv(pair_definition_csv, index=False)

    print("\nChromophore-pair definitions:")
    print(pair_definition_df.to_string(index=False))

    frame_pair_rows = []
    overlap_rows = []

    reference_pos_id_to_index = {
        int(atom_id): index for index, atom_id in enumerate(pos_ids)
    }
    reference_neg_id_to_index = {
        int(atom_id): index for index, atom_id in enumerate(neg_ids)
    }

    for frame_number, frame in enumerate(frames, start=1):
        oxygen_mask = frame.types == args.oxygen_type
        water_xyz = frame.xyz[oxygen_mask]
        water_mols = frame.mols[oxygen_mask]

        id_to_xyz = {
            int(atom_id): frame.xyz[index]
            for index, atom_id in enumerate(frame.ids)
        }

        contact_matrix = np.zeros(
            (len(water_xyz), len(pairing)),
            dtype=bool,
        )

        for pair_index, (ipos, ineg, _) in enumerate(pairing):
            positive_id = int(pos_ids[ipos])
            negative_id = int(neg_ids[ineg])

            site_a = id_to_xyz[positive_id]
            site_b = id_to_xyz[negative_id]

            distances = minimum_distance_to_two_sites(
                water_xyz,
                site_a,
                site_b,
            )
            contacts = distances <= args.cutoff_A
            contact_matrix[:, pair_index] = contacts

            frame_pair_rows.append(
                {
                    "label": args.label,
                    "frame": frame_number,
                    "timestep": frame.timestep,
                    "time_from_first_frame_fs": (
                        frame_number - 1
                    ) * args.frame_spacing_fs,
                    "pair_index": pair_index + 1,
                    "positive_atom_id": positive_id,
                    "negative_atom_id": negative_id,
                    "cutoff_A": args.cutoff_A,
                    "contact_count": int(contacts.sum()),
                    "contact_fraction": float(contacts.mean()),
                    "minimum_distance_A": float(distances.min()),
                }
            )

        number_of_pairs_per_water = contact_matrix.sum(axis=1)

        overlap_rows.append(
            {
                "label": args.label,
                "frame": frame_number,
                "timestep": frame.timestep,
                "waters_contacting_any_pair": int(
                    np.any(contact_matrix, axis=1).sum()
                ),
                "waters_contacting_multiple_pairs": int(
                    (number_of_pairs_per_water > 1).sum()
                ),
                "maximum_pairs_contacted_by_one_water": int(
                    number_of_pairs_per_water.max()
                ),
            }
        )

        print(
            f"[{args.label}] frame {frame_number}/{len(frames)}, "
            f"step {frame.timestep}, "
            f"any-pair waters = "
            f"{np.any(contact_matrix, axis=1).sum()}"
        )

    frame_pair_df = pd.DataFrame(frame_pair_rows)
    overlap_df = pd.DataFrame(overlap_rows)

    pair_summary_df = (
        frame_pair_df.groupby(
            [
                "label",
                "pair_index",
                "positive_atom_id",
                "negative_atom_id",
                "cutoff_A",
            ],
            as_index=False,
        )
        .agg(
            n_frames=("frame", "count"),
            mean_contact_count=("contact_count", "mean"),
            std_contact_count=("contact_count", "std"),
            minimum_contact_count=("contact_count", "min"),
            maximum_contact_count=("contact_count", "max"),
            mean_minimum_distance_A=("minimum_distance_A", "mean"),
            minimum_observed_distance_A=("minimum_distance_A", "min"),
        )
    )

    pair_summary_df["std_contact_count"] = (
        pair_summary_df["std_contact_count"].fillna(0.0)
    )

    frame_pair_csv = (
        args.output_dir / f"{args.label}_contact_by_frame_and_pair.csv"
    )
    pair_summary_csv = (
        args.output_dir / f"{args.label}_contact_pair_summary.csv"
    )
    overlap_csv = (
        args.output_dir / f"{args.label}_contact_overlap_summary.csv"
    )

    frame_pair_df.to_csv(frame_pair_csv, index=False)
    pair_summary_df.to_csv(pair_summary_csv, index=False)
    overlap_df.to_csv(overlap_csv, index=False)

    matrix = frame_pair_df.pivot(
        index="pair_index",
        columns="timestep",
        values="contact_count",
    ).sort_index()

    fig, axis = plt.subplots(figsize=(12, 6))
    image = axis.imshow(
        matrix.to_numpy(),
        aspect="auto",
        interpolation="nearest",
        origin="lower",
    )
    axis.set_xlabel("Saved frame")
    axis.set_ylabel("Chromophore dipole-pair index")
    axis.set_title(
        f"Water occupancy within {args.cutoff_A:g} Å by chromophore pair\n"
        f"{args.label}"
    )

    x_positions = np.arange(matrix.shape[1])
    axis.set_xticks(x_positions)
    axis.set_xticklabels(
        [str(value) for value in matrix.columns],
        rotation=45,
        ha="right",
        fontsize=8,
    )

    y_positions = np.arange(matrix.shape[0])
    axis.set_yticks(y_positions)
    axis.set_yticklabels(
        [str(value) for value in matrix.index],
    )

    colorbar = fig.colorbar(image, ax=axis)
    colorbar.set_label("Water oxygen count")
    fig.tight_layout()

    heatmap_png = (
        args.output_dir / f"{args.label}_contact_pair_heatmap.png"
    )
    fig.savefig(heatmap_png, dpi=200)
    plt.close(fig)

    print("\nPair-resolved contact analysis complete.")
    print(f"Pair definitions: {pair_definition_csv}")
    print(f"Frame/pair records: {frame_pair_csv}")
    print(f"Pair summary: {pair_summary_csv}")
    print(f"Overlap summary: {overlap_csv}")
    print(f"Heatmap: {heatmap_png}")

    print("\nPair summary:")
    print(
        pair_summary_df[
            [
                "pair_index",
                "positive_atom_id",
                "negative_atom_id",
                "mean_contact_count",
                "std_contact_count",
                "minimum_contact_count",
                "maximum_contact_count",
                "minimum_observed_distance_A",
            ]
        ].to_string(index=False)
    )

    print("\nOverlap diagnostics:")
    print(overlap_df.to_string(index=False))


if __name__ == "__main__":
    main()
