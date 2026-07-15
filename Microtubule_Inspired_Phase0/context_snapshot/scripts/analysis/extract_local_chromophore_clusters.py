#!/usr/bin/env python3

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd


SCAFFOLD_TYPES = {1, 2}
OXYGEN_TYPE = 3
HYDROGEN_TYPE = 4
POSITIVE_SITE_TYPE = 5
NEGATIVE_SITE_TYPE = 6

ELEMENTS = {
    1: "B",
    2: "N",
    3: "O",
    4: "H",
    5: "X",
    6: "X",
}


@dataclass
class Frame:
    timestep: int
    ids: np.ndarray
    mols: np.ndarray
    types: np.ndarray
    xyz: np.ndarray


def read_single_frame(path: Path) -> Frame:
    with path.open("r", encoding="utf-8") as handle:
        line = handle.readline()

        if not line.startswith("ITEM: TIMESTEP"):
            raise ValueError(f"{path}: missing TIMESTEP header")

        timestep = int(handle.readline().strip())

        if not handle.readline().startswith("ITEM: NUMBER OF ATOMS"):
            raise ValueError(f"{path}: missing NUMBER OF ATOMS")

        natoms = int(handle.readline().strip())

        if not handle.readline().startswith("ITEM: BOX BOUNDS"):
            raise ValueError(f"{path}: missing BOX BOUNDS")

        for _ in range(3):
            handle.readline()

        atom_header = handle.readline().strip()

        if not atom_header.startswith("ITEM: ATOMS"):
            raise ValueError(f"{path}: missing ATOMS section")

        columns = atom_header.split()[2:]
        required = {"id", "mol", "type", "x", "y", "z"}
        missing = required.difference(columns)

        if missing:
            raise ValueError(
                f"{path}: missing columns {sorted(missing)}"
            )

        col = {name: columns.index(name) for name in required}

        ids = np.empty(natoms, dtype=np.int64)
        mols = np.empty(natoms, dtype=np.int64)
        types = np.empty(natoms, dtype=np.int32)
        xyz = np.empty((natoms, 3), dtype=np.float64)

        for index in range(natoms):
            fields = handle.readline().split()

            ids[index] = int(fields[col["id"]])
            mols[index] = int(fields[col["mol"]])
            types[index] = int(fields[col["type"]])

            xyz[index, 0] = float(fields[col["x"]])
            xyz[index, 1] = float(fields[col["y"]])
            xyz[index, 2] = float(fields[col["z"]])

    return Frame(
        timestep=timestep,
        ids=ids,
        mols=mols,
        types=types,
        xyz=xyz,
    )


def define_pairs(frame: Frame) -> list[dict[str, object]]:
    positive_mask = frame.types == POSITIVE_SITE_TYPE
    negative_mask = frame.types == NEGATIVE_SITE_TYPE

    positive_ids = frame.ids[positive_mask]
    negative_ids = frame.ids[negative_mask]
    positive_xyz = frame.xyz[positive_mask]
    negative_xyz = frame.xyz[negative_mask]

    if len(positive_ids) != 12 or len(negative_ids) != 12:
        raise ValueError(
            f"Expected 12 positive and 12 negative sites; found "
            f"{len(positive_ids)} and {len(negative_ids)}"
        )

    distances = np.linalg.norm(
        positive_xyz[:, None, :] - negative_xyz[None, :, :],
        axis=2,
    )

    unused_positive = set(range(12))
    unused_negative = set(range(12))
    pairs = []

    while unused_positive:
        best = None

        for ipos in unused_positive:
            for ineg in unused_negative:
                candidate = (
                    float(distances[ipos, ineg]),
                    ipos,
                    ineg,
                )

                if best is None or candidate < best:
                    best = candidate

        if best is None:
            raise RuntimeError("Chromophore pairing failed")

        distance, ipos, ineg = best
        center = 0.5 * (
            positive_xyz[ipos] + negative_xyz[ineg]
        )

        pairs.append(
            {
                "positive_id": int(positive_ids[ipos]),
                "negative_id": int(negative_ids[ineg]),
                "positive_xyz": positive_xyz[ipos],
                "negative_xyz": negative_xyz[ineg],
                "center": center,
                "pair_distance_A": distance,
            }
        )

        unused_positive.remove(ipos)
        unused_negative.remove(ineg)

    pairs.sort(
        key=lambda item: (
            float(item["center"][2]),
            float(
                np.arctan2(
                    item["center"][1],
                    item["center"][0],
                )
            ),
        )
    )

    for pair_index, pair in enumerate(pairs, start=1):
        pair["pair_index"] = pair_index

    return pairs


def minimum_distance_to_pair(
    xyz: np.ndarray,
    positive_xyz: np.ndarray,
    negative_xyz: np.ndarray,
) -> np.ndarray:
    distance_positive = np.linalg.norm(
        xyz - positive_xyz[None, :],
        axis=1,
    )
    distance_negative = np.linalg.norm(
        xyz - negative_xyz[None, :],
        axis=1,
    )

    return np.minimum(distance_positive, distance_negative)


def select_cluster(
    frame: Frame,
    pair: dict[str, object],
    water_cutoff_A: float,
    scaffold_cutoff_A: float,
) -> tuple[np.ndarray, dict[str, object]]:
    positive_id = int(pair["positive_id"])
    negative_id = int(pair["negative_id"])
    positive_xyz = np.asarray(pair["positive_xyz"])
    negative_xyz = np.asarray(pair["negative_xyz"])

    selected = np.zeros(len(frame.ids), dtype=bool)

    # Always include the two chromophore pseudo-sites.
    selected |= frame.ids == positive_id
    selected |= frame.ids == negative_id

    # Include nearby BN scaffold atoms.
    scaffold_mask = np.isin(
        frame.types,
        list(SCAFFOLD_TYPES),
    )

    scaffold_distance = minimum_distance_to_pair(
        frame.xyz,
        positive_xyz,
        negative_xyz,
    )

    selected |= (
        scaffold_mask
        & (scaffold_distance <= scaffold_cutoff_A)
    )

    # Select complete water molecules using oxygen distance.
    oxygen_mask = frame.types == OXYGEN_TYPE
    oxygen_xyz = frame.xyz[oxygen_mask]
    oxygen_mols = frame.mols[oxygen_mask]

    oxygen_distance = minimum_distance_to_pair(
        oxygen_xyz,
        positive_xyz,
        negative_xyz,
    )

    selected_water_mols = oxygen_mols[
        oxygen_distance <= water_cutoff_A
    ]

    if len(selected_water_mols):
        selected |= np.isin(frame.mols, selected_water_mols) & np.isin(
            frame.types,
            [OXYGEN_TYPE, HYDROGEN_TYPE],
        )

    selected_indices = np.where(selected)[0]

    selected_types = frame.types[selected_indices]

    summary = {
        "pair_index": int(pair["pair_index"]),
        "positive_atom_id": positive_id,
        "negative_atom_id": negative_id,
        "pair_distance_A": float(pair["pair_distance_A"]),
        "center_x_A": float(pair["center"][0]),
        "center_y_A": float(pair["center"][1]),
        "center_z_A": float(pair["center"][2]),
        "water_cutoff_A": water_cutoff_A,
        "scaffold_cutoff_A": scaffold_cutoff_A,
        "n_atoms": len(selected_indices),
        "n_scaffold_atoms": int(
            np.isin(selected_types, list(SCAFFOLD_TYPES)).sum()
        ),
        "n_chromophore_sites": int(
            np.isin(
                selected_types,
                [POSITIVE_SITE_TYPE, NEGATIVE_SITE_TYPE],
            ).sum()
        ),
        "n_water_oxygen": int(
            (selected_types == OXYGEN_TYPE).sum()
        ),
        "n_water_hydrogen": int(
            (selected_types == HYDROGEN_TYPE).sum()
        ),
        "n_water_molecules": int(
            (selected_types == OXYGEN_TYPE).sum()
        ),
    }

    return selected_indices, summary


def write_xyz(
    frame: Frame,
    indices: np.ndarray,
    output: Path,
    comment: str,
) -> None:
    order = indices[np.argsort(frame.ids[indices])]

    with output.open("w", encoding="utf-8") as handle:
        handle.write(f"{len(order)}\n")
        handle.write(comment + "\n")

        for index in order:
            atom_type = int(frame.types[index])
            element = ELEMENTS.get(atom_type, "X")
            x, y, z = frame.xyz[index]

            handle.write(
                f"{element:2s} "
                f"{x:16.8f} "
                f"{y:16.8f} "
                f"{z:16.8f}\n"
            )


def write_mapping(
    frame: Frame,
    indices: np.ndarray,
    output: Path,
) -> None:
    order = indices[np.argsort(frame.ids[indices])]

    mapping = pd.DataFrame(
        {
            "cluster_atom_index": np.arange(1, len(order) + 1),
            "original_atom_id": frame.ids[order],
            "original_molecule_id": frame.mols[order],
            "atom_type": frame.types[order],
            "element": [
                ELEMENTS.get(int(atom_type), "X")
                for atom_type in frame.types[order]
            ],
            "x_A": frame.xyz[order, 0],
            "y_A": frame.xyz[order, 1],
            "z_A": frame.xyz[order, 2],
        }
    )

    mapping.to_csv(output, index=False)


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--frame",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--label",
        required=True,
    )
    parser.add_argument(
        "--pairs",
        nargs="+",
        type=int,
        default=[4, 5, 7, 8],
    )
    parser.add_argument(
        "--water-cutoff-A",
        type=float,
        default=6.0,
    )
    parser.add_argument(
        "--scaffold-cutoff-A",
        type=float,
        default=8.0,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("systems/local_clusters"),
    )

    args = parser.parse_args()

    if not args.frame.exists():
        raise FileNotFoundError(args.frame)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    frame = read_single_frame(args.frame)
    pairs = define_pairs(frame)

    available_pairs = {
        int(pair["pair_index"]): pair
        for pair in pairs
    }

    invalid = [
        pair_index
        for pair_index in args.pairs
        if pair_index not in available_pairs
    ]

    if invalid:
        raise ValueError(f"Invalid pair indices: {invalid}")

    summaries = []

    for pair_index in args.pairs:
        pair = available_pairs[pair_index]

        indices, summary = select_cluster(
            frame=frame,
            pair=pair,
            water_cutoff_A=args.water_cutoff_A,
            scaffold_cutoff_A=args.scaffold_cutoff_A,
        )

        cluster_label = (
            f"{args.label}_pair_{pair_index:02d}_"
            f"water{args.water_cutoff_A:g}A_"
            f"scaffold{args.scaffold_cutoff_A:g}A"
        )

        xyz_path = args.output_dir / f"{cluster_label}.xyz"
        mapping_path = (
            args.output_dir / f"{cluster_label}_mapping.csv"
        )

        comment = (
            f"{cluster_label}; source_step={frame.timestep}; "
            f"pair={pair_index}; uncapped_structural_cluster"
        )

        write_xyz(
            frame=frame,
            indices=indices,
            output=xyz_path,
            comment=comment,
        )

        write_mapping(
            frame=frame,
            indices=indices,
            output=mapping_path,
        )

        summary.update(
            {
                "source_label": args.label,
                "source_frame": str(args.frame),
                "source_timestep": frame.timestep,
                "cluster_label": cluster_label,
                "xyz": str(xyz_path),
                "mapping_csv": str(mapping_path),
                "qm_ready": False,
                "status": "uncapped_structural_cluster",
            }
        )

        summaries.append(summary)

        print(
            f"[{args.label}] pair {pair_index:02d}: "
            f"{summary['n_atoms']} atoms, "
            f"{summary['n_scaffold_atoms']} BN atoms, "
            f"{summary['n_water_molecules']} waters"
        )

        print(f"  XYZ: {xyz_path}")
        print(f"  Map: {mapping_path}")

    summary_df = pd.DataFrame(summaries)

    summary_path = (
        args.output_dir
        / f"{args.label}_local_cluster_summary.csv"
    )

    summary_df.to_csv(summary_path, index=False)

    print(f"\nWrote {summary_path}")
    print("\nSummary:")
    print(
        summary_df[
            [
                "pair_index",
                "n_atoms",
                "n_scaffold_atoms",
                "n_water_molecules",
                "n_chromophore_sites",
                "qm_ready",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
