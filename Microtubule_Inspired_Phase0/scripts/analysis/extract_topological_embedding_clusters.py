#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import argparse

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
        assert handle.readline().startswith("ITEM: TIMESTEP")
        timestep = int(handle.readline().strip())

        assert handle.readline().startswith("ITEM: NUMBER OF ATOMS")
        natoms = int(handle.readline().strip())

        assert handle.readline().startswith("ITEM: BOX BOUNDS")
        for _ in range(3):
            handle.readline()

        header = handle.readline().split()
        columns = header[2:]
        col = {name: columns.index(name) for name in columns}

        required = {"id", "mol", "type", "x", "y", "z"}
        missing = required.difference(columns)
        if missing:
            raise ValueError(f"Missing columns: {sorted(missing)}")

        ids = np.empty(natoms, dtype=np.int64)
        mols = np.empty(natoms, dtype=np.int64)
        types = np.empty(natoms, dtype=np.int32)
        xyz = np.empty((natoms, 3), dtype=np.float64)

        for i in range(natoms):
            fields = handle.readline().split()
            ids[i] = int(fields[col["id"]])
            mols[i] = int(fields[col["mol"]])
            types[i] = int(fields[col["type"]])
            xyz[i] = [
                float(fields[col["x"]]),
                float(fields[col["y"]]),
                float(fields[col["z"]]),
            ]

    return Frame(timestep=timestep, ids=ids, mols=mols, types=types, xyz=xyz)


def define_pairs(frame: Frame) -> list[dict[str, object]]:
    positive_mask = frame.types == POSITIVE_SITE_TYPE
    negative_mask = frame.types == NEGATIVE_SITE_TYPE

    positive_ids = frame.ids[positive_mask]
    negative_ids = frame.ids[negative_mask]
    positive_xyz = frame.xyz[positive_mask]
    negative_xyz = frame.xyz[negative_mask]

    distances = np.linalg.norm(
        positive_xyz[:, None, :] - negative_xyz[None, :, :],
        axis=2,
    )

    unused_positive = set(range(len(positive_ids)))
    unused_negative = set(range(len(negative_ids)))
    pairs = []

    while unused_positive:
        best = None
        for ipos in unused_positive:
            for ineg in unused_negative:
                candidate = (distances[ipos, ineg], ipos, ineg)
                if best is None or candidate < best:
                    best = candidate

        _, ipos, ineg = best
        center = 0.5 * (positive_xyz[ipos] + negative_xyz[ineg])

        pairs.append(
            {
                "positive_id": int(positive_ids[ipos]),
                "negative_id": int(negative_ids[ineg]),
                "positive_xyz": positive_xyz[ipos],
                "negative_xyz": negative_xyz[ineg],
                "center": center,
            }
        )

        unused_positive.remove(ipos)
        unused_negative.remove(ineg)

    pairs.sort(
        key=lambda item: (
            float(item["center"][2]),
            float(np.arctan2(item["center"][1], item["center"][0])),
        )
    )

    for index, pair in enumerate(pairs, start=1):
        pair["pair_index"] = index

    return pairs


def load_component_atoms(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)

    atoms = pd.read_csv(path)
    atoms = atoms.sort_values(["component_index", "z_A"]).reset_index(drop=True)

    atoms["chain_position"] = atoms.groupby("component_index").cumcount()

    return atoms


def minimum_distance_to_pair(xyz, positive_xyz, negative_xyz):
    dpos = np.linalg.norm(xyz - positive_xyz[None, :], axis=1)
    dneg = np.linalg.norm(xyz - negative_xyz[None, :], axis=1)
    return np.minimum(dpos, dneg)


def write_xyz(frame: Frame, indices: np.ndarray, output: Path, comment: str) -> None:
    order = indices[np.argsort(frame.ids[indices])]

    with output.open("w", encoding="utf-8") as handle:
        handle.write(f"{len(order)}\n")
        handle.write(comment + "\n")

        for index in order:
            atom_type = int(frame.types[index])
            element = ELEMENTS.get(atom_type, "X")
            x, y, z = frame.xyz[index]
            handle.write(f"{element:2s} {x:16.8f} {y:16.8f} {z:16.8f}\n")


def write_mapping(frame: Frame, indices: np.ndarray, output: Path) -> None:
    order = indices[np.argsort(frame.ids[indices])]

    mapping = pd.DataFrame(
        {
            "cluster_atom_index": np.arange(1, len(order) + 1),
            "original_atom_id": frame.ids[order],
            "original_molecule_id": frame.mols[order],
            "atom_type": frame.types[order],
            "element": [ELEMENTS.get(int(t), "X") for t in frame.types[order]],
            "x_A": frame.xyz[order, 0],
            "y_A": frame.xyz[order, 1],
            "z_A": frame.xyz[order, 2],
        }
    )

    mapping.to_csv(output, index=False)


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument("--frame", type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--pairs", nargs="+", type=int, default=[4, 5, 7, 8])
    parser.add_argument("--n-chains", type=int, default=3)
    parser.add_argument("--half-window-sites", type=int, default=5)
    parser.add_argument("--water-cutoff-A", type=float, default=6.0)
    parser.add_argument(
        "--component-atoms",
        type=Path,
        default=Path("results/scaffold_topology/day009_scaffold_component_atoms.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("systems/hybrid/bn_like_chromophore_scaffold_water/topological_embedding_clusters"),
    )

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    frame = read_single_frame(args.frame)
    pairs = define_pairs(frame)
    pair_by_index = {int(pair["pair_index"]): pair for pair in pairs}

    component_atoms = load_component_atoms(args.component_atoms)

    id_to_frame_index = {int(atom_id): i for i, atom_id in enumerate(frame.ids)}

    summaries = []

    for pair_index in args.pairs:
        pair = pair_by_index[pair_index]

        positive_xyz = np.asarray(pair["positive_xyz"])
        negative_xyz = np.asarray(pair["negative_xyz"])

        chain_rows = []

        for component_index, group in component_atoms.groupby("component_index"):
            xyz = group[["x_A", "y_A", "z_A"]].to_numpy()
            distances = minimum_distance_to_pair(xyz, positive_xyz, negative_xyz)
            local_min_index = int(np.argmin(distances))

            chain_rows.append(
                {
                    "component_index": int(component_index),
                    "minimum_distance_A": float(distances[local_min_index]),
                    "nearest_chain_position": int(group.iloc[local_min_index]["chain_position"]),
                    "nearest_atom_id": int(group.iloc[local_min_index]["atom_id"]),
                }
            )

        chain_rank = (
            pd.DataFrame(chain_rows)
            .sort_values(["minimum_distance_A", "component_index"])
            .head(args.n_chains)
            .reset_index(drop=True)
        )

        selected_atom_ids = set([int(pair["positive_id"]), int(pair["negative_id"])])

        selected_chain_components = []

        for _, row in chain_rank.iterrows():
            component_index = int(row["component_index"])
            center_position = int(row["nearest_chain_position"])

            chain = component_atoms[
                component_atoms["component_index"].eq(component_index)
            ].copy()

            lo = max(0, center_position - args.half_window_sites)
            hi = min(len(chain) - 1, center_position + args.half_window_sites)

            segment = chain[
                chain["chain_position"].between(lo, hi)
            ]

            selected_atom_ids.update(segment["atom_id"].astype(int).tolist())
            selected_chain_components.append(component_index)

        oxygen_mask = frame.types == OXYGEN_TYPE
        oxygen_xyz = frame.xyz[oxygen_mask]
        oxygen_mols = frame.mols[oxygen_mask]

        oxygen_distance = minimum_distance_to_pair(
            oxygen_xyz,
            positive_xyz,
            negative_xyz,
        )

        selected_water_mols = set(
            oxygen_mols[oxygen_distance <= args.water_cutoff_A].astype(int)
        )

        if selected_water_mols:
            water_indices = np.where(
                np.isin(frame.mols, list(selected_water_mols))
                & np.isin(frame.types, [OXYGEN_TYPE, HYDROGEN_TYPE])
            )[0]

            selected_atom_ids.update(frame.ids[water_indices].astype(int).tolist())

        selected_indices = np.asarray(
            sorted(id_to_frame_index[atom_id] for atom_id in selected_atom_ids),
            dtype=int,
        )

        selected_types = frame.types[selected_indices]

        cluster_label = (
            f"{args.label}_pair_{pair_index:02d}_"
            f"chains{args.n_chains}_half{args.half_window_sites}_"
            f"water{args.water_cutoff_A:g}A"
        )

        xyz_path = args.output_dir / f"{cluster_label}.xyz"
        mapping_path = args.output_dir / f"{cluster_label}_mapping.csv"
        chain_rank_path = args.output_dir / f"{cluster_label}_chain_rank.csv"

        comment = (
            f"{cluster_label}; source_step={frame.timestep}; "
            f"coarse_grained_embedding_cluster; not_qm_ready"
        )

        write_xyz(frame, selected_indices, xyz_path, comment)
        write_mapping(frame, selected_indices, mapping_path)
        chain_rank.to_csv(chain_rank_path, index=False)

        summary = {
            "source_label": args.label,
            "source_frame": str(args.frame),
            "source_timestep": frame.timestep,
            "pair_index": pair_index,
            "cluster_label": cluster_label,
            "n_selected_chains": args.n_chains,
            "selected_chain_components": " ".join(str(x) for x in selected_chain_components),
            "half_window_sites": args.half_window_sites,
            "water_cutoff_A": args.water_cutoff_A,
            "n_atoms": len(selected_indices),
            "n_scaffold_atoms": int(np.isin(selected_types, [1, 2]).sum()),
            "n_water_molecules": int((selected_types == OXYGEN_TYPE).sum()),
            "n_chromophore_sites": int(np.isin(selected_types, [5, 6]).sum()),
            "xyz": str(xyz_path),
            "mapping_csv": str(mapping_path),
            "chain_rank_csv": str(chain_rank_path),
            "qm_ready": False,
            "status": "coarse_grained_embedding_cluster_not_qm_ready",
        }

        summaries.append(summary)

        print(
            f"[{args.label}] pair {pair_index:02d}: "
            f"chains={summary['selected_chain_components']}, "
            f"atoms={summary['n_atoms']}, "
            f"scaffold={summary['n_scaffold_atoms']}, "
            f"waters={summary['n_water_molecules']}"
        )

    summary_df = pd.DataFrame(summaries)
    output = args.output_dir / f"{args.label}_topological_embedding_cluster_summary.csv"
    summary_df.to_csv(output, index=False)

    print(f"\nWrote {output}")
    print(
        summary_df[
            [
                "source_label",
                "pair_index",
                "selected_chain_components",
                "n_atoms",
                "n_scaffold_atoms",
                "n_water_molecules",
                "status",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
