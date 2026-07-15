#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


SCAFFOLD_TYPES = {1, 2}


def read_lammpstrj_frame(path: Path) -> dict[str, np.ndarray | int]:
    with path.open("r", encoding="utf-8") as handle:
        if not handle.readline().startswith("ITEM: TIMESTEP"):
            raise ValueError(f"{path}: missing TIMESTEP")

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
            raise ValueError(f"{path}: missing ATOMS")

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

    return {
        "timestep": timestep,
        "ids": ids,
        "mols": mols,
        "types": types,
        "xyz": xyz,
    }


def connected_components(
    atom_ids: list[int],
    adjacency: dict[int, set[int]],
) -> list[list[int]]:
    unvisited = set(atom_ids)
    components: list[list[int]] = []

    while unvisited:
        start = min(unvisited)
        stack = [start]
        component: list[int] = []

        while stack:
            current = stack.pop()

            if current not in unvisited:
                continue

            unvisited.remove(current)
            component.append(current)

            for neighbor in adjacency.get(current, set()):
                if neighbor in unvisited:
                    stack.append(neighbor)

        components.append(sorted(component))

    components.sort(key=lambda component: (-len(component), component[0]))

    return components


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--inventory",
        type=Path,
        default=Path(
            "systems/hybrid/bn_like_chromophore_scaffold_water/"
            "local_clusters/day009_local_cluster_inventory.csv"
        ),
    )
    parser.add_argument(
        "--bond-cutoff-A",
        type=float,
        default=1.90,
        help="Maximum distance used to identify a B-N scaffold bond",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "results/local_cluster_connectivity"
        ),
    )

    args = parser.parse_args()

    if not args.inventory.exists():
        raise FileNotFoundError(args.inventory)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    inventory = pd.read_csv(args.inventory)

    required = {
        "cluster_label",
        "source_frame",
        "mapping_csv",
        "source_label",
        "pair_index",
    }
    missing = required.difference(inventory.columns)

    if missing:
        raise ValueError(
            f"Inventory is missing columns: {sorted(missing)}"
        )

    all_cluster_rows: list[dict[str, object]] = []
    all_atom_rows: list[dict[str, object]] = []

    frame_cache: dict[str, dict[str, np.ndarray | int]] = {}

    for _, inventory_row in inventory.iterrows():
        cluster_label = str(inventory_row["cluster_label"])
        source_frame_path = Path(str(inventory_row["source_frame"]))
        mapping_path = Path(str(inventory_row["mapping_csv"]))

        if not source_frame_path.exists():
            raise FileNotFoundError(source_frame_path)

        if not mapping_path.exists():
            raise FileNotFoundError(mapping_path)

        source_key = str(source_frame_path)

        if source_key not in frame_cache:
            frame_cache[source_key] = read_lammpstrj_frame(
                source_frame_path
            )

        frame = frame_cache[source_key]
        mapping = pd.read_csv(mapping_path)

        cluster_bn = mapping[
            mapping["atom_type"].isin(SCAFFOLD_TYPES)
        ].copy()

        cluster_bn_ids = set(
            cluster_bn["original_atom_id"].astype(int)
        )

        frame_ids = np.asarray(frame["ids"])
        frame_types = np.asarray(frame["types"])
        frame_xyz = np.asarray(frame["xyz"])

        full_bn_mask = np.isin(
            frame_types,
            list(SCAFFOLD_TYPES),
        )

        full_bn_ids = frame_ids[full_bn_mask]
        full_bn_types = frame_types[full_bn_mask]
        full_bn_xyz = frame_xyz[full_bn_mask]

        id_to_full_index = {
            int(atom_id): index
            for index, atom_id in enumerate(full_bn_ids)
        }

        tree = cKDTree(full_bn_xyz)
        neighbor_pairs = tree.query_pairs(
            r=args.bond_cutoff_A,
            output_type="set",
        )

        full_adjacency: dict[int, set[int]] = {
            int(atom_id): set()
            for atom_id in full_bn_ids
        }

        bond_lengths: dict[tuple[int, int], float] = {}

        for index_a, index_b in neighbor_pairs:
            atom_a = int(full_bn_ids[index_a])
            atom_b = int(full_bn_ids[index_b])

            type_a = int(full_bn_types[index_a])
            type_b = int(full_bn_types[index_b])

            # BN scaffold bonds must connect unlike atom types.
            if type_a == type_b:
                continue

            distance = float(
                np.linalg.norm(
                    full_bn_xyz[index_a] - full_bn_xyz[index_b]
                )
            )

            full_adjacency[atom_a].add(atom_b)
            full_adjacency[atom_b].add(atom_a)

            key = tuple(sorted((atom_a, atom_b)))
            bond_lengths[key] = distance

        internal_adjacency: dict[int, set[int]] = {
            atom_id: set()
            for atom_id in cluster_bn_ids
        }

        internal_lengths: list[float] = []
        cut_lengths: list[float] = []

        total_cut_bonds = 0
        n_undercoordinated_atoms = 0
        n_full_coordination_anomalies = 0

        for atom_id in sorted(cluster_bn_ids):
            if atom_id not in id_to_full_index:
                raise ValueError(
                    f"{cluster_label}: atom {atom_id} is absent "
                    f"from the source scaffold"
                )

            full_neighbors = full_adjacency.get(atom_id, set())
            internal_neighbors = (
                full_neighbors & cluster_bn_ids
            )
            omitted_neighbors = (
                full_neighbors - cluster_bn_ids
            )

            internal_adjacency[atom_id] = set(
                internal_neighbors
            )

            full_coordination = len(full_neighbors)
            internal_coordination = len(internal_neighbors)
            cut_bond_count = len(omitted_neighbors)

            if cut_bond_count > 0:
                n_undercoordinated_atoms += 1

            if full_coordination != 3:
                n_full_coordination_anomalies += 1

            total_cut_bonds += cut_bond_count

            for neighbor_id in internal_neighbors:
                if atom_id < neighbor_id:
                    internal_lengths.append(
                        bond_lengths[
                            tuple(sorted((atom_id, neighbor_id)))
                        ]
                    )

            for neighbor_id in omitted_neighbors:
                cut_lengths.append(
                    bond_lengths[
                        tuple(sorted((atom_id, neighbor_id)))
                    ]
                )

            atom_type = int(
                full_bn_types[id_to_full_index[atom_id]]
            )

            all_atom_rows.append(
                {
                    "cluster_label": cluster_label,
                    "source_label": inventory_row["source_label"],
                    "pair_index": int(inventory_row["pair_index"]),
                    "atom_id": atom_id,
                    "atom_type": atom_type,
                    "element": "B" if atom_type == 1 else "N",
                    "full_coordination": full_coordination,
                    "internal_coordination": internal_coordination,
                    "n_cut_bonds": cut_bond_count,
                    "internal_neighbor_ids": " ".join(
                        str(value)
                        for value in sorted(internal_neighbors)
                    ),
                    "omitted_neighbor_ids": " ".join(
                        str(value)
                        for value in sorted(omitted_neighbors)
                    ),
                }
            )

        components = connected_components(
            sorted(cluster_bn_ids),
            internal_adjacency,
        )

        n_internal_bonds = sum(
            len(neighbors)
            for neighbors in internal_adjacency.values()
        ) // 2

        cluster_status = "passivation_required"

        if len(components) != 1:
            cluster_status = "disconnected_bn_patch"

        if n_full_coordination_anomalies > 0:
            cluster_status = (
                cluster_status
                + ";source_coordination_anomaly"
            )

        all_cluster_rows.append(
            {
                "cluster_label": cluster_label,
                "source_label": inventory_row["source_label"],
                "pair_index": int(inventory_row["pair_index"]),
                "source_timestep": int(
                    inventory_row["source_timestep"]
                ),
                "n_bn_atoms": len(cluster_bn_ids),
                "n_internal_bn_bonds": n_internal_bonds,
                "n_connected_components": len(components),
                "component_sizes": " ".join(
                    str(len(component))
                    for component in components
                ),
                "n_undercoordinated_bn_atoms": (
                    n_undercoordinated_atoms
                ),
                "n_cut_bn_bonds": total_cut_bonds,
                "n_full_coordination_anomalies": (
                    n_full_coordination_anomalies
                ),
                "minimum_internal_bn_bond_A": (
                    min(internal_lengths)
                    if internal_lengths else np.nan
                ),
                "maximum_internal_bn_bond_A": (
                    max(internal_lengths)
                    if internal_lengths else np.nan
                ),
                "mean_internal_bn_bond_A": (
                    np.mean(internal_lengths)
                    if internal_lengths else np.nan
                ),
                "mean_cut_bn_bond_A": (
                    np.mean(cut_lengths)
                    if cut_lengths else np.nan
                ),
                "proposed_number_of_H_caps": total_cut_bonds,
                "status": cluster_status,
            }
        )

        print(
            f"{cluster_label}: "
            f"BN={len(cluster_bn_ids)}, "
            f"bonds={n_internal_bonds}, "
            f"components={len(components)}, "
            f"cut bonds={total_cut_bonds}, "
            f"coordination anomalies="
            f"{n_full_coordination_anomalies}"
        )

    cluster_summary = pd.DataFrame(all_cluster_rows)
    atom_summary = pd.DataFrame(all_atom_rows)

    cluster_output = (
        args.output_dir
        / "day009_bn_cluster_connectivity_summary.csv"
    )
    atom_output = (
        args.output_dir
        / "day009_bn_cluster_atom_connectivity.csv"
    )

    cluster_summary.to_csv(cluster_output, index=False)
    atom_summary.to_csv(atom_output, index=False)

    print("\nCluster-level connectivity summary:")
    print(
        cluster_summary[
            [
                "source_label",
                "pair_index",
                "n_bn_atoms",
                "n_internal_bn_bonds",
                "n_connected_components",
                "n_undercoordinated_bn_atoms",
                "n_cut_bn_bonds",
                "proposed_number_of_H_caps",
                "status",
            ]
        ].to_string(index=False)
    )

    print(f"\nWrote {cluster_output}")
    print(f"Wrote {atom_output}")


if __name__ == "__main__":
    main()
