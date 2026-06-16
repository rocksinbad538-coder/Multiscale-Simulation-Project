#!/usr/bin/env python3

from __future__ import annotations

from collections import deque
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


SOURCE = Path(
    "systems/hybrid/bn_like_chromophore_scaffold_water/"
    "representative_frames/fieldfree_representative_step_14000.lammpstrj"
)

OUTPUT_DIR = Path("results/scaffold_topology")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CUTOFF_A = 3.50


def read_frame(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        if not handle.readline().startswith("ITEM: TIMESTEP"):
            raise ValueError("Missing TIMESTEP")
        timestep = int(handle.readline().strip())

        if not handle.readline().startswith("ITEM: NUMBER OF ATOMS"):
            raise ValueError("Missing NUMBER OF ATOMS")
        natoms = int(handle.readline().strip())

        box_header = handle.readline()
        if not box_header.startswith("ITEM: BOX BOUNDS"):
            raise ValueError("Missing BOX BOUNDS")

        for _ in range(3):
            handle.readline()

        header = handle.readline().split()
        columns = header[2:]
        col = {name: columns.index(name) for name in columns}

        required = {"id", "type", "x", "y", "z"}
        missing = required.difference(columns)

        if missing:
            raise ValueError(f"Missing columns: {sorted(missing)}")

        ids = np.empty(natoms, dtype=np.int64)
        types = np.empty(natoms, dtype=np.int32)
        xyz = np.empty((natoms, 3), dtype=np.float64)

        for index in range(natoms):
            fields = handle.readline().split()

            ids[index] = int(fields[col["id"]])
            types[index] = int(fields[col["type"]])
            xyz[index] = [
                float(fields[col["x"]]),
                float(fields[col["y"]]),
                float(fields[col["z"]]),
            ]

    return timestep, ids, types, xyz


def connected_components(
    adjacency: list[set[int]],
) -> list[list[int]]:
    unvisited = set(range(len(adjacency)))
    components: list[list[int]] = []

    while unvisited:
        start = min(unvisited)
        queue = deque([start])
        component: list[int] = []

        while queue:
            current = queue.popleft()

            if current not in unvisited:
                continue

            unvisited.remove(current)
            component.append(current)

            for neighbor in adjacency[current]:
                if neighbor in unvisited:
                    queue.append(neighbor)

        components.append(sorted(component))

    components.sort(key=lambda values: (-len(values), values[0]))
    return components


def main() -> None:
    timestep, ids, types, xyz = read_frame(SOURCE)

    mask = np.isin(types, [1, 2])
    ids = ids[mask]
    types = types[mask]
    xyz = xyz[mask]

    tree = cKDTree(xyz)
    candidate_pairs = tree.query_pairs(CUTOFF_A)

    adjacency = [set() for _ in range(len(ids))]
    edge_rows = []

    for index_a, index_b in candidate_pairs:
        if types[index_a] == types[index_b]:
            continue

        distance = float(
            np.linalg.norm(xyz[index_a] - xyz[index_b])
        )

        adjacency[index_a].add(index_b)
        adjacency[index_b].add(index_a)

        edge_rows.append(
            {
                "atom_id_i": int(ids[index_a]),
                "atom_id_j": int(ids[index_b]),
                "type_i": int(types[index_a]),
                "type_j": int(types[index_b]),
                "distance_A": distance,
            }
        )

    components = connected_components(adjacency)

    component_rows = []
    atom_rows = []

    for component_index, component in enumerate(
        components,
        start=1,
    ):
        component_xyz = xyz[component]
        component_ids = ids[component]
        component_types = types[component]

        center = component_xyz.mean(axis=0)
        center_r = float(np.hypot(center[0], center[1]))
        center_phi = float(np.arctan2(center[1], center[0]))

        internal_edges = (
            sum(len(adjacency[index]) for index in component) // 2
        )

        degrees = np.asarray(
            [len(adjacency[index]) for index in component],
            dtype=int,
        )

        component_rows.append(
            {
                "component_index": component_index,
                "n_atoms": len(component),
                "n_edges": internal_edges,
                "is_tree": internal_edges == len(component) - 1,
                "n_B": int((component_types == 1).sum()),
                "n_N": int((component_types == 2).sum()),
                "n_degree_1": int((degrees == 1).sum()),
                "n_degree_2": int((degrees == 2).sum()),
                "n_degree_other": int(
                    ((degrees != 1) & (degrees != 2)).sum()
                ),
                "minimum_atom_id": int(component_ids.min()),
                "maximum_atom_id": int(component_ids.max()),
                "center_x_A": float(center[0]),
                "center_y_A": float(center[1]),
                "center_z_A": float(center[2]),
                "center_r_A": center_r,
                "center_phi_rad": center_phi,
                "minimum_z_A": float(component_xyz[:, 2].min()),
                "maximum_z_A": float(component_xyz[:, 2].max()),
                "z_span_A": float(
                    component_xyz[:, 2].max()
                    - component_xyz[:, 2].min()
                ),
            }
        )

        for local_index in component:
            atom_rows.append(
                {
                    "component_index": component_index,
                    "atom_id": int(ids[local_index]),
                    "atom_type": int(types[local_index]),
                    "element": (
                        "B" if types[local_index] == 1 else "N"
                    ),
                    "degree": len(adjacency[local_index]),
                    "x_A": xyz[local_index, 0],
                    "y_A": xyz[local_index, 1],
                    "z_A": xyz[local_index, 2],
                }
            )

    component_df = pd.DataFrame(component_rows)
    atom_df = pd.DataFrame(atom_rows)
    edge_df = pd.DataFrame(edge_rows)

    component_path = (
        OUTPUT_DIR / "day009_scaffold_component_summary.csv"
    )
    atom_path = (
        OUTPUT_DIR / "day009_scaffold_component_atoms.csv"
    )
    edge_path = (
        OUTPUT_DIR / "day009_scaffold_edges.csv"
    )

    component_df.to_csv(component_path, index=False)
    atom_df.to_csv(atom_path, index=False)
    edge_df.to_csv(edge_path, index=False)

    print(f"Timestep: {timestep}")
    print(f"BN-like sites: {len(ids)}")
    print(f"Edges: {len(edge_df)}")
    print(f"Connected components: {len(component_df)}")

    print("\nComponent-size distribution:")
    print(
        component_df.groupby(
            [
                "n_atoms",
                "n_edges",
                "is_tree",
                "n_degree_1",
                "n_degree_2",
                "n_degree_other",
            ]
        )
        .size()
        .reset_index(name="n_components")
        .to_string(index=False)
    )

    print("\nCoordinate ranges over components:")
    print(
        component_df[
            [
                "n_atoms",
                "center_r_A",
                "minimum_z_A",
                "maximum_z_A",
                "z_span_A",
            ]
        ].describe().to_string()
    )

    print("\nFirst 12 components:")
    print(component_df.head(12).to_string(index=False))

    print(f"\nWrote {component_path}")
    print(f"Wrote {atom_path}")
    print(f"Wrote {edge_path}")


if __name__ == "__main__":
    main()
