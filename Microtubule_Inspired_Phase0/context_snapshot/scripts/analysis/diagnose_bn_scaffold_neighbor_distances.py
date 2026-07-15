#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


SOURCE = Path(
    "systems/hybrid/bn_like_chromophore_scaffold_water/"
    "representative_frames/fieldfree_representative_step_14000.lammpstrj"
)

OUTPUT_DIR = Path("results/local_cluster_connectivity")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def read_frame(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        if not handle.readline().startswith("ITEM: TIMESTEP"):
            raise ValueError("Missing TIMESTEP")
        timestep = int(handle.readline().strip())

        if not handle.readline().startswith("ITEM: NUMBER OF ATOMS"):
            raise ValueError("Missing NUMBER OF ATOMS")
        natoms = int(handle.readline().strip())

        box_header = handle.readline().strip()
        if not box_header.startswith("ITEM: BOX BOUNDS"):
            raise ValueError("Missing BOX BOUNDS")

        bounds = []
        for _ in range(3):
            values = handle.readline().split()
            bounds.append((float(values[0]), float(values[1])))

        atom_header = handle.readline().strip()
        columns = atom_header.split()[2:]
        col = {name: columns.index(name) for name in columns}

        required = {"id", "type", "x", "y", "z"}
        missing = required.difference(columns)
        if missing:
            raise ValueError(f"Missing columns: {sorted(missing)}")

        ids = np.empty(natoms, dtype=np.int64)
        types = np.empty(natoms, dtype=np.int32)
        xyz = np.empty((natoms, 3), dtype=np.float64)

        for i in range(natoms):
            fields = handle.readline().split()
            ids[i] = int(fields[col["id"]])
            types[i] = int(fields[col["type"]])
            xyz[i] = [
                float(fields[col["x"]]),
                float(fields[col["y"]]),
                float(fields[col["z"]]),
            ]

    return timestep, ids, types, xyz, bounds


timestep, ids, types, xyz, bounds = read_frame(SOURCE)

b_mask = types == 1
n_mask = types == 2

b_ids = ids[b_mask]
n_ids = ids[n_mask]
b_xyz = xyz[b_mask]
n_xyz = xyz[n_mask]

print(f"Timestep: {timestep}")
print(f"B atoms: {len(b_ids)}")
print(f"N atoms: {len(n_ids)}")

tree_n = cKDTree(n_xyz)
tree_b = cKDTree(b_xyz)

# Several nearest opposite-type neighbors for every atom.
b_dist, b_idx = tree_n.query(b_xyz, k=6)
n_dist, n_idx = tree_b.query(n_xyz, k=6)

rows = []

for local_index, atom_id in enumerate(b_ids):
    for rank in range(b_dist.shape[1]):
        rows.append(
            {
                "central_atom_id": int(atom_id),
                "central_type": 1,
                "central_element": "B",
                "neighbor_rank": rank + 1,
                "neighbor_atom_id": int(n_ids[b_idx[local_index, rank]]),
                "neighbor_type": 2,
                "neighbor_element": "N",
                "distance_A": float(b_dist[local_index, rank]),
            }
        )

for local_index, atom_id in enumerate(n_ids):
    for rank in range(n_dist.shape[1]):
        rows.append(
            {
                "central_atom_id": int(atom_id),
                "central_type": 2,
                "central_element": "N",
                "neighbor_rank": rank + 1,
                "neighbor_atom_id": int(b_ids[n_idx[local_index, rank]]),
                "neighbor_type": 1,
                "neighbor_element": "B",
                "distance_A": float(n_dist[local_index, rank]),
            }
        )

records = pd.DataFrame(rows)

records_path = OUTPUT_DIR / "day009_bn_opposite_type_neighbor_distances.csv"
records.to_csv(records_path, index=False)

summary = (
    records.groupby(
        ["central_element", "neighbor_rank"],
        as_index=False,
    )
    .agg(
        minimum_A=("distance_A", "min"),
        q01_A=("distance_A", lambda x: x.quantile(0.01)),
        q05_A=("distance_A", lambda x: x.quantile(0.05)),
        median_A=("distance_A", "median"),
        q95_A=("distance_A", lambda x: x.quantile(0.95)),
        maximum_A=("distance_A", "max"),
    )
)

summary_path = OUTPUT_DIR / "day009_bn_neighbor_rank_summary.csv"
summary.to_csv(summary_path, index=False)

print("\nOpposite-type neighbor distances by rank:")
print(summary.to_string(index=False))

print("\nNearest-neighbor global quantiles:")
nearest = records[records["neighbor_rank"] == 1]["distance_A"]
print(nearest.describe(percentiles=[0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]))

# Coordination scan over candidate cutoffs.
scan_rows = []

for cutoff in np.arange(1.5, 5.01, 0.05):
    b_coord = (b_dist <= cutoff).sum(axis=1)
    n_coord = (n_dist <= cutoff).sum(axis=1)
    coordination = np.concatenate([b_coord, n_coord])

    scan_rows.append(
        {
            "cutoff_A": cutoff,
            "mean_coordination": coordination.mean(),
            "minimum_coordination": coordination.min(),
            "maximum_coordination": coordination.max(),
            "fraction_coordination_0": np.mean(coordination == 0),
            "fraction_coordination_1": np.mean(coordination == 1),
            "fraction_coordination_2": np.mean(coordination == 2),
            "fraction_coordination_3": np.mean(coordination == 3),
            "fraction_coordination_gt3": np.mean(coordination > 3),
        }
    )

scan = pd.DataFrame(scan_rows)
scan_path = OUTPUT_DIR / "day009_bn_coordination_cutoff_scan.csv"
scan.to_csv(scan_path, index=False)

print("\nBest candidate cutoffs by fraction with coordination 3:")
print(
    scan.sort_values(
        ["fraction_coordination_3", "fraction_coordination_gt3"],
        ascending=[False, True],
    )
    .head(15)
    .to_string(index=False)
)

print(f"\nWrote {records_path}")
print(f"Wrote {summary_path}")
print(f"Wrote {scan_path}")
