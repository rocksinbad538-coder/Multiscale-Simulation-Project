#!/usr/bin/env python3

from pathlib import Path
import json
import shutil

import numpy as np
import pandas as pd

from sklearn.decomposition import PCA


ROOT = Path(__file__).resolve().parents[2]

BASE = (
    ROOT
    / "runs"
    / "phase2"
    / "campaign_phase5_corrected"
)

DESCRIPTOR = (
    BASE
    / "snapshot_descriptors"
    / "snapshot_descriptor_matrix.csv"
)

ASSIGNMENTS = (
    BASE
    / "cluster_characterization"
    / "k3_assignments.csv"
)

SOURCE = (
    BASE
    / "representative_ensemble"
)

OUT = (
    BASE
    / "tddft_core_selection"
)

SELECTED = (
    OUT
    / "selected_dumps"
)

N_PER_CLUSTER = 5


if OUT.exists():
    shutil.rmtree(OUT)

SELECTED.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 1. LOAD FINAL STRUCTURAL MODEL
# ============================================================

desc = pd.read_csv(
    DESCRIPTOR
)

assign = pd.read_csv(
    ASSIGNMENTS
)

feature_cols = [
    c for c in desc.columns
    if c.startswith("f")
]

if len(desc) != 250:
    raise RuntimeError(
        f"Expected 250 descriptors, found {len(desc)}"
    )

if len(assign) != 250:
    raise RuntimeError(
        f"Expected 250 cluster assignments, found {len(assign)}"
    )


X = (
    desc[
        feature_cols
    ]
    .to_numpy(dtype=float)
)

pca = PCA()

PC_full = pca.fit_transform(
    X
)

cumvar = np.cumsum(
    pca.explained_variance_ratio_
)

n95 = int(
    np.searchsorted(
        cumvar,
        0.95
    ) + 1
)

if n95 != 11:
    raise RuntimeError(
        f"Expected PCA95=11 components, found {n95}"
    )

PC = PC_full[:, :n95]


# Preserve descriptor row identity.
meta = desc[
    [
        "temperature_K",
        "snapshot",
    ]
].copy()

for i in range(n95):
    meta[
        f"PC{i+1}"
    ] = PC[:, i]


df = meta.merge(
    assign[
        [
            "temperature_K",
            "snapshot",
            "cluster",
        ]
    ],
    on=[
        "temperature_K",
        "snapshot",
    ],
    how="inner",
    validate="one_to_one",
)

if len(df) != 250:
    raise RuntimeError(
        "Descriptor/cluster merge did not preserve 250 structures."
    )


# ============================================================
# 2. GREEDY DIVERSITY SELECTION
# ============================================================

pc_cols = [
    f"PC{i+1}"
    for i in range(n95)
]

selected_records = []


for cluster in (0, 1, 2):

    sub = (
        df[
            df["cluster"] == cluster
        ]
        .copy()
        .reset_index(drop=True)
    )

    coords = (
        sub[
            pc_cols
        ]
        .to_numpy(dtype=float)
    )

    centroid = (
        coords.mean(axis=0)
    )

    distance_to_centroid = (
        np.linalg.norm(
            coords-centroid,
            axis=1
        )
    )

    # Start with observed structure nearest cluster centroid.
    first = int(
        np.argmin(
            distance_to_centroid
        )
    )

    chosen = [first]

    available_temperatures = sorted(
        sub[
            "temperature_K"
        ]
        .astype(int)
        .unique()
        .tolist()
    )

    while len(chosen) < N_PER_CLUSTER:

        remaining = [
            i
            for i in range(len(sub))
            if i not in chosen
        ]

        represented_T = {
            int(
                sub.iloc[i][
                    "temperature_K"
                ]
            )
            for i in chosen
        }

        missing_T = [
            T
            for T in available_temperatures
            if T not in represented_T
        ]

        slots_left = (
            N_PER_CLUSTER
            - len(chosen)
        )

        candidates = remaining

        # If enough slots remain, force temperature coverage
        # before spending additional selections within an
        # already represented temperature.
        if (
            missing_T
            and slots_left <= len(missing_T)
        ):
            candidates = [
                i
                for i in remaining
                if int(
                    sub.iloc[i][
                        "temperature_K"
                    ]
                ) in missing_T
            ]

        elif missing_T:
            restricted = [
                i
                for i in remaining
                if int(
                    sub.iloc[i][
                        "temperature_K"
                    ]
                ) in missing_T
            ]

            if restricted:
                candidates = restricted

        best = None
        best_score = -np.inf

        for i in candidates:

            distances = [
                np.linalg.norm(
                    coords[i]
                    - coords[j]
                )
                for j in chosen
            ]

            score = min(
                distances
            )

            if score > best_score:
                best_score = score
                best = i

        chosen.append(
            int(best)
        )


    # --------------------------------------------------------
    # Record / copy selected structures
    # --------------------------------------------------------

    for rank, idx in enumerate(
        chosen,
        start=1
    ):

        row = sub.iloc[idx]

        T = int(
            row["temperature_K"]
        )

        snapshot = str(
            row["snapshot"]
        )

        src = (
            SOURCE
            / f"{T}K"
            / f"{snapshot}.dump"
        )

        if not src.exists():
            raise RuntimeError(
                f"Missing source snapshot: {src}"
            )

        dst_name = (
            f"cluster{cluster}_"
            f"rank{rank:02d}_"
            f"{T}K_"
            f"{snapshot}.dump"
        )

        dst = (
            SELECTED
            / dst_name
        )

        shutil.copy2(
            src,
            dst
        )

        selected_records.append({
            "cluster":
                cluster,
            "selection_rank":
                rank,
            "temperature_K":
                T,
            "snapshot":
                snapshot,
            "source":
                str(
                    src.relative_to(ROOT)
                ),
            "selected_file":
                str(
                    dst.relative_to(ROOT)
                ),
            "distance_to_cluster_centroid":
                float(
                    distance_to_centroid[idx]
                ),
            "selection_role":
                (
                    "MEDOID_LIKE"
                    if rank == 1
                    else
                    "DIVERSITY_REPRESENTATIVE"
                ),
        })


selection = pd.DataFrame(
    selected_records
)


# ============================================================
# 3. AUDIT
# ============================================================

cluster_counts = (
    selection[
        "cluster"
    ]
    .value_counts()
    .sort_index()
)

temperature_by_cluster = (
    pd.crosstab(
        selection[
            "cluster"
        ],
        selection[
            "temperature_K"
        ]
    )
)

population = (
    assign[
        "cluster"
    ]
    .value_counts(
        normalize=True
    )
    .sort_index()
)

population_counts = (
    assign[
        "cluster"
    ]
    .value_counts()
    .sort_index()
)


selection.to_csv(
    OUT
    / "PHASE5_TDDFT_CORE_SELECTION.csv",
    index=False
)


manifest = {
    "selection_purpose":
        "STATE_BALANCED_ELECTRONIC_CHARACTERIZATION",
    "not_population_weighted":
        True,
    "final_structural_model":
        {
            "k":
                3,
            "PCA_components":
                n95,
            "PCA_retained_variance":
                float(
                    cumvar[n95-1]
                ),
        },
    "selection_strategy":
        (
            "cluster medoid-like structure plus "
            "temperature-aware farthest-point sampling "
            "in PCA95 structural space"
        ),
    "structures_per_cluster":
        N_PER_CLUSTER,
    "total_selected":
        int(
            len(selection)
        ),
    "original_cluster_population_counts":
        {
            str(int(k)):
                int(v)
            for k,v in
            population_counts.items()
        },
    "original_cluster_population_fractions":
        {
            str(int(k)):
                float(v)
            for k,v in
            population.items()
        },
    "selected_structures":
        selected_records,
}


(
    OUT
    / "PHASE5_TDDFT_CORE_SELECTION.json"
).write_text(
    json.dumps(
        manifest,
        indent=2
    )
    + "\n"
)


checks = [
    len(selection) == 15,
    all(
        cluster_counts.get(c,0) == 5
        for c in (0,1,2)
    ),
    selection[
        "selected_file"
    ].is_unique,
]

all_pass = all(
    checks
)


print("="*90)
print("PHASE5-E34")
print("STATE-BALANCED TDDFT CORE ENSEMBLE")
print("="*90)

print()
print(
    f"PCA_COMPONENTS={n95}"
)

print(
    f"PCA_RETAINED_VARIANCE="
    f"{cumvar[n95-1]:.6f}"
)

print()
print("ORIGINAL CLUSTER POPULATIONS:")

for c in (0,1,2):

    print(
        f"cluster={c}  "
        f"N={population_counts[c]}  "
        f"fraction={population[c]:.4f}"
    )

print()
print("SELECTED COUNTS:")
print(
    cluster_counts.to_string()
)

print()
print(
    "SELECTED TEMPERATURE x CLUSTER:"
)

print(
    temperature_by_cluster.to_string()
)

print()
print("SELECTION:")

print(
    selection[
        [
            "cluster",
            "selection_rank",
            "temperature_K",
            "snapshot",
            "selection_role",
            "distance_to_cluster_centroid",
        ]
    ].to_string(
        index=False
    )
)

print()
print(
    "TDDFT_CORE_SELECTION="
    + (
        "PASS"
        if all_pass
        else "FAIL"
    )
)

print(
    OUT
    / "PHASE5_TDDFT_CORE_SELECTION.json"
)

if not all_pass:
    raise SystemExit(1)
