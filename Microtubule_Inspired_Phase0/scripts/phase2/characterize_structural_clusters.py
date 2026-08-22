#!/usr/bin/env python3

from pathlib import Path
import json
import numpy as np
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import (
    silhouette_samples,
    silhouette_score,
)

ROOT = Path(__file__).resolve().parents[2]

BASE = (
    ROOT
    / "runs"
    / "phase2"
    / "campaign_phase5_corrected"
)

DATA = (
    BASE
    / "snapshot_descriptors"
    / "snapshot_descriptor_matrix.csv"
)

OUT = (
    BASE
    / "cluster_characterization"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)

df = pd.read_csv(DATA)

feature_cols = [
    c for c in df.columns
    if c.startswith("f")
]

X = df[
    feature_cols
].to_numpy(dtype=float)

# ------------------------------------------------------------
# PCA retaining >=95% variance
# ------------------------------------------------------------

pca = PCA()
PC_full = pca.fit_transform(X)

cumvar = np.cumsum(
    pca.explained_variance_ratio_
)

n95 = int(
    np.searchsorted(
        cumvar,
        0.95
    ) + 1
)

PC = PC_full[:, :n95]

print("="*90)
print("PHASE5-E23")
print("STRUCTURAL CLUSTER CHARACTERIZATION")
print("="*90)

print(f"PCA_COMPONENTS={n95}")
print(
    f"PCA_RETAINED_VARIANCE="
    f"{cumvar[n95-1]:.6f}"
)

summary = {}

for k in (2, 3, 4):

    model = KMeans(
        n_clusters=k,
        random_state=42,
        n_init=100,
    )

    labels = model.fit_predict(PC)

    sil_samples = silhouette_samples(
        PC,
        labels
    )

    sil_global = silhouette_score(
        PC,
        labels
    )

    local = df[
        [
            "temperature_K",
            "snapshot",
        ]
    ].copy()

    local["cluster"] = labels
    local["silhouette"] = sil_samples

    # PC coordinates useful for diagnostics.
    local["PC1"] = PC_full[:,0]
    local["PC2"] = PC_full[:,1]

    local.to_csv(
        OUT / f"k{k}_assignments.csv",
        index=False
    )

    counts = (
        local["cluster"]
        .value_counts()
        .sort_index()
    )

    contingency = pd.crosstab(
        local["temperature_K"],
        local["cluster"],
    )

    # Normalize within temperature.
    fractions = contingency.div(
        contingency.sum(axis=1),
        axis=0,
    )

    # Distance of each point to assigned centroid.
    distances = np.linalg.norm(
        PC
        - model.cluster_centers_[labels],
        axis=1,
    )

    # Distance between cluster centroids.
    centroid_distances = {}

    for i in range(k):
        for j in range(i+1, k):

            d = float(
                np.linalg.norm(
                    model.cluster_centers_[i]
                    - model.cluster_centers_[j]
                )
            )

            centroid_distances[
                f"{i}-{j}"
            ] = d

    # Medoid-like representative:
    # observed snapshot closest to centroid.
    representatives = []

    for c in range(k):

        idx = np.where(
            labels == c
        )[0]

        subdist = distances[idx]

        best = idx[
            np.argmin(subdist)
        ]

        representatives.append({
            "cluster":
                int(c),
            "temperature_K":
                int(
                    df.iloc[best][
                        "temperature_K"
                    ]
                ),
            "snapshot":
                str(
                    df.iloc[best][
                        "snapshot"
                    ]
                ),
            "distance_to_centroid":
                float(
                    distances[best]
                ),
            "PC1":
                float(
                    PC_full[best,0]
                ),
            "PC2":
                float(
                    PC_full[best,1]
                ),
        })

    # Cluster-wise silhouette statistics.
    cluster_silhouette = {}

    for c in range(k):

        vals = sil_samples[
            labels == c
        ]

        cluster_silhouette[
            str(c)
        ] = {
            "n":
                int(len(vals)),
            "mean":
                float(np.mean(vals)),
            "median":
                float(np.median(vals)),
            "min":
                float(np.min(vals)),
            "fraction_negative":
                float(
                    np.mean(vals < 0)
                ),
        }

    result = {
        "k":
            k,
        "global_silhouette":
            float(sil_global),
        "cluster_counts":
            {
                str(int(i)): int(v)
                for i, v
                in counts.items()
            },
        "cluster_silhouette":
            cluster_silhouette,
        "centroid_distances":
            centroid_distances,
        "representatives":
            representatives,
    }

    summary[str(k)] = result

    print()
    print("="*72)
    print(f"k = {k}")
    print("="*72)

    print(
        f"GLOBAL_SILHOUETTE="
        f"{sil_global:.6f}"
    )

    print()
    print("CLUSTER_COUNTS:")
    print(counts.to_string())

    print()
    print("TEMPERATURE x CLUSTER COUNTS:")
    print(contingency.to_string())

    print()
    print("TEMPERATURE x CLUSTER FRACTIONS:")
    print(
        fractions
        .round(3)
        .to_string()
    )

    print()
    print("CLUSTER SILHOUETTE:")

    for c, rec in cluster_silhouette.items():

        print(
            f"cluster={c} "
            f"n={rec['n']} "
            f"mean={rec['mean']:.4f} "
            f"median={rec['median']:.4f} "
            f"negative_fraction="
            f"{rec['fraction_negative']:.4f}"
        )

    print()
    print("REPRESENTATIVES:")

    for r in representatives:

        print(
            f"cluster={r['cluster']} "
            f"T={r['temperature_K']}K "
            f"snapshot={r['snapshot']} "
            f"distance={r['distance_to_centroid']:.6f}"
        )


(
    OUT
    / "cluster_characterization.json"
).write_text(
    json.dumps(
        {
            "PCA_components":
                n95,
            "PCA_retained_variance":
                float(
                    cumvar[n95-1]
                ),
            "models":
                summary,
        },
        indent=2,
    )
    + "\n"
)

print()
print("="*90)
print(
    OUT
    / "cluster_characterization.json"
)
print("PHASE5-E23 COMPLETE")
print("="*90)
