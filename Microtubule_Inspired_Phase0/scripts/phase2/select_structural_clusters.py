#!/usr/bin/env python3

from pathlib import Path
import json

import numpy as np
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import (
    silhouette_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    adjusted_rand_score,
)

ROOT = Path(__file__).resolve().parents[2]

DATA = (
    ROOT
    / "runs"
    / "phase2"
    / "campaign_phase5_corrected"
    / "snapshot_descriptors"
    / "snapshot_descriptor_matrix.csv"
)

OUT = (
    ROOT
    / "runs"
    / "phase2"
    / "campaign_phase5_corrected"
    / "cluster_selection"
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
# PCA
# ------------------------------------------------------------

pca_full = PCA()

PC_full = pca_full.fit_transform(
    X
)

cumvar = np.cumsum(
    pca_full.explained_variance_ratio_
)

n95 = int(
    np.searchsorted(
        cumvar,
        0.95
    ) + 1
)

n90 = int(
    np.searchsorted(
        cumvar,
        0.90
    ) + 1
)

n99 = int(
    np.searchsorted(
        cumvar,
        0.99
    ) + 1
)

# Cluster using the 95%-variance PCA subspace.
PC = PC_full[:, :n95]


# ------------------------------------------------------------
# Candidate k evaluation
# ------------------------------------------------------------

K_MIN = 2
K_MAX = 15

SEEDS = [
    11,
    23,
    42,
    71,
    101,
]

rows = []

for k in range(
    K_MIN,
    K_MAX + 1
):

    reference_labels = None

    silhouettes = []
    calinskis = []
    davies = []
    ari_to_reference = []
    inertias = []

    for seed in SEEDS:

        model = KMeans(
            n_clusters=k,
            random_state=seed,
            n_init=50,
        )

        labels = model.fit_predict(
            PC
        )

        silhouettes.append(
            silhouette_score(
                PC,
                labels
            )
        )

        calinskis.append(
            calinski_harabasz_score(
                PC,
                labels
            )
        )

        davies.append(
            davies_bouldin_score(
                PC,
                labels
            )
        )

        inertias.append(
            model.inertia_
        )

        if reference_labels is None:

            reference_labels = (
                labels.copy()
            )

        else:

            ari_to_reference.append(
                adjusted_rand_score(
                    reference_labels,
                    labels
                )
            )

    rows.append({
        "k":
            k,
        "silhouette_mean":
            float(
                np.mean(silhouettes)
            ),
        "silhouette_std":
            float(
                np.std(silhouettes)
            ),
        "calinski_harabasz_mean":
            float(
                np.mean(calinskis)
            ),
        "davies_bouldin_mean":
            float(
                np.mean(davies)
            ),
        "kmeans_stability_ARI_mean":
            float(
                np.mean(
                    ari_to_reference
                )
            )
            if ari_to_reference
            else 1.0,
        "inertia_mean":
            float(
                np.mean(inertias)
            ),
    })


metrics = pd.DataFrame(
    rows
)

metrics.to_csv(
    OUT
    / "cluster_selection_metrics.csv",
    index=False
)


# ------------------------------------------------------------
# Rankings
# ------------------------------------------------------------

metrics[
    "rank_silhouette"
] = metrics[
    "silhouette_mean"
].rank(
    ascending=False,
    method="min",
)

metrics[
    "rank_calinski"
] = metrics[
    "calinski_harabasz_mean"
].rank(
    ascending=False,
    method="min",
)

metrics[
    "rank_davies"
] = metrics[
    "davies_bouldin_mean"
].rank(
    ascending=True,
    method="min",
)

metrics[
    "rank_stability"
] = metrics[
    "kmeans_stability_ARI_mean"
].rank(
    ascending=False,
    method="min",
)

metrics[
    "consensus_rank"
] = (
    metrics["rank_silhouette"]
    +
    metrics["rank_calinski"]
    +
    metrics["rank_davies"]
    +
    metrics["rank_stability"]
)

metrics = metrics.sort_values(
    [
        "consensus_rank",
        "rank_silhouette",
    ]
)

recommended_k = int(
    metrics.iloc[0]["k"]
)


summary = {
    "n_snapshots":
        int(len(df)),
    "n_cartesian_features":
        int(len(feature_cols)),
    "pca_components_90pct":
        n90,
    "pca_components_95pct":
        n95,
    "pca_components_99pct":
        n99,
    "PC1_variance_fraction":
        float(
            pca_full.explained_variance_ratio_[0]
        ),
    "PC2_variance_fraction":
        float(
            pca_full.explained_variance_ratio_[1]
        ),
    "PC1_PC2_variance_fraction":
        float(
            pca_full.explained_variance_ratio_[:2].sum()
        ),
    "recommended_k_by_consensus":
        recommended_k,
    "candidate_range":
        [K_MIN, K_MAX],
    "seeds":
        SEEDS,
    "cluster_space":
        f"PCA_{n95}_components_95pct_variance",
}

(
    OUT
    / "cluster_selection_summary.json"
).write_text(
    json.dumps(
        summary,
        indent=2,
    )
    + "\n"
)


print("="*90)
print("PHASE5-E22")
print("PCA + CLUSTER MODEL SELECTION")
print("="*90)

print(
    f"SNAPSHOTS={len(df)}"
)

print(
    f"CARTESIAN_FEATURES="
    f"{len(feature_cols)}"
)

print()
print(
    f"PCA_90_COMPONENTS={n90}"
)

print(
    f"PCA_95_COMPONENTS={n95}"
)

print(
    f"PCA_99_COMPONENTS={n99}"
)

print()

print(
    f"PC1_VARIANCE="
    f"{pca_full.explained_variance_ratio_[0]:.6f}"
)

print(
    f"PC2_VARIANCE="
    f"{pca_full.explained_variance_ratio_[1]:.6f}"
)

print(
    f"PC1_PC2_VARIANCE="
    f"{pca_full.explained_variance_ratio_[:2].sum():.6f}"
)

print()
print("CLUSTER METRICS:")
print(
    metrics[
        [
            "k",
            "silhouette_mean",
            "calinski_harabasz_mean",
            "davies_bouldin_mean",
            "kmeans_stability_ARI_mean",
            "consensus_rank",
        ]
    ].to_string(
        index=False
    )
)

print()
print(
    f"RECOMMENDED_K="
    f"{recommended_k}"
)

print(
    f"CLUSTER_SPACE="
    f"PCA_{n95}_COMPONENTS"
)

print()
print(
    OUT
    / "cluster_selection_summary.json"
)
