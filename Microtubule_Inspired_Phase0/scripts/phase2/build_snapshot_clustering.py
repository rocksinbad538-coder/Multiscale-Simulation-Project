#!/usr/bin/env python3

from pathlib import Path
import json

import matplotlib.pyplot as plt
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

ROOT = Path(__file__).resolve().parents[2]

DATA = ROOT/"runs"/"phase2"/"campaign_phase5_corrected"/"snapshot_descriptors"/"snapshot_descriptor_matrix.csv"

OUT = ROOT/"runs"/"phase2"/"campaign_phase5_corrected"/"snapshot_clustering"

OUT.mkdir(exist_ok=True)

df = pd.read_csv(DATA)

X = df.iloc[:,2:].values

pca = PCA(n_components=2)

PC = pca.fit_transform(X)

kmeans = KMeans(
    n_clusters=10,
    random_state=42,
    n_init=20
)

labels = kmeans.fit_predict(PC)

df["PC1"] = PC[:,0]
df["PC2"] = PC[:,1]
df["cluster"] = labels

df.to_csv(
    OUT/"snapshot_clusters.csv",
    index=False
)

centers = kmeans.cluster_centers_

selected = []

for i in range(10):

    sub = df[df.cluster==i].copy()

    dx = sub.PC1-centers[i,0]
    dy = sub.PC2-centers[i,1]

    d = dx*dx + dy*dy

    idx = d.idxmin()

    selected.append(df.loc[idx])

selected = pd.DataFrame(selected)

selected.to_csv(
    OUT/"cluster_representatives.csv",
    index=False
)

plt.figure(figsize=(8,6))

plt.scatter(
    df.PC1,
    df.PC2,
    c=df.cluster,
    s=18
)

plt.scatter(
    centers[:,0],
    centers[:,1],
    marker="X",
    s=250,
    color="black"
)

plt.xlabel("PC1")
plt.ylabel("PC2")

plt.title("Snapshot clustering")

plt.tight_layout()

plt.savefig(
    OUT/"SnapshotClusters.png",
    dpi=300
)

plt.close()

summary = {
    "n_snapshots": int(len(df)),
    "n_clusters": 10,
    "explained_variance": pca.explained_variance_ratio_.tolist()
}

with open(
    OUT/"cluster_summary.json",
    "w"
) as f:

    json.dump(
        summary,
        f,
        indent=2
    )

print("="*90)
print("PHASE3-A02")
print("PCA + KMEANS")
print("="*90)

print(summary)

print()

print(OUT)
