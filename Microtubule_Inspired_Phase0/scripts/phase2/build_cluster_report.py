#!/usr/bin/env python3

from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

CLUSTERS = ROOT/"runs"/"phase2"/"campaign"/"snapshot_clustering"/"snapshot_clusters.csv"
REP = ROOT/"runs"/"phase2"/"campaign"/"snapshot_clustering"/"cluster_representatives.csv"

OUT = ROOT/"runs"/"phase2"/"campaign"/"snapshot_clustering"

df = pd.read_csv(CLUSTERS)
rep = pd.read_csv(REP)

rows = []

for cluster in sorted(df.cluster.unique()):

    sub = df[df.cluster == cluster]

    rows.append({

        "cluster": int(cluster),

        "population": int(len(sub)),

        "fraction": float(len(sub)/len(df)),

        "temperatures":

            sorted(sub.temperature_K.unique().tolist()),

        "representative_snapshot":

            rep.loc[rep.cluster==cluster,"snapshot"].iloc[0]

    })

report = pd.DataFrame(rows)

report.to_csv(
    OUT/"cluster_population_report.csv",
    index=False
)

with open(
    OUT/"cluster_population_report.json",
    "w"
) as f:

    json.dump(
        rows,
        f,
        indent=2
    )

print("="*90)
print("PHASE3-A03")
print("CLUSTER POPULATION REPORT")
print("="*90)

print(report)

print()

print(OUT)
