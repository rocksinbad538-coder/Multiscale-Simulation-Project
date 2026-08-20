#!/usr/bin/env python3

from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

BASE = (
    ROOT
    / "runs"
    / "phase2"
    / "campaign_phase5_corrected"
)

SEL = (
    BASE
    / "tddft_selected_structures"
    / "TDDFT_SELECTION.json"
)

CLUSTERS = (
    BASE
    / "snapshot_clustering"
    / "snapshot_clusters.csv"
)

REP = (
    BASE
    / "snapshot_clustering"
    / "cluster_representatives.csv"
)

selection = json.loads(
    SEL.read_text()
)

cluster_df = pd.read_csv(
    CLUSTERS
)

rep_df = pd.read_csv(
    REP
)

rows = []

for item in selection:

    T = int(
        item["temperature_K"]
    )

    snap = str(
        item["snapshot"]
    )

    expected_cluster = int(
        item["cluster"]
    )

    cluster_match = cluster_df[
        (cluster_df["temperature_K"] == T)
        &
        (cluster_df["snapshot"].astype(str) == snap)
    ]

    representative_match = rep_df[
        (rep_df["temperature_K"] == T)
        &
        (rep_df["snapshot"].astype(str) == snap)
        &
        (rep_df["cluster"] == expected_cluster)
    ]

    valid_cluster_record = (
        len(cluster_match) == 1
    )

    cluster_consistent = (
        valid_cluster_record
        and
        int(cluster_match["cluster"].iloc[0])
        == expected_cluster
    )

    unique_representative = (
        len(representative_match) == 1
    )

    provenance_pass = (
        valid_cluster_record
        and
        cluster_consistent
        and
        unique_representative
    )

    rows.append({
        "temperature_K": T,
        "snapshot": snap,
        "expected_cluster":
            expected_cluster,
        "cluster_record_unique":
            valid_cluster_record,
        "cluster_consistent":
            cluster_consistent,
        "representative_unique":
            unique_representative,
        "selected_for_TDDFT":
            True,
        "provenance_pass":
            provenance_pass,
    })


audit = pd.DataFrame(rows)

OUT = BASE / "audit"
OUT.mkdir(
    parents=True,
    exist_ok=True
)

csv = (
    OUT
    / "pipeline_provenance.csv"
)

audit.to_csv(
    csv,
    index=False
)

all_pass = (
    len(audit) > 0
    and audit["provenance_pass"].all()
)

print("="*90)
print("PHASE5-D52")
print("PIPELINE PROVENANCE AUDIT")
print("="*90)

print(audit)

print()
print(
    "PIPELINE_PROVENANCE="
    f"{'PASS' if all_pass else 'FAIL'}"
)

print(csv)

if not all_pass:
    raise SystemExit(1)
