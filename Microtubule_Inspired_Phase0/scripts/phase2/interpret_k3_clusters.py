#!/usr/bin/env python3

from pathlib import Path
import json

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

ENSEMBLE = (
    BASE
    / "representative_ensemble"
    / "ensemble_manifest.json"
)

SELECTION = (
    BASE
    / "cluster_selection"
    / "cluster_selection_summary.json"
)

CHARACTERIZATION = (
    BASE
    / "cluster_characterization"
    / "cluster_characterization.json"
)

OUT = (
    BASE
    / "cluster_interpretation_k3"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 1. INPUT CONTRACTS
# ============================================================

required = [
    DESCRIPTOR,
    ASSIGNMENTS,
    ENSEMBLE,
    SELECTION,
    CHARACTERIZATION,
]

missing = [
    str(p)
    for p in required
    if not p.exists()
]

if missing:
    raise RuntimeError(
        "Missing required files:\n"
        + "\n".join(missing)
    )


selection = json.loads(
    SELECTION.read_text()
)

characterization = json.loads(
    CHARACTERIZATION.read_text()
)

recommended_k = int(
    selection[
        "recommended_k_by_consensus"
    ]
)

pca_components = int(
    characterization[
        "PCA_components"
    ]
)

pca_retained = float(
    characterization[
        "PCA_retained_variance"
    ]
)

if recommended_k != 3:
    raise RuntimeError(
        f"Expected final k=3, found k={recommended_k}"
    )

if "3" not in characterization["models"]:
    raise RuntimeError(
        "k=3 model missing from characterization JSON."
    )


# ============================================================
# 2. LOAD FINAL ASSIGNMENTS
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

if len(feature_cols) != 111:
    raise RuntimeError(
        f"Expected 111 Cartesian features, "
        f"found {len(feature_cols)}"
    )

needed_assignment_cols = {
    "temperature_K",
    "snapshot",
    "cluster",
    "PC1",
    "PC2",
}

if not needed_assignment_cols.issubset(
    assign.columns
):
    raise RuntimeError(
        "k3 assignments missing required columns."
    )

if len(assign) != 250:
    raise RuntimeError(
        f"Expected 250 assignments, found {len(assign)}"
    )

clusters = sorted(
    assign["cluster"]
    .astype(int)
    .unique()
    .tolist()
)

if clusters != [0, 1, 2]:
    raise RuntimeError(
        f"Unexpected cluster labels: {clusters}"
    )


# ============================================================
# 3. MANIFEST -> PHYSICAL TIME
# ============================================================

manifest = json.loads(
    ENSEMBLE.read_text()
)

manifest_rows = []

for record in manifest:

    T = int(
        record["temperature_K"]
    )

    for snap in record["snapshots"]:

        manifest_rows.append({
            "temperature_K":
                T,
            "snapshot":
                Path(
                    snap["file"]
                ).stem,
            "step":
                int(
                    snap["step"]
                ),
        })


manifest_df = pd.DataFrame(
    manifest_rows
)

df = assign.merge(
    manifest_df,
    on=[
        "temperature_K",
        "snapshot",
    ],
    how="left",
    validate="one_to_one",
)

if df["step"].isna().any():
    raise RuntimeError(
        "Could not map all assignments to trajectory timesteps."
    )

df["step"] = (
    df["step"]
    .astype(int)
)

df["time_ns"] = (
    df["step"]
    * 0.25
    / 1.0e6
)


# ============================================================
# 4. JOIN PHYSICAL OBSERVABLES
# ============================================================

joined = []

for T in (
    150,
    200,
    250,
    300,
    350,
):

    sub = df[
        df["temperature_K"] == T
    ].copy()

    analysis = (
        BASE
        / f"{T}K"
        / "analysis"
    )

    geom = pd.read_csv(
        analysis
        / "trajectory_summary.csv"
    ).rename(
        columns={
            "timestep":
                "step"
        }
    )

    rmsd = pd.read_csv(
        analysis
        / "aligned_rmsd.csv"
    ).rename(
        columns={
            "timestep":
                "step"
        }
    )

    sub = sub.merge(
        geom[
            [
                "step",
                "Rg",
                "Lx",
                "Ly",
                "Lz",
            ]
        ],
        on="step",
        how="left",
        validate="many_to_one",
    )

    sub = sub.merge(
        rmsd[
            [
                "step",
                "AlignedRMSD",
            ]
        ],
        on="step",
        how="left",
        validate="many_to_one",
    )

    joined.append(
        sub
    )


df = pd.concat(
    joined,
    ignore_index=True
)

physical_cols = [
    "Rg",
    "AlignedRMSD",
    "Lx",
    "Ly",
    "Lz",
    "PC1",
    "PC2",
]

if df[
    physical_cols
].isna().any().any():
    raise RuntimeError(
        "Missing structural observables after join."
    )


# ============================================================
# 5. CLUSTER PHYSICAL SUMMARY
# ============================================================

summary_rows = []

for c in clusters:

    sub = df[
        df["cluster"] == c
    ]

    row = {
        "cluster":
            int(c),
        "n":
            int(len(sub)),
    }

    for metric in physical_cols:

        row[
            f"{metric}_mean"
        ] = float(
            sub[metric].mean()
        )

        row[
            f"{metric}_std"
        ] = float(
            sub[metric].std(ddof=1)
        )

        row[
            f"{metric}_median"
        ] = float(
            sub[metric].median()
        )

    summary_rows.append(
        row
    )


summary = pd.DataFrame(
    summary_rows
)

summary.to_csv(
    OUT
    / "k3_cluster_physical_summary.csv",
    index=False
)


# ============================================================
# 6. TEMPERATURE COMPOSITION
# ============================================================

contingency = pd.crosstab(
    df["temperature_K"],
    df["cluster"],
)

fractions = contingency.div(
    contingency.sum(axis=1),
    axis=0,
)

contingency.to_csv(
    OUT
    / "k3_temperature_cluster_counts.csv"
)

fractions.to_csv(
    OUT
    / "k3_temperature_cluster_fractions.csv"
)


# ============================================================
# 7. PAIRWISE STANDARDIZED DIFFERENCES
# ============================================================

effect_rows = []

for metric in physical_cols:

    for i in range(
        len(clusters)
    ):

        for j in range(
            i+1,
            len(clusters)
        ):

            c1 = clusters[i]
            c2 = clusters[j]

            a = (
                df[
                    df["cluster"] == c1
                ][metric]
                .to_numpy(dtype=float)
            )

            b = (
                df[
                    df["cluster"] == c2
                ][metric]
                .to_numpy(dtype=float)
            )

            pooled = np.sqrt(
                (
                    (len(a)-1)
                    * np.var(
                        a,
                        ddof=1
                    )
                    +
                    (len(b)-1)
                    * np.var(
                        b,
                        ddof=1
                    )
                )
                /
                (
                    len(a)
                    + len(b)
                    - 2
                )
            )

            if pooled > 0.0:

                d = (
                    np.mean(b)
                    - np.mean(a)
                ) / pooled

            else:

                d = np.nan

            effect_rows.append({
                "metric":
                    metric,
                "cluster_A":
                    int(c1),
                "cluster_B":
                    int(c2),
                "mean_A":
                    float(
                        np.mean(a)
                    ),
                "mean_B":
                    float(
                        np.mean(b)
                    ),
                "standardized_difference_B_minus_A":
                    float(d),
            })


effects = pd.DataFrame(
    effect_rows
)

effects.to_csv(
    OUT
    / "k3_pairwise_standardized_differences.csv",
    index=False
)


# ============================================================
# 8. TEMPORAL MEMBERSHIP / INTERCONVERSIONS
# ============================================================

timeline_summary = []

for T in (
    150,
    200,
    250,
    300,
    350,
):

    sub = (
        df[
            df["temperature_K"] == T
        ]
        .sort_values("step")
        .copy()
    )

    labels = (
        sub["cluster"]
        .astype(int)
        .to_numpy()
    )

    transitions = int(
        np.sum(
            labels[1:]
            != labels[:-1]
        )
    )

    cluster_counts = {
        str(c):
            int(
                np.sum(
                    labels == c
                )
            )
        for c in clusters
    }

    timeline_summary.append({
        "temperature_K":
            int(T),
        "n_snapshots":
            int(len(sub)),
        "cluster_transitions":
            transitions,
        "cluster_counts":
            cluster_counts,
        "first_time_ns":
            float(
                sub["time_ns"].min()
            ),
        "last_time_ns":
            float(
                sub["time_ns"].max()
            ),
    })

    sub[
        [
            "snapshot",
            "step",
            "time_ns",
            "cluster",
            "Rg",
            "AlignedRMSD",
            "Lx",
            "Ly",
            "Lz",
            "PC1",
            "PC2",
        ]
    ].to_csv(
        OUT
        / f"{T}K_k3_timeline.csv",
        index=False
    )


# ============================================================
# 9. PCA LOADINGS
# ============================================================

X = (
    desc[
        feature_cols
    ]
    .to_numpy(dtype=float)
)

pca = PCA()

pca.fit(
    X
)

cumvar = np.cumsum(
    pca.explained_variance_ratio_
)

n95_recomputed = int(
    np.searchsorted(
        cumvar,
        0.95
    ) + 1
)

if n95_recomputed != pca_components:
    raise RuntimeError(
        "PCA dimensionality mismatch: "
        f"characterization={pca_components}, "
        f"recomputed={n95_recomputed}"
    )


# Features correspond to aligned Cartesian coordinates:
# [atom1_x, atom1_y, atom1_z, atom2_x, ...].
# We report descriptor atom position 1..37 rather than
# asserting a chemical identity that is not stored in the
# descriptor CSV itself.

loading_rows = []

for pc_index in (
    0,
    1,
):

    vector = (
        pca.components_[
            pc_index
        ]
        .reshape(
            37,
            3
        )
    )

    magnitude = np.linalg.norm(
        vector,
        axis=1
    )

    order = np.argsort(
        magnitude
    )[::-1]

    for rank, atom_index in enumerate(
        order,
        start=1
    ):

        loading_rows.append({
            "PC":
                pc_index + 1,
            "rank":
                rank,
            "descriptor_atom_index":
                int(
                    atom_index + 1
                ),
            "loading_x":
                float(
                    vector[
                        atom_index,
                        0
                    ]
                ),
            "loading_y":
                float(
                    vector[
                        atom_index,
                        1
                    ]
                ),
            "loading_z":
                float(
                    vector[
                        atom_index,
                        2
                    ]
                ),
            "loading_magnitude":
                float(
                    magnitude[
                        atom_index
                    ]
                ),
        })


loadings = pd.DataFrame(
    loading_rows
)

loadings.to_csv(
    OUT
    / "PC1_PC2_atom_loading_ranking.csv",
    index=False
)


# ============================================================
# 10. FINAL OUTPUTS
# ============================================================

df.to_csv(
    OUT
    / "k3_snapshot_assignments_with_physics.csv",
    index=False
)


result = {
    "final_model":
        {
            "k":
                recommended_k,
            "PCA_components":
                pca_components,
            "PCA_retained_variance":
                pca_retained,
            "PC1_variance_fraction":
                float(
                    selection[
                        "PC1_variance_fraction"
                    ]
                ),
            "PC2_variance_fraction":
                float(
                    selection[
                        "PC2_variance_fraction"
                    ]
                ),
            "PC1_PC2_variance_fraction":
                float(
                    selection[
                        "PC1_PC2_variance_fraction"
                    ]
                ),
        },
    "temperature_cluster_counts":
        {
            str(int(T)):
                {
                    str(int(c)):
                        int(
                            contingency
                            .loc[T, c]
                        )
                    for c in contingency.columns
                }
            for T in contingency.index
        },
    "temperature_cluster_fractions":
        {
            str(int(T)):
                {
                    str(int(c)):
                        float(
                            fractions
                            .loc[T, c]
                        )
                    for c in fractions.columns
                }
            for T in fractions.index
        },
    "temporal_membership":
        timeline_summary,
    "PC1_top_descriptor_atoms":
        (
            loadings[
                loadings["PC"] == 1
            ]
            .head(10)
            .to_dict(
                orient="records"
            )
        ),
    "PC2_top_descriptor_atoms":
        (
            loadings[
                loadings["PC"] == 2
            ]
            .head(10)
            .to_dict(
                orient="records"
            )
        ),
}


outfile = (
    OUT
    / "PHASE5_FINAL_K3_PHYSICAL_INTERPRETATION.json"
)

outfile.write_text(
    json.dumps(
        result,
        indent=2
    )
    + "\n"
)


print("="*90)
print("PHASE5-E33")
print("FINAL k=3 PHYSICAL CHARACTERIZATION")
print("="*90)

print()
print(
    f"FINAL_K={recommended_k}"
)

print(
    f"PCA_COMPONENTS={pca_components}"
)

print(
    f"PCA_RETAINED_VARIANCE="
    f"{pca_retained:.6f}"
)

print(
    f"PC1_PC2_VARIANCE="
    f"{selection['PC1_PC2_variance_fraction']:.6f}"
)

print()
print("CLUSTER PHYSICAL SUMMARY:")
print(
    summary.to_string(
        index=False
    )
)

print()
print("TEMPERATURE x CLUSTER COUNTS:")
print(
    contingency.to_string()
)

print()
print("TEMPERATURE x CLUSTER FRACTIONS:")
print(
    fractions
    .round(3)
    .to_string()
)

print()
print("PAIRWISE STANDARDIZED DIFFERENCES:")
print(
    effects.to_string(
        index=False
    )
)

print()
print("TEMPORAL INTERCONVERSIONS:")

for record in timeline_summary:

    print(
        f"{record['temperature_K']}K  "
        f"transitions="
        f"{record['cluster_transitions']}  "
        f"counts="
        f"{record['cluster_counts']}"
    )

print()
print("TOP PC1 DESCRIPTOR-ATOM LOADINGS:")
print(
    loadings[
        loadings["PC"] == 1
    ]
    .head(10)
    [
        [
            "rank",
            "descriptor_atom_index",
            "loading_magnitude",
        ]
    ]
    .to_string(
        index=False
    )
)

print()
print("TOP PC2 DESCRIPTOR-ATOM LOADINGS:")
print(
    loadings[
        loadings["PC"] == 2
    ]
    .head(10)
    [
        [
            "rank",
            "descriptor_atom_index",
            "loading_magnitude",
        ]
    ]
    .to_string(
        index=False
    )
)

print()
print(outfile)

print()
print("FINAL_K3_PHYSICAL_CHARACTERIZATION=PASS")
