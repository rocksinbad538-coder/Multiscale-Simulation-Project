#!/usr/bin/env python3

from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

META = ROOT/"runs"/"phase2"/"campaign"/"campaign_meta_summary.csv"

OUT = ROOT/"runs"/"phase2"/"campaign"/"descriptor_correlations"

OUT.mkdir(exist_ok=True)

df = pd.read_csv(META)

cols = [
    "final_Rg_A",
    "mean_RMSF_A",
    "maximum_RMSF_A",
    "mean_aligned_rmsd_A",
    "mean_relative_shape_anisotropy",
    "maximum_atomic_displacement_A",
]

corr = df[cols].corr()

corr.to_csv(OUT/"descriptor_correlation_matrix.csv")

plt.figure(figsize=(7,6))

im = plt.imshow(
    corr,
    vmin=-1,
    vmax=1,
    interpolation="nearest"
)

plt.xticks(range(len(cols)), cols, rotation=90, fontsize=8)
plt.yticks(range(len(cols)), cols, fontsize=8)

plt.colorbar(im, label="Pearson correlation")

plt.tight_layout()

plt.savefig(
    OUT/"CorrelationHeatmap.png",
    dpi=300
)

plt.close()

print("="*90)
print("DAY049 / PHASE2-B11")
print("DESCRIPTOR CORRELATIONS")
print("="*90)

print(corr)

print()

print(OUT)
