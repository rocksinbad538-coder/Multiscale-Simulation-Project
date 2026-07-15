#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

root = Path(".")
outdir = root / "figures" / "day003"
outdir.mkdir(parents=True, exist_ok=True)

figs = [
    ("figures/thermo/bn_like_30000w_with_hold_water_temperature.png", "Water temperature"),
    ("figures/thermo/bn_like_30000w_with_hold_water_msd.png", "Water MSD"),
    ("figures/confinement/nvt300_hold/bn_like_30000w_fraction_inside_lumen_segment.png", "Fraction inside lumen"),
    ("results/profiles/nvt300_hold_radial_number_density.png", "Volume-corrected radial density"),
]

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

for ax, (path, title) in zip(axes.flatten(), figs):
    p = root / path
    if not p.exists():
        ax.text(0.5, 0.5, f"Missing:\n{path}", ha="center", va="center")
        ax.axis("off")
        continue
    img = mpimg.imread(p)
    ax.imshow(img)
    ax.set_title(title, fontsize=12)
    ax.axis("off")

fig.suptitle(
    "Day 003 — BN-like scaffold + 30,000 confined waters: 300 K stability evidence",
    fontsize=16,
)
plt.tight_layout()
outfile = outdir / "day003_bn_like_30000w_stability_dashboard.png"
plt.savefig(outfile, dpi=200)
plt.close()

print(f"Wrote {outfile}")
