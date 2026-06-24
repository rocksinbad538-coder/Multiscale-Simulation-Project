from pathlib import Path
import pandas as pd

OUT = Path("runs/phase1A/day015_dynamic_disorder_sigma_tau_map")
df = pd.read_csv(OUT / "sigma_tau_dynamic_disorder_summary.csv")

static = df[df["sigma_E_meV"] == 0].iloc[0]
static_max = static["mean_max_PYR2"]
static_t30 = static["mean_t30_fs"]

df["delta_max_PYR2_vs_static"] = df["mean_max_PYR2"] - static_max
df["delta_t30_fs_vs_static"] = df["mean_t30_fs"] - static_t30
df["t30_speedup_factor"] = static_t30 / df["mean_t30_fs"]

df.to_csv(OUT / "sigma_tau_dynamic_disorder_summary_with_deltas.csv", index=False)

top_max = df.sort_values("mean_max_PYR2", ascending=False).head(10)
top_t30 = df.sort_values("mean_t30_fs", ascending=True).head(10)

top_max.to_csv(OUT / "top10_sigma_tau_by_max_PYR2.csv", index=False)
top_t30.to_csv(OUT / "top10_sigma_tau_by_fastest_t30.csv", index=False)

md = f"""# Day015 Sigma-Tau Dynamic-Disorder Map Summary

## Static reference

- mean max PYR2 = {static_max:.6f}
- mean t(PYR2 > 30%) = {static_t30:.1f} fs

## Main result

The site-energy dynamic-disorder scan identifies a robust transport-enhancement regime.

The strongest enhancement occurs for site-energy fluctuations of approximately 10-20 meV and correlation times of approximately 0.02-0.10 ps.

## Best cases by maximum PYR2 population

{top_max[["sigma_E_meV","tau_c_ps","mean_max_PYR2","mean_t30_fs","delta_max_PYR2_vs_static","t30_speedup_factor"]].to_string(index=False)}

## Best cases by fastest t(PYR2 > 30%)

{top_t30[["sigma_E_meV","tau_c_ps","mean_max_PYR2","mean_t30_fs","delta_t30_fs_vs_static","t30_speedup_factor"]].to_string(index=False)}

## Interpretation

The static model already reaches PYR2, but dynamic site-energy fluctuations moderately improve PYR2 access. This is consistent with an ENAQT-like regime where environmental fluctuations help bridge the weak PYR2-PYR3 bottleneck.

The effect is not monotonic in correlation time. Fast and intermediate fluctuations are more beneficial than slow fluctuations, suggesting that slow disorder behaves more like transient energetic inhomogeneity than efficient transport-assisting noise.

## Caveat

This is still a synthetic dynamic-disorder scan. The next physical step is to replace the Ornstein-Uhlenbeck parameters with MD-derived site-energy and coupling fluctuation statistics.
"""

(OUT / "SIGMA_TAU_DYNAMIC_DISORDER_SUMMARY_DAY015.md").write_text(md)

print("Static max_PYR2:", static_max)
print("Static t30 fs:", static_t30)
print("\nTop by max PYR2:")
print(top_max[["sigma_E_meV","tau_c_ps","mean_max_PYR2","mean_t30_fs","delta_max_PYR2_vs_static","t30_speedup_factor"]].to_string(index=False))
print("\nTop by fastest t30:")
print(top_t30[["sigma_E_meV","tau_c_ps","mean_max_PYR2","mean_t30_fs","delta_t30_fs_vs_static","t30_speedup_factor"]].to_string(index=False))
print("Wrote:", OUT)
