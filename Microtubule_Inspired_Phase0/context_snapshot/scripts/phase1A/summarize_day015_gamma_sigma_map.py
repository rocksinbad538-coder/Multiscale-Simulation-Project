from pathlib import Path
import pandas as pd

OUT = Path("runs/phase1A/day015_dynamic_disorder_gamma_sigma_map")
df = pd.read_csv(OUT / "gamma_sigma_dynamic_disorder_summary.csv")

static = df[(df["gamma_phi_ps"] == 1.0) & (df["sigma_E_meV"] == 0)].iloc[0]
static_max = static["mean_max_PYR2"]
static_t30 = static["mean_t30_fs"]

df["delta_max_PYR2_vs_static"] = df["mean_max_PYR2"] - static_max
df["delta_t30_fs_vs_static"] = df["mean_t30_fs"] - static_t30
df["t30_speedup_factor"] = static_t30 / df["mean_t30_fs"]

df.to_csv(OUT / "gamma_sigma_dynamic_disorder_summary_with_deltas.csv", index=False)

top_max = df.sort_values("mean_max_PYR2", ascending=False).head(10)
top_t30 = df.sort_values("mean_t30_fs", ascending=True).head(10)

top_max.to_csv(OUT / "top10_gamma_sigma_by_max_PYR2.csv", index=False)
top_t30.to_csv(OUT / "top10_gamma_sigma_by_fastest_t30.csv", index=False)

md = f"""# Day015 Gamma-Sigma Dynamic-Disorder Map Summary

## Static reference

Reference case: gamma_phi = 1 ps^-1, sigma_E = 0 meV.

- mean max PYR2 = {static_max:.6f}
- mean t(PYR2 > 30%) = {static_t30:.1f} fs

## Main result

The map shows that transport enhancement is driven primarily by site-energy dynamic disorder amplitude rather than by additional Lindblad pure dephasing.

The strongest enhancement occurs for sigma_E around 20 meV. Enhancement is already present even when gamma_phi = 0 ps^-1, indicating that dynamic Hamiltonian fluctuations themselves provide an effective decoherence/modulation pathway.

## Best cases by maximum PYR2 population

{top_max[["gamma_phi_ps","sigma_E_meV","mean_max_PYR2","mean_t30_fs","delta_max_PYR2_vs_static","t30_speedup_factor"]].to_string(index=False)}

## Best cases by fastest t(PYR2 > 30%)

{top_t30[["gamma_phi_ps","sigma_E_meV","mean_max_PYR2","mean_t30_fs","delta_t30_fs_vs_static","t30_speedup_factor"]].to_string(index=False)}

## Interpretation

The result supports an ENAQT-like mechanism, but specifically one dominated by site-energy fluctuations. In this regime, moderate dynamic energetic disorder helps bridge the weak PYR2-PYR3 bottleneck.

## Caveat

The disorder remains synthetic. The next physical step is to determine whether real MD-derived site-energy fluctuations fall in the 10-20 meV amplitude range with sub-100 fs to few-hundred-fs correlation times.
"""
(OUT / "GAMMA_SIGMA_DYNAMIC_DISORDER_SUMMARY_DAY015.md").write_text(md)

print("Static max_PYR2:", static_max)
print("Static t30 fs:", static_t30)
print("\nTop by max PYR2:")
print(top_max[["gamma_phi_ps","sigma_E_meV","mean_max_PYR2","mean_t30_fs","delta_max_PYR2_vs_static","t30_speedup_factor"]].to_string(index=False))
print("\nTop by fastest t30:")
print(top_t30[["gamma_phi_ps","sigma_E_meV","mean_max_PYR2","mean_t30_fs","delta_t30_fs_vs_static","t30_speedup_factor"]].to_string(index=False))
print("Wrote:", OUT)
