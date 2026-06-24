from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

OUT = Path("runs/phase1A/day015_open_system_comparison")
OUT.mkdir(parents=True, exist_ok=True)

HS = Path("runs/phase1A/day015_exciton_model/haken_strobl_audit_summary.csv")
LB = Path("runs/phase1A/day015_lindblad_relaxation/lindblad_relaxation_summary.csv")

hs = pd.read_csv(HS)
lb = pd.read_csv(LB)

# Representative Haken-Strobl case: gamma = 10 ps^-1
hs10 = hs.loc[hs["gamma_ps"] == 10].iloc[0]

rows = []

rows.append({
    "model": "Haken-Strobl",
    "parameter": "gamma=10 ps^-1",
    "final_PYR2": hs10["final_PYR2"],
    "final_PYR3": hs10["final_PYR3"],
    "final_PYR4": hs10["final_PYR4"],
    "final_PYR5": hs10["final_PYR5"],
    "max_PYR2": hs10["max_PYR2"],
    "t_PYR2_1pct_fs": hs10["t_PYR2_1pct_fs"],
    "t_PYR2_5pct_fs": hs10["t_PYR2_5pct_fs"],
    "t_PYR2_10pct_fs": hs10["t_PYR2_10pct_fs"],
    "interpretation": "pure dephasing; long-time uniform population is a model artifact",
})

for _, r in lb.iterrows():
    rows.append({
        "model": "Lindblad exciton relaxation",
        "parameter": f"k_down={r['k_down_ps']:g} ps^-1, gamma_phi={r['gamma_phi_ps']:g} ps^-1",
        "final_PYR2": r["final_PYR2"],
        "final_PYR3": r["final_PYR3"],
        "final_PYR4": r["final_PYR4"],
        "final_PYR5": r["final_PYR5"],
        "max_PYR2": r["max_PYR2"],
        "t_PYR2_1pct_fs": r["t_PYR2_1pct_fs"],
        "t_PYR2_5pct_fs": r["t_PYR2_5pct_fs"],
        "t_PYR2_10pct_fs": r["t_PYR2_10pct_fs"],
        "interpretation": "exciton-basis relaxation with detailed-balance structure",
    })

cmp = pd.DataFrame(rows)
cmp.to_csv(OUT / "open_system_model_comparison_day015.csv", index=False)

plt.figure(figsize=(7,4))
plt.bar(range(len(cmp)), cmp["final_PYR2"])
plt.xticks(range(len(cmp)), cmp["parameter"], rotation=45, ha="right")
plt.ylabel("Final PYR2 population")
plt.tight_layout()
plt.savefig(OUT / "comparison_final_PYR2.png", dpi=300)
plt.close()

plt.figure(figsize=(7,4))
plt.bar(range(len(cmp)), cmp["t_PYR2_10pct_fs"])
plt.xticks(range(len(cmp)), cmp["parameter"], rotation=45, ha="right")
plt.ylabel("t(PYR2 > 10%) (fs)")
plt.tight_layout()
plt.savefig(OUT / "comparison_PYR2_10pct_arrival.png", dpi=300)
plt.close()

summary = """# Day015 Open-System Model Comparison

## Purpose

Compare pure-dephasing Haken-Strobl dynamics against an exciton-basis Lindblad relaxation model for the hydrated four-site PYR2-PYR5 Hamiltonian.

## Main conclusion

The Haken-Strobl model is mathematically valid but produces long-time uniform populations because it includes pure dephasing without energetic relaxation. This long-time limit is therefore not interpreted as physical equilibration.

The Lindblad exciton-relaxation model changes the long-time behavior by introducing energetic relaxation and detailed-balance-like upward rates. This produces nonuniform final site populations and stronger accumulation on PYR2 than Haken-Strobl pure dephasing.

## Physical interpretation

The combined analysis supports three conclusions:

1. PYR2 access is limited by the weak PYR2-PYR3 coupling.
2. Dephasing assists access to PYR2 across this bottleneck.
3. Energetic relaxation further redistributes population toward lower-energy excitonic structure, increasing PYR2 occupation.

## Caveat

The Lindblad rates remain phenomenological. Quantitative prediction requires rate calibration from MD-derived fluctuations, spectral-density assumptions, or literature-based relaxation times.

## Files

- `open_system_model_comparison_day015.csv`
- `comparison_final_PYR2.png`
- `comparison_PYR2_10pct_arrival.png`
"""
(OUT / "OPEN_SYSTEM_MODEL_COMPARISON_DAY015.md").write_text(summary)

print(cmp.to_string(index=False))
print("Wrote:", OUT)
