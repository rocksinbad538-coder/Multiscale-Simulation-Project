from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT = Path("runs/phase1A/day015_lindblad_relaxation")

summary = pd.read_csv(OUT / "lindblad_relaxation_summary.csv")
thermal = pd.read_csv(OUT / "lindblad_expected_thermal_equilibrium.csv")
site_thermal = thermal[thermal["basis"] == "site"].set_index("state")["thermal_population"].to_dict()

site_labels = ["PYR2", "PYR3", "PYR4", "PYR5"]

def first_time_within_tolerance(df, expected, tol):
    for _, row in df.iterrows():
        ok = True
        for site in site_labels:
            if abs(row[site] - expected[site]) > tol:
                ok = False
                break
        if ok:
            return float(row["time_fs"])
    return np.nan

rows = []

for _, s in summary.iterrows():
    k = s["k_down_ps"]
    tag = f"kdown_{k:g}ps".replace(".", "p")
    csv = OUT / f"lindblad_relaxation_{tag}.csv"

    df = pd.read_csv(csv)

    t_thermal_5pct = first_time_within_tolerance(df, site_thermal, 0.05)
    t_thermal_2pct = first_time_within_tolerance(df, site_thermal, 0.02)
    t_thermal_1pct = first_time_within_tolerance(df, site_thermal, 0.01)

    rows.append({
        "k_down_ps": k,
        "gamma_phi_ps": s["gamma_phi_ps"],
        "max_PYR2": s["max_PYR2"],
        "final_PYR2": s["final_PYR2"],
        "t_PYR2_1pct_fs": s["t_PYR2_1pct_fs"],
        "t_PYR2_5pct_fs": s["t_PYR2_5pct_fs"],
        "t_PYR2_10pct_fs": s["t_PYR2_10pct_fs"],
        "t_thermal_within_5pct_fs": t_thermal_5pct,
        "t_thermal_within_2pct_fs": t_thermal_2pct,
        "t_thermal_within_1pct_fs": t_thermal_1pct,
    })

out = pd.DataFrame(rows)
out.to_csv(OUT / "lindblad_timescale_analysis.csv", index=False)

plt.figure(figsize=(6,4))
plt.plot(out["k_down_ps"], out["t_thermal_within_5pct_fs"], "o-", label="within 5%")
plt.plot(out["k_down_ps"], out["t_thermal_within_2pct_fs"], "o-", label="within 2%")
plt.plot(out["k_down_ps"], out["t_thermal_within_1pct_fs"], "o-", label="within 1%")
plt.xlabel("k_down (ps^-1)")
plt.ylabel("Thermalization time (fs)")
plt.legend()
plt.tight_layout()
plt.savefig(OUT / "lindblad_thermalization_time_vs_kdown.png", dpi=300)
plt.close()

interpretation = """# Day015 Lindblad Timescale Analysis

## Purpose

Quantify PYR2 access times and approximate thermalization times for the exciton-basis Lindblad relaxation model.

## Definition

Thermalization time is estimated as the first time at which all site populations are within a specified absolute tolerance of the expected 300 K thermal site populations.

## Interpretation

This analysis separates early transport into PYR2 from full relaxation toward the imposed excitonic thermal equilibrium.

## Caveat

The relaxation rates remain phenomenological. Therefore, these timescales should be interpreted as model-regime diagnostics rather than final physical times.
"""
(OUT / "LINDBLAD_TIMESCALE_ANALYSIS_DAY015.md").write_text(interpretation)

print(out.to_string(index=False))
print("Wrote:", OUT / "lindblad_timescale_analysis.csv")
