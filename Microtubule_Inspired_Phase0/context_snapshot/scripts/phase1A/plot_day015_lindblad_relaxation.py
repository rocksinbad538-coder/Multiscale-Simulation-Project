from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

OUT = Path("runs/phase1A/day015_lindblad_relaxation")
df = pd.read_csv(OUT / "lindblad_relaxation_summary.csv")

plt.figure(figsize=(6,4))
plt.plot(df["k_down_ps"], df["max_PYR2"], "o-")
plt.xlabel("Exciton relaxation rate k_down (ps^-1)")
plt.ylabel("Maximum PYR2 population")
plt.tight_layout()
plt.savefig(OUT / "lindblad_max_PYR2_vs_kdown.png", dpi=300)
plt.close()

plt.figure(figsize=(6,4))
for col in ["t_PYR2_1pct_fs", "t_PYR2_5pct_fs", "t_PYR2_10pct_fs"]:
    plt.plot(df["k_down_ps"], df[col], "o-", label=col)
plt.xlabel("Exciton relaxation rate k_down (ps^-1)")
plt.ylabel("First arrival time (fs)")
plt.legend()
plt.tight_layout()
plt.savefig(OUT / "lindblad_PYR2_arrival_times_vs_kdown.png", dpi=300)
plt.close()

summary = """# Day015 Lindblad Exciton-Relaxation Interpretation

## Purpose

Introduce energetic relaxation in the exciton eigenbasis to move beyond pure-dephasing Haken-Strobl dynamics.

## Mathematical audit

The Lindblad implementation conserves trace and Hermiticity within numerical tolerance and maintains density-matrix positivity.

## Main result

Adding exciton-basis relaxation changes the long-time behavior relative to pure Haken-Strobl dephasing. The system no longer relaxes to uniform site population. Instead, the final population is biased toward lower-energy excitonic structure, with substantial PYR2 population.

## Physical interpretation

The current Lindblad model suggests that energetic relaxation can enhance access to PYR2, despite weak PYR2-PYR3 coupling. However, the relaxation rates used here are phenomenological. These results should be interpreted as a controlled model scan, not yet as a final quantitative prediction.

## Key caveat

Future refinement should estimate relaxation and dephasing rates from MD-derived fluctuations, literature-based spectral densities, or a calibrated open-system model.

## Files

- `lindblad_relaxation_summary.csv`
- `lindblad_max_PYR2_vs_kdown.png`
- `lindblad_PYR2_arrival_times_vs_kdown.png`
- `lindblad_input_hamiltonian_site_basis.csv`
- `lindblad_exciton_energies.csv`
- `lindblad_exciton_eigenvectors.csv`
"""
(OUT / "LINDBLAD_RELAXATION_INTERPRETATION_DAY015.md").write_text(summary)

print("Wrote Lindblad interpretation and figures to", OUT)
