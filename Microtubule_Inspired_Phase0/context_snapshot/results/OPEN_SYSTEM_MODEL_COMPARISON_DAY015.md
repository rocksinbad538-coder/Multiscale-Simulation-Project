# Day015 Open-System Model Comparison

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
