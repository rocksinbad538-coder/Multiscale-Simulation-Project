# Day015 Haken-Strobl Physical Interpretation

## Mathematical audit

The Haken-Strobl implementation is numerically consistent for the hydrated four-site Hamiltonian. Trace, Hermiticity, and population are conserved to numerical precision. Positivity is maintained within numerical tolerance.

## Physical interpretation

The long-time uniform population produced by the pure-dephasing Haken-Strobl model is not interpreted as thermodynamic equilibration. It is a model artifact arising from pure dephasing without energetic relaxation or detailed balance.

The early-time transport metrics are physically more informative.

## Key result

Starting from an initial excitation on PYR5, access to PYR2 is strongly dephasing-dependent:

- gamma = 0 ps^-1: maximum PYR2 population is only ~0.022.
- gamma = 5 ps^-1: maximum PYR2 population increases to ~0.194.
- gamma = 10 ps^-1: maximum PYR2 population increases to ~0.234.
- gamma = 50 ps^-1: PYR2 approaches the uniform-population limit.

This indicates noise-assisted access to PYR2 in the current four-site hydrated Hamiltonian.

## Scientific caveat

This should not yet be interpreted as a final physical transport prediction. A more physical model should include energetic relaxation, temperature or detailed balance, possible recombination/trapping channels, and environmental spectral information from MD or literature-based assumptions.

## Files

- haken_strobl_audit_summary.csv
- HAKEN_STROBL_AUDIT_DAY015.md
- audit_max_PYR2_vs_gamma.png
- audit_PYR2_arrival_times_vs_gamma.png
