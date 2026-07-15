# Day019 point-transition-dipole coupling baseline

## Scope

- Character-indexed eight-state local basis.
- Six intersite pairs Ã four state-family combinations Ã 21 frames = 504 couplings.
- Local-state transition-dipole signs were gauge-aligned to the corresponding frame000 vector before evaluating signed couplings.
- Same-site alternate/bright electronic couplings are zero in the local adiabatic basis used here. Time-derivative nonadiabatic couplings are not included.

## Numerical controls

- Transition-dipole observations: 168/168
- Intersite coupling observations: 504/504
- Hamiltonian snapshots: 21/21
- Relative-permittivity divisor: 1.000000
- Gauge-stability groups passing cosine >= 0.80: 4/8
- Minimum aligned cosine to frame000: 0.011738
- Minimum site-centroid distance: 26.894353 Ã
- Minimum interatomic distance: 20.686104 Ã

## Coupling magnitudes

- Overall mean |J|: 0.228796 meV
- Overall maximum |J|: 1.424239 meV
- Overall signed-J SD: 0.461670 meV
- Bright-bright mean |J|: 0.733840 meV
- Bright-bright maximum |J|: 1.424239 meV
- Mixed-family mean |J|: 0.085615 meV
- Alternate-alternate mean |J|: 0.010112 meV

## Scale comparison

- Minimum local S1-S2 gap: 53.000 meV
- Maximum diagonal energy SD: 16.034 meV
- max|J| / minimum local gap: 0.026872
- max|J| / maximum diagonal SD: 0.088826

## Pair-resolved statistics

| Site A | Family A | Site B | Family B | Distance (Ã) | Mean J (meV) | SD (meV) | Mean |J| (meV) | Max |J| (meV) |
|---|---|---|---|---:|---:|---:|---:|---:|
| PYR2 | alternate_like | PYR3 | alternate_like | 26.896125 | 0.014388 | 0.019151 | 0.016449 | 0.070391 |
| PYR2 | alternate_like | PYR3 | bright_like | 26.896125 | -0.109243 | 0.083983 | 0.112048 | 0.344134 |
| PYR2 | alternate_like | PYR4 | alternate_like | 42.490545 | -0.001566 | 0.002683 | 0.002144 | 0.008874 |
| PYR2 | alternate_like | PYR4 | bright_like | 42.490545 | 0.021267 | 0.025384 | 0.027323 | 0.085100 |
| PYR2 | alternate_like | PYR5 | alternate_like | 46.893453 | -0.001047 | 0.004903 | 0.004146 | 0.010049 |
| PYR2 | alternate_like | PYR5 | bright_like | 46.893453 | 0.003645 | 0.018946 | 0.017220 | 0.036782 |
| PYR2 | bright_like | PYR3 | alternate_like | 26.896125 | 0.138906 | 0.087414 | 0.141239 | 0.352074 |
| PYR2 | bright_like | PYR3 | bright_like | 26.896125 | -1.169367 | 0.026070 | 1.169367 | 1.201115 |
| PYR2 | bright_like | PYR4 | alternate_like | 42.490545 | -0.018231 | 0.016652 | 0.019684 | 0.049329 |
| PYR2 | bright_like | PYR4 | bright_like | 42.490545 | 0.256270 | 0.005039 | 0.256270 | 0.260616 |
| PYR2 | bright_like | PYR5 | alternate_like | 46.893453 | -0.022890 | 0.019738 | 0.028417 | 0.045972 |
| PYR2 | bright_like | PYR5 | bright_like | 46.893453 | 0.075082 | 0.004044 | 0.075082 | 0.081270 |
| PYR3 | alternate_like | PYR4 | alternate_like | 26.894353 | 0.013871 | 0.015422 | 0.014792 | 0.045465 |
| PYR3 | alternate_like | PYR4 | bright_like | 26.894353 | -0.151984 | 0.096394 | 0.157116 | 0.385945 |
| PYR3 | alternate_like | PYR5 | alternate_like | 43.421130 | 0.002945 | 0.003768 | 0.003517 | 0.010097 |
| PYR3 | alternate_like | PYR5 | bright_like | 43.421130 | -0.028229 | 0.020785 | 0.028662 | 0.075924 |
| PYR3 | bright_like | PYR4 | alternate_like | 26.894353 | -0.137367 | 0.090635 | 0.140787 | 0.334273 |
| PYR3 | bright_like | PYR4 | bright_like | 26.894353 | 1.228828 | 0.017191 | 1.228828 | 1.266366 |
| PYR3 | bright_like | PYR5 | alternate_like | 43.421130 | -0.028138 | 0.026426 | 0.032425 | 0.077104 |
| PYR3 | bright_like | PYR5 | bright_like | 43.421130 | 0.285290 | 0.003780 | 0.285290 | 0.291032 |
| PYR4 | alternate_like | PYR5 | alternate_like | 27.639866 | -0.009983 | 0.021004 | 0.019625 | 0.053120 |
| PYR4 | alternate_like | PYR5 | bright_like | 27.639866 | 0.160185 | 0.105929 | 0.166146 | 0.380607 |
| PYR4 | bright_like | PYR5 | alternate_like | 27.639866 | 0.061097 | 0.173255 | 0.156318 | 0.406456 |
| PYR4 | bright_like | PYR5 | bright_like | 27.639866 | -1.388204 | 0.026260 | 1.388204 | 1.424239 |

## Interpretation boundary

These are point-transition-dipole couplings and constitute a controlled baseline, not a transition-density benchmark. The minimum interatomic separation exceeds 20 Ã, which supports using the dipolar term as a first approximation; nevertheless, finite-size multipolar and dielectric-screening effects remain unquantified. Coupling signs are reported in the explicit frame000 local-state gauge and are not independently observable under arbitrary local basis-phase changes.
