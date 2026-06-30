# Day019 full two-state embedded-manifold analysis

## Validation

- Embedded outputs selected: 84/84
- S1/S2 state observations parsed: 168/168
- Bright-root mapping preserved: 84/84
- Bright-root switches detected: 0

## Electronic-character convention

| Site | Bright-like local state | Alternate-like local state |
|---|---:|---:|
| PYR2 | S2 | S1 |
| PYR3 | S2 | S1 |
| PYR4 | S2 | S1 |
| PYR5 | S1 | S2 |

The family labels are defined by oscillator-strength ordering within S1/S2 and are consistent with the Day019 NTO occupation and cross-site subspace analyses.

## Aggregate results

- Full S1-S2 gap range: 53.000 to 102.000 meV
- Mean S1-S2 gap: 74.881 meV
- Mean bright-state share of S1/S2 oscillator strength: 0.978472
- Minimum bright-state share of S1/S2 oscillator strength: 0.878722

## Site-resolved statistics

| Site | Family | Root | Mean energy (eV) | SD (meV) | Energy range (meV) | Mean fosc |
|---|---|---:|---:|---:|---:|---:|
| PYR2 | bright_like | S2 | 4.088571 | 10.817 | 36.000 | 0.493290 |
| PYR2 | alternate_like | S1 | 4.006143 | 8.061 | 27.000 | 0.011173 |
| PYR3 | bright_like | S2 | 4.078048 | 11.615 | 44.000 | 0.501778 |
| PYR3 | alternate_like | S1 | 4.003190 | 13.037 | 50.000 | 0.009893 |
| PYR4 | bright_like | S2 | 4.090619 | 7.403 | 35.000 | 0.510859 |
| PYR4 | alternate_like | S1 | 4.019238 | 11.820 | 57.000 | 0.008602 |
| PYR5 | bright_like | S1 | 3.776952 | 16.034 | 47.000 | 0.648089 |
| PYR5 | alternate_like | S2 | 3.847810 | 14.949 | 56.000 | 0.018624 |

## Model-space conclusion

**The two local electronic families are preserved across all 84 embedded calculations.** The eight-state basis can therefore be indexed by electronic character rather than by a globally fixed root number. The four-state model is the bright-like subset of this basis.

This analysis provides the complete diagonal-energy and oscillator-strength time series for the two-state-per-site manifold. Interstate and intersite couplings remain the next required ingredient before constructing the final excitonic Hamiltonian.
