# Day015 Open-System Dynamics Report

## Scope

This report consolidates the first open-system dynamics analysis for the hydrated four-site PYR2-PYR5 excitonic Hamiltonian derived from explicit-water TDDFT calculations.

## Hydrated Hamiltonian

The analysis uses the validated Day014 hydrated Hamiltonian with nearest-neighbor effective couplings:

| Pair | Coupling |
|---|---:|
| PYR2-PYR3 | 2.37 meV |
| PYR3-PYR4 | 17.00 meV |
| PYR4-PYR5 | 33.00 meV |

The weak PYR2-PYR3 coupling creates the primary transport bottleneck.

## Haken-Strobl dynamics

Pure-dephasing Haken-Strobl simulations were used as a diagnostic model.

For gamma = 10 ps^-1 and initial excitation on PYR5:

| Metric | Value |
|---|---:|
| t(PYR2 > 1%) | 80 fs |
| t(PYR2 > 5%) | 760 fs |
| t(PYR2 > 10%) | 1810 fs |
| max PYR2 | 0.234 |

The long-time uniform population produced by pure dephasing is interpreted as a model artifact, not as thermodynamic equilibration.

## J23 sensitivity

Sensitivity analysis confirmed that PYR2 access is strongly controlled by the PYR2-PYR3 coupling.

At gamma = 10 ps^-1:

| J23 | max PYR2 | t(PYR2 > 10%) |
|---:|---:|---:|
| 1.00 meV | 0.098 | not reached |
| 2.00 meV | 0.215 | 2550 fs |
| 2.37 meV | 0.234 | 1810 fs |
| 5.00 meV | 0.250 | 400 fs |
| 10.00 meV | 0.250 | 70 fs |
| 20.00 meV | 0.333 | 50 fs |

## PYR2 detuning sensitivity

With fixed J23 = 2.37 meV, moderate energetic detuning of PYR2 changes arrival times but does not dominate transport as strongly as J23.

## Lindblad exciton-relaxation model

A Lindblad model was implemented in the exciton eigenbasis with detailed-balance-compatible upward rates at 300 K.

Expected thermal exciton populations:

| Exciton | Population |
|---|---:|
| X1 | 0.536772 |
| X2 | 0.306374 |
| X3 | 0.126535 |
| X4 | 0.030319 |

Projected thermal site populations:

| Site | Population |
|---|---:|
| PYR2 | 0.305432 |
| PYR3 | 0.160919 |
| PYR4 | 0.283158 |
| PYR5 | 0.250492 |

For k_down >= 1 ps^-1, simulated final populations matched the expected thermal populations with errors near 1e-9 or smaller.

## Lindblad timescales

| k_down | t(PYR2 > 10%) | thermalization within 1% |
|---:|---:|---:|
| 0.1 ps^-1 | 1340 fs | not reached in 10 ps |
| 1.0 ps^-1 | 130 fs | 1200 fs |
| 5.0 ps^-1 | 30 fs | 240 fs |
| 10.0 ps^-1 | 20 fs | 190 fs |

## Initial-condition robustness

For k_down = 1 ps^-1 and gamma_phi = 1 ps^-1, all initial site-localized excitations converge to the same thermal site distribution.

Access to PYR2 from non-PYR2 initial states:

| Initial site | t(PYR2 > 10%) | t(PYR2 > 20%) | t(PYR2 > 30%) |
|---|---:|---:|---:|
| PYR3 | 110 fs | 310 fs | 950 fs |
| PYR4 | 150 fs | 410 fs | 1530 fs |
| PYR5 | 130 fs | 390 fs | 1350 fs |

## Main physical conclusions

1. The hydrated excitonic chain is not homogeneous.
2. PYR2-PYR3 is the dominant transport bottleneck.
3. Pure dephasing assists access to PYR2 but produces an artificial uniform long-time limit.
4. Exciton-basis Lindblad relaxation produces a nonuniform thermal equilibrium.
5. PYR2 becomes significantly populated at equilibrium because the low-energy excitonic manifold has substantial PYR2 character.
6. The open-system dynamics are internally consistent and ready for rate calibration from MD fluctuations or spectral-density modeling.

## Caveat

The current Lindblad relaxation rates are phenomenological. Quantitative physical prediction requires future calibration from MD-derived site-energy/coupling fluctuations, spectral-density assumptions, or literature-based relaxation/dephasing timescales.
