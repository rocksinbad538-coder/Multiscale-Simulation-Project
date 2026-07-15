# Day015 Dynamic-Disorder Synthesis

## Scope

This document consolidates the synthetic dynamic-disorder analysis performed on the hydrated four-site PYR2-PYR5 excitonic Hamiltonian.

## Static reference

The reference model uses the validated hydrated Hamiltonian with Lindblad exciton relaxation.

Reference metrics:

- max PYR2 = 0.3056
- t(PYR2 > 30%) = 1350 fs

## Ensemble disorder scan

The first ensemble scan compared static dynamics, site-energy disorder, coupling disorder, and combined site/coupling disorder.

| Case | max PYR2 | final PYR2 | t(PYR2 > 30%) |
|---|---:|---:|---:|
| static | 0.306 | 0.305 | 1350 fs |
| site disorder 5 meV | 0.318 | 0.311 | 1204 fs |
| site disorder 10 meV | 0.332 | 0.315 | 935 fs |
| coupling disorder 2 meV | 0.318 | 0.304 | 1390 fs |
| site 10 + coupling 2 meV | 0.345 | 0.309 | 981 fs |

Main result: site-energy disorder enhances PYR2 access more strongly than coupling disorder alone.

## Sigma-tau map

A 2D map over site-energy disorder amplitude and correlation time identified a robust enhancement region:

- sigma_E ≈ 10-20 meV
- tau_c ≈ 0.02-0.10 ps

The best sigma-tau case was:

- sigma_E = 20 meV
- tau_c = 0.02 ps
- max PYR2 = 0.3418
- t(PYR2 > 30%) = 811 fs
- speedup = 1.66x relative to the static model

## Gamma-sigma map

A second map varied additional Lindblad pure dephasing and site-energy disorder amplitude at fixed tau_c = 0.05 ps.

The best case by fastest t(PYR2 > 30%) was:

- gamma_phi = 0.1 ps^-1
- sigma_E = 20 meV
- max PYR2 = 0.3428
- t(PYR2 > 30%) = 717.5 fs
- speedup = 1.88x relative to the static model

Importantly, enhancement was also present at gamma_phi = 0 ps^-1:

- gamma_phi = 0 ps^-1
- sigma_E = 20 meV
- max PYR2 = 0.3433
- t(PYR2 > 30%) = 732.5 fs
- speedup = 1.84x

## Physical interpretation

The hydrated excitonic chain appears to support an ENAQT-like regime. However, the dominant mechanism is not additional Markovian pure dephasing. Instead, explicit dynamic site-energy fluctuations provide the main transport-assisting mechanism.

The weak PYR2-PYR3 coupling remains the bottleneck, but moderate fluctuations of approximately 10-20 meV can transiently reduce energetic mismatch and improve access to PYR2.

## Caveat

All disorder in this block is synthetic Ornstein-Uhlenbeck disorder. These results should be interpreted as regime identification, not final quantitative prediction.

## Next physical step

Replace synthetic disorder parameters with MD-derived site-energy and coupling fluctuation statistics:

1. Extract site-energy/coupling time series from MD-sampled chromophore-water configurations.
2. Compute variance, autocorrelation time and spectral-density proxies.
3. Determine whether real fluctuations fall inside the identified enhancement regime.
4. Re-run the dynamic-disorder Lindblad model using MD-derived parameters.
