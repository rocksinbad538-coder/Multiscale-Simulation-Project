# Day019 bright-model readiness audit

## Input model

- Four bright local states: PYR2, PYR3, PYR4, and PYR5.
- TDC-AC finite-size-corrected bright-bright couplings.
- Snapshots: 21/21.
- Sampling interval: 5.000 ps.
- Chromophore geometries are frozen; the time dependence originates from the solvent embedding.

## Diagonal disorder

| State | Mean energy (eV) | SD (meV) | Range (meV) |
|---|---:|---:|---:|
| PYR2_bright | 0.338571 | 10.817 | 36.000 |
| PYR3_bright | 0.328048 | 11.615 | 44.000 |
| PYR4_bright | 0.340619 | 7.403 | 35.000 |
| PYR5_bright | 0.026952 | 16.034 | 47.000 |

## Coupling and detuning scales

- Overall maximum |J|: 1.492201 meV.
- Diagonal-energy SD range: 7.403-16.034 meV.
- Minimum sampled absolute site detuning: 1.000000 meV at frame 002 for PYR2_bright-PYR4_bright.
- Maximum sampled |J|/|Delta|: 1.228157 at frame 018 for PYR3_bright-PYR4_bright.

## Eigenstate localization

- Participation-ratio range: 1.000034-2.320444.
- Mean participation ratio: 1.093888.
- Minimum PYR5 population in the lowest eigenstate: 0.99997235.
- Minimum adjacent eigenvalue spacing: 0.986980 meV.
- Full sampled eigenvalue range: 354.095407 meV.

## Time-resolution audit

- Electronic coupling timescale hbar/max|J|: 0.441101 ps.
- Snapshot interval / coupling timescale: 11.335.
- Lag-1 (5.0 ps) autocorrelation for PYR2_bright: -0.088596.
- Lag-1 (5.0 ps) autocorrelation for PYR3_bright: -0.417922.
- Lag-1 (5.0 ps) autocorrelation for PYR4_bright: -0.286243.
- Lag-1 (5.0 ps) autocorrelation for PYR5_bright: +0.071819.

- The autocorrelations are descriptive only: 21 points are insufficient for a converged spectral density.
- The 5 ps sampling interval is substantially longer than the electronic coupling timescale.

## Readiness decision

- Static-disorder ensemble: READY.
- Continuous stochastic Hamiltonian propagation: NOT READY.
- Bath autocorrelation/spectral-density extraction: NOT READY.

## Accepted use of the present dataset

The 21 Hamiltonians can be used as a quasi-static disorder ensemble: each snapshot defines one frozen Hamiltonian realization for independent coherent propagation, ensemble averaging, eigenstate analysis, or comparison with phenomenological dephasing models. They should not be connected sequentially as a continuous 5 ps-resolved stochastic trajectory.

## Required data for dynamical bath models

A defensible time-dependent bath treatment requires a substantially finer embedding sampling interval and a longer trajectory. The required interval must resolve both the sub-picosecond electronic coupling scale and the relevant solvent fluctuations. The current 21-point series cannot determine a reliable spectral density, memory kernel, or microscopic pure-dephasing rate.
