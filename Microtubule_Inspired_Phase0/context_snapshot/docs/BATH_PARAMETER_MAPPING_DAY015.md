# Day015 Bath-Parameter Mapping

## Purpose

Prepare the bridge between the current phenomenological open-system model and physically calibrated bath parameters.

## Context

The current Lindblad model is internally consistent and converges to the imposed 300 K thermal equilibrium. However, the dephasing and relaxation rates are still phenomenological. The next model-improvement step is to derive or constrain these rates from MD fluctuations and spectral-density analysis.

## Current physical scales

At 300 K:

- kBT = 25.85 meV
- J23 = 2.37 meV
- J34 = 17.00 meV
- J45 = 33.00 meV

Therefore:

- J23 << kBT
- J34 is comparable to kBT
- J45 is slightly larger than kBT

## Interpretation

The weak J23 coupling explains the PYR2 access bottleneck. The stronger J34 and J45 couplings support faster redistribution within the PYR3-PYR5 subnetwork.

The current Lindblad relaxation scan should be interpreted as a controlled phenomenological model. To make it physically predictive, the following quantities should be extracted from MD or calibrated from literature:

1. Site-energy fluctuation time series.
2. Coupling fluctuation time series.
3. Autocorrelation functions.
4. Spectral densities.
5. Reorganization energies.
6. Dephasing rates.
7. Relaxation rates satisfying detailed balance.

## Immediate next computational target

Build an MD-to-bath-parameter workflow template:

- read chromophore-resolved site-energy proxy time series;
- compute fluctuation autocorrelation C(t);
- estimate correlation time;
- estimate dephasing scale;
- generate a first spectral-density proxy.

This will connect the current Lindblad model to physically motivated bath parameters.
