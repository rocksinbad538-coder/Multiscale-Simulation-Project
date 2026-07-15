# Real MD-to-Bath Input Specification

## Purpose

Define the required input format for extracting bath parameters from real MD-derived fluctuation trajectories.

## Required CSV format

The bath-parameter pipeline expects a CSV file with the following structure:

```csv
time_ps,PYR2,PYR3,PYR4,PYR5,J23,J34,J45
0.000,...
0.010,...
0.020,...
```

## Required columns

| Column | Unit | Meaning |
|---|---|---|
| time_ps | ps | MD or sampled electronic-structure time. |
| PYR2 | eV | Site excitation energy of PYR2. |
| PYR3 | eV | Site excitation energy of PYR3. |
| PYR4 | eV | Site excitation energy of PYR4. |
| PYR5 | eV | Site excitation energy of PYR5. |
| J23 | eV | Coupling between PYR2 and PYR3. |
| J34 | eV | Coupling between PYR3 and PYR4. |
| J45 | eV | Coupling between PYR4 and PYR5. |

## How values should be generated

Preferred route:

1. Sample MD snapshots from the hydrated Phase 1A trajectory.
2. Extract local chromophore plus water environments.
3. Compute monomer site excitation energies for PYR2-PYR5.
4. Compute nearest-neighbor dimer couplings for PYR2-PYR3, PYR3-PYR4, and PYR4-PYR5.
5. Assemble the time series into the required CSV format.

## Minimum practical sampling

For a first-pass estimate:

- Electronic-structure sampling interval: 10-50 fs.
- Total sampled window: 5-20 ps.
- Number of samples: preferably at least 200.

## Quantities extracted by the pipeline

For each selected variable:

- Mean value.
- Standard deviation.
- Variance.
- Normalized autocorrelation function.
- Correlation time.
- Spectral-density proxy.
- Motional-narrowing dephasing proxy.

## Caveat

The current proxy does not replace a full quantum spectral-density model. It is intended as a controlled bridge between MD fluctuations and phenomenological open-system rates.
