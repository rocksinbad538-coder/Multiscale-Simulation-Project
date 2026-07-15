# MD Fluctuation to Bath-Parameter Proxy

## Purpose

Convert MD-derived site-energy or coupling time series into first-pass bath parameters for open-system exciton dynamics.

## Required input format

CSV with:

- `time_ps`
- one or more fluctuation columns, for example `PYR2`, `PYR3`, `PYR4`, `PYR5`, `J23`, `J34`, `J45`

## Outputs

For each selected column:

- fluctuation autocorrelation
- normalized autocorrelation
- correlation time
- spectral-density proxy
- dephasing-rate proxy

## Caveat

The dephasing estimate is a motional-narrowing proxy, not a final quantum bath model. It should be used as a diagnostic bridge between MD and Lindblad/Haken-Strobl parameters.
