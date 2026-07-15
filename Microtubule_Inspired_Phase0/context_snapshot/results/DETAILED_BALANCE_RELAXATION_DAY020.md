# Day020 Detailed-Balance Relaxation Sensitivity

## Scope

This calculation introduces population relaxation separately from the previous pure-dephasing model.

The jump operators are constructed in the instantaneous energy eigenbasis of each bright-state Hamiltonian. Upward and downward rates satisfy detailed balance at the selected temperature.

The relaxation amplitude `kappa_ref` is phenomenological and is not derived from a microscopic spectral density.

## Parameter sweep

- Temperatures: 150, 200, 250, 300 K.
- Kappa values: 0.1, 1, 10 ps^-1.
- Propagation interval: 0-100 ps.
- Output interval: 0.05 ps.

## Numerical validation

- Maximum trace error: 9.208e-13.
- Maximum Hermiticity error: 1.262e-13.
- Minimum sampled density eigenvalue: 0.000e+00.
- Maximum detailed-balance error: 4.441e-16.
- Maximum Gibbs-stationarity error: 2.166e-13.
- Overall validation: PASS.

## Interpretation limits

This model tests whether a thermally consistent population-relaxation bath can transfer excitation toward the low-energy PYR5 state. It does not provide a microscopic relaxation time because no bath spectral density has been derived from the current trajectory.
