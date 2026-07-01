# Day020 Combined Dephasing and Relaxation Sensitivity

## Scope

This calculation combines coherent dynamics, local site-basis pure dephasing, and global energy-basis detailed-balance relaxation.

The dephasing and relaxation terms represent separate phenomenological environments. For gamma_phi > 0, the combined stationary state is not assumed to equal the Gibbs state and is solved explicitly.

## Parameter sweep

- Temperatures: 150, 300 K.
- Kappa values: 0.1, 1, 10 ps^-1.
- Gamma values: 0, 1, 20, 100 ps^-1.
- Propagation interval: 0-100 ps.

## Numerical validation

- Maximum trace error: 1.830e-12.
- Maximum Hermiticity error: 1.151e-13.
- Minimum sampled density eigenvalue: 0.000e+00.
- Maximum combined stationary residual: 5.460e-13.
- Gamma=0 reference error: 0.000e+00.
- Overall validation: PASS.

## Aggregate diagnostic ranges

- PYR4 thermal gateway fraction: 0.842996 to 0.914134.
- Mean stationary PYR5 population: 0.241117 to 0.999977.

## Interpretation limit

Absolute relaxation times remain proportional to the phenomenological kappa scale. Local pure dephasing is not a thermally balanced population bath, so deviations of the combined stationary state from Gibbs are expected and must be interpreted as model sensitivity rather than a microscopic equilibrium prediction.
## Stationary-state solver correction

For `gamma_phi = 0`, the stationary density is now assigned directly to the analytical Gibbs state. This removes the small loss of positivity produced by the ill-conditioned complex least-squares nullspace solve for nearly pure Gibbs states. No propagated trajectories were modified.
