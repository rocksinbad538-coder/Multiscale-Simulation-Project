# Day020 High-Dephasing and Zeno-Scaling Audit

## Model

The effective Haken-Strobl pair rate was evaluated as

\[
k_{ij} = \frac{2(J_{ij}/\hbar)^2\gamma_\phi}{\gamma_\phi^2 + (\Delta E_{ij}/\hbar)^2}.
\]

The rate is maximal when `gamma_phi = |Delta E|/hbar`.

## Validation

- Hamiltonians: 21/21.
- Maximum Hamiltonian symmetry error: 0.000e+00 eV.
- Maximum classical probability-norm error: 3.997e-15.
- Selected gamma values: 20, 50, 100, 200 ps^-1.

## Characteristic turnover scales

- Mean PYR2-PYR4 predicted turnover scale: 20.884 ps^-1.
- Mean turnover scale for pairs involving PYR5: 469.140 ps^-1.
- Mean high-manifold k(100)/k(200): 1.917808.
- Mean PYR5-pair k(100)/k(200): 0.566727.

## Interpretation

The PYR2-PYR4 pairs have detuning frequencies within or below the tested 20-50 ps^-1 turnover region. Their effective rates therefore decrease at sufficiently large gamma_phi, consistent with high-dephasing or Zeno suppression.

Pairs involving PYR5 have substantially larger detuning frequencies because PYR5 lies roughly 300 meV below the high-energy manifold. The tested range up to 200 ps^-1 does not necessarily exceed their pair-specific turnover scales. Their rates and PYR5 populations can therefore continue to increase while PYR2-PYR4 transfer is already being suppressed.

The PYR5 increase remains a phenomenological pure-dephasing result. It is not thermal downhill relaxation because the model has no detailed balance or population-relaxation bath.

## Classical-reduction comparison

| gamma_phi (ps^-1) | Start time (ps) | RMSE all | RMSE PYR2-PYR4 | RMSE PYR5 |
|---:|---:|---:|---:|---:|
| 20 | 0.2500 | 4.778686e-03 | 6.371262e-03 | 7.149916e-05 |
| 50 | 0.1000 | 1.514381e-03 | 2.017977e-03 | 6.902209e-05 |
| 100 | 0.0500 | 5.433656e-04 | 7.191808e-04 | 7.744768e-05 |
| 200 | 0.0500 | 2.189372e-04 | 2.711989e-04 | 9.100258e-05 |
