# Day019 Haken-Strobl pure-dephasing sensitivity

## Protocol

- Primary Hamiltonian: four-state bright TDC-AC-corrected model.
- Twenty-one solvent snapshots were propagated independently.
- Lindblad operators: `L_i = sqrt(gamma_phi) |i><i|` for each local bright state.
- Therefore each off-diagonal density-matrix element decays directly at `gamma_phi` in the absence of the Hamiltonian.
- Gamma sweep: 0, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200 ps^-1.
- Propagation interval: 0-20.000 ps; output interval: 0.010 ps.
- No population-relaxation operators or detailed-balance constraints were included.

## Numerical validation

- Maximum trace error: 2.296e-13.
- Maximum Hermiticity error: 1.091e-13.
- Minimum sampled density-matrix eigenvalue: -2.801e-14.
- Gamma=0 maximum population error versus the accepted coherent trajectories: 5.778e-13.
- Overall numerical validation: PASS.

## Aggregate high-energy-manifold response

| gamma_phi (ps^-1) | T_phi (ps) | Mean minimum survival | Mean time-averaged survival | Mean maximum population on another PYR2-PYR4 site | Maximum PYR5 population | Integrated l1 coherence |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | inf | 0.878506 | 0.920923 | 0.121480 | 8.200093e-05 | 9.943157 |
| 0.05 | 20 | 0.867697 | 0.902208 | 0.132275 | 8.683199e-05 | 8.014375 |
| 0.1 | 10 | 0.848203 | 0.885868 | 0.151752 | 1.211120e-04 | 7.014756 |
| 0.2 | 5 | 0.806033 | 0.858231 | 0.193892 | 1.940558e-04 | 5.921089 |
| 0.5 | 2 | 0.718431 | 0.799117 | 0.281402 | 4.030322e-04 | 4.507755 |
| 1 | 1 | 0.639616 | 0.738835 | 0.360067 | 7.223806e-04 | 3.523819 |
| 2 | 0.5 | 0.557209 | 0.670646 | 0.442171 | 1.300114e-03 | 2.664842 |
| 5 | 0.2 | 0.453022 | 0.577273 | 0.545455 | 2.797083e-03 | 1.728744 |
| 10 | 0.1 | 0.396392 | 0.518148 | 0.600588 | 5.000637e-03 | 1.203172 |
| 20 | 0.05 | 0.365498 | 0.486837 | 0.628527 | 9.266735e-03 | 0.860303 |
| 50 | 0.02 | 0.365377 | 0.514756 | 0.620170 | 2.363147e-02 | 0.623800 |
| 100 | 0.01 | 0.414702 | 0.594556 | 0.558350 | 5.071455e-02 | 0.504068 |
| 200 | 0.005 | 0.515815 | 0.700213 | 0.439289 | 9.813771e-02 | 0.373591 |

## Sensitivity extrema

- Gamma=0 mean maximum high-manifold transfer: 0.121480.
- Gamma maximizing mean high-manifold transfer: 20 ps^-1, with value 0.628527.
- Gamma minimizing mean time-averaged survival: 20 ps^-1, with value 0.486837.
- Largest PYR5 population anywhere in the phenomenological gamma sweep: 9.813771e-02 at gamma=200 ps^-1.

## Interpretation boundary

This calculation is a phenomenological robustness analysis. The gamma values are not extracted from the 21 solvent snapshots and are not microscopic dephasing rates. Pure dephasing can convert coherent coupling into population redistribution, but it contains no energy-selective downhill relaxation and does not enforce thermal detailed balance. Any increase of PYR5 population must therefore not be interpreted as physical bath-assisted relaxation. A subsequent relaxation model requires explicitly declared rates or a separately justified bath spectral density.
