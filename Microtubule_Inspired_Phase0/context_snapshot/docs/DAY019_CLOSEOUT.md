# Day019 Closeout — Bright-State Couplings and Dynamics

## Scope completed

Day019 established, validated, and dynamically tested the four-state bright excitonic model derived from the accepted frozen-solute MD ensemble and embedded TDDFT calculations.

Completed components:

1. Transition-density cube grid convergence.
2. Bright-state transition-density production at 80 x 80 x 80.
3. ORCA transition-dipole normalization validation.
4. Transition-density-derived atom-centered charge couplings.
5. Point-dipole versus finite-size coupling comparison.
6. Finite-size-corrected bright Hamiltonians for 21 solvent snapshots.
7. Model-readiness and temporal-resolution audit.
8. Quasi-static coherent ensemble dynamics.
9. Haken-Strobl pure-dephasing sensitivity analysis.

## Transition-density normalization

The raw ORCA transition-density cubes reproduced approximately 1/sqrt(2) of the independently printed ORCA transition dipole.

After sqrt(2) normalization:

- Four bright states passed.
- Production grid: 80 x 80 x 80.
- Maximum dipole-magnitude error: 0.1540%.
- Minimum directional cosine: 0.99999666.
- Maximum absolute transition-charge residual: below 1.0e-4 e.
- Boundary densities were negligible.

Accepted bright roots:

- PYR2: S2
- PYR3: S2
- PYR4: S2
- PYR5: S1

## Finite-size coupling benchmark

Transition-density-derived atom-centered charge couplings were constructed using a nearest-atom Voronoi partition followed by a minimum-norm correction imposing:

- zero total transition charge;
- the exact ORCA transition dipole.

The largest absolute correction relative to the point-dipole model was:

- 0.066408 meV at frame000;
- 0.067962 meV over the full 21-frame ensemble.

All coupling signs were preserved.

The maximum relative difference was 7.21%, occurring for the smallest long-range coupling. The corrected and point-dipole Hamiltonians remained spectrally and structurally equivalent within the relevant energy scales.

## Corrected Hamiltonian sensitivity

Across the 21 snapshots:

- Maximum bright-state eigenvalue shift: 0.042074 meV.
- Minimum matched eigenvector overlap: 0.99985570.
- Maximum correction / maximum diagonal-energy SD: 0.004239.
- Maximum correction / minimum local S1-S2 gap: 0.001282.

The TDC-AC-corrected four-state bright Hamiltonian was accepted as the primary corrected control model.

## Model-readiness audit

The ensemble contains:

- 21 snapshots;
- 5 ps spacing;
- diagonal-energy SD range: 7.403-16.034 meV;
- maximum corrected bright coupling: 1.492201 meV;
- electronic coupling timescale hbar/max|J|: 0.441101 ps.

The snapshot interval is 11.335 times the electronic coupling timescale.

Accepted use:

- quasi-static disorder ensemble;
- independent coherent propagation;
- eigenstate and localization analysis;
- phenomenological dephasing sensitivity.

Not supported by the current sampling:

- sequential stochastic Hamiltonian propagation;
- microscopic bath autocorrelation;
- spectral-density extraction;
- memory-kernel extraction;
- microscopic pure-dephasing rates.

## Quasi-static coherent dynamics

Exact unitary propagation was performed for:

- point-dipole and TDC-AC models;
- four localized initial states;
- 21 independent solvent snapshots;
- 0-20 ps;
- 0.005 ps output interval.

Maximum norm error:

- 3.553e-15.

Corrected-model minimum ensemble survival:

- PYR2: 0.873483
- PYR3: 0.856748
- PYR4: 0.905243
- PYR5: 0.999914

Maximum ensemble-mean PYR5 population from PYR2-PYR4:

- 8.200093e-05.

The large PYR5 energy offset suppresses direct coherent transfer in the closed bright-state model.

## Haken-Strobl sensitivity

A phenomenological pure-dephasing sweep was performed for:

gamma_phi = 0, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, and 10 ps^-1.

Numerical validation:

- Maximum trace error: 2.005e-13.
- Maximum Hermiticity error: 1.091e-13.
- Minimum sampled density-matrix eigenvalue: -2.801e-14.
- Gamma=0 error versus accepted coherent trajectories: 5.778e-13.
- Overall status: PASS.

Within the tested range:

- high-manifold population redistribution increased monotonically with gamma_phi;
- the Zeno turnover was not yet identified;
- gamma_phi = 10 ps^-1 produced the largest tested redistribution;
- maximum PYR5 population reached 5.000637e-03.

The PYR5 population generated under pure dephasing must not be interpreted as physical downhill relaxation. The Haken-Strobl model contains no detailed balance or energy-selective relaxation.

## Accepted models

Primary corrected model:

- Four-state bright TDC-AC-corrected Hamiltonian ensemble.

Controls:

- Four-state bright point-dipole ensemble.
- Full eight-state hybrid Hamiltonian for sensitivity only.

The alternate-state and mixed-family signed couplings remain insufficiently benchmarked for primary dynamical use.

## Next technical steps

1. Extend the pure-dephasing sweep beyond 10 ps^-1 to locate the Zeno turnover.
2. Quantify transfer-rate scaling in the high-dephasing regime.
3. Decide whether a phenomenological energy-relaxation model is required.
4. If relaxation is introduced, enforce explicit detailed balance and clearly declare all assumed rates.
5. Generate final figures comparing:
   - coherent ensemble dynamics;
   - point-dipole versus TDC-AC dynamics;
   - dephasing sensitivity;
   - PYR5 accessibility.
6. Determine whether a new, finer-resolution embedded trajectory is required for microscopic bath modeling.

## Interpretation boundary

The current results support a quasi-static bright-state ensemble and phenomenological open-system sensitivity analysis. They do not provide microscopic solvent dephasing rates, a spectral density, or a physically parameterized energy-relaxation model.
