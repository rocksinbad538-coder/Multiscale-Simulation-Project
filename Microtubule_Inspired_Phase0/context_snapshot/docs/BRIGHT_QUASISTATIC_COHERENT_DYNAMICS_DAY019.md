# Day019 quasi-static coherent bright-state dynamics

## Protocol

- Twenty-one solvent snapshots were propagated independently as frozen Hamiltonian realizations.
- Exact unitary propagation was used; no numerical time integrator or interpolation between solvent frames was introduced.
- Propagation interval: 0-20.000 ps.
- Output interval: 0.005 ps.
- Initial conditions: one localized excitation on each of the four bright local states.
- Compared models: point-transition-dipole and TDC-AC finite-size-corrected bright Hamiltonians.

## Numerical validation

- Maximum population-norm error: 3.553e-15.
- Frames per model and initial condition: 21/21.
- Maximum absolute difference between ensemble-mean point and corrected populations: 5.124317e-02.
- Difference maximum context: initial PYR2_bright, target PYR2_bright, time 16.915 ps.

## Corrected-model ensemble summary

| Initial state | Minimum ensemble survival | Time (ps) | Mean diagonal-ensemble survival |
|---|---:|---:|---:|
| PYR2_bright | 0.873483 | 5.805 | 0.919967 |
| PYR3_bright | 0.856748 | 17.145 | 0.903037 |
| PYR4_bright | 0.905243 | 8.415 | 0.939477 |
| PYR5_bright | 0.999914 | 0.020 | 0.999955 |

## PYR5 coherent-accessibility audit

- Maximum ensemble-mean PYR5 population from a PYR2/PYR3/PYR4 localized initial excitation: 8.200093e-05.
- Maximum mean diagonal-ensemble PYR5 population from those three initial conditions: 4.179469e-05.
- Maximum PYR5 population in any individual snapshot from a high-energy initial state: 1.083852e-04 (initial PYR4_bright, frame 014).

The large PYR5 energy offset suppresses direct coherent transfer into PYR5 in the present closed-system model. This result does not exclude bath-assisted downhill relaxation, which is absent from unitary propagation.

## Interpretation boundary

Damping of ensemble-averaged oscillations in this analysis is inhomogeneous dephasing caused by averaging over static Hamiltonian realizations. It is not a microscopic decoherence rate. The calculation contains no population relaxation, pure-dephasing operator, spectral density, or sequential 5 ps solvent dynamics. Near-resonant PYR2-PYR3-PYR4 frames can support transient coherent mixing, whereas PYR5 remains energetically isolated in the closed bright-state model.
