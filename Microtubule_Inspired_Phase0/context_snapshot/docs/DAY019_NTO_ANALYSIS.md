# Day019 NTO analysis and model-space assessment

## Validation

- Representative jobs analyzed: 8/8
- NTO states analyzed: 16/16
- Tracked root equals brightest S1/S2 root: 8/8 jobs
- Every output passed normal termination, SCF, and TDDFT checks.

## State-resolved metrics

| Type | Frame | Site | Root | Tracked | Energy (eV) | fosc | n1 | n2 | n1+n2 | PR | Character |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| embedded | 3 | PYR5 | S1 | True | 3.750 | 0.615413 | 0.875815 | 0.093929 | 0.969744 | 1.2874 | single_pair_dominated |
| embedded | 3 | PYR5 | S2 | False | 3.852 | 0.042652 | 0.630689 | 0.338557 | 0.969245 | 1.9478 | two_pair_mixed |
| embedded | 5 | PYR2 | S1 | False | 4.006 | 0.058606 | 0.630607 | 0.341388 | 0.971995 | 1.9409 | two_pair_mixed |
| embedded | 5 | PYR2 | S2 | True | 4.068 | 0.424630 | 0.837745 | 0.138190 | 0.975935 | 1.3858 | single_pair_dominated |
| embedded | 5 | PYR5 | S1 | True | 3.797 | 0.647308 | 0.877578 | 0.091586 | 0.969164 | 1.2829 | single_pair_dominated |
| embedded | 5 | PYR5 | S2 | False | 3.871 | 0.020591 | 0.613507 | 0.355037 | 0.968544 | 1.9858 | two_pair_mixed |
| embedded | 13 | PYR5 | S1 | True | 3.783 | 0.661990 | 0.879568 | 0.089522 | 0.969090 | 1.2778 | single_pair_dominated |
| embedded | 13 | PYR5 | S2 | False | 3.836 | 0.010825 | 0.611480 | 0.357359 | 0.968839 | 1.9887 | two_pair_mixed |
| vacuum_reference | 0 | PYR2 | S1 | False | 4.001 | 0.000794 | 0.523468 | 0.447566 | 0.971034 | 2.1034 | two_pair_mixed |
| vacuum_reference | 0 | PYR2 | S2 | True | 4.090 | 0.505315 | 0.846859 | 0.127692 | 0.974551 | 1.3618 | single_pair_dominated |
| vacuum_reference | 0 | PYR3 | S1 | False | 4.001 | 0.000675 | 0.521070 | 0.449910 | 0.970980 | 2.1052 | two_pair_mixed |
| vacuum_reference | 0 | PYR3 | S2 | True | 4.081 | 0.514411 | 0.849181 | 0.125225 | 0.974405 | 1.3556 | single_pair_dominated |
| vacuum_reference | 0 | PYR4 | S1 | False | 4.015 | 0.001167 | 0.526102 | 0.444762 | 0.970864 | 2.1022 | two_pair_mixed |
| vacuum_reference | 0 | PYR4 | S2 | True | 4.089 | 0.518148 | 0.850711 | 0.123702 | 0.974413 | 1.3516 | single_pair_dominated |
| vacuum_reference | 0 | PYR5 | S1 | True | 3.780 | 0.664486 | 0.879565 | 0.089235 | 0.968800 | 1.2778 | single_pair_dominated |
| vacuum_reference | 0 | PYR5 | S2 | False | 3.846 | 0.005932 | 0.577559 | 0.390787 | 0.968345 | 2.0511 | two_pair_mixed |

## Aggregate findings

- Representative S1-S2 gaps span 53.0-102.0 meV.
- Tracked roots are single-pair dominated in 100.0% of cases.
- Alternate roots are two-pair mixed in 100.0% of cases.
- Mean tracked-root dominant occupation: 0.862128 (SD 0.016431).
- Mean alternate-root dominant occupation: 0.579310 (SD 0.045853).

## Physical interpretation

For PYR2-PYR4 vacuum references, S2 is the bright tracked root and is dominated by the 52a->53a NTO pair, while S1 contains two comparably weighted NTO pairs. For PYR5, the ordering is reversed: S1 is the bright single-pair-dominated state and S2 is the two-pair-mixed state. Representative embedded cases preserve this distinction while changing the degree of mixing.

## Preliminary model-space decision

**Primary recommendation: 8-state.**

Use the 8-state manifold as the primary low-energy excitonic model because both S1 and S2 remain within 120 meV and have systematically distinct NTO composition. Retain the 4-state tracked-bright manifold as a reduced control model.

This decision uses energies, oscillator strengths, and NTO occupation spectra. Orbital-shape visualization remains the next validation step because matching orbital indices do not by themselves prove spatial equivalence across sites and environments.
