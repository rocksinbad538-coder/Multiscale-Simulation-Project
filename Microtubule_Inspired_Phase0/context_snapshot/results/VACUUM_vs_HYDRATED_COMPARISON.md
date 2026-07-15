# Day013 Vacuum vs Hydrated Exciton Hamiltonian Comparison

## Purpose

Assess whether local explicit hydration preserves or suppresses excitonic delocalization in the Phase 1A four-pyrene network.

## Model

Two Hamiltonians were compared:

1. Vacuum Hamiltonian:
   - isolated pyrene TDDFT site energies;
   - first-pass splitting-derived dimer couplings.

2. Hydrated Hamiltonian:
   - pyrene + local water 0.50 nm TDDFT site energies;
   - same first-pass splitting-derived dimer couplings.

## Hydrated site energies

| Site | Vacuum S1 (eV) | Hydrated S1 (eV) | Shift (meV) |
|---|---:|---:|---:|
| PYR2 | 3.779 | 3.754 | -25.0 |
| PYR3 | 3.774 | 3.765 | -9.0 |
| PYR4 | 3.782 | 3.779 | -3.0 |
| PYR5 | 3.767 | 3.746 | -21.0 |

## Site-energy disorder

Vacuum site-energy range:

15.0 meV

Hydrated site-energy range:

33.0 meV

## Hydrated exciton eigenstates

| State | Energy (eV) | Shift from hydrated mean site energy (meV) | PR | Dominant site |
|---|---:|---:|---:|---|
| X1 | 3.738122 | -22.878 | 1.913 | PYR5 |
| X2 | 3.755444 | -5.556 | 2.468 | PYR2 |
| X3 | 3.761773 | 0.773 | 2.313 | PYR3 |
| X4 | 3.788661 | 27.661 | 1.903 | PYR4 |

## Vacuum vs hydrated participation

| State | Vacuum PR | Hydrated PR | Change |
|---|---:|---:|---:|
| X1 | 2.658 | 1.913 | -0.745 |
| X2 | 1.909 | 2.468 | 0.560 |
| X3 | 1.776 | 2.313 | 0.537 |
| X4 | 3.255 | 1.903 | -1.352 |

## Interpretation

Explicit local hydration red-shifts all pyrene S1 energies and increases energetic disorder relative to the isolated-chromophore Hamiltonian.

The key question is whether this disorder suppresses excitonic delocalization. The participation-ratio comparison quantifies this directly.

These results remain preliminary because the off-diagonal couplings are still splitting-derived and were not recalculated in the hydrated environment. The next refinement should compare splitting-derived couplings against transition-dipole or transition-density-based couplings.
