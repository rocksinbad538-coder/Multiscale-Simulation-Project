# Day013 Local-Water Site-Energy Interpretation

ORCA TDDFT calculations were completed for each pyrene chromophore with its local water shell within 0.50 nm.

## Main result

Local hydration red-shifts all four pyrene S1 energies.

| Site | Vacuum S1 (eV) | Water S1 (eV) | Shift (meV) |
|---|---:|---:|---:|
| PYR2 | 3.779 | 3.754 | -25 |
| PYR3 | 3.774 | 3.765 | -9 |
| PYR4 | 3.782 | 3.779 | -3 |
| PYR5 | 3.767 | 3.746 | -21 |

## Interpretation

The local water environment increases the site-energy spread from approximately 15 meV in vacuum to approximately 33 meV in the explicit local-water clusters.

This means hydration introduces additional energetic disorder, but does not completely destroy the near-homogeneity of the four-pyrene network.

The hydrated site-energy disorder is now comparable to, and somewhat larger than, the preliminary excitonic couplings of approximately 1–11 meV.

## Physical implication

The Phase 1A system remains compatible with weak excitonic delocalization, but the explicit water environment may partially localize the exciton states relative to the vacuum Hamiltonian.

## Next step

Build and diagonalize a water-shifted Hamiltonian using:

- hydrated site energies on the diagonal;
- current first-pass TDDFT splitting-derived couplings off-diagonal.

This will estimate whether explicit local hydration preserves or suppresses the partially delocalized exciton states.
