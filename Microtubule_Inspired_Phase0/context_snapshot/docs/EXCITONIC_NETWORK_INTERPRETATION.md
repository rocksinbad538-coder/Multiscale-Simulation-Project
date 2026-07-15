# Day013 Excitonic Network Interpretation

## System

Phase 1A hydrated h-BN scaffold with four pyrene chromophores.

Input electronic data:

- ORCA 6.1.1
- PBE0/def2-SVP
- TDDFT
- Isolated pyrene monomers
- All six pyrene dimers

## Site Energies

| Site | S1 energy (eV) |
|---|---:|
| PYR2 | 3.779 |
| PYR3 | 3.774 |
| PYR4 | 3.782 |
| PYR5 | 3.767 |

Average site energy:

- 3.7755 eV

Maximum S1 spread:

- 0.015 eV

Interpretation:

The pyrene chromophores form a nearly homogeneous site-energy landscape.

## Effective Couplings

Preliminary couplings were estimated from dimer excited-state splittings:

J_eff = (E2 - E1) / 2

| Pair | J_eff (meV) | Geometry class |
|---|---:|---|
| PYR4-PYR5 | 11.0 | near-orthogonal |
| PYR3-PYR4 | 10.5 | near-orthogonal |
| PYR2-PYR5 | 8.5 | near-orthogonal |
| PYR2-PYR3 | 6.5 | near-orthogonal |
| PYR2-PYR4 | 2.5 | near-parallel long-range |
| PYR3-PYR5 | 1.0 | near-parallel long-range |

## Main Physical Interpretation

The strongest preliminary couplings are not found in the long-range near-parallel pairs. Instead, they occur mostly among closer near-neighbor pairs.

This suggests that the current Phase 1A excitonic network is dominated by local chromophore connectivity rather than by same-orientation long-range pairing.

The system is not a close pi-stacked aggregate. Minimum pyrene-pyrene atom distances are approximately 22–42 Å, so the coupling regime is expected to be weak and primarily Coulombic/dipolar rather than exchange-dominated.

## Preliminary Hamiltonian

The first 4x4 excitonic Hamiltonian has been generated in:

- `exciton_hamiltonian_4x4_eV.csv`
- `exciton_hamiltonian_4x4_meV_relative.csv`

Diagonal entries are isolated monomer S1 site energies.

Off-diagonal entries are first-pass splitting-derived couplings.

## Limitations

The current coupling estimates are preliminary.

They should be refined using at least one of the following:

1. transition-dipole coupling analysis;
2. transition-density cube analysis;
3. fragment excitation difference analysis;
4. electrostatic embedding with local water/scaffold environment;
5. repeated calculations on additional MD frames.

## Next Technical Step

Diagonalize the preliminary 4x4 Hamiltonian to obtain:

- exciton eigenenergies;
- eigenvectors;
- participation ratios;
- localization/delocalization character;
- dominant chromophore contributions to each exciton state.

## Hamiltonian Diagonalization Result

The preliminary 4-site Hamiltonian was diagonalized.

| Exciton | Energy (eV) | Shift from mean site energy (meV) | Participation ratio | Dominant site |
|---|---:|---:|---:|---|
| X1 | 3.756985 | -18.515 | 2.658 | PYR5 |
| X2 | 3.770327 | -5.173 | 1.909 | PYR3 |
| X3 | 3.777730 | 2.230 | 1.776 | PYR2 |
| X4 | 3.796958 | 21.458 | 3.255 | PYR4 |

Interpretation:

The preliminary excitonic eigenstates are not fully localized on single chromophores. Participation ratios range from approximately 1.8 to 3.3, indicating partial delocalization across the four-pyrene network.

The lowest exciton, X1, is dominated by PYR5 but includes non-negligible contributions from PYR2, PYR3, and PYR4. The highest exciton, X4, has the largest participation ratio and is the most delocalized state in this preliminary Hamiltonian.

Scientific implication:

The Phase 1A chromophore network supports weak but non-negligible collective excitonic states. The current model suggests local-network-mediated excitonic delocalization rather than isolated chromophore behavior.

Important limitation:

These eigenstates depend on first-pass splitting-derived couplings. The next refinement should use transition dipole moments, oscillator strengths, and/or transition-density-based couplings.
