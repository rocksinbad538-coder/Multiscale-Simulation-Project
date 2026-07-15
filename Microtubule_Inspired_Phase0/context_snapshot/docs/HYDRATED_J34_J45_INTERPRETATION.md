# Hydrated Hamiltonian Sensitivity: J34 = 17 meV, J45 = 33 meV

## Input

This model updates the 4-site PYR2-PYR5 excitonic Hamiltonian using explicit-water dimer-derived effective couplings:

- PYR3-PYR4: J34 ≈ 17 meV
- PYR4-PYR5: J45 ≈ 33 meV

Site energies used:

- PYR2 = 3.755 eV
- PYR3 = 3.765 eV
- PYR4 = 3.778 eV
- PYR5 = 3.812 eV

## Result

The updated hydrated Hamiltonian produces strong mixing across the PYR3-PYR4-PYR5 block. PYR2 remains isolated in this specific model because its couplings were set to zero pending explicit hydrated-dimer calculations involving PYR2.

## Key interpretation

The hydrated coupling pattern is no longer a weak perturbation of the vacuum Hamiltonian. Explicit hydration increases the effective PYR4-PYR5 coupling and introduces substantial redistribution of excitonic weights across PYR3, PYR4, and PYR5.

The lowest excitonic state X1 is distributed mainly over PYR3/PYR4/PYR5, while the highest state X4 is dominated by PYR5 with significant PYR4 contribution. This supports the emerging conclusion that local hydration modifies both the state ordering and the delocalization structure of the chromophore aggregate.

## Limitation

This is still a partial hydrated Hamiltonian because only J34 and J45 have been updated from explicit hydrated dimer TDDFT diagnostics. The remaining pair couplings should be refined with additional explicit-water dimer calculations or a transition-density coupling protocol.
