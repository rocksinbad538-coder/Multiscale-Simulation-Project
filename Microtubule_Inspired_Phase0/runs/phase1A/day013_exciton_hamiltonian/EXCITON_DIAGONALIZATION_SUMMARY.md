# Day013 Exciton Hamiltonian Diagonalization

Basis:

PYR2, PYR3, PYR4, PYR5

Input Hamiltonian:

- `exciton_hamiltonian_4x4_eV.csv`

Outputs:

- `exciton_eigenstates_day013.csv`
- `exciton_eigenvectors_day013.csv`

## Exciton energies

| State | Energy (eV) | Shift from mean site energy (meV) | Participation ratio | Dominant site |
|---|---:|---:|---:|---|
| X1 | 3.756985 | -18.515 | 2.658 | PYR5 |
| X2 | 3.770327 | -5.173 | 1.909 | PYR3 |
| X3 | 3.777730 | 2.230 | 1.776 | PYR2 |
| X4 | 3.796958 | 21.458 | 3.255 | PYR4 |

## Interpretation

The participation ratio estimates how many pyrene sites contribute appreciably to each excitonic eigenstate.

- PR ≈ 1: localized mostly on one chromophore.
- PR ≈ 2: delocalized over approximately two chromophores.
- PR ≈ 4: delocalized over the full four-site network.

These results are based on a preliminary splitting-derived Hamiltonian and should be refined after transition-dipole or transition-density coupling analysis.
