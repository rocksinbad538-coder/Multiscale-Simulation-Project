# Day013 Preliminary Excitonic Hamiltonian

Hamiltonian basis:

PYR2, PYR3, PYR4, PYR5

Diagonal entries are isolated pyrene S1 site energies from ORCA TDDFT.

Off-diagonal entries are preliminary effective couplings estimated as:

J_eff = (E2 - E1) / 2

from the two lowest TDDFT dimer excited states.

Average site energy:

3.775500 eV

Files:

- exciton_hamiltonian_4x4_eV.csv
- exciton_hamiltonian_4x4_meV_relative.csv

Important limitation:

These are first-pass effective couplings from excited-state splittings. They are useful for a preliminary Hamiltonian, but should later be refined using transition-density, transition-dipole, or fragment-based coupling analysis.
