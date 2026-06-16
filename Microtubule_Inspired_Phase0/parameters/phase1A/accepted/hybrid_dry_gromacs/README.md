# Accepted dry h-BN + four-pyrene GROMACS baseline

## Status

Accepted as Phase 1A dry hybrid GROMACS baseline.

System:

- h-BN scaffold: fixed/dummy/non-interacting placeholder
- pyrene molecules: 4
- pyrene force field: GAFF2
- pyrene charges: AM1-BCC neutralized
- solvent: absent
- water: absent

## Accepted counts

- h-BN atoms: 1680
- pyrene molecules: 4
- pyrene atoms total: 104
- total atoms: 1784

## Accepted topology gate

- GROMACS TPR generated: true
- fatal/error detected: false
- warnings: 0
- notes: 2

The two notes are expected:

1. h-BN has no internal bonded potentials because it is intentionally used as a fixed/dummy scaffold at this stage.
2. PME mesh performance note.

## Accepted minimization

- minimization converged: true
- fatal/NaN detected: false
- final potential energy: 184.69263 kJ/mol
- final maximum force: 47.837498 kJ/mol/nm

## Use constraints

Accepted labels:

- dry h-BN + four-pyrene GROMACS assembly
- preliminary GROMACS-compatible hybrid baseline
- topology/minimization sanity-tested dry system

Rejected labels:

- physical h-BN force-field model
- hydrated system
- production MD-ready system
- QM-ready full hybrid
- TDDFT-ready full hybrid
- validated excitonic Hamiltonian
