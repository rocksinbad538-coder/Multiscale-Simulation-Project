# Accepted local TIP4P/2005 water block

## Status

Accepted as local TIP4P/2005 water model block for Phase 1A hydration tests.

## Model

- Water model: TIP4P/2005
- Reference: Abascal & Vega, J. Chem. Phys. 123, 234505 (2005)
- Geometry:
  - O-H distance: 0.09572 nm
  - H-H distance: 0.15139 nm
  - M-site defined as virtual site
- Charges:
  - OW: 0.0000 e
  - HW1: +0.5564 e
  - HW2: +0.5564 e
  - MW: -1.1128 e
- Lennard-Jones site:
  - OW only
  - sigma = 0.315890 nm
  - epsilon = 0.774900 kJ/mol

## Accepted syntax test

- GROMACS version: 2025.4
- TPR generated: true
- fatal/error detected: false
- warnings: 0
- notes: 1

The single note is PME performance-related and does not invalidate the topology.

## Use constraints

Accepted label:

- local TIP4P/2005 GROMACS-compatible water block

Rejected labels:

- hydrated production system
- validated confined-water simulation
- final MD production setup
