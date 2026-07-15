# Day021 Mobile Bonded-Topology Transition Audit

## Preserved model components

- Defaults equal: True.
- B0 atom type equal: True.
- N0 atom type equal: True.
- Mobile HBN LJ nonzero: True.
- HBN atom table equal: True.
- Pyrene model equal: True.
- TIP4P/2005 model equal: True.
- Molecular composition equal: True.

## Intended HBN changes

- Frozen HBN bonded counts: {'atoms': 1680, 'bonds': 0, 'angles': 0, 'dihedrals': 0, 'pairs': 0, 'constraints': 0, 'exclusions': 0}.
- Mobile HBN bonded counts: {'atoms': 1680, 'bonds': 2460, 'angles': 4860, 'dihedrals': 1620, 'pairs': 0, 'constraints': 0, 'exclusions': 0}.
- Expected bonded counts pass: True.

## GROMACS reconstruction

- `grompp`: PASS.
- Rebuilt atoms: 68320.
- Rebuilt interaction parameter types: 902.

## Decision

- Intended frozen-to-bonded model transition: PASS.
- This test validates model continuity and topology construction only. It does not yet authorize unrestrained production MD.
