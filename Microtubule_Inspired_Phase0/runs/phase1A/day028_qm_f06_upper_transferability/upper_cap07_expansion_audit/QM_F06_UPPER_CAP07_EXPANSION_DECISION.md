# QM_F06 UPPER HCAP07 Expansion Decision — Day028

## Correction to the automated inventory

The initial expansion summary counted three heavy-atom additions and zero
real hydrogen additions.

The boundary-edge table demonstrates that restored atom
`A:UPPER:14:2` is covalently bonded to the existing real R2 passivant:

`H4:UPPER:0046:0`

through edge:

`E:3113`

This real hydrogen must be included in the expanded QM fragment.

## Selected UPPER Boundary V2 construction

Remove:

- `HCAP:UPPER:07`

Add real R2 atoms:

- `A:UPPER:11:3`
- `A:UPPER:13:3`
- `A:UPPER:14:2`
- `H4:UPPER:0046:0`

Add artificial peripheral caps along original cut-bond vectors:

- `E:2897`: `A:UPPER:11:3 — A:UPPER:10:2`
- `E:2898`: `A:UPPER:11:3 — A:UPPER:10:4`
- `E:2913`: `A:UPPER:14:2 — A:UPPER:13:1`

Retain all other atoms and caps from the repaired 22-atom UPPER
fragment.

## Expected composition

- Initial atoms: 22
- Removed artificial caps: 1
- Added real heavy atoms: 3
- Added real R2 hydrogens: 1
- Added new artificial caps: 3
- Expected final atoms: 28

## Authorization state

- UPPER Boundary V2 construction: **AUTHORIZED**
- Pre-QM audit: **REQUIRED**
- ORCA input preparation: **NOT AUTHORIZED**
- ORCA execution: **NOT AUTHORIZED**
