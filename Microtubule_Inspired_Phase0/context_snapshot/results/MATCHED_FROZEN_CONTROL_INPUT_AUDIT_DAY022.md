# Matched Frozen-Control Input Audit

## Stage02 authoritative initial state

- TPR: `runs/phase1A/day021_mobile_restraint_protocol/execution/02_nvt_k10000_1ps/02_nvt_k10000_1ps.tpr`
- Atoms: **68320**
- Initial coordinate SHA256: `7c5c78f35a458ab26170df2b2eb7f618cde99e98763800c58e7c62441de90e86`
- Initial water-velocity SHA256: `d11ac84fa2d4a2a9a594a91ac0a6e0714dc3caa9da33f0ed12617371feff722d`

## Initial velocity populations

- HBN nonzero fraction: 1.000000
- PYR nonzero fraction: 1.000000
- Water nonzero fraction: 0.750000
- Water speed mean/std: 1.41348138/1.35258344 nm ps^-1

## Closest GRO representation

- File: `runs/phase1A/day021_mobile_restraint_protocol/execution/01_em_k10000/01_em_k10000.gro`
- All-atom RMS/max difference from TPR coordinates: 0.00000000/0.00000000 nm

## Stage02 MDP settings

- `integrator`: `md`
- `dt`: `0.00025`
- `nsteps`: `4000`
- `continuation`: `no`
- `gen-vel`: `yes`
- `gen-seed`: `20260706`
- `tcoupl`: `V-rescale`
- `tc-grps`: `System`
- `tau-t`: `1.0`
- `ref-t`: `300`
- `freezegrps`: `MISSING`
- `freezedim`: `MISSING`
- `comm-mode`: `None`
- `comm-grps`: `MISSING`
- `nstcomm`: `MISSING`
- `constraints`: `none`
- `constraint-algorithm`: `MISSING`
- `pbc`: `xyz`
- `nstxout-compressed`: `400`

## Index groups

- `System`, `HBN`, `PYR`, `HBN_PYR`, `SOL`

## Decision

- Matched frozen-control execution authorized: **NO**
- Next requirement: construct the frozen-control TPR from the authoritative Stage02 initial state and verify that its water-velocity SHA256 is identical before running MD.
