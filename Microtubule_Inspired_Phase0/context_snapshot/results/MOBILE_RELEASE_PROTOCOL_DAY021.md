# Day021 Mobile-Release Protocol

## Provenance

- Source topology: `parameters/phase1A/accepted/hybrid_hbnBonded_kang2000_improperGeo100_validated/hbn_pyrene_4_hydratable_gap45_pyr5shift_clean032_hbnBonded_kang2000_improperGeo100.top`
- Starting GRO: `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.gro`
- Accepted source files modified: no.

## Index groups

- System: 68320 atoms.
- HBN: 1680 atoms.
- PYR: 104 atoms.
- HBN_PYR: 1784 atoms.
- SOL: 66536 atoms.
- Index validation: PASS.

## Static GROMACS validation

- 00_em_k100000: grompp=PASS, position restraints=1706/1706.
- 01_em_k10000: grompp=PASS, position restraints=1706/1706.
- 02_nvt_k10000_1ps: grompp=PASS, position restraints=1706/1706.
- 03_nvt_k1000_2ps: grompp=PASS, position restraints=1706/1706.
- 04_nvt_k100_2ps: grompp=PASS, position restraints=1706/1706.
- 05_nvt_unrestrained_2ps: grompp=PASS, position restraints=0/0.

## Decision

- Protocol static validation: PASS.
- Scientific calculation started: no.
