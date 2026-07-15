# Day021 Original Frozen-Topology Identity Verification

## Candidate

- TOP: `parameters/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032/hbn_pyrene_4_hydratable_gap45_pyr5shift_clean032.top`
- HBN ITP: `parameters/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032/hbn_fixed_dummy.itp`
- Explicit reconstructed `ld-seed`: 2099210158

## Exact TPR-section comparison

- Atom count: 68320 versus 68320.
- `inputrec` exact equality: True.
- Differing `inputrec` parameters: 0.
- `topology` exact equality: True.
- `box` exact equality: True.
- First coordinate difference: atom index 1784.
- Frozen HBN+PYR coordinate prefix: True.

## Interpretation

The unbonded HBN topology stored under the accepted parameter set is the exact topology embedded in the accepted frozen-solute TPR. The first coordinate difference occurs at atom 1784, the first water atom, because the reconstruction uses the final 100 ps GRO while the accepted TPR stores the initial production state.

## Decision

- Exact original-topology identity: PASS.
- Provenance validation: PASS.
- Mobile bonded-HBN topology remains a separate model and must be validated independently before mobile production.
