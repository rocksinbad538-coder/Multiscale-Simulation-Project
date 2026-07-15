# Day019 NTO input audit

## Scope

Eight representative calculations were prepared for natural transition orbital analysis. Both S1 and S2 are requested for every case to support comparison of the two low-lying local states.

## NTO settings

- `DoNTO true`
- `NTOStates 1,2`
- `NTOThresh 1e-4`

## Validation

- Inputs generated: 8/8
- Embedded cases: 4
- Vacuum-reference cases: 4
- Embedded point-charge files copied and hash-verified: 4/4
- Unexpected source-input changes: 0
- States requested per case: S1 and S2

## Cases

| Type | Frame | Site | Tracked root | Reason | Target job |
|---|---:|---|---:|---|---|
| embedded | 3 | PYR5 | S1 | minimum_PYR5_solvent_shift | `frame003_PYR5_embedding_nto_s1_s2` |
| embedded | 5 | PYR2 | S2 | minimum_tracked_character_weight | `frame005_PYR2_embedding_nto_s1_s2` |
| embedded | 5 | PYR5 | S1 | maximum_PYR5_solvent_shift | `frame005_PYR5_embedding_nto_s1_s2` |
| embedded | 13 | PYR5 | S1 | minimum_S1_S2_separation_meV | `frame013_PYR5_embedding_nto_s1_s2` |
| vacuum_reference | 0 | PYR2 | S2 | vacuum_reference_for_each_chromophore | `frame000_PYR2_vacuum_reference_nto_s1_s2` |
| vacuum_reference | 0 | PYR3 | S2 | vacuum_reference_for_each_chromophore | `frame000_PYR3_vacuum_reference_nto_s1_s2` |
| vacuum_reference | 0 | PYR4 | S2 | vacuum_reference_for_each_chromophore | `frame000_PYR4_vacuum_reference_nto_s1_s2` |
| vacuum_reference | 0 | PYR5 | S1 | vacuum_reference_for_each_chromophore | `frame000_PYR5_vacuum_reference_nto_s1_s2` |
