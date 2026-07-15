# R0 t=0 Hydrated Reference Audit

## Authoritative sources

- Accepted trajectory:
  `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.xtc`
- Accepted run input:
  `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.tpr`
- Accepted GRO:
  `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.gro`
- Validated hydration time series:
  `runs/phase1A/day021_mobile_restraint_protocol/execution/08_nvt_mobile_100ps/mobile_vs_frozen_water/mobile_frozen_water_timeseries.csv`

## Extracted R1 starting state

- Extracted GRO:
  `runs/phase1A/day023_confinement_design/01_r0_t0_reference/r0_accepted_t0_hydrated_system.gro`
- SHA256:
  `3e2f207361765c7448099591664f17bcecd4e2f53c516520d2ebcfb512028754`
- Parsed time:
  **0.0 ps**
- Atoms:
  **68320**

## Composition

- HBN atoms: **1680**
- PYR atoms: **104**
- Water atoms: **66536**
- TIP4P/2005 waters: **16634**
- Water-residue chunk consistency:
  **True**
- First water-site names:
  `OW:16634`

## Consistency with accepted R0

- Solute RMS difference:
  **0.000000000000 nm**
- Solute maximum difference:
  **0.000000000000 nm**
- Box maximum difference:
  **0.000000000000 nm**

The accepted GRO may represent the final frozen trajectory state.
Solute identity is still expected because HBN and PYR were frozen.

## Tube geometry

- PCA axis:
  **(0.00000000,
  0.00000000,
  1.00000000)**
- Robust axial planes:
  **-3.008714/3.009286 nm**
- p98 tube length:
  **6.018000 nm**
- Wall radius mean/median:
  **1.199126/1.199111 nm**
- Wall-radius q01/q99:
  **1.198667/1.199602 nm**
- Provisional accessible radius:
  **0.949111 nm**
- Provisional accessible volume:
  **17.030831 nm³**

## Authoritative hydration state

- t=0 lumen occupancy:
  **428 waters**
- Maximum accepted-trajectory occupancy:
  **437 waters**
- t=0 fraction of maximum:
  **0.979405**
- Accepted endpoint occupancy:
  **23 waters**
- t=0 lumen density:
  **15.744292 nm^-3**
- Independent provisional geometric count:
  **428 waters**

The validated time-series occupancy remains authoritative.
The independent geometric value is retained only as a cross-check.

## End-zone audit

### Lower end

- Ring atoms: **60**
- Water oxygens inside 0.30 nm:
  **22**
- Water oxygens outside 0.30 nm:
  **31**
- Minimum water-O/end-ring distance:
  **0.333138 nm**

### Upper end

- Ring atoms: **60**
- Water oxygens inside 0.30 nm:
  **22**
- Water oxygens outside 0.30 nm:
  **24**
- Minimum water-O/end-ring distance:
  **0.313480 nm**

## Gates

- `extracted_frame_has_68320_atoms`: **PASS**
- `ordered_HBN_count_is_1680`: **PASS**
- `ordered_PYR_count_is_104`: **PASS**
- `ordered_water_atom_count_is_66536`: **PASS**
- `TIP4P_water_count_is_16634`: **PASS**
- `water_chunks_are_residue_consistent`: **PASS**
- `solute_matches_accepted_frozen_GRO`: **PASS**
- `box_matches_accepted_frozen_GRO`: **PASS**
- `authoritative_t0_occupancy_is_positive`: **PASS**
- `authoritative_t0_occupancy_at_least_100`: **PASS**
- `authoritative_t0_is_at_least_90pct_of_max`: **PASS**
- `tube_length_is_positive`: **PASS**
- `wall_radius_is_positive`: **PASS**

## Decision

- Audit decision: **PASS**
- Authoritative R1 start accepted:
  **YES**
- Failed gates:
  **NONE**
- Required next step:
  `DEFINE_AND_GENERATE_R1_STERIC_CAP_PROTOTYPE`

No MD, topology modification, cap generation, or QM calculation was
performed by this audit.
