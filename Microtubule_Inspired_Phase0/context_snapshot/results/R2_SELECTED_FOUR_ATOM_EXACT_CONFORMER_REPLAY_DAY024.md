# R2 Selected Four-Atom BN Bridge Exact Conformer Replay

## Scope

This gate replays the deterministic Gate 3L selection using the original
source code, random seed, conformer library and global coordinate-descent
selection.

## Recovered selection

- Random seed: **35004**
- Selected conformers: **30**
- Lower/upper conformers:
  **15/15**
- Internal bridge coordinates:
  **120**
- Lower/upper internal coordinates:
  **60/60**

## Recovered geometry

- Minimum angle:
  **70.899839788 degrees**
- Minimum local clearance:
  **0.159747397 nm**
- Maximum bond deviation:
  **4.410696146e-05 nm**
- Maximum library distance error:
  **2.318709098e-04 nm**
- Local clashes:
  **0**
- Local angle violations:
  **0**

## Reproducibility

- Selected mappings match original Gate 3L:
  **True**
- Differences:
  **NONE**

## Gates

- `Gate3L_source_script_checksum_matches`: **PASS**
- `Gate3L_source_decision_is_accepted`: **PASS**
- `instrumented_replay_decision_is_accepted`: **PASS**
- `selected_mapping_rows_match_original_Gate3L`: **PASS**
- `30_selected_conformers_were_recovered`: **PASS**
- `15_selected_conformers_were_recovered_per_end`: **PASS**
- `120_internal_bridge_coordinates_were_recovered`: **PASS**
- `60_internal_bridge_coordinates_were_recovered_per_end`: **PASS**
- `all_bridge_path_identifiers_are_unique`: **PASS**
- `coordinate_paths_match_selected_conformer_paths`: **PASS**
- `all_bridge_node_identifiers_are_unique`: **PASS**
- `all_selected_candidates_have_zero_local_clashes`: **PASS**
- `all_selected_candidates_have_zero_angle_violations`: **PASS**
- `minimum_selected_angle_is_at_least70deg`: **PASS**
- `minimum_selected_local_clearance_is_at_least0p120nm`: **PASS**
- `maximum_selected_bond_deviation_is_at_most0p003nm`: **PASS**
- `library_seed_is_35004_for_every_selected_conformer`: **PASS**
- `coordinates_are_marked_as_exact_Gate3L_selection`: **PASS**

## Decision

- Decision:
  **R2_SELECTED_FOUR_ATOM_BN_BRIDGE_EXACT_CONFORMERS_RECOVERED**
- Failed gates:
  **NONE**
- Coordinates applied to complete structure:
  **NO**
- Hydrogen coordinates generated:
  **NO**
- Molecular topology generation authorized:
  **NO**
- Energy minimization authorized:
  **NO**
- MD authorized:
  **NO**
- QM authorized:
  **NO**
- Required next step:
  `BUILD_AND_VALIDATE_R2_SELECTED_FOUR_ATOM_BN_BRIDGE_STATIC_COORDINATE_EMBEDDING`
