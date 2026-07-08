# R2 TPR-Extension Difference Classification

## Finding

The source and extended `gmx dump` outputs contain one residual
difference after normalization of `nsteps`.

That difference is the first-line filename header:

- Source:
  `/Users/alejandro/projects/Multiscale-Simulation-Project/Microtubule_Inspired_Phase0/runs/phase1A/day023_confinement_design/15_r2_frozen_solute_nvt_20ps_preparation/r2_frozen_solute_nvt_20ps.tpr:`
- Extended:
  `/Users/alejandro/projects/Multiscale-Simulation-Project/Microtubule_Inspired_Phase0/runs/phase1A/day023_confinement_design/18_r2_nvt_50ps_checkpoint_continuation_preparation/r2_frozen_solute_nvt_to_50ps.tpr:`

This line identifies the file supplied to `gmx dump`; it is not part of
the simulation input record, coordinates, velocities, box state, force
field, thermostat state, or integration state.

## State comparison

- Coordinate mismatches:
  **0**
- Velocity mismatches:
  **0**
- Box-state mismatches:
  **0**
- Scalar-parameter mismatches:
  **0**
- Changed parameters:
  **NONE**
- Residual body mismatches after removing the filename header and
  normalizing `nsteps`:
  **0**

## Run-length change

- Source `nsteps`:
  **40000**
- Extended `nsteps`:
  **100000**
- Source atoms:
  **68332**
- Extended atoms:
  **68332**

## Gates

- `original_preparation_is_expected_review_state`: **PASS**
- `original_preparation_failed_exactly_one_gate`: **PASS**
- `audit_found_exactly_one_semantic_mismatch`: **PASS**
- `audit_found_zero_coordinate_mismatches`: **PASS**
- `audit_found_zero_velocity_mismatches`: **PASS**
- `audit_found_zero_box_state_mismatches`: **PASS**
- `audit_found_zero_scalar_parameter_mismatches`: **PASS**
- `audit_found_no_changed_parameters`: **PASS**
- `audit_mismatch_family_is_only_other`: **PASS**
- `mismatch_record_is_only_dump_filename_header`: **PASS**
- `source_first_line_is_dump_filename_header`: **PASS**
- `extended_first_line_is_dump_filename_header`: **PASS**
- `normalized_TPR_bodies_are_exactly_equal`: **PASS**
- `source_nsteps_is_40000`: **PASS**
- `extended_nsteps_is_100000`: **PASS**
- `source_and_extended_atom_counts_are_68332`: **PASS**
- `checkpoint_copy_is_bitwise_identical`: **PASS**
- `checkpoint_step_is_40000`: **PASS**
- `checkpoint_time_is_20ps`: **PASS**
- `original_contract_did_not_authorize_execution`: **PASS**
- `original_contract_forbids_velocity_regeneration`: **PASS**
- `original_contract_forbids_source_rerun`: **PASS**

## Decision

- Decision:
  **R2_CHECKPOINT_CONTINUATION_TO_50PS_AUTHORIZED**
- Difference classification:
  **BENIGN_GMX_DUMP_FILENAME_HEADER**
- Physical TPR differences beyond `nsteps`:
  **0**
- Checkpoint-continuation execution authorized:
  **YES**
- Velocity regeneration authorized:
  **NO**
- Source 0–20 ps rerun authorized:
  **NO**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `RUN_R2_20_TO_50PS_CHECKPOINT_CONTINUATION`
