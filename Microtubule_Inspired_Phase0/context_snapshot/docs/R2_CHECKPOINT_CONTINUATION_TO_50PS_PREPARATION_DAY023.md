# R2 Checkpoint Continuation to 50 ps — Preparation

## Purpose

The R2 20 ps frozen-solute screen remains formally nonstationary.
The final trajectory windows show flattening, so an exact checkpoint
continuation is prepared to test plateau formation versus continued
depletion.

The completed 0–20 ps trajectory will not be repeated.

## Source state

- Source atoms:
  **68332**
- Source steps:
  **40000**
- Source time:
  **20.000000 ps**
- Source checkpoint step/time:
  **40000/20.000000 ps**
- Source failed gate:
  **second_half_occupancy_slope_is_acceptable**

## Extended TPR

- Extended atoms:
  **68332**
- Extended steps:
  **100000**
- Target time:
  **50.000000 ps**
- Remaining steps:
  **60000**
- Remaining time:
  **30.000000 ps**
- XTC stride:
  **1000**
- Expected continuation frames:
  **61**
- TPR differences beyond normalized `nsteps`:
  **8**

## Checkpoint integrity

- Coordinate entries:
  **204996**
- Velocity entries:
  **204996**
- Source checkpoint SHA256:
  `59dd34a9a4d073c77dff5d51e54b3962446a2a788ac954127e336df1b7e5ca4f`
- Copied checkpoint SHA256:
  `59dd34a9a4d073c77dff5d51e54b3962446a2a788ac954127e336df1b7e5ca4f`
- Bitwise-identical copy:
  **YES**

## Gates

- `source_preparation_is_valid`: **PASS**
- `source_run_is_nonstationary_review_state`: **PASS**
- `transient_audit_authorized_extension_preparation`: **PASS**
- `source_failed_only_stationarity_gate`: **PASS**
- `source_run_completed_without_instability`: **PASS**
- `convert_tpr_return_code_zero`: **PASS**
- `source_TPR_dump_return_code_zero`: **PASS**
- `extended_TPR_dump_return_code_zero`: **PASS**
- `checkpoint_dump_return_code_zero`: **PASS**
- `source_TPR_has_68332_atoms`: **PASS**
- `extended_TPR_has_68332_atoms`: **PASS**
- `source_TPR_nsteps_is_40000`: **PASS**
- `extended_TPR_nsteps_is_100000`: **PASS**
- `source_and_extended_dt_are_0p0005ps`: **PASS**
- `source_and_extended_integrators_are_md`: **PASS**
- `continuation_flag_is_preserved`: **PASS**
- `XTC_stride_is_preserved_at_1000`: **PASS**
- `TPRs_differ_only_by_nsteps`: **FAIL**
- `checkpoint_step_is_40000`: **PASS**
- `checkpoint_time_is_20ps`: **PASS**
- `checkpoint_coordinate_entries_are_complete`: **PASS**
- `checkpoint_velocity_entries_are_complete`: **PASS**
- `checkpoint_coordinate_and_velocity_counts_match`: **PASS**
- `checkpoint_copy_is_bitwise_identical`: **PASS**
- `remaining_steps_are_60000`: **PASS**
- `remaining_time_is_30ps`: **PASS**
- `expected_continuation_frames_are_61`: **PASS**

## Decision

- Decision:
  **R2_CHECKPOINT_CONTINUATION_PREPARATION_REQUIRES_REVIEW**
- Failed gates:
  **TPRs_differ_only_by_nsteps**
- Checkpoint-continuation execution authorized:
  **NO**
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
  `REVIEW_R2_CHECKPOINT_CONTINUATION_PREPARATION_FAILURES`
