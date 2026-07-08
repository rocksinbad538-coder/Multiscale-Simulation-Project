# R1 NVT Checkpoint Continuation to 50 ps

## Purpose

The 20 ps R1 frozen-solute trajectory showed a short initial
thermalization transient followed by stable behavior:

- 5–20 ps temperature:
  **300.1487
  ± 1.3613 K**
- 10–20 ps temperature:
  **300.1393
  ± 1.3364 K**
- Endpoint lumen occupancy:
  **428 waters**
- Endpoint initially luminal waters retained:
  **428 waters**

A 30 ps checkpoint continuation was therefore prepared to extend the
same trajectory from 20 to 50 ps.

## Continuation method

The original TPR was extended using `gmx convert-tpr`.

No new `grompp` operation was performed. No velocities were generated.

The future run must use:

- the extended TPR;
- the exact 20 ps checkpoint;
- `mdrun -cpi`;
- `-noappend`, preserving the original 20 ps result.

This retains the checkpoint coordinates, velocities, thermostat state,
random state, and current integration step.

## Source state

- Source atoms:
  **68314**
- Source steps:
  **40000**
- Source nominal duration:
  **20.000000 ps**
- Checkpoint step:
  **40000**
- Checkpoint time:
  **20.000000 ps**
- Final GRO atoms:
  **68314**
- Checkpoint SHA256:
  `5d105645d6037544147b4f3546dcad77c340d206d14a0b8c543fd484aef982cb`

## Extended TPR

- Extended atoms:
  **68314**
- Extended total steps:
  **100000**
- Extended target time:
  **50.000000 ps**
- Remaining steps:
  **60000**
- Remaining duration:
  **30.000000 ps**
- Compressed trajectory stride:
  **1000 steps**
- Expected continuation frames:
  **61**
- TPR differences beyond `nsteps`:
  **0**

## Gates

- `convert_tpr_return_code_zero`: **PASS**
- `source_TPR_dump_return_code_zero`: **PASS**
- `extended_TPR_dump_return_code_zero`: **PASS**
- `checkpoint_dump_return_code_zero`: **PASS**
- `source_TPR_has_68314_atoms`: **PASS**
- `extended_TPR_has_68314_atoms`: **PASS**
- `source_final_GRO_has_68314_atoms`: **PASS**
- `source_TPR_has_40000_steps`: **PASS**
- `extended_TPR_has_100000_steps`: **PASS**
- `source_dt_is_0p0005ps`: **PASS**
- `extended_dt_is_0p0005ps`: **PASS**
- `source_total_time_is_20ps`: **PASS**
- `extended_total_time_is_50ps`: **PASS**
- `checkpoint_step_is_40000`: **PASS**
- `checkpoint_time_is_20ps`: **PASS**
- `remaining_steps_are_60000`: **PASS**
- `remaining_time_is_30ps`: **PASS**
- `trajectory_stride_is_preserved`: **PASS**
- `expected_continuation_frames_are_61`: **PASS**
- `no_TPR_changes_beyond_nsteps`: **PASS**
- `checkpoint_copy_is_bitwise_identical`: **PASS**
- `thermal_review_authorized_continuation`: **PASS**
- `no_unreviewed_source_gate_failures`: **PASS**

## Decision

- Decision:
  **R1_CHECKPOINT_CONTINUATION_TO_50PS_PREPARED**
- Failed gates:
  **NONE**
- Checkpoint-continuation execution authorized:
  **YES**
- Velocity regeneration authorized:
  **NO**
- New independent trajectory authorized:
  **NO**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `RUN_R1_20_TO_50PS_CHECKPOINT_CONTINUATION`

The continuation remains part of the R1 frozen steric positive-control
screening. It does not establish chemical realizability.
