# R2 Source-versus-Extended TPR Difference Audit

## Scope

This audit compares the complete `gmx dump` representations of the
20 ps source TPR and the 50 ps extended TPR.

The expected `nsteps` change is normalized before comparison. No TPR,
checkpoint, trajectory, or molecular-dynamics result is modified.

## Diff accounting

- Unified-diff total lines:
  **8**
- Actual removed records:
  **1**
- Actual added records:
  **1**
- Semantic mismatch records after normalizing `nsteps`:
  **1**

## State integrity

- Coordinate mismatches:
  **0**
- Velocity mismatches:
  **0**
- Box-state mismatches:
  **0**
- Scalar-parameter mismatches:
  **0**

## Changed parameters

- NONE

## Exact mismatch records

- Line 1 [other] ``
  - source: `/Users/alejandro/projects/Multiscale-Simulation-Project/Microtubule_Inspired_Phase0/runs/phase1A/day023_confinement_design/15_r2_frozen_solute_nvt_20ps_preparation/r2_frozen_solute_nvt_20ps.tpr:`
  - extended: `/Users/alejandro/projects/Multiscale-Simulation-Project/Microtubule_Inspired_Phase0/runs/phase1A/day023_confinement_design/18_r2_nvt_50ps_checkpoint_continuation_preparation/r2_frozen_solute_nvt_to_50ps.tpr:`

## Current authorization

- Checkpoint-continuation execution authorized:
  **NO**
- Required next step:
  `CLASSIFY_TPR_EXTENSION_RESIDUAL_DIFFERENCES`
