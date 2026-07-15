# R2 20 ps Occupancy-Transient Audit

## Scope

This audit evaluates the completed R2 20 ps frozen-solute trajectory.
It does not repeat molecular dynamics and does not change the original
stationarity threshold.

The original 20 ps gate remains unpassed because the second-half
occupancy slope exceeded the predefined absolute limit of
0.50 waters/ps.

## Source result

- Source decision:
  **R2_FROZEN_SOLUTE_NVT_20PS_REQUIRES_REVIEW**
- Failed source gates:
  **second_half_occupancy_slope_is_acceptable**
- Initial/minimum/endpoint occupancy:
  **428/414/416**
- Endpoint occupancy fraction:
  **0.971963**
- Endpoint retained initial identities:
  **414/428**
- Endpoint identity-retention fraction:
  **0.967290**
- Endpoint noninitial luminal waters:
  **2**
- Maximum noninitial luminal waters:
  **4**

## Stationarity metric

- Second-half occupancy slope:
  **-0.509091 waters/ps**
- Original absolute limit:
  **0.500000 waters/ps**
- Excess beyond limit:
  **0.009091 waters/ps**

The slope is not reclassified as passing. It is used only to determine
whether a checkpoint extension is scientifically justified as a
stationarity diagnostic.

## Cumulative windows

- 0.0–20.0 ps: mean=420.8537; min/max=414/429; change=-12; slope=-0.767247 waters/ps
- 5.0–20.0 ps: mean=418.8065; min/max=414/426; change=-10; slope=-0.722581 waters/ps
- 10.0–20.0 ps: mean=416.7619; min/max=414/422; change=-5; slope=-0.509091 waters/ps
- 12.5–20.0 ps: mean=415.8125; min/max=414/417; change=-1; slope=-0.150000 waters/ps
- 15.0–20.0 ps: mean=415.6364; min/max=414/417; change=+0; slope=-0.163636 waters/ps

## Five-picosecond blocks

- 0.0–5.0 ps: mean=427.2000; change=-1; slope=-0.121212 waters/ps; maximum noninitial occupancy=4
- 5.0–10.0 ps: mean=423.1000; change=-4; slope=-0.836364 waters/ps; maximum noninitial occupancy=4
- 10.0–15.0 ps: mean=418.0000; change=-4; slope=-1.381818 waters/ps; maximum noninitial occupancy=2
- 15.0–20.0 ps: mean=415.6364; change=+0; slope=-0.163636 waters/ps; maximum noninitial occupancy=2

## Safety and execution gates

- `source_decision_requires_review`: **PASS**
- `exactly_one_source_gate_failed`: **PASS**
- `only_stationarity_slope_gate_failed`: **PASS**
- `all_other_source_gates_passed`: **PASS**
- `mdrun_completed_successfully`: **PASS**
- `trajectory_check_completed_successfully`: **PASS**
- `checkpoint_was_written`: **PASS**
- `no_instability_signatures`: **PASS**
- `trajectory_contains_41_frames`: **PASS**
- `initial_occupancy_is_428`: **PASS**
- `minimum_occupancy_remains_above_80_percent`: **PASS**
- `endpoint_occupancy_remains_above_90_percent`: **PASS**
- `endpoint_identity_retention_remains_above_90_percent`: **PASS**
- `noninitial_lumen_waters_demonstrate_exchange`: **PASS**
- `minimum_CAP_OW_distance_is_safe`: **PASS**
- `post5_temperature_mean_is_stable`: **PASS**
- `post5_temperature_std_is_canonical`: **PASS**
- `post5_temperature_slope_is_small`: **PASS**
- `HBN_remained_frozen`: **PASS**
- `PYR_remained_frozen`: **PASS**
- `CAPS_remained_frozen`: **PASS**
- `CAP_SOL_LJ_remained_below_100kJmol`: **PASS**

## Decision

- Audit decision:
  **R2_OCCUPANCY_TRANSIENT_CHECKPOINT_EXTENSION_JUSTIFIED**
- R2 20 ps validation status:
  **NOT VALIDATED — NONSTATIONARY**
- MD rerun required:
  **NO**
- Checkpoint-continuation preparation authorized:
  **YES**
- Checkpoint-continuation execution authorized:
  **NO**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `PREPARE_R2_30PS_CHECKPOINT_CONTINUATION_TO_50PS`

The proposed extension must use the exact 20 ps checkpoint and must not
regenerate velocities. Its purpose is to determine whether occupancy
approaches a plateau or continues to decline.
