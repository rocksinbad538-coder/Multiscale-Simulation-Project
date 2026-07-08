# R2 Frozen-Solute NVT — Combined 0–50 ps Validation

## Execution scope

The accepted 0–20 ps trajectory was not repeated. The simulation was
continued from the exact 20 ps checkpoint to the 50 ps target using the
extended TPR.

- Source checkpoint:
  **step 40000, 20.0 ps**
- Final checkpoint:
  **step 100000, 50.000000 ps**
- Velocity regeneration:
  **NO**
- Thermostat-state regeneration:
  **NO**
- Source 0–20 ps rerun:
  **NO**
- Continuation instability signatures:
  **NONE**

## Trajectory integrity

- Source raw frames:
  **41**
- Continuation raw frames:
  **61**
- Continuation time range:
  **20.000–
  50.000 ps**
- Combined unique frames:
  **101**
- Combined time range:
  **0.000–50.000 ps**
- Combined frame interval:
  **0.500000 ps**

## Thermal behavior

- Continuation mean/std/min/max:
  **299.8727/
  1.3326/
  295.8250/
  303.3355 K**
- Continuation slope:
  **-0.013063 K/ps**
- Final 15 ps mean/std/slope:
  **299.8449/
  1.3533/
  -0.044881 K, K, K/ps**
- Theoretical canonical temperature standard deviation:
  **1.345772 K**

## Frozen-group integrity

Maximum displacement over the combined trajectory:

- HBN:
  **0.000000000000 nm**
- PYR:
  **0.000000000000 nm**
- CAPS:
  **0.000000000000 nm**

Water-O final RMS/max displacement:

- **0.840909/3.275839 nm**

## Lumen-water behavior

- Initial/mean/minimum/maximum/endpoint occupancy:
  **428/
  416.6337/
  409/
  429/
  411**
- Endpoint retained initial identities:
  **409/428**
- Endpoint identity-retention fraction:
  **0.955607**
- Minimum CAP–OW distance:
  **0.166949 nm**

### Cumulative windows

- 0.0–50.0 ps: mean=416.6337; min/max=409/429; change=-17; slope=-0.305859 waters/ps
- 10.0–50.0 ps: mean=414.5309; min/max=409/422; change=-10; slope=-0.193451 waters/ps
- 20.0–50.0 ps: mean=413.7869; min/max=409/420; change=-5; slope=-0.238710 waters/ps
- 25.0–50.0 ps: mean=413.0980; min/max=409/418; change=-5; slope=-0.217376 waters/ps
- 30.0–50.0 ps: mean=412.4390; min/max=409/417; change=-5; slope=-0.173868 waters/ps
- 35.0–50.0 ps: mean=411.6452; min/max=409/415; change=-2; slope=0.025806 waters/ps
- 40.0–50.0 ps: mean=411.7143; min/max=409/415; change=+0; slope=0.064935 waters/ps

### Ten-picosecond blocks

- 0.0–10.0 ps: mean=425.1500; min/max=421/429; change=-6; slope=-0.735338 waters/ps
- 10.0–20.0 ps: mean=416.8000; min/max=414/422; change=-5; slope=-0.565414 waters/ps
- 20.0–30.0 ps: mean=416.5500; min/max=414/420; change=+0; slope=-0.218045 waters/ps
- 30.0–40.0 ps: mean=413.2000; min/max=409/417; change=-7; slope=-0.637594 waters/ps
- 40.0–50.0 ps: mean=411.7143; min/max=409/415; change=+0; slope=0.064935 waters/ps

## Continuation energetics

- Potential initial/final/change:
  **-767635.250000/
  -768991.062500/
  -1355.812500 kJ/mol**
- CAP–SOL LJ initial/final/maximum:
  **30.687618/
  25.050507/
  49.542679 kJ/mol**

## Gates

- `authorized_contract_is_valid`: **PASS**
- `source_checkpoint_hash_matches_contract`: **PASS**
- `velocity_regeneration_remained_disabled`: **PASS**
- `thermostat_state_regeneration_remained_disabled`: **PASS**
- `source_0_to_20ps_was_not_rerun`: **PASS**
- `mdrun_return_code_zero`: **PASS**
- `mdrun_completion_confirmed`: **PASS**
- `checkpoint_continuation_confirmed_in_logs`: **PASS**
- `no_instability_signatures`: **PASS**
- `continuation_trajectory_check_return_code_zero`: **PASS**
- `final_checkpoint_dump_return_code_zero`: **PASS**
- `final_checkpoint_step_is_100000`: **PASS**
- `final_checkpoint_time_is_50ps`: **PASS**
- `source_trajectory_has_41_frames`: **PASS**
- `continuation_trajectory_has_60_or_61_frames`: **PASS**
- `continuation_starts_at_20_or_20p5ps`: **PASS**
- `continuation_ends_at_50ps`: **PASS**
- `combined_trajectory_has_101_unique_frames`: **PASS**
- `combined_trajectory_starts_at_0ps`: **PASS**
- `combined_trajectory_ends_at_50ps`: **PASS**
- `combined_trajectory_interval_is_0p5ps`: **PASS**
- `HBN_remained_frozen`: **PASS**
- `PYR_remained_frozen`: **PASS**
- `CAPS_remained_frozen`: **PASS**
- `water_displacement_is_finite`: **PASS**
- `minimum_CAP_OW_distance_is_safe`: **PASS**
- `continuation_temperature_mean_is_295_to_305K`: **PASS**
- `continuation_temperature_std_is_acceptable`: **PASS**
- `continuation_temperature_slope_is_acceptable`: **PASS**
- `final15_temperature_mean_is_295_to_305K`: **PASS**
- `final15_temperature_std_is_acceptable`: **PASS**
- `final15_temperature_slope_is_acceptable`: **PASS**
- `initial_lumen_occupancy_is_428`: **PASS**
- `combined_minimum_occupancy_remains_above_80_percent`: **PASS**
- `final20_mean_occupancy_is_at_least_90_percent`: **PASS**
- `endpoint_occupancy_is_at_least_90_percent`: **PASS**
- `endpoint_initial_identity_retention_is_at_least_50_percent`: **PASS**
- `final20_occupancy_slope_is_stationary`: **PASS**
- `final15_occupancy_slope_is_stationary`: **PASS**
- `final10_occupancy_slope_is_stationary`: **PASS**
- `final10_net_occupancy_change_is_at_most_5_waters`: **PASS**
- `temperature_series_is_finite`: **PASS**
- `potential_series_is_finite`: **PASS**
- `total_energy_series_is_finite`: **PASS**
- `CAP_SOL_LJ_series_is_finite`: **PASS**
- `CAP_SOL_LJ_remains_below_100kJmol`: **PASS**

## Decision

- Decision:
  **R2_FROZEN_SOLUTE_NVT_50PS_VALIDATED**
- Failed gates:
  **NONE**
- R2-versus-R1 comparison authorized:
  **YES**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `COMPARE_R2_WITH_R1_AND_DECIDE_NEXT_ARCHITECTURE_GATE`

R2 remains a neutral frozen steric screening architecture. Even if the
50 ps gate passes, this result does not establish chemical realizability
or authorize long mobile, multitemperature, or quantum calculations.
