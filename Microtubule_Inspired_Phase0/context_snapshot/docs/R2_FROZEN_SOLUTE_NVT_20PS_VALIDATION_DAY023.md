# R2 Frozen-Solute NVT — 20 ps Validation

## Scope

R2 was screened for 20 ps with HBN, all four pyrenes, and both
partial-cap assemblies frozen. Only TIP4P/2005 water was mobile and
thermostatted.

## Execution

- Mdrun return code:
  **0**
- Trajectory-check return code:
  **0**
- Completion confirmed:
  **YES**
- Checkpoint written:
  **YES**
- Instability signatures:
  **NONE**
- Frames:
  **41**
- Time range:
  **0.000–20.000 ps**
- Mean frame interval:
  **0.500000 ps**

## Thermal behavior

- Full mean/std/min/max:
  **298.1977/
  10.9859/
  179.1389/
  303.7430 K**
- 5–20 ps mean/std/min/max:
  **299.7600/
  1.2273/
  296.2640/
  303.7430 K**
- 5–20 ps slope:
  **-0.046683 K/ps**
- Theoretical canonical standard deviation:
  **1.345772 K**
- 10–20 ps mean/std/slope:
  **299.6069/
  1.2091/
  -0.018416 K, K, K/ps**

## Frozen-group integrity

Maximum displacement over all frames:

- HBN:
  **0.000000000000 nm**
- PYR:
  **0.000000000000 nm**
- CAPS:
  **0.000000000000 nm**

Water-O final RMS/max displacement:

- **0.563030/2.998818 nm**

## Lumen-water behavior

- Initial/mean/minimum/maximum/endpoint occupancy:
  **428/
  420.8537/
  414/
  429/
  416**
- Second-half mean:
  **416.7619 waters**
- Second-half slope:
  **-0.509091 waters/ps**
- Endpoint retention of initially luminal waters:
  **414/428**
- Endpoint initial-identity retention fraction:
  **0.967290**
- Minimum CAP–OW distance:
  **0.169759 nm**

## Energetics

- Potential initial/final/change:
  **-904577.625000/
  -767611.562500/
  136966.062500 kJ/mol**
- CAP–SOL LJ initial/final/maximum:
  **14.025614/
  30.687595/
  74.219818 kJ/mol**

## Gates

- `preparation_decision_is_valid`: **PASS**
- `run_contract_authorized_execution`: **PASS**
- `velocity_regeneration_remained_disabled`: **PASS**
- `mdrun_return_code_zero`: **PASS**
- `mdrun_completion_confirmed`: **PASS**
- `checkpoint_was_written`: **PASS**
- `no_instability_signatures`: **PASS**
- `trajectory_check_return_code_zero`: **PASS**
- `trajectory_has_41_frames`: **PASS**
- `trajectory_starts_at_0ps`: **PASS**
- `trajectory_ends_at_20ps`: **PASS**
- `trajectory_interval_is_0p5ps`: **PASS**
- `HBN_remained_frozen`: **PASS**
- `PYR_remained_frozen`: **PASS**
- `CAPS_remained_frozen`: **PASS**
- `water_is_mobile`: **PASS**
- `water_displacement_is_finite`: **PASS**
- `post5_temperature_mean_is_295_to_305K`: **PASS**
- `post5_temperature_std_is_acceptable`: **PASS**
- `post5_temperature_slope_is_acceptable`: **PASS**
- `initial_lumen_occupancy_is_428`: **PASS**
- `trajectory_does_not_approach_complete_drying`: **PASS**
- `second_half_mean_occupancy_is_at_least_90_percent`: **PASS**
- `endpoint_occupancy_is_at_least_90_percent`: **PASS**
- `endpoint_identity_retention_is_at_least_50_percent`: **PASS**
- `second_half_occupancy_slope_is_acceptable`: **FAIL**
- `minimum_CAP_OW_distance_is_safe`: **PASS**
- `temperature_series_is_finite`: **PASS**
- `potential_series_is_finite`: **PASS**
- `total_energy_series_is_finite`: **PASS**
- `CAP_SOL_LJ_series_is_finite`: **PASS**
- `CAP_SOL_LJ_remains_below_100kJmol`: **PASS**

## Decision

- Decision:
  **R2_FROZEN_SOLUTE_NVT_20PS_REQUIRES_REVIEW**
- Failed gates:
  **second_half_occupancy_slope_is_acceptable**
- Checkpoint-continuation preparation authorized:
  **NO**
- Checkpoint-continuation execution authorized:
  **NO**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `REVIEW_R2_FROZEN_SOLUTE_NVT_20PS_GATE_FAILURES`

R2 remains a neutral frozen steric screening architecture. Passing
this gate would justify only preparation of the matched checkpoint
continuation from 20 to 50 ps.
