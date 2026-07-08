# R1 Frozen-Solute NVT 20 ps Validation

## Scope

The R1 fully capped steric positive control was simulated for
20 ps at 300 K.

Frozen throughout the trajectory:

- HBN;
- all pyrene chromophores;
- both steric caps.

Mobile and thermostatted:

- 16551 TIP4P/2005 water molecules.

This is a confinement-methodology screening, not a chemically
realizable final device model.

## Execution

- Mdrun return code:
  **0**
- Finished mdrun:
  **YES**
- Frames:
  **41/41**
- Time range:
  **0.0–20.0 ps**
- Instability signatures:
  **NONE**

## Temperature and energy

- Temperature mean ± standard deviation:
  **298.4674 ± 10.9928 K**
- Temperature minimum/maximum:
  **179.0396/303.4521 K**
- Potential initial/final/change:
  **-899823.062500/
  -765130.250000/
  134692.812500 kJ/mol**
- Total energy initial/final:
  **-775952.187500/
  -641436.062500 kJ/mol**
- CAP–SOL LJ initial/final/maximum:
  **14.305803/
  15.390732/
  45.470341 kJ/mol**

## Frozen-group integrity

Final RMS/max displacement:

- HBN:
  **0.000000000000/0.000000000000 nm**
- PYR:
  **0.000000000000/0.000000000000 nm**
- CAPS:
  **0.000000000000/0.000000000000 nm**

Maximum displacement observed in compressed trajectory:

- HBN:
  **0.000000 nm**
- PYR:
  **0.000000 nm**
- CAPS:
  **0.000000 nm**

## Water motion and confinement

- Final water-O RMS/max displacement:
  **0.544019/2.019448 nm**
- Initial lumen occupancy:
  **428 waters**
- Mean ± standard deviation:
  **427.8293 ± 0.4951 waters**
- Minimum/maximum occupancy:
  **426/428 waters**
- Endpoint occupancy:
  **428 waters**
- Zero-occupancy fraction:
  **0.000000**
- Initial lumen waters retained at endpoint:
  **428/428
  (1.000000)**
- Second-half mean occupancy:
  **427.8571 waters**
- Second-half occupancy slope:
  **-0.046753 waters/ps**
- Minimum CAP–OW distance:
  **0.170798 nm**

## Validation gates

- `mdrun_return_code_zero`: **PASS**
- `mdrun_finished`: **PASS**
- `no_instability_signatures`: **PASS**
- `trajectory_check_return_code_zero`: **PASS**
- `final_atom_count_is_68314`: **PASS**
- `trajectory_has_41_frames`: **PASS**
- `trajectory_starts_at_0ps`: **PASS**
- `trajectory_ends_at_20ps`: **PASS**
- `box_is_unchanged`: **PASS**
- `HBN_final_coordinates_are_frozen`: **PASS**
- `PYR_final_coordinates_are_frozen`: **PASS**
- `CAPS_final_coordinates_are_frozen`: **PASS**
- `HBN_trajectory_is_frozen`: **PASS**
- `PYR_trajectory_is_frozen`: **PASS**
- `CAPS_trajectory_is_frozen`: **PASS**
- `temperature_series_is_finite`: **PASS**
- `temperature_mean_is_295_to_305K`: **PASS**
- `temperature_standard_deviation_is_acceptable`: **FAIL**
- `potential_series_is_finite`: **PASS**
- `total_energy_series_is_finite`: **PASS**
- `CAP_SOL_energy_series_is_finite`: **PASS**
- `initial_lumen_occupancy_is_428`: **PASS**
- `no_complete_lumen_drying`: **PASS**
- `minimum_lumen_occupancy_is_at_least_90_percent`: **PASS**
- `endpoint_lumen_occupancy_is_at_least_90_percent`: **PASS**
- `endpoint_initial_lumen_retention_is_at_least_90_percent`: **PASS**
- `CAP_OW_distance_remains_above_0p15nm`: **PASS**
- `water_coordinates_are_mobile`: **PASS**

## Decision

- Decision:
  **R1_FROZEN_SOLUTE_NVT_20PS_REQUIRES_REVIEW**
- Failed gates:
  **temperature_standard_deviation_is_acceptable**
- Preparation of a 50 ps frozen-solute extension authorized:
  **NO**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `REVIEW_R1_NVT_20PS_GATE_FAILURES`

The R1 caps remain a neutral frozen steric positive control. A positive
20 ps result demonstrates that the retention-analysis methodology can
detect a closed-boundary confinement condition; it does not establish
chemical realizability or justify a final device architecture.
