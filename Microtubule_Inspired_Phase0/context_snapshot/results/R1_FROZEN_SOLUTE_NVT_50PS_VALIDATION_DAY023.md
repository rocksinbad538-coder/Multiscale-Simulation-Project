# R1 Frozen-Solute NVT 50 ps Validation

## Scope

The R1 fully capped neutral steric positive control was continued from
20 to 50 ps using the exact 20 ps checkpoint.

No velocities were regenerated. HBN, all pyrenes, and both caps
remained frozen. Only TIP4P/2005 water was mobile and thermostatted.

## Execution

- Mdrun return code:
  **0**
- Finished mdrun:
  **YES**
- Checkpoint continuation reported:
  **YES**
- Continuation frames:
  **61**
- Continuation time range:
  **20.000–50.000 ps**
- Combined unique frames:
  **101**
- Combined time range:
  **0.000–50.000 ps**
- Instability signatures:
  **NONE**

## Temperature

Continuation 20–50 ps:

- Mean ± standard deviation:
  **299.9643
  ± 1.3557 K**
- Minimum/maximum:
  **296.5693/
  303.7985 K**
- Linear slope:
  **-0.013563 K/ps**

Final 15 ps:

- Mean ± standard deviation:
  **299.7803
  ± 1.2221 K**
- Minimum/maximum:
  **296.5693/
  303.4604 K**
- Linear slope:
  **-0.020140 K/ps**

## Frozen-group integrity

Final RMS/max displacement:

- HBN:
  **0.000000000000/0.000000000000 nm**
- PYR:
  **0.000000000000/0.000000000000 nm**
- CAPS:
  **0.000000000000/0.000000000000 nm**

Maximum displacement in the continuation trajectory:

- HBN:
  **0.000000 nm**
- PYR:
  **0.000000 nm**
- CAPS:
  **0.000000 nm**

## Confinement over 0–50 ps

- Initial lumen occupancy:
  **428 waters**
- Mean ± standard deviation:
  **427.9109
  ± 0.3493 waters**
- Minimum/maximum:
  **426/
  428 waters**
- Endpoint occupancy:
  **428 waters**
- Endpoint initially luminal waters retained:
  **428/428
  (1.000000)**
- Zero-occupancy fraction:
  **0.000000**
- Occupancy slope over 25–50 ps:
  **0.004706 waters/ps**
- Minimum CAP–OW distance:
  **0.170798 nm**

## Energy

- Potential initial/final/change:
  **-765130.250000/
  -762399.375000/
  2730.875000 kJ/mol**
- Total energy initial/final:
  **-641436.062500/
  -638832.625000 kJ/mol**
- CAP–SOL LJ initial/final/maximum:
  **15.390732/
  15.371712/
  36.945129 kJ/mol**

## Gates

- `mdrun_return_code_zero`: **PASS**
- `mdrun_finished`: **PASS**
- `checkpoint_continuation_was_reported`: **PASS**
- `no_instability_signatures`: **PASS**
- `trajectory_check_return_code_zero`: **PASS**
- `continuation_frame_count_is_60_or_61`: **PASS**
- `continuation_last_time_is_50ps`: **PASS**
- `combined_trajectory_has_101_unique_frames`: **PASS**
- `combined_trajectory_starts_at_0ps`: **PASS**
- `combined_trajectory_ends_at_50ps`: **PASS**
- `combined_frame_spacing_is_0p5ps`: **PASS**
- `final_atom_count_is_68314`: **PASS**
- `box_is_unchanged`: **PASS**
- `HBN_final_coordinates_are_frozen`: **PASS**
- `PYR_final_coordinates_are_frozen`: **PASS**
- `CAPS_final_coordinates_are_frozen`: **PASS**
- `HBN_continuation_trajectory_is_frozen`: **PASS**
- `PYR_continuation_trajectory_is_frozen`: **PASS**
- `CAPS_continuation_trajectory_is_frozen`: **PASS**
- `water_coordinates_are_mobile`: **PASS**
- `continuation_temperature_mean_is_295_to_305K`: **PASS**
- `continuation_temperature_std_is_at_most_5K`: **PASS**
- `continuation_temperature_range_is_280_to_320K`: **PASS**
- `continuation_temperature_slope_is_small`: **PASS**
- `last15ps_temperature_mean_is_295_to_305K`: **PASS**
- `last15ps_temperature_std_is_at_most_5K`: **PASS**
- `last15ps_temperature_slope_is_small`: **PASS**
- `potential_series_is_finite`: **PASS**
- `total_energy_series_is_finite`: **PASS**
- `CAP_SOL_energy_series_is_finite`: **PASS**
- `CAP_SOL_energy_remains_below_500kJmol`: **PASS**
- `initial_lumen_occupancy_is_428`: **PASS**
- `no_complete_lumen_drying_over_0_to_50ps`: **PASS**
- `combined_minimum_occupancy_is_at_least_90_percent`: **PASS**
- `continuation_endpoint_occupancy_is_at_least_90_percent`: **PASS**
- `continuation_endpoint_initial_retention_is_at_least_90_percent`: **PASS**
- `combined_second_half_occupancy_slope_is_small`: **PASS**
- `CAP_OW_distance_remains_above_0p15nm`: **PASS**

## Decision

- Decision:
  **R1_FROZEN_SOLUTE_50PS_POSITIVE_CONTROL_VALIDATED**
- Failed gates:
  **NONE**
- R2 partial-cap static design authorized:
  **YES**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `BEGIN_R2_PARTIAL_CAP_DESIGN_AND_STATIC_GATE`

A passing result validates R1 only as a frozen neutral steric positive
control for the confinement-analysis methodology. It does not establish
that R1 is chemically realizable or appropriate as the final device.
