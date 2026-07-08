# R1 NVT Grompp Warning Authorization

## Warning reviewed

The only GROMACS warning is:

> Some temperature coupling groups do not use temperature coupling. We will assume their temperature is not more than 300.000 K. If their temperature is higher, the energy error and the Verlet buffer might be underestimated.

This warning is authorized only for this R1 frozen-solute screening
protocol because:

- HBN_PYR has 0.0 degrees of freedom;
- CAPS has 0.0 degrees of freedom;
- both groups are frozen in all Cartesian dimensions;
- both groups use tau-t = -1 and are not thermostatted;
- SOL is the only mobile and thermostatted group;
- no additional warning or fatal error is present.

The use of `-maxwarn 1` is therefore restricted to this exact,
programmatically audited warning.

## Coupling state

- T-coupling groups:
  **SOL HBN_PYR CAPS**
- tau-t:
  **0.1 -1.0 -1.0 ps**
- ref-t:
  **300.0 300.0 300.0 K**
- Freeze groups:
  **HBN_PYR CAPS**
- Freeze dimensions:
  **Y Y Y Y Y Y**

## Degrees of freedom

- SOL:
  **99303.0**
- HBN_PYR:
  **0.0**
- CAPS:
  **0.0**

## Gates

- `TPR_exists_and_is_nonempty`: **PASS**
- `exactly_one_grompp_warning`: **PASS**
- `warning_text_is_expected`: **PASS**
- `no_fatal_error_after_authorization`: **PASS**
- `grompp_used_exactly_maxwarn_1`: **PASS**
- `temperature_groups_are_correct`: **PASS**
- `tau_t_values_are_correct`: **PASS**
- `reference_temperatures_are_300K`: **PASS**
- `freeze_groups_are_correct`: **PASS**
- `all_six_freeze_dimensions_are_active`: **PASS**
- `SOL_has_positive_degrees_of_freedom`: **PASS**
- `HBN_PYR_has_zero_degrees_of_freedom`: **PASS**
- `CAPS_has_zero_degrees_of_freedom`: **PASS**
- `expected_CAPL_unbound_note_is_present`: **PASS**
- `expected_CAPU_unbound_note_is_present`: **PASS**
- `expected_VCM_note_is_present`: **PASS**
- `preparation_decision_is_pass`: **PASS**
- `preparation_authorized_execution`: **PASS**

## Decision

- Decision:
  **R1_NVT_GROMPP_WARNING_AUTHORIZED**
- Failed gates:
  **NONE**
- NVT execution authorized:
  **YES**
- Required next step:
  `RUN_R1_FROZEN_SOLUTE_NVT_20PS`
