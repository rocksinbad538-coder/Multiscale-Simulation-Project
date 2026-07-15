# R2 Water-Only Energy Minimization

## Scope

The R2 symmetric partial-cap system was minimized with HBN, all four
pyrenes, and both cap assemblies frozen in all dimensions.

Only the 16565 TIP4P/2005 water molecules were allowed to move.

## Execution

- Grompp return code:
  **0**
- Grompp warnings:
  **0**
- Mdrun return code:
  **0**
- Steepest-descents convergence:
  **YES**
- EM steps:
  **1**
- Maximum force:
  **382.908780 kJ mol^-1 nm^-1**
- Norm of force:
  **16.557585 kJ mol^-1 nm^-1**

## Frozen-group integrity

RMS/max displacement:

- HBN:
  **0.000000000000/0.000000000000 nm**
- PYR:
  **0.000000000000/0.000000000000 nm**
- CAPS:
  **0.000000000000/0.000000000000 nm**
- Water oxygen:
  **0.000029404/0.000530000 nm**

## Water confinement

- Initial lumen occupancy:
  **428**
- Final lumen occupancy:
  **428**
- Initially luminal waters retained:
  **428/428**
- Initial retention fraction:
  **1.000000**
- Initial CAP–OW distance:
  **0.220189 nm**
- Final CAP–OW distance:
  **0.220189 nm**

## Energy

- Energy records:
  **1**
- Potential initial/final/change:
  **-904585.312500/
  -904585.312500/
  0.000000 kJ/mol**
- CAP–SOL LJ initial/final/maximum:
  **14.033027/
  14.033027/
  14.033027 kJ/mol**
- Maximum absolute CAP–SOL Coulomb:
  **0.000000000 kJ/mol**

## Gates

- `R2_static_gate_is_validated`: **PASS**
- `R2_strict_instability_audit_is_clean`: **PASS**
- `grompp_return_code_zero`: **PASS**
- `grompp_warning_count_zero`: **PASS**
- `mdrun_return_code_zero`: **PASS**
- `no_instability_signatures`: **PASS**
- `steepest_descents_converged`: **PASS**
- `maximum_force_is_finite`: **PASS**
- `maximum_force_is_at_most_emtol`: **PASS**
- `norm_force_is_finite`: **PASS**
- `initial_atom_count_is_68332`: **PASS**
- `final_atom_count_is_68332`: **PASS**
- `box_is_unchanged`: **PASS**
- `HBN_is_exactly_frozen`: **PASS**
- `PYR_is_exactly_frozen`: **PASS**
- `CAPS_are_exactly_frozen`: **PASS**
- `water_displacement_is_finite`: **PASS**
- `water_displacement_is_local`: **PASS**
- `initial_lumen_occupancy_is_428`: **PASS**
- `final_lumen_occupancy_retains_at_least_98_percent`: **PASS**
- `initial_lumen_identity_retains_at_least_98_percent`: **PASS**
- `initial_CAP_OW_distance_is_safe`: **PASS**
- `final_CAP_OW_distance_is_safe`: **PASS**
- `potential_series_is_finite`: **PASS**
- `CAP_SOL_LJ_series_is_finite`: **PASS**
- `CAP_SOL_LJ_remains_below_100kJmol`: **PASS**
- `CAP_SOL_Coulomb_is_zero`: **PASS**

## Decision

- Decision:
  **R2_WATER_ONLY_ENERGY_MINIMIZATION_VALIDATED**
- Failed gates:
  **NONE**
- Short frozen-solute NVT preparation authorized:
  **YES**
- Short frozen-solute NVT execution authorized:
  **NO**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `PREPARE_R2_FROZEN_SOLUTE_NVT_20PS`

This result applies only to the neutral frozen steric R2 screening
model. It does not establish chemical realizability or long-time water
retention.
