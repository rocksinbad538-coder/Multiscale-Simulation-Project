# R2 Frozen-Solute NVT 20 ps Preparation

## Scope

This stage prepares, but does not execute, a 20 ps frozen-solute NVT
screen for R2.

The input state is the validated R2 water-only minimized structure.
HBN, all four pyrenes, and both partial-cap assemblies will remain
frozen. Only TIP4P/2005 water has nonzero dynamical degrees of freedom.

## Protocol

- Atoms:
  **68332**
- Waters:
  **16565**
- CAPL/CAPU:
  **144/144 beads**
- Integrator:
  **md**
- Time step:
  **0.0005000 ps**
- Steps:
  **40000**
- Duration:
  **20.000 ps**
- Temperature:
  **300.0 K**
- Thermostat:
  **V-rescale**
- Effective thermostat group:
  **SOL**
- Frozen groups:
  **HBN_PYR, CAPS**
- Velocity seed:
  **20260708**
- XTC interval:
  **0.500 ps**
- Expected XTC frames:
  **41**

## Grompp audit

- Probe return code:
  **1**
- Final return code:
  **0**
- Warning count:
  **1**
- Controlled `-maxwarn 1` used:
  **YES**
- SOL degrees of freedom:
  **99387.000**
- Expected SOL degrees of freedom:
  **99387**
- HBN_PYR degrees of freedom:
  **0.000**
- CAPS degrees of freedom:
  **0.000**

Any accepted warning is restricted to the known GROMACS
Verlet-buffer/frozen-particle warning. No unrelated warning is
authorized.

## Generated velocities

- OW RMS speed:
  **0.684687 nm/ps**
- OW maximum speed:
  **1.865740 nm/ps**

Velocity generation has already been encoded in the TPR. Mdrun must
not regenerate velocities.

## Gates

- `R2_water_only_EM_is_validated`: **PASS**
- `R2_EM_authorized_NVT_preparation`: **PASS**
- `source_GRO_has_68332_atoms`: **PASS**
- `source_topology_molecule_counts_are_correct`: **PASS**
- `source_index_group_counts_are_correct`: **PASS**
- `grompp_return_code_zero`: **PASS**
- `grompp_warning_policy_is_valid`: **PASS**
- `TPR_dump_return_code_zero`: **PASS**
- `TPR_has_68332_atoms`: **PASS**
- `TPR_integrator_is_md`: **PASS**
- `TPR_dt_is_0p0005ps`: **PASS**
- `TPR_nsteps_is_40000`: **PASS**
- `TPR_total_time_is_20ps`: **PASS**
- `TPR_continuation_is_false`: **PASS**
- `TPR_nstlog_is_100`: **PASS**
- `TPR_nstenergy_is_100`: **PASS**
- `TPR_XTC_stride_is_1000`: **PASS**
- `TPR_has_all_velocity_entries`: **PASS**
- `water_oxygen_velocities_are_nonzero`: **PASS**
- `SOL_degrees_of_freedom_are_expected`: **PASS**
- `HBN_PYR_degrees_of_freedom_are_zero`: **PASS**
- `CAPS_degrees_of_freedom_are_zero`: **PASS**
- `processed_MDP_contract_is_correct`: **PASS**

## Decision

- Decision:
  **R2_FROZEN_SOLUTE_NVT_20PS_PREPARED**
- Failed gates:
  **NONE**
- NVT execution authorized:
  **YES**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `RUN_R2_FROZEN_SOLUTE_NVT_20PS`

R2 remains a frozen neutral steric screening design. Preparation of
this short trajectory does not establish chemical realizability,
long-time retention, or a final device architecture.
