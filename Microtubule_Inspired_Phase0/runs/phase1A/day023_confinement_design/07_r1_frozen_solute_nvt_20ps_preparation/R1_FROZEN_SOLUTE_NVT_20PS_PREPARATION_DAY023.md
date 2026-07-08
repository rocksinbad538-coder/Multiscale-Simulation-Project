# R1 Frozen-Solute NVT 20 ps Preparation

## Scope

This gate prepared and statically validated the first short R1
confinement screening trajectory.

No molecular dynamics was executed.

## Inputs

- Minimized coordinates:
  `runs/phase1A/day023_confinement_design/07_r1_frozen_solute_nvt_20ps_preparation/r1_frozen_solute_nvt_20ps_input.gro`
- Selected topology:
  `runs/phase1A/day023_confinement_design/07_r1_frozen_solute_nvt_20ps_preparation/r1_frozen_solute_nvt_20ps.top`
- Index:
  `runs/phase1A/day023_confinement_design/07_r1_frozen_solute_nvt_20ps_preparation/r1_frozen_solute_nvt_20ps.ndx`
- NVT template:
  `runs/phase1A/day021_mobile_restraint_protocol/protocol_inputs/mdp/02_nvt_k10000_1ps.mdp`

## Protocol

- Integrator: **MD**
- Temperature: **300.0 K**
- Time step: **0.0005000 ps**
- Steps: **40000**
- Duration: **20.0 ps**
- Trajectory interval:
  **0.5 ps**
- Expected frames:
  **41**
- Velocity seed:
  **20260708**
- Frozen groups:
  **HBN_PYR CAPS**
- T-coupling partition:
  **SOL HBN_PYR CAPS**
- Mobile and thermostatted group:
  **SOL, tau-t = 0.100 ps**
- Frozen and unthermostatted groups:
  **HBN_PYR and CAPS, tau-t = -1 ps**
- Pressure coupling:
  **disabled**

## TPR and topology audit

- GRO atoms:
  **68314**
- TPR atoms:
  **68314**
- Active position-restraint sections:
  **0**
- Processed SOL/CAPL/CAPU counts:
  **16551/
  1/
  1**
- CAP–OW override preserved:
  **YES**

## Velocity audit

- `HBN`: atoms=1680; nonzero fraction=1.000000; RMS speed=0.775970 nm/ps
- `PYR`: atoms=104; nonzero fraction=1.000000; RMS speed=1.569049 nm/ps
- `SOL`: atoms=66204; nonzero fraction=0.750000; RMS speed=1.956395 nm/ps
- `CAPS`: atoms=326; nonzero fraction=1.000000; RMS speed=0.804803 nm/ps

The expected nonzero fraction for four-site TIP4P/2005 water is
approximately 0.75 because the virtual M site has no independently
generated velocity.

Frozen-group velocities are recorded for provenance. Their coordinates
remain fixed by the six active freeze dimensions.

## Energy-minimization reporting note

- Potential-energy records in the EM XVG:
  **1**
- Initial/final EM energy change independently resolvable:
  **NO**

If only one energy record is present, the previously reported zero
energy change is a sampling limitation, not evidence that the
minimization failed. EM acceptance remains based on explicit
convergence, final force, structural integrity, and preserved hydration.

## Gates

- `grompp_return_code_zero`: **PASS**
- `TPR_dump_return_code_zero`: **PASS**
- `input_GRO_atom_count_is_68314`: **PASS**
- `TPR_atom_count_is_68314`: **PASS**
- `all_TPR_velocities_parsed`: **PASS**
- `all_TPR_velocities_are_finite`: **PASS**
- `water_velocity_fraction_is_TIP4P_consistent`: **PASS**
- `integrator_is_md`: **PASS**
- `dt_is_0p0005_ps`: **PASS**
- `nsteps_is_40000`: **PASS**
- `continuation_is_no`: **PASS**
- `generation_temperature_is_300K`: **PASS**
- `generation_seed_is_fixed`: **PASS**
- `freeze_groups_are_HBN_PYR_and_CAPS`: **PASS**
- `all_freeze_dimensions_are_enabled`: **PASS**
- `temperature_groups_partition_mobile_and_frozen_atoms`: **PASS**
- `mobile_group_tau_t_is_0p1_ps`: **PASS**
- `frozen_groups_are_not_thermostatted`: **PASS**
- `reference_temperatures_have_three_300K_entries`: **PASS**
- `pressure_coupling_is_disabled`: **PASS**
- `trajectory_interval_is_0p5ps`: **PASS**
- `no_active_position_restraints`: **PASS**
- `SOL_count_is_16551`: **PASS**
- `CAPL_count_is_one`: **PASS**
- `CAPU_count_is_one`: **PASS**
- `selected_CAP_OW_override_is_preserved`: **PASS**

## Decision

- Decision:
  **R1_FROZEN_SOLUTE_NVT_20PS_PREPARED**
- Failed gates:
  **NONE**
- NVT execution authorized:
  **YES**
- Required next step:
  `RUN_R1_FROZEN_SOLUTE_NVT_20PS`

The 20 ps run is an initial positive-control screening. Extension to a
longer frozen-solute trajectory will require validation of temperature,
cap integrity, lumen occupancy, zero-occupancy fraction, and axial
water retention.
