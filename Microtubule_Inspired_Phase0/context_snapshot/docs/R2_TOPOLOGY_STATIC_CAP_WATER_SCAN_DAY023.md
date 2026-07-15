# R2 Topology and Static CAP–Water Scan

## Scope

The R2 symmetric partial-cap geometry was converted into a complete
GROMACS topology by reusing the validated R1 neutral CAP model.

No molecular dynamics or energy minimization was performed. The
GROMACS run contained zero integration steps and served only as an
independent static energy evaluation.

## System

- Total atoms:
  **68332**
- Waters:
  **16565**
- CAPL/CAPU beads:
  **144/144**
- Total cap beads:
  **288**

## Reused CAP–OW model

- Sigma:
  **-0.170000000000 nm**
- Epsilon:
  **3.117923481807 kJ/mol**
- C12:
  **7.266286217926e-09 kJ mol^-1 nm^12**

The negative sigma suppresses the C6 term under combination-rule 2,
leaving the validated purely repulsive C12/r12 interaction.

## Static initial-state results

- Minimum CAP–OW distance:
  **0.220189 nm**
- Analytic CAP–SOL LJ:
  **14.024857060 kJ/mol**
- GROMACS CAP–SOL LJ:
  **14.024806976 kJ/mol**
- Absolute difference:
  **5.008329091e-05 kJ/mol**
- Relative difference:
  **3.571037530e-06**
- Maximum total CAP force on one water oxygen:
  **34.234493
  kJ mol^-1 nm^-1**
- CAP–SOL Coulomb:
  **0.000000000
  kJ/mol**
- CAP–HBN/PYR LJ:
  **0.000000000
  kJ/mol**
- CAP–HBN/PYR Coulomb:
  **0.000000000
  kJ/mol**

## Actual potential-defined aperture

At the selected 5 kBT CAP–OW model:

- Lower aperture radius:
  **0.419703 nm**
- Upper aperture radius:
  **0.419703 nm**
- Mean aperture diameter:
  **0.839406
  nm**
- Lower maximum centerline barrier:
  **0.000084
  kBT**
- Upper maximum centerline barrier:
  **0.000084
  kBT**

The potential-defined aperture incorporates the summed interaction of
all cap beads and is therefore stricter than the nearest-bead geometric
estimate.

## Gates

- `R2_geometry_gate_passed`: **PASS**
- `R2_geometry_authorized_topology`: **PASS**
- `R2_GRO_has_68332_atoms`: **PASS**
- `R2_caps_GRO_has_288_beads`: **PASS**
- `CAPL_topology_has_144_atoms`: **PASS**
- `CAPU_topology_has_144_atoms`: **PASS**
- `SOL_count_is_16565`: **PASS**
- `CAPL_and_CAPU_counts_are_one`: **PASS**
- `CAP_OW_function_is_one`: **PASS**
- `CAP_OW_sigma_is_negative_0p17nm`: **PASS**
- `CAP_OW_epsilon_matches_R1_selected_model`: **PASS**
- `grompp_return_code_zero`: **PASS**
- `TPR_dump_return_code_zero`: **PASS**
- `TPR_has_68332_atoms`: **PASS**
- `mdrun_return_code_zero`: **PASS**
- `mdrun_finished`: **PASS**
- `no_instability_signatures`: **PASS**
- `minimum_CAP_OW_distance_is_valid`: **PASS**
- `selected_CAP_SOL_energy_is_finite`: **PASS**
- `selected_CAP_SOL_energy_is_below_100kJmol`: **PASS**
- `selected_maximum_water_force_is_below_250`: **PASS**
- `GROMACS_and_analytic_CAP_SOL_agree`: **PASS**
- `CAP_SOL_Coulomb_is_zero`: **PASS**
- `CAP_HBNPYR_LJ_is_zero`: **PASS**
- `CAP_HBNPYR_Coulomb_is_zero`: **PASS**
- `lower_actual_aperture_radius_is_0p20_to_0p45nm`: **PASS**
- `upper_actual_aperture_radius_is_0p20_to_0p45nm`: **PASS**
- `lower_upper_aperture_radii_are_symmetric`: **PASS**
- `central_axial_path_barrier_is_at_most_5kBT`: **PASS**

## Decision

- Decision:
  **R2_TOPOLOGY_AND_STATIC_CAP_WATER_MODEL_VALIDATED**
- Failed gates:
  **NONE**
- Water-only minimization authorized:
  **YES**
- Short frozen-solute NVT authorized:
  **NO**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `RUN_R2_WATER_ONLY_ENERGY_MINIMIZATION`

R2 remains a frozen neutral steric screening design. Static acceptance
does not establish water retention, exchange kinetics, chemical
realizability, or device suitability.
