# R2 Inner-H Reflected-Direction Refinement

## Applied refinement

- Scenario:
  **UPPER_DRIVES_LOWER**
- Modified coordinates:
  **12 LOWER inner-rim H**
- Preserved heavy coordinates:
  **True**
- Preserved remaining H coordinates:
  **True**

## Bond and angle geometry

- B-H/N-H bonds:
  **102/102**
- Maximum X-H deviation:
  **8.881784197001e-16 nm**
- H-angle range:
  **70.026550401–170.273738102 degrees**
- H-angle violations:
  **0**

## Nonbonded clearance

- Minimum H-heavy / clashes:
  **0.154841883 nm / 0**
- Minimum H-H / clashes:
  **0.230044781 nm / 0**

## Inner-H nuclear aperture proxy

- LOWER:
  **1.028648811 nm**
- UPPER:
  **1.031124689 nm**
- Asymmetry:
  **0.002475879 nm**

This remains a nuclear geometric proxy and is not the effective
water-accessible aperture at 5 kBT.

## Gates

- `diagnostic_selected_UPPER_DRIVES_LOWER`: **PASS**
- `2316_nodes_have_coordinates`: **PASS**
- `2112_heavy_coordinates_are_unchanged`: **PASS**
- `exactly_12_lower_inner_H_coordinates_were_modified`: **PASS**
- `remaining_192_H_coordinates_are_unchanged`: **PASS**
- `102_BH_and_102_NH_bonds_are_present`: **PASS**
- `all_XH_bonds_are_within_0p002nm`: **PASS**
- `all_H_angles_are_within70to175deg`: **PASS**
- `no_global_H_heavy_clashes`: **PASS**
- `no_global_H_H_clashes`: **PASS**
- `lower_upper_inner_H_aperture_asymmetry_is_at_most0p010nm`: **PASS**
- `inner_H_aperture_is_recorded_as_geometric_proxy_only`: **PASS**
- `no_topology_charges_parameterization_minimization_MD_or_QM`: **PASS**

## Decision

- Decision:
  **R2_SELECTED_FOUR_ATOM_BN_BRIDGE_HYDROGEN_COORDINATES_VALIDATED_AFTER_SYMMETRY_REFINEMENT**
- Failed gates:
  **NONE**
- Molecular topology generation authorized:
  **NO**
- Formal charge assignment authorized:
  **NO**
- Force-field parameterization authorized:
  **NO**
- Energy minimization authorized:
  **NO**
- MD authorized:
  **NO**
- QM authorized:
  **NO**
- Required next step:
  `AUDIT_R2_SELECTED_FOUR_ATOM_BN_BRIDGE_CHEMICAL_REALIZABILITY_AND_PARAMETERIZATION_SCOPE`
