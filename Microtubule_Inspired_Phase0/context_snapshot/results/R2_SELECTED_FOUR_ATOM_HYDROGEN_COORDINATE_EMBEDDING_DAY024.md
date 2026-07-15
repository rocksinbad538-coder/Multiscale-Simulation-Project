# R2 Four-Atom BN Bridge Hydrogen Coordinate Embedding

## Inventory

- Total coordinate nodes: **2316**
- Heavy atoms: **2112**
- Hydrogen atoms: **204**
- Lower/upper H: **102/102**
- B-H/N-H: **102/102**

## X-H geometry

- Maximum X-H deviation:
  **0.000000000 nm**
- Local H-angle violations:
  **0**

## Nonbonded clearance

- Minimum H-heavy / clashes:
  **0.154841883/0**
- Minimum H-H / clashes:
  **0.230044781/0**

## Inner-H aperture proxy

- Lower:
  **1.008734672 nm**
- Upper:
  **1.031124689 nm**
- Asymmetry:
  **0.022390018 nm**

This is a nuclear geometric proxy, not the effective 5 kBT aperture.

## Gates

- `Gate3M_graph_is_accepted`: **PASS**
- `Gate3O2_semantics_resolution_is_accepted`: **PASS**
- `2316_nodes_received_coordinates`: **PASS**
- `2112_heavy_coordinates_are_unchanged`: **PASS**
- `204_H_coordinates_were_generated`: **PASS**
- `102_H_were_generated_per_end`: **PASS**
- `120_bridge_H_were_generated`: **PASS**
- `30_seed_H_were_generated`: **PASS**
- `30_outer_H_were_generated`: **PASS**
- `24_inner_H_were_generated`: **PASS**
- `102_BH_and_102_NH_bonds_are_present`: **PASS**
- `all_XH_bonds_are_within_0p002nm`: **PASS**
- `all_local_H_angles_are_within70to175deg`: **PASS**
- `all_local_H_heavy_clearance_checks_pass`: **PASS**
- `no_global_H_heavy_clashes`: **PASS**
- `no_global_H_H_clashes`: **PASS**
- `inner_H_aperture_is_lower_upper_symmetric`: **FAIL**
- `inner_H_aperture_is_recorded_as_geometric_proxy_only`: **PASS**
- `no_energy_minimization_or_MD_was_performed`: **PASS**

## Decision

- Decision:
  **R2_SELECTED_FOUR_ATOM_BN_BRIDGE_HYDROGEN_COORDINATES_REQUIRE_REVIEW**
- Failed gates:
  **inner_H_aperture_is_lower_upper_symmetric**
- Candidate is final chemistry:
  **NO**
- Molecular topology generation authorized:
  **NO**
- Energy minimization authorized:
  **NO**
- MD authorized:
  **NO**
- QM authorized:
  **NO**
- Required next step:
  `REFINE_R2_SELECTED_FOUR_ATOM_HYDROGEN_ORIENTATIONS`
