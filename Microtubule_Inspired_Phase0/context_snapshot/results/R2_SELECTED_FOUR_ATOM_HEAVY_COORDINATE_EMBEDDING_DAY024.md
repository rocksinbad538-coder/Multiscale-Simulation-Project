# R2 Selected Four-Atom BN Bridge Heavy Coordinate Embedding

## Scope

This gate combines the fixed parent/seed/annulus coordinates with the
exact 120 bridge coordinates recovered from Gate 3L.

Hydrogen coordinates are not generated in this gate.

## Coordinate inventory

- Fixed heavy atoms: **1992**
- Four-atom bridge atoms: **120**
- Total heavy atoms: **2112**
- Bridge paths: **30**

## Bond geometry

- Maximum B-N deviation:
  **0.001195394 nm**
- Maximum bridge B-N deviation:
  **0.000044107 nm**

## Critical angles

- Minimum/mean/maximum:
  **70.899840/
  117.353534/
  168.181473 degrees**
- RMS deviation from 120 degrees:
  **16.110536 degrees**

## Heavy-atom clearance

- Minimum nonbonded heavy-heavy distance:
  **0.159747397 nm**
- Heavy-heavy clashes:
  **0**

## Aperture and symmetry

- Lower/upper aperture:
  **1.045415170/
  1.045415170 nm**
- Maximum lower-upper asymmetry:
  **0.000000000 nm**

## Gates

- `Gate3M_graph_is_accepted`: **PASS**
- `Gate3N_exact_replay_is_accepted`: **PASS**
- `2112_heavy_nodes_received_coordinates`: **PASS**
- `1992_fixed_heavy_nodes_were_preserved`: **PASS**
- `120_exact_bridge_coordinates_were_applied`: **PASS**
- `30_graph_paths_match_coordinate_paths`: **PASS**
- `3066_heavy_edges_were_audited`: **PASS**
- `all_BN_bonds_are_within_0p003nm`: **PASS**
- `all_four_atom_bridge_BN_bonds_are_within_0p003nm`: **PASS**
- `critical_angle_minimum_is_at_least70deg`: **PASS**
- `critical_angle_maximum_is_at_most175deg`: **PASS**
- `critical_angle_RMS_deviation_is_at_most30deg`: **PASS**
- `no_nonbonded_heavy_heavy_clashes`: **PASS**
- `aperture_errors_are_within10percent`: **FAIL**
- `outer_radius_errors_are_within15percent`: **PASS**
- `lower_upper_asymmetry_is_within0p010nm`: **PASS**
- `no_H_coordinates_were_generated`: **PASS**

## Decision

- Decision: **R2_SELECTED_FOUR_ATOM_BN_BRIDGE_HEAVY_COORDINATE_EMBEDDING_REQUIRES_REVIEW**
- Failed gates:
  **aperture_errors_are_within10percent**
- Hydrogen coordinates generated: **NO**
- Candidate is final chemistry: **NO**
- Molecular topology generation authorized: **NO**
- Energy minimization authorized: **NO**
- MD authorized: **NO**
- QM authorized: **NO**
- Required next step:
  `REVIEW_R2_SELECTED_FOUR_ATOM_HEAVY_COORDINATE_EMBEDDING_FAILURES`
