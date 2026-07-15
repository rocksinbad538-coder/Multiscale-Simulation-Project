# R2 Graded Heteropolar BN Collar Connectivity Blueprint

## Scope

This gate constructs and validates a coordinate-free graph blueprint
for the selected C2 graded heteropolar BN collar–annulus candidate.

No three-dimensional coordinates, molecular topology, bond-order
model, partial charges, force-field parameters, minimization, MD, or
QM calculation were generated.

## Selected ring-population sequence

- Sequence:
  **30-30-30-22-20-12**
- Layers:
  **6**
- Added B/N atoms per end:
  **144**
- Screening estimate:
  **145.133 atoms/end**
- Relative heavy-atom error:
  **0.007806**
- Inner-boundary population:
  **12**
- Inner-boundary spacing proxy:
  **0.219756 nm**
- Parent terminal-site spacing:
  **0.251144 nm**
- Relative spacing error:
  **0.124981**

The selected graph contains exactly 144 added B/N atoms per end,
matching the 144 steric R2 beads per end.

## End-specific composition

### Lower B-terminated end

- First added layer:
  **N**
- Added B/N/H:
  **64/
  80/
  12**
- Parent-to-collar edges:
  **60**
- Same-element heavy edges:
  **0**
- Connected:
  **True**

### Upper N-terminated end

- First added layer:
  **B**
- Added B/N/H:
  **80/
  64/
  12**
- Parent-to-collar edges:
  **60**
- Same-element heavy edges:
  **0**
- Connected:
  **True**

### Combined added structure

- B/N/H:
  **144/
  144/
  24**
- Added heavy atoms:
  **288**
- Parent-to-collar bonds:
  **120**
- Coordination failures:
  **0**
- Same-element heavy edges:
  **0**

## Interfaces

- LOWER `PARENT_TO_L1`: 30 → 30 nodes; 60 edges; maximum angular span 0.033333 turns
- LOWER `L1_TO_L2`: 30 → 30 nodes; 30 edges; maximum angular span 0.000000 turns
- LOWER `L2_TO_L3`: 30 → 30 nodes; 60 edges; maximum angular span 0.033333 turns
- LOWER `L3_TO_L4`: 30 → 22 nodes; 30 edges; maximum angular span 0.027273 turns
- LOWER `L4_TO_L5`: 22 → 20 nodes; 36 edges; maximum angular span 0.063636 turns
- LOWER `L5_TO_L6`: 20 → 12 nodes; 24 edges; maximum angular span 0.066667 turns
- LOWER `INNER_BOUNDARY_TO_H`: 12 → 12 nodes; 12 edges; maximum angular span 0.000000 turns
- UPPER `PARENT_TO_L1`: 30 → 30 nodes; 60 edges; maximum angular span 0.033333 turns
- UPPER `L1_TO_L2`: 30 → 30 nodes; 30 edges; maximum angular span 0.000000 turns
- UPPER `L2_TO_L3`: 30 → 30 nodes; 60 edges; maximum angular span 0.033333 turns
- UPPER `L3_TO_L4`: 30 → 22 nodes; 30 edges; maximum angular span 0.027273 turns
- UPPER `L4_TO_L5`: 22 → 20 nodes; 36 edges; maximum angular span 0.063636 turns
- UPPER `L5_TO_L6`: 20 → 12 nodes; 24 edges; maximum angular span 0.066667 turns
- UPPER `INNER_BOUNDARY_TO_H`: 12 → 12 nodes; 12 edges; maximum angular span 0.000000 turns

## Gates

- `Gate3A_parent_audit_is_accepted`: **PASS**
- `Gate3B_candidate_ranking_is_accepted`: **PASS**
- `primary_candidate_is_C2_graded_heteropolar_BN_collar`: **PASS**
- `selected_ring_population_sequence_is_expected`: **PASS**
- `selected_blueprint_has_six_layers`: **PASS**
- `selected_blueprint_has_144_added_heavy_atoms_per_end`: **PASS**
- `total_added_heavy_atoms_match_288_R2_steric_beads`: **PASS**
- `heavy_atom_estimate_relative_error_is_within_2_percent`: **PASS**
- `inner_boundary_spacing_error_is_within_15_percent`: **PASS**
- `lower_parent_first_layer_is_B_to_N`: **PASS**
- `upper_parent_first_layer_is_N_to_B`: **PASS**
- `parent_to_collar_bond_count_is_60_per_end`: **PASS**
- `parent_to_collar_bond_count_is_120_total`: **PASS**
- `all_parent_and_added_atom_coordination_targets_are_met`: **PASS**
- `all_heavy_atom_edges_are_heteropolar_BN`: **PASS**
- `lower_end_graph_is_connected`: **PASS**
- `upper_end_graph_is_connected`: **PASS**
- `blueprint_has_no_duplicate_edges`: **PASS**
- `blueprint_has_no_self_edges`: **PASS**
- `interface_angular_locality_is_within_threshold`: **PASS**
- `inner_boundary_has_12_H_passivants_per_end`: **PASS**
- `total_inner_passivants_are_24`: **PASS**
- `combined_added_BN_composition_is_balanced`: **PASS**
- `no_formal_charges_were_assigned`: **PASS**
- `no_coordinates_were_assigned`: **PASS**

## Decision

- Decision:
  **R2_GRADED_HETEROPOLAR_BN_COLLAR_CONNECTIVITY_BLUEPRINT_VALIDATED**
- Failed gates:
  **NONE**
- Candidate is final chemistry:
  **NO**
- Static coordinate embedding authorized:
  **YES**
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
  `BUILD_AND_VALIDATE_R2_GRADED_HETEROPOLAR_COLLAR_STATIC_COORDINATE_EMBEDDING`

## Interpretation limitation

The graph demonstrates that the parent valence deficits, end
asymmetry, heavy-atom population, connectivity, and heteropolar
bonding constraints can be satisfied simultaneously at the abstract
topological level. It does not establish that the graph can be
embedded in three dimensions with chemically acceptable B–N and X–H
bond lengths, bond angles, strain, planarity, aperture size, or
energetic stability.
