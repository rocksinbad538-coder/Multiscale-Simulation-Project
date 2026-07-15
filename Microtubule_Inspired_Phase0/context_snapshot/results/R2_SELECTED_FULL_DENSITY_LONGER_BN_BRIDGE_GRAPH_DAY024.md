# R2 Selected Full-Density Four-Atom BN Bridge Graph

## Scope

This gate replaces the previous three-atom BN bridges with the shortest
full-density longer bridge class selected by Gate 3L.

No coordinates, molecular topology, formal charges, force-field
parameters, minimization, MD or QM calculation were generated.

## Architecture

- Bridge atoms per path:
  **4**
- Bonds per path:
  **5**
- Attachments per end:
  **15**
- Total paths:
  **30**
- Total bridge atoms:
  **120**

### LOWER

- Mapping: **M4:LOWER:P0:O-1:R14**
- Bridge sequence: **B-N-B-N**
- Paths/bridge atoms: **15/60**
- Bridge B/N: **30/30**
- H seed/outer/inner/bridge/other/total: **15/15/12/60/0/102**
- Bridge-containing cycles: **16–16**
### UPPER

- Mapping: **M4:UPPER:P0:O-1:R0**
- Bridge sequence: **N-B-N-B**
- Paths/bridge atoms: **15/60**
- Bridge B/N: **30/30**
- H seed/outer/inner/bridge/other/total: **15/15/12/60/0/102**
- Bridge-containing cycles: **16–16**

## Complete graph

- Heavy atoms:
  **2112**
- H atoms:
  **204**
- Total nodes:
  **2316**
- Heavy/H edges:
  **3066/204**
- Heavy/H degree failures:
  **0/
  0**
- Nonheteropolar heavy edges:
  **0**
- Full/heavy connected components:
  **1/1**
- Bipartite:
  **True**
- Four-member cycles:
  **0**
- Bridge cycle range:
  **16–16**

## Gates

- `Gate3I_source_graph_is_accepted`: **PASS**
- `Gate3L_longer_bridge_screen_is_accepted`: **PASS**
- `selected_bridge_class_has_four_atoms`: **PASS**
- `selected_architecture_has_15_attachments_per_end`: **PASS**
- `30_old_trimer_paths_were_replaced`: **PASS**
- `30_new_four_atom_bridge_paths_were_built`: **PASS**
- `15_bridge_paths_were_built_per_end`: **PASS**
- `60_bridge_atoms_were_added_per_end`: **PASS**
- `120_bridge_atoms_were_added_total`: **PASS**
- `2112_heavy_atoms_are_present`: **PASS**
- `204_H_atoms_are_present`: **PASS**
- `2316_total_nodes_are_present`: **PASS**
- `3066_heavy_edges_are_present`: **PASS**
- `204_H_edges_are_present`: **PASS**
- `all_heavy_atoms_have_total_coordination3`: **PASS**
- `all_H_atoms_have_coordination1`: **PASS**
- `no_heavy_atom_has_heavy_degree_above3`: **PASS**
- `all_heavy_edges_are_heteropolar_BN`: **PASS**
- `full_graph_is_connected`: **PASS**
- `heavy_graph_is_connected`: **PASS**
- `heavy_graph_is_bipartite`: **PASS**
- `heavy_graph_contains_no_four_member_cycles`: **PASS**
- `every_bridge_path_participates_in_a_cycle`: **PASS**
- `all_bridge_containing_cycles_have_length16`: **PASS**
- `no_coordinates_were_assigned`: **PASS**
- `no_formal_charges_were_assigned`: **PASS**
- `no_force_field_types_were_assigned`: **PASS**

## Decision

- Decision:
  **R2_SELECTED_FULL_DENSITY_FOUR_ATOM_BN_BRIDGE_GRAPH_VALIDATED**
- Failed gates:
  **NONE**
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
  `BUILD_AND_VALIDATE_R2_SELECTED_FOUR_ATOM_BN_BRIDGE_STATIC_COORDINATE_EMBEDDING`
