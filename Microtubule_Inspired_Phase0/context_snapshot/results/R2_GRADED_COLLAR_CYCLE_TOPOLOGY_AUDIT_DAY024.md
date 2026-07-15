# R2 Graded-Collar Cycle-Topology Audit

## Scope

This stage audits the heavy-atom graph produced by Gate 3C.

No coordinates, formal bond orders, charges, force-field parameters,
minimization, MD, or QM calculation were generated.

## Lower end

- Heavy nodes/edges:
  **174/240**
- Connected components:
  **1**
- Bipartite:
  **True**
- Cycle rank:
  **67**
- Girth:
  **4**
- Cycles of length 3/4/5/6/7/8:
  **0/
  39/
  0/
  77/
  0/
  32**
- Added BN nodes participating in four-member cycles:
  **98**

## Upper end

- Heavy nodes/edges:
  **174/240**
- Connected components:
  **1**
- Bipartite:
  **True**
- Cycle rank:
  **67**
- Girth:
  **4**
- Cycles of length 3/4/5/6/7/8:
  **0/
  39/
  0/
  77/
  0/
  32**
- Added BN nodes participating in four-member cycles:
  **98**

## Gates

- `LOWER_heavy_graph_is_connected`: **PASS**
- `LOWER_heavy_graph_is_bipartite`: **PASS**
- `LOWER_has_no_self_edges`: **PASS**
- `LOWER_has_no_duplicate_edges`: **PASS**
- `LOWER_parent_blueprint_degree_is_two`: **PASS**
- `LOWER_added_BN_degree_is_three`: **FAIL**
- `LOWER_four_cycle_enumerators_agree`: **PASS**
- `LOWER_contains_no_triangles`: **PASS**
- `LOWER_contains_no_four_member_cycles`: **FAIL**
- `LOWER_contains_no_five_member_cycles`: **PASS**
- `LOWER_contains_no_odd_cycles_up_to_length10`: **PASS**
- `LOWER_girth_is_at_least_six`: **FAIL**
- `LOWER_contains_hexagonal_cycles`: **PASS**
- `LOWER_added_BN_nodes_avoid_four_cycles`: **FAIL**
- `LOWER_parent_collar_edges_avoid_four_cycles`: **FAIL**
- `UPPER_heavy_graph_is_connected`: **PASS**
- `UPPER_heavy_graph_is_bipartite`: **PASS**
- `UPPER_has_no_self_edges`: **PASS**
- `UPPER_has_no_duplicate_edges`: **PASS**
- `UPPER_parent_blueprint_degree_is_two`: **PASS**
- `UPPER_added_BN_degree_is_three`: **FAIL**
- `UPPER_four_cycle_enumerators_agree`: **PASS**
- `UPPER_contains_no_triangles`: **PASS**
- `UPPER_contains_no_four_member_cycles`: **FAIL**
- `UPPER_contains_no_five_member_cycles`: **PASS**
- `UPPER_contains_no_odd_cycles_up_to_length10`: **PASS**
- `UPPER_girth_is_at_least_six`: **FAIL**
- `UPPER_contains_hexagonal_cycles`: **PASS**
- `UPPER_added_BN_nodes_avoid_four_cycles`: **FAIL**
- `UPPER_parent_collar_edges_avoid_four_cycles`: **FAIL**
- `Gate3C_blueprint_was_previously_accepted`: **PASS**
- `Gate3C_had_no_failed_upstream_gates`: **PASS**
- `lower_and_upper_cycle_statistics_are_symmetric`: **PASS**

## Decision

- Decision:
  **R2_GRADED_HETEROPOLAR_BN_COLLAR_REQUIRES_HEXAGONAL_GRAPH_REDESIGN**
- Failed gates:
  **LOWER_added_BN_degree_is_three | LOWER_contains_no_four_member_cycles | LOWER_girth_is_at_least_six | LOWER_added_BN_nodes_avoid_four_cycles | LOWER_parent_collar_edges_avoid_four_cycles | UPPER_added_BN_degree_is_three | UPPER_contains_no_four_member_cycles | UPPER_girth_is_at_least_six | UPPER_added_BN_nodes_avoid_four_cycles | UPPER_parent_collar_edges_avoid_four_cycles**
- Previous coordinate-embedding authorization superseded:
  **YES**
- Static coordinate embedding authorized:
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
  `REDESIGN_R2_COLLAR_GRAPH_USING_HEXAGONAL_LATTICE_TEMPLATE_AND_REAUDIT`

## Interpretation

Gate 3C established an abstract trivalent, connected and heteropolar
graph. This audit determines whether that graph also has a
low-strain hexagonal-network cycle topology.

The presence of four-member B-N-B-N cycles blocks direct promotion to
the coordinate-embedding stage. In that event, the graph must be
rebuilt from an explicit hexagonal-lattice or nanotube-junction
template rather than from unconstrained bipartite interface flows.
