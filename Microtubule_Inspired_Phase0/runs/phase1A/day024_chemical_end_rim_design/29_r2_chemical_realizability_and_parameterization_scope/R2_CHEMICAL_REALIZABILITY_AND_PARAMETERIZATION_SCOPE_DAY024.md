# R2 Chemical Realizability and Parameterization Scope

## Static chemical inventory

- Total nodes: **2316**
- Heavy/H: **2112/204**
- B/N: **1056/1056**
- Graph edges: **3270**
- B-N/B-H/N-H bonds:
  **3066/102/102**

## Bond geometry

- Maximum B-N deviation:
  **0.001195394 nm**
- Maximum B-H deviation:
  **8.881784197001e-16 nm**
- Maximum N-H deviation:
  **8.604228440845e-16 nm**

## Valence and graph chemistry

- Heavy valence failures:
  **0**
- H valence failures:
  **0**
- Nonheteropolar heavy-heavy edges:
  **0**
- Invalid H edges:
  **0**

## Chemical environments

- Unique local environments:
  **40**
- Parameterization-critical centers:
  **468**

## Parameterization scope

- **PARENT_HBN_BULK_LIKE**: 1620 nodes; QM priority=MEDIUM; novel parameterization=POSSIBLY_NO_FOR_BULK_INTERIOR.
- **PARENT_AND_ANNULUS_ATTACHMENTS**: 60 nodes; QM priority=HIGH; novel parameterization=LIKELY.
- **ANNULUS_AND_SEED_EDGE**: 312 nodes; QM priority=HIGH; novel parameterization=LIKELY.
- **FOUR_ATOM_BRIDGE**: 120 nodes; QM priority=HIGHEST; novel parameterization=YES.
- **HYDROGEN_PASSIVANTS**: 204 nodes; QM priority=HIGH; novel parameterization=LIKELY.

## Gates

- `Gate3M_graph_is_accepted`: **PASS**
- `Gate3P2_refined_coordinates_are_accepted`: **PASS**
- `graph_contains_2316_nodes`: **PASS**
- `coordinates_exist_for_all_2316_nodes`: **PASS**
- `graph_contains_2112_heavy_and_204_H_nodes`: **PASS**
- `graph_has_no_missing_edge_nodes`: **PASS**
- `graph_has_no_self_edges`: **PASS**
- `graph_has_no_duplicate_edges`: **PASS**
- `graph_is_connected`: **PASS**
- `all_heavy_atoms_are_three_coordinate`: **PASS**
- `all_H_atoms_are_one_coordinate`: **PASS**
- `all_heavy_heavy_edges_are_BN`: **PASS**
- `all_H_edges_are_BH_or_NH`: **PASS**
- `102_BH_and_102_NH_bonds_are_present`: **PASS**
- `BN_bond_deviation_is_at_most0p003nm`: **PASS**
- `BH_and_NH_deviations_are_at_most0p002nm`: **PASS**
- `local_chemical_environments_are_fully_enumerated`: **PASS**
- `parameterization_critical_centers_are_explicitly_identified`: **PASS**
- `force_field_coverage_is_not_assumed`: **PASS**
- `no_topology_charges_parameters_minimization_MD_or_QM_generated`: **PASS**

## Decision

- Decision:
  **R2_STATIC_CHEMICAL_REALIZABILITY_VALIDATED_PARAMETERIZATION_SCOPE_DEFINED**
- Failed gates:
  **NONE**
- Static chemical graph realizable:
  **YES**
- Existing force-field coverage established:
  **NO**
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
- Force-field and QM-reference audit authorized:
  **YES**
- Required next step:
  `AUDIT_BN_H_FORCE_FIELD_COVERAGE_AND_DEFINE_QM_REFERENCE_FRAGMENT_SET`
