# R2 Aperture Target Semantics Resolution

## Resolved definition

The value **0.839406210 nm**
is an effective water-accessible aperture obtained from the R2 neutral
steric-control model at the 5 kBT interaction threshold.

It is not a direct heavy-atom nuclear diameter.

## Heavy embedding

- Lower/upper heavy nuclear aperture:
  **1.045415170/1.045415170 nm**
- Heavy atoms:
  **2112**
- Exact bridge atoms:
  **120**
- Heavy-heavy clashes:
  **0**
- Minimum critical angle:
  **70.899839788 degrees**
- Maximum B-N deviation:
  **0.001195394 nm**

## Interpretation

The heavy-only aperture is retained as a structural descriptor. It is
not compared directly against the effective 5 kBT target.

The next static geometric proxy is the aperture defined by the inner
hydrogen nuclei after the 204 passivants are placed. Final validation
of the effective water-accessible aperture remains deferred until an
explicit nonbonded model is authorized.

## Gates

- `Gate3O_failed_only_the_non_equivalent_aperture_comparison`: **PASS**
- `Gate3O1_confirmed_all_nuclear_aperture_metrics_exceed_effective_target`: **PASS**
- `heavy_embedding_contains_2112_heavy_nodes`: **PASS**
- `heavy_embedding_contains_120_exact_bridge_nodes`: **PASS**
- `heavy_embedding_contains_3066_heavy_edges`: **PASS**
- `heavy_embedding_has_zero_clashes`: **PASS**
- `heavy_embedding_minimum_angle_is_at_least70deg`: **PASS**
- `heavy_embedding_BN_deviation_is_at_most0p003nm`: **PASS**
- `fixed_coordinates_are_unchanged`: **PASS**
- `lower_upper_heavy_geometry_is_symmetric`: **PASS**
- `target_is_explicitly_classified_as_effective_5kBT_aperture`: **PASS**
- `heavy_nuclear_aperture_is_not_relabelled_as_effective_5kBT_aperture`: **PASS**
- `204_H_coordinates_are_still_pending`: **PASS**

## Decision

- Decision:
  **R2_HEAVY_COORDINATE_EMBEDDING_VALIDATED_APERTURE_FUNCTIONAL_GATE_DEFERRED**
- Failed gates:
  **NONE**
- Heavy coordinate embedding validated:
  **YES**
- Hydrogen coordinate generation authorized:
  **YES**
- Molecular topology generation authorized:
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
  `GENERATE_AND_VALIDATE_R2_SELECTED_FOUR_ATOM_BN_BRIDGE_HYDROGEN_COORDINATES`
