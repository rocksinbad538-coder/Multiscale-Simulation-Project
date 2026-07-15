# R2 Partial-Attachment Annulus Static Coordinate Embedding

## Scope

This stage generated an analytic, non-minimized coordinate embedding
for the Gate 3F graph.

No molecular topology, formal charges, force-field parameters,
minimization, MD, or QM calculation was generated.

## Coordinate inventory

- Parent HBN atoms:
  **1680**
- Added B/N atoms:
  **312**
- Added H atoms:
  **84**
- Total coordinate nodes:
  **2076**
- Parent coordinates modified:
  **NO**

## Attachment geometry

- Seed–annulus mean/minimum/maximum bond:
  **0.310185/
  0.177563/
  0.441245 nm**
- Seed–annulus RMS/max target deviation:
  **0.198283/
  0.296272 nm**
- Lower/upper plane gap:
  **0.000000/
  0.000020 nm**
- Lower/upper annulus-center offset:
  **0.000000/
  0.000000 nm**

## Junction angles

- Minimum/mean/maximum:
  **0.519/
  113.338/
  167.746 degrees**
- RMS deviation from 120 degrees:
  **28.811 degrees**
- Junction degree failures:
  **0**

## Aperture proxy

- Target diameter:
  **0.839406 nm**
- Lower nuclear H-defined diameter:
  **0.807415 nm**
- Upper nuclear H-defined diameter:
  **0.807415 nm**
- Lower/upper relative error:
  **0.038112/
  0.038112**

This is a nucleus-position proxy. It is not a hydrated free-energy
aperture and cannot replace an excluded-volume or PMF calculation.

## Nonbonded contacts

- Heavy–heavy minimum/count below threshold:
  **0.076175/
  24**
- H–heavy minimum/count below threshold:
  **0.047737/
  21**
- H–H minimum/count below threshold:
  **0.118919/
  0**

## Gates

- `Gate3A_parent_audit_is_accepted`: **PASS**
- `Gate3F_graph_design_is_accepted`: **PASS**
- `Gate3F_has_no_failed_gates`: **PASS**
- `all_2076_nodes_received_finite_coordinates`: **PASS**
- `parent_coordinates_are_unchanged`: **PASS**
- `parent_seed_bonds_match_target_within_0p003nm`: **PASS**
- `annulus_internal_bonds_match_target_within_0p002nm`: **PASS**
- `seed_annulus_attachment_RMS_deviation_within_0p020nm`: **FAIL**
- `seed_annulus_attachment_max_deviation_within_0p035nm`: **FAIL**
- `all_XH_bonds_match_provisional_targets_within_0p002nm`: **PASS**
- `junction_nodes_all_have_three_neighbors`: **PASS**
- `junction_angle_minimum_is_at_least_70deg`: **FAIL**
- `junction_angle_maximum_is_at_most_170deg`: **PASS**
- `junction_angle_RMS_deviation_is_at_most_35deg`: **PASS**
- `annulus_center_offset_is_at_most_0p050nm`: **PASS**
- `attachment_plane_gap_is_positive_and_at_most_0p250nm`: **FAIL**
- `annulus_is_planar_to_numerical_precision`: **PASS**
- `nuclear_aperture_error_is_within10percent`: **PASS**
- `outer_radius_error_is_within15percent`: **PASS**
- `no_nonbonded_heavy_heavy_clashes`: **FAIL**
- `no_nonbonded_H_heavy_clashes`: **FAIL**
- `no_nonbonded_H_H_clashes`: **PASS**
- `lower_and_upper_embeddings_are_symmetric`: **FAIL**

## Decision

- Decision:
  **R2_PARTIAL_ATTACHMENT_ANNULUS_STATIC_COORDINATE_EMBEDDING_REQUIRES_CONSTRAINED_GEOMETRIC_OPTIMIZATION**
- Failed gates:
  **seed_annulus_attachment_RMS_deviation_within_0p020nm | seed_annulus_attachment_max_deviation_within_0p035nm | junction_angle_minimum_is_at_least_70deg | attachment_plane_gap_is_positive_and_at_most_0p250nm | no_nonbonded_heavy_heavy_clashes | no_nonbonded_H_heavy_clashes | lower_and_upper_embeddings_are_symmetric**
- Candidate is final chemistry:
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
- Required next step:
  `OPTIMIZE_R2_PARTIAL_ATTACHMENT_STATIC_EMBEDDING_WITH_CONSTRAINED_GEOMETRY`
