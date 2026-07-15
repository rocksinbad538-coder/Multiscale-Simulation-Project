# R2 Alternating BN Trimer-Bridge Static Coordinate Embedding

## Scope

This gate embeds all 30 alternating BN trimer bridges and all 174 H
passivants while retaining the accepted parent coordinates.

No topology, formal charges, force-field parameters, minimization, MD,
or QM calculation was generated.

## Coordinate inventory

- Parent atoms: **1680**
- Total heavy atoms: **2082**
- H atoms: **174**
- Total nodes: **2256**
- Explicit bridge conformers: **30**
- Parent coordinates changed:
  **NO**

## Bond geometry

- Maximum B-N deviation:
  **0.001195 nm**
- Maximum bridge B-N deviation:
  **0.000029 nm**
- Maximum X-H deviation:
  **0.000000 nm**
- Maximum conformer-library distance mismatch:
  **0.000134 nm**

## Critical valence angles

- Minimum/mean/maximum:
  **23.689/117.331/171.568 degrees**
- RMS deviation from 120 degrees:
  **15.628 degrees**

## Nonbonded contacts

- Heavy-heavy minimum/clashes:
  **0.059513/6**
- H-heavy minimum/clashes:
  **0.029396/19**
- H-H minimum/clashes:
  **0.060473/0**

## Aperture

- Target diameter:
  **0.839406 nm**
- Lower/upper H-defined diameter:
  **0.816386/
  0.816386 nm**

## Gates

- `Gate3A_parent_audit_is_accepted`: **PASS**
- `Gate3G1_direct_junction_is_rejected`: **PASS**
- `Gate3I_trimer_bridge_graph_is_accepted`: **PASS**
- `Gate3I_has_no_failed_gates`: **PASS**
- `all_2256_nodes_received_finite_coordinates`: **PASS**
- `parent_coordinates_are_unchanged`: **PASS**
- `30_explicit_trimer_conformers_were_embedded`: **PASS**
- `library_endpoint_distance_error_is_within_0p0005nm`: **PASS**
- `all_BN_bonds_are_within_0p003nm_of_target`: **PASS**
- `all_bridge_BN_bonds_are_within_0p003nm_of_target`: **PASS**
- `all_XH_bonds_are_within_0p002nm_of_target`: **PASS**
- `critical_valence_angle_minimum_is_at_least_70deg`: **FAIL**
- `critical_valence_angle_maximum_is_at_most_175deg`: **PASS**
- `critical_valence_angle_RMS_deviation_is_at_most_30deg`: **PASS**
- `annulus_center_offsets_are_within_0p050nm`: **PASS**
- `aperture_errors_are_within10percent`: **PASS**
- `outer_radius_errors_are_within15percent`: **PASS**
- `no_nonbonded_heavy_heavy_clashes`: **FAIL**
- `no_nonbonded_H_heavy_clashes`: **FAIL**
- `no_nonbonded_H_H_clashes`: **PASS**
- `lower_and_upper_embeddings_are_symmetric_within_0p010nm`: **PASS**

## Decision

- Decision:
  **R2_ALTERNATING_BN_TRIMER_BRIDGE_STATIC_COORDINATE_EMBEDDING_REQUIRES_CONFORMER_REFINEMENT**
- Failed gates:
  **critical_valence_angle_minimum_is_at_least_70deg | no_nonbonded_heavy_heavy_clashes | no_nonbonded_H_heavy_clashes**
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
  `REFINE_R2_TRIMER_BRIDGE_CONFORMERS_AND_H_PASSIVANT_ORIENTATIONS`

## Interpretation

This is a deterministic static conformer embedding. Passing the gate
would establish geometric consistency only, not energetic stability,
synthetic feasibility, or force-field validity.
