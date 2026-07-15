# R2 Trimer-Bridge Conformer and H Refinement

## Scope

This gate globally refines the discrete conformers of all 30 BN trimer
bridges and the orientations of all 174 H passivants.

The validated graph and the parent, seed and annulus coordinates were
not changed.

No topology, formal charges, force-field parameters, minimization, MD
or QM calculation was generated.

## Refined inventory

- Coordinate nodes:
  **2256**
- Bridge conformers:
  **30**
- H orientations:
  **174**
- Parent coordinates unchanged:
  **True**
- Seed and annulus coordinates unchanged:
  **True**

## Bond geometry

- Maximum B-N deviation:
  **0.001195 nm**
- Maximum bridge B-N deviation:
  **0.000199 nm**
- Maximum X-H deviation:
  **0.000000 nm**
- Maximum library endpoint mismatch:
  **0.000922 nm**

## Critical angles

- Minimum/mean/maximum:
  **57.220/
  117.222/
  154.795 degrees**
- RMS deviation from 120 degrees:
  **12.206 degrees**

## Nonbonded contacts

- Heavy-heavy minimum/clashes:
  **0.137562/
  0**
- H-heavy minimum/clashes:
  **0.090154/
  0**
- H-H minimum/clashes:
  **0.081619/
  0**

## Aperture

- Target diameter:
  **0.839406 nm**
- Lower/upper refined diameter:
  **0.816386/
  0.816386 nm**

## Gates

- `Gate3A_parent_audit_is_accepted`: **PASS**
- `Gate3I_trimer_bridge_graph_is_accepted`: **PASS**
- `Gate3I_has_no_failed_gates`: **PASS**
- `Gate3J_requires_conformer_refinement`: **PASS**
- `all_2256_nodes_received_finite_coordinates`: **PASS**
- `parent_coordinates_are_unchanged`: **PASS**
- `seed_and_annulus_coordinates_are_unchanged`: **PASS**
- `30_bridge_conformers_were_globally_refined`: **PASS**
- `174_H_orientations_were_globally_refined`: **PASS**
- `library_endpoint_errors_are_within_0p0005nm`: **FAIL**
- `all_BN_bonds_are_within_0p003nm`: **PASS**
- `all_bridge_BN_bonds_are_within_0p003nm`: **PASS**
- `all_XH_bonds_are_within_0p002nm`: **PASS**
- `critical_angle_minimum_is_at_least70deg`: **FAIL**
- `critical_angle_maximum_is_at_most175deg`: **PASS**
- `critical_angle_RMS_deviation_is_at_most30deg`: **PASS**
- `no_nonbonded_heavy_heavy_clashes`: **PASS**
- `no_nonbonded_H_heavy_clashes`: **PASS**
- `no_nonbonded_H_H_clashes`: **PASS**
- `aperture_errors_are_within10percent`: **PASS**
- `outer_radius_errors_are_within15percent`: **PASS**
- `lower_upper_asymmetry_is_within0p010nm`: **PASS**

## Decision

- Decision:
  **R2_ALTERNATING_BN_TRIMER_BRIDGE_REFINEMENT_REQUIRES_BRIDGE_TOPOLOGY_REDESIGN**
- Failed gates:
  **library_endpoint_errors_are_within_0p0005nm | critical_angle_minimum_is_at_least70deg**
- Heavy geometry pass:
  **False**
- Hydrogen geometry pass:
  **True**
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
  `EVALUATE_R2_LONGER_OR_TOPOLOGICALLY_REVISED_BRIDGE_ARCHITECTURE`
