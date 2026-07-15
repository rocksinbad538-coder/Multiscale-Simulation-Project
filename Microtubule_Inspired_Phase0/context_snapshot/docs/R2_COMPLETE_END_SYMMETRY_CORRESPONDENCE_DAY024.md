# R2 Complete End-Symmetry Correspondence

## Complete correspondence

- Heavy pairs: **216**
- Seed-heavy pairs: **30**
- H pairs: **102**
- Seed-H pairs: **15**
- Inner-H pairs: **12**

## Transform comparison

- ALL_END_HEAVY / PROPER_ROTATION_ONLY: det=1.000000; RMSD=3.009299227340e-01 nm; max=1.042387923012e+00 nm
- ALL_END_HEAVY / ORTHOGONAL_REFLECTION_ALLOWED: det=-1.000000; RMSD=2.662866243847e-01 nm; max=1.041179193167e+00 nm
- RIGID_RIM_HEAVY / PROPER_ROTATION_ONLY: det=1.000000; RMSD=2.290270549693e-01 nm; max=4.488905673313e-01 nm
- RIGID_RIM_HEAVY / ORTHOGONAL_REFLECTION_ALLOWED: det=-1.000000; RMSD=1.815919844246e-01 nm; max=3.455694579138e-01 nm
- ANNULUS_HEAVY / PROPER_ROTATION_ONLY: det=1.000000; RMSD=5.688217534143e-15 nm; max=6.466036496704e-15 nm
- ANNULUS_HEAVY / ORTHOGONAL_REFLECTION_ALLOWED: det=-1.000000; RMSD=2.666883902359e-15 nm; max=4.070144838902e-15 nm
- INNER_BOUNDARY_HEAVY / PROPER_ROTATION_ONLY: det=1.000000; RMSD=1.905793032954e-15 nm; max=2.035072419451e-15 nm
- INNER_BOUNDARY_HEAVY / ORTHOGONAL_REFLECTION_ALLOWED: det=-1.000000; RMSD=1.905793032954e-15 nm; max=2.035072419451e-15 nm
- SEED_HEAVY / PROPER_ROTATION_ONLY: det=1.000000; RMSD=6.858855589371e-04 nm; max=1.242621274189e-03 nm
- SEED_HEAVY / ORTHOGONAL_REFLECTION_ALLOWED: det=-1.000000; RMSD=6.580852468392e-04 nm; max=1.121191285701e-03 nm
- FOUR_ATOM_BRIDGE_HEAVY / PROPER_ROTATION_ONLY: det=1.000000; RMSD=3.190158074630e-01 nm; max=8.318395985537e-01 nm
- FOUR_ATOM_BRIDGE_HEAVY / ORTHOGONAL_REFLECTION_ALLOWED: det=-1.000000; RMSD=2.920967229912e-01 nm; max=8.319616047124e-01 nm

## Selected rigid-rim transformation

- Subset: **INNER_BOUNDARY_HEAVY**
- Type: **ORTHOGONAL_REFLECTION_ALLOWED**
- Determinant:
  **-1.000000000000**
- RMSD:
  **1.905793032954e-15 nm**
- Maximum deviation:
  **2.035072419451e-15 nm**

## Existing inner-H asymmetry

- RMSD:
  **2.160444102433e-01 nm**
- Maximum:
  **2.181886234714e-01 nm**

## Gates

- `216_lower_heavy_nodes_are_paired`: **PASS**
- `216_upper_heavy_nodes_are_used_once`: **PASS**
- `30_seed_pairs_are_resolved_by_circumferential_index`: **PASS**
- `102_lower_H_nodes_are_paired`: **PASS**
- `102_upper_H_nodes_are_used_once`: **PASS**
- `15_seed_H_pairs_are_resolved_through_seed_correspondence`: **PASS**
- `12_inner_H_pairs_are_resolved`: **PASS**
- `selected_rigid_rim_transform_is_exact`: **PASS**
- `coordinates_are_not_modified`: **PASS**
- `no_topology_minimization_MD_or_QM_is_generated`: **PASS**

## Decision

- Decision: **R2_COMPLETE_END_SYMMETRY_CORRESPONDENCE_VALIDATED**
- Failed gates:
  **NONE**
- Coordinates modified: **NO**
- Molecular topology generated: **NO**
- Energy minimization performed: **NO**
- MD performed: **NO**
- QM performed: **NO**
- Required next step:
  `REFINE_R2_HYDROGEN_ORIENTATIONS_WITH_VALIDATED_END_SYMMETRY_TRANSFORM`
