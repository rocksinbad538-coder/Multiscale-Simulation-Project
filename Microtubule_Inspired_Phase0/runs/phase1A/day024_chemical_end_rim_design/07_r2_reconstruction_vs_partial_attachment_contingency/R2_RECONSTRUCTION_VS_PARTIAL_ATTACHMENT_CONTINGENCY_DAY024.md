# R2 Reconstruction versus Partial-Attachment Contingency

## Scope

This gate compares two responses to the negative Gate 3E result:

1. removal of complete graph shells from the BNNT end in an attempt to
   expose a mixed 15 B / 15 N coordination-two rim;
2. partial heteropolar attachment of the validated hexagonal edge seed
   to the n=5, m=2 annulus, with explicit hydrogen passivation of every
   unmatched coordination-two site.

No three-dimensional molecular coordinates, partial charges,
force-field parameters, minimization, MD, or QM calculation were
generated.

## Full-shell parent reconstruction screen

- LOWER depth 0: frontier=30; B/N=30/0; d1/d2/d3=30/0/0; mixed-30-site-rim=False
- LOWER depth 1: frontier=30; B/N=0/30; d1/d2/d3=0/30/0; mixed-30-site-rim=False
- LOWER depth 2: frontier=30; B/N=30/0; d1/d2/d3=30/0/0; mixed-30-site-rim=False
- LOWER depth 3: frontier=30; B/N=0/30; d1/d2/d3=0/30/0; mixed-30-site-rim=False
- LOWER depth 4: frontier=30; B/N=30/0; d1/d2/d3=30/0/0; mixed-30-site-rim=False
- UPPER depth 0: frontier=30; B/N=0/30; d1/d2/d3=30/0/0; mixed-30-site-rim=False
- UPPER depth 1: frontier=30; B/N=30/0; d1/d2/d3=0/30/0; mixed-30-site-rim=False
- UPPER depth 2: frontier=30; B/N=0/30; d1/d2/d3=30/0/0; mixed-30-site-rim=False
- UPPER depth 3: frontier=30; B/N=30/0; d1/d2/d3=0/30/0; mixed-30-site-rim=False
- UPPER depth 4: frontier=30; B/N=0/30; d1/d2/d3=30/0/0; mixed-30-site-rim=False

Simple full-shell reconstruction candidates:

- **0**

## Selected partial-attachment architecture

- Annulus:
  **n=5, m=2**
- Annulus B/N atoms:
  **126 per end**
- Edge-completion seed:
  **30 atoms per end**
- Total added heavy atoms:
  **156 per end**
- Target heavy atoms:
  **145.133 per end**
- Relative heavy-atom error:
  **0.074877**
- Direct heteropolar attachments:
  **15 per end**
- H-passivated seed sites:
  **15 per end**
- H-passivated outer-annulus sites:
  **15 per end**
- H-passivated inner-annulus sites:
  **12 per end**
- Total H passivants:
  **42 per end**

## Final graph audit

- Heavy-graph components:
  **1**
- Heavy-graph bipartite:
  **True**
- Heavy-graph girth:
  **6**
- Four-member heavy cycles:
  **0**
- B/N coordination failures:
  **0**
- H coordination failures:
  **0**
- Nonheteropolar heavy edges:
  **0**

## Geometric proxies

- Outer radius:
  **1.074041 nm**
- Parent-rim target:
  **1.199126 nm**
- Outer-radius relative error:
  **0.104313**
- Inner radius:
  **0.522709 nm**
- Aperture-radius target:
  **0.419703 nm**
- Inner-radius relative error:
  **0.245426**

The inner-radius value is a nucleus-position lattice proxy. The
effective steric aperture after inward H termination must be measured
during the static coordinate-embedding gate.

## Gates

- `Gate3A_parent_audit_is_accepted`: **PASS**
- `Gate3D_hexagonal_edge_completion_seed_is_accepted`: **PASS**
- `Gate3E_direct_30_site_attachment_was_rejected`: **PASS**
- `simple_full_shell_reconstruction_does_not_create_30_site_15B_15N_rim`: **PASS**
- `selected_annulus_is_n5_m2`: **PASS**
- `selected_annulus_has_126_heavy_atoms`: **PASS**
- `selected_annulus_has_30_outer_boundary_sites`: **PASS**
- `selected_annulus_is_connected`: **PASS**
- `selected_annulus_is_bipartite`: **PASS**
- `selected_annulus_has_no_four_member_cycles`: **PASS**
- `selected_annulus_has_only_degree2_or_degree3_atoms`: **PASS**
- `15_alternating_seed_sites_are_attached_per_end`: **PASS**
- `15_remaining_seed_sites_are_H_passivated_per_end`: **PASS**
- `15_noncomplementary_outer_sites_are_H_passivated_per_end`: **PASS**
- `all_inner_boundary_sites_are_H_passivated`: **PASS**
- `all_BN_atoms_reach_total_coordination3`: **PASS**
- `all_H_atoms_have_coordination1`: **PASS**
- `all_heavy_edges_are_heteropolar_BN`: **PASS**
- `combined_heavy_graph_is_connected`: **PASS**
- `combined_heavy_graph_is_bipartite`: **PASS**
- `combined_heavy_graph_has_no_four_member_cycles`: **PASS**
- `combined_heavy_graph_girth_is_at_least6`: **PASS**
- `attachment_mapping_is_angularly_local`: **PASS**
- `attachment_created_cycles_have_length_at_least6`: **PASS**
- `heavy_atom_population_error_is_within15_percent`: **PASS**
- `outer_radius_error_is_within15_percent`: **PASS**
- `inner_radius_proxy_error_is_within30_percent`: **PASS**
- `no_coordinates_were_assigned`: **PASS**
- `no_formal_charges_were_assigned`: **PASS**
- `no_force_field_types_were_assigned`: **PASS**

## Decision

- Decision:
  **R2_PARTIAL_HETEROPOLAR_ANNULUS_ATTACHMENT_AND_COMPLEMENTARY_PASSIVATION_GRAPH_VALIDATED**
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
  `BUILD_AND_VALIDATE_R2_PARTIAL_ATTACHMENT_ANNULUS_STATIC_COORDINATE_EMBEDDING`
