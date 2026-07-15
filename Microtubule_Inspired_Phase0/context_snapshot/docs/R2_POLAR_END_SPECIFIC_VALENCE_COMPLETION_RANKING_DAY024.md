# R2 Polar End-Specific Valence-Completion Candidate Ranking

## Scope

This stage ranks chemical-architecture classes for replacing the
validated neutral steric R2 cap.

No coordinates, molecular topology, force-field parameters, partial
charges, minimization, MD, or QM calculation were generated.

The ranking is a deterministic screening instrument, not a formation
energy, synthetic-yield prediction, or proof of chemical stability.

## Parent terminal requirements

- Degree-1 terminal sites:
  **60**
- Lower end:
  **30 B-terminated sites**
- Upper end:
  **30 N-terminated sites**
- Missing bond incidences per parent terminal site:
  **2**
- Required new parent–junction bonds per end:
  **60**
- Required new parent–junction bonds total:
  **120**

The lower B-terminated end requires an N-complementary first junction
layer. The upper N-terminated end requires a B-complementary first
junction layer. Same-element primary parent-junction bonds are not
authorized by this gate.

## Geometric target

- Aperture diameter:
  **0.839406 nm**
- Aperture radius:
  **0.419703 nm**
- Open-area fraction:
  **0.142928**
- Parent-rim radius:
  **1.199126 nm**
- Required radial occlusion:
  **0.779423 nm**
- Annular area per end:
  **3.963910 nm²**
- Screening estimate:
  **145.133 BN atoms/end**

## Candidate ranking

- Rank 1: `C2_GRADED_HETEROPOLAR_BN_COLLAR_ANNULUS` — 77.70; **ADVANCE_TO_CONNECTIVITY_BLUEPRINT_GATE**
- Rank 2: `C3_RECONSTRUCTED_EDGE_PLUS_GRADED_BN_ANNULUS` — 68.70; **RETAIN_AS_PARENT_RECONSTRUCTION_CONTINGENCY**
- Rank 3: `C1_DIRECT_PLANAR_BN_ANNULUS_DUAL_LINK` — 66.80; **DEFER_DIRECT_SEAM_HIGH_STRAIN**
- Rank 4: `C4_END_SPECIFIC_ORGANIC_MACROCYCLE` — 47.50; **FALLBACK_AFTER_INORGANIC_CANDIDATE_FAILURE**
- Rank 5: `C6_METAL_OR_SPIN_ACTIVE_COORDINATED_CAP` — 40.70; **REJECT_AT_CURRENT_GATE**
- Rank 6: `C0_SMALL_GROUP_PASSIVATION_ONLY` — 35.10; **REJECT_AS_R2_CAP_REPLACEMENT**
- Rank 7: `C5_INWARD_FLEXIBLE_TETHER_CORONA` — 27.00; **DEFER_PORE_NOT_STRUCTURALLY_FIXED**

## Primary hypothesis

`C2_GRADED_HETEROPOLAR_BN_COLLAR_ANNULUS`

This candidate uses a graded, end-specific, heteropolar BN collar
between the accepted tubular parent and the planar BN annulus:

- lower B parent → N-rich first collar layer;
- upper N parent → B-rich first collar layer;
- two additional parent-junction incidences per terminal site;
- 60 new parent-junction bonds per end;
- 120 new parent-junction bonds total;
- target central aperture preserved.

The candidate is not final chemistry. The next gate must determine
whether a graph with these requirements can be constructed without
invalid coordination, same-element primary junctions, disconnected
components, or an impossible ring topology.

## Contingency

`C3_RECONSTRUCTED_EDGE_PLUS_GRADED_BN_ANNULUS`

Parent-rim reconstruction is retained only as a contingency because it
would modify the accepted scaffold and require a separate structural
comparability assessment.

## Gates

- `Gate3A_parent_audit_is_accepted`: **PASS**
- `Gate3A_has_no_failed_gates`: **PASS**
- `explicit_topology_classification_is_accepted`: **PASS**
- `parent_has_60_degree1_terminal_atoms`: **PASS**
- `parent_has_zero_degree2_atoms`: **PASS**
- `parent_has_1620_degree3_atoms`: **PASS**
- `lower_end_is_30B_0N`: **PASS**
- `upper_end_is_0B_30N`: **PASS**
- `terminal_coordination_deficit_is_120`: **PASS**
- `lower_end_requires_60_new_parent_bonds`: **PASS**
- `upper_end_requires_60_new_parent_bonds`: **PASS**
- `target_aperture_is_finite_and_open`: **PASS**
- `open_area_fraction_is_valid`: **PASS**
- `candidate_set_contains_at_least_six_classes`: **PASS**
- `all_candidate_scores_are_finite`: **PASS**
- `primary_candidate_is_highest_ranked`: **PASS**
- `primary_candidate_requires_two_bonds_per_site`: **PASS**
- `primary_candidate_accounts_for_120_parent_bonds`: **PASS**
- `primary_candidate_preserves_target_aperture`: **PASS**
- `primary_candidate_resolves_end_asymmetry`: **PASS**
- `primary_candidate_does_not_require_parent_reconstruction`: **PASS**
- `parent_reconstruction_is_retained_as_contingency`: **PASS**
- `small_group_passivation_is_not_selected_as_cap`: **PASS**
- `metal_spin_active_candidate_is_rejected`: **PASS**

## Decision

- Decision:
  **R2_POLAR_END_SPECIFIC_VALENCE_COMPLETION_CANDIDATES_RANKED**
- Failed gates:
  **NONE**
- Primary candidate:
  **C2_GRADED_HETEROPOLAR_BN_COLLAR_ANNULUS**
- Primary screening score:
  **77.70**
- Primary candidate is final chemistry:
  **NO**
- Connectivity-blueprint gate authorized:
  **YES**
- Explicit geometry generation authorized:
  **NO**
- Topology generation authorized:
  **NO**
- Energy minimization authorized:
  **NO**
- MD authorized:
  **NO**
- QM authorized:
  **NO**
- Required next step:
  `BUILD_AND_VALIDATE_R2_GRADED_HETEROPOLAR_COLLAR_CONNECTIVITY_BLUEPRINT`

## Literature-status limitation

Open BNNT edges and BNNT functionalization provide chemical context,
but no cited work is treated as direct evidence that this exact
graded collar–annulus junction is stable or synthesizable. Its first
test is therefore graph-level valence and connectivity feasibility.
