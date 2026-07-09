# R2 Parent Rim and Chemical-Constraint Audit

## Scope

This audit characterizes the accepted R2 parent BNNT and translates
the validated steric-cap geometry into chemical and atomistic design
constraints.

No coordinates, topology, trajectory, checkpoint, or accepted
simulation result were modified.

No minimization, molecular dynamics, DFT, TDDFT, or other quantum
calculation was executed.

## Parent system

- Total atoms:
  **68332**
- HBN atoms:
  **1680**
- B/N atoms:
  **840/840**
- Pyrene atoms:
  **104**
- Water atoms:
  **66260**
- Cap beads:
  **288**

The accepted HBN screening topology assigns zero charge to all
1680 B/N atoms. This is retained as a record of the
screening model and is not interpreted as final chemical
electrostatics.

## Geometry-derived BN connectivity

- Bonds:
  **2460**
- Expected bonds:
  **2460**
- Mean/median BN distance:
  **0.144973/0.145000 nm**
- Minimum/maximum BN distance:
  **0.143778/0.146034 nm**
- Degree-2 terminal atoms:
  **60**
- Degree-3 interior atoms:
  **1620**
- Coordination anomalies:
  **0**

The accepted topology and the independent geometric graph identify
60 degree-1 terminal sites and 1620 degree-3 interior sites. The lower
end contains 30 B-only terminal atoms, whereas the upper end contains
30 N-only terminal atoms. These are strongly undercoordinated polar
termini in the current parent model. Frozen-coordinate stability does
not establish their chemical stability.

## Terminal rims

### Lower end

- Atoms:
  **30**
- B/N:
  **30/0**
- Mean radius:
  **1.199126 nm**
- Axial standard deviation:
  **0.000000 nm**
- Dominant element/purity:
  **B/
  1.000000**
- Element-alternation fraction:
  **0.000000**
- Cap-plane offset:
  **0.200000 nm**

### Upper end

- Atoms:
  **30**
- B/N:
  **0/30**
- Mean radius:
  **1.199126 nm**
- Axial standard deviation:
  **0.000000 nm**
- Dominant element/purity:
  **N/
  1.000000**
- Element-alternation fraction:
  **0.000000**
- Cap-plane offset:
  **0.200000 nm**

## Validated R2 target

- Effective aperture diameter:
  **0.839406 nm**
- Effective aperture radius:
  **0.419703 nm**
- Open-area fraction:
  **0.142928**
- Mean parent-rim radius:
  **1.199126 nm**
- Required radial occlusion:
  **0.779423 nm**
- Annular area per end:
  **3.963910 nm²**
- Estimated one-layer h-BN annular population:
  **145.133 atoms/end**
- Validated steric beads:
  **144 beads/end**
- Validated minimum CAP–OW separation:
  **0.166949 nm**

Simple H/OH/NHx termination cannot be promoted directly as the R2
chemical replacement. It cannot reproduce the validated aperture, and
the parent atoms are degree-1 polar termini rather than conventional
degree-2 edge sites. The chemical design must therefore resolve an
aggregate parent-side coordination deficit of
**120 missing neighbor
incidences**, while also providing an inward radial closure of
**0.779423 nm**.

## Preliminary candidate classes

- `C0_SIMPLE_EDGE_PASSIVATION_ONLY` — priority 0: **REJECT_AS_STANDALONE_CAP_REPLACEMENT**
- `C1_LINKED_BN_ANNULAR_NANOFLAKE` — priority 1: **ADVANCE_TO_EXPLICIT_GEOMETRY_AND_JUNCTION_DESIGN**
- `C2_RIGID_ORGANIC_OR_HYBRID_MACROCYCLE` — priority 2: **ADVANCE_AS_FALLBACK_CANDIDATE_CLASS**
- `C3_INWARD_TETHERED_FUNCTIONAL_CORONA` — priority 3: **DEFER_UNLESS_RIGID_ANNULAR_CANDIDATES_FAIL**
- `C4_METAL_OR_SPIN_ACTIVE_CAP` — priority 4: **REJECT_AT_CURRENT_GATE**

These are candidate classes, not selected final chemistries.

The leading geometric analogue remains a separately defined annular
nanoflake or rigid macrocycle positioned at the validated cap plane.
However, the lower B-terminated rim and upper N-terminated rim require
distinct junction chemistry or an explicitly justified compensation
strategy. A symmetric identical linker assignment and a direct
unstrained 90-degree seamless sp2 junction are not assumed.

## Gates

- `R2_architecture_selection_is_valid`: **PASS**
- `explicit_topology_terminal_classification_is_valid`: **PASS**
- `system_has_68332_atoms`: **PASS**
- `HBN_has_1680_atoms`: **PASS**
- `HBN_atoms_are_first_and_contiguous`: **PASS**
- `pyrene_atom_count_is_104`: **PASS**
- `water_atom_count_is_66260`: **PASS**
- `cap_atom_count_is_288`: **PASS**
- `cap_split_is_144_per_end`: **PASS**
- `TPR_dump_return_code_zero`: **PASS**
- `TPR_HBN_atom_count_is_1680`: **PASS**
- `HBN_contains_840_B_and_840_N`: **PASS**
- `HBN_contains_only_B_and_N`: **PASS**
- `current_HBN_screening_charges_are_zero`: **PASS**
- `geometry_bond_count_is_expected`: **PASS**
- `geometry_has_60_degree1_terminal_atoms`: **PASS**
- `geometry_has_zero_degree2_atoms`: **PASS**
- `geometry_has_1620_degree3_interior_atoms`: **PASS**
- `geometry_has_no_coordination_anomalies`: **PASS**
- `geometry_BN_bond_lengths_are_plausible`: **PASS**
- `lower_end_has_60_atoms`: **PASS**
- `upper_end_has_60_atoms`: **PASS**
- `lower_end_is_30B_and_0N`: **PASS**
- `upper_end_is_0B_and_30N`: **PASS**
- `lower_end_is_planar`: **PASS**
- `upper_end_is_planar`: **PASS**
- `terminal_radii_are_symmetric`: **PASS**
- `lower_end_is_elementally_pure_B`: **PASS**
- `upper_end_is_elementally_pure_N`: **PASS**
- `terminal_to_terminal_bond_count_is_zero`: **PASS**
- `terminal_to_interior_bond_count_is_60`: **PASS**
- `validated_aperture_is_finite_and_open`: **PASS**
- `validated_open_area_fraction_is_valid`: **PASS**
- `required_radial_occlusion_is_positive`: **PASS**
- `monolayer_annulus_population_is_comparable_to_R2_beads`: **PASS**

## Decision

- Decision:
  **R2_PARENT_RIM_AND_CHEMICAL_CONSTRAINTS_AUDITED**
- Failed gates:
  **NONE**
- Static candidate-design work authorized:
  **YES**
- Explicit geometry generation authorized:
  **NO**
- Energy minimization authorized:
  **NO**
- New MD authorized:
  **NO**
- Mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `DEFINE_AND_RANK_R2_POLAR_END_SPECIFIC_VALENCE_COMPLETION_CANDIDATES`
