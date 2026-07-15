# HBN Explicit Topology and Terminal Coordination Classification

## Purpose

This stage determines whether the 60 degree-1 atoms identified from
geometry are also present in the accepted explicit HBN topology.

No coordinates or topology were modified. No minimization, MD, or QM
calculation was executed.

## Accepted topology

- Selected topology:
  `parameters/phase1A/accepted/hybrid_hbnBonded_kang2000_improperGeo100_validated/hbn_pyrene_4_hydratable_gap45_pyr5shift_clean032_hbnBonded_kang2000_improperGeo100.top`
- Local files in include tree:
  **4**
- HBN atom records:
  **1680**
- Explicit HBN connections:
  **2460**

## Explicit versus geometric connectivity

- Geometric connections:
  **2460**
- Explicit-only connections:
  **0**
- Geometry-only connections:
  **0**
- Bond sets identical:
  **YES**

## Explicit coordination distribution

- Degree 0:
  **0**
- Degree 1:
  **60**
- Degree 2:
  **0**
- Degree 3:
  **1620**
- Degree 4 or greater:
  **0**

The accepted scaffold therefore contains 60 singly coordinated
terminal sites and 1620 three-coordinate interior sites. The prior
assumption of 120 two-coordinate terminal sites is not applicable to
this structure.

## Terminal-site distribution

### Lower end

- Terminal atoms:
  **30**
- B/N:
  **30/0**
- Axial mean/std:
  **-3.008714/
  0.000000 nm**
- Mean radius:
  **1.199126 nm**

### Upper end

- Terminal atoms:
  **30**
- B/N:
  **0/30**
- Axial mean/std:
  **3.009286/
  0.000000 nm**
- Mean radius:
  **1.199126 nm**

### Total terminal composition

- B/N:
  **30/30**

## Bond distances

- Mean/median:
  **0.144973/
  0.145000 nm**
- Minimum/maximum:
  **0.143778/
  0.146034 nm**

## Gates

- `accepted_topology_contains_1680_HBN_atoms`: **PASS**
- `accepted_topology_contains_2460_HBN_bonds`: **PASS**
- `geometry_contains_2460_HBN_bonds`: **PASS**
- `explicit_and_geometry_bond_sets_are_identical`: **PASS**
- `explicit_and_geometry_degree_arrays_are_identical`: **PASS**
- `explicit_topology_has_60_degree1_atoms`: **PASS**
- `explicit_topology_has_zero_degree2_atoms`: **PASS**
- `explicit_topology_has_1620_degree3_atoms`: **PASS**
- `explicit_topology_has_no_degree0_atoms`: **PASS**
- `explicit_topology_has_no_degree4plus_atoms`: **PASS**
- `lower_end_has_30_terminal_atoms`: **PASS**
- `upper_end_has_30_terminal_atoms`: **PASS**
- `terminal_population_contains_30B_and_30N`: **PASS**
- `lower_terminal_cluster_is_at_axial_minimum`: **PASS**
- `upper_terminal_cluster_is_at_axial_maximum`: **PASS**
- `terminal_end_radii_are_symmetric`: **PASS**
- `all_explicit_connectivity_is_BN`: **PASS**
- `explicit_bond_lengths_are_chemically_plausible`: **PASS**

## Decision

- Decision:
  **HBN_EXPLICIT_TOPOLOGY_AND_TERMINAL_COORDINATION_CLASSIFIED**
- Failed gates:
  **NONE**
- Parent-rim auditor repair authorized:
  **YES**
- Explicit geometry generation authorized:
  **NO**
- Energy minimization authorized:
  **NO**
- MD authorized:
  **NO**
- QM authorized:
  **NO**
- Required next step:
  `REPAIR_R2_PARENT_RIM_AUDITOR_FOR_EXPLICIT_DEGREE1_TERMINI`
