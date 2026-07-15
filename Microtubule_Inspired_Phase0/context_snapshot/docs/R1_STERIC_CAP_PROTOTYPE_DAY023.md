# R1 Steric Cap Prototype

## Scientific role

R1 is a **neutral, frozen steric positive control** designed to test
whether blocking the two axial exits prevents the progressive water
depletion observed in R0.

R1 is not yet a chemically realizable final device architecture.
No physical cap atom type, bonded model, or nonbonded parameter set has
been assigned at this stage.

## R0 reference

- Source:
  `runs/phase1A/day023_confinement_design/01_r0_t0_reference/r0_accepted_t0_hydrated_system.gro`
- Source SHA256:
  `3e2f207361765c7448099591664f17bcecd4e2f53c516520d2ebcfb512028754`
- Initial lumen waters:
  **428**
- Tube axis:
  **(0.00000000,
  0.00000000,
  1.00000000)**
- Robust axial planes:
  **-3.008725/3.009275 nm**
- Wall radius q99:
  **1.199616 nm**
- Provisional accessible radius:
  **0.949111 nm**

## Candidate scan

- Axial offsets:
  **0.20, 0.25, 0.30 nm**
- Triangular-lattice spacings:
  **0.18, 0.20, 0.22 nm**
- Disk radius:
  **1.359616 nm**
- Candidates evaluated:
  **9**
- Candidates passing:
  **9**

- `offset_0.20_spacing_0.18`: PASS; beads/end=211; coverage hole=0.1034 nm; removed waters=91; retained lumen=428/428
- `offset_0.20_spacing_0.20`: PASS; beads/end=163; coverage hole=0.1149 nm; removed waters=83; retained lumen=428/428
- `offset_0.20_spacing_0.22`: PASS; beads/end=139; coverage hole=0.1260 nm; removed waters=90; retained lumen=428/428
- `offset_0.25_spacing_0.18`: PASS; beads/end=211; coverage hole=0.1034 nm; removed waters=106; retained lumen=428/428
- `offset_0.25_spacing_0.20`: PASS; beads/end=163; coverage hole=0.1149 nm; removed waters=99; retained lumen=428/428
- `offset_0.25_spacing_0.22`: PASS; beads/end=139; coverage hole=0.1260 nm; removed waters=99; retained lumen=428/428
- `offset_0.30_spacing_0.18`: PASS; beads/end=211; coverage hole=0.1034 nm; removed waters=126; retained lumen=428/428
- `offset_0.30_spacing_0.20`: PASS; beads/end=163; coverage hole=0.1149 nm; removed waters=117; retained lumen=428/428
- `offset_0.30_spacing_0.22`: PASS; beads/end=139; coverage hole=0.1260 nm; removed waters=120; retained lumen=428/428

## Selected candidate

- Candidate:
  **offset_0.20_spacing_0.20**
- Offset:
  **0.200 nm**
- Lattice spacing:
  **0.200 nm**
- Beads per cap:
  **163**
- Total cap beads:
  **326**
- Maximum planar coverage hole:
  **0.114878 nm**
- Provisional axial steric overlap:
  **0.140000 nm**
- Minimum cap-HBN distance:
  **0.200001 nm**
- Minimum cap-PYR distance:
  **0.959497 nm**

## Water preservation

- Initial waters:
  **16634**
- Removed waters:
  **83**
- Retained waters:
  **16551**
- Initial lumen waters:
  **428**
- Removed lumen waters:
  **0**
- Retained lumen waters:
  **428**
- Retained lumen fraction:
  **1.000000**

Only complete TIP4P/2005 water molecules were removed.

## Derived files

- Cap-only geometry:
  `runs/phase1A/day023_confinement_design/02_r1_steric_cap_prototype/r1_selected_steric_caps_only.gro`
- Geometry-only full system:
  `runs/phase1A/day023_confinement_design/02_r1_steric_cap_prototype/r1_t0_hydrated_with_steric_caps_geometry_only.gro`
- Removed-water audit:
  `runs/phase1A/day023_confinement_design/02_r1_steric_cap_prototype/r1_removed_water_molecules_due_to_cap_overlap.csv`
- Cap coordinate table:
  `runs/phase1A/day023_confinement_design/02_r1_steric_cap_prototype/r1_selected_cap_atom_coordinates.csv`
- Machine-readable definition:
  `runs/phase1A/day023_confinement_design/02_r1_steric_cap_prototype/r1_selected_steric_cap_definition.json`

## Decision

- Geometry prototype accepted: **YES**
- Topology construction authorized: **YES**
- Energy minimization authorized: **NO**
- MD execution authorized: **NO**
- QM execution authorized: **NO**
- Required next step:
  `DEFINE_R1_CAP_NONBONDED_MODEL_AND_BUILD_TOPOLOGY`

Before simulation, the cap model must explicitly define:

1. zero net charge;
2. cap-cap exclusions or zero cap-cap interaction;
3. a water-oxygen steric interaction;
4. interactions with HBN and PYR;
5. frozen coordinate groups;
6. static energy and overlap validation.
