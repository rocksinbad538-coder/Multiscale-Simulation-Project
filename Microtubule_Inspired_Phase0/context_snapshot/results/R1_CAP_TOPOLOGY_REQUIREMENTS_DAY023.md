# R1 Cap Topology Requirements

## Purpose

This audit identifies the exact force-field conventions required to
add the two R1 steric caps without introducing uncontrolled attraction,
charge, or cap–solute interactions.

No topology was modified and no GROMACS preprocessing or dynamics was
performed.

## Baseline nonbonded convention

- `nbfunc`: **1**
- `comb-rule`: **2**
- parameter semantics: **SIGMA_EPSILON**
- `gen-pairs`: **yes**
- active atom types: **7**
- explicit `[ nonbond_params ]` entries:
  **0**

## Baseline molecule inventory

- `HBN`: 1 × 1680 atoms = 1680
- `PYR`: 4 × 26 atoms = 104
- `SOL`: 16634 × 4 atoms = 66536 — identified water

- Baseline atom count:
  **68320/68320**

## Water model identification

- Molecule type:
  **SOL**
- Sites per molecule:
  **4**
- Oxygen atom:
  **OW**
- Oxygen atom type:
  **OW**
- Oxygen charge:
  **0.00000000 e**

## R1 coordinate ordering

- R1 GRO atoms:
  **68314**
- Retained waters:
  **16551**
- Lower-cap atoms:
  **163**
- Upper-cap atoms:
  **163**
- Water four-site ordering:
  **PASS**

## Proposed R1 molecule inventory

- `HBN`: 1 × 1680 = 1680 (RETAIN)
- `PYR`: 4 × 26 = 104 (RETAIN)
- `SOL`: 16551 × 4 = 66204 (REPLACE_WATER_COUNT)
- `CAPL`: 1 × 163 = 163 (ADD_FROZEN_CAP)
- `CAPU`: 1 × 163 = 163 (ADD_FROZEN_CAP)

- Proposed topology atom count:
  **68314**
- R1 coordinate atom count:
  **68314**

## Proposed cap contract

The provisional cap model will contain:

- atom type: `CAP`;
- charge: **0 e**;
- mass: **12.011 u**;
- base nonbonded parameters: **zero**;
- cap–cap interaction: **zero**;
- cap–HBN interaction: **zero**;
- cap–PYR interaction: **zero**;
- cap–water H/M interaction: **zero**;
- cap–water O interaction:
  **requires an alternative implementation**;
- two molecule types: `CAPL` and `CAPU`;
- all cap coordinates frozen during the screening.

No C12 value is authorized yet. Candidate values are recorded in:

`runs/phase1A/day023_confinement_design/03_r1_topology_model/r1_cap_water_repulsion_calibration.csv`

They must be tested with static energy and force scans before energy
minimization or MD.

## Gates

- `baseline_atom_count_is_68320`: **PASS**
- `R1_GRO_atom_count_is_68314`: **PASS**
- `proposed_topology_atom_count_matches_R1_GRO`: **PASS**
- `water_molecule_has_four_sites`: **PASS**
- `water_oxygen_type_is_unique`: **PASS**
- `R1_water_chunks_are_consistent`: **PASS**
- `lower_cap_has_expected_resname`: **PASS**
- `upper_cap_has_expected_resname`: **PASS**
- `no_existing_cap_atomtype_collision`: **PASS**
- `no_existing_cap_nonbond_override_collision`: **PASS**

## Decision

- Decision: **TOPOLOGY_CONTRACT_VALID_BUT_PURE_R12_REQUIRES_ALTERNATIVE_IMPLEMENTATION**
- Topology contract valid:
  **YES**
- Standard pure-r12 override feasible:
  **NO**
- Energy minimization authorized: **NO**
- MD execution authorized: **NO**
- Required next step:
  `SELECT_TABULATED_OR_VALIDATED_ALTERNATIVE_CAP_INTERACTION`
