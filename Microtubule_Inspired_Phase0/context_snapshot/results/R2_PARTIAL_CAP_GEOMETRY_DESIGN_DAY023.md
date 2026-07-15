# R2 Partial-Cap Geometry Design

## Scope

R2 introduces symmetric central axial apertures into the validated R1
steric-cap lattice.

R2 remains a neutral frozen steric screening architecture. This gate
does not claim chemical realizability and does not authorize molecular
dynamics.

## Design basis

- Authoritative hydrated source:
  `runs/phase1A/day023_confinement_design/01_r0_t0_reference/r0_accepted_t0_hydrated_system.gro`
- Validated R1 cap lattice:
  `runs/phase1A/day023_confinement_design/02_r1_steric_cap_prototype/r1_selected_steric_caps_only.gro`
- Initial lumen occupancy:
  **428 waters**
- R1 cap beads:
  **163 per end**
- CAP–OW 5 kBT distance:
  **0.170 nm**
- Initial-overlap cutoff:
  **0.220 nm**
- Target effective aperture radius:
  **0.300 nm**

## Candidate scan

- `remove_r0p20`: remove radius=0.200 nm; caps=160/160; effective aperture radius=0.030219 nm; open area=0.001014; retained lumen waters=428/428; status=FAIL
- `remove_r0p25`: remove radius=0.250 nm; caps=156/156; effective aperture radius=0.176182 nm; open area=0.034458; retained lumen waters=428/428; status=FAIL
- `remove_r0p30`: remove radius=0.300 nm; caps=156/156; effective aperture radius=0.176182 nm; open area=0.034458; retained lumen waters=428/428; status=FAIL
- `remove_r0p35`: remove radius=0.350 nm; caps=150/150; effective aperture radius=0.229656 nm; open area=0.058549; retained lumen waters=428/428; status=PASS
- `remove_r0p40`: remove radius=0.400 nm; caps=148/148; effective aperture radius=0.230041 nm; open area=0.058746; retained lumen waters=428/428; status=PASS
- `remove_r0p45`: remove radius=0.450 nm; caps=144/144; effective aperture radius=0.358819 nm; open area=0.142928; retained lumen waters=428/428; status=PASS
- `remove_r0p50`: remove radius=0.500 nm; caps=144/144; effective aperture radius=0.358819 nm; open area=0.142928; retained lumen waters=428/428; status=PASS
- `remove_r0p55`: remove radius=0.550 nm; caps=132/132; effective aperture radius=0.429593 nm; open area=0.204871; retained lumen waters=428/428; status=FAIL
- `remove_r0p60`: remove radius=0.600 nm; caps=128/128; effective aperture radius=0.430207 nm; open area=0.205457; retained lumen waters=428/428; status=FAIL
- `remove_r0p65`: remove radius=0.650 nm; caps=126/126; effective aperture radius=0.522344 nm; open area=0.302887; retained lumen waters=428/428; status=FAIL

## Selected R2 candidate

- Candidate:
  **remove_r0p45**
- Nominal cap-bead removal radius:
  **0.450000 nm**
- Lower/upper cap beads:
  **144/
  144**
- Removed cap beads:
  **38**
- Effective 5 kBT aperture radius/diameter:
  **0.358819/
  0.717638 nm**
- Conservative aperture radius/diameter:
  **0.308819/
  0.617638 nm**
- Open-area fraction:
  **0.142928**
- Removed/retained water molecules:
  **69/
  16565**
- Retained lumen waters:
  **428/
  428**
- Lumen-water retention fraction:
  **1.000000**
- Minimum CAP–HBN distance:
  **0.200002 nm**
- Minimum CAP–PYR distance:
  **0.959393 nm**
- Minimum CAP–OW distance:
  **0.220189 nm**
- R2 atom count:
  **68332**

## Static gates

- `R1_positive_control_is_validated`: **PASS**
- `R1_authorized_R2_static_design`: **PASS**
- `R0_atom_count_is_68320`: **PASS**
- `R1_cap_lattice_has_326_beads`: **PASS**
- `R0_initial_lumen_occupancy_is_428`: **PASS**
- `at_least_one_candidate_passed`: **PASS**
- `selected_candidate_passed_all_static_gates`: **PASS**
- `selected_caps_are_symmetric`: **PASS**
- `selected_system_atom_count_is_consistent`: **PASS**

## Decision

- Decision:
  **R2_PARTIAL_CAP_GEOMETRY_STATIC_GATE_PASSED**
- Failed gates:
  **NONE**
- Topology generation authorized:
  **YES**
- Energy minimization authorized:
  **NO**
- Molecular dynamics authorized:
  **NO**
- Required next step:
  `BUILD_R2_TOPOLOGY_AND_RUN_STATIC_CAP_WATER_SCAN`

The selected aperture is a geometry-screening target. Actual water
exchange and retention must be established through a subsequent
validated topology, static interaction scan, water-only minimization,
and short frozen-solute trajectory.
