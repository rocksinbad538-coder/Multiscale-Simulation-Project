# Day023 R1 Design Input Audit

## Repository state

- Project root: `.`
- Git root: `/Users/alejandro/projects/Multiscale-Simulation-Project`
- Branch: `main`
- HEAD: `b3af4af0568279876523a8316bda56b6056b09b5`

## Known authoritative inputs

- Validated topology:
  **PASS**
  `parameters/phase1A/accepted/hybrid_hbnBonded_kang2000_improperGeo100_validated/hbn_pyrene_4_hydratable_gap45_pyr5shift_clean032_hbnBonded_kang2000_improperGeo100.top`
- Stage02 starting coordinates:
  **PASS**
  `runs/phase1A/day021_mobile_restraint_protocol/execution/01_em_k10000/01_em_k10000.gro`
- Stage02 TPR:
  **PASS**
- Stage08 mobile TPR/XTC:
  **PASS**
- Matched frozen-control TPR/XTC:
  **PASS**
- Index file:
  **PASS**

## Accepted frozen R0 directory

- Directory:
  `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute`
- Directory status:
  **PASS**
- GRO files: **1**
- GRO files with 68320 atoms:
  **1**
- TPR files: **1**
- XTC files: **1**
- CPT files: **1**

### Accepted R0 files

- `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/audit_contacts_nvt_100ps_frozenSolute.csv` — analysis_or_report
- `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/audit_pyrene_hbn_nvt_100ps_frozenSolute_NOPBC.csv` — analysis_or_report
- `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.cpt` — checkpoint
- `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.gro` — coordinate_candidate; atoms=68320
- `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.log` — supporting_file
- `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.mdp` — md_parameters
- `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.tpr` — run_input
- `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.xtc` — trajectory

## Accepted trajectory metadata

- `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.xtc`: return code=0; atoms=68320; frames=201; interval=0.5 ps; duration=100.0 ps

## Existing reusable workflow candidates

- `scripts/phase1A/analyze_day020_confined_water_profile_guided_regions.py`
- `scripts/phase1A/analyze_day020_confined_water_regions.py`
- `scripts/phase1A/analyze_day022_matched_mobile_vs_frozen_water.py`
- `scripts/phase1A/analyze_day022_mobile_vs_frozen_water.py`
- `scripts/phase1A/audit_day017_pyrene_structure_and_state_identity.py`
- `scripts/phase1A/audit_day018_pyrene_structure_and_state_identity.py`
- `scripts/phase1A/audit_day020_confined_water_region_sensitivity.py`
- `scripts/phase1A/audit_day020_hbn_axial_architecture.py`
- `scripts/phase1A/audit_day020_md_confined_water_inputs.py`
- `scripts/phase1A/audit_day021_accepted_hydrated_topology.py`
- `scripts/phase1A/audit_day021_mobile_bonded_topology_transition.py`
- `scripts/phase1A/audit_day022_water_depletion_provenance.py`
- `scripts/phase1A/compare_day021_stage05_stage06_hbn_improper_persistence.py`
- `scripts/phase1A/diagnose_day021_hbn_stage_transition.py`
- `scripts/phase1A/diagnose_day021_stage05_hbn_improper_phase.py`
- `scripts/phase1A/diagnose_day021_stage06_hbn_improper_phase.py`
- `scripts/phase1A/diagnose_day022_stage07_hbn_improper_phase.py`
- `scripts/phase1A/diagnose_day022_stage08_hbn_improper_phase.py`
- `scripts/phase1A/finalize_day020_confined_water_axial_radial_density.py`
- `scripts/phase1A/inventory_day021_hydrated_topology_candidates.py`
- `scripts/phase1A/repair_day020_confined_water_density_dat_parser.py`
- `scripts/phase1A/repair_day020_confined_water_density_grid_geometry.py`
- `scripts/phase1A/review_day022_matched_water_effects.py`
- `scripts/phase1A/run_day020_confined_water_axial_radial_density.py`
- `scripts/phase1A/summarize_day022_water_depletion_kinetics.py`
- `scripts/phase1A/verify_day021_original_frozen_topology_identity.py`
- `scripts/phase1A/verify_day021_working_frozen_topology_semantic_identity.py`

## Cap-generation status

- Existing cap-related scripts:
  **0**
- Authoritative hydrated R1 starting state resolved:
  **YES**

## Decision

- R0 remains the validated open-tube reference.
- No existing R0 file is modified by this audit.
- No MD or QM calculation was executed.
- Next required step:
  **EXTRACT_AND_AUDIT_ACCEPTED_R0_T0_HYDRATED_STATE**

R1 must be generated from the earliest hydrated state of the accepted
frozen-solute R0 trajectory, not from the depleted mobile branch state.
