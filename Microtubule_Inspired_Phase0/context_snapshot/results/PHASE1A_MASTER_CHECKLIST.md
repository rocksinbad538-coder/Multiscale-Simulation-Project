# Phase 1A Master Checklist

Date: 2026-06-19

## Scope reference

Phase 1A was initiated to move from coarse Phase 0 screening to a chemically valid atomistic route including:

- real chromophores
- explicit water
- DFT/TDDFT parameterization
- excitonic Hamiltonian construction

## Checklist

| Item | Status | Evidence |
|---|---|---|
| Chemically explicit h-BN scaffold | COMPLETE | `structures/phase1A/accepted/PHASE1A_ACCEPTED_MANIFEST.md` |
| Pyrene chromophore model | COMPLETE | `structures/phase1A/accepted/pyrene/` |
| Dry h-BN + 4 pyrene hybrid | COMPLETE | `structures/phase1A/accepted/hybrid_dry/` |
| GROMACS dry/minimized hybrid | COMPLETE | `runs/phase1A/accepted/hybrid_dry_gromacs_min/` |
| TIP4P/2005 water model | COMPLETE | `parameters/phase1A/accepted/water_tip4p2005/` |
| Hydrated atomistic model | COMPLETE | `runs/phase1A/accepted/hybrid_hbnBonded_kang2000_improperGeo100_hydrated_baseline_10ps100K_postmin/` |
| Scaffold mechanics validation | COMPLETE | `runs/phase1A/accepted/hybrid_hbnBonded_kang2000_improperGeo100_validated_100ps_100K/` |
| Confined-water retention | COMPLETE | `internal_water_timeseries.csv`; final internal OW = 423 |
| Hydration spatial profile | COMPLETE | `water_radial_profile_postmin.csv`, `water_axial_profile_postmin.csv` |
| Pyrene hydration shell audit | COMPLETE | `pyrene_hydration_audit.csv` |
| Pyrene pair geometry | COMPLETE | `runs/phase1A/day012_electronic_prep/excitonic_geometry/` |
| Excitonic dimer prioritization | COMPLETE | `runs/phase1A/day012_exciton_network/` |
| QM monomer extracts | COMPLETE | `qm_extracts_pyr_only/` |
| QM hydrated pyrene extracts | COMPLETE | `qm_extracts/` |
| QM pyrene-pair extracts | COMPLETE | `qm_extracts_pyr_pairs/` |
| ORCA TDDFT input generation | COMPLETE | `orca_inputs/`; audited 8/8 valid inputs |
| First TDDFT calculations | PENDING | ORCA not yet installed locally |
| Site energies | PENDING | Requires TDDFT outputs |
| Transition dipoles | PENDING | Requires TDDFT outputs |
| Pairwise excitonic couplings | PENDING | Requires TDDFT / coupling model |
| Excitonic Hamiltonian | PENDING | Requires site energies + couplings |
| Open quantum dynamics | NOT STARTED | Downstream of Hamiltonian |
| Optical / EM response | NOT STARTED | Downstream of electronic parameters |
| THz / microwave response | NOT STARTED | Downstream of dipole/electronic response |
| Spin module | NOT STARTED / CONDITIONAL | No spin-active center yet defined |

## Current technical conclusion

Phase 1A structural and hydration baseline is effectively closed.

The project is now ready to begin the electronic-structure stage, starting with local ORCA installation and pilot TDDFT calculations.

## Immediate next step

Install ORCA locally on the Mac and run the first pilot:

`PYR2_only.inp`

Then parse:

- normal termination
- excitation energies
- oscillator strengths
- transition information
