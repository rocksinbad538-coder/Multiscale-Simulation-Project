# Phase 1A Electronic-Structure Preparation: QM Extracts Manifest

Source structure:

- Accepted hydrated Phase 1A baseline
- File: `runs/phase1A/accepted/hybrid_hbnBonded_kang2000_improperGeo100_hydrated_baseline_10ps100K_postmin/phase1A_hydrated_baseline.gro`

## Extract classes

### 1. PYR-only monomers

Directory:

`runs/phase1A/day012_electronic_prep/qm_extracts_pyr_only/`

Files:

- `PYR2_only.xyz`
- `PYR3_only.xyz`
- `PYR4_only.xyz`
- `PYR5_only.xyz`

Purpose:

- Isolated chromophore reference calculations.
- Baseline DFT/TDDFT site energies.
- Reference transition dipoles.
- Comparison against hydrated local environments.

### 2. PYR + local water shell

Directory:

`runs/phase1A/day012_electronic_prep/qm_extracts/`

Files:

- `PYR2_water0p50.xyz`: 102 atoms, 19 waters
- `PYR3_water0p50.xyz`: 98 atoms, 18 waters
- `PYR4_water0p50.xyz`: 102 atoms, 19 waters
- `PYR5_water0p50.xyz`: 118 atoms, 23 waters

Purpose:

- Local solvatochromic shift estimation.
- Water-induced site-energy perturbation.
- Local electrostatic/environmental sensitivity tests.

### 3. PYR–PYR pairs

Directory:

`runs/phase1A/day012_electronic_prep/qm_extracts_pyr_pairs/`

Files:

- `PYR2_PYR3_pair.xyz`
- `PYR2_PYR4_pair.xyz`
- `PYR2_PYR5_pair.xyz`
- `PYR3_PYR4_pair.xyz`
- `PYR3_PYR5_pair.xyz`
- `PYR4_PYR5_pair.xyz`

Purpose:

- Pairwise electronic-structure inputs.
- Excitonic-coupling preparation.
- Geometry-dependent chromophore-pair comparison.

## Current interpretation

The pyrene array is not a close pi-stacked aggregate. It forms an alternating quasi-orthogonal arrangement:

- Nearest-neighbor pairs are spatially closer but nearly orthogonal.
- Same-orientation pairs are nearly parallel but farther apart.

This implies that direct excitonic coupling is expected to be pair-dependent and likely weak to moderate rather than strongly aggregate-like.

## Next recommended electronic-structure workflow

1. Run single-point DFT/TDDFT on PYR-only monomers.
2. Run TDDFT on PYR + water0.50 local environments.
3. Compare excitation energies and transition dipoles.
4. Use pair geometries to estimate or compute pairwise couplings.
5. Assemble preliminary site-energy/coupling table for an effective excitonic Hamiltonian.
