# Phase 1A ORCA TDDFT Pilot Inputs

## Purpose

This folder contains the first pilot ORCA TDDFT inputs generated from the accepted hydrated Phase 1A baseline.

These inputs are not yet production-level electronic-structure calculations. Their purpose is to validate:

- coordinate extraction,
- charge and multiplicity assumptions,
- ORCA input formatting,
- approximate computational cost,
- suitability of the selected QM model for later TDDFT/site-energy analysis.

## Method

Current pilot method:

- Functional: PBE0
- Basis set: def2-SVP
- Dispersion: D3BJ
- Approximation: RIJCOSX
- SCF: TightSCF
- TDDFT roots: 10
- Charge: 0
- Multiplicity: 1
- Parallel cores: 4

## Input classes

### PYR-only TDDFT

Directory:

`runs/phase1A/day012_electronic_prep/orca_inputs/pyr_only_tddft/`

Files:

- `PYR2_only.inp`
- `PYR3_only.inp`
- `PYR4_only.inp`
- `PYR5_only.inp`

Purpose:

- Isolated pyrene site-energy references.
- Baseline transition dipoles.
- Comparison against hydrated local environments.

### PYR + water0.50 TDDFT

Directory:

`runs/phase1A/day012_electronic_prep/orca_inputs/pyr_water0p50_tddft/`

Files:

- `PYR2_water0p50.inp`
- `PYR3_water0p50.inp`
- `PYR4_water0p50.inp`
- `PYR5_water0p50.inp`

Purpose:

- Local hydration-induced site-energy shifts.
- Local solvent perturbation of TDDFT transitions.
- Comparison with isolated chromophore calculations.

## Important technical note

The water-shell inputs currently treat the selected water molecules quantum mechanically. This is acceptable for pilot tests but may be too expensive or noisy for production.

Recommended next alternatives:

1. full QM PYR-only;
2. PYR + explicit water shell;
3. PYR with water represented as fixed point charges;
4. PYR embedded in h-BN/water electrostatic environment.

## Operational status

The ORCA input-generation workflow is now scripted in:

`scripts/electronic_prep/generate_orca_tddft_inputs.py`

This improves reproducibility and allows the input set to be regenerated from the current QM extracts.
