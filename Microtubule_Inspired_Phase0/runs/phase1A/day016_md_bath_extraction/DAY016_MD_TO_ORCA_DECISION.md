# Day016 MD-to-ORCA Extraction Decision

## Current trajectory

Accepted trajectory:

- `runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.xtc`
- 68320 atoms
- 201 frames
- 0–100 ps
- 0.5 ps sampling interval

## Mapping

MD labels are mapped to the exciton model as:

- MD_PYR1 -> PYR2
- MD_PYR2 -> PYR3
- MD_PYR3 -> PYR4
- MD_PYR4 -> PYR5

## Geometry validation

Nearest-neighbor chromophore distances in the final accepted GRO are consistent with the excitonic chain:

- PYR2-PYR3 ≈ 26.9 Å
- PYR3-PYR4 ≈ 26.9 Å
- PYR4-PYR5 ≈ 27.6 Å

## Extracted frames

21 frames were extracted every 5 ps from 0 to 100 ps.

Each frame contains:

- HBN: 1680 atoms
- PYR: 104 atoms
- SOL: 66536 atoms

## Local QM clusters

For each frame, local clusters were generated with a 5 Å water shell around each monomer and nearest-neighbor dimer.

Typical sizes:

- monomers: ~40–44 waters, ~185–202 atoms
- dimers: ~82–84 waters, ~378–388 atoms

## Technical decision

The explicit-water clusters are useful for geometry and electrostatic screening, but are too large for systematic TDDFT production.

Production extraction should use:

- PYR-only QM region
- nearby waters represented as electrostatic point-charge embedding
- optional limited explicit-water pilots for calibration

## Important limitation

The accepted 100 ps trajectory is sampled every 0.5 ps. This is too coarse to resolve sub-100 fs bath correlation times.

Therefore, this trajectory is appropriate for:

- MD-to-ORCA pipeline validation
- structural ensemble sampling
- slow disorder estimates
- preliminary site-energy/coupling distributions

It is not sufficient for final bath correlation-time extraction in the ENAQT window identified on Day015.

## Next required physical step

Run or obtain a short high-time-resolution MD trajectory, ideally:

- 10–20 ps total length
- 5–10 fs output interval
- same accepted hydrated structural model

This is required to estimate MD-derived autocorrelation times in the 20–100 fs range.
