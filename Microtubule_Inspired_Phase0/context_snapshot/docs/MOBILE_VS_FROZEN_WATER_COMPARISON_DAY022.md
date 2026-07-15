# Mobile versus Frozen Water Comparison

## Scope

- Mobile trajectory: **201 frames / 100 ps**
- Frozen-solute trajectory: **201 frames / 100 ps**
- Water molecules: **16634**
- Geometric lumen definition: water oxygen inside the instantaneous
  HBN wall radius and between the 1st and 99th percentiles of the
  HBN axial coordinates.

## Lumen water

- Mobile occupancy mean/std:
  1.9055/
  2.2341
- Frozen occupancy mean/std:
  281.9005/
  128.7457
- Mobile versus frozen occupancy difference:
  -99.3241 %

- Mobile density mean/std:
  0.068249/
  0.080076 nm^-3
- Frozen density mean/std:
  10.369915/
  4.736004 nm^-3
- Mobile versus frozen density difference:
  -99.3419 %

## Spatial distributions

- Radial Jensen-Shannon divergence:
  0.01471366
- Axial Jensen-Shannon divergence:
  0.46555245
- Mobile/frozen mean normalized radial position:
  0.525039/
  0.471811

## Temporal behavior

- Mobile occupancy slope:
  0.00815132
  water molecules ps^-1
- Frozen occupancy slope:
  -4.31993202
  water molecules ps^-1

## PYR hydration

- Maximum absolute difference among all PYR/cutoff comparisons:
  15.5134 %

## Screening interpretation

- Analysis status: **COMPLETE**
- Interpretation: **MOBILITY_DEPENDENT_REORGANIZATION_CANDIDATE**
- Screening flags:
  mean lumen occupancy differs by more than 5% | axial Jensen-Shannon divergence exceeds 0.02 | local PYR hydration differs by more than 10%

These screening thresholds are comparative diagnostics and are not
force-field acceptance criteria.

## Snapshot candidates

- Representative matched mobile/frozen snapshot pairs:
  **21**
- Electronic recalculation authorized by this analysis:
  **NO**

The candidate set can be reviewed before recalculating electronic
properties and comparing against the existing time-dependent
solvent-induced site energies under frozen-solute conditions.
