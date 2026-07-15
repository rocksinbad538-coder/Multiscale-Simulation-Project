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
  0.4577/
  0.8342
- Mobile versus frozen occupancy difference:
  316.3043 %

- Mobile density mean/std:
  0.068249/
  0.080076 nm^-3
- Frozen density mean/std:
  0.016829/
  0.030673 nm^-3
- Mobile versus frozen density difference:
  305.5465 %

## Spatial distributions

- Radial Jensen-Shannon divergence:
  0.06134384
- Axial Jensen-Shannon divergence:
  0.02234335
- Mobile/frozen mean normalized radial position:
  0.525040/
  0.537610

## Temporal behavior

- Mobile occupancy slope:
  0.00815132
  water molecules ps^-1
- Frozen occupancy slope:
  0.00627457
  water molecules ps^-1

## PYR hydration

- Maximum absolute difference among all PYR/cutoff comparisons:
  32.8102 %

## Screening interpretation

- Analysis status: **COMPLETE**
- Interpretation: **MOBILITY_DEPENDENT_REORGANIZATION_CANDIDATE**
- Screening flags:
  mean lumen occupancy differs by more than 5% | radial Jensen-Shannon divergence exceeds 0.02 | axial Jensen-Shannon divergence exceeds 0.02 | local PYR hydration differs by more than 10%

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
