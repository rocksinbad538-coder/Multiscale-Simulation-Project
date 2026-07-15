# Matched Water Scientific Review

## Scope

- Matched branch interval: **44-144 ps**
- Frames per trajectory: **201**
- Coordinate interval: **0.5 ps**
- Temporal blocks: **10 × 10 ps**
- Block-bootstrap replicates: **50000**

This is a screening-level paired time-series analysis. The ten
temporal blocks reduce frame-level autocorrelation but do not replace
independent trajectory replicas.

## Lumen occupancy

- Mobile mean: 1.905473
- Matched frozen mean: 0.457711
- Absolute mobile-minus-frozen difference:
  1.447761 waters
- Block-bootstrap 95% interval:
  [0.515714,
  2.650000]
- Exact block sign-flip p value:
  0.00585938
- Fraction of blocks with positive mobile-minus-frozen difference:
  0.900000
- Screening effect supported:
  **YES**

The previously reported 316% relative difference is secondary because
the matched-frozen mean is below one water molecule.

## Reentry-event kinetics

- Mobile zero-occupancy fraction:
  0.318408
- Matched-frozen zero-occupancy fraction:
  0.716418
- Mobile/matched-frozen positive episodes:
  20/
  18
- Longest mobile/matched-frozen positive episode:
  18.500/
  5.000 ps

## Spatial distributions

- Radial Jensen-Shannon divergence:
  0.06134384
- Axial Jensen-Shannon divergence:
  0.02234335

These divergences remain descriptive because the distributions contain
few lumen-water observations and no independent trajectory replicas.

## Local pyrene hydration

- Supported PYR/cutoff metrics:
  NONE
- Suggestive PYR/cutoff metrics:
  NONE

## Decision

- Screening interpretation:
  **MOBILITY_ASSOCIATED_TRANSIENT_REHYDRATION_SUPPORTED**
- Publication-level causal claim authorized: **NO**
- Electronic recalculation authorized: **NO**
- Longer mobile production authorized: **NO**
- Authorized next step:
  `REVIEW_21_MATCHED_SNAPSHOT_PAIRS_FOR_LIMITED_MOBILE_QM_PILOT`

The scientifically appropriate interpretation remains transient
rehydration or local hydration reorganization, not stable confined
water retention.
