# R2 Aperture Metric Consistency Audit

## Scope

This gate investigates the isolated aperture failure reported by Gate 3O.
No coordinates are modified.

## Coordinate consistency

- Maximum fixed-coordinate difference between Gate 3K and Gate 3O:
  **0.000000000000e+00 nm**

## Compared aperture metrics

- LOWER / TWO_TIMES_MINIMUM_RADIAL_DISTANCE: 1.045415170 nm; relative error=0.245423; pass=False
- LOWER / MINIMUM_OPPOSITE_PAIR_DISTANCE: 1.045415170 nm; relative error=0.245423; pass=False
- LOWER / MEAN_OPPOSITE_PAIR_DISTANCE: 1.045415170 nm; relative error=0.245423; pass=False
- LOWER / MINIMUM_PROJECTED_WIDTH: 1.004403541 nm; relative error=0.196565; pass=False
- UPPER / TWO_TIMES_MINIMUM_RADIAL_DISTANCE: 1.045415170 nm; relative error=0.245423; pass=False
- UPPER / MINIMUM_OPPOSITE_PAIR_DISTANCE: 1.045415170 nm; relative error=0.245423; pass=False
- UPPER / MEAN_OPPOSITE_PAIR_DISTANCE: 1.045415170 nm; relative error=0.245423; pass=False
- UPPER / MINIMUM_PROJECTED_WIDTH: 1.004403541 nm; relative error=0.196565; pass=False

## Preferred metric

- Preferred metric:
  **NONE**
- Lower value:
  ** nm**
- Upper value:
  ** nm**

## Gates

- `Gate3O_has_expected_aperture_only_review_decision`: **PASS**
- `fixed_annulus_coordinates_are_unchanged_from_Gate3K`: **PASS**
- `Gate3O_radial_metric_fails_both_ends`: **PASS**
- `at_least_one_geometric_free_aperture_metric_passes_both_ends`: **FAIL**
- `all_aperture_metrics_are_lower_upper_symmetric`: **PASS**
- `heavy_geometry_other_than_aperture_remains_accepted`: **PASS**

## Decision

- Decision: **R2_APERTURE_METRIC_AUDIT_REQUIRES_FURTHER_REVIEW**
- Failed gates:
  **at_least_one_geometric_free_aperture_metric_passes_both_ends**
- Heavy coordinate embedding retained:
  **NO**
- Hydrogen coordinate generation authorized:
  **NO**
- Molecular topology generation authorized:
  **NO**
- Energy minimization authorized:
  **NO**
- MD authorized:
  **NO**
- QM authorized:
  **NO**
- Required next step:
  `REVIEW_R2_APERTURE_METRIC_DEFINITION_AND_TARGET`
