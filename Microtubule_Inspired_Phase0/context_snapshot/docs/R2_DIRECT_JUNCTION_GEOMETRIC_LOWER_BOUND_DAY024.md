# R2 Direct-Junction Geometric Lower-Bound Audit

## Scope

This stage determines whether any cyclic mapping between the
15 alternating seed sites and the 15 complementary annulus sites can
support a direct B-N junction with moderate deformation.

No coordinates were replaced. No molecular topology, force-field
parameters, minimization, MD, or QM calculation was generated.

## Search

- Fits per end:
  **300**
- Total fits:
  **600**
- Mapping variables:
  seed parity, circumferential orientation, discrete rotation and
  transformation chirality.
- Transformation models:
  rigid, similarity and affine.

### LOWER

- Fits screened: **300**
- Passing local-preoptimization fits: **0**
- Best rigid RMS/max deviation: **0.198197/0.296148 nm**
- Best similarity RMS/max deviation: **0.199928/0.307035 nm**
- Best similarity principal strain: **0.068369**
- Best affine RMS/max deviation: **0.199928/0.307033 nm**
- Best affine principal strain/anisotropy: **0.068395/1.000049**
- Local constrained optimization feasible: **False**
### UPPER

- Fits screened: **300**
- Passing local-preoptimization fits: **0**
- Best rigid RMS/max deviation: **0.198197/0.296148 nm**
- Best similarity RMS/max deviation: **0.199928/0.307035 nm**
- Best similarity principal strain: **0.068369**
- Best affine RMS/max deviation: **0.199928/0.307033 nm**
- Best affine principal strain/anisotropy: **0.068395/1.000049**
- Local constrained optimization feasible: **False**

## Acceptance thresholds

- Attachment RMS deviation:
  **≤ 0.020 nm**
- Attachment maximum deviation:
  **≤ 0.035 nm**
- Annulus internal-bond maximum deviation:
  **≤ 0.010 nm**
- Maximum principal strain:
  **≤ 0.050**
- Maximum affine anisotropy:
  **≤ 1.100**
- Annulus-center offset:
  **≤ 0.050 nm**
- Axial gap:
  **positive and ≤ 0.250 nm**

## Audit gates

- `Gate3F_graph_design_is_accepted`: **PASS**
- `Gate3G_embedding_has_expected_review_decision`: **PASS**
- `600_mapping_and_transform_fits_were_screened`: **PASS**
- `300_fits_were_screened_per_end`: **PASS**
- `all_required_fit_modes_are_present`: **PASS**
- `all_core_fit_metrics_are_finite`: **PASS**

## Decision

- Decision:
  **R2_PARTIAL_ATTACHMENT_DIRECT_BN_JUNCTION_GEOMETRIC_LOWER_BOUND_FAILED**
- Failed audit-integrity gates:
  **NONE**
- Local constrained optimization authorized:
  **NO**
- Current direct B-N attachment graph retained:
  **NO**
- Coordinate update authorized:
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
  `REDESIGN_R2_ANNULUS_JUNCTION_WITH_EXPLICIT_BRIDGING_LINKER_OR_REVISED_ATTACHMENT_TOPOLOGY`

## Interpretation

A negative result means that the current graph cannot be repaired by a
local geometric relaxation without requiring excessive annulus strain,
anisotropy, junction displacement, or nonlocal bond deformation.

In that case, the direct seed-to-annulus B-N edges must be replaced by
an explicit bridging junction or by a different attachment topology.
