# R2 Inner-H Reflected-Direction Diagnostic

## Scope

The validated inner-boundary reflection is applied to X-H directions,
not directly to hydrogen coordinates. Each rebuilt hydrogen uses the
B-H or N-H length required by its receiving heavy atom.

## Scenarios

- LOWER_DRIVES_UPPER: aperture L/U/asymmetry=1.008734672/1.002561666/0.006173005 nm; H-heavy clashes=0; H-H clashes=0; angle violations=0; pass=True
- UPPER_DRIVES_LOWER: aperture L/U/asymmetry=1.028648811/1.031124689/0.002475879 nm; H-heavy clashes=0; H-H clashes=0; angle violations=0; pass=True

## Decision

- Decision: **R2_INNER_H_REFLECTED_DIRECTION_REFINEMENT_PATH_IDENTIFIED**
- Passing scenarios: **2**
- Selected scenario: **UPPER_DRIVES_LOWER**
- Coordinates modified: **NO**
- Molecular topology generated: **NO**
- Energy minimization performed: **NO**
- MD performed: **NO**
- QM performed: **NO**
- Required next step:
  `APPLY_AND_VALIDATE_SELECTED_INNER_H_REFLECTED_DIRECTION_SCENARIO`
