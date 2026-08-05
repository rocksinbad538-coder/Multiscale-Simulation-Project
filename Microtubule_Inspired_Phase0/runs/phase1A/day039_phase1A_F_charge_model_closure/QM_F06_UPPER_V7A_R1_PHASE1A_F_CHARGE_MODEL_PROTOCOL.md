# QM_F06_UPPER_V7A_R1 — Phase 1A-F working charge model

## Formal status

Phase 1A-F is closed.

The model documented here is adopted as the working charge model for
force-field integration and validation in Phase 1A-G.

It is not yet an adopted production force field.

## Source system

- Source QM model: QM_F06_UPPER_V7A_R1
- Source atoms: 52
- Retained real atoms: 37
- Removed artificial QM boundary caps: 15
- Target total charge: 0 e
- Charge unit: elementary charge

## Adopted fitting method

The selected real-atom charges minimize the ESP residual plus the
lambda-weighted squared deviation from the real-atom RESP Stage 1
charges.

The adopted settings are:

- lambda = 4
- exact neutrality: sum(q) = 0
- local inequality: q[A:UPPER:8:4] >= 0
- active optimum: q[A:UPPER:8:4] = 0 e

## Adopted working vector

The full-grid A12 solution is stored in:

- QM_F06_UPPER_V7A_R1_PHASE1A_F_ADOPTED_WORKING_CHARGES.csv
- QM_F06_UPPER_V7A_R1_PHASE1A_F_ADOPTED_WORKING_CHARGES.dat
- QM_F06_UPPER_V7A_R1_PHASE1A_F_ADOPTED_WORKING_CHARGES.json

Atom ordering follows the retained real-atom order from the original
52-atom transferability table. Reordering is prohibited during
force-field integration.

## Internal validation completed

1. ORCA VPOT source and unit validation.
2. Amber ESP format and round-trip validation.
3. RESP Stage 1 execution and charge audit.
4. Artificial-cap partition and transferability analysis.
5. Real37 constrained-refit feasibility.
6. Lambda-path and lambda=4 selection.
7. Local nonnegative-B constraint review.
8. KKT optimality validation.
9. Deterministic interleaved holdout.
10. Genuine train-only holdout.
11. Six-fold blocked spatial cross-validation.

## Blocked spatial cross-validation

All six blocked folds satisfied the KKT conditions, activated the
nonnegative-B boundary, improved RMSE relative to the unmodified real37
model and produced highly correlated charge vectors.

Spatial extrapolation was nevertheless anisotropic.

Mean validation RMSE:
0.0144114284289 a.u.

Maximum validation RMSE:
0.0213414332078 a.u.

Minimum Pearson correlation:
-0.721256673608

Minimum same-sign fraction:
0.245621099255

Minimum charge-vector correlation against the full-grid candidate:
0.993000500173

Accordingly, the model is adopted as a working transferable charge
model, not as evidence of uniformly accurate regional ESP
extrapolation.

## Phase 1A-G authorization

Authorized activities:

- map the 37 charges into the target topology
- verify atom IDs and atom order
- audit bonded and nonbonded parameter coverage
- verify topology net charge
- perform single-point energy checks
- perform minimization
- perform short controlled validation MD

Still blocked:

- final force-field adoption
- production MD
- Phase 1A closure
- Phase 1B execution
- Phase 2 execution

## Scientific limitation

All current ESP evidence comes from a single geometry and a single
authorized ORCA ESP calculation. Conformational transferability must be
tested during Phase 1A-G and, where needed, using additional QM
snapshots.
