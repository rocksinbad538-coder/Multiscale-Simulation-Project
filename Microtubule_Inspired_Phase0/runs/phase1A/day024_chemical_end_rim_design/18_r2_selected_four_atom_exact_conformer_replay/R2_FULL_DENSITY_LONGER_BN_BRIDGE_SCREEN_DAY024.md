# R2 Full-Density Longer BN Bridge Screen

## Scope

This gate screens alternating BN bridges containing four, five or six
heavy atoms while retaining 15 attachments per end.

No accepted coordinates were replaced. No molecular topology, formal
charges, force-field parameters, minimization, MD or QM calculation was
generated.

## Results

- m=4: lower/upper feasible mappings=1/1; feasible=True
- m=5: lower/upper feasible mappings=1/1; feasible=True
- m=6: lower/upper feasible mappings=3/3; feasible=True

## Selection

- Selected bridge atoms per attachment: **4**
- Selected attachments per end: **15**
- Selected minimum angle: **70.89983978777892**
- Selected minimum clearance: **0.15974739665991966**
- Selected graph checks pass: **True**

## Audit gates

- `Gate3I_graph_is_accepted`: **PASS**
- `Gate3K1_trimer_redesign_is_confirmed`: **PASS**
- `bridge_classes_4_5_6_were_screened`: **PASS**
- `60_mappings_per_end_and_class_were_screened`: **PASS**
- `all_pair_metrics_are_finite`: **PASS**

## Decision

- Decision: **R2_FULL_DENSITY_LONGER_BN_BRIDGE_CLASS_IDENTIFIED**
- Failed audit-integrity gates: **NONE**
- Selected graph generation authorized: **YES**
- Coordinate generation authorized: **NO**
- Molecular topology generation authorized: **NO**
- Energy minimization authorized: **NO**
- MD authorized: **NO**
- QM authorized: **NO**
- Required next step: `BUILD_AND_VALIDATE_R2_SELECTED_FULL_DENSITY_LONGER_BN_BRIDGE_GRAPH`
