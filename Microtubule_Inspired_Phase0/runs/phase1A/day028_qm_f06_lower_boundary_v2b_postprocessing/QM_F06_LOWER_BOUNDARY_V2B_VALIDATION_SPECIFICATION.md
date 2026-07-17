# QM_F06 LOWER Boundary V2-B Validation Specification — Day028

## Purpose

Determine whether joint relaxation of all hydrogen atoms and the restored
real boundary heavy atoms resolves the residual artificial-boundary
contact while preserving the validated B–N–B–N bridge/scaffold core.

## Required execution gates

- ORCA return code: 0
- Normal ORCA termination: YES
- SCF convergence: YES
- Geometry convergence: YES
- Final optimized XYZ generated: YES

## Constraint-integrity gates

- Expected fixed atoms: 10
- Expected mobile atoms: 18
- Maximum displacement of fixed atoms: <= 1.0e-6 Å
- Constraint indices must match the V2-B atom-role map exactly

## Chemical-topology gates

- Fragment graph remains connected
- Intended bonds: 27
- Degree/valence failures: 0
- Bond-range failures: 0
- Artificial-cap bond failures: 0
- Unintended covalent contacts: 0

## Boundary-artifact gates

The residual V2-A contact:

`HCAP:LOWER:03 — A:LOWER:13:-3`

must either:

1. rise above the hard-contact threshold of 0.70 of the vdW-radius sum; or
2. be explicitly reclassified only if chemical and geometric evidence
   demonstrates that it is nonblocking.

Required final conditions:

- Hard bridge–cap contacts: 0
- Unintended cap-induced covalent contacts: 0
- No new blocking artificial-cap contacts

## Downstream authorization

Only if all structural and chemical gates pass:

- V2-B geometry acceptance: AUTHORIZED
- Final electronic-property input preparation: AUTHORIZED
- ESP/RESP execution: still requires a separately defined protocol
- Force-field parameter adoption: NOT AUTHORIZED by this gate
