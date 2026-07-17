# QM_F06 UPPER Workflow Decision — Day028

## Scientific objective

Construct and independently validate the QM_F06 UPPER bridge reference
fragment before deriving transferable electrostatic or bonded parameters.

## LOWER reference status

`QM_F06_LOWER_BOUNDARY_V2B` is accepted at the structural and electronic
levels.

The LOWER result must be used as a methodological reference, not copied
as evidence that the UPPER environment is equivalent.

## UPPER workflow requirements

The UPPER fragment must independently demonstrate:

- chemically complete boundary construction;
- preservation of the B–N–B–N bridge;
- valid covalent connectivity and valence;
- absence of cap-induced artifacts;
- converged QM geometry;
- electronic validation of any compressed real contacts;
- quantitative comparison with LOWER.

## Immediate task

Audit the existing UPPER 22-atom repaired fragment and identify the
minimum chemically complete boundary expansion required for an
independent UPPER reference calculation.

## Authorization state

- UPPER topology/boundary audit: AUTHORIZED
- UPPER boundary construction: PENDING AUDIT
- UPPER ORCA execution: NOT AUTHORIZED
- ESP/RESP execution: NOT AUTHORIZED
- Charge adoption: NOT AUTHORIZED
- Force-field parameter adoption: NOT AUTHORIZED
