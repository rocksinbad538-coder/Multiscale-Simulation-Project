# QM_F06 ORCA Input Technical Audit — Day026

## Decision: **QM_F06_ORCA_INPUTS_PASS_TECHNICAL_AUDIT**

Stage-1 and Stage-2 input files were checked against the repaired atom manifests and the atom-role/constraint map.

## QM_F06_LOWER_CAPPED_REPAIRED

- Stage-1 expected fixed atoms: **11**
- Stage-2 expected fixed atoms: **4**
- Bridge/attachment atoms incorrectly fixed: **0 expected**

## QM_F06_UPPER_CAPPED_REPAIRED

- Stage-1 expected fixed atoms: **11**
- Stage-2 expected fixed atoms: **4**
- Bridge/attachment atoms incorrectly fixed: **0 expected**

## Workflow dependency

The current Stage-2 inputs contain the repaired initial coordinates. They are templates only. Before execution, their coordinate blocks must be replaced with the optimized Stage-1 geometries.

Likewise, Stage 3 must be generated from the optimized Stage-2 geometry rather than from the current initial coordinates.

## Authorization state

- Static input audit: **PASSED**
- Sequential workflow preparation: **PENDING**
- QM execution: **NOT AUTHORIZED**

