# QM_F06 Pair-Specific Steric Audit — Day026

## Rationale

The previous 0.70 Å universal threshold only excluded near-coincident atoms. The present audit uses interatomic distances normalized by the sum of elemental van der Waals radii.

## Decision: **QM_F06_CAPPED_FRAGMENTS_REQUIRE_CAP_GEOMETRY_REPAIR**

## QM_F06_LOWER_CAPPED

- Nonbonded contacts below vdW sum: **61**
- Classification counts: `{'ACCEPTABLE': 15, 'CLOSE_CONTACT': 19, 'SEVERE_CLASH': 18, 'STRONG_COMPRESSION': 9}`
- Most compressed contact: `BR4:LOWER:00:1 — P:48`
- Distance: **1.886000 Å**
- Distance/vdW-sum ratio: **0.4911**
- Pair-specific steric gate: **FAIL**

## QM_F06_UPPER_CAPPED

- Nonbonded contacts below vdW sum: **62**
- Classification counts: `{'ACCEPTABLE': 16, 'CLOSE_CONTACT': 16, 'SEVERE_CLASH': 17, 'STRONG_COMPRESSION': 13}`
- Most compressed contact: `A:UPPER:14:4 — BR4:UPPER:00:4`
- Distance: **1.701843 Å**
- Distance/vdW-sum ratio: **0.4432**
- Pair-specific steric gate: **FAIL**

## Interpretation

SEVERE_CLASH and STRONG_COMPRESSION contacts must be corrected before electronic-structure input preparation. CLOSE_CONTACT geometries may be retained as initial structures but must be monitored during optimization.

## Authorization state

- Pair-specific steric audit: **COMPLETED**
- QM input preparation: **NOT AUTHORIZED**
- QM calculation executed: **NO**

