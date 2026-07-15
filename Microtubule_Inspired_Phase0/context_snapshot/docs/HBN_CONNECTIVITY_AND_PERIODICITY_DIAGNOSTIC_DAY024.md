# HBN Connectivity and Periodicity Diagnostic

## Scope

This diagnostic investigates why the original Day024 auditor found no
degree-2 terminal atoms.

No minimization, molecular dynamics, topology generation, or quantum
calculation was executed.

## Geometry

- Box:
  **7.394530, 7.394530, 10.017500 nm**
- HBN coordinate span:
  **2.398000, 2.398000,
  6.018000 nm**
- PCA tube axis:
  **0.00000000, 0.00000000,
  1.00000000**
- Axial tube span:
  **6.018000 nm**
- Box extent projected on tube axis:
  **10.017500 nm**

## Original geometric bond rule

Original search band:

- minimum: **0.115000 nm**
- maximum: **0.175000 nm**

Raw, nonperiodic bond count in this band:

- **2460**

Minimum-image XYZ bond count in this band:

- **2460**

Degree graph using all B–N distances up to 0.175 nm:

- Raw:
  bonds=2460,
  d0/d1/d2/d3/d4+ =
  0/
  60/
  0/
  1620/
  0
- Minimum-image XYZ:
  bonds=2460,
  d0/d1/d2/d3/d4+ =
  0/
  60/
  0/
  1620/
  0

Potential PBC-created pairs below 0.175 nm:

- **0**

B–N pairs below 0.100 nm:

- Raw: **0**
- Minimum-image XYZ: **0**

## Exact expected graph matches

- NONE

## Ten best cutoff/mode candidates

- RAW_NO_PBC: cutoff=0.1475 nm; score=480; bonds=2460; d0/d1/d2/d3/d4+=0/60/0/1620/0
- MINIMUM_IMAGE_XYZ: cutoff=0.1475 nm; score=480; bonds=2460; d0/d1/d2/d3/d4+=0/60/0/1620/0
- RAW_NO_PBC: cutoff=0.1500 nm; score=480; bonds=2460; d0/d1/d2/d3/d4+=0/60/0/1620/0
- MINIMUM_IMAGE_XYZ: cutoff=0.1500 nm; score=480; bonds=2460; d0/d1/d2/d3/d4+=0/60/0/1620/0
- RAW_NO_PBC: cutoff=0.1525 nm; score=480; bonds=2460; d0/d1/d2/d3/d4+=0/60/0/1620/0
- MINIMUM_IMAGE_XYZ: cutoff=0.1525 nm; score=480; bonds=2460; d0/d1/d2/d3/d4+=0/60/0/1620/0
- RAW_NO_PBC: cutoff=0.1550 nm; score=480; bonds=2460; d0/d1/d2/d3/d4+=0/60/0/1620/0
- MINIMUM_IMAGE_XYZ: cutoff=0.1550 nm; score=480; bonds=2460; d0/d1/d2/d3/d4+=0/60/0/1620/0
- RAW_NO_PBC: cutoff=0.1575 nm; score=480; bonds=2460; d0/d1/d2/d3/d4+=0/60/0/1620/0
- MINIMUM_IMAGE_XYZ: cutoff=0.1575 nm; score=480; bonds=2460; d0/d1/d2/d3/d4+=0/60/0/1620/0

## Status

- Auditor repair authorized: **NO**
- MD executed: **NO**
- QM executed: **NO**
- Required next step:
  `CLASSIFY_CONNECTIVITY_FAILURE_AND_REPAIR_AUDITOR`
