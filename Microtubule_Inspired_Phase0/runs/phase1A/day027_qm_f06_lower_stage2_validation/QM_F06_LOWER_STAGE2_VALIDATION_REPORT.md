# QM_F06 LOWER Stage-2 Validation — Day027

## Decision: **QM_F06_LOWER_STAGE2_OR_STAGE3_STATIC_GATE_FAILED**

## Stage-2 calculation

- ORCA termination: **NORMAL**
- Geometry optimization: **CONVERGED**
- Method: **PBE0-D4/def2-TZVP**
- Charge/multiplicity: **0 / 1**

## Stage 1 → Stage 2 displacement

- Direct Cartesian RMSD: **0.448949 Å**
- Maximum atomic displacement: **1.208706 Å**
- Maximum-displacement atom: `HCAP:LOWER:03`

## Bonded structure

- Reconstructed bonds: **21**
- Bond-range failures: **0**
- Artificial-cap bond failures: **0**

## Stage-1 relaxation targets

- Target contacts evaluated: **2**
- Both contacts improved: **YES**
- `BR4:LOWER:00:3 — P:48`: 2.6558 Å → 2.7538 Å; improved = **True**
- `BR4:LOWER:00:3 — H4:LOWER:0017:0`: 2.1278 Å → 2.2602 Å; improved = **True**

## Long-range contacts

- Contacts below 0.90 of vdW sum: **11**
- Hard contacts below 0.70: **1**
- Hard contacts involving artificial caps: **1**

## Stage-3 static input

- Static syntax/content gate: **PASS**
- Exact Stage-2 coordinate propagation: **PASS**
- Maximum coordinate-transfer error: **4.656e-11 Å**

## Authorization state

- Stage-2 structural validation: **FAIL**
- Stage-3 static input readiness: **PASS**
- Stage-3 execution: **NOT YET AUTHORIZED**
- RESP/ESP charge protocol: **NOT YET DEFINED**
- Force-field parameter adoption: **NOT AUTHORIZED**

## Required scientific decision

Define the Stage-3 electronic-property protocol. The current input is a valid energy single point, but an energy-only calculation is insufficient for RESP or electrostatic-potential-derived charges.

