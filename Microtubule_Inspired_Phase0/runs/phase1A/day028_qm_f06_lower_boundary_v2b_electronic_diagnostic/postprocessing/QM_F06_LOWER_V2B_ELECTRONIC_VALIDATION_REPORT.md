# QM_F06 LOWER V2-B Electronic Validation — Day028

## Decision: **QM_F06_LOWER_V2B_RESIDUAL_CONTACT_HAS_NO_SIGNIFICANT_COVALENT_CHARACTER_ELECTRONIC_GATE_PASS**

## Execution

- Return code: **0**
- SCF convergence: **YES**
- Normal ORCA termination: **YES**
- Final single-point energy: **-541.318700348194 Eh**

## Target contact

- B atom: `BR4:LOWER:00:3` — ORCA index `5`
- H atom: `H4:LOWER:0017:0` — ORCA index `20`
- Mayer print threshold: **0.01**
- Target Mayer bond order printed: **NO**
- Target Mayer bond-order result: **< 0.01**

## Covalent-reference bonds

- Intended B5–H9 Mayer bond order: **0.9964**
- Intended N2–H20 Mayer bond order: **0.8695**

## Charge analyses

- HIRSHFELD: sum = **+2.000000e-06 e**; B5 = **+0.106987 e**; H20 = **+0.140788 e**
- MBIS: sum = **+5.000000e-06 e**; B5 = **+0.979552 e**; H20 = **+0.408799 e**
- CHELPG: sum = **-1.000000e-06 e**; B5 = **+0.501677 e**; H20 = **+0.303981 e**

## Interpretation

The B5···H20 target pair is absent from the Mayer table printed at a 0.01 threshold, while the intended local B–H and N–H bonds have Mayer bond orders close to unity. The compressed 1–5 contact therefore has no significant covalent-bond character at the PBE0-D4/def2-TZVP level.

Hirshfeld, MBIS and CHELPG values remain diagnostic only. They are not adopted as force-field charges by this validation.

## Authorization state

- LOWER V2-B geometry formally accepted: **YES**
- ESP/RESP protocol definition: **AUTHORIZED**
- ESP/RESP execution: **NOT AUTHORIZED**
- Charge adoption: **NOT AUTHORIZED**
- Force-field parameter adoption: **NOT AUTHORIZED**

