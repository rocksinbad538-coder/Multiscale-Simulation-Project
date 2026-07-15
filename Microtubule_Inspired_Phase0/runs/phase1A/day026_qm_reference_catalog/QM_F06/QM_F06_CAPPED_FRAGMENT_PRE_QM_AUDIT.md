# QM_F06 Capped-Fragment Pre-QM Audit — Day026

## Scope

The chemically capped LOWER and UPPER bridge fragments were audited for atom counts, graph valence, bonded distances, nonbonded clashes and electronic-state parity.

## Decision: **QM_F06_CAPPED_FRAGMENTS_PASS_PRE_QM_GATE**

## QM_F06_LOWER_CAPPED

- Formula: **B5N5H12**
- Atoms: **22**
- Reconstructed bonds: **21**
- Degree/valence failures: **0**
- Bond-distance failures: **0**
- Minimum nonbonded distance: **1.218614 Å** (`H4:LOWER:0017:0 — HCAP:LOWER:07`)
- Steric clashes: **0**
- Neutral valence-electron count: **52**
- Provisional electronic state: **charge 0, multiplicity 1**
- Pre-QM gate: **PASS**

## QM_F06_UPPER_CAPPED

- Formula: **B5N5H12**
- Atoms: **22**
- Reconstructed bonds: **21**
- Degree/valence failures: **0**
- Bond-distance failures: **0**
- Minimum nonbonded distance: **1.389046 Å** (`H4:UPPER:0047:0 — HCAP:UPPER:07`)
- Steric clashes: **0**
- Neutral valence-electron count: **52**
- Provisional electronic state: **charge 0, multiplicity 1**
- Pre-QM gate: **PASS**

## Charge and multiplicity interpretation

Both fragments have formula B5N5H12 and 52 neutral valence electrons. A neutral closed-shell singlet (charge 0, multiplicity 1) is therefore the provisional initial electronic state.

This is an electron-count and stoichiometric assignment, not yet an electronic-structure validation. The eventual QM workflow must confirm SCF stability and absence of a lower-energy open-shell solution.

## Authorization state

- Geometry construction: **COMPLETED**
- Pre-QM audit: **COMPLETED**
- QM input preparation: **PENDING THIS GATE**
- QM calculation executed: **NO**
