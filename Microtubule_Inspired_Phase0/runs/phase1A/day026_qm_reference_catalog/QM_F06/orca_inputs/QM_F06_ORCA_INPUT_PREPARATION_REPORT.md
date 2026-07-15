# QM_F06 ORCA Input Preparation — Day026

## Decision

**QM_F06_ORCA_INPUTS_PREPARED_EXECUTION_NOT_AUTHORIZED**

## Electronic-structure model

- Functional: **PBE0**
- Dispersion: **D4**
- Basis: **def2-TZVP**
- Approximation: **RIJCOSX / def2-J**
- SCF: **TightSCF, maximum 500 iterations**
- Charge: **0**
- Multiplicity: **1**

## Optimization hierarchy

1. Stage 1 keeps artificial caps and peripheral scaffold atoms fixed.
2. Stage 2 releases all hydrogen atoms while retaining peripheral heavy-atom constraints.
3. Stage 3 is a single-point template and must receive the optimized Stage-2 coordinates.

## QM_F06_LOWER_CAPPED_REPAIRED

- Atoms: **22**
- Stage-1 fixed atoms: **11**
- Stage-1 mobile atoms: **11**
- Stage-2 fixed atoms: **4**
- Stage-2 mobile atoms: **18**
- Charge/multiplicity: **0 / 1**

## QM_F06_UPPER_CAPPED_REPAIRED

- Atoms: **22**
- Stage-1 fixed atoms: **11**
- Stage-1 mobile atoms: **11**
- Stage-2 fixed atoms: **4**
- Stage-2 mobile atoms: **18**
- Charge/multiplicity: **0 / 1**

## Executable detection

- ORCA detected in current shell: **YES**
- Detected path: `/Users/alejandro/projects/orca_6_1_1_macosx_intel_openmpi411/orca `

## Authorization state

- Input preparation: **COMPLETED**
- QM execution: **NOT AUTHORIZED**
- Charge fitting: **NOT YET DEFINED**
- Force-field fitting: **NOT YET EXECUTED**

