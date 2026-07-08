# R2 Static Instability-Signature Audit

## Scope

This audit evaluates the completed R2 zero-step single-point run. It
does not rerun `grompp`, `mdrun`, minimization, or molecular dynamics.

The earlier terminal wrapper used the unrestricted expression `nan`,
which can match ordinary text containing those three characters. The
strict audit uses the bounded numerical token `\bnan\b` and explicit
GROMACS failure signatures.

## Core Gate 2B state

- Decision:
  **R2_TOPOLOGY_AND_STATIC_CAP_WATER_MODEL_VALIDATED**
- Expected decision:
  **R2_TOPOLOGY_AND_STATIC_CAP_WATER_MODEL_VALIDATED**
- Core gates:
  **29**
- All core gates passed:
  **YES**

## Pattern results

- Legacy broad-pattern matches:
  **1**
- Strict instability matches:
  **0**
- Broad-pattern false positives:
  **1**

### Broad-pattern-only matches

- `r2_static_single_point.log:230` — `nan` in: `   lincs-warnangle                = 30`

### Strict instability matches

- NONE

## Decision

- Audit decision:
  **R2_STATIC_INSTABILITY_SCAN_CONFIRMED_CLEAN**
- Static or MD rerun required:
  **NO**
- Water-only energy minimization authorized:
  **YES**
- Short frozen-solute NVT authorized:
  **NO**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `RUN_R2_WATER_ONLY_ENERGY_MINIMIZATION`
