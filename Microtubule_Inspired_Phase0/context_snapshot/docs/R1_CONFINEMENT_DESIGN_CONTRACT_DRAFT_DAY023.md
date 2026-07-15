# R1 Confinement Design Contract — Draft

## Purpose

R1 is a fully capped positive control for persistent lumen-water
confinement. It is not assumed to be the final device architecture.

## Invariants inherited from R0

R1 must preserve, unless an explicit audit demonstrates otherwise:

- the validated h-BN scaffold;
- all four pyrene chromophores;
- atom ordering for the inherited R0 atoms;
- the accepted force-field parameters;
- TIP4P/2005 water;
- the simulation box;
- the 300 K reference condition;
- the established lumen-axis convention;
- the existing structural-analysis definitions.

## Starting-state requirement

R1 must start from the earliest hydrated frame of the accepted
frozen-solute R0 trajectory.

The depleted Stage02/mobile-branch state must not be used as the
primary R1 positive-control starting state.

## Cap requirements

The initial R1 caps must:

- close both axial exits;
- avoid atomic overlaps with HBN, PYR, and water;
- preserve the accessible lumen interior;
- introduce no uncontrolled net charge;
- have explicitly documented composition and bonding;
- be generated reproducibly by a script under `scripts/phase1A/`.

## Gate sequence

1. Identify and extract the authoritative hydrated R0 starting state.
2. Quantify lumen axis, end planes, radius, and accessible volume.
3. Generate an R1 cap prototype.
4. Audit composition, bonding, overlaps, charge, and geometry.
5. Minimize and perform static preprocessing.
6. Run only a short frozen-solute confinement screening.
7. Authorize longer or mobile simulations only after screening passes.

## Current prohibitions

- No long mobile MD.
- No multitemperature production.
- No new QM calculations.
- No replacement of the accepted R0 baseline.
- No interpretation of R1 as a final experimental architecture before
  R2 and R3 are evaluated.

## Provisional scientific objective

R1 should demonstrate that eliminating axial escape prevents the
progressive depletion observed in R0. Its role is to validate the
confinement methodology and establish a positive-control reference.
