# QM_F06 ESP/RESP Protocol Definition — Day028

## Status

Protocol definition authorized.

ESP/RESP calculation execution remains unauthorized.

## Scientific objective

Derive transferable electrostatic parameters for the real R2 chemical
environments represented by QM_F06 without transferring artificial
boundary-cap artifacts into the full scaffold force field.

## Accepted reference geometry

`QM_F06_LOWER_BOUNDARY_V2B`

Method used for structural and electronic validation:

`PBE0-D4/def2-TZVP`, neutral singlet.

## Atom classes requiring explicit treatment

1. Alternating B–N–B–N bridge atoms.
2. Bridge passivant B–H and N–H atoms.
3. Annulus bridge-attachment atoms.
4. Parent/seed bridge-attachment atoms.
5. Real annulus edge hydrogens.
6. Artificial QM boundary caps.

## Required protocol decisions

### 1. Charge-fitting target

Determine whether the production charge model will use:

- RESP;
- restrained ESP with a non-AMBER implementation;
- MBIS-derived charges;
- a hybrid transferable charge model.

No model is selected solely from the diagnostic Hirshfeld, MBIS or
CHELPG outputs.

### 2. Treatment of artificial caps

Artificial cap charges must not be transferred directly into the full
R2 scaffold.

The protocol must specify whether cap charges are:

- excluded from the transferable atom-type averages;
- constrained to predefined values;
- absorbed into adjacent real-atom charges;
- removed through a charge-renormalization scheme.

### 3. Charge conservation

The transferable real-atom charge set must satisfy:

- the intended formal charge of the corresponding full-system region;
- elemental and environment-specific equivalence constraints;
- reproducibility across LOWER and UPPER fragments.

### 4. Symmetry and equivalence restraints

Candidate equivalent environments include:

- chemically equivalent bridge B atoms;
- chemically equivalent bridge N atoms;
- equivalent B–H passivants;
- equivalent N–H passivants;
- symmetry-related LOWER/UPPER attachment environments.

Equivalence must be demonstrated from topology and local geometry rather
than assumed from elemental identity alone.

### 5. Conformational sampling

A single optimized geometry may be insufficient for transferable charge
derivation.

The protocol must determine whether to include:

- LOWER V2-B geometry;
- independently optimized UPPER geometry;
- additional constrained conformers;
- representative geometries from later thermal sampling.

### 6. Validation criteria

Any adopted charge set must be checked for:

- exact total-charge reproduction;
- stability with respect to conformer selection;
- LOWER/UPPER consistency;
- electrostatic-potential reproduction;
- absence of cap-dominated artifacts;
- compatibility with the eventual force-field implementation.

## Current authorization state

- Protocol design: **IN PROGRESS**
- Additional diagnostic analysis: **AUTHORIZED**
- LOWER RESP/ESP execution: **NOT AUTHORIZED**
- UPPER fragment construction/optimization: **PENDING**
- Charge adoption: **NOT AUTHORIZED**
- Force-field parameter adoption: **NOT AUTHORIZED**
