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

## Status update — 2026-07-31

The UPPER fragment is no longer pending.

`QM_F06_UPPER_V7A_R1` has completed:

- ORCA geometry optimization;
- normal ORCA termination;
- independent post-QM structural acceptance;
- final-coordinate consistency validation;
- coordinate adoption with SHA256 provenance;
- RESP input-design validation;
- Day036 RESP preparation.

The accepted UPPER reference contains 52 atoms with composition
B17N14H21. It includes 37 real atoms and 15 artificial QM boundary
caps.

Artificial caps remain part of the finite QM electrostatic model but
are excluded from direct transfer into the full scaffold and from
transferable atom-type averages.

No non-singleton charge-equivalence restraints have yet been enforced.
Candidate equivalence groups require joint topology and local-geometry
validation across LOWER and UPPER fragments.

The independent ORCA CHELPG input preflight has passed. The ESP
single-point execution is formally authorized.

RESP input generation, RESP fitting, charge adoption, force-field
adoption and molecular dynamics remain unauthorized until their
respective downstream gates pass.

## ESP result update — 2026-07-31

The authorized QM_F06 UPPER V7-A R1 ORCA CHELPG single-point
calculation completed successfully.

- SCF converged after 17 cycles.
- Final energy: -1201.445434386447 Eh.
- ORCA terminated normally.
- Shell exit status: 0.
- Standard error output: empty.
- The CHELPG charge block was generated and independently audited.
- Input identity and geometry provenance remained intact.

CHELPG is itself an ESP-fitted charge model. Its fitted atomic charges
are retained as diagnostic electronic-property results.

Before AmberTools RESP input generation can be authorized, a separate
stage must demonstrate how the raw ORCA electrostatic-potential data
will be extracted or converted into the ESP format consumed by RESP.

RESP input generation, RESP execution, charge adoption, force-field
adoption and molecular dynamics remain unauthorized.
