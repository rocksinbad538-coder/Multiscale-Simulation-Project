# Phase 1A.5 — Force-field and charge-assignment plan

## Objective

Define a defensible force-field and charge-assignment strategy before adding water or preparing molecular dynamics input files.

This checkpoint is required because the accepted Phase 1A structures are geometrically valid but not yet MD-ready.

Accepted structures so far:

- atomistic non-passivated h-BN nanotube;
- standalone pyrene;
- dry h-BN + four-pyrene hybrid.

None of these structures is yet labeled MD-production-ready, QM-ready, or TDDFT-ready.

## Scientific constraint

No bonded, Lennard-Jones, or charge parameters should be invented.

All parameters must come from one of the following:

1. established force-field literature;
2. force-field tool assignment with audit;
3. QM-derived charges;
4. clearly labeled provisional exploratory parameters.

## Component-level strategy

### 1. Pyrene

Pyrene is an organic aromatic molecule.

Recommended initial force-field route:

- GAFF/GAFF2 or OPLS-AA as the first classical MD route;
- AM1-BCC as a practical first-pass charge model;
- RESP charges from DFT as the more rigorous later option.

Minimum required checks:

- atom typing succeeds;
- total molecular charge = 0;
- aromatic carbon and hydrogen types are assigned consistently;
- bonded terms exist for all bonds, angles, dihedrals;
- no missing parameters.

Accepted status after parameter assignment:

- MD-ready preliminary pyrene if all parameters are assigned and charge neutrality passes.

Rejected status:

- TDDFT-ready pyrene. TDDFT requires separate electronic-structure preparation and likely geometry reoptimization.

### 2. h-BN nanotube

The h-BN scaffold is inorganic and cannot be blindly parameterized with generic organic force fields.

Possible routes:

#### Route A — Fixed rigid or restrained scaffold

Use the h-BN nanotube as a rigid or strongly restrained structural scaffold during early MD.

Advantages:

- minimizes risk from uncertain B-N bonded parameters;
- preserves the validated geometry;
- allows pyrene and water environment testing;
- appropriate for early hydration and placement screening.

Limitations:

- does not test intrinsic h-BN flexibility;
- not sufficient for final material mechanical response.

Recommended for the first short MD controls.

#### Route B — Literature h-BN force-field parameters

Use B-N bonded and nonbonded parameters from literature for h-BN sheets/nanotubes.

Advantages:

- physically better if validated;
- allows scaffold vibrations and thermal stability checks.

Limitations:

- requires source verification;
- must audit B-N equilibrium distances, angle terms, partial charges, and LJ parameters;
- cannot proceed until parameters are cited and tested.

Recommended after Route A if a reliable parameter set is identified.

#### Route C — QM-derived or fitted parameters

Derive parameters from quantum calculations or fit to reference data.

Advantages:

- highest rigor.

Limitations:

- too expensive and premature for Phase 1A geometry screening.

Not recommended as the immediate next step.

### 3. Water

Recommended water model for initial MD:

- TIP3P as first-pass compatibility model if using GAFF/OPLS workflows;
- SPC/E or TIP4P variants can be considered later for dipolar/dielectric sensitivity tests.

Minimum required checks:

- water model declared explicitly;
- no mixing of incompatible water parameters without justification;
- water geometry and charge model documented.

### 4. h-BN / pyrene / water nonbonded interactions

Initial cross-interactions may use Lorentz-Berthelot mixing only if the underlying force field supports it.

Minimum required checks:

- all atom types have epsilon/sigma or equivalent LJ parameters;
- no zeroed nonbonded interactions unless deliberately documented;
- no artificial bonding between pyrene and h-BN;
- pyrene remains physisorbed, not covalently attached.

## Recommended immediate MD strategy

Use a staged route:

### Stage 1 — Dry hybrid sanity test

System:

- h-BN nanotube fixed or strongly restrained;
- four pyrenes mobile or weakly restrained;
- no water.

Purpose:

- check pyrene placement;
- check nonbonded interactions;
- check absence of violent forces;
- perform minimization only.

### Stage 2 — Hydrated hybrid minimization

System:

- h-BN nanotube fixed or restrained;
- four pyrenes;
- explicit water.

Purpose:

- remove bad contacts;
- audit water placement;
- check pyrene displacement;
- check water occupancy inside/outside nanotube.

### Stage 3 — Very short NVT

System:

- same as Stage 2.

Purpose:

- 10–50 ps sanity dynamics;
- no production claims;
- monitor temperature, energy, pyrene positions, h-BN geometry, and water confinement.

## Current decision

Immediate next technical path:

1. keep h-BN scaffold fixed or strongly restrained for first MD-preparation tests;
2. parameterize pyrene using an established organic force-field route;
3. add water only after dry force-field sanity checks;
4. do not run production MD until minimization and short NVT audits pass.

## Non-negotiable constraints

Do not claim:

- QM-ready full hybrid;
- TDDFT-ready full hybrid;
- exciton-ready Hamiltonian;
- validated optical response;
- spin-active system.

until each corresponding module is explicitly parameterized and audited.

