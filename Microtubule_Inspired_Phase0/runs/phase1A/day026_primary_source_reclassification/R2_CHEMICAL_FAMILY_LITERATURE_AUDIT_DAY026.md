# R2 Chemical Family Literature Audit — Day026

## Objective

Determine which chemical families of the R2 system are already supported by the audited literature (Rajan 2018, Bamane 2023 and Ghorai 2025) and identify the minimum set of chemical environments that still require dedicated QM parametrization.

---

## Chemical Family 1 — Pristine h-BN lattice

### Included node types

- PARENT_HBN

### Literature evidence

#### Rajan (2018)

Direct evidence.

Provides the bonded and nonbonded atomistic force field for pristine h-BN, validated against lattice structure, elastic constants and phonon dispersion.

#### Bamane (2023)

Direct evidence.

Extends the same fixed-topology philosophy to BN nanotubes and confirms transferability to curved BN structures and molecular interfaces.

#### Ghorai (2025)

Direct evidence.

Explicitly reuses the Rajan h-BN parameters as the foundation for the functionalized system.

### Coverage assessment

DIRECTLY SUPPORTED

### Decision

REUSE WITHOUT MODIFICATION

### Scientific justification

The pristine h-BN framework is consistently supported by all three primary references and does not require additional parametrization.

---

## Chemical Family 2 — Curved h-BN Annulus Interior

### Included node types

- ANNULUS_INTERIOR

### Scientific description

These environments correspond to fully coordinated boron and nitrogen atoms
located inside the curved annular scaffold. They preserve the pristine
three-fold BN coordination but experience geometric curvature introduced by
the nanotubular architecture.

No new covalent chemistry is introduced in these environments.

### Literature evidence

#### Rajan (2018)

Provides the atomistic bonded and nonbonded description of pristine
hexagonal boron nitride.

Supports the local chemical identity of these atoms.

#### Bamane (2023)

Provides direct evidence that the same fixed-topology philosophy remains
valid for curved BN nanotubes and reproduces structural and interfacial
properties.

This is the strongest literature support for these environments.

#### Ghorai (2025)

Uses the same h-BN force-field foundation while extending only the
functionalized edge chemistry.

No modification of the annulus interior is introduced.

### Coverage assessment

DIRECT LITERATURE SUPPORT

### Decision

REUSE WITHOUT MODIFICATION

### Scientific justification

The annulus interior does not introduce novel local chemistry.
Its only distinguishing feature relative to pristine h-BN is the
curvature of the scaffold, which is already represented within the
validated BNNT framework described by Bamane.

No dedicated QM parametrization is currently justified.

---


## Chemical Family 3 — Edge BN Environments

### Included node types

- ANNULUS_INNER_BOUNDARY
- ANNULUS_OUTER_BOUNDARY
- HEXAGONAL_EDGE_COMPLETION_SEED

### Scientific description

These environments correspond to boron and nitrogen atoms located at the
structural boundary of the annulus.

Unlike the pristine lattice, these atoms experience coordination
interruption and represent the transition between the basal BN network
and the terminal passivation chemistry.

No chromophore or bridge chemistry is introduced at this stage.

### Literature evidence

#### Rajan (2018)

Provides the reference force field for pristine h-BN but does not
explicitly validate edge chemistry.

#### Bamane (2023)

Validates curved BNNTs and molecular interfaces but does not introduce
dedicated edge parametrization.

#### Ghorai (2025)

Provides the strongest evidence for hydrogen-functionalized edge
environments and explicitly extends the Rajan force field to treat
functionalized BN edge chemistry.

### Coverage assessment

PARTIALLY DIRECTLY SUPPORTED

### Decision

REUSE WITH LOCAL VERIFICATION

### Scientific justification

The literature supports the use of a fixed-topology description for
functionalized BN edges.

However, the exact local coordination present in the R2 annulus should
be verified against the chemical environments parametrized by Ghorai
before complete parameter transfer is authorized.

No evidence currently suggests that a completely new force field is
required for these environments.

---


## Chemical Family 4 — Hydrogen-Terminated BN Edge Environments

### Included node types

- ANNULUS_INNER_PASSIVANT_H
- ANNULUS_OUTER_PASSIVANT_H
- BRIDGE_PASSIVANT_H
- SEED_PASSIVANT_H

### Scientific description

These environments correspond to hydrogen atoms covalently bound to
terminal boron or nitrogen atoms at the annulus boundary.

Their primary role is to saturate dangling bonds generated during the
construction of the finite BN scaffold.

No chromophore attachment is present in these environments.

### Literature evidence

#### Rajan (2018)

Does not include hydrogen-passivated edge chemistry.

The force field is restricted to pristine h-BN.

#### Bamane (2023)

Uses a fixed-topology BNNT description but does not develop dedicated
bonded parameters for hydrogen-terminated BN edge atoms.

#### Ghorai (2025)

Provides direct parametrization of hydrogen-functionalized h-BN edge
chemistry.

The supplementary force-field explicitly introduces bonded and
nonbonded parameters for hydrogen-functionalized BN environments while
retaining the pristine h-BN framework.

### Coverage assessment

DIRECT LITERATURE SUPPORT

### Decision

REUSE WITH LOCAL VERIFICATION

### Scientific justification

Hydrogen-passivated BN environments are explicitly represented within
the Ghorai force field.

Before parameter transfer, the local coordination of the R2 terminal
atoms should be verified to ensure chemical equivalence with the
functionalized environments described in the publication.

No evidence currently justifies a completely new parametrization for
simple BN–H termination.

---


## Chemical Family 5 — Alternating Four-Atom B–N–B–N Bridge

### Included node types

- ALTERNATING_BN_FOUR_ATOM_BRIDGE

### Scientific description

This family corresponds to the alternating B–N–B–N bridge that connects
the reconstructed annular scaffold.

Although each atom preserves the expected elemental valence, the bridge
creates a collective local environment involving new bonded,
angular, torsional and bridge-attachment interactions that are not
explicitly validated in the audited literature.

### Literature evidence

#### Rajan (2018)

No direct representation of the four-atom bridge.

Provides transferable pristine h-BN bonded interactions only.

#### Bamane (2023)

Validates curved BN nanotubes but does not introduce reconstructed
alternating bridge motifs.

#### Ghorai (2025)

Introduces functionalized edge chemistry while preserving the underlying
BN lattice.

No alternating four-atom bridge or equivalent coupled environment is
parametrized.

### Coverage assessment

NO DIRECT LITERATURE SUPPORT

### Decision

QM EXTENSION REQUIRED

### Scientific justification

The individual B–N, B–H and N–H bonds are represented in the literature.

However, the R2 bridge introduces a coupled chemical environment with
new bond, angle, torsional and attachment combinations that are not
validated by any audited force field.

Consequently, transferability cannot be assumed from elemental
similarity alone.

Dedicated QM reference fragments are therefore required to validate
the bridge region before force-field parameter assignment.

---

