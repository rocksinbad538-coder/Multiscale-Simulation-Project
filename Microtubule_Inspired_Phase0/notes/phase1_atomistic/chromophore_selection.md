# Phase 1A.1 — Chromophore selection freeze

## Objective

Freeze the first real chromophore to be used in the atomistic Phase 1 model before generating hybrid h-BN/chromophore/water systems.

The selected chromophore must support the master multiscale workflow:

Structural model → MD → representative fragment extraction → DFT/TDDFT → site energies and transition dipoles → excitonic Hamiltonian → open-system dynamics → optical/electromagnetic response.

## Decision

Pyrene is selected as the first methodological chromophore for Phase 1A.

This does not mean pyrene is the final biological chromophore of the project. It is selected as the first controlled aromatic chromophore to validate the electronic-structure-to-exciton pipeline.

## Rationale

Pyrene is suitable for the first Phase 1A implementation because it is:

- small enough for repeated DFT/TDDFT calculations;
- rigid and planar, reducing conformational ambiguity;
- aromatic and electronically active;
- compatible with Frenkel-exciton parameterization;
- photophysically well studied;
- useful for testing chromophore–chromophore coupling;
- easier to audit geometrically than larger chromophores;
- suitable for controlled placement near an h-BN nanotube scaffold.

## Scientific role in this project

Pyrene will be used to validate:

1. real chromophore placement around the atomistic scaffold;
2. chromophore–chromophore distance and orientation control;
3. extraction of chemically valid representative fragments;
4. TDDFT calculation of excitation energies and transition dipoles;
5. construction of a minimal four-site excitonic Hamiltonian;
6. sensitivity of site energies and couplings to hydration and local environment.

## Alternatives considered

### Indole / tryptophan-like chromophore

Advantages:

- stronger biomimetic connection to microtubule aromatic residues;
- relevant to tryptophan-network hypotheses.

Limitations for first implementation:

- smaller oscillator strengths and potentially weaker excitonic coupling;
- more difficult biological interpretation if introduced before the inorganic scaffold route is stable;
- better suited as a biomimetic control after the pyrene pipeline is validated.

Decision:

Retain as a later biomimetic control.

### Perylene

Advantages:

- stronger extended π-system;
- potentially stronger excitonic coupling.

Limitations:

- larger system;
- stronger aggregation/excimer complications;
- higher TDDFT cost.

Decision:

Retain as a later stronger-coupling control.

### Anthracene

Advantages:

- small PAH;
- electronically active.

Limitations:

- less robust than pyrene as the first methodological chromophore;
- potentially less useful for controlled pyrene-like excimer/exciton benchmarking.

Decision:

Not selected for Phase 1A.

### Porphyrin

Advantages:

- strong optical transitions;
- biologically and optically relevant.

Limitations:

- much larger;
- more complex electronic structure;
- possible metal/spin complications depending on substitution;
- not appropriate before the basic pipeline is validated.

Decision:

Not selected for Phase 1A.

## Accepted Phase 1A chromophore status

Accepted:

- pyrene as first methodological chromophore;
- four pyrene molecules as minimal excitonic test network;
- dry and hydrated controls;
- field-free and field-on controls after geometry validation.

Not accepted yet:

- pyrene as final biological chromophore;
- TDDFT-ready hybrid fragments;
- exciton-ready Hamiltonian;
- spin-active model.

## Required validation before Phase 1B

Before adding water or running MD, the following must be generated and audited:

- `pyrene.xyz`
- `pyrene.pdb`
- `hbn_pyrene_4_dry.xyz`
- `hbn_pyrene_4_dry.pdb`
- `chromophore_mapping.csv`
- `audit_chromophore_distances.csv`
- `audit_overlaps.csv`

Minimum checks:

- pyrene bond lengths chemically plausible;
- four pyrenes uniquely indexed;
- no pyrene–pyrene overlaps;
- no pyrene–h-BN overlaps;
- controlled pyrene–pyrene distances;
- controlled pyrene orientation relative to nanotube axis;
- clear mapping from atoms to chromophore IDs.

## Current conclusion

Proceed with pyrene as the first real chromophore for Phase 1A.

Next technical block:

Generate and audit a standalone pyrene molecule before placing four pyrenes around the h-BN nanotube.
