# Phase 1A accepted structural baselines

## Status

This directory freezes the accepted Phase 1A structural baselines generated during Day 009.

These structures are accepted as geometry-audited starting points for subsequent force-field assignment, hydration, short MD preparation, representative fragment extraction, and later DFT/TDDFT parameterization.

## Accepted structures

### 1. Atomistic h-BN nanotube, non-passivated

Directory:

- `hbn_noH/`

Files:

- `hbn_nt_initial_noH.xyz`
- `hbn_nt_initial_noH.pdb`
- `audit_hbn_summary_noH.csv`
- `audit_hbn_distances_noH.txt`
- `audit_hbn_overlaps_noH.txt`
- `audit_hbn_connectivity_noH.txt`
- `audit_hbn_bonds_noH.txt`

Accepted audit values:

- B atoms: 840
- N atoms: 840
- H atoms: 0
- total atoms: 1680
- B-N bonds: 2460
- B-N mean distance: approximately 1.44967 Å
- overlap count: 0

Interpretation:

Accepted as an atomistic h-BN structural scaffold.

Not labeled QM-ready because finite edges are not passivated.

### 2. Standalone pyrene

Directory:

- `pyrene/`

Files:

- `pyrene.xyz`
- `pyrene.pdb`
- `audit_pyrene_summary.csv`
- `audit_pyrene_bonds.csv`
- `audit_pyrene_overlaps.csv`

Accepted audit values:

- formula: C16H10
- total atoms: 26
- C atoms: 16
- H atoms: 10
- nonbonded overlap count: 0
- carbon-plane RMS deviation: approximately 0.000005 Å

Interpretation:

Accepted as the first standalone methodological chromophore geometry.

Not labeled TDDFT-final; later electronic-structure optimization may still be required.

### 3. Dry h-BN + four-pyrene hybrid

Directory:

- `hybrid_dry/`

Files:

- `hbn_pyrene_4_dry.xyz`
- `hbn_pyrene_4_dry.pdb`
- `chromophore_mapping.csv`
- `audit_chromophore_distances.csv`
- `audit_overlaps.csv`
- `audit_summary.csv`

Accepted audit values:

- h-BN atoms: 1680
- pyrene molecules: 4
- pyrene atoms total: 104
- total atoms: 1784
- intercomponent overlap count: 0
- pyrene-hBN minimum distances: approximately 3.30-3.36 Å
- nearest pyrene-pyrene minimum atom-atom distance: approximately 19.30 Å

Interpretation:

Accepted as the initial dry atomistic h-BN + four-pyrene hybrid geometry.

Not yet MD-ready final because force-field types, charges, and simulation parameters are still pending.

## Script versions

Directory:

- `scripts/`

Frozen scripts:

- `build_hbn_nanotube.py`
- `build_pyrene.py`
- `place_pyrene_4_on_hbn.py`

## Logs

Directory:

- `logs/`

Available build logs are copied when present.

## Use constraints

Accepted labels:

- atomistic h-BN initial scaffold
- standalone pyrene initial geometry
- dry h-BN + four-pyrene hybrid initial geometry

Rejected labels at this stage:

- QM-ready full system
- TDDFT-ready full system
- MD-production-ready system
- final passivated h-BN nanotube
- validated excitonic Hamiltonian
