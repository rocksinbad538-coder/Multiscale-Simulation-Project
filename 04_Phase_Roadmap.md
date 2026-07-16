# 04_PHASE_ROADMAP.md

# Multiscale Simulation – Microtubule Inspired Artificial System
## Phase Roadmap

Version: Day 026

---

# Purpose

This document tracks the scientific progress of every project phase described in the Master Project Document.

Unlike the daily logs, this roadmap should evolve slowly.

Its purpose is to answer one question:

"Where are we in the complete project?"

---

# Overall project roadmap

Phase 0

↓

Phase 1A

↓

Phase 1B

↓

Phase 2

↓

Phase 3

↓

Phase 4

↓

Phase 5

↓

Final validation

↓

Final deliverables

---

# PHASE 0
## Structural system definition

Objective

Design the complete atomistic model of the artificial microtubule-inspired system.

Major tasks

✓ Selection of scaffold

✓ Chromophore placement

✓ Structural validation

✓ Initial molecular model

Status

COMPLETED

Deliverables

✓ Atomistic structure

✓ Initial coordinates

✓ Functionalized model

Dependencies

None

---

# PHASE 1A
## Classical molecular simulation

Objective

Generate a physically realistic molecular dynamics trajectory suitable for electronic-structure calculations.

Major tasks

✓ Force-field implementation

✓ Molecular dynamics

✓ Frame extraction

✓ Chromophore mapping

✓ Local QM clusters

✓ Electrostatic embedding

✓ Embedded TDDFT workflow

✓ Site-energy extraction

✓ Diagonal Hamiltonians

Current work

• Force-field comparison

• Parametrization decision

• Remaining production TDDFT calculations

Remaining before closing Phase 1A

□ Complete all 84 embedded TDDFT calculations

□ Extract complete site-energy trajectory

□ Final statistical validation

□ Finalize diagonal Hamiltonian dataset

Current status

ACTIVE

Estimated completion

Near-term

---

# PHASE 1B
## Electronic couplings

Objective

Construct the complete excitonic Hamiltonian.

Major tasks

□ Electronic couplings

□ Coupling validation

□ Time-dependent Hamiltonians

Expected outputs

H(t)

including

diagonal terms

+

off-diagonal couplings

Dependencies

Requires completion of Phase 1A.

Status

NOT STARTED

---

# PHASE 2
## Open quantum dynamics

Objective

Simulate exciton dynamics.

Candidate methodologies

• Lindblad

• Haken–Strobl

• Redfield (if justified)

Expected outputs

Population dynamics

Coherence dynamics

Relaxation

Transfer efficiencies

Dependencies

Complete Hamiltonian.

Status

NOT STARTED

---

# PHASE 3
## Optical properties

Objective

Predict optical observables.

Expected calculations

Absorption spectra

Emission

Linear response

Electronic transitions

Excitonic spectra

Dependencies

Quantum dynamics.

Status

NOT STARTED

---

# PHASE 4
## Electromagnetic response

Objective

Evaluate electromagnetic properties of the complete artificial system.

Topics

THz response

Microwave response

Field interactions

Collective behavior

Dependencies

Electronic structure

+

Quantum dynamics

Status

NOT STARTED

---

# PHASE 5
## Integration and validation

Objective

Integrate every computational level into one coherent framework.

Tasks

□ Cross-validation

□ Internal consistency

□ Literature comparison

□ Sensitivity analyses

□ Final scientific interpretation

Status

NOT STARTED

---

# Final deliverables

The project is expected to deliver

✓ Reproducible repository

✓ Complete computational workflow

✓ Molecular dynamics

✓ Embedded TDDFT

✓ Time-dependent site energies

✓ Excitonic Hamiltonians

✓ Quantum dynamics

✓ Optical response

✓ Electromagnetic response

✓ Documentation

✓ Scientific reports

Potential publications

The project architecture is intended to support future journal publications after complete validation.

---

# Current priorities

Priority 1

Complete force-field assessment.

Priority 2

Finish embedded TDDFT production.

Priority 3

Complete site-energy trajectory.

Priority 4

Complete diagonal Hamiltonians.

Priority 5

Begin electronic coupling calculations.

Priority 6

Start Phase 1B.

---

# Current completion estimate

Phase 0

██████████

100%

Phase 1A

████████░░

approximately 80%

Phase 1B

░░░░░░░░░░

0%

Phase 2

░░░░░░░░░░

0%

Phase 3

░░░░░░░░░░

0%

Phase 4

░░░░░░░░░░

0%

Phase 5

░░░░░░░░░░

0%

---

# Roadmap rule

No phase should begin before the scientific objectives of the previous phase have been satisfied.

Exceptions should only occur if explicitly justified by the Master Project Document.

The roadmap should be reviewed periodically to ensure that daily work remains aligned with the overall scientific objectives of the project.