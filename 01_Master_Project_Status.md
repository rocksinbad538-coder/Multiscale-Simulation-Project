# Multiscale Simulation – Microtubule Inspired Artificial System
# MASTER PROJECT STATUS
Version: Day 026
Status: ACTIVE DEVELOPMENT

---

# 1. Project objective

The objective of this project is to develop a complete multiscale simulation framework for a microtubule-inspired artificial system capable of describing:

- atomistic structure;
- thermal stability;
- confined water;
- chromophore environment;
- excited electronic states;
- excitonic Hamiltonians;
- open quantum dynamics;
- optical response;
- THz / microwave response.

The entire workflow follows the roadmap established in the Project Master Document and associated supporting documentation.

The philosophy adopted throughout the project is:

MD
→ Local QM
→ TDDFT
→ Excitonic Hamiltonian
→ Quantum Dynamics
→ Optical Response

No stage should be implemented unless it directly contributes to the final deliverables promised in the project proposal.

---

# 2. Global project status

Current overall status:

Phase 0
COMPLETED

Phase 1A
ACTIVE

Phase 1B
NOT STARTED

Phase 2
NOT STARTED

Phase 3
NOT STARTED

Phase 4
NOT STARTED

Phase 5
NOT STARTED

Overall completion:
approximately Phase 1A.

The project is currently focused on establishing a fully validated production workflow for extracting site energies from molecular dynamics snapshots using embedded TDDFT calculations.

---

# 3. Scientific philosophy adopted

Several architectural alternatives were evaluated during the first weeks of the project.

The final decision was:

✓ chemically realistic atomistic scaffold

instead of

✗ simplified phenomenological scaffold.

Specifically:

• h-BN nanotube adopted as structural analogue.
• Explicit chromophores.
• Explicit confined water.
• Explicit molecular dynamics.
• Explicit electrostatic embedding.
• Excited states computed from TDDFT.

This decision aligns the project with publishable computational chemistry standards instead of a toy model.

---

# 4. Current scientific workflow

The current validated workflow is

MD trajectory

↓

Frame extraction

↓

Chromophore identification

↓

Local QM cluster extraction

↓

Electrostatic environment generation

↓

ORCA embedded TDDFT

↓

Site energies

↓

Diagonal excitonic Hamiltonian

↓

(time-dependent trajectory)

↓

Future quantum dynamics

Every block above already exists in the repository.

Only the production calculations are still running.

---

# 5. Completed milestones

The following milestones have been completed.

## Molecular system

✓ Functionalized nanotube constructed.

✓ Chromophore positions validated.

✓ Force-field implementation established.

✓ MD production trajectory generated.

---

## MD analysis

✓ MD trajectory extracted.

✓ Representative frames identified.

✓ Automatic frame extraction pipeline implemented.

---

## Chromophore mapping

✓ Automatic chromophore identification.

✓ Residue mapping.

✓ Model numbering.

✓ Distance validation.

---

## Local QM clusters

Automatic extraction implemented.

Each cluster contains

chromophore

+

nearby water

+

electrostatic environment.

---

## Electrostatic embedding

Validated.

Automatic generation of

.xyz

.pc

.inp

files.

Automatic ORCA directories.

Automatic manifests.

Automatic bookkeeping.

---

## Embedded TDDFT

Production workflow validated.

Pilot calculations:

Frames

000

010

020

completed successfully.

Later production expanded to

001

002

003

All calculations:

Normal ORCA termination.

TDDFT finished.

Automatic parsing.

Automatic validation.

---

## Site energies

Automatic extraction implemented.

Current outputs include

raw S1 energies

centered energies

global-centered energies

trajectory files

statistics

Hamiltonian diagonals

All automatically generated.

---

## Hamiltonian generation

Automatic generation implemented.

Outputs include

Hdiag_frameXXX.csv

Hdiag_frameXXX.npy

summary tables

metadata

The current Hamiltonians contain

diagonal elements only.

Off-diagonal couplings will be introduced in future phases.

---

# 6. Current production status

Current MD frames:

21

(frame000 through frame020)

Chromophores per frame:

4

(PYR2

PYR3

PYR4

PYR5)

Total embedded TDDFT jobs planned:

84

Completed:

24

Remaining:

60

The production workflow has already demonstrated that the calculations can be executed automatically without manual intervention.

This is considered the principal validation milestone before scaling to the remaining trajectory.

---

# 7. Current observations

Several physically meaningful observations have already emerged.

PYR5 consistently exhibits lower excitation energy than the remaining chromophores.

The energetic separation remains approximately stable along the sampled trajectory.

Frame-to-frame fluctuations remain modest.

Hamiltonian diagonal ranges remain physically reasonable.

No numerical instabilities have been detected.

These observations are preliminary.

Scientific interpretation will only be performed after the complete trajectory has been processed.

---

# 8. Current work (Day 026)

The current objective is NOT to generate more physics.

The objective is to finalize the computational infrastructure required for large-scale production.

Specifically,

complete comparison of candidate force fields,

evaluate whether one previously validated parametrization can be reused,

produce an evidence-based recommendation,

avoid unnecessary parametrization work,

prepare the remaining production calculations.

This comparison is expected to determine the force-field strategy adopted for the remainder of Phase 1A.

---

# 9. Immediate next milestones

Immediate priorities are

1.
Complete force-field comparison matrix.

2.
Finalize parametrization decision.

3.
Launch remaining embedded TDDFT production.

4.
Extract complete site-energy trajectory.

5.
Generate complete diagonal Hamiltonian trajectory.

Only after these are completed should the project move toward excitonic coupling calculations.

---

# 10. Medium-term roadmap

After completion of Phase 1A the project will proceed toward

electronic couplings,

time-dependent excitonic Hamiltonian,

spectral densities,

open quantum dynamics,

optical spectra,

electromagnetic response,

comparison with experimental observables.

---

# 11. Long-term roadmap

The final project should deliver a fully reproducible multiscale simulation pipeline capable of connecting

atomistic structure

↓

molecular dynamics

↓

electronic structure

↓

excitonic physics

↓

quantum dynamics

↓

optical response

within a single coherent computational framework.

This remains fully aligned with the Master Project Document.

---

# 12. Repository status

Repository:

active

Git:

clean

GitHub:

synchronized

Automation:

high

Reproducibility:

high

Documentation:

active

Daily logs:

maintained

Scientific traceability:

maintained

---

# 13. General status assessment

Project health:

GOOD

Scientific risk:

LOW

Computational workflow:

VALIDATED

Production calculations:

ONGOING

Current bottleneck:

Computational throughput (remaining ORCA calculations), not methodology.

No scientific blockers are currently identified.