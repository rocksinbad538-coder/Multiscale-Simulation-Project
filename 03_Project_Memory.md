# 03_PROJECT_MEMORY.md

# Multiscale Simulation – Microtubule Inspired Artificial System
## Project Memory

Version: Day 026

---

# Purpose

This document records the persistent knowledge accumulated during the project.

It intentionally stores decisions that should not be rediscovered.

Whenever development resumes after several days (or months), this document should be reviewed before writing code or launching simulations.

It is the "long-term memory" of the project.

---

# 1. General philosophy

The project is not a collection of isolated simulations.

Everything must contribute to one single scientific objective:

A fully reproducible multiscale simulation framework connecting

atomistic structure

↓

molecular dynamics

↓

electronic structure

↓

excitonic Hamiltonian

↓

open quantum dynamics

↓

optical response.

Every new script must support this workflow.

---

# 2. Guiding rule

Whenever multiple possible developments exist, always prioritize the one that moves the project closer to the final objectives described in the Master Project Document.

Avoid interesting but unnecessary calculations.

Scientific value is always preferred over computational complexity.

---

# 3. Repository philosophy

The repository should remain reproducible.

Manual workflows should be converted into scripts.

Every important calculation should leave

inputs

outputs

metadata

logs

summary tables

documentation.

No important result should exist only inside a terminal session.

---

# 4. Automation philosophy

Whenever a task is repeated twice,

automate it.

Whenever a parser is written,

make it reusable.

Whenever a calculation is validated,

convert it into production.

---

# 5. Molecular dynamics

MD is not the final objective.

It provides the structural ensemble from which electronic properties are extracted.

The project therefore uses MD snapshots rather than a single optimized geometry.

---

# 6. Local QM philosophy

Never perform TDDFT on the complete nanotube.

Always isolate

one chromophore

+

local water

+

electrostatic environment.

This was chosen because it provides a physically meaningful compromise between accuracy and computational cost.

---

# 7. Electrostatic embedding

The electrostatic embedding workflow has already been validated.

Previously solved problems include

point-charge ordering,

ORCA point-charge formatting,

automatic generation of embedding files,

automatic manifests,

automatic execution,

automatic parsing.

These issues should not be revisited unless new evidence appears.

---

# 8. Embedded TDDFT

The embedded TDDFT workflow is considered validated.

Successful production calculations already exist.

The remaining calculations represent computational work rather than methodological development.

Therefore,

future effort should focus on production throughput instead of redesigning the workflow.

---

# 9. Hamiltonian philosophy

Current Hamiltonians contain only diagonal terms.

This is intentional.

The diagonal energies must first be validated over the complete MD trajectory.

Only afterwards should electronic couplings be incorporated.

Never attempt to estimate couplings before completing the site-energy trajectory.

---

# 10. Scientific interpretation

Avoid interpreting isolated frames.

Scientific conclusions should only be drawn after analyzing the complete trajectory.

Single-frame observations are useful only as consistency checks.

---

# 11. Force-field philosophy

The project should not generate a new parametrization unless necessary.

Before creating new parameters,

always verify whether an existing validated parametrization can be reused.

The current work (Day 026) follows exactly this philosophy.

---

# 12. Current bottleneck

The limiting factor is no longer methodology.

The limiting factor is computational throughput.

The production workflow already works.

Remaining work consists primarily of executing and analyzing the remaining TDDFT calculations.

---

# 13. Interaction with Vitalii

Communication philosophy:

One concise Slack update every working hour.

Updates should describe

completed work,

current status,

next objective.

Avoid discussing speculative ideas.

Communicate completed scientific progress.

---

# 14. Daily workflow

Each workday begins with

1.

Slack message to Vitalii.

2.

Update Excel planning sheet.

3.

Review repository status.

4.

Review current phase.

5.

Resume technical work.

During the day

hourly Slack updates.

End of day

README updates if needed,

daily notes,

Git commit,

GitHub push,

technical summary.

---

# 15. Documentation philosophy

The documentation hierarchy is

Master Project Document

↓

README

↓

Daily Notes

↓

Generated Reports

↓

Automatic Metadata

This hierarchy should always remain consistent.

---

# 16. Coding philosophy

Prefer

small,

modular,

reusable,

well-documented

Python scripts.

Avoid large monolithic scripts.

Whenever practical,

one script should perform one scientific task.

---

# 17. Validation philosophy

Every workflow should be validated in three stages.

Pilot

↓

Controlled production

↓

Full production

Never jump directly to production without successful pilot validation.

---

# 18. Scientific priorities

Current priority order:

1.

Force-field decision.

2.

Remaining embedded TDDFT calculations.

3.

Complete site-energy trajectory.

4.

Complete diagonal Hamiltonian trajectory.

5.

Electronic couplings.

6.

Quantum dynamics.

7.

Optical calculations.

---

# 19. Important observations

Current recurring observations include

PYR5 consistently exhibits the lowest excitation energy.

Energy fluctuations remain moderate.

Hamiltonian traces remain stable.

No numerical instabilities have appeared.

These observations should continue to be monitored as production progresses.

---

# 20. Things already solved

Do NOT spend time re-solving

ORCA point-charge formatting.

Embedding input generation.

Automatic manifests.

Site-energy extraction.

Hamiltonian generation.

Batch execution.

Production bookkeeping.

These components are already operational.

---

# 21. Things still missing

Complete production trajectory.

Complete site-energy statistics.

Electronic couplings.

Time-dependent Hamiltonian.

Spectral densities.

Quantum dynamics.

Optical spectra.

Electromagnetic calculations.

Final validation against the scientific objectives of the project.

---

# 22. Golden rule

Whenever uncertainty appears, always ask:

"Does this task directly contribute to the final multiscale simulation framework promised in the Master Project Document?"

If the answer is no,

it should probably not be done.