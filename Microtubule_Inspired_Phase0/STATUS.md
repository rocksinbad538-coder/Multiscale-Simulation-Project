# Microtubule Inspired Phase 0 — Current Status

Last updated: 2026-07-27

---

## Current active model

**QM_F06 UPPER V7A**

Status:

- ORCA geometry optimization running
- Formal construction completed
- Global pre-QM audit passed
- ORCA input preparation completed
- ORCA execution authorized and in progress

---

## Previous model

### QM_F06 UPPER V6B

Status:

Rejected after post-QM structural audit.

Reason for rejection:

- artificial boundary cap migrated
- canonical B-H bond lost
- noncanonical H-S interaction formed
- secondary B-P reconnectivity detected
- two atoms became overcoordinated

Consequences:

- RESP prohibited
- force-field generation prohibited
- molecular dynamics prohibited

---

## Current workflow

Completed:

- canonical topology repair
- V7A formal construction
- boundary provenance audit
- coordinate preflight
- cap generation
- global structural pre-QM audit
- ORCA input preparation
- execution packaging
- supervised ORCA execution

Running:

- ORCA geometry optimization (V7A)

Pending:

- post-QM structural audit
- RESP authorization review
- RESP calculation
- force-field parameter generation
- molecular dynamics

---

## Reproducibility policy

This repository preserves complete scientific provenance.

Successful and unsuccessful intermediate models are intentionally retained, together with:

- generating scripts
- intermediate reports
- audits
- ORCA inputs
- ORCA outputs
- execution manifests
- decision reports

No intermediate scientific artifacts are intentionally removed.

