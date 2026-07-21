# QM_F06 UPPER V4 Pipeline (Day030)

## Objective

Repair the defective V3 boundary while preserving the validated V3 core,
eliminating the unstable geminal cap topology, restoring canonical hBN
connectivity, and generating a QM-ready V4 model.

---

## Pipeline

### Phase 1 — Expansion analysis

audit_day030_qm_f06_upper_v4_expansion_scope.py

Purpose

• determine the minimum chemically necessary expansion
• identify newly exposed heavy-heavy cuts

Output

QM_F06_UPPER_V4_EXPANSION_SCOPE.json

Decision

Expansion required.

---

### Phase 2 — Boundary delta

audit_day030_qm_f06_upper_v4_boundary_delta.py

Purpose

Compare V3 and proposed V4 boundaries.

Output

QM_F06_UPPER_V4_BOUNDARY_DELTA.json

Decision

Two previous cuts repaired.
Six genuinely new cuts identified.

---

### Phase 3 — Second-shell audit

audit_day030_qm_f06_upper_v4_second_shell.py

Purpose

Determine whether a second heavy shell is sufficient.

Decision

Second shell still leaves external cuts.

---

### Phase 4 — Final closure

audit_day030_qm_f06_upper_v4_final_closure.py

Purpose

Verify direct B–H closure.

Decision

Topologically valid.

---

### Phase 5 — Selective expansion

audit_day030_qm_f06_upper_v4_p1523_selective_expansion.py

Purpose

Evaluate replacing a double-H cap by restoring P:1523.

Decision

Accepted.

Artificial H reduced by one.

---

### Phase 6 — Coordinate alignment

audit_day030_qm_f06_upper_v4_coordinate_alignment.py

Purpose

Compare raw coordinates.

Decision

Translation detected.

Rigid alignment required.

---

### Phase 7 — Rigid alignment

audit_day030_qm_f06_upper_v4_rigid_alignment.py

Purpose

Compute optimal rigid-body transform.

Decision

Large residuals originate only from the intentionally mobile region.

---

### Phase 8 — Alignment subset analysis

audit_day030_qm_f06_upper_v4_alignment_subsets.py

Purpose

Determine the optimal alignment anchor.

Decision

FIXED_SHARED_REAL selected.

---

### Phase 9 — Fixed-anchor transformation

finalize_day030_qm_f06_upper_v4_fixed_anchor.py

Purpose

Generate transformed restoration coordinates.

Decision

PASS.

---

### Phase 10 — V4 construction

build_day030_qm_f06_upper_v4.py

Output

46 atoms

15 B

11 N

20 H

Decision

PASS

---

### Phase 11 — Structural pre-QM audit

audit_day030_qm_f06_upper_v4_pre_qm.py

Checks

• connectivity

• valence

• cap geometry

• hard contacts

Decision

PASS

---

### Phase 12 — Constraint design

design_day030_qm_f06_upper_v4_constraints.py

Result

15 fixed atoms

31 mobile atoms

Decision

PASS

---

### Phase 13 — ORCA input generation

prepare_day030_qm_f06_upper_v4_orca_input.py

Result

Fresh SCF

46 atoms

PBE0-D4

Decision

PASS

---

### Phase 14 — Execution preflight

prepare_and_preflight_day030_qm_f06_upper_v4_execution.py

Decision

PASS

Execution directory generated.

---

### Phase 15 — ORCA execution

run_v4_supervised.sh

Status

Running.

---

## Current state

ORCA optimization is in progress.

No additional geometry modifications are permitted until
the optimization completes.
