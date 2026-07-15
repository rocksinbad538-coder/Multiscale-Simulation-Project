# QM_F06 Repaired-Fragment Final Audit — Day026

## Decision

**QM_F06_REPAIRED_FRAGMENTS_READY_FOR_QM_INPUT_PREPARATION**

## Gate logic

Hard long-range contacts involving artificial caps are blocking. Hard contacts composed only of original R2 atoms are retained as inherited conformational strain and must be monitored during QM optimization.

## QM_F06_LOWER_CAPPED_REPAIRED

- Excluded 1–2 pairs: **21**
- Excluded 1–3 pairs: **30**
- Excluded 1–4 pairs: **36**
- Long-range contacts below vdW sum: **10**
- Classification counts: `{'ACCEPTABLE': 4, 'CLOSE_CONTACT': 3, 'STRONG_COMPRESSION': 3}`
- Hard contacts inherited from R2: **3**
- Hard contacts involving artificial caps: **0**
- Repaired cap-distance failures: **0**
- Artificial-cap steric gate: **PASS**

## QM_F06_UPPER_CAPPED_REPAIRED

- Excluded 1–2 pairs: **21**
- Excluded 1–3 pairs: **30**
- Excluded 1–4 pairs: **36**
- Long-range contacts below vdW sum: **11**
- Classification counts: `{'ACCEPTABLE': 6, 'CLOSE_CONTACT': 3, 'SEVERE_CLASH': 1, 'STRONG_COMPRESSION': 1}`
- Hard contacts inherited from R2: **2**
- Hard contacts involving artificial caps: **0**
- Repaired cap-distance failures: **0**
- Artificial-cap steric gate: **PASS**

## Authorization state

- Artificial-cap geometry gate: **PASSED**
- Repaired X–H distance gate: **PASSED**
- QM input preparation: **AUTHORIZED**
- QM calculation execution: **NOT AUTHORIZED**

## Required next step

Prepare reproducible electronic-structure input files for the LOWER and UPPER repaired fragments if both boundary gates pass. Execution remains a separate authorization decision.

