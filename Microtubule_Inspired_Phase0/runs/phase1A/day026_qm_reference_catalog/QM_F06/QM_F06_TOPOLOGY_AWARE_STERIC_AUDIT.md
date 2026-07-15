# QM_F06 Topology-Aware Steric Audit — Day026

## Correction to the previous audit

The previous pair-specific audit incorrectly applied a van der Waals overlap criterion to bonded 1–2, angular 1–3 and torsional 1–4 intramolecular pairs.

The present audit excludes all pairs with graph separation ≤3 and applies the hard steric gate only to longer-range or disconnected atom pairs.

## Decision: **QM_F06_CAPPED_FRAGMENTS_REQUIRE_TARGETED_GEOMETRY_REPAIR**

## QM_F06_LOWER_CAPPED

- Excluded 1–2 pairs: **21**
- Excluded 1–3 pairs: **30**
- Excluded 1–4 pairs: **36**
- Long-range contacts below vdW sum: **16**
- Classification counts: `{'ACCEPTABLE': 6, 'CLOSE_CONTACT': 5, 'SEVERE_CLASH': 1, 'STRONG_COMPRESSION': 4}`
- Hard steric failures: **5**
- Hard failures involving artificial caps: **2**
- Topology-aware steric gate: **FAIL**

## QM_F06_UPPER_CAPPED

- Excluded 1–2 pairs: **21**
- Excluded 1–3 pairs: **30**
- Excluded 1–4 pairs: **36**
- Long-range contacts below vdW sum: **14**
- Classification counts: `{'ACCEPTABLE': 6, 'CLOSE_CONTACT': 3, 'SEVERE_CLASH': 2, 'STRONG_COMPRESSION': 3}`
- Hard steric failures: **5**
- Hard failures involving artificial caps: **3**
- Topology-aware steric gate: **FAIL**

## Authorization state

- Graph and valence gate: **PASSED**
- Topology-aware steric gate: **FAILED**
- QM input preparation: **NOT AUTHORIZED**
- QM calculation executed: **NO**

