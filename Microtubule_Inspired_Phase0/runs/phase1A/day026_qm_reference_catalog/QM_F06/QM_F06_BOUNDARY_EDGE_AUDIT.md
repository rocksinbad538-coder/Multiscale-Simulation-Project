# QM_F06 Boundary-Edge Audit — Day026

## Scope

All graph edges cut during extraction of the LOWER and UPPER QM_F06 fragments were chemically classified.

No atoms were added, removed or moved. Hydrogen capping remains unauthorized pending this audit.

## Combined result

- Total cut edges: **16**
- Element-pair counts: `{('B', 'H'): 1, ('B', 'N'): 7, ('N', 'B'): 7, ('N', 'H'): 1}`
- Preliminary actions: `{'CANDIDATE_BN_CUT_FOR_HYDROGEN_CAPPING': 14, 'INCLUDE_EXISTING_HYDROGEN': 2}`
- Decision: **FRAGMENT_EXPANSION_REQUIRED_BEFORE_CAPPING**

## QM_F06_LOWER

- Boundary edges audited: **8**
- Element-pair counts: `{('B', 'N'): 4, ('N', 'B'): 3, ('N', 'H'): 1}`
- Preliminary actions: `{'CANDIDATE_BN_CUT_FOR_HYDROGEN_CAPPING': 7, 'INCLUDE_EXISTING_HYDROGEN': 1}`

## QM_F06_UPPER

- Boundary edges audited: **8**
- Element-pair counts: `{('B', 'H'): 1, ('B', 'N'): 3, ('N', 'B'): 4}`
- Preliminary actions: `{'CANDIDATE_BN_CUT_FOR_HYDROGEN_CAPPING': 7, 'INCLUDE_EXISTING_HYDROGEN': 1}`

## Authorization state

- Artificial capping authorized: **NO**
- Geometry optimization authorized: **NO**
- QM calculation executed: **NO**

## Required next step

Review every preliminary action and either enlarge the graph fragment or define chemically valid B–H/N–H cap placement along the original cut-bond vectors.
