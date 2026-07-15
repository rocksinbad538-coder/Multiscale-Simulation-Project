# QM_F06 Chemically Capped Fragments — Day026

## Decision

The original 14-atom bridge fragments required only minimal boundary completion. Existing R2 hydrogen atoms were restored, and peripheral B-N cuts were saturated with artificial H caps along the original bond vectors.

No additional full coordination shell was required because none of the cut edges touched the bridge core or either attachment center.

## Capping geometry

- B-H target distance: **1.19 Å**
- N-H target distance: **1.01 Å**
- Placement: along each original inside-to-outside B-N bond vector.

## QM_F06_LOWER_CAPPED

- Original extracted atoms: **14**
- Existing R2 hydrogen atoms incorporated: **1**
- Artificial caps added: **7**
- Final atoms: **22**
- Element counts: `{'B': 5, 'H': 12, 'N': 5}`
- Cap-distance validation failures: **0**
- Geometry optimized: **NO**
- QM calculation executed: **NO**

## QM_F06_UPPER_CAPPED

- Original extracted atoms: **14**
- Existing R2 hydrogen atoms incorporated: **1**
- Artificial caps added: **7**
- Final atoms: **22**
- Element counts: `{'B': 5, 'H': 12, 'N': 5}`
- Cap-distance validation failures: **0**
- Geometry optimized: **NO**
- QM calculation executed: **NO**

## Authorization state

- Fragment construction: **COMPLETED**
- Artificial caps: **GEOMETRIC INITIAL GUESSES ONLY**
- Geometry optimization authorized: **NO**
- QM calculation executed: **NO**

## Required next step

Audit capped-fragment valence, atom count, minimum interatomic distances and net charge/multiplicity requirements before preparing electronic-structure inputs.
