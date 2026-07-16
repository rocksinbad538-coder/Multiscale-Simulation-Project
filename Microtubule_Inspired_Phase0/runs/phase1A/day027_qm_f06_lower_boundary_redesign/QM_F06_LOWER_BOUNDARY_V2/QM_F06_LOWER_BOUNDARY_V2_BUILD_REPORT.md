# QM_F06 LOWER Boundary V2 Build — Day027

## Correction

- Removed artificial cap: `HCAP:LOWER:07`
- Restored real atom: `A:LOWER:13:-3`
- Restored first-shell atoms: `A:LOWER:11:-3`, `A:LOWER:14:-2`
- Restored real passivant: `H4:LOWER:0016:0`
- New peripheral caps: **3**

## Result

- Total atoms: **28**
- Element counts: `{'N': 7, 'B': 6, 'H': 15}`
- Role counts: `{'ORIGINAL_FRAGMENT_ATOM': 14, 'ARTIFICIAL_BOUNDARY_CAP': 6, 'EXISTING_R2_HYDROGEN_ADDED': 1, 'REAL_R2_BOUNDARY_EXPANSION_ATOM': 4, 'ARTIFICIAL_BOUNDARY_CAP_V2': 3}`
- Geometry optimized: **NO**
- QM calculation executed: **NO**

## Coordinate policy

Atoms retained from V1 use the converged Stage-2 geometry. Newly restored R2 atoms use validated Day024 coordinates. New caps are placed along the original R2 cut-bond vectors.

## Required next step

Reconstruct connectivity and audit valence, bond distances, boundary contacts and bridge proximity before preparing a V2 QM optimization.

