# Water-Depletion Provenance Audit

## Temporal provenance

- Existing frozen trajectory relative to the mobile branch point:
  **-100 to 0 ps**
- Stage08 mobile trajectory relative to the branch point:
  **44 to 144 ps**
- Temporally matched windows: **NO**
- Direct attribution of the water difference to solute mobility:
  **NOT SUPPORTED**

## Count continuity

- Accepted frozen endpoint / frozen trajectory last frame:
  23 / 23
- Stage07 endpoint / Stage08 first frame:
  0 / 0
- Stage08 endpoint / Stage08 last frame:
  0 / 0
- Geometry/count consistency: **PASS**

## Required control

A matched frozen continuation must begin from the same branch
coordinates and use the same water-velocity initialization as the
mobile protocol. It should span 144 ps. The scientifically matched
comparison window is frozen 44-144 ps versus mobile Stage08 0-100 ps.

No electronic snapshot selection or recalculation is authorized from
the current unmatched comparison.
