# R1 Authoritative Starting State

The authoritative starting state for the R1 fully capped positive
control is:

`runs/phase1A/day023_confinement_design/01_r0_t0_reference/r0_accepted_t0_hydrated_system.gro`

SHA256:

`3e2f207361765c7448099591664f17bcecd4e2f53c516520d2ebcfb512028754`

Source trajectory:

`runs/phase1A/accepted/hybrid_hydrated_gap45_pyr5shift_clean032_nvt_100ps_frozenSolute/nvt_100ps_frozenSolute.xtc`

Source time:

**0.0 ps**

Accepted starting hydration:

- lumen occupancy:
  **428 waters**
- maximum occupancy observed in the accepted R0 trajectory:
  **437 waters**
- t=0 fraction of trajectory maximum:
  **0.979405**
- endpoint occupancy:
  **23 waters**

This state preserves the accepted R0 atom ordering:

- HBN: atoms 1-1680
- PYR: atoms 1681-1784
- TIP4P/2005 water:
  atoms 1785-68320

R1 cap construction must not modify this file in place. A derived
structure must be created in a new R1 design directory and must retain
a provenance link to this SHA256.
