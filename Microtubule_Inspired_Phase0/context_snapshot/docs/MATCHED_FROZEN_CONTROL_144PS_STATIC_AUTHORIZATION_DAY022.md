# Matched Frozen-Control 144 ps Static Authorization

## Control definition

- Initial coordinates:
  `runs/phase1A/day021_mobile_restraint_protocol/execution/01_em_k10000/01_em_k10000.gro`
- Reference mobile initial state:
  `runs/phase1A/day021_mobile_restraint_protocol/execution/02_nvt_k10000_1ps/02_nvt_k10000_1ps.tpr`
- HBN/PYR treatment: **frozen in X, Y, and Z**
- Position restraints: **0**
- Duration: **144.0 ps**
- Time step: **0.0005 ps**
- Expected trajectory frames: **289**

## Initial-state identity

- Coordinate exact match: **True**
- Coordinate maximum absolute difference:
  **0.000000000000 nm**
- Stage02 water-velocity SHA256:
  `d11ac84fa2d4a2a9a594a91ac0a6e0714dc3caa9da33f0ed12617371feff722d`
- Control water-velocity SHA256:
  `d11ac84fa2d4a2a9a594a91ac0a6e0714dc3caa9da33f0ed12617371feff722d`
- Water-velocity exact match:
  **True**
- Water-velocity maximum absolute difference:
  **0.000000000000 nm ps^-1**

## Decision

- Static decision: **PASS**
- Matched frozen-control execution authorized:
  **YES**
- Failure reasons:
  **NONE**

No MD execution was performed by this preparation workflow.
