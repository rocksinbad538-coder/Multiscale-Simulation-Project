# QM_F06 LOWER Final Acceptance — Day028

## Accepted model

Fragment:

`QM_F06_LOWER_BOUNDARY_V2B`

Accepted geometry source:

`runs/phase1A/day028_qm_f06_lower_boundary_v2b_executions/v2b_20260717T101410/v2b_partial_relax.xyz`

Formula:

`B6N7H15`

Charge and multiplicity:

`0 / 1`

Electronic-structure level:

`PBE0-D4/def2-TZVP`, with `RIJCOSX/def2-J`

## Structural validation

The accepted V2-B geometry preserves:

- all 27 intended covalent edges;
- complete graph connectivity;
- expected atomic valences;
- valid B–N, B–H and N–H bond distances;
- all artificial-cap bonds;
- the intended alternating B–N–B–N bridge topology.

No unintended covalent contacts or hard bridge–cap contacts were detected.

## Boundary-artifact resolution

The artificial boundary contact that motivated Boundary V2 was removed by:

1. replacing the problematic artificial termination with real R2
   coordination atoms;
2. performing a constrained V2-A boundary relaxation;
3. jointly relaxing all hydrogens and restored boundary atoms in V2-B.

The original blocking cap-related contact no longer remains a hard
boundary artifact.

## Residual real 1–5 contact

Pair:

`BR4:LOWER:00:3 — H4:LOWER:0017:0`

Properties:

- graph separation: 4 bonds;
- V2-B distance: 2.088349 Å;
- distance/vdW-sum ratio: 0.669343;
- artificial-cap involvement: NO;
- unintended covalent-distance criterion: NO.

Electronic diagnostic:

- target Mayer bond order: `< 0.01`;
- intended local B–H Mayer bond order: `0.9964`;
- intended local N–H Mayer bond order: `0.8695`.

Decision:

The residual contact has no significant covalent-bond character and does
not invalidate the accepted LOWER geometry.

## Diagnostic charge analyses

Target-atom charges:

| Scheme | B5 | H20 |
|---|---:|---:|
| Hirshfeld | +0.106987 | +0.140788 |
| MBIS | +0.979552 | +0.408799 |
| CHELPG | +0.501677 | +0.303981 |

All three charge sums reproduce the neutral molecular charge within
numerical precision.

These values are diagnostic only. Their substantial scheme dependence
precludes direct force-field charge adoption without a separately
defined charge-fitting protocol.

## Final decision

`QM_F06_LOWER_GEOMETRIC_AND_ELECTRONIC_REFERENCE_ACCEPTED`

## Authorization state

- LOWER geometry acceptance: **AUTHORIZED**
- LOWER electronic diagnostic: **COMPLETED**
- ESP/RESP protocol definition: **AUTHORIZED**
- ESP/RESP execution: **NOT AUTHORIZED**
- Charge adoption: **NOT AUTHORIZED**
- Bonded-parameter fitting: **NOT AUTHORIZED**
- Force-field parameter adoption: **NOT AUTHORIZED**
- MD execution: **NOT AUTHORIZED**
