# Day019 cube-grid and direct vacuum/embedding overlap audit

## Scope

- Parsed and numerically integrated all 48 selected NTO cubes, including ORCA molecular-orbital dataset-ID records.
- Verified the cube-grid and frozen-geometry identity for the 24 vacuum/embedding orbital comparisons available for PYR2 and PYR5.
- Computed sign-invariant normalized orbital overlaps without interpolation only where grids and atom records were identical.

## Integrity results

- Cubes parsed successfully: 48/48
- Finite positive orbital norms: 48/48
- Orbital-norm range: 0.9999752506 to 1.0000632503
- Direct vacuum/embedding comparisons: 24/24
- Identical grids: 24/24
- Identical frozen atom records: 24/24

## Direct orbital-shape results

- Absolute-overlap range across 24 orbitals: 0.75925303 to 0.99964404
- Tracked-pair geometric-mean overlap range: 0.97202136 to 0.99951137
- Alternate-pair geometric-mean overlap range: 0.83210328 to 0.99675213

## Pair-resolved comparison

| Frame | Site | Root | Tracked | Pair | Hole | Particle | Geometric mean | Minimum |
|---:|---|---:|---:|---|---:|---:|---:|---:|
| 3 | PYR5 | S1 | True | `52a->53a` | 0.98494947 | 0.96182357 | 0.97331784 | 0.96182357 |
| 3 | PYR5 | S2 | False | `52a->53a` | 0.92695109 | 0.77900980 | 0.84976702 | 0.77900980 |
| 3 | PYR5 | S2 | False | `51a->54a` | 0.92455845 | 0.77937191 | 0.84886682 | 0.77937191 |
| 5 | PYR2 | S1 | False | `52a->53a` | 0.91194350 | 0.75925303 | 0.83210328 | 0.75925303 |
| 5 | PYR2 | S1 | False | `51a->54a` | 0.91344425 | 0.76262650 | 0.83463572 | 0.76262650 |
| 5 | PYR2 | S2 | True | `52a->53a` | 0.97843696 | 0.96564783 | 0.97202136 | 0.96564783 |
| 5 | PYR5 | S1 | True | `52a->53a` | 0.99440839 | 0.98824872 | 0.99132377 | 0.98824872 |
| 5 | PYR5 | S2 | False | `52a->53a` | 0.96055057 | 0.89063927 | 0.92493462 | 0.89063927 |
| 5 | PYR5 | S2 | False | `51a->54a` | 0.95946683 | 0.89053745 | 0.92435985 | 0.89053745 |
| 13 | PYR5 | S1 | True | `52a->53a` | 0.99964404 | 0.99937872 | 0.99951137 | 0.99937872 |
| 13 | PYR5 | S2 | False | `52a->53a` | 0.99786445 | 0.99530854 | 0.99658568 | 0.99530854 |
| 13 | PYR5 | S2 | False | `51a->54a` | 0.99772224 | 0.99578297 | 0.99675213 | 0.99578297 |

## Acceptance

**Day019 cube-grid and direct-overlap audit: PASS.**

Cross-site comparisons are deliberately excluded here because PYR2-PYR5 occupy different Cartesian frames. Those comparisons require atom-based rigid alignment and field interpolation; direct voxel-wise overlap would be invalid.
