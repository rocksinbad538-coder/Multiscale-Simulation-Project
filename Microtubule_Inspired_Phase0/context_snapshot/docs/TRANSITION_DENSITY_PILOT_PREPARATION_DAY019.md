# Day019 transition-density pilot preparation

## Scope

- Prepared the four frame000 bright-like monomer states.
- Verified matching ORCA `.gbw` and `.cis` files.
- Created standardized per-site symlinks named `pilot.gbw` and `pilot.cis`.
- Ranked all six site pairs by molecular-size/separation geometry for transition-density benchmarking.

## Pilot states

| Site | Bright root | GBW bytes | CIS bytes | All-atom radius (A) | Status |
|---|---:|---:|---:|---:|---|
| PYR2 | S2 | 2507310 | 654596 | 4.604400 | READY |
| PYR3 | S2 | 2507310 | 654596 | 4.605037 | READY |
| PYR4 | S2 | 2507310 | 654596 | 4.604559 | READY |
| PYR5 | S1 | 2507310 | 654596 | 4.756961 | READY |

## Geometry diagnostic

| Rank | Pair | Centroid distance (A) | Minimum atom distance (A) | (a_i+a_j)/R |
|---:|---|---:|---:|---:|
| 1 | PYR3-PYR4 | 26.894353 | 20.686104 | 0.342436 |
| 2 | PYR2-PYR3 | 26.896125 | 20.692322 | 0.342408 |
| 3 | PYR4-PYR5 | 27.639866 | 21.356100 | 0.338696 |
| 4 | PYR2-PYR4 | 42.490545 | 38.445491 | 0.216730 |
| 5 | PYR3-PYR5 | 43.421130 | 39.484764 | 0.215609 |
| 6 | PYR2-PYR5 | 46.893453 | 40.251493 | 0.199630 |

## Priority pairs

- PYR3-PYR4: R=26.894 A, minimum atom distance=20.686 A.
- PYR2-PYR3: R=26.896 A, minimum atom distance=20.692 A.
- PYR4-PYR5: R=27.640 A, minimum atom distance=21.356 A.

## Next controlled action

Run `launch_orca_plot_probe.sh PYR2` and select the CIS/TD-DFT transition-density option. The purpose of the first probe is to capture the exact ORCA 6.1.1 interactive prompt sequence before automating all four cubes. Production generation must use a single documented grid convention and must preserve the transition-density sign by matching the cube-derived dipole to the audited ORCA transition dipole.
