# Day019 bright TDC-derived atomic-charge coupling benchmark

## Method

- Input: four phase-aligned, sqrt(2)-normalized N80 bright transition-density cubes.
- Each voxel charge was assigned to its nearest atom (atom-centered Voronoi partition).
- A minimum-norm correction imposed exactly zero total transition charge and the independently printed ORCA transition dipole.
- Coulomb couplings were evaluated between the resulting 26 atom-centered transition charges on each monomer.
- This is a transition-density-derived atomic-charge (TDC-AC) model, not the exact double integral over the two continuous three-dimensional densities.

## Numerical validation

- Sites partitioned and constrained: 4/4
- Maximum independent point-dipole reproduction error: 2.837537e-08 meV
- Maximum RMS atomic-charge correction: 1.547956e-02 e
- Maximum absolute atomic-charge correction: 2.723205e-02 e

## Site charge validation

| Site | Root | Cube Q | Raw dipole error | Corrected dipole error | RMS dq (e) | Max |dq| (e) |
|---|---:|---:|---:|---:|---:|---:|
| PYR2 | S2 | +1.308e-04 | 1.4916% | 8.988e-15 | 8.767e-03 | 1.898e-02 |
| PYR3 | S2 | -1.397e-04 | 1.4113% | 4.902e-15 | 1.548e-02 | 2.723e-02 |
| PYR4 | S2 | +5.976e-06 | 1.4400% | 5.210e-16 | 2.483e-03 | 5.022e-03 |
| PYR5 | S1 | -1.033e-04 | 1.1281% | 6.527e-15 | 1.085e-02 | 2.315e-02 |

## Coupling comparison

| Pair | Point dipole (meV) | Raw Voronoi (meV) | Corrected TDC-AC (meV) | Delta (meV) | TDC-AC / point |
|---|---:|---:|---:|---:|---:|
| PYR2-PYR3 | -1.170805 | -1.220059 | -1.208873 | -0.038067 | +1.032514 |
| PYR2-PYR4 | +0.254942 | +0.260059 | +0.253160 | -0.001782 | +0.993009 |
| PYR2-PYR5 | +0.073324 | +0.076158 | +0.078609 | +0.005285 | +1.072071 |
| PYR3-PYR4 | +1.221499 | +1.268898 | +1.200548 | -0.020951 | +0.982848 |
| PYR3-PYR5 | +0.286177 | +0.291318 | +0.282232 | -0.003945 | +0.986214 |
| PYR4-PYR5 | -1.391671 | -1.440054 | -1.458079 | -0.066408 | +1.047718 |

## Aggregate deviations

- Maximum |TDC-AC - point|: 0.066408 meV
- Maximum absolute relative deviation: 7.2071%
- Mean absolute relative deviation, three nearest pairs: 3.2461%
- Mean absolute relative deviation, three distant pairs: 3.0950%

## Interpretation boundary

The corrected TDC-AC model reproduces the exact transition charge and dipole by construction, so differences from the point-dipole result originate from the atom-resolved finite extent of the transition density. Agreement at the more distant pairs supports the far-field limit. Deviations at neighboring pairs quantify finite-size corrections within this atom-centered discretization. Because the Voronoi partition is not unique, these values should be treated as a finite-size benchmark rather than as an exact continuous-density Coulomb integral. No dielectric screening is included.
