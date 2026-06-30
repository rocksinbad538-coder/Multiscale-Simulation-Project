# Day019 bright finite-size-corrected Hamiltonians

## Construction

- Six static correction factors were defined as `frame000 TDC-AC / frame000 point-dipole`.
- Each factor was applied to the corresponding bright-bright point-dipole coupling in all 21 solvent frames.
- The four-state bright model is the primary finite-size-corrected coupling model.
- A full eight-state hybrid sensitivity model was also written; only its bright-bright elements are corrected. Mixed and alternate-state couplings remain the original point-dipole values and are not transition-density benchmarked.

## Pair correction factors

| Pair | Point frame000 (meV) | TDC-AC frame000 (meV) | Factor | Change | Sign preserved |
|---|---:|---:|---:|---:|---|
| PYR2-PYR3 | -1.170805 | -1.208873 | 1.032514 | +3.251% | yes |
| PYR2-PYR4 | +0.254942 | +0.253160 | 0.993009 | -0.699% | yes |
| PYR2-PYR5 | +0.073324 | +0.078609 | 1.072071 | +7.207% | yes |
| PYR3-PYR4 | +1.221499 | +1.200548 | 0.982848 | -1.715% | yes |
| PYR3-PYR5 | +0.286177 | +0.282232 | 0.986214 | -1.379% | yes |
| PYR4-PYR5 | -1.391671 | -1.458079 | 1.047718 | +4.772% | yes |

## Dynamic impact across 21 frames

- Hamiltonian snapshots generated: 21/21
- Maximum |coupling correction|: 0.067962 meV
- Maximum full-eight-state eigenvalue shift: 0.042075 meV
- Maximum four-state bright eigenvalue shift: 0.042074 meV
- Minimum matched bright eigenvector overlap: 0.9998556977
- Maximum correction / maximum diagonal-energy SD: 0.004239
- Maximum correction / minimum local S1-S2 gap: 0.001282

## Pair-class summary

- Mean absolute factor change, nearest pairs: 3.2461%
- Mean absolute factor change, distant pairs: 3.0950%

## Interpretation boundary

The pair-specific factors transfer a static frame000 finite-size correction to all solvent frames. This is consistent with the frozen chromophore geometries and the observed stability of the bright transition-dipole directions, but it does not constitute an embedded transition-density calculation for every frame. The four-state bright Hamiltonian is therefore the defensible corrected control model. The eight-state hybrid is retained only for sensitivity analysis until alternate-state transition densities and their gauge are treated explicitly. No dielectric screening is included.
