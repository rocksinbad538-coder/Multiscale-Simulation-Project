# Day019 cross-site NTO aligned-subspace analysis

## Method

- The four vacuum-reference chromophores were compared pairwise.
- Moving geometries were aligned to reference geometries by a proper-rotation Kabsch fit over the 16 carbon atoms.
- Moving orbital fields were evaluated on the reference cube grid by trilinear interpolation.
- Both directions were calculated for every site pair and averaged to expose interpolation asymmetry.
- The tracked bright state was compared as a one-dimensional hole/particle pair.
- The alternate low state was compared as separate two-dimensional hole and particle subspaces. This is invariant to arbitrary rotations between near-degenerate NTO pairs.

## Numerical controls

- Pairwise site comparisons: 6/6
- Directional calculations: 12/12
- Minimum interpolated-orbital captured norm: 0.98688434
- Minimum reference-grid coverage fraction: 0.95077930
- Maximum directional similarity asymmetry: 0.00039128

## Symmetric cross-site results

| Site A | Site B | Heavy RMSD (Ã) | Tracked hole | Tracked particle | Tracked pair | Alternate hole subspace | Alternate particle subspace | Alternate transition subspace | Minimum principal cosine |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| PYR2 | PYR3 | 0.007352 | 0.99881719 | 0.99933308 | 0.99907510 | 0.99963714 | 0.99962106 | 0.99962910 | 0.99944304 |
| PYR2 | PYR4 | 0.008186 | 0.99773430 | 0.99930534 | 0.99851951 | 0.99961113 | 0.99964719 | 0.99962916 | 0.99955728 |
| PYR2 | PYR5 | 0.090479 | 0.99292263 | 0.98938616 | 0.99115282 | 0.99198938 | 0.98960231 | 0.99079513 | 0.98807321 |
| PYR3 | PYR4 | 0.006289 | 0.99966844 | 0.99986163 | 0.99976503 | 0.99981239 | 0.99985247 | 0.99983243 | 0.99972029 |
| PYR3 | PYR5 | 0.089509 | 0.99333905 | 0.99119032 | 0.99226410 | 0.99218371 | 0.98978644 | 0.99098435 | 0.98812104 |
| PYR4 | PYR5 | 0.089396 | 0.99264394 | 0.99126618 | 0.99195482 | 0.99227303 | 0.98995237 | 0.99111202 | 0.98838416 |

## Aggregate ranges

- Tracked-pair cross-site similarity range: 0.99115282 to 0.99976503
- Alternate-transition-subspace similarity range: 0.99079513 to 0.99983243

## Interpretation boundary

This analysis establishes spatial similarity after rigid alignment. It does not by itself provide diabatic state phases or interstate couplings. The tracked-root and alternate-root spaces should only be frozen into a Hamiltonian after these overlap results are reviewed together with the S1-S2 energy splittings, oscillator strengths, and subsequent coupling calculations.
