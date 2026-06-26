# Day017 production QC summary

## Objective

Promote the validated Day016 electrostatic-embedding TDDFT workflow from a controlled pilot to a scalable production workflow for MD-derived site-energy extraction.

## Production status

- Full extracted MD set: 21 frames
- Chromophores per frame: 4
- Planned embedded TDDFT site-energy jobs: 84
- Completed jobs so far: 24
- Completed frames: 000, 001, 002, 003, 010, 020
- Completion rate among launched jobs: 24/24 = 100%
- Completion rate relative to full production set: 24/84 = 28.6%

## Workflow advances

1. Reconstructed a full-production ORCA manifest from local QM cluster XYZ files.
2. Expanded the embedding input set from 12 pilot calculations to 84 production inputs.
3. Preserved the original pilot manifest as backup.
4. Implemented a resumable ORCA batch runner that skips completed calculations and only launches missing jobs.
5. Successfully executed the first controlled production batch beyond the pilot: frames 001–003.
6. Rebuilt the site-energy trajectory and diagonal excitonic Hamiltonians using all completed frames.

## Completed site-energy trajectory

| Frame | PYR2 | PYR3 | PYR4 | PYR5 |
|---:|---:|---:|---:|---:|
| 000 | 3.999 | 4.000 | 4.013 | 3.779 |
| 001 | 3.994 | 4.010 | 4.017 | 3.756 |
| 002 | 4.009 | 4.013 | 4.021 | 3.779 |
| 003 | 4.015 | 4.004 | 4.019 | 3.750 |
| 010 | 3.998 | 3.996 | 4.020 | 3.791 |
| 020 | 4.010 | 4.000 | 4.016 | 3.768 |

## Main physical observation

The first completed production subset confirms the qualitative pattern observed in the pilot: PYR2–PYR4 remain clustered near 4.00–4.02 eV, while PYR5 remains consistently red-shifted by approximately 0.21–0.27 eV relative to the higher-energy group.

## Hamiltonian status

Diagonal Hamiltonians were regenerated for all completed frames:

- `Hdiag_frame000`
- `Hdiag_frame001`
- `Hdiag_frame002`
- `Hdiag_frame003`
- `Hdiag_frame010`
- `Hdiag_frame020`

The current diagonal energy spread per frame ranges from 229 to 269 meV.

## Next step

Do not launch additional ORCA calculations until this production subset is committed. The next production step should continue with the remaining frames using the resumable runner, followed by automatic extraction of S1 site energies and regeneration of the complete diagonal Hamiltonian trajectory.
