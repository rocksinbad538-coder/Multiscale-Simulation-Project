# Phase 2 Multi-temperature Molecular Dynamics Campaign


## Simulation summary


- Number of simulations: 5
- Temperature range: 150–350 K
- Production length: 1 ns
- Force field: Phase1B
- Ensemble: NVT

## Comparative results


| Temperature_K | Final_Rg_A | Mean_RMSD_A | Final_PE | Mean_Temperature_K | Mean_RMSF_A | Shape_kappa2 | Aligned_RMSD_A |
|---|---|---|---|---|---|---|---|
| 150.000000 | 4.133315 | 5.351213 | 30.808597 | 150.173997 | 2.064502 | 0.257895 | 4.494224 |
| 200.000000 | 4.170315 | 5.659469 | 38.250705 | 199.853453 | 2.143727 | 0.253342 | 4.550294 |
| 250.000000 | 4.185672 | 5.770567 | 46.632321 | 249.714420 | 2.066536 | 0.258305 | 4.294216 |
| 300.000000 | 4.094927 | 5.674974 | 43.320885 | 299.611265 | 1.813078 | 0.256949 | 4.131270 |
| 350.000000 | 3.938007 | 6.011266 | 52.021063 | 350.118176 | 2.118131 | 0.259094 | 4.460192 |

## Preliminary observations


- Temperature control remained close to the target value for all simulations.
- Radius of gyration shows moderate temperature dependence.
- RMSD varies smoothly across temperatures.
- RMSF remains comparatively stable.
- Relative shape anisotropy remains nearly constant.

## Figures


- figures/Rg_vs_Temperature.png
- figures/Mean_RMSD_vs_Temperature.png
- figures/Mean_RMSF_vs_Temperature.png
- figures/PotentialEnergy_vs_Temperature.png
- figures/ShapeAnisotropy_vs_Temperature.png
- figures/AlignedRMSD_vs_Temperature.png