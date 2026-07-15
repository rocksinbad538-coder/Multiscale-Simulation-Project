# Day019 transition-dipole and geometry audit

## Validation

- Embedded ORCA outputs audited: 84/84
- S1/S2 transition-dipole records: 168/168
- QM geometry records: 84/84
- Frozen site geometries: 4/4 PASS
- Maximum D2 reconstruction error: 2.43145000e-05 auÂ² (tolerance 5.00000000e-05)
- Maximum oscillator-strength reconstruction error: 2.32073333e-06 (tolerance 2.00000000e-05)
- Maximum state/table energy discrepancy: 4.96000000e-04 eV

## Transition-dipole magnitudes

- Bright-like |mu| range: 2.064105 to 2.694929 au
- Alternate-like |mu| range: 0.087293 to 0.772748 au

## Site-pair geometry

| Site A | Site B | Carbon-centroid distance (Ã) | Minimum atom distance (Ã) | Minimum C-C distance (Ã) |
|---|---|---:|---:|---:|
| PYR2 | PYR3 | 26.896125 | 20.692322 | 22.224367 |
| PYR2 | PYR4 | 42.490545 | 38.445491 | 39.494165 |
| PYR2 | PYR5 | 46.893453 | 40.251493 | 42.233354 |
| PYR3 | PYR4 | 26.894353 | 20.686104 | 22.219496 |
| PYR3 | PYR5 | 43.421130 | 39.484764 | 40.512028 |
| PYR4 | PYR5 | 27.639866 | 21.356100 | 22.885884 |

## Coupling-model boundary

The audited electric transition dipoles and fixed site geometry are sufficient to construct a point transition-dipole coupling baseline for all 24 intersite state pairs per frame. That baseline is not automatically a final coupling model: its validity must be judged against chromophore size and separation, and it should be superseded by a transition-charge or transition-density treatment when short-range effects are material.
