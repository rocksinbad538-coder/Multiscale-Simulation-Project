# Day019 bright transition-density N80 validation

## Production convention

- Frame: 000
- Grid: 80 x 80 x 80
- Bright roots: PYR2=S2, PYR3=S2, PYR4=S2, PYR5=S1
- Each raw transition density is gauge-aligned to the independently printed ORCA transition dipole.
- The aligned density is multiplied by sqrt(2), following the normalization established by the PYR2 grid-convergence study.

## Site validation

| Site | Root | Gauge sign | Raw |mu| (au) | ORCA |mu| (au) | sqrt(2)-scaled |mu| (au) | Scaled error | Cosine | Qtr | Boundary ratio | Status |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| PYR2 | S2 | -1 | 1.570542822 | 2.224508208 | 2.221082960 | 0.1540% | 0.999996741 | -9.247e-05 | 1.672e-07 | PASS |
| PYR3 | S2 | -1 | 1.586682848 | 2.246803007 | 2.243908403 | 0.1288% | 0.999996662 | 9.881e-05 | 1.665e-07 | PASS |
| PYR4 | S2 | -1 | 1.583542244 | 2.240876457 | 2.239466918 | 0.0629% | 0.999997224 | -4.211e-06 | 1.821e-07 | PASS |
| PYR5 | S1 | -1 | 1.877980600 | 2.655189142 | 2.655865634 | 0.0255% | 0.999999930 | 7.307e-05 | 2.907e-07 | PASS |

## Global controls

- Valid sites: 4/4
- Maximum sqrt(2)-normalized dipole error: 0.153978%
- Minimum directional cosine: 0.9999966616
- Maximum |raw transition charge|: 9.881414e-05
- Maximum boundary ratio: 2.906788e-07
- Empirical scale range: 1.413853339 to 1.416394495
- Overall production status: PASS

## Interpretation boundary

Passing this validation establishes a common, phase-consistent and far-field-normalized bright-state transition-density set. It does not by itself establish dielectric screening. The first transition-density coupling benchmark should therefore be computed as an unscreened Coulomb integral and compared directly with the unscreened point-transition-dipole baseline.
