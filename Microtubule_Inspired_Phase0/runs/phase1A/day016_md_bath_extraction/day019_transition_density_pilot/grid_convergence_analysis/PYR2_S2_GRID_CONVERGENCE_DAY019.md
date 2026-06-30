# Day019 PYR2 S2 transition-density grid convergence

## Resolution study

| Grid | Values | Raw |mu| (au) | Raw/ORCA | sqrt(2)-scaled |mu| (au) | Scaled/ORCA | Empirical scale | Qtr | Boundary ratio |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 40^3 | 64000 | 1.561270470 | 0.701849723 | 2.207969874 | 0.992565397 | 1.424806432 | 5.677e-04 | 1.798e-07 |
| 80^3 | 512000 | 1.570542822 | 0.706017994 | 2.221082960 | 0.998460222 | 1.416394495 | -9.247e-05 | 1.672e-07 |
| 120^3 | 1728000 | 1.571985567 | 0.706666562 | 2.223123309 | 0.999377436 | 1.415094550 | 8.692e-05 | 1.048e-07 |

## Finest-grid validation

- Finest grid: 120^3
- Directional cosine with ORCA: 0.9999999755
- Raw cube/ORCA dipole ratio: 0.7066665618
- sqrt(2)-scaled cube/ORCA ratio: 0.9993774358
- sqrt(2)-scaled magnitude error: 0.062256%
- Empirical scale factor: 1.4150945496
- Empirical scale / sqrt(2): 1.0006229520
- Previous-to-finest raw-dipole change: 0.091779%

## Decision controls

- Net transition charge: PASS
- Boundary containment: PASS
- Dipole direction: PASS
- sqrt(2)-scaled magnitude within 2.0%: PASS
- Grid convergence within 0.5%: PASS
- Overall convergence status: PASS

## Interpretation boundary

The near-1/sqrt(2) raw dipole ratio is treated here as a normalization hypothesis, not as an assumed ORCA convention. It is accepted only if the ratio is resolution-independent and the sqrt(2)-scaled moment converges to the independently printed ORCA transition dipole. Before production transition-density couplings, the same normalization test must also pass for PYR3, PYR4, and PYR5. State-specific empirical dipole normalization remains the fallback because it guarantees the correct far-field limit of each transition density.
