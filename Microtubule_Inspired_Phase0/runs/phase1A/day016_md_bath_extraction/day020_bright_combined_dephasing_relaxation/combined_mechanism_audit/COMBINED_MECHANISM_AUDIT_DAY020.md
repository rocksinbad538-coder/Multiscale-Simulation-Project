# Day020 Combined-Mechanism Audit

## Numerical status

- All numerical validation criteria: PASS.
- Conditions audited: 24/24.

## Robust pathway result

- PYR4 gateway fraction: 0.842995742 to 0.914133516.
- PYR4 remains the dominant thermal gateway in every tested condition.

## Transient PYR5 population

- Full 100 ps final-population range: 0.000160392 to 0.110107150.
- At gamma_phi=100 ps^-1: 0.105185501 to 0.110107150.
- Maximum final population occurs at T=300 K, kappa=10, gamma=100.

## Stationary-state competition

- Thermal-only stationary PYR5 range: 0.999950643 to 0.999977286.
- Nonzero-dephasing stationary PYR5 range: 0.241116708 to 0.749992076.
- Overall stationary range: 0.241116708 to 0.999977286.
- Minimum stationary PYR5 occurs at T=150 K, kappa=10, gamma=100.

## Temperature sensitivity

- Maximum absolute temperature effect on the 100 ps PYR5 population: 0.004921649.
- Maximum absolute temperature effect on stationary PYR5: 0.042051749.

## Rate-ratio assessment

- At T=150 K and gamma/kappa=10, the final-PYR5 spread across absolute scales is 0.103528292, while the stationary spread is 0.023655126.
- At T=300 K and gamma/kappa=10, the final-PYR5 spread across absolute scales is 0.108449595, while the stationary spread is 0.013819357.

The ratio gamma_phi/kappa identifies the dominant dissipative channel but does not fully determine either the transient or stationary dynamics. Absolute rates relative to the coherent Hamiltonian remain relevant.

## Accepted interpretation

The upper PYR2-PYR4 manifold undergoes rapid redistribution, followed by kinetically slow capture into PYR5 predominantly through PYR4. Local dephasing enhances transient access to PYR5 but competes with thermal relaxation in the stationary state. Absolute relaxation times remain phenomenological because no microscopic spectral density has been derived.
