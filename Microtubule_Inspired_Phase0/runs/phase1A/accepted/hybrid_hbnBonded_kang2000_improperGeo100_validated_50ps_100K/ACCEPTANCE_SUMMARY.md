# Phase 1A Accepted Scaffold Model

Accepted model:

- h-BN bonded scaffold
- angular stiffness: kang2000
- geometry-referenced improper stiffness: improperGeo100
- validation run: NVT, 50 ps, 100 K
- PYR and SOL frozen during scaffold-mechanics validation

## Acceptance criteria

The model is accepted as the current Phase 1A scaffold-mechanics baseline because:

- The 50 ps / 100 K NVT run completed without crash.
- No PYR-HBN contacts below 0.30 nm were observed.
- PYR4, previously the problematic chromophore, relaxed to 0.47406 nm minimum PYR-HBN distance.
- Mean radial deformation remained small: +0.01313 nm.
- The remaining PYR2/PYR3 distances around 0.336–0.342 nm are interpreted as localized adsorption-like contacts, not nanotube collapse.

## Key audit values

PYR-HBN minimum distances after 50 ps:

- PYR2: 0.33603 nm
- PYR3: 0.34234 nm
- PYR4: 0.47406 nm
- PYR5: 0.69622 nm

Radial deformation after 50 ps:

- min Δr: -0.17872 nm
- max Δr: +0.24672 nm
- mean Δr: +0.01313 nm

## Operational conclusion

This model is the accepted scaffold-mechanics baseline for the next Phase 1A stage.
