# Day019 PYR2 S2 transition-density cube validation

## Cube structure

- Atoms: 26
- Grid: 40 Ã 40 Ã 40
- Values: 64000
- Voxel volume: 2.2456736481e-01 bohrÂ³
- Grid lengths: 8.378476 Ã 16.558477 Ã 14.228478 Ã

## Integral diagnostics

- Net transition charge: 5.6771959687e-04
- Integral of |rho_tr|: 1.0267740893e+00
- Integral of rho_trÂ²: 3.6240489718e-03
- Boundary maximum / global maximum: 1.7975850542e-07
- Boundary |rho| integral fraction: 6.9487079254e-08

## Dipole reconstruction

- ORCA transition dipole: (0.2289900000, 2.2126900000, -0.0018200000) au
- ORCA magnitude: 2.2245066869 au
- Best cube convention: `+centered_first_moment`
- Best cube dipole: (-0.1604205461, -1.5530069604, -0.0003329862) au
- Best cube magnitude: 1.5612704703 au
- Absolute cosine with ORCA: 0.9999994499
- Relative magnitude error: 29.814980%

## Pilot acceptance

- Net charge target |Q| <= 1.0e-03: PASS
- Dipole cosine target >= 0.99: PASS
- Dipole magnitude error <= 5.0%: REVIEW
- Boundary maximum ratio <= 1.0e-03: PASS
- Overall pilot status: REVIEW

## Interpretation boundary

The generic cube comment `Total electron density` is an ORCA header string and does not override the transition-density identity established by the `cistp02` filename and the `orca_plot` generation path. A production grid should only be frozen after the integral, dipole, and boundary diagnostics are reviewed together.
