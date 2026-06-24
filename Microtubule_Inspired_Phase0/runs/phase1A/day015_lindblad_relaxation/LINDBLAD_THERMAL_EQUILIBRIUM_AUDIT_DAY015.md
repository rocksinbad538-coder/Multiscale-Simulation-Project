# Day015 Lindblad Thermal-Equilibrium Audit

## Purpose

Validate that the exciton-basis Lindblad relaxation model converges to the intended 300 K thermal equilibrium distribution.

## Result

The expected thermal exciton populations are:

- X1: 0.536772
- X2: 0.306374
- X3: 0.126535
- X4: 0.030319

Projected into the site basis, this gives:

- PYR2: 0.305432
- PYR3: 0.160919
- PYR4: 0.283158
- PYR5: 0.250492

For k_down >= 1 ps^-1, the simulated final site populations match the expected thermal populations with errors near 1e-9 or smaller.

## Interpretation

This confirms that the Lindblad detailed-balance implementation is internally consistent. Unlike the pure-dephasing Haken-Strobl model, the long-time Lindblad state is not a uniform-population artifact; it corresponds to the imposed thermal equilibrium in the exciton basis.

## Caveat

The equilibrium structure is physically meaningful within the current Hamiltonian and temperature assumptions. The absolute relaxation rates remain phenomenological and should later be calibrated using MD-derived fluctuations, spectral-density modeling, or literature-based timescales.

## Files

- `lindblad_expected_thermal_equilibrium.csv`
- `lindblad_thermal_equilibrium_audit.csv`
- `lindblad_expected_exciton_thermal_population.png`
- `lindblad_expected_site_thermal_population.png`
- `lindblad_thermal_equilibrium_error_vs_kdown.png`
