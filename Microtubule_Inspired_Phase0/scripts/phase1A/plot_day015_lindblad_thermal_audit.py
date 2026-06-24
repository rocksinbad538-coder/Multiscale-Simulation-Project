from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

OUT = Path("runs/phase1A/day015_lindblad_relaxation")

eq = pd.read_csv(OUT / "lindblad_expected_thermal_equilibrium.csv")
audit = pd.read_csv(OUT / "lindblad_thermal_equilibrium_audit.csv")

site = eq[eq["basis"] == "site"]
exc = eq[eq["basis"] == "exciton"]

plt.figure(figsize=(6,4))
plt.bar(exc["state"], exc["thermal_population"])
plt.xlabel("Exciton state")
plt.ylabel("Thermal population at 300 K")
plt.tight_layout()
plt.savefig(OUT / "lindblad_expected_exciton_thermal_population.png", dpi=300)
plt.close()

plt.figure(figsize=(6,4))
plt.bar(site["state"], site["thermal_population"])
plt.xlabel("Site")
plt.ylabel("Thermal site population at 300 K")
plt.tight_layout()
plt.savefig(OUT / "lindblad_expected_site_thermal_population.png", dpi=300)
plt.close()

plt.figure(figsize=(6,4))
plt.plot(audit["k_down_ps"], audit["abs_error_PYR2"], "o-", label="PYR2")
plt.plot(audit["k_down_ps"], audit["abs_error_PYR3"], "o-", label="PYR3")
plt.plot(audit["k_down_ps"], audit["abs_error_PYR4"], "o-", label="PYR4")
plt.plot(audit["k_down_ps"], audit["abs_error_PYR5"], "o-", label="PYR5")
plt.xlabel("k_down (ps^-1)")
plt.ylabel("Absolute final-vs-thermal error")
plt.legend()
plt.tight_layout()
plt.savefig(OUT / "lindblad_thermal_equilibrium_error_vs_kdown.png", dpi=300)
plt.close()

summary = """# Day015 Lindblad Thermal-Equilibrium Audit

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
"""
(OUT / "LINDBLAD_THERMAL_EQUILIBRIUM_AUDIT_DAY015.md").write_text(summary)

print("Wrote Lindblad thermal-equilibrium audit figures and summary.")
