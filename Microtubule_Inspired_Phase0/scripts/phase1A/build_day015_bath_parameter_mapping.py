from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path("runs/phase1A/day015_bath_parameter_mapping")
OUT.mkdir(parents=True, exist_ok=True)

kBT_eV = 8.617333262145e-5 * 300.0

data = [
    ["kBT_300K", kBT_eV * 1000, "meV", "Thermal energy at 300 K"],
    ["J23", 2.37, "meV", "Weak PYR2-PYR3 bottleneck coupling"],
    ["J34", 17.00, "meV", "Intermediate PYR3-PYR4 coupling"],
    ["J45", 33.00, "meV", "Dominant PYR4-PYR5 coupling"],
    ["gamma_phi_tested", 1.0, "ps^-1", "Phenomenological exciton-basis pure dephasing used in Lindblad scan"],
    ["k_down_tested_low", 0.1, "ps^-1", "Slow relaxation scan point"],
    ["k_down_tested_reference", 1.0, "ps^-1", "Reference relaxation scan point"],
    ["k_down_tested_high", 10.0, "ps^-1", "Fast relaxation scan point"],
]

df = pd.DataFrame(data, columns=["quantity", "value", "unit", "interpretation"])
df.to_csv(OUT / "day015_energy_rate_scale_mapping.csv", index=False)

# Convert couplings to characteristic coherent times tau ~ hbar/J
HBAR_eV_fs = 0.6582119514
couplings = {
    "J23": 2.37,
    "J34": 17.00,
    "J45": 33.00,
}

rows = []
for name, J_meV in couplings.items():
    J_eV = J_meV / 1000
    tau_fs = HBAR_eV_fs / J_eV
    period_fs = 2 * np.pi * HBAR_eV_fs / J_eV
    rows.append({
        "coupling": name,
        "J_meV": J_meV,
        "hbar_over_J_fs": tau_fs,
        "2pi_hbar_over_J_fs": period_fs,
        "J_over_kBT": J_meV / (kBT_eV * 1000),
    })

pd.DataFrame(rows).to_csv(OUT / "day015_coupling_timescale_mapping.csv", index=False)

summary = """# Day015 Bath-Parameter Mapping

## Purpose

Prepare the bridge between the current phenomenological open-system model and physically calibrated bath parameters.

## Context

The current Lindblad model is internally consistent and converges to the imposed 300 K thermal equilibrium. However, the dephasing and relaxation rates are still phenomenological. The next model-improvement step is to derive or constrain these rates from MD fluctuations and spectral-density analysis.

## Current physical scales

At 300 K:

- kBT = 25.85 meV
- J23 = 2.37 meV
- J34 = 17.00 meV
- J45 = 33.00 meV

Therefore:

- J23 << kBT
- J34 is comparable to kBT
- J45 is slightly larger than kBT

## Interpretation

The weak J23 coupling explains the PYR2 access bottleneck. The stronger J34 and J45 couplings support faster redistribution within the PYR3-PYR5 subnetwork.

The current Lindblad relaxation scan should be interpreted as a controlled phenomenological model. To make it physically predictive, the following quantities should be extracted from MD or calibrated from literature:

1. Site-energy fluctuation time series.
2. Coupling fluctuation time series.
3. Autocorrelation functions.
4. Spectral densities.
5. Reorganization energies.
6. Dephasing rates.
7. Relaxation rates satisfying detailed balance.

## Immediate next computational target

Build an MD-to-bath-parameter workflow template:

- read chromophore-resolved site-energy proxy time series;
- compute fluctuation autocorrelation C(t);
- estimate correlation time;
- estimate dephasing scale;
- generate a first spectral-density proxy.

This will connect the current Lindblad model to physically motivated bath parameters.
"""

(OUT / "BATH_PARAMETER_MAPPING_DAY015.md").write_text(summary)

print(df.to_string(index=False))
print()
print(pd.read_csv(OUT / "day015_coupling_timescale_mapping.csv").to_string(index=False))
print("Wrote:", OUT)
