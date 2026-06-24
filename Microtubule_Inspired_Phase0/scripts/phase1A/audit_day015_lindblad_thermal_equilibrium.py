from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path("runs/phase1A/day015_lindblad_relaxation")

KB = 8.617333262145e-5
T = 300.0
kBT = KB * T

labels_site = ["PYR2", "PYR3", "PYR4", "PYR5"]
labels_exc = ["X1", "X2", "X3", "X4"]

H = pd.read_csv(OUT / "lindblad_input_hamiltonian_site_basis.csv", index_col=0).values
E = pd.read_csv(OUT / "lindblad_exciton_energies.csv")["energy_eV"].values
U = pd.read_csv(OUT / "lindblad_exciton_eigenvectors.csv", index_col=0).values

E0 = E.min()
boltz = np.exp(-(E - E0) / kBT)
boltz = boltz / boltz.sum()

rho_exc_eq = np.diag(boltz)
rho_site_eq = U @ rho_exc_eq @ U.T
site_eq = np.diag(rho_site_eq)

rows = []
for i, lab in enumerate(labels_exc):
    rows.append({
        "basis": "exciton",
        "state": lab,
        "energy_eV": E[i],
        "thermal_population": boltz[i],
    })

for i, lab in enumerate(labels_site):
    rows.append({
        "basis": "site",
        "state": lab,
        "energy_eV": np.nan,
        "thermal_population": site_eq[i],
    })

pd.DataFrame(rows).to_csv(OUT / "lindblad_expected_thermal_equilibrium.csv", index=False)

summary = pd.read_csv(OUT / "lindblad_relaxation_summary.csv")

cmp_rows = []
for _, r in summary.iterrows():
    cmp_rows.append({
        "k_down_ps": r["k_down_ps"],
        "final_PYR2": r["final_PYR2"],
        "expected_PYR2": site_eq[0],
        "abs_error_PYR2": abs(r["final_PYR2"] - site_eq[0]),
        "final_PYR3": r["final_PYR3"],
        "expected_PYR3": site_eq[1],
        "abs_error_PYR3": abs(r["final_PYR3"] - site_eq[1]),
        "final_PYR4": r["final_PYR4"],
        "expected_PYR4": site_eq[2],
        "abs_error_PYR4": abs(r["final_PYR4"] - site_eq[2]),
        "final_PYR5": r["final_PYR5"],
        "expected_PYR5": site_eq[3],
        "abs_error_PYR5": abs(r["final_PYR5"] - site_eq[3]),
    })

cmp = pd.DataFrame(cmp_rows)
cmp.to_csv(OUT / "lindblad_thermal_equilibrium_audit.csv", index=False)

print("Expected exciton thermal populations:")
for lab, pop in zip(labels_exc, boltz):
    print(f"{lab}: {pop:.6f}")

print("\nExpected site populations:")
for lab, pop in zip(labels_site, site_eq):
    print(f"{lab}: {pop:.6f}")

print("\nFinal-vs-expected audit:")
print(cmp.to_string(index=False))
