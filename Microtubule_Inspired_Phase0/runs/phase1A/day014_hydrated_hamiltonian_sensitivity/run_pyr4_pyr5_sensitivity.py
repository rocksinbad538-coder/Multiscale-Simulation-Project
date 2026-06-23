import numpy as np
import pandas as pd
from pathlib import Path

outdir = Path("runs/phase1A/day014_hydrated_hamiltonian_sensitivity")
basis = ["PYR2","PYR3","PYR4","PYR5"]

H_base = pd.read_csv("runs/phase1A/day013_hydrated_exciton_hamiltonian/hydrated_exciton_hamiltonian_4x4_eV.csv", index_col=0)

cases = {
    "baseline_J45_11meV": 0.011,
    "hydrated_J45_33meV": 0.033,
}

summary = []

for name, j45 in cases.items():
    H = H_base.copy()
    H.loc["PYR4","PYR5"] = j45
    H.loc["PYR5","PYR4"] = j45

    vals, vecs = np.linalg.eigh(H.values)
    mean_site = np.mean(np.diag(H.values))

    weights = vecs**2
    pr = 1.0 / np.sum(weights**2, axis=0)

    H.to_csv(outdir / f"{name}_hamiltonian_eV.csv")

    eig_rows = []
    vec_rows = []

    for i in range(4):
        state = f"X{i+1}"
        dominant = basis[int(np.argmax(weights[:, i]))]
        eig_rows.append({
            "case": name,
            "exciton_state": state,
            "energy_eV": vals[i],
            "shift_meV_from_mean_site": (vals[i] - mean_site)*1000,
            "participation_ratio": pr[i],
            "dominant_site": dominant,
            **{f"weight_{b}": weights[j, i] for j, b in enumerate(basis)}
        })
        vec_rows.append({"site": basis[i], **{f"X{k+1}": vecs[i,k] for k in range(4)}})

    pd.DataFrame(eig_rows).to_csv(outdir / f"{name}_eigenstates.csv", index=False)
    pd.DataFrame(vec_rows).to_csv(outdir / f"{name}_eigenvectors.csv", index=False)
    summary.extend(eig_rows)

df = pd.DataFrame(summary)
df.to_csv(outdir / "PYR4_PYR5_J45_sensitivity_summary.csv", index=False)

pivot = df.pivot(index="exciton_state", columns="case", values=["energy_eV","participation_ratio","dominant_site"])
print(df.to_string(index=False))
print()
print("Wrote:", outdir)
