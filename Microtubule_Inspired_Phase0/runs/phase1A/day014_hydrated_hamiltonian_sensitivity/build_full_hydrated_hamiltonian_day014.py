from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path("runs/phase1A/day014_hydrated_hamiltonian_sensitivity")
OUT.mkdir(parents=True, exist_ok=True)

# Site energies from the current Phase 1A working model.
# Keep these aligned with the Day013/Day014 hydrated Hamiltonian convention.
sites = {
    "PYR2": 3.755,
    "PYR3": 3.7775,
    "PYR4": 3.7775,
    "PYR5": 3.7775,
}

# Hydrated nearest-neighbor effective couplings from Day014 dimer analysis.
# Units: eV
couplings = {
    ("PYR2", "PYR3"): 0.00237,
    ("PYR3", "PYR4"): 0.01700,
    ("PYR4", "PYR5"): 0.03300,
}

labels = list(sites)
H = np.zeros((len(labels), len(labels)))

for i, a in enumerate(labels):
    H[i, i] = sites[a]

for (a, b), J in couplings.items():
    i = labels.index(a)
    j = labels.index(b)
    H[i, j] = J
    H[j, i] = J

evals, evecs = np.linalg.eigh(H)

rows = []
for k, E in enumerate(evals, start=1):
    weights = evecs[:, k-1] ** 2
    pr = 1.0 / np.sum(weights**2)
    dom = labels[int(np.argmax(weights))]
    rows.append({
        "exciton_state": f"X{k}",
        "energy_eV": E,
        "shift_meV_from_mean_site": (E - np.mean(list(sites.values()))) * 1000,
        "participation_ratio": pr,
        "dominant_site": dom,
        **{f"weight_{lab}": weights[i] for i, lab in enumerate(labels)}
    })

pd.DataFrame(H, index=labels, columns=labels).to_csv(OUT / "full_hydrated_hamiltonian_4x4_eV.csv")
pd.DataFrame(H * 1000, index=labels, columns=labels).to_csv(OUT / "full_hydrated_hamiltonian_4x4_meV.csv")
pd.DataFrame(rows).to_csv(OUT / "full_hydrated_exciton_eigenstates_day014.csv", index=False)
pd.DataFrame(evecs, index=labels, columns=[f"X{i}" for i in range(1, len(labels)+1)]).to_csv(
    OUT / "full_hydrated_exciton_eigenvectors_day014.csv"
)

summary = f"""# Full Hydrated Excitonic Hamiltonian — Day014

## Model

Four-site PYR2–PYR5 Hamiltonian with hydrated nearest-neighbor effective couplings:

| Pair | J_eff (meV) | Assignment |
|---|---:|---|
| PYR2-PYR3 | 2.37 | weak/mixed Frenkel-like effective coupling |
| PYR3-PYR4 | 17.00 | bright Frenkel-like effective coupling |
| PYR4-PYR5 | 33.00 | bright Frenkel-like effective coupling |

## Hamiltonian files

- `full_hydrated_hamiltonian_4x4_eV.csv`
- `full_hydrated_hamiltonian_4x4_meV.csv`
- `full_hydrated_exciton_eigenstates_day014.csv`
- `full_hydrated_exciton_eigenvectors_day014.csv`

## Interpretation

The hydrated network is strongly asymmetric. PYR4-PYR5 remains the dominant hydrated coupling, PYR3-PYR4 is intermediate, and PYR2-PYR3 is weak/mixed. This produces a nonuniform excitonic chain rather than a homogeneous nearest-neighbor Frenkel model.

The low-energy CT-like dark states detected in the explicit hydrated dimers were excluded from the Frenkel Hamiltonian construction. Only bright or operationally Frenkel-like states were used to define effective couplings.
"""

(OUT / "FULL_HYDRATED_HAMILTONIAN_DAY014.md").write_text(summary)

print("Wrote:", OUT)
print(pd.DataFrame(H, index=labels, columns=labels))
print(pd.DataFrame(rows).to_string(index=False))
