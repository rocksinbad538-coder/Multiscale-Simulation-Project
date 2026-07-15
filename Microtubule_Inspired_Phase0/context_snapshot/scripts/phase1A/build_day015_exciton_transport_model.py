from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT = Path("runs/phase1A/day015_exciton_model")
OUT.mkdir(parents=True, exist_ok=True)

labels = ["PYR2", "PYR3", "PYR4", "PYR5"]

# Day014 validated hydrated Hamiltonian, eV
H = np.array([
    [3.75500, 0.00237, 0.00000, 0.00000],
    [0.00237, 3.77750, 0.01700, 0.00000],
    [0.00000, 0.01700, 3.77750, 0.03300],
    [0.00000, 0.00000, 0.03300, 3.77750],
])

evals, evecs = np.linalg.eigh(H)
weights = evecs**2
pr = 1.0 / np.sum(weights**2, axis=0)

pd.DataFrame(H, index=labels, columns=labels).to_csv(OUT / "hydrated_hamiltonian.csv")
pd.DataFrame({
    "pair": ["PYR2-PYR3", "PYR3-PYR4", "PYR4-PYR5"],
    "J_meV": [2.37, 17.00, 33.00],
    "assignment": [
        "weak/mixed Frenkel-like effective coupling",
        "bright Frenkel-like effective coupling",
        "bright Frenkel-like effective coupling",
    ],
}).to_csv(OUT / "hydrated_couplings.csv", index=False)

rows = []
for k, E in enumerate(evals):
    w = weights[:, k]
    rows.append({
        "exciton_state": f"X{k+1}",
        "energy_eV": E,
        "participation_ratio": pr[k],
        "dominant_site": labels[int(np.argmax(w))],
        **{f"weight_{lab}": w[i] for i, lab in enumerate(labels)}
    })

pd.DataFrame(rows).to_csv(OUT / "hydrated_exciton_states.csv", index=False)

plt.figure(figsize=(6,4))
plt.plot(range(1, 5), np.diag(H), "o-", label="site energies")
plt.plot(range(1, 5), evals, "s-", label="exciton energies")
plt.xlabel("state/site index")
plt.ylabel("Energy (eV)")
plt.legend()
plt.tight_layout()
plt.savefig(OUT / "hydrated_exciton_spectrum.png", dpi=300)

summary = """# Day015 Hydrated Exciton Model

This package consolidates the Day014 explicit-hydration Hamiltonian into a reusable exciton-model input.

## Main conclusion

The hydrated four-site PYR2-PYR5 network is strongly asymmetric:

- PYR2-PYR3: 2.37 meV
- PYR3-PYR4: 17.00 meV
- PYR4-PYR5: 33.00 meV

This implies a weakly connected PYR2 site and a more strongly coupled PYR3-PYR5 subnetwork.

## Files

- hydrated_hamiltonian.csv
- hydrated_couplings.csv
- hydrated_exciton_states.csv
- hydrated_exciton_spectrum.png
"""
(OUT / "DAY015_EXCITON_MODEL_SUMMARY.md").write_text(summary)

print("Wrote:", OUT)
print(pd.DataFrame(rows).to_string(index=False))
