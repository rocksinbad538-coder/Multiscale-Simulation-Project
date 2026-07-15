from pathlib import Path
import csv
import numpy as np

sites = ["PYR2", "PYR3", "PYR4", "PYR5"]

H = np.array([
    [3.7790, 0.0065, 0.0025, 0.0085],
    [0.0065, 3.7740, 0.0105, 0.0010],
    [0.0025, 0.0105, 3.7820, 0.0110],
    [0.0085, 0.0010, 0.0110, 3.7670],
], dtype=float)

outdir = Path("runs/phase1A/day013_exciton_hamiltonian")
outdir.mkdir(parents=True, exist_ok=True)

evals, evecs = np.linalg.eigh(H)

# Ensure eigenvectors are columns; participation ratio = 1/sum(|c_i|^4)
rows = []
for k, E in enumerate(evals):
    coeffs = evecs[:, k]
    weights = coeffs**2
    pr = 1.0 / np.sum(weights**2)
    dominant = sites[int(np.argmax(weights))]
    rows.append({
        "exciton_state": f"X{k+1}",
        "energy_eV": f"{E:.6f}",
        "energy_shift_meV_from_mean_site": f"{(E - np.mean(np.diag(H))) * 1000:.3f}",
        "participation_ratio": f"{pr:.3f}",
        "dominant_site": dominant,
        **{f"weight_{s}": f"{w:.6f}" for s, w in zip(sites, weights)}
    })

csvout = outdir / "exciton_eigenstates_day013.csv"
with csvout.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)

vecout = outdir / "exciton_eigenvectors_day013.csv"
with vecout.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["site"] + [f"X{k+1}" for k in range(len(evals))])
    for i, site in enumerate(sites):
        w.writerow([site] + [f"{evecs[i,k]:.8f}" for k in range(len(evals))])

summary = outdir / "EXCITON_DIAGONALIZATION_SUMMARY.md"
summary.write_text(f"""# Day013 Exciton Hamiltonian Diagonalization

Basis:

{", ".join(sites)}

Input Hamiltonian:

- `exciton_hamiltonian_4x4_eV.csv`

Outputs:

- `exciton_eigenstates_day013.csv`
- `exciton_eigenvectors_day013.csv`

## Exciton energies

| State | Energy (eV) | Shift from mean site energy (meV) | Participation ratio | Dominant site |
|---|---:|---:|---:|---|
""" + "\n".join(
    f"| {r['exciton_state']} | {r['energy_eV']} | {r['energy_shift_meV_from_mean_site']} | {r['participation_ratio']} | {r['dominant_site']} |"
    for r in rows
) + """

## Interpretation

The participation ratio estimates how many pyrene sites contribute appreciably to each excitonic eigenstate.

- PR ≈ 1: localized mostly on one chromophore.
- PR ≈ 2: delocalized over approximately two chromophores.
- PR ≈ 4: delocalized over the full four-site network.

These results are based on a preliminary splitting-derived Hamiltonian and should be refined after transition-dipole or transition-density coupling analysis.
""")

print(csvout.read_text())
print(vecout.read_text())
print(summary.read_text())
