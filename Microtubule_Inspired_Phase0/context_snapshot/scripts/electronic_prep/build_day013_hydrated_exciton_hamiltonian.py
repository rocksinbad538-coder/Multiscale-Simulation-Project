from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt

sites = ["PYR2", "PYR3", "PYR4", "PYR5"]

# Hydrated site energies from pyrene + local water 0.50 nm TDDFT
site_energies_hyd = {
    "PYR2": 3.754,
    "PYR3": 3.765,
    "PYR4": 3.779,
    "PYR5": 3.746,
}

# Vacuum site energies from isolated pyrene TDDFT
site_energies_vac = {
    "PYR2": 3.779,
    "PYR3": 3.774,
    "PYR4": 3.782,
    "PYR5": 3.767,
}

# First-pass splitting-derived couplings, eV
couplings = {
    ("PYR2","PYR3"): 0.0065,
    ("PYR2","PYR4"): 0.0025,
    ("PYR2","PYR5"): 0.0085,
    ("PYR3","PYR4"): 0.0105,
    ("PYR3","PYR5"): 0.0010,
    ("PYR4","PYR5"): 0.0110,
}

outdir = Path("runs/phase1A/day013_hydrated_exciton_hamiltonian")
outdir.mkdir(parents=True, exist_ok=True)

def build_H(site_energies):
    H = np.zeros((4,4), dtype=float)
    for i,a in enumerate(sites):
        for j,b in enumerate(sites):
            if i == j:
                H[i,j] = site_energies[a]
            else:
                H[i,j] = couplings[tuple(sorted((a,b)))]
    return H

def diagonalize(H):
    evals, evecs = np.linalg.eigh(H)
    rows = []
    for k,E in enumerate(evals):
        coeffs = evecs[:,k]
        weights = coeffs**2
        pr = 1.0 / np.sum(weights**2)
        rows.append({
            "exciton_state": f"X{k+1}",
            "energy_eV": E,
            "shift_meV_from_mean_site": (E - np.mean(np.diag(H))) * 1000.0,
            "participation_ratio": pr,
            "dominant_site": sites[int(np.argmax(weights))],
            **{f"weight_{s}": weights[i] for i,s in enumerate(sites)}
        })
    return evals, evecs, rows

H_hyd = build_H(site_energies_hyd)
H_vac = build_H(site_energies_vac)

evals_hyd, evecs_hyd, rows_hyd = diagonalize(H_hyd)
evals_vac, evecs_vac, rows_vac = diagonalize(H_vac)

# Write hydrated Hamiltonian
with (outdir / "hydrated_exciton_hamiltonian_4x4_eV.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["site"] + sites)
    for i,s in enumerate(sites):
        w.writerow([s] + [f"{H_hyd[i,j]:.6f}" for j in range(4)])

# Write hydrated eigenstates
with (outdir / "hydrated_exciton_eigenstates.csv").open("w", newline="") as f:
    fieldnames = list(rows_hyd[0].keys())
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for r in rows_hyd:
        rr = r.copy()
        for k,v in rr.items():
            if isinstance(v, float):
                rr[k] = f"{v:.6f}"
        w.writerow(rr)

# Write hydrated eigenvectors
with (outdir / "hydrated_exciton_eigenvectors.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["site"] + [f"X{k+1}" for k in range(4)])
    for i,s in enumerate(sites):
        w.writerow([s] + [f"{evecs_hyd[i,k]:.8f}" for k in range(4)])

# Vacuum vs hydrated comparison
with (outdir / "vacuum_vs_hydrated_exciton_comparison.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow([
        "exciton_state",
        "vacuum_energy_eV",
        "hydrated_energy_eV",
        "energy_shift_meV",
        "vacuum_PR",
        "hydrated_PR",
        "PR_change"
    ])
    for rv, rh in zip(rows_vac, rows_hyd):
        w.writerow([
            rv["exciton_state"],
            f"{rv['energy_eV']:.6f}",
            f"{rh['energy_eV']:.6f}",
            f"{(rh['energy_eV'] - rv['energy_eV']) * 1000.0:.3f}",
            f"{rv['participation_ratio']:.3f}",
            f"{rh['participation_ratio']:.3f}",
            f"{rh['participation_ratio'] - rv['participation_ratio']:.3f}",
        ])

# Site-energy shifts
with (outdir / "hydrated_site_energy_shifts.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["site", "vacuum_S1_eV", "hydrated_S1_eV", "shift_meV"])
    for s in sites:
        w.writerow([
            s,
            f"{site_energies_vac[s]:.6f}",
            f"{site_energies_hyd[s]:.6f}",
            f"{(site_energies_hyd[s] - site_energies_vac[s])*1000.0:.3f}"
        ])

# Plot PR comparison
plt.figure()
x = np.arange(1,5)
plt.plot(x, [r["participation_ratio"] for r in rows_vac], marker="o", label="vacuum")
plt.plot(x, [r["participation_ratio"] for r in rows_hyd], marker="o", label="hydrated")
plt.xticks(x, [f"X{i}" for i in x])
plt.xlabel("Exciton state")
plt.ylabel("Participation ratio")
plt.title("Vacuum vs hydrated exciton delocalization")
plt.legend()
plt.tight_layout()
plt.savefig(outdir / "vacuum_vs_hydrated_participation_ratio.png", dpi=300)
plt.close()

# Plot exciton energies
plt.figure()
plt.plot(x, [r["energy_eV"] for r in rows_vac], marker="o", label="vacuum")
plt.plot(x, [r["energy_eV"] for r in rows_hyd], marker="o", label="hydrated")
plt.xticks(x, [f"X{i}" for i in x])
plt.xlabel("Exciton state")
plt.ylabel("Energy (eV)")
plt.title("Vacuum vs hydrated exciton energies")
plt.legend()
plt.tight_layout()
plt.savefig(outdir / "vacuum_vs_hydrated_exciton_energies.png", dpi=300)
plt.close()

summary = outdir / "VACUUM_vs_HYDRATED_COMPARISON.md"
summary.write_text("""# Day013 Vacuum vs Hydrated Exciton Hamiltonian Comparison

## Purpose

Assess whether local explicit hydration preserves or suppresses excitonic delocalization in the Phase 1A four-pyrene network.

## Model

Two Hamiltonians were compared:

1. Vacuum Hamiltonian:
   - isolated pyrene TDDFT site energies;
   - first-pass splitting-derived dimer couplings.

2. Hydrated Hamiltonian:
   - pyrene + local water 0.50 nm TDDFT site energies;
   - same first-pass splitting-derived dimer couplings.

## Hydrated site energies

| Site | Vacuum S1 (eV) | Hydrated S1 (eV) | Shift (meV) |
|---|---:|---:|---:|
""" + "\n".join(
    f"| {s} | {site_energies_vac[s]:.3f} | {site_energies_hyd[s]:.3f} | {(site_energies_hyd[s]-site_energies_vac[s])*1000.0:.1f} |"
    for s in sites
) + f"""

## Site-energy disorder

Vacuum site-energy range:

{(max(site_energies_vac.values()) - min(site_energies_vac.values()))*1000.0:.1f} meV

Hydrated site-energy range:

{(max(site_energies_hyd.values()) - min(site_energies_hyd.values()))*1000.0:.1f} meV

## Hydrated exciton eigenstates

| State | Energy (eV) | Shift from hydrated mean site energy (meV) | PR | Dominant site |
|---|---:|---:|---:|---|
""" + "\n".join(
    f"| {r['exciton_state']} | {r['energy_eV']:.6f} | {r['shift_meV_from_mean_site']:.3f} | {r['participation_ratio']:.3f} | {r['dominant_site']} |"
    for r in rows_hyd
) + """

## Vacuum vs hydrated participation

| State | Vacuum PR | Hydrated PR | Change |
|---|---:|---:|---:|
""" + "\n".join(
    f"| {rv['exciton_state']} | {rv['participation_ratio']:.3f} | {rh['participation_ratio']:.3f} | {rh['participation_ratio']-rv['participation_ratio']:.3f} |"
    for rv, rh in zip(rows_vac, rows_hyd)
) + """

## Interpretation

Explicit local hydration red-shifts all pyrene S1 energies and increases energetic disorder relative to the isolated-chromophore Hamiltonian.

The key question is whether this disorder suppresses excitonic delocalization. The participation-ratio comparison quantifies this directly.

These results remain preliminary because the off-diagonal couplings are still splitting-derived and were not recalculated in the hydrated environment. The next refinement should compare splitting-derived couplings against transition-dipole or transition-density-based couplings.
""")

print((outdir / "hydrated_exciton_hamiltonian_4x4_eV.csv").read_text())
print((outdir / "hydrated_exciton_eigenstates.csv").read_text())
print((outdir / "vacuum_vs_hydrated_exciton_comparison.csv").read_text())
print(summary.read_text())
print("Generated figures:")
for p in sorted(outdir.glob("*.png")):
    print(p)
