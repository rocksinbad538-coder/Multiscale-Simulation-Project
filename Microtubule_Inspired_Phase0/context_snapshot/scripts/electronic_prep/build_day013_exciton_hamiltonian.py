from pathlib import Path
import csv

site_energies = {
    "PYR2": 3.779,
    "PYR3": 3.774,
    "PYR4": 3.782,
    "PYR5": 3.767,
}

couplings = {
    ("PYR2","PYR3"): 0.0065,
    ("PYR2","PYR4"): 0.0025,
    ("PYR2","PYR5"): 0.0085,
    ("PYR3","PYR4"): 0.0105,
    ("PYR3","PYR5"): 0.0010,
    ("PYR4","PYR5"): 0.0110,
}

sites = ["PYR2","PYR3","PYR4","PYR5"]

outdir = Path("runs/phase1A/day013_exciton_hamiltonian")
outdir.mkdir(parents=True, exist_ok=True)

matrix = []
for i,a in enumerate(sites):
    row = []
    for j,b in enumerate(sites):
        if i == j:
            row.append(site_energies[a])
        else:
            key = tuple(sorted((a,b)))
            row.append(couplings[key])
    matrix.append(row)

csvout = outdir / "exciton_hamiltonian_4x4_eV.csv"
with csvout.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["site"] + sites)
    for site,row in zip(sites,matrix):
        w.writerow([site] + [f"{x:.6f}" for x in row])

mevout = outdir / "exciton_hamiltonian_4x4_meV_relative.csv"
ref = sum(site_energies.values()) / len(site_energies)

with mevout.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["site"] + sites)
    for i,site in enumerate(sites):
        row = []
        for j,x in enumerate(matrix[i]):
            if i == j:
                row.append((x-ref)*1000)
            else:
                row.append(x*1000)
        w.writerow([site] + [f"{x:.3f}" for x in row])

summary = outdir / "EXCITON_HAMILTONIAN_SUMMARY.md"
summary.write_text(f"""# Day013 Preliminary Excitonic Hamiltonian

Hamiltonian basis:

{", ".join(sites)}

Diagonal entries are isolated pyrene S1 site energies from ORCA TDDFT.

Off-diagonal entries are preliminary effective couplings estimated as:

J_eff = (E2 - E1) / 2

from the two lowest TDDFT dimer excited states.

Average site energy:

{ref:.6f} eV

Files:

- exciton_hamiltonian_4x4_eV.csv
- exciton_hamiltonian_4x4_meV_relative.csv

Important limitation:

These are first-pass effective couplings from excited-state splittings. They are useful for a preliminary Hamiltonian, but should later be refined using transition-density, transition-dipole, or fragment-based coupling analysis.
""")

print(csvout.read_text())
print(mevout.read_text())
print(summary.read_text())
