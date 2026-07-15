from pathlib import Path
import csv
import matplotlib.pyplot as plt

outdir = Path("runs/phase1A/day013_exciton_hamiltonian")
eig_csv = outdir / "exciton_eigenstates_day013.csv"
geo_csv = outdir / "geometry_vs_coupling_day013.csv"

states = []
energies = []
prs = []

with eig_csv.open() as f:
    r = csv.DictReader(f)
    for row in r:
        states.append(row["exciton_state"])
        energies.append(float(row["energy_eV"]))
        prs.append(float(row["participation_ratio"]))

pairs = []
J = []
dist = []
angle = []

with geo_csv.open() as f:
    r = csv.DictReader(f)
    for row in r:
        pairs.append(row["pair"])
        J.append(float(row["J_eff_meV"]))
        dist.append(float(row["center_distance_A"]))
        angle.append(float(row["plane_angle_deg"]))

plt.figure()
plt.bar(states, energies)
plt.ylabel("Exciton energy (eV)")
plt.xlabel("Exciton state")
plt.title("Day013 preliminary exciton energies")
plt.tight_layout()
plt.savefig(outdir / "day013_exciton_energies.png", dpi=300)
plt.close()

plt.figure()
plt.bar(states, prs)
plt.ylabel("Participation ratio")
plt.xlabel("Exciton state")
plt.title("Day013 exciton participation ratios")
plt.tight_layout()
plt.savefig(outdir / "day013_exciton_participation_ratios.png", dpi=300)
plt.close()

plt.figure()
plt.scatter(dist, J)
for x, y, label in zip(dist, J, pairs):
    plt.annotate(label, (x, y), fontsize=8)
plt.xlabel("Center-center distance (Å)")
plt.ylabel("J_eff (meV)")
plt.title("Geometry vs coupling: distance")
plt.tight_layout()
plt.savefig(outdir / "day013_J_vs_distance.png", dpi=300)
plt.close()

plt.figure()
plt.scatter(angle, J)
for x, y, label in zip(angle, J, pairs):
    plt.annotate(label, (x, y), fontsize=8)
plt.xlabel("Plane angle (deg)")
plt.ylabel("J_eff (meV)")
plt.title("Geometry vs coupling: orientation")
plt.tight_layout()
plt.savefig(outdir / "day013_J_vs_plane_angle.png", dpi=300)
plt.close()

print("Generated:")
for p in sorted(outdir.glob("day013_*.png")):
    print(p)
