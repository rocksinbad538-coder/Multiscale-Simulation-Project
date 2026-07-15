from pathlib import Path
import csv

# Geometry from Day012 excitonic geometry audit
geometry = {
    "PYR2-PYR3": {"center_distance_nm": 2.69177, "min_atom_distance_nm": 2.22428, "plane_angle_deg": 89.997},
    "PYR3-PYR4": {"center_distance_nm": 2.69074, "min_atom_distance_nm": 2.22258, "plane_angle_deg": 89.923},
    "PYR4-PYR5": {"center_distance_nm": 2.76521, "min_atom_distance_nm": 2.29794, "plane_angle_deg": 89.460},
    "PYR2-PYR4": {"center_distance_nm": 4.25169, "min_atom_distance_nm": 3.95167, "plane_angle_deg": 1.149},
    "PYR3-PYR5": {"center_distance_nm": 4.34339, "min_atom_distance_nm": 4.05089, "plane_angle_deg": 0.993},
    "PYR2-PYR5": {"center_distance_nm": 4.69020, "min_atom_distance_nm": 4.22280, "plane_angle_deg": 89.518},
}

couplings = {
    "PYR2-PYR3": 6.5,
    "PYR2-PYR4": 2.5,
    "PYR2-PYR5": 8.5,
    "PYR3-PYR4": 10.5,
    "PYR3-PYR5": 1.0,
    "PYR4-PYR5": 11.0,
}

outdir = Path("runs/phase1A/day013_exciton_hamiltonian")
outdir.mkdir(parents=True, exist_ok=True)

rows = []
for pair in sorted(couplings):
    g = geometry[pair]
    rows.append({
        "pair": pair,
        "center_distance_nm": f"{g['center_distance_nm']:.5f}",
        "center_distance_A": f"{g['center_distance_nm']*10:.3f}",
        "min_atom_distance_nm": f"{g['min_atom_distance_nm']:.5f}",
        "min_atom_distance_A": f"{g['min_atom_distance_nm']*10:.3f}",
        "plane_angle_deg": f"{g['plane_angle_deg']:.3f}",
        "J_eff_meV": f"{couplings[pair]:.3f}",
        "geometry_class": (
            "near-parallel_long-range" if g["plane_angle_deg"] < 10
            else "near-orthogonal"
        )
    })

csvout = outdir / "geometry_vs_coupling_day013.csv"
with csvout.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)

mdout = outdir / "GEOMETRY_COUPLING_INTERPRETATION.md"
mdout.write_text("""# Day013 Geometry–Coupling Interpretation

This table combines MD-derived pyrene pair geometry with preliminary TDDFT-derived effective couplings.

Effective couplings were estimated as:

J_eff = (E2 - E1) / 2

from the two lowest dimer excited states.

Important interpretation:

- The pyrenes are not in close pi-stacked contact.
- Minimum atom-atom separations are approximately 22–42 Å.
- Couplings are weak, in the 1–11 meV range.
- The strongest preliminary couplings occur in near-neighbor, near-orthogonal geometries.
- The same-orientation long-range pairs are not necessarily the strongest coupled pairs.

This indicates that the excitonic network is geometry-sensitive and should not be inferred from plane-angle classification alone.

Limitations:

These J_eff values are first-pass splitting-derived couplings. They should be refined with transition-dipole, transition-density, or fragment-based coupling analysis.
""")

print(csvout.read_text())
print(mdout.read_text())
