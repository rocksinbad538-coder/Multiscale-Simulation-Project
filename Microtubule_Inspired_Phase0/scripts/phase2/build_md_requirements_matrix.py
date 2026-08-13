#!/usr/bin/env python3

from pathlib import Path
import csv
import json

ROOT = Path(__file__).resolve().parents[2]

requirements = [

    ["Coordinates","YES","Current MD analysis","YES","dump custom"],

    ["Simulation box","YES","Periodic geometry","YES","dump"],

    ["Potential energy","YES","Thermodynamics","YES","thermo_style"],

    ["Temperature","YES","Thermodynamics","YES","thermo_style"],

    ["Pressure","YES","Thermodynamics","YES","thermo_style"],

    ["Radius of gyration","Derived","Structure","YES","post-processing"],

    ["RMSD","Derived","Structure","YES","post-processing"],

    ["RMSF","Derived","Dynamics","YES","post-processing"],

    ["Shape descriptors","Derived","Morphology","YES","post-processing"],

    ["Velocities","Candidate","VACF / diffusion","NO","dump custom vx vy vz"],

    ["Forces","Candidate","Mechanical analysis","NO","dump custom fx fy fz"],

    ["Per-atom energy","Candidate","Energy localization","NO","compute pe/atom"],

    ["Stress tensor","Candidate","Elastic response","NO","compute stress/atom"],

    ["Dipole moment","Pending model review","Dielectric response","NO","compute dipole"],

    ["Per-atom charge","Pending model review","Electrostatics","NO","dump q"]

]

OUT = ROOT / "runs" / "phase2"

csvfile = OUT / "MD_REQUIREMENTS_MATRIX.csv"

with open(csvfile, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(
        [
            "Observable",
            "Priority",
            "Scientific_use",
            "Currently_available",
            "Implementation"
        ]
    )
    writer.writerows(requirements)

jsonfile = OUT / "MD_REQUIREMENTS_MATRIX.json"

jsonfile.write_text(
    json.dumps(requirements, indent=2)
)

print("=" * 90)
print("DAY046 / PHASE2-A32")
print("MD REQUIREMENTS MATRIX")
print("=" * 90)
print(csvfile)
print(jsonfile)
