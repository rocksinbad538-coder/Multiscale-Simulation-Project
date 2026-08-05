#!/usr/bin/env python3

from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]

PIPELINE = [

"build_internal_topology.py",

"build_bond_connectivity.py",

"build_molecular_graph.py",

"build_angles.py",

"type_angles.py",

"map_angle_parameters.py",

"build_impropers.py",

"build_dihedrals.py",

"validate_dihedrals.py",

"build_molecular_system.py",

"assign_bond_parameters.py",

"assign_angle_parameters.py",

"assign_planarity_parameters.py",

"forcefield_coverage_audit.py",

"validate_molecular_system.py",

"build_lammps_model.py",

"validate_lammps_model.py",

"build_lammps_indices.py",

"validate_indexed_model.py",

"build_lammps_sections.py",

"validate_lammps_sections.py",

"write_lammps_data.py"

]

print("="*100)
print("PHASE1B REBUILD")
print("="*100)

for script in PIPELINE:

    print()
    print("-"*80)
    print(script)

    result = subprocess.run(
        [
            sys.executable,
            str(ROOT/"scripts"/"phase1B"/script)
        ]
    )

    if result.returncode != 0:

        print()
        print("="*100)
        print("PIPELINE STOPPED")
        print(script)
        print("="*100)

        sys.exit(result.returncode)

print()
print("="*100)
print("PHASE1B REBUILD COMPLETE")
print("="*100)
