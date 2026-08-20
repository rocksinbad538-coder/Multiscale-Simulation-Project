#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import json
import collections

HARMONIC_K_KJMOL_TO_LAMMPS_REAL = 1.0 / (2.0 * 4.184)
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

SYSTEM = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_parameterized_system"
    / "PHASE1B_MOLECULAR_SYSTEM.json"
)

RUN = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_parameterized_system"
)


def utc():
    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00", "Z")
    )


system = json.loads(SYSTEM.read_text())

atom_type = {
    a["atom_id"]: a["forcefield_type"]
    for a in system["atoms"]
}

improper_hist = collections.Counter()

for imp in system["impropers"]:

    center = imp["center"]

    ctype = atom_type[center]

    if ctype == "B":

        imp["parameter_type"] = "B_PLANAR"
        imp["theta0_deg"] = 0.0
        imp["k_source_kJmol"] = 258.0
        imp["k_lammps_real_kcalmol"] = (
            258.0 * HARMONIC_K_KJMOL_TO_LAMMPS_REAL
        )

    elif ctype == "N":

        imp["parameter_type"] = "N_PLANAR"
        imp["theta0_deg"] = 0.0
        imp["k_source_kJmol"] = 657.5
        imp["k_lammps_real_kcalmol"] = (
            657.5 * HARMONIC_K_KJMOL_TO_LAMMPS_REAL
        )

    improper_hist[imp["parameter_type"]] += 1


dihedral_hist = collections.Counter()

for dih in system["dihedrals"]:

    types = [
        atom_type[dih["atom1"]],
        atom_type[dih["atom2"]],
        atom_type[dih["atom3"]],
        atom_type[dih["atom4"]],
    ]

    label = "-".join(types)

    if label in ("N-B-N-HN", "HN-N-B-N"):

        ff = "N-B-N-HN"

    elif label in ("B-N-B-HB", "HB-B-N-B"):

        ff = "B-N-B-HB"

    else:

        ff = "BASAL"

    dih["parameter_type"] = ff
    dih["opls"] = [0.0, 0.0, 0.0, 0.0]

    dihedral_hist[ff] += 1


system["metadata"]["improper_parameter_types"] = dict(improper_hist)
system["metadata"]["dihedral_parameter_types"] = dict(dihedral_hist)

system["status"] = "FULL_FORCEFIELD_ASSIGNED"

system["timestamp_parameterization"] = utc()

outfile = RUN / "PHASE1B_MOLECULAR_SYSTEM.json"

outfile.write_text(
    json.dumps(
        system,
        indent=2
    )
)

print("="*100)
print("DAY041 / PHASE1B-A13")
print("PLANARITY PARAMETER ASSIGNMENT")
print("="*100)
print()

print("[1] IMPROPERS")

for k in sorted(improper_hist):
    print(f"{k:12s} {improper_hist[k]}")

print()

print("[2] DIHEDRALS")

for k in sorted(dihedral_hist):
    print(f"{k:12s} {dihedral_hist[k]}")

print()

print("[3] OUTPUT")
print(outfile)

print()

print("[4] DECISION")
print("PHASE1B_FULL_FORCEFIELD_ASSIGNED")
