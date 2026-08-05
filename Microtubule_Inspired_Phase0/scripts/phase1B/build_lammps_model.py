#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import json
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
    / "day041_lammps_model"
)

RUN.mkdir(parents=True, exist_ok=True)


def utc():
    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )


system = json.loads(SYSTEM.read_text())

atom_type_ids = {
    "B":1,
    "N":2,
    "HB":3,
    "HN":4
}

bond_type_ids = {
    "B-N":1,
    "B-HB":2,
    "N-HN":3
}

angle_type_ids = {
    "B-N-B":1,
    "N-B-N":2,
    "B-N-HN":3,
    "N-B-HB":4
}

improper_type_ids = {
    "B_PLANAR":1,
    "N_PLANAR":2
}

dihedral_type_ids = {
    "BASAL":1,
    "B-N-B-HB":2,
    "N-B-N-HN":3
}

model = {

    "timestamp":utc(),

    "counts":system["metadata"],

    "type_tables":{

        "atoms":atom_type_ids,

        "bonds":bond_type_ids,

        "angles":angle_type_ids,

        "impropers":improper_type_ids,

        "dihedrals":dihedral_type_ids

    },

    "system":system

}

outfile = RUN/"PHASE1B_LAMMPS_MODEL.json"

outfile.write_text(
    json.dumps(
        model,
        indent=2
    )
)

print("="*100)
print("DAY041 / PHASE1B-A16")
print("LAMMPS MODEL BUILDER")
print("="*100)
print()

print("Atom types      :",len(atom_type_ids))
print("Bond types      :",len(bond_type_ids))
print("Angle types     :",len(angle_type_ids))
print("Improper types  :",len(improper_type_ids))
print("Dihedral types  :",len(dihedral_type_ids))

print()

print(outfile)

print()

print("PHASE1B_LAMMPS_MODEL_CREATED")
