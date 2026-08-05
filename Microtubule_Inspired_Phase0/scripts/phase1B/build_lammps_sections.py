#!/usr/bin/env python3

from __future__ import annotations

import json
import pathlib
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

MODEL = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_lammps_indexed_model"
    / "PHASE1B_LAMMPS_INDEXED_MODEL.json"
)

RUN = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_lammps_sections"
)

RUN.mkdir(parents=True, exist_ok=True)


def utc():
    return (
        datetime.datetime.now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )


model=json.loads(MODEL.read_text())

system=model["system"]

sections={

"Header":{

"atoms":len(system["atoms"]),

"bonds":len(system["bonds"]),

"angles":len(system["angles"]),

"impropers":len(system["impropers"]),

"dihedrals":len(system["dihedrals"])

},

"Masses":model["type_tables"]["atoms"],

"BondCoeffs":model["type_tables"]["bonds"],

"AngleCoeffs":model["type_tables"]["angles"],

"ImproperCoeffs":model["type_tables"]["impropers"],

"DihedralCoeffs":model["type_tables"]["dihedrals"],

"Atoms":system["atoms"],

"Bonds":system["bonds"],

"Angles":system["angles"],

"Impropers":system["impropers"],

"Dihedrals":system["dihedrals"]

}

outfile=RUN/"PHASE1B_LAMMPS_SECTIONS.json"

outfile.write_text(

json.dumps(

{

"timestamp":utc(),

"sections":sections

},

indent=2

)

)

print("="*100)
print("DAY041 / PHASE1B-A20")
print("LAMMPS SECTION BUILDER")
print("="*100)
print()

for k,v in sections["Header"].items():

    print(f"{k:12s}{v}")

print()

print(outfile)

print()

print("PHASE1B_LAMMPS_SECTIONS_CREATED")
