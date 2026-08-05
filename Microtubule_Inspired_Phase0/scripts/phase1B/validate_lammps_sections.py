#!/usr/bin/env python3

from __future__ import annotations

import json
import pathlib
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

SECTIONS=(
ROOT
/"runs"
/
"phase1B"
/
"day041_lammps_sections"
/
"PHASE1B_LAMMPS_SECTIONS.json"
)

RUN=(
ROOT
/"runs"
/
"phase1B"
/
"day041_lammps_section_validation"
)

RUN.mkdir(parents=True,exist_ok=True)

def utc():

    return (
        datetime.datetime.now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )

data=json.loads(SECTIONS.read_text())

s=data["sections"]

tests={

"Header":True,

"Masses":len(s["Masses"])>0,

"BondCoeffs":len(s["BondCoeffs"])>0,

"AngleCoeffs":len(s["AngleCoeffs"])>0,

"ImproperCoeffs":len(s["ImproperCoeffs"])>0,

"DihedralCoeffs":len(s["DihedralCoeffs"])>0,

"Atoms":len(s["Atoms"])==s["Header"]["atoms"],

"Bonds":len(s["Bonds"])==s["Header"]["bonds"],

"Angles":len(s["Angles"])==s["Header"]["angles"],

"Impropers":len(s["Impropers"])==s["Header"]["impropers"],

"Dihedrals":len(s["Dihedrals"])==s["Header"]["dihedrals"]

}

decision="PASS"

if not all(tests.values()):

    decision="FAIL"

outfile=RUN/"LAMMPS_SECTION_VALIDATION.json"

outfile.write_text(

json.dumps(

{

"timestamp":utc(),

"tests":tests,

"decision":decision

},

indent=2

)

)

print("="*100)
print("DAY041 / PHASE1B-A21")
print("LAMMPS SECTION VALIDATION")
print("="*100)
print()

for k,v in tests.items():

    print(f"{k:20s}{v}")

print()

print(outfile)

print()

print(decision)
