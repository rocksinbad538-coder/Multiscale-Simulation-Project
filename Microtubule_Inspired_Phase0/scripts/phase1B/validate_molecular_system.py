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
    / "day041_system_validation"
)

RUN.mkdir(parents=True,exist_ok=True)


def utc():
    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )


system=json.loads(SYSTEM.read_text())

ids={a["atom_id"] for a in system["atoms"]}

errors=[]


def check(sequence,keys):

    for obj in sequence:

        for k in keys:

            if obj[k] not in ids:

                errors.append(obj)


check(system["bonds"],["atom1","atom2"])
check(system["angles"],["atom1","atom2","atom3"])
check(system["dihedrals"],["atom1","atom2","atom3","atom4"])

for imp in system["impropers"]:

    if imp["center"] not in ids:

        errors.append(imp)

    for n in imp["neighbors"]:

        if n not in ids:

            errors.append(imp)

report={

"timestamp":utc(),

"atom_count":len(ids),

"errors":len(errors),

"decision":

"PASS"

if len(errors)==0

else

"FAIL"

}

outfile=RUN/"MOLECULAR_SYSTEM_VALIDATION.json"

outfile.write_text(json.dumps(report,indent=2))

print("="*100)
print("DAY041 / PHASE1B-A15")
print("MOLECULAR SYSTEM VALIDATION")
print("="*100)
print()

print("Atoms :",len(ids))
print("Errors:",len(errors))

print()

print(outfile)

print()

print(report["decision"])
