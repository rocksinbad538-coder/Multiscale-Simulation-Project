#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import json
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
    / "day041_lammps_index_validation"
)

RUN.mkdir(parents=True,exist_ok=True)


def utc():
    return (
        datetime.datetime.now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )


model=json.loads(MODEL.read_text())

indices=sorted(

a["lammps_index"]

for a in model["system"]["atoms"]

)

expected=list(range(1,len(indices)+1))

ok=(indices==expected)

report={

"timestamp":utc(),

"atom_count":len(indices),

"indices_contiguous":ok,

"first_index":indices[0],

"last_index":indices[-1],

"decision":"PASS" if ok else "FAIL"

}

outfile=RUN/"LAMMPS_INDEX_VALIDATION.json"

outfile.write_text(

json.dumps(report,indent=2)

)

print("="*100)
print("DAY041 / PHASE1B-A19")
print("LAMMPS INDEX VALIDATION")
print("="*100)
print()

for k,v in report.items():

    if k in ("timestamp","decision"):
        continue

    print(f"{k:20s} {v}")

print()

print(outfile)

print()

print(report["decision"])
