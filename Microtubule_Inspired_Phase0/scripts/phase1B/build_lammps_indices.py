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
    / "day041_lammps_model"
    / "PHASE1B_LAMMPS_MODEL.json"
)

RUN = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_lammps_indexed_model"
)

RUN.mkdir(parents=True, exist_ok=True)


def utc():
    return (
        datetime.datetime.now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )


model = json.loads(MODEL.read_text())

atoms = model["system"]["atoms"]

index = {}

for i,atom in enumerate(atoms,start=1):

    atom["lammps_index"]=i

    index[atom["atom_id"]]=i


for collection,keylist in [

("bonds",("atom1","atom2")),

("angles",("atom1","atom2","atom3")),

("dihedrals",("atom1","atom2","atom3","atom4"))

]:

    for obj in model["system"][collection]:

        obj["indices"]=[

            index[obj[k]]

            for k in keylist

        ]


for imp in model["system"]["impropers"]:

    imp["indices"]=[

        index[imp["center"]],

        *[index[x] for x in imp["neighbors"]]

    ]


model["index_table"]=index

model["timestamp_indexing"]=utc()

outfile=RUN/"PHASE1B_LAMMPS_INDEXED_MODEL.json"

outfile.write_text(

json.dumps(

model,

indent=2

)

)

print("="*100)
print("DAY041 / PHASE1B-A18")
print("LAMMPS INDEX ASSIGNMENT")
print("="*100)
print()

print("Atoms indexed :",len(index))

print()

print(outfile)

print()

print("PHASE1B_INDEXING_COMPLETE")
