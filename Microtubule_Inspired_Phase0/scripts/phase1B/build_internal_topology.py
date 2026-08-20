#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import csv
import json
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

INPUT = (
ROOT
/"runs"
/
"phase1A"
/
"day040_phase1A_parameter_mapping"
/
"PHASE1B_PARAMETER_MAPPING.csv"
)

ADOPTED = (
ROOT
/"runs"
/
"phase1A"
/
"day039_phase1A_F_charge_model_closure"
/
"QM_F06_UPPER_V7A_R1_PHASE1A_F_ADOPTED_WORKING_CHARGES.csv"
)

RUN = (
ROOT
/"runs"
/
"phase1B"
/
"day041_internal_topology"
)

RUN.mkdir(
parents=True,
exist_ok=True
)

def utc():

    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )

adopted_charge={}

with open(ADOPTED) as f:

    reader=csv.DictReader(f)

    for row in reader:

        adopted_charge[
            row["atom_id"]
        ]=float(
            row["adopted_working_charge_e"]
        )

atoms=[]

with open(INPUT) as f:

    reader=csv.DictReader(f)

    for row in reader:

        if row["proposed_forcefield_type"]=="H":

            continue

        atom=dict(row)

        atom["forcefield_type"]=row["proposed_forcefield_type"]

        atom["RESP_charge_e"]=adopted_charge[row["atom_id"]]

        atom["role"]=row["atom_role"]

        atoms.append(atom)

print("="*100)
print("DAY041 / PHASE1B-A1")
print("INTERNAL TOPOLOGY CONSTRUCTION")
print("="*100)
print()

print("[1] PHYSICAL ATOMS")

print(len(atoms))

types={}

for a in atoms:

    t=a["forcefield_type"]

    types[t]=types.get(t,0)+1

print()

print("[2] FORCE FIELD TYPES")

for k in sorted(types):

    print(f"{k:>3s} {types[k]}")

topology={

"timestamp":utc(),

"physical_atom_count":len(atoms),

"atom_types":types,

"atoms":atoms,

"status":"INTERNAL_TOPOLOGY_INITIALIZED"

}

outfile=RUN/"PHASE1B_INTERNAL_TOPOLOGY.json"

outfile.write_text(

json.dumps(

topology,

indent=2

)

)

print()

print("[3] OUTPUT")

print(outfile)

print()

print("[4] DECISION")

print("PHASE1B_INTERNAL_TOPOLOGY_CREATED")
