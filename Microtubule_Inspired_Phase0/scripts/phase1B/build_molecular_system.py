#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import json
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

INPUTS = {

"atoms":
ROOT/"runs"/"phase1B"/"day041_internal_topology"/"PHASE1B_INTERNAL_TOPOLOGY.json",

"bonds":
ROOT/"runs"/"phase1B"/"day041_bond_connectivity"/"PHASE1B_BOND_CONNECTIVITY.json",

"angles":
ROOT/"runs"/"phase1B"/"day041_angle_parameter_mapping"/"PHASE1B_CANONICAL_ANGLES.json",

"impropers":
ROOT/"runs"/"phase1B"/"day041_impropers"/"PHASE1B_IMPROPERS.json",

"dihedrals":
ROOT/"runs"/"phase1B"/"day041_dihedrals"/"PHASE1B_DIHEDRALS.json"

}

RUN = (
ROOT
/"runs"
/
"phase1B"
/
"day041_molecular_system"
)

RUN.mkdir(parents=True,exist_ok=True)

def utc():

    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )

system = {

"timestamp":utc(),

"phase":"Phase1B",

"status":"TOPOLOGY_ASSEMBLED",

"metadata":{},


"atoms":None,

"bonds":None,

"angles":None,

"impropers":None,

"dihedrals":None

}

print("="*100)
print("DAY041 / PHASE1B-A10")
print("MOLECULAR SYSTEM ASSEMBLY")
print("="*100)
print()

for name,path in INPUTS.items():

    data=json.loads(path.read_text())

    if name=="atoms":

        system["atoms"]=data["atoms"]

        system["metadata"]["atom_count"]=len(data["atoms"])

    elif name=="bonds":

        system["bonds"]=data["bonds"]

        system["metadata"]["bond_count"]=len(data["bonds"])

    elif name=="angles":

        system["angles"]=data["angles"]

        system["metadata"]["angle_count"]=len(data["angles"])

    elif name=="impropers":

        system["impropers"]=data["impropers"]

        system["metadata"]["improper_count"]=len(data["impropers"])

    elif name=="dihedrals":

        system["dihedrals"]=data["dihedrals"]

        system["metadata"]["dihedral_count"]=len(data["dihedrals"])

    print(f"{name:12s} PASS")

outfile=RUN/"PHASE1B_MOLECULAR_SYSTEM.json"

outfile.write_text(

json.dumps(

system,

indent=2

)

)

print()

print("[SUMMARY]")

for k,v in system["metadata"].items():

    print(f"{k:20s} {v}")

print()

print("[OUTPUT]")

print(outfile)

print()

print("[DECISION]")

print("PHASE1B_MOLECULAR_SYSTEM_CREATED")
