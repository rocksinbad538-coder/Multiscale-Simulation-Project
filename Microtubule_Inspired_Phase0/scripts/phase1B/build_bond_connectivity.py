#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import csv
import json
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

TOPOLOGY = (
ROOT
/"runs"
/
"phase1B"
/
"day041_internal_topology"
/
"PHASE1B_INTERNAL_TOPOLOGY.json"
)

EDGES = (
ROOT
/"runs"
/
"phase1A"
/
"day035_qm_f06_upper_v7a_r1_coordinate_adoption"
/
"QM_F06_UPPER_V7A_ADOPTED_nominal_edges.csv"
)

RUN = (
ROOT
/"runs"
/
"phase1B"
/
"day041_bond_connectivity"
)

RUN.mkdir(parents=True,exist_ok=True)

def utc():

    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )

topology=json.loads(TOPOLOGY.read_text())

physical_ids=set()

for atom in topology["atoms"]:

    physical_ids.add(atom["atom_id"])

bonds=[]

with open(EDGES) as f:

    reader=csv.DictReader(f)

    for row in reader:

        a=row["first_atom"]

        b=row["second_atom"]

        if a not in physical_ids:

            continue

        if b not in physical_ids:

            continue

        bonds.append({

            "atom1":a,

            "atom2":b,

            "distance_A":float(row["distance_A"]),

            "edge_type":row["edge_type"],

            "provenance":row["provenance"]

        })

print("="*100)
print("DAY041 / PHASE1B-A2")
print("BONDED CONNECTIVITY")
print("="*100)
print()

print("[1] PHYSICAL BONDS")

print(len(bonds))

report={

"timestamp":utc(),

"bond_count":len(bonds),

"bonds":bonds,

"status":"PHYSICAL_CONNECTIVITY_CREATED"

}

outfile=RUN/"PHASE1B_BOND_CONNECTIVITY.json"

outfile.write_text(

json.dumps(

report,

indent=2

)

)

print()

print("[2] OUTPUT")

print(outfile)

print()

print("[3] DECISION")

print("PHASE1B_BOND_CONNECTIVITY_CREATED")
