#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import json
import itertools
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

GRAPH = (
ROOT
/"runs"
/
"phase1B"
/
"day041_molecular_graph"
/
"PHASE1B_MOLECULAR_GRAPH.json"
)

RUN = (
ROOT
/"runs"
/
"phase1B"
/
"day041_angles"
)

RUN.mkdir(parents=True,exist_ok=True)

def utc():

    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )

graph=json.loads(GRAPH.read_text())["graph"]

angles=[]

seen=set()

for center,data in graph.items():

    neigh=data["neighbors"]

    if len(neigh)<2:

        continue

    for a,c in itertools.combinations(sorted(neigh),2):

        key=(a,center,c)

        if key in seen:

            continue

        seen.add(key)

        angles.append({

            "atom1":a,

            "atom2":center,

            "atom3":c

        })

print("="*100)
print("DAY041 / PHASE1B-A4")
print("ANGLE ENUMERATION")
print("="*100)
print()

print("[1] ANGLES")

print(len(angles))

outfile=RUN/"PHASE1B_ANGLES.json"

outfile.write_text(

json.dumps({

"timestamp":utc(),

"angle_count":len(angles),

"angles":angles

},indent=2)

)

print()

print("[2] OUTPUT")

print(outfile)

print()

print("[3] DECISION")

print("PHASE1B_ANGLES_CREATED")
