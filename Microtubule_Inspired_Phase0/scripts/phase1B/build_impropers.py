#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import json
import itertools
import collections
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
"day041_impropers"
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

impropers=[]

hist=collections.Counter()

for center,data in graph.items():

    if len(data["neighbors"])!=3:

        continue

    neigh=sorted(data["neighbors"])

    atomtypes=[graph[x]["type"] for x in neigh]

    center_type=data["type"]

    improper_type=f"{center_type}"

    hist[improper_type]+=1

    impropers.append({

        "center":center,

        "neighbors":neigh,

        "improper_type":improper_type

    })

print("="*100)
print("DAY041 / PHASE1B-A7")
print("IMPROPER ENUMERATION")
print("="*100)
print()

print("[1] IMPROPERS")

print(len(impropers))

print()

print("[2] CENTRAL TYPES")

for k in sorted(hist):

    print(f"{k:3s} {hist[k]}")

outfile=RUN/"PHASE1B_IMPROPERS.json"

outfile.write_text(

json.dumps({

"timestamp":utc(),

"impropers":impropers,

"histogram":dict(hist)

},indent=2)

)

print()

print("[3] OUTPUT")

print(outfile)

print()

print("[4] DECISION")

print("PHASE1B_IMPROPERS_CREATED")
