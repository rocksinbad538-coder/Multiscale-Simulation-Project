#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import json
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

GRAPH = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_molecular_graph"
    / "PHASE1B_MOLECULAR_GRAPH.json"
)

RUN = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_dihedrals"
)

RUN.mkdir(parents=True, exist_ok=True)


def utc():
    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )


graph = json.loads(GRAPH.read_text())["graph"]

dihedrals = []

seen = set()

for b,data_b in graph.items():

    for c in data_b["neighbors"]:

        data_c = graph[c]

        for a in data_b["neighbors"]:

            if a == c:
                continue

            for d in data_c["neighbors"]:

                if d == b:
                    continue

                key = (a,b,c,d)

                rev = (d,c,b,a)

                if key in seen or rev in seen:
                    continue

                seen.add(key)

                dihedrals.append({

                    "atom1":a,

                    "atom2":b,

                    "atom3":c,

                    "atom4":d

                })

print("="*100)
print("DAY041 / PHASE1B-A8")
print("DIHEDRAL ENUMERATION")
print("="*100)
print()

print("[1] DIHEDRALS")
print(len(dihedrals))

outfile = RUN/"PHASE1B_DIHEDRALS.json"

outfile.write_text(

    json.dumps({

        "timestamp":utc(),

        "dihedral_count":len(dihedrals),

        "dihedrals":dihedrals

    },indent=2)

)

print()

print("[2] OUTPUT")
print(outfile)

print()

print("[3] DECISION")
print("PHASE1B_DIHEDRALS_CREATED")
