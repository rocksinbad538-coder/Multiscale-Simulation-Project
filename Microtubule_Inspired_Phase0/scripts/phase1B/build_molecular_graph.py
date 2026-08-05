#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import json
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

TOPOLOGY = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_internal_topology"
    / "PHASE1B_INTERNAL_TOPOLOGY.json"
)

BONDS = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_bond_connectivity"
    / "PHASE1B_BOND_CONNECTIVITY.json"
)

RUN = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_molecular_graph"
)

RUN.mkdir(parents=True, exist_ok=True)


def utc():
    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00", "Z")
    )


topology = json.loads(TOPOLOGY.read_text())
bond_data = json.loads(BONDS.read_text())

graph = {}

for atom in topology["atoms"]:

    graph[atom["atom_id"]] = {
        "element": atom["element"],
        "type": atom["forcefield_type"],
        "neighbors": []
    }

for bond in bond_data["bonds"]:

    a = bond["atom1"]
    b = bond["atom2"]

    graph[a]["neighbors"].append(b)
    graph[b]["neighbors"].append(a)

coordination = {}

for atom_id in sorted(graph):

    coordination[atom_id] = len(graph[atom_id]["neighbors"])

print("="*100)
print("DAY041 / PHASE1B-A3")
print("MOLECULAR GRAPH")
print("="*100)
print()

print("[1] ATOMS")
print(len(graph))
print()

print("[2] MAX COORDINATION")
print(max(coordination.values()))
print()

print("[3] MIN COORDINATION")
print(min(coordination.values()))
print()

print("[4] COORDINATION HISTOGRAM")

hist = {}

for c in coordination.values():
    hist[c] = hist.get(c,0)+1

for k in sorted(hist):
    print(f"{k:2d} : {hist[k]}")

report = {

    "timestamp": utc(),

    "graph": graph,

    "coordination_histogram": hist,

    "status": "MOLECULAR_GRAPH_CREATED"

}

outfile = RUN / "PHASE1B_MOLECULAR_GRAPH.json"

outfile.write_text(
    json.dumps(
        report,
        indent=2
    )
)

print()
print("[5] OUTPUT")
print(outfile)
print()

print("[6] DECISION")
print("PHASE1B_GRAPH_CREATED")
