#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import json
import collections
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

GRAPH = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_molecular_graph"
    / "PHASE1B_MOLECULAR_GRAPH.json"
)

DIH = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_dihedrals"
    / "PHASE1B_DIHEDRALS.json"
)

RUN = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_dihedral_validation"
)

RUN.mkdir(parents=True, exist_ok=True)


def utc():
    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00", "Z")
    )

graph = json.loads(GRAPH.read_text())["graph"]
dihedrals = json.loads(DIH.read_text())["dihedrals"]

duplicates = 0
seen = set()

hist = collections.Counter()

for d in dihedrals:

    key = (
        d["atom1"],
        d["atom2"],
        d["atom3"],
        d["atom4"]
    )

    rev = tuple(reversed(key))

    if key in seen or rev in seen:
        duplicates += 1

    seen.add(key)

    t = (
        graph[d["atom1"]]["type"],
        graph[d["atom2"]]["type"],
        graph[d["atom3"]]["type"],
        graph[d["atom4"]]["type"]
    )

    hist["-".join(t)] += 1

print("="*100)
print("DAY041 / PHASE1B-A9")
print("DIHEDRAL VALIDATION")
print("="*100)
print()

print("[1] DIHEDRALS")
print(len(dihedrals))

print()

print("[2] DUPLICATES")
print(duplicates)

print()

print("[3] UNIQUE CHEMICAL TYPES")

for k,v in sorted(hist.items()):
    print(f"{k:24s} {v}")

report = {

    "timestamp": utc(),

    "dihedral_count": len(dihedrals),

    "duplicates": duplicates,

    "chemical_histogram": dict(hist),

    "decision":
        "PASS" if duplicates == 0 else "REVIEW"

}

outfile = RUN / "PHASE1B_DIHEDRAL_VALIDATION.json"

outfile.write_text(
    json.dumps(report, indent=2)
)

print()

print("[4] OUTPUT")
print(outfile)

print()

print("[5] DECISION")
print(report["decision"])
