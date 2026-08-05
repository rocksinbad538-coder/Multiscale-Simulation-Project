#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import json
import collections
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

TOPOLOGY = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_internal_topology"
    / "PHASE1B_INTERNAL_TOPOLOGY.json"
)

ANGLES = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_angles"
    / "PHASE1B_ANGLES.json"
)

RUN = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_angle_typing"
)

RUN.mkdir(parents=True, exist_ok=True)


def utc():
    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )


atoms = {}

for a in json.loads(TOPOLOGY.read_text())["atoms"]:
    atoms[a["atom_id"]] = a["forcefield_type"]

angles = json.loads(ANGLES.read_text())["angles"]

typed = []

hist = collections.Counter()

for ang in angles:

    t1 = atoms[ang["atom1"]]
    t2 = atoms[ang["atom2"]]
    t3 = atoms[ang["atom3"]]

    label = f"{t1}-{t2}-{t3}"

    hist[label] += 1

    typed.append({

        **ang,

        "angle_type":label

    })

print("="*100)
print("DAY041 / PHASE1B-A5")
print("CHEMICAL ANGLE TYPING")
print("="*100)
print()

print("[1] UNIQUE ANGLE TYPES")

for k in sorted(hist):
    print(f"{k:18s} {hist[k]}")

outfile = RUN/"PHASE1B_TYPED_ANGLES.json"

outfile.write_text(

json.dumps({

"timestamp":utc(),

"typed_angles":typed,

"histogram":dict(hist)

},indent=2)

)

print()

print("[2] OUTPUT")

print(outfile)

print()

print("[3] DECISION")

print("PHASE1B_ANGLE_TYPES_CREATED")
