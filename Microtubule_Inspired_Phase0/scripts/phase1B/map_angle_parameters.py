#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import json
import collections
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

ANGLES = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_angle_typing"
    / "PHASE1B_TYPED_ANGLES.json"
)

RUN = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_angle_parameter_mapping"
)

RUN.mkdir(parents=True, exist_ok=True)


def utc():
    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00", "Z")
    )


canonical = {

    "B-N-B":"B-N-B",

    "N-B-N":"N-B-N",

    "B-N-HN":"B-N-HN",

    "HN-N-B":"B-N-HN",

    "N-B-HB":"N-B-HB",

    "HB-B-N":"N-B-HB"

}

angles = json.loads(ANGLES.read_text())["typed_angles"]

mapped=[]

hist=collections.Counter()

for a in angles:

    fftype=canonical[a["angle_type"]]

    hist[fftype]+=1

    mapped.append({

        **a,

        "forcefield_parameter":fftype

    })

print("="*100)
print("DAY041 / PHASE1B-A6")
print("CANONICAL ANGLE PARAMETER MAPPING")
print("="*100)
print()

print("[1] FORCE-FIELD ANGLE TYPES")

for k in sorted(hist):

    print(f"{k:12s} {hist[k]}")

outfile=RUN/"PHASE1B_CANONICAL_ANGLES.json"

outfile.write_text(

json.dumps({

"timestamp":utc(),

"angles":mapped,

"histogram":dict(hist)

},indent=2)

)

print()

print("[2] OUTPUT")

print(outfile)

print()

print("[3] DECISION")

print("PHASE1B_CANONICAL_ANGLE_MAPPING_COMPLETE")
