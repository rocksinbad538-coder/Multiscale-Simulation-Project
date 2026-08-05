#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import csv
import json
import hashlib
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

RUN = (
ROOT
/"runs"
/
"phase1A"
/
"day040_phase1A_topology_mapping"
)

RUN.mkdir(
parents=True,
exist_ok=True,
)

def utc():

    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )

print("="*100)
print("DAY040 / D040-A9")
print("TOPOLOGY MAPPING SPECIFICATION")
print("="*100)
print()
rows=[]

with open(INPUT) as f:

    reader=csv.DictReader(f)

    for r in reader:

        keep=False

        reason=""

        ff=r["proposed_forcefield_type"]

        role=r["atom_role"]

        if ff in ("HB","HN"):

            keep=True
            reason="FUNCTIONAL_EDGE"

        elif r["element"]!="H":

            keep=True
            reason="FRAMEWORK"

        else:

            keep=False
            reason="NON_FUNCTIONAL_H"

        r["topology_keep"]=keep
        r["topology_reason"]=reason

        rows.append(r)

print("[1] TOPOLOGY CLASSIFICATION")

print("atoms =",len(rows))

print()

keep=sum(
r["topology_keep"]
for r in rows
)

print("kept =",keep)

print("discarded =",len(rows)-keep)

print()
csv_file=RUN/"PHASE1B_TOPOLOGY_MAPPING.csv"

with open(csv_file,"w",newline="") as f:

    writer=csv.DictWriter(
        f,
        fieldnames=rows[0].keys()
    )

    writer.writeheader()

    writer.writerows(rows)

report={

"timestamp":utc(),

"total_atoms":len(rows),

"topology_atoms":keep,

"mapping_complete":True,

"decision":
"D040_A9_TOPOLOGY_MAPPING_COMPLETE",

"phase1B_topology_ready":True,

}

json_file=RUN/"PHASE1B_TOPOLOGY_MAPPING.json"

json_file.write_text(

json.dumps(
report,
indent=2,
)

)

print("[2] OUTPUTS")

print(csv_file)

print(json_file)

print()

print("[3] DECISION")

print(report["decision"])
