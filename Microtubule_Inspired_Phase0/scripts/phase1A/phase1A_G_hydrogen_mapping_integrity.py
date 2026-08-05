#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import csv
import json
import hashlib
import datetime
from collections import Counter

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
"day040_phase1A_hydrogen_integrity_audit"
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
print("DAY040 / D040-A10A")
print("HYDROGEN MAPPING INTEGRITY AUDIT")
print("="*100)
print()
rows=[]

hydrogen_ids=[]

with open(INPUT) as f:

    reader=csv.DictReader(f)

    for row in reader:

        rows.append(row)

        if row["element"]=="H":

            hydrogen_ids.append(row["atom_id"])

counter=Counter(hydrogen_ids)

duplicates={

k:v

for k,v in counter.items()

if v>1

}

print("[1] HYDROGEN COUNTS")

print()

print("total_H_rows =",len(hydrogen_ids))

print("unique_H =",len(counter))

print("duplicate_atom_ids =",len(duplicates))

print()

if duplicates:

    print("[2] DUPLICATED HYDROGENS")

    print()

    for atom_id,n in sorted(duplicates.items()):

        print(atom_id,"count=",n)

print()
report={

"timestamp":utc(),

"total_hydrogen_rows":len(hydrogen_ids),

"unique_hydrogen_atom_ids":len(counter),

"duplicate_count":len(duplicates),

"duplicates":duplicates,

"decision":

"D040_A10A_HYDROGEN_MAPPING_INTEGRITY_COMPLETE"

}

json_file=RUN/"HYDROGEN_MAPPING_INTEGRITY.json"

json_file.write_text(

json.dumps(

report,

indent=2,

)

)

print("[3] OUTPUT")

print(json_file)

print()

print("[4] DECISION")

print(report["decision"])
