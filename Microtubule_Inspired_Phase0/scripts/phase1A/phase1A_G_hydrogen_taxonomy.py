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
"day040_phase1A_hydrogen_taxonomy"
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
print("DAY040 / D040-A10")
print("HYDROGEN TAXONOMY FINALIZATION")
print("="*100)
print()
rows=[]

counts={

"FUNCTIONAL_EDGE":0,

"ORIGINAL_FRAGMENT":0,

"RESP_TEMPORARY":0,

}

with open(INPUT) as f:

    reader=csv.DictReader(f)

    for row in reader:

        if row["element"]!="H":

            continue

        ff=row["proposed_forcefield_type"]

        role=row["atom_role"].upper()

        if ff in ("HB","HN"):

            taxonomy="FUNCTIONAL_EDGE"

        elif "ORIGINAL_FRAGMENT" in role:

            taxonomy="ORIGINAL_FRAGMENT"

        else:

            taxonomy="RESP_TEMPORARY"

        row["hydrogen_taxonomy"]=taxonomy

        counts[taxonomy]+=1

        rows.append(row)

print("[1] HYDROGEN TAXONOMY")

print()

for k,v in counts.items():

    print(f"{k:20s} {v}")

print()

print("total =",len(rows))

print()
csv_file=RUN/"HYDROGEN_TAXONOMY.csv"

with open(csv_file,"w",newline="") as f:

    writer=csv.DictWriter(

        f,

        fieldnames=rows[0].keys()

    )

    writer.writeheader()

    writer.writerows(rows)

report={

"timestamp":utc(),

"functional_edge":counts["FUNCTIONAL_EDGE"],

"original_fragment":counts["ORIGINAL_FRAGMENT"],

"resp_temporary":counts["RESP_TEMPORARY"],

"decision":

"D040_A10_HYDROGEN_TAXONOMY_COMPLETE",

"phase1A_scientific_taxonomy_complete":True

}

json_file=RUN/"HYDROGEN_TAXONOMY.json"

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
