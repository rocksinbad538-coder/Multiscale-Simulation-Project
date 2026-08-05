#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import csv
import json
import datetime
from collections import Counter

ROOT = pathlib.Path(__file__).resolve().parents[2]

TRANSFER = (
ROOT
/"runs"
/
"phase1A"
/
"day038_resp_stage1_executions"
/
"resp_stage1_upper_v7a_r1_20260803T202335Z"
/
"QM_F06_UPPER_V7A_R1_RESP_STAGE1_TRANSFERABILITY.csv"
)

RUN = (
ROOT
/"runs"
/
"phase1A"
/
"day040_phase1A_RESP_hydrogen_origin_audit"
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
print("DAY040 / D040-A10C")
print("RESP HYDROGEN ORIGIN AUDIT")
print("="*100)
print()
rows=[]

role_counter=Counter()

node_counter=Counter()

artificial_counter=Counter()

transfer_counter=Counter()

with open(TRANSFER) as f:

    reader=csv.DictReader(f)

    for row in reader:

        if row["element"]!="H":

            continue

        rows.append(row)

        role_counter[row["atom_role"]]+=1

        node_counter[row["node_type"]]+=1

        artificial_counter[row["artificial_cap"]]+=1

        transfer_counter[row["transfer_status"]]+=1

print("[1] TOTAL HYDROGENS")

print(len(rows))

print()

print("[2] ROLE DISTRIBUTION")

for k,v in sorted(role_counter.items()):

    print(f"{v:3d}   {k}")

print()

print("[3] NODE TYPES")

for k,v in sorted(node_counter.items()):

    print(f"{v:3d}   {k}")

print()

print("[4] ARTIFICIAL CAP")

for k,v in sorted(artificial_counter.items()):

    print(f"{k:8s} {v}")

print()

print("[5] TRANSFER STATUS")

for k,v in sorted(transfer_counter.items()):

    print(f"{v:3d}   {k}")

print()
csv_file=RUN/"RESP_HYDROGEN_ORIGIN_TABLE.csv"

with open(csv_file,"w",newline="") as f:

    writer=csv.DictWriter(
        f,
        fieldnames=rows[0].keys()
    )

    writer.writeheader()

    writer.writerows(rows)

report={

"timestamp":utc(),

"hydrogen_count":len(rows),

"role_distribution":dict(role_counter),

"node_distribution":dict(node_counter),

"artificial_cap_distribution":dict(artificial_counter),

"transfer_status_distribution":dict(transfer_counter),

"decision":"D040_A10C_RESP_HYDROGEN_ORIGIN_AUDIT_COMPLETE"

}

json_file=RUN/"RESP_HYDROGEN_ORIGIN_AUDIT.json"

json_file.write_text(

json.dumps(

report,

indent=2,

)

)

print("[6] OUTPUTS")

print(csv_file)

print(json_file)

print()

print("[7] DECISION")

print(report["decision"])
