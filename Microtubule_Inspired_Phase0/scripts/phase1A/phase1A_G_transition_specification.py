#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import csv
import json
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

INPUT = (
ROOT
/"runs"
/
"phase1A"
/
"day040_phase1A_RESP_hydrogen_origin_audit"
/
"RESP_HYDROGEN_ORIGIN_TABLE.csv"
)

RUN = (
ROOT
/"runs"
/
"phase1A"
/
"day040_phase1A_transition_specification"
)

RUN.mkdir(parents=True, exist_ok=True)

def utc():
    return (
        datetime.datetime.now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )

print("="*100)
print("DAY040 / D040-A11")
print("RESP TO FORCE-FIELD TRANSITION SPECIFICATION")
print("="*100)
print()
physical=[]

qm_caps=[]

with open(INPUT) as f:

    reader=csv.DictReader(f)

    for row in reader:

        if row["transfer_status"]=="TRANSFERABLE_REAL_ATOM":

            physical.append(row)

        else:

            qm_caps.append(row)

print("[1] TRANSITION SUMMARY")

print()

print("physical hydrogens =",len(physical))

print("QM caps =",len(qm_caps))

print()
phys_csv=RUN/"PHASE1B_PHYSICAL_HYDROGENS.csv"

with open(phys_csv,"w",newline="") as f:

    writer=csv.DictWriter(
        f,
        fieldnames=physical[0].keys()
    )

    writer.writeheader()
    writer.writerows(physical)

caps_csv=RUN/"PHASE1B_QM_CAPS.csv"

with open(caps_csv,"w",newline="") as f:

    writer=csv.DictWriter(
        f,
        fieldnames=qm_caps[0].keys()
    )

    writer.writeheader()
    writer.writerows(qm_caps)

report={

"timestamp":utc(),

"physical_hydrogens":len(physical),

"qm_caps":len(qm_caps),

"phase1A_complete":True,

"phase1B_ready":True,

"decision":"D040_A11_PHASE1B_TRANSITION_READY"

}

json_file=RUN/"PHASE1B_TRANSITION_SPECIFICATION.json"

json_file.write_text(

json.dumps(report,indent=2)

)

print("[2] OUTPUTS")

print(phys_csv)
print(caps_csv)
print(json_file)

print()

print("[3] DECISION")

print(report["decision"])
