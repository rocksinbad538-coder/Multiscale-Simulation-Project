#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import hashlib
import json
import csv
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

INPUT = (
ROOT
/"runs"
/"phase1A"
/"day040_phase1A_semantic_forcefield_audit"
/"PHASE1B_PARAMETER_COMPLETENESS_CERTIFICATE.json"
)

assert INPUT.exists()

TRANSFERABILITY = (
ROOT
/"runs"
/"phase1A"
/"day038_resp_stage1_executions"
/"resp_stage1_upper_v7a_r1_20260803T202335Z"
/
"QM_F06_UPPER_V7A_R1_RESP_STAGE1_TRANSFERABILITY.csv"
)

assert TRANSFERABILITY.exists()


PARENT_MAPPING = (
ROOT
/"runs"
/"phase1A"
/"day040_phase1A_G_hydrogen_augmentation_design"
/
"QM_F06_UPPER_V7A_R1_HYDROGEN_PARENT_MAPPING.csv"
)

assert PARENT_MAPPING.exists()


RUN = (
ROOT
/"runs"
/"phase1A"
/"day040_phase1A_parameter_mapping"
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
print("DAY040 / D040-A8")
print("PARAMETER MAPPING DESIGN")
print("="*100)
print()
mapping = []

parent_lookup = {}

with open(PARENT_MAPPING) as f:

    reader = csv.DictReader(f)

    for row in reader:

        parent_lookup[
            row["hydrogen_atom_id"]
        ] = row

with open(TRANSFERABILITY) as f:

    reader = csv.DictReader(f)

    for row in reader:

        atom_id = row["atom_id"]

        element = row["element"]

        role = row["atom_role"]

        parent_element = ""

        proposed_type = element

        if atom_id in parent_lookup:

            parent_element = parent_lookup[
                atom_id
            ]["parent_element"]

            if parent_element == "B":

                proposed_type = "HB"

            elif parent_element == "N":

                proposed_type = "HN"

        record = dict(row)

        record["parent_element"] = parent_element

        record["proposed_forcefield_type"] = proposed_type

        record["RESP_stage1_charge_e"] = row["RESP_stage1_charge_e_float"]

        mapping.append(record)

print("[1] TRACEABILITY TABLE")

print("rows =", len(mapping))

print()

counts = {}

for m in mapping:

    t = m["proposed_forcefield_type"]

    counts[t] = counts.get(t, 0) + 1

print("[2] PROPOSED FORCE-FIELD TYPES")

for k in sorted(counts):

    print(f"{k:>4s} : {counts[k]}")

print()
csv_file=RUN/"PHASE1B_PARAMETER_MAPPING.csv"

with open(csv_file,"w",newline="") as f:

    writer=csv.DictWriter(
        f,
        fieldnames=list(mapping[0].keys())
    )

    writer.writeheader()

    writer.writerows(mapping)

report={

"timestamp":utc(),

"mapped_atoms":len(mapping),

"mapping_generated":True,

"HB_detected":any(
x["proposed_forcefield_type"]=="HB"
for x in mapping
),

"HN_detected":any(
x["proposed_forcefield_type"]=="HN"
for x in mapping
),

"decision":
"D040_A8_PARAMETER_MAPPING_COMPLETE",

"phase1B_mapping_ready":True

}

json_file=RUN/"PHASE1B_PARAMETER_MAPPING.json"

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
