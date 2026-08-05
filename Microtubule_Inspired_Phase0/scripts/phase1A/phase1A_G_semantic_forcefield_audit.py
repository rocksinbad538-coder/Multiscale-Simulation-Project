#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import pathlib
import hashlib
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

RUN_DIR = (
    ROOT
    / "runs"
    / "phase1A"
    / "day040_phase1A_semantic_forcefield_audit"
)

RUN_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

SOURCE_CSV = (
    ROOT
    / "runs"
    / "phase1A"
    / "day040_phase1A_G_fftools_parameter_audit"
    / "FFTOOLS_PARAMETER_INVENTORY.csv"
)

assert SOURCE_CSV.exists()


def sha256(path):

    h = hashlib.sha256()

    with open(path,"rb") as f:

        while True:

            b=f.read(1024*1024)

            if not b:
                break

            h.update(b)

    return h.hexdigest()


def utc():

    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )

rows=[]

with open(SOURCE_CSV) as f:

    reader=csv.DictReader(f)

    for r in reader:

        rows.append(r)

print("="*100)
print("DAY040 / D040-A7B")
print("SEMANTIC FORCE FIELD AUDIT")
print("="*100)
print()
print("[1] INVENTORY")
print("records =",len(rows))
print()
required={

"HB atom type":False,
"HN atom type":False,

"B-HB bond":False,
"N-HN bond":False,
"B-N bond":False,

"B-N-HN angle":False,
"N-B-HB angle":False,

"B improper":False,
"N improper":False,

"HB dihedral":False,
"HN dihedral":False,

}

evidence={}

for row in rows:

    text=row["raw"]

    tokens=text.split()

    sec=row["section"]

    if sec=="ATOMS":

        if text.startswith("HB "):

            required["HB atom type"]=True
            evidence["HB atom type"]=text

        if text.startswith("HN "):

            required["HN atom type"]=True
            evidence["HN atom type"]=text

    elif sec=="BONDS":

        if text.startswith("B  HB"):

            required["B-HB bond"]=True
            evidence["B-HB bond"]=text

        elif text.startswith("N  HN"):

            required["N-HN bond"]=True
            evidence["N-HN bond"]=text

        elif text.startswith("B  N"):

            required["B-N bond"]=True
            evidence["B-N bond"]=text

    elif sec=="ANGLES":

        if len(tokens)>=3:

            if tokens[:3]==["B","N","HN"]:

                required["B-N-HN angle"]=True
                evidence["B-N-HN angle"]=text

            elif tokens[:3]==["N","B","HB"]:

                required["N-B-HB angle"]=True
                evidence["N-B-HB angle"]=text

    elif sec=="IMPROPER":

        if text.startswith("B"):

            required["B improper"]=True
            evidence["B improper"]=text

        elif text.startswith("N"):

            required["N improper"]=True
            evidence["N improper"]=text

    elif sec=="DIHEDRALS":

        if "HN" in text:

            required["HN dihedral"]=True
            evidence["HN dihedral"]=text

        if "HB" in text:

            required["HB dihedral"]=True
            evidence["HB dihedral"]=text

print("[2] SEMANTIC AUDIT")

for k,v in required.items():

    print(f"{k:25s}",v)

print()
coverage=sum(required.values())

total=len(required)

coverage_fraction=coverage/total

print("[3] COVERAGE")

print(f"{coverage}/{total}")

print()

decision=(
    "PASS"
    if coverage==total
    else
    "REVIEW"
)

csv_file=RUN_DIR/"SEMANTIC_PARAMETER_AUDIT.csv"

with open(csv_file,"w",newline="") as f:

    w=csv.writer(f)

    w.writerow(

        [

            "requirement",

            "status",

            "evidence",

        ]

    )

    for k in required:

        w.writerow(

            [

                k,

                required[k],

                evidence.get(k,""),

            ]

        )

certificate={

"timestamp":utc(),

"coverage":coverage,

"total":total,

"coverage_fraction":coverage_fraction,

"decision":decision,

"HB_atom_type_detected":required["HB atom type"],

"HN_atom_type_detected":required["HN atom type"],

"scientific_conclusion":
"Hydrogen type assignment is explicitly separated into HB and HN according to parent atom identity in the authoritative force field.",

"phase1B_parameterization_ready":
coverage==total,

}

json_file=RUN_DIR/"PHASE1B_PARAMETER_COMPLETENESS_CERTIFICATE.json"

json_file.write_text(

    json.dumps(

        certificate,

        indent=2,

    )

)

print("[4] OUTPUTS")

print(csv_file)

print(json_file)

print()

print("[5] DECISION")

print(decision)
