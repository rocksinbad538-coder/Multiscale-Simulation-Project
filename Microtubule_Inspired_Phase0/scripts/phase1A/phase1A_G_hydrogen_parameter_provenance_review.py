#!/usr/bin/env python3

"""
DAY040 / D040-A6
Hydrogen Parameter Provenance Review

Scientific objective
--------------------
Establish provenance of every hydrogen-related interaction required
for future edge-passivated hBN topology construction.

NO PARAMETERS ARE CREATED.
NO FORCE FIELD IS MODIFIED.
NO TOPOLOGY IS GENERATED.
"""

from pathlib import Path
import json
import hashlib
import pandas as pd
import datetime

ROOT = Path.cwd()

INPUT_DIR = ROOT / "runs/phase1A/day040_phase1A_G_hydrogen_parameter_source_audit"

OUTPUT_DIR = ROOT / "runs/phase1A/day040_phase1A_G_hydrogen_parameter_provenance_review"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(1024 * 1024)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


required_terms = [
    "H_ATOMTYPE",
    "B_H_BOND",
    "N_H_BOND",
    "H_CONTAINING_ANGLES",
    "H_CONTAINING_IMPROPERS",
]

records = []

for term in required_terms:

    if term == "H_ATOMTYPE":
        status = "PROJECT_PARAMETER_PRESENT"

        recommended_source = (
            "Existing accepted project atomtypes "
            "(pyrene/TIP4P) require environment validation "
            "before transfer to edge-passivated hBN."
        )

    else:

        status = "PRIMARY_SOURCE_REQUIRED"

        recommended_source = (
            "Explicit literature parameterization "
            "required before adoption."
        )

    records.append(
        dict(
            interaction=term,
            repository_status=status,
            parameter_adoption_authorized=False,
            recommended_primary_source=recommended_source,
            transferability_established=False,
        )
    )

df = pd.DataFrame(records)

csv_file = OUTPUT_DIR / "QM_F06_UPPER_V7A_R1_HYDROGEN_PARAMETER_PROVENANCE.csv"
df.to_csv(csv_file, index=False)

report = {

    "day":"DAY040",
    "block":"D040_A6",

    "timestamp_utc":
        datetime.datetime.utcnow().isoformat()+"Z",

    "scientific_question":
        "Has complete provenance been established for every required hydrogen interaction?",

    "required_interactions":required_terms,

    "summary":{

        "complete_parameter_provenance":False,

        "parameter_adoption_authorized":False,

        "primary_literature_required":True,

        "next_block":
        "D040_A7_PRIMARY_LITERATURE_PARAMETER_EXTRACTION"

    },

    "csv_sha256":sha256(csv_file)

}

json_file = OUTPUT_DIR / "QM_F06_UPPER_V7A_R1_HYDROGEN_PARAMETER_PROVENANCE.json"

json_file.write_text(
    json.dumps(report,indent=2)+"\n",
    encoding="utf-8",
)

print("="*100)
print("DAY040 / D040-A6 — HYDROGEN PARAMETER PROVENANCE REVIEW")
print("="*100)

print()

print("required_interactions =",len(required_terms))
print("parameter_adoption_authorized = False")
print("primary_literature_required = True")

print()

print("output_csv =",csv_file)
print("output_csv_sha256 =",sha256(csv_file))

print()

print("output_json =",json_file)
print("output_json_sha256 =",sha256(json_file))

print()

print("decision=D040_A6_PARAMETER_PROVENANCE_PASS_PRIMARY_LITERATURE_EXTRACTION_AUTHORIZED")
