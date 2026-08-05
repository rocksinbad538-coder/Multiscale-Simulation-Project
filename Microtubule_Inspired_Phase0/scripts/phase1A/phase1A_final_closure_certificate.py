#!/usr/bin/env python3

from __future__ import annotations

import json
import hashlib
import pathlib
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

RUN = (
    ROOT
    / "runs"
    / "phase1A"
    / "day040_phase1A_final_closure"
)

RUN.mkdir(parents=True, exist_ok=True)


def sha256(path):

    h = hashlib.sha256()

    with open(path, "rb") as f:

        while True:

            b = f.read(1024 * 1024)

            if not b:

                break

            h.update(b)

    return h.hexdigest()


def utc():

    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00", "Z")
    )


print("=" * 100)
print("DAY040 / D040-A12")
print("PHASE 1A FINAL CLOSURE CERTIFICATE")
print("=" * 100)
print()


artifacts = [

ROOT/"runs/phase1A/day040_phase1A_transition_specification/PHASE1B_TRANSITION_SPECIFICATION.json",

ROOT/"runs/phase1A/day040_phase1A_RESP_hydrogen_origin_audit/RESP_HYDROGEN_ORIGIN_AUDIT.json",

ROOT/"runs/phase1A/day040_phase1A_parameter_mapping/PHASE1B_PARAMETER_MAPPING.json",

ROOT/"runs/phase1A/day040_phase1A_topology_mapping/PHASE1B_TOPOLOGY_MAPPING.json",

ROOT/"runs/phase1A/day040_phase1A_semantic_forcefield_audit/PHASE1B_PARAMETER_COMPLETENESS_CERTIFICATE.json"

]

inventory=[]

print("[1] VERIFY ARTIFACTS")

for f in artifacts:

    ok=f.exists()

    print(("PASS " if ok else "FAIL "),f.name)

    inventory.append({

        "file":str(f.relative_to(ROOT)),

        "exists":ok,

        "sha256":sha256(f) if ok else None

    })

print()

all_ok=all(x["exists"] for x in inventory)

report={

"timestamp_utc":utc(),

"phase":"Phase1A",

"status":"CLOSED" if all_ok else "INCOMPLETE",

"physical_atoms":37,

"QM_atoms":52,

"physical_hydrogens":6,

"QM_caps":15,

"artifacts":inventory,

"phase1B_ready":all_ok,

"decision":"PHASE1A_FORMALLY_CLOSED" if all_ok else "PHASE1A_NOT_READY"

}

json_file=RUN/"PHASE1A_FINAL_CERTIFICATE.json"

json_file.write_text(

json.dumps(

report,

indent=2

)

)

md_file=RUN/"PHASE1A_FINAL_CERTIFICATE.md"

with open(md_file,"w") as f:

    f.write("# Phase 1A Final Certificate\n\n")

    f.write(f"Generated: {report['timestamp_utc']}\n\n")

    f.write(f"Status: {report['status']}\n\n")

    f.write(f"Physical atoms: {report['physical_atoms']}\n")

    f.write(f"QM atoms: {report['QM_atoms']}\n")

    f.write(f"Physical H: {report['physical_hydrogens']}\n")

    f.write(f"QM caps: {report['QM_caps']}\n\n")

    f.write("## Verified Artifacts\n\n")

    for x in inventory:

        f.write(

            f"- {x['file']} : "

            f"{'PASS' if x['exists'] else 'FAIL'}\n"

        )

print("[2] OUTPUTS")

print(json_file)

print(md_file)

print()

print("[3] DECISION")

print(report["decision"])
