#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import json
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

RUN = (
    ROOT
    / "runs"
    / "phase1A"
    / "day040_phase1_readiness_verification"
)

RUN.mkdir(parents=True, exist_ok=True)


def utc():

    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )


required = {

"Phase1A certificate":
ROOT/"runs/phase1A/day040_phase1A_final_closure/PHASE1A_FINAL_CERTIFICATE.json",

"Phase1B manifest":
ROOT/"runs/phase1A/day040_phase1B_input_manifest/PHASE1B_INPUT_MANIFEST.json",

"Transition specification":
ROOT/"runs/phase1A/day040_phase1A_transition_specification/PHASE1B_TRANSITION_SPECIFICATION.json",

"Parameter mapping":
ROOT/"runs/phase1A/day040_phase1A_parameter_mapping/PHASE1B_PARAMETER_MAPPING.csv",

"Topology mapping":
ROOT/"runs/phase1A/day040_phase1A_topology_mapping/PHASE1B_TOPOLOGY_MAPPING.csv",

"Hydrogen taxonomy":
ROOT/"runs/phase1A/day040_phase1A_hydrogen_taxonomy/HYDROGEN_TAXONOMY.csv",

"RESP hydrogen audit":
ROOT/"runs/phase1A/day040_phase1A_RESP_hydrogen_origin_audit/RESP_HYDROGEN_ORIGIN_AUDIT.json",

"Force-field semantic audit":
ROOT/"runs/phase1A/day040_phase1A_semantic_forcefield_audit/PHASE1B_PARAMETER_COMPLETENESS_CERTIFICATE.json",

"Force field":
ROOT/"references/force_fields/Rajan_JPCL_2018_hBN_Functionalized/hBN_functionalized-FFTOOLS.ff"

}

print("="*100)
print("DAY040 / D040-A14")
print("PHASE 1 READINESS VERIFICATION")
print("="*100)
print()

results=[]

for name,path in required.items():

    ok=path.exists()

    results.append(ok)

    print(f"{name:35s} {'PASS' if ok else 'FAIL'}")

ready=all(results)

report={

"timestamp_utc":utc(),

"phase1A_closed":ready,

"phase1B_authorized":ready,

"verification":{

k:v.exists()

for k,v in required.items()

},

"decision":

"BEGIN_PHASE1B_AUTHORIZED"

if ready

else

"PHASE1B_BLOCKED"

}

json_file=RUN/"PHASE1_READINESS_REPORT.json"

json_file.write_text(

json.dumps(

report,

indent=2

)

)

print()

print("[OUTPUT]")

print(json_file)

print()

print("[DECISION]")

print(report["decision"])
