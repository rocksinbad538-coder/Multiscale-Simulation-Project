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
    / "day040_phase1B_input_manifest"
)

RUN.mkdir(parents=True, exist_ok=True)


def utc():
    return (
        datetime.datetime.now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )

manifest = {

    "timestamp_utc": utc(),

    "phase": "Phase1B",

    "coordinates":
    "runs/phase1A/day035_qm_f06_upper_v7a_r1_coordinate_adoption/",

    "RESP_transferability":
    "runs/phase1A/day038_resp_stage1_executions/resp_stage1_upper_v7a_r1_20260803T202335Z/QM_F06_UPPER_V7A_R1_RESP_STAGE1_TRANSFERABILITY.csv",

    "parameter_mapping":
    "runs/phase1A/day040_phase1A_parameter_mapping/PHASE1B_PARAMETER_MAPPING.csv",

    "topology_mapping":
    "runs/phase1A/day040_phase1A_topology_mapping/PHASE1B_TOPOLOGY_MAPPING.csv",

    "hydrogen_taxonomy":
    "runs/phase1A/day040_phase1A_hydrogen_taxonomy/HYDROGEN_TAXONOMY.csv",

    "transition_specification":
    "runs/phase1A/day040_phase1A_transition_specification/PHASE1B_TRANSITION_SPECIFICATION.json",

    "force_field":
    "references/force_fields/Rajan_JPCL_2018_hBN_Functionalized/hBN_functionalized-FFTOOLS.ff",

    "physical_atoms":37,

    "QM_atoms":52,

    "physical_hydrogens":6,

    "QM_caps":15,

    "status":"READY_FOR_PHASE1B"

}

json_file = RUN/"PHASE1B_INPUT_MANIFEST.json"

json_file.write_text(

    json.dumps(
        manifest,
        indent=2
    )

)

md_file = RUN/"PHASE1B_INPUT_MANIFEST.md"

with open(md_file,"w") as f:

    f.write("# Phase 1B Input Manifest\n\n")

    for k,v in manifest.items():

        f.write(f"- **{k}** : {v}\n")

print("="*100)
print("DAY040 / D040-A13")
print("PHASE1B INPUT MANIFEST")
print("="*100)
print()

print(json_file)
print(md_file)

print()

print("DECISION")

print("PHASE1B_INPUT_MANIFEST_COMPLETE")
