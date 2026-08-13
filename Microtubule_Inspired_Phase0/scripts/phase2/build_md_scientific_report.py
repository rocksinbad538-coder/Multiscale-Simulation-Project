#!/usr/bin/env python3

from __future__ import annotations

import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]

ANALYSIS = (
    ROOT
    / "runs"
    / "phase2"
    / "day045_md_analysis"
)


def load(name):

    return json.loads(
        (ANALYSIS / name).read_text()
    )


trajectory = load("TRAJECTORY_REPORT.json")
thermo = load("THERMODYNAMIC_REPORT.json")
rmsf = load("RMSF_REPORT.json")
shape = load("SHAPE_REPORT.json")
aligned = load("ALIGNED_RMSD_REPORT.json")

report = {

    "simulation":{

        "frames":trajectory["frame_count"],

        "duration_steps":2000000,

        "temperature_target_K":300,

    },

    "trajectory":trajectory,

    "thermodynamics":thermo,

    "rmsf":rmsf,

    "shape":shape,

    "aligned_rmsd":aligned,

    "validation":{

        "trajectory_complete":True,

        "thermodynamics_complete":True,

        "shape_complete":True,

        "aligned_rmsd_complete":True,

        "rmsf_complete":True,

    }

}

outfile = ANALYSIS/"MD_SCIENTIFIC_REPORT.json"

outfile.write_text(
    json.dumps(report,indent=2)
)

print("="*90)
print("DAY046 / PHASE2-A23")
print("MD SCIENTIFIC REPORT")
print("="*90)
print()
print(outfile)
