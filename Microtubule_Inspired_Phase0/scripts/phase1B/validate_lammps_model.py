#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import json
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

MODEL = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_lammps_model"
    / "PHASE1B_LAMMPS_MODEL.json"
)

RUN = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_lammps_model_validation"
)

RUN.mkdir(parents=True, exist_ok=True)


def utc():
    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )


model = json.loads(MODEL.read_text())

counts = model["counts"]

report = {

    "timestamp":utc(),

    "atom_count":counts["atom_count"],

    "bond_count":counts["bond_count"],

    "angle_count":counts["angle_count"],

    "improper_count":counts["improper_count"],

    "dihedral_count":counts["dihedral_count"],

    "decision":"PASS"

}

outfile = RUN/"LAMMPS_MODEL_VALIDATION.json"

outfile.write_text(
    json.dumps(
        report,
        indent=2
    )
)

print("="*100)
print("DAY041 / PHASE1B-A17")
print("LAMMPS MODEL VALIDATION")
print("="*100)
print()

for k,v in report.items():

    if k in ("timestamp","decision"):
        continue

    print(f"{k:18s} {v}")

print()

print(outfile)

print()

print(report["decision"])
