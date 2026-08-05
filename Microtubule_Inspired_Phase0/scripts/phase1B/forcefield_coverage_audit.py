#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import json
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

SYSTEM = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_parameterized_system"
    / "PHASE1B_MOLECULAR_SYSTEM.json"
)

RUN = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_forcefield_coverage"
)

RUN.mkdir(parents=True, exist_ok=True)


def utc():
    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )


system = json.loads(SYSTEM.read_text())


def coverage(items, field):

    assigned = sum(field in x for x in items)

    return assigned, len(items)


bond_ok,bond_total = coverage(system["bonds"],"parameter_type")
angle_ok,angle_total = coverage(system["angles"],"parameter_type")
improper_ok,improper_total = coverage(system["impropers"],"parameter_type")
dihedral_ok,dihedral_total = coverage(system["dihedrals"],"parameter_type")

overall_ok = bond_ok+angle_ok+improper_ok+dihedral_ok
overall_total = bond_total+angle_total+improper_total+dihedral_total

coverage_percent = 100.0*overall_ok/overall_total

report = {

"timestamp":utc(),

"bond":

{"assigned":bond_ok,"total":bond_total},

"angle":

{"assigned":angle_ok,"total":angle_total},

"improper":

{"assigned":improper_ok,"total":improper_total},

"dihedral":

{"assigned":dihedral_ok,"total":dihedral_total},

"coverage_percent":coverage_percent,

"decision":

"PASS"

if coverage_percent==100.0

else

"FAIL"

}

outfile=RUN/"FORCEFIELD_COVERAGE_REPORT.json"

outfile.write_text(json.dumps(report,indent=2))

print("="*100)
print("DAY041 / PHASE1B-A14")
print("FORCE FIELD COVERAGE")
print("="*100)
print()

print(f"Bonds      {bond_ok}/{bond_total}")
print(f"Angles     {angle_ok}/{angle_total}")
print(f"Impropers  {improper_ok}/{improper_total}")
print(f"Dihedrals  {dihedral_ok}/{dihedral_total}")

print()

print(f"Coverage = {coverage_percent:.1f}%")

print()

print(outfile)

print()

print(report["decision"])
