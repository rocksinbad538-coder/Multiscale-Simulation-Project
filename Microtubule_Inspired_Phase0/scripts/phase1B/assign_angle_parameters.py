#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import json
import collections

HARMONIC_K_KJMOL_TO_LAMMPS_REAL = 1.0 / (2.0 * 4.184)
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

SYSTEM = (
ROOT
/"runs"
/
"phase1B"
/
"day041_parameterized_system"
/
"PHASE1B_MOLECULAR_SYSTEM.json"
)

ANGLES = (
ROOT
/"runs"
/
"phase1B"
/
"day041_angle_parameter_mapping"
/
"PHASE1B_CANONICAL_ANGLES.json"
)

RUN = (
ROOT
/"runs"
/
"phase1B"
/
"day041_parameterized_system"
)

def utc():

    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )

system=json.loads(SYSTEM.read_text())

angle_lookup={}

for a in json.loads(ANGLES.read_text())["angles"]:

    key=(
        a["atom1"],
        a["atom2"],
        a["atom3"]
    )

    angle_lookup[key]=a["forcefield_parameter"]

parameters={

"B-N-B":(120.0,588.70),

"N-B-N":(120.0,230.90),

"B-N-HN":(119.0,197.25),

"N-B-HB":(119.9,250.89)

}

hist=collections.Counter()

new_angles=[]

for a in system["angles"]:

    key=(
        a["atom1"],
        a["atom2"],
        a["atom3"]
    )

    ff=angle_lookup[key]

    theta,k_source=parameters[ff]

    aa=dict(a)

    aa["parameter_type"]=ff
    aa["theta0_deg"]=theta

    # Canonical literature value before engine-specific conversion.
    aa["k_source_kJmol_rad2"]=k_source

    # LAMMPS harmonic uses E = K(theta-theta0)^2 and units real.
    aa["k_lammps_real_kcalmol_rad2"] = (
        k_source * HARMONIC_K_KJMOL_TO_LAMMPS_REAL
    )

    hist[ff]+=1

    new_angles.append(aa)

system["angles"]=new_angles

system["metadata"]["angle_parameter_types"]=dict(hist)

system["status"]="ANGLE_PARAMETERS_ASSIGNED"

outfile=RUN/"PHASE1B_MOLECULAR_SYSTEM.json"

outfile.write_text(

json.dumps(

system,

indent=2

)

)

print("="*100)
print("DAY041 / PHASE1B-A12")
print("ANGLE PARAMETER ASSIGNMENT")
print("="*100)
print()

print("[1] PARAMETER TYPES")

for k in sorted(hist):

    print(f"{k:10s} {hist[k]}")

print()

print("[2] OUTPUT")

print(outfile)

print()

print("[3] DECISION")

print("PHASE1B_ANGLE_PARAMETERS_ASSIGNED")
