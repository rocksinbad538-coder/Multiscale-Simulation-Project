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
"day041_molecular_system"
/
"PHASE1B_MOLECULAR_SYSTEM.json"
)

RUN = (
ROOT
/"runs"
/
"phase1B"
/
"day041_parameterized_system"
)

RUN.mkdir(parents=True,exist_ok=True)

def utc():

    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )

system=json.loads(SYSTEM.read_text())

atom_types={}

for atom in system["atoms"]:

    atom_types[atom["atom_id"]]=atom["forcefield_type"]

bond_ff={

frozenset(("B","N")):
("B-N",1.32,2488.80),

frozenset(("B","HB")):
("B-HB",1.20,2112.85),

frozenset(("N","HN")):
("N-HN",1.02,4315.23)

}

hist=collections.Counter()

for bond in system["bonds"]:

    t1=atom_types[bond["atom1"]]

    t2=atom_types[bond["atom2"]]

    key=frozenset((t1,t2))

    ff,r0,k_source=bond_ff[key]

    bond["parameter_type"]=ff
    bond["r0_A"]=r0

    # Canonical literature value before engine-specific conversion.
    bond["k_source_kJmol_A2"]=k_source

    # LAMMPS harmonic uses E = K(r-r0)^2, i.e. the conventional
    # 1/2 factor is already included in K.  units real => kcal/mol.
    bond["k_lammps_real_kcalmol_A2"] = (
        k_source * HARMONIC_K_KJMOL_TO_LAMMPS_REAL
    )

    hist[ff]+=1

system["metadata"]["bond_parameter_types"]=dict(hist)

system["status"]="BOND_PARAMETERS_ASSIGNED"

outfile=RUN/"PHASE1B_MOLECULAR_SYSTEM.json"

outfile.write_text(

json.dumps(

system,

indent=2

)

)

print("="*100)
print("DAY041 / PHASE1B-A11")
print("BOND PARAMETER ASSIGNMENT")
print("="*100)
print()

print("[1] PARAMETER TYPES")

for k in sorted(hist):

    print(f"{k:8s} {hist[k]}")

print()

print("[2] OUTPUT")

print(outfile)

print()

print("[3] DECISION")

print("PHASE1B_BOND_PARAMETERS_ASSIGNED")
