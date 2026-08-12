#!/usr/bin/env python3

from __future__ import annotations

import json
import csv
import math
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]

MODEL = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_lammps_indexed_model"
    / "PHASE1B_LAMMPS_INDEXED_MODEL.json"
)

DUMP = (
    ROOT
    / "runs"
    / "phase2"
    / "day042_import_test"
    / "minimized.xyz"
)

RUN = (
    ROOT
    / "runs"
    / "phase2"
    / "day043_bond_relaxation"
)

RUN.mkdir(parents=True, exist_ok=True)


def dist(a,b):

    return math.sqrt(
        (a[0]-b[0])**2+
        (a[1]-b[1])**2+
        (a[2]-b[2])**2
    )


model=json.loads(MODEL.read_text())

initial={}
atom_info={}

for atom in model["system"]["atoms"]:

    idx=str(atom["lammps_index"])

    initial[idx]=(

        float(atom["x_A"]),
        float(atom["y_A"]),
        float(atom["z_A"])
    )

    atom_info[idx]=atom

final={}

lines=DUMP.read_text().splitlines()

start=None

for i,line in enumerate(lines):

    if line.startswith("ITEM: ATOMS"):

        start=i+1

        break

for line in lines[start:]:

    s=line.split()

    final[s[0]]=(

        float(s[2]),
        float(s[3]),
        float(s[4])

    )

rows=[]

for bond in model["system"]["bonds"]:

    i,j=bond["indices"]

    i=str(i)
    j=str(j)

    d0=dist(initial[i],initial[j])

    d1=dist(final[i],final[j])

    
    
    delta=d1-d0

    ai=atom_info[i]
    aj=atom_info[j]

    rows.append({

        "atom1":i,
        "atom2":j,

        "atom_id_1":ai["atom_id"],
        "atom_id_2":aj["atom_id"],

        "element1":ai["element"],
        "element2":aj["element"],

        "role1":ai["role"],
        "role2":aj["role"],

        "node1":ai["node_type"],
        "node2":aj["node_type"],

        "parameter":bond["parameter_type"],

        "initial_A":d0,

        "final_A":d1,

        "delta_A":delta,

        "percent_change":100.0*delta/d0

    })


rows.sort(
    key=lambda r: abs(r["delta_A"]),
    reverse=True
)

for rank,row in enumerate(rows,1):

    row["rank"]=rank

csvfile=RUN/"bond_relaxation.csv"

with open(csvfile,"w",newline="") as f:

    writer=csv.DictWriter(f,fieldnames=rows[0].keys())

    writer.writeheader()

    writer.writerows(rows)

deltas=[abs(x["delta_A"]) for x in rows]

report={

    "bond_count":len(rows),

    "mean_delta_A":sum(deltas)/len(deltas),

    "max_delta_A":max(deltas),

    "min_delta_A":min(deltas)

}

jsonfile=RUN/"BOND_RELAXATION_REPORT.json"

jsonfile.write_text(

    json.dumps(report,indent=2)

)

print("="*90)
print("DAY043 / PHASE2-A9")
print("BOND RELAXATION ANALYSIS")
print("="*90)
print()

print("Bonds:",report["bond_count"])

print(f"Mean ΔL : {report['mean_delta_A']:.6f} Å")

print(f"Max  ΔL : {report['max_delta_A']:.6f} Å")

print()

print(csvfile)

print(jsonfile)
