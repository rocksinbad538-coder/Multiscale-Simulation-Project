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
    / "day043_angle_relaxation"
)

RUN.mkdir(parents=True, exist_ok=True)


def angle(a,b,c):

    ba=[a[i]-b[i] for i in range(3)]
    bc=[c[i]-b[i] for i in range(3)]

    nba=math.sqrt(sum(x*x for x in ba))
    nbc=math.sqrt(sum(x*x for x in bc))

    cosang=sum(ba[i]*bc[i] for i in range(3))/(nba*nbc)

    cosang=max(-1.0,min(1.0,cosang))

    return math.degrees(math.acos(cosang))


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

for ang in model["system"]["angles"]:

    i,j,k=[str(x) for x in ang["indices"]]

    a0=angle(initial[i],initial[j],initial[k])
    a1=angle(final[i],final[j],final[k])

    delta=a1-a0

    rows.append({

        "atom1":i,
        "atom2":j,
        "atom3":k,

        "atom_id_1":atom_info[i]["atom_id"],
        "atom_id_2":atom_info[j]["atom_id"],
        "atom_id_3":atom_info[k]["atom_id"],

        "parameter":ang["parameter_type"],

        "initial_deg":a0,
        "final_deg":a1,

        "delta_deg":delta,

        "percent_change":100.0*delta/a0

    })


rows.sort(
    key=lambda r: abs(r["delta_deg"]),
    reverse=True
)

for rank,row in enumerate(rows,1):

    row["rank"]=rank


csvfile=RUN/"angle_relaxation.csv"

with open(csvfile,"w",newline="") as f:

    writer=csv.DictWriter(
        f,
        fieldnames=rows[0].keys()
    )

    writer.writeheader()

    writer.writerows(rows)


vals=[abs(r["delta_deg"]) for r in rows]

report={

    "angle_count":len(rows),

    "mean_delta_deg":sum(vals)/len(vals),

    "max_delta_deg":max(vals),

    "min_delta_deg":min(vals)

}

jsonfile=RUN/"ANGLE_RELAXATION_REPORT.json"

jsonfile.write_text(
    json.dumps(report,indent=2)
)

print("="*90)
print("DAY043 / PHASE2-A10")
print("ANGLE RELAXATION ANALYSIS")
print("="*90)
print()

print("Angles:",report["angle_count"])
print(f"Mean Δθ : {report['mean_delta_deg']:.6f} deg")
print(f"Max  Δθ : {report['max_delta_deg']:.6f} deg")
print()

print(csvfile)
print(jsonfile)
