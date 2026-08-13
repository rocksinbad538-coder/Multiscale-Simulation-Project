#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import pathlib
import csv

from md_analysis.paths import get_paths
from md_analysis.io import read_lammps_dump

from md_analysis.geometry import (
    centroid,
    radius_of_gyration,
    bounding_box,
)

ROOT = pathlib.Path(__file__).resolve().parents[2]

PATHS = get_paths()

OUT = PATHS["OUT"]

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

def _legacy_read_lammps_dump():

    with open(PATHS["XYZ"]) as f:
        lines = [line.rstrip() for line in f]

    frames = []

    i = 0

    while i < len(lines):

        if lines[i] != "ITEM: TIMESTEP":
            raise RuntimeError(
                f"Unexpected format at line {i+1}"
            )

        timestep = int(lines[i+1])

        natoms = int(lines[i+3])

        atom_start = i + 9

        atoms = []

        for line in lines[atom_start:atom_start+natoms]:

            s = line.split()

            atoms.append({

                "id": int(s[0]),
                "type": int(s[1]),
                "x": float(s[2]),
                "y": float(s[3]),
                "z": float(s[4])

            })

        frames.append({

            "timestep": timestep,
            "natoms": natoms,
            "atoms": atoms

        })

        i = atom_start + natoms

    return frames


def _legacy_center_of_mass(frame):

    n = len(frame["atoms"])

    cx = sum(a["x"] for a in frame["atoms"]) / n
    cy = sum(a["y"] for a in frame["atoms"]) / n
    cz = sum(a["z"] for a in frame["atoms"]) / n

    return cx, cy, cz


def _legacy_radius_of_gyration(frame):

    cx, cy, cz = center_of_mass(frame)

    s = 0.0

    for a in frame["atoms"]:

        dx = a["x"] - cx
        dy = a["y"] - cy
        dz = a["z"] - cz

        s += dx*dx + dy*dy + dz*dz

    return (s / len(frame["atoms"])) ** 0.5




def rmsd(frame, reference):

    s = 0.0

    for a,b in zip(frame["atoms"], reference["atoms"]):

        dx = a["x"]-b["x"]
        dy = a["y"]-b["y"]
        dz = a["z"]-b["z"]

        s += dx*dx + dy*dy + dz*dz

    return (s/len(frame["atoms"]))**0.5


def max_displacement(frame, reference):

    m = 0.0

    for a,b in zip(frame["atoms"], reference["atoms"]):

        dx = a["x"]-b["x"]
        dy = a["y"]-b["y"]
        dz = a["z"]-b["z"]

        d = (dx*dx+dy*dy+dz*dz)**0.5

        if d>m:
            m=d

    return m


def _legacy_bounding_box(frame):

    xs = [a["x"] for a in frame["atoms"]]
    ys = [a["y"] for a in frame["atoms"]]
    zs = [a["z"] for a in frame["atoms"]]

    return {

        "xmin": min(xs),
        "xmax": max(xs),

        "ymin": min(ys),
        "ymax": max(ys),

        "zmin": min(zs),
        "zmax": max(zs)

    }


frames = read_lammps_dump(PATHS["XYZ"])

print("=" * 90)
print("DAY045 / PHASE2-A17")
print("TRAJECTORY PARSER")
print("=" * 90)
print()

print("Frames :", len(frames))
print("Atoms  :", frames[0]["natoms"])
print("First timestep :", frames[0]["timestep"])
print("Last timestep  :", frames[-1]["timestep"])


print()

print("Calculating geometric descriptors...")

reference = frames[0]

summary = []

for frame in frames:

    com = centroid(frame['atoms'])

    rg = radius_of_gyration(frame['atoms'])

    box = bounding_box(frame['atoms'])

    summary.append({

        "timestep": frame["timestep"],

        "COMx": com[0],
        "COMy": com[1],
        "COMz": com[2],

        "Rg": rg,

        "RMSD": rmsd(frame,reference),

        "MaxDisplacement": max_displacement(frame,reference),

        "Lx": box["xmax"]-box["xmin"],
        "Ly": box["ymax"]-box["ymin"],
        "Lz": box["zmax"]-box["zmin"]

    })

print("Frames analyzed :", len(summary))
print("Final Rg        :", f"{summary[-1]['Rg']:.6f} Å")


csvfile = OUT / "trajectory_summary.csv"

with open(csvfile,"w",newline="") as f:

    writer = csv.DictWriter(

        f,

        fieldnames=summary[0].keys()

    )

    writer.writeheader()

    writer.writerows(summary)

print(csvfile)


report={

    "frame_count":len(summary),

    "final_Rg_A":summary[-1]["Rg"],

    "maximum_RMSD_A":max(x["RMSD"] for x in summary),

    "mean_RMSD_A":sum(x["RMSD"] for x in summary)/len(summary),

    "maximum_atomic_displacement_A":max(
        x["MaxDisplacement"] for x in summary
    )

}

outfile=OUT/"TRAJECTORY_REPORT.json"

outfile.write_text(
    json.dumps(report,indent=2)
)

print(outfile)



