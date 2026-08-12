#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]

INITIAL = (
    ROOT
    / "runs"
    / "phase1B"
    / "day041_lammps_indexed_model"
    / "PHASE1B_LAMMPS_INDEXED_MODEL.json"
)

FINAL = (
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
    / "day042_relaxation_analysis"
)

RUN.mkdir(
    parents=True,
    exist_ok=True,
)


def read_initial():

    obj = json.loads(INITIAL.read_text())

    atoms = {}

    for atom in obj["system"]["atoms"]:

        idx = str(atom["lammps_index"])

        atoms[idx] = (
            float(atom["x_A"]),
            float(atom["y_A"]),
            float(atom["z_A"]),
        )

    return atoms


def read_dump():

    atoms = {}

    with open(FINAL) as f:
        lines = f.readlines()

    start = None

    for i, line in enumerate(lines):

        if line.startswith("ITEM: ATOMS"):

            start = i + 1
            break

    if start is None:
        raise RuntimeError("ATOMS section not found.")

    for line in lines[start:]:

        s = line.split()

        atoms[s[0]] = (
            float(s[2]),
            float(s[3]),
            float(s[4]),
        )

    return atoms


def main():

    initial = read_initial()
    final = read_dump()

    rows = []
    displacements = []

    for idx in sorted(initial, key=int):

        xi, yi, zi = initial[idx]
        xf, yf, zf = final[idx]

        d = math.sqrt(
            (xf - xi) ** 2
            + (yf - yi) ** 2
            + (zf - zi) ** 2
        )

        displacements.append(d)

        rows.append(
            {
                "lammps_index": int(idx),
                "displacement_A": d,
            }
        )

    rmsd = math.sqrt(
        sum(d * d for d in displacements)
        / len(displacements)
    )

    mean_disp = (
        sum(displacements)
        / len(displacements)
    )

    max_disp = max(displacements)

    csvfile = RUN / "relaxation.csv"

    with open(csvfile, "w", newline="") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=[
                "lammps_index",
                "displacement_A",
            ],
        )

        writer.writeheader()
        writer.writerows(rows)

    report = {

        "atom_count": len(displacements),

        "RMSD_A": rmsd,

        "mean_displacement_A": mean_disp,

        "maximum_shift_A": max_disp,

    }

    jsonfile = RUN / "RELAXATION_REPORT.json"

    jsonfile.write_text(
        json.dumps(
            report,
            indent=2,
        )
    )

    print("=" * 90)
    print("DAY042 / PHASE2-A8")
    print("STRUCTURAL RELAXATION ANALYSIS")
    print("=" * 90)
    print()

    print("Atoms:", len(displacements))
    print(f"RMSD               {rmsd:.6f} Å")
    print(f"Mean displacement  {mean_disp:.6f} Å")
    print(f"Maximum shift      {max_disp:.6f} Å")
    print()

    print(csvfile)
    print(jsonfile)


if __name__ == "__main__":

    main()
