#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import csv
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[2]

sys.path.insert(
    0,
    str(ROOT / "scripts" / "phase2")
)

from md_analysis.io import iter_lammps_dump

from md_analysis.geometry import (
    centroid,
    radius_of_gyration,
    bounding_box,
)

from md_analysis.alignment import (
    aligned_rmsd,
)


if len(sys.argv) != 2:
    raise SystemExit(
        "Usage: analyze_phase5_streaming.py "
        "<temperature-folder>"
    )


WORK = Path(
    sys.argv[1]
).resolve()

TRAJECTORY = WORK / "production.xyz"
LOG = WORK / "production.log"
OUT = WORK / "analysis"

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


if not TRAJECTORY.exists():
    raise RuntimeError(
        f"Missing trajectory: {TRAJECTORY}"
    )

if not LOG.exists():
    raise RuntimeError(
        f"Missing production log: {LOG}"
    )


# ============================================================
# 1. STREAMING STRUCTURAL ANALYSIS
# ============================================================

trajectory_csv = (
    OUT / "trajectory_summary.csv"
)

rmsd_csv = (
    OUT / "aligned_rmsd.csv"
)

start_time = time.time()

frame_count = 0
first_step = None
last_step = None
reference_atoms = None

rg_sum = 0.0
rmsd_sum = 0.0
rmsd_max = 0.0
final_rg = None
final_rmsd = None


with (
    trajectory_csv.open(
        "w",
        newline="",
    ) as tf,
    rmsd_csv.open(
        "w",
        newline="",
    ) as rf,
):

    trajectory_writer = csv.DictWriter(
        tf,
        fieldnames=[
            "timestep",
            "COMx",
            "COMy",
            "COMz",
            "Rg",
            "Lx",
            "Ly",
            "Lz",
        ],
    )

    rmsd_writer = csv.DictWriter(
        rf,
        fieldnames=[
            "timestep",
            "AlignedRMSD",
        ],
    )

    trajectory_writer.writeheader()
    rmsd_writer.writeheader()

    for frame in iter_lammps_dump(
        TRAJECTORY
    ):

        if reference_atoms is None:

            reference_atoms = [
                dict(a)
                for a in frame["atoms"]
            ]

            first_step = int(
                frame["timestep"]
            )

        com = centroid(
            frame["atoms"]
        )

        rg = radius_of_gyration(
            frame["atoms"]
        )

        box = bounding_box(
            frame["atoms"]
        )

        rmsd = aligned_rmsd(
            reference_atoms,
            frame["atoms"],
        )

        trajectory_writer.writerow({
            "timestep":
                frame["timestep"],
            "COMx":
                com[0],
            "COMy":
                com[1],
            "COMz":
                com[2],
            "Rg":
                rg,
            "Lx":
                box["xmax"]-box["xmin"],
            "Ly":
                box["ymax"]-box["ymin"],
            "Lz":
                box["zmax"]-box["zmin"],
        })

        rmsd_writer.writerow({
            "timestep":
                frame["timestep"],
            "AlignedRMSD":
                rmsd,
        })

        frame_count += 1

        last_step = int(
            frame["timestep"]
        )

        rg_sum += rg
        rmsd_sum += rmsd

        rmsd_max = max(
            rmsd_max,
            rmsd,
        )

        final_rg = rg
        final_rmsd = rmsd

        if frame_count % 100000 == 0:

            elapsed = (
                time.time()-start_time
            )

            print(
                f"frames={frame_count} "
                f"step={last_step} "
                f"elapsed_s={elapsed:.1f}",
                flush=True,
            )


if frame_count == 0:
    raise RuntimeError(
        "No trajectory frames parsed."
    )


# ============================================================
# 2. STREAM THERMODYNAMICS FROM PRODUCTION.LOG
# ============================================================

thermo_csv = (
    OUT / "thermodynamics.csv"
)

thermo_rows = 0
thermo_header = None
active = False


with (
    LOG.open() as source,
    thermo_csv.open(
        "w",
        newline="",
    ) as destination,
):

    writer = None

    for raw in source:

        line = raw.strip()

        if not line:
            continue

        fields = line.split()

        if (
            fields
            and fields[0] == "Step"
            and "PotEng" in fields
        ):

            thermo_header = fields

            writer = csv.DictWriter(
                destination,
                fieldnames=thermo_header,
            )

            writer.writeheader()

            active = True

            continue

        if not active:
            continue

        if line.startswith(
            "Loop time of"
        ):
            break

        fields = line.split()

        if (
            len(fields)
            != len(thermo_header)
        ):
            continue

        try:
            float(fields[0])
        except ValueError:
            continue

        record = dict(
            zip(
                thermo_header,
                fields,
            )
        )

        writer.writerow(record)

        thermo_rows += 1


if thermo_rows == 0:
    raise RuntimeError(
        "No thermodynamic records parsed."
    )


# ============================================================
# 3. REPORT
# ============================================================

elapsed_total = (
    time.time()-start_time
)

report = {
    "trajectory":
        str(TRAJECTORY),
    "frame_count":
        int(frame_count),
    "first_timestep":
        int(first_step),
    "last_timestep":
        int(last_step),
    "thermodynamic_records":
        int(thermo_rows),
    "mean_Rg_A":
        float(rg_sum/frame_count),
    "final_Rg_A":
        float(final_rg),
    "mean_aligned_rmsd_A":
        float(rmsd_sum/frame_count),
    "maximum_aligned_rmsd_A":
        float(rmsd_max),
    "final_aligned_rmsd_A":
        float(final_rmsd),
    "elapsed_seconds":
        float(elapsed_total),
    "memory_model":
        "STREAMING_ONE_FRAME_AT_A_TIME",
}


(
    OUT
    / "PHASE5_STREAMING_ANALYSIS_REPORT.json"
).write_text(
    json.dumps(
        report,
        indent=2,
    )
    + "\n"
)


print()
print("="*90)
print("PHASE5 STREAMING ANALYSIS")
print("="*90)

print(
    f"FRAMES={frame_count}"
)

print(
    f"FIRST_STEP={first_step}"
)

print(
    f"LAST_STEP={last_step}"
)

print(
    f"THERMO_ROWS={thermo_rows}"
)

print(
    f"ELAPSED_SECONDS={elapsed_total:.2f}"
)

print(
    "STATUS=PHASE5_STREAMING_ANALYSIS_PASS"
)
