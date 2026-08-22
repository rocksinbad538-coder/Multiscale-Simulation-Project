#!/usr/bin/env python3

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[2]

DEFAULT_CAMPAIGN = (
    ROOT
    / "runs"
    / "phase2"
    / "campaign"
)

CAMPAIGN = pathlib.Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else DEFAULT_CAMPAIGN

LAMMPS = "/Users/alejandro/projects/lammps/src/lmp_mpi"

STAGES = [
    "in.minimize",
    "in.heating",
    "in.nvt",
    "in.production",
]

# Only physical temperature directories are campaign targets.
# This prevents audit/, analysis/, or other downstream directories
# from being interpreted as MD jobs.
temperatures = sorted(
    d for d in CAMPAIGN.glob("*K")
    if d.is_dir()
)

if len(sys.argv) > 2:
    target = sys.argv[2]
    temperatures = [
        d for d in temperatures
        if d.name == target
    ]



LOGFILE = CAMPAIGN / "campaign.log"

LOG = open(LOGFILE,"w",buffering=1)

def log(msg):

    print(msg)

    LOG.write(msg + "\n")

    LOG.flush()


print("=" * 90)
print("DAY046 / PHASE2-A25")
print("MD CAMPAIGN LAUNCHER")
print("=" * 90)

for folder in temperatures:

    print()
    log(f"=== {folder.name} ===")

    report = {}

    for stage in STAGES:

        start = time.time()

        result = subprocess.run(
            [LAMMPS, "-in", stage],
            cwd=folder
        )

        elapsed = time.time() - start

        report[stage] = {
            "returncode": result.returncode,
            "elapsed_seconds": elapsed,
            "pass": result.returncode == 0
        }

        status = "PASS" if result.returncode == 0 else "FAIL"

        log(
            f"   {stage:15s} {status} ({elapsed:.2f} s)"
        )

        if result.returncode != 0:

            log("Campaign aborted.")

            LOG.close()

            raise SystemExit(result.returncode)

    outfile = folder / "campaign_report.json"

    outfile.write_text(
        json.dumps(report, indent=2)
    )

    log(str(outfile))


LOG.close()
