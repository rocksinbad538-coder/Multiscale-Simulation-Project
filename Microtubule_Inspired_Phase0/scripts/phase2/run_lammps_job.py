#!/usr/bin/env python3

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[2]

LAMMPS = pathlib.Path(
    "/Users/alejandro/projects/lammps/src/lmp_mpi"
)

WORKDIR = (
    ROOT
    / "runs"
    / "phase2"
    / "day044_md_protocol"
)

OUTDIR = (
    ROOT
    / "runs"
    / "phase2"
    / "day045_md_runs"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)


def run_job(script_name):

    t0 = time.time()

    proc = subprocess.run(

        [
            str(LAMMPS),
            "-in",
            script_name
        ],

        cwd=WORKDIR,

        capture_output=True,

        text=True

    )

    elapsed = time.time() - t0

    logfile = OUTDIR / f"{script_name}.log"

    logfile.write_text(
        proc.stdout +
        "\n" +
        proc.stderr
    )

    result = {

        "script": script_name,

        "returncode": proc.returncode,

        "elapsed_seconds": elapsed,

        "pass": proc.returncode == 0

    }

    outfile = OUTDIR / f"{script_name}.json"

    outfile.write_text(

        json.dumps(

            result,

            indent=2

        )

    )

    print("="*80)

    print(script_name)

    print(result["pass"])

    print(f"{elapsed:.3f} s")

    return proc.returncode


if __name__ == "__main__":

    if len(sys.argv) != 2:

        print(
            "usage: run_lammps_job.py input_script"
        )

        sys.exit(1)

    sys.exit(

        run_job(

            sys.argv[1]

        )

    )
