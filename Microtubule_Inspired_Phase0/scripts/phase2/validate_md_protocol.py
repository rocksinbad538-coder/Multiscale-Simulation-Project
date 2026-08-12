#!/usr/bin/env python3

from __future__ import annotations

import json
import pathlib
import subprocess

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

OUTPUT = (
    ROOT
    / "runs"
    / "phase2"
    / "day044_md_validation"
)

OUTPUT.mkdir(
    parents=True,
    exist_ok=True,
)

tests = [

    ("Minimize", "in.minimize"),

    ("Heating", "in.heating"),

    ("NVT", "in.nvt"),

    ("Production", "in.production")

]

results = {}

print("=" * 90)
print("DAY044 / PHASE2-A14")
print("MD PROTOCOL VALIDATION")
print("=" * 90)
print()

for name, script in tests:

    print(f"Testing {script} ...")

    proc = subprocess.run(

        [
            str(LAMMPS),
            "-in",
            script
        ],

        cwd=WORKDIR,

        capture_output=True,

        text=True

    )

    ok = proc.returncode == 0

    results[name] = {

        "pass": ok,

        "returncode": proc.returncode

    }

    logfile = OUTPUT / f"{script}.log"

    logfile.write_text(
        proc.stdout + "\n" + proc.stderr
    )

    print("PASS" if ok else "FAIL")

summary = OUTPUT / "MD_PROTOCOL_VALIDATION.json"

summary.write_text(
    json.dumps(
        results,
        indent=2
    )
)

print()
print(summary)
