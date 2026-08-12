#!/usr/bin/env python3

import subprocess
import sys

scripts = [
    "in.minimize",
    "in.heating",
    "in.nvt",
    "in.production",
]

for s in scripts:

    print("="*80)
    print("RUNNING", s)
    print("="*80)

    rc = subprocess.run(
        [
            sys.executable,
            "scripts/phase2/run_lammps_job.py",
            s
        ]
    )

    if rc.returncode != 0:

        raise SystemExit(f"\nFAILED: {s}")

print()
print("="*80)
print("MOLECULAR DYNAMICS PIPELINE COMPLETED")
print("="*80)
