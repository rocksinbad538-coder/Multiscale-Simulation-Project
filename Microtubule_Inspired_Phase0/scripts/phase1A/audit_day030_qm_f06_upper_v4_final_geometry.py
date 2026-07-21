#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
import re

from collections import Counter
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

LATEST = (
    ROOT
    / (
        Path(
            (
                ROOT
                / "runs/phase1A/"
                  "day030_qm_f06_upper_v4_executions/"
                  "LATEST_V4_EXECUTION.txt"
            ).read_text(
                encoding="utf-8",
            ).strip()
        )
    )
)

OUT = (
    ROOT
    / "runs/phase1A/"
      "day030_qm_f06_upper_v4_final_geometry"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

V4_OUT = LATEST / "v4.out"
START_XYZ = LATEST / "v4_start.xyz"

REPORT = (
    OUT
    / "QM_F06_UPPER_V4_FINAL_GEOMETRY_AUDIT.json"
)

CSV_GEOM = (
    OUT
    / "QM_F06_UPPER_V4_FINAL_GEOMETRY.csv"
)

FINAL_XYZ = (
    OUT
    / "QM_F06_UPPER_V4_FINAL.xyz"
)

LAST_FRAME = (
    OUT
    / "QM_F06_UPPER_V4_LAST_FRAME.xyz"
)

BOND_CSV = (
    OUT
    / "QM_F06_UPPER_V4_BOND_CHECK.csv"
)

CAP_CSV = (
    OUT
    / "QM_F06_UPPER_V4_CAP_CHECK.csv"
)


def require(path: Path):
    if (not path.exists()) or path.stat().st_size == 0:
        raise RuntimeError(path)


def main():

    require(V4_OUT)
    require(START_XYZ)

    print("=" * 78)
    print("QM_F06 UPPER V4 FINAL GEOMETRY AUDIT")
    print("=" * 78)

    text = V4_OUT.read_text(
        encoding="utf-8",
        errors="ignore",
    )

    optimization_converged = (
        "OPTIMIZATION CONVERGED"
        in text
    )

    scf_converged = (
        "FINAL SINGLE POINT ENERGY"
        in text
    )

    failed = any(
        token in text.upper()
        for token in (
            "FAILED SCF",
            "ERROR TERMINATION",
            "ABORTING",
        )
    )

    maxiter = (
        "MAX ITER"
        in text.upper()
    )

    energy = None

    match = re.findall(
        r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)",
        text,
    )

    if match:
        energy = float(match[-1])

    report = {
        "optimization_converged":
            optimization_converged,
        "scf_converged":
            scf_converged,
        "failed":
            failed,
        "maxiter":
            maxiter,
        "final_energy_Eh":
            energy,
        "authorization": {
            "RESP_authorized": False,
            "force_field_authorized": False,
            "MD_authorized": False,
        },
    }

    REPORT.write_text(
        json.dumps(
            report,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        "Optimization converged:",
        optimization_converged,
    )

    print(
        "SCF converged:",
        scf_converged,
    )

    print(
        "Failed:",
        failed,
    )

    print(
        "MaxIter reached:",
        maxiter,
    )

    print(
        "Final energy:",
        energy,
    )

    if (
        optimization_converged
        and scf_converged
        and not failed
        and not maxiter
    ):
        print()
        print(
            "Decision:",
            "POST_QM_GEOMETRY_VALIDATION_READY"
        )
    else:
        print()
        print(
            "Decision:",
            "WAITING_FOR_ORCA_COMPLETION"
        )


if __name__ == "__main__":
    main()

