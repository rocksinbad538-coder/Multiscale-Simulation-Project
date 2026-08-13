#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[3]


def get_paths():

    if len(sys.argv) > 1:

        work = pathlib.Path(sys.argv[1]).resolve()

        return {

            "WORK": work,

            "XYZ": work / "production.xyz",

            "LOG": work / "production.log",

            "OUT": work

        }

    work = (
        ROOT
        / "runs"
        / "phase2"
        / "day044_md_protocol"
    )

    out = (
        ROOT
        / "runs"
        / "phase2"
        / "day045_md_analysis"
    )

    return {

        "WORK": work,

        "XYZ": work / "production.xyz",

        "LOG": work / "in.production.log",

        "OUT": out

    }
