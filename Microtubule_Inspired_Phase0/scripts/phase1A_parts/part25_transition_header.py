#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import csv
import json
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

INPUT = (
ROOT
/"runs"
/
"phase1A"
/
"day040_phase1A_RESP_hydrogen_origin_audit"
/
"RESP_HYDROGEN_ORIGIN_TABLE.csv"
)

RUN = (
ROOT
/"runs"
/
"phase1A"
/
"day040_phase1A_transition_specification"
)

RUN.mkdir(parents=True, exist_ok=True)

def utc():
    return (
        datetime.datetime.now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )

print("="*100)
print("DAY040 / D040-A11")
print("RESP TO FORCE-FIELD TRANSITION SPECIFICATION")
print("="*100)
print()
