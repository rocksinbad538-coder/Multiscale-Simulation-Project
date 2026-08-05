#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import csv
import json
import hashlib
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

INPUT = (
ROOT
/"runs"
/
"phase1A"
/
"day040_phase1A_parameter_mapping"
/
"PHASE1B_PARAMETER_MAPPING.csv"
)

RUN = (
ROOT
/"runs"
/
"phase1A"
/
"day040_phase1A_topology_mapping"
)

RUN.mkdir(
parents=True,
exist_ok=True,
)

def utc():

    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )

print("="*100)
print("DAY040 / D040-A9")
print("TOPOLOGY MAPPING SPECIFICATION")
print("="*100)
print()
