#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import pathlib
import hashlib
import datetime

ROOT = pathlib.Path(__file__).resolve().parents[2]

RUN_DIR = (
    ROOT
    / "runs"
    / "phase1A"
    / "day040_phase1A_semantic_forcefield_audit"
)

RUN_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

SOURCE_CSV = (
    ROOT
    / "runs"
    / "phase1A"
    / "day040_phase1A_G_fftools_parameter_audit"
    / "FFTOOLS_PARAMETER_INVENTORY.csv"
)

assert SOURCE_CSV.exists()


def sha256(path):

    h = hashlib.sha256()

    with open(path,"rb") as f:

        while True:

            b=f.read(1024*1024)

            if not b:
                break

            h.update(b)

    return h.hexdigest()


def utc():

    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00","Z")
    )

rows=[]

with open(SOURCE_CSV) as f:

    reader=csv.DictReader(f)

    for r in reader:

        rows.append(r)

print("="*100)
print("DAY040 / D040-A7B")
print("SEMANTIC FORCE FIELD AUDIT")
print("="*100)
print()
print("[1] INVENTORY")
print("records =",len(rows))
print()
