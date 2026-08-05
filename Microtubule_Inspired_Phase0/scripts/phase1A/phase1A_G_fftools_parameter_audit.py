#!/usr/bin/env python3
"""
DAY040 / D040-A7

Authoritative FFTOOLS Force Field Audit

PART 1 / 5
"""

from __future__ import annotations

import csv
import json
import hashlib
import pathlib
import datetime

from dataclasses import dataclass

ROOT = pathlib.Path(__file__).resolve().parents[2]

FF_FILE = (
    ROOT
    / "references"
    / "force_fields"
    / "Rajan_JPCL_2018_hBN_Functionalized"
    / "hBN_functionalized-FFTOOLS.ff"
)

RUN_DIR = (
    ROOT
    / "runs"
    / "phase1A"
    / "day040_phase1A_G_fftools_parameter_audit"
)

RUN_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


def sha256(path):

    h = hashlib.sha256()

    with open(path, "rb") as f:

        while True:

            block = f.read(1024 * 1024)

            if not block:

                break

            h.update(block)

    return h.hexdigest()


def utc():

    return (
        datetime.datetime
        .now(datetime.UTC)
        .isoformat()
        .replace("+00:00", "Z")
    )


@dataclass
class FFEntry:

    section: str

    tokens: list[str]

    raw: str


print("=" * 100)
print("DAY040 / D040-A7")
print("AUTHORITATIVE FFTOOLS PARAMETER AUDIT")
print("=" * 100)
print()

assert FF_FILE.exists()

print("[1] SOURCE FILE")

print(FF_FILE)

print("SHA256 =", sha256(FF_FILE))

print()

entries = []

current_section = None
sections = {

    "ATOMS",

    "BONDS",

    "ANGLES",

    "IMPROPER",

    "DIHEDRALS",

}

for line in FF_FILE.read_text().splitlines():

    line = line.strip()

    if not line:

        continue

    if line.startswith("#"):

        continue

    upper = line.upper()

    if upper in sections:

        current_section = upper

        continue

    if current_section is None:

        continue

    entries.append(

        FFEntry(

            section=current_section,

            tokens=line.split(),

            raw=line,

        )

    )

print("[2] TOTAL ENTRIES")

print(len(entries))

print()
summary = {}

for e in entries:

    summary.setdefault(

        e.section,

        0,

    )

    summary[e.section] += 1

print("[3] SECTION COUNTS")

for k in (

    "ATOMS",

    "BONDS",

    "ANGLES",

    "IMPROPER",

    "DIHEDRALS",

):

    print(

        k,

        summary.get(k,0),

    )

print()
csv_file = RUN_DIR / "FFTOOLS_PARAMETER_INVENTORY.csv"

with open(

    csv_file,

    "w",

    newline="",

) as f:

    w = csv.writer(f)

    w.writerow(

        [

            "section",

            "raw",

        ]

    )

    for e in entries:

        w.writerow(

            [

                e.section,

                e.raw,

            ]

        )

report = {

    "timestamp": utc(),

    "sections": summary,

    "entry_count": len(entries),

    "source_sha256": sha256(FF_FILE),

}

json_file = RUN_DIR / "FFTOOLS_PARAMETER_AUDIT.json"

json_file.write_text(

    json.dumps(

        report,

        indent=2,

    )

)

print("[4] OUTPUTS")

print(csv_file)

print(json_file)

print()

print("DECISION")

print(

    "AUTHORITATIVE_FORCE_FIELD_PARSED"

)
