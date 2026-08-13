#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import shutil
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]

PROTOCOL = ROOT / "runs" / "phase2" / "day047_md_protocol_v2"
CAMPAIGN = ROOT / "runs" / "phase2" / "campaign"

TEMPERATURES = [150, 200, 250, 300, 350]

CAMPAIGN.mkdir(parents=True, exist_ok=True)

for T in TEMPERATURES:

    folder = CAMPAIGN / f"{T}K"
    folder.mkdir(exist_ok=True)

    shutil.copy2(
        ROOT / "runs" / "phase1B" / "day041_lammps_export" / "data.lammps",
        folder / "data.lammps"
    )

    for infile in (
        "in.minimize",
        "in.heating",
        "in.nvt",
        "in.production",
    ):

        src = PROTOCOL / infile
        dst = folder / infile

        shutil.copy2(src, dst)

        text = dst.read_text()

        if infile == "in.heating":

            text = re.sub(
                r"velocity all create\s+\S+",
                f"velocity all create {T}",
                text,
            )

            text = re.sub(
                r"fix 1 all nvt temp\s+\S+\s+\S+",
                f"fix 1 all nvt temp {T} {T}",
                text,
            )

        elif infile in ("in.nvt", "in.production"):

            text = re.sub(
                r"fix 1 all nvt temp\s+\S+\s+\S+",
                f"fix 1 all nvt temp {T} {T}",
                text,
            )

        
        text = text.replace(
            "read_data ../../phase1B/day041_lammps_export/data.lammps",
            "read_data data.lammps"
        )

        dst.write_text(text)


print("="*90)
print("DAY046 / PHASE2-A24")
print("MD CAMPAIGN CREATED")
print("="*90)

for T in TEMPERATURES:
    print(CAMPAIGN / f"{T}K")
