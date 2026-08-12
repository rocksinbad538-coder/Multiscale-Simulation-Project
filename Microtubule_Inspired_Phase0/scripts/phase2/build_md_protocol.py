#!/usr/bin/env python3

from __future__ import annotations

import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]

RUN = (
    ROOT
    / "runs"
    / "phase2"
    / "day044_md_protocol"
)

RUN.mkdir(
    parents=True,
    exist_ok=True,
)

protocol = {

    "protocol_version": "1.0",

    "forcefield_version": "Phase1B",

    "units": "real",

    "atom_style": "full",

    "pair_style": "lj/cut",

    "pair_cutoff_A": 12.0,

    "bond_style": "harmonic",

    "angle_style": "harmonic",

    "improper_style": "harmonic",

    "dihedral_style": "zero",

    "neighbor_skin_A": 2.0,

    "neighbor_style": "bin",

    "timestep_fs": 0.25,

    "integrator": "velocity-verlet",

    "ensemble": "NVT",

    "thermostat": "Nose-Hoover",

    "temperature_schedule_K": [

        10,

        50,

        100,

        300

    ],

    "thermo_every": 100,

    "dump_every": 100,

    "restart_every": 1000

}

outfile = RUN / "MD_PROTOCOL.json"

outfile.write_text(
    json.dumps(
        protocol,
        indent=2
    )
)

print("=" * 90)
print("DAY044 / PHASE2-A12")
print("MD PROTOCOL")
print("=" * 90)
print()
print(outfile)
print()
print("MD_PROTOCOL_CREATED")
