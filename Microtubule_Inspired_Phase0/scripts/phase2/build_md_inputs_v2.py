#!/usr/bin/env python3

"""
Builds the definitive Phase-2 MD production protocol.

This script upgrades the production input generated previously
without modifying the original build_md_inputs.py.
"""

from pathlib import Path
import shutil
import json
import re

ROOT = Path(__file__).resolve().parents[2]

SRC = ROOT/"runs"/"phase2"/"day044_md_protocol"
DST = ROOT/"runs"/"phase2"/"day047_md_protocol_v2"

DST.mkdir(parents=True, exist_ok=True)

for f in SRC.iterdir():
    if f.is_file():
        shutil.copy2(f, DST/f.name)

prod = DST/"in.production"

txt = prod.read_text()

# ------------------------------------------------------------------
# Better thermo output
# ------------------------------------------------------------------

txt = re.sub(
    r"thermo_style.*",
    "thermo_style custom step temp pe ke etotal ebond eangle edihed eimp evdwl press",
    txt,
)

# ------------------------------------------------------------------
# Extended dump
# ------------------------------------------------------------------

txt = re.sub(
    r"dump 1 all custom .*",
    (
        "dump 1 all custom 100 "
        "production.xyz "
        "id mol type q x y z ix iy iz vx vy vz"
    ),
    txt,
)

# ------------------------------------------------------------------
# Stable ordering
# ------------------------------------------------------------------

if "dump_modify 1 sort id" not in txt:
    txt += "\n\ndump_modify 1 sort id\n"

prod.write_text(txt)

protocol = {

    "version":"2.0",

    "coordinates":True,

    "charges":True,

    "velocities":True,

    "forces":False,

    "stress":False,

    "per_atom_energy":False,

    "future_ready_for_dipole":True

}

(DST/"MD_PROTOCOL_V2.json").write_text(
    json.dumps(protocol,indent=2)
)

print("="*90)
print("DAY046 / PHASE2-A33")
print("MD PROTOCOL V2 GENERATED")
print("="*90)
print(DST)
