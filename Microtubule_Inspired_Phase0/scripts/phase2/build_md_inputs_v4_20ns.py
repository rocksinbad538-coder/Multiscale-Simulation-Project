#!/usr/bin/env python3

from pathlib import Path
import shutil
import json

ROOT = Path(__file__).resolve().parents[2]

# Canonical Phase5-validated MD protocol.
SRC = ROOT / "runs" / "phase2" / "day044_md_protocol"

DST = ROOT / "runs" / "phase2" / "day049_md_protocol_v4_20ns"

if DST.exists():
    shutil.rmtree(DST)

shutil.copytree(SRC, DST)

prod = DST / "in.production"
txt = prod.read_text()

# Canonical protocol is 1 ns at dt = 0.25 fs:
# 4,000,000 steps.
# Extend production to 20 ns:
# 80,000,000 steps.
if "run 4000000" not in txt:
    raise RuntimeError(
        "Expected canonical 1 ns production command not found."
    )

txt = txt.replace(
    "run 4000000",
    "run 80000000",
    1,
)

# Original restart interval:
# 1000 steps * 0.25 fs = 250 fs.
#
# For the 20 ns campaign use:
# 2,000,000 steps * 0.25 fs = 0.5 ns.
if "restart 1000 production.restart" not in txt:
    raise RuntimeError(
        "Expected canonical restart command not found."
    )

txt = txt.replace(
    "restart 1000 production.restart",
    "restart 2000000 production.restart",
    1,
)

prod.write_text(txt)

protocol = {
    "version": "4.1_PHASE5_CORRECTED",
    "production_length_ns": 20.0,
    "timesteps": 80000000,
    "timestep_fs": 0.25,
    "dump_every_steps": 100,
    "restart_every_steps": 2000000,
    "restart_interval_ns": 0.5,
    "analysis_ready": True,
    "derived_from": "day044_md_protocol_PHASE5_validated",
    "forcefield_checkpoint": "PHASE5_MD_CORRECTION_VALIDATED",
    "electrostatics": {
        "pair_style": "lj/cut/coul/long",
        "pair_cutoff_A": 12.0,
        "kspace_style": "pppm",
        "kspace_accuracy": 1.0e-5,
        "pair_mixing_rule": "geometric",
        "special_bonds_lj_coul": [0.0, 0.0, 0.0],
    },
}

(DST / "MD_PROTOCOL_V4_20NS.json").write_text(
    json.dumps(protocol, indent=2) + "\n"
)

print("=" * 90)
print("PHASE5-D08")
print("CORRECTED 20 ns PROTOCOL READY")
print("=" * 90)
print(DST)
