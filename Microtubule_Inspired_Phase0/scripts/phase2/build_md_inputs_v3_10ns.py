#!/usr/bin/env python3

from pathlib import Path
import shutil
import json

ROOT = Path(__file__).resolve().parents[2]

SRC = ROOT/"runs"/"phase2"/"day047_md_protocol_v2"

DST = ROOT/"runs"/"phase2"/"day048_md_protocol_v3_10ns"

if DST.exists():
    shutil.rmtree(DST)

shutil.copytree(SRC, DST)

prod = DST/"in.production"

txt = prod.read_text()

txt = txt.replace(

    "run 4000000",

    "run 40000000"

)

prod.write_text(txt)

protocol = {

    "version":"3.0",

    "production_length_ns":10,

    "timesteps":40000000,

    "timestep_fs":0.25,

    "dump_every_steps":100,

    "analysis_ready":True,

    "derived_from":"day047_md_protocol_v2"

}

(DST/"MD_PROTOCOL_V3_10NS.json").write_text(

    json.dumps(protocol,indent=2)

)

print("="*90)
print("DAY047 / PHASE2-A48")
print("10 ns PROTOCOL READY")
print("="*90)
print(DST)
