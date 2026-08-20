#!/usr/bin/env python3

from pathlib import Path
import json
import shutil

ROOT = Path(__file__).resolve().parents[2]

XYZ = ROOT/"runs"/"phase2"/"campaign_phase5_corrected"/"tddft_xyz"

OUT = ROOT/"runs"/"phase2"/"campaign_phase5_corrected"/"orca_campaign"

OUT.mkdir(parents=True,exist_ok=True)

manifest=[]

orca_template = """! CAM-B3LYP def2-SVP TightSCF RIJCOSX D3BJ TDDFT

%pal
 nprocs 8
end

%tddft
 nroots 20
end

* xyz 0 1
{xyz}
*
"""

for xyzfile in sorted(XYZ.glob("*.xyz")):

    name = xyzfile.stem

    work = OUT/name

    work.mkdir(exist_ok=True)

    shutil.copy2(
        xyzfile,
        work/xyzfile.name
    )

    coords = "\n".join(
        xyzfile.read_text().splitlines()[2:]
    )

    inp = orca_template.format(
        xyz=coords
    )

    (work/"job.inp").write_text(inp)

    manifest.append({

        "structure":name,

        "xyz": str((work / xyzfile.name).relative_to(ROOT)),

        "input":str((work/"job.inp").relative_to(ROOT)),

        "status":"ready"

    })

with open(
    OUT/"ORCA_CAMPAIGN_MANIFEST.json",
    "w"
) as f:

    json.dump(
        manifest,
        f,
        indent=2
    )

print("="*90)
print("PHASE4-A04")
print("ORCA CAMPAIGN GENERATED")
print("="*90)
print("Jobs :",len(manifest))
print(OUT)
