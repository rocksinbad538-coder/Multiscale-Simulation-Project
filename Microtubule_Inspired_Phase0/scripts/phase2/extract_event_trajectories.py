#!/usr/bin/env python3

from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[2]

traj = ROOT/"runs"/"phase2"/"campaign"/"300K"/"production.xyz"

analysis = ROOT/"runs"/"phase2"/"day045_md_analysis"

frames_to_keep = set(
    json.loads(
        (analysis/"event_windows.json").read_text()
    )
)

outdir = analysis/"event_trajectories"
outdir.mkdir(exist_ok=True)

clusters = {

    "104k": [],
    "583k": [],
    "645k": [],
    "3123k": [],
    "3195k": [],
    "3410k": [],

}

with traj.open() as f:

    while True:

        line = f.readline()

        if not line:
            break

        if not line.startswith("ITEM: TIMESTEP"):
            continue

        block = [line]

        timestep = int(f.readline())
        block.append(str(timestep)+"\n")

        for _ in range(7):
            block.append(f.readline())

        natoms = int(block[3])

        for _ in range(natoms):
            block.append(f.readline())

        if timestep in frames_to_keep:

            if timestep < 200000:
                clusters["104k"].append(block)

            elif timestep < 640000:
                clusters["583k"].append(block)

            elif timestep < 700000:
                clusters["645k"].append(block)

            elif timestep < 3130000:
                clusters["3123k"].append(block)

            elif timestep < 3200000:
                clusters["3195k"].append(block)

            else:
                clusters["3410k"].append(block)

for name, frames in clusters.items():

    outfile = outdir/f"event_{name}.xyz"

    with outfile.open("w") as out:

        for frame in frames:

            out.writelines(frame)

print("="*90)
print("DAY047 / PHASE2-A41")
print("EVENT TRAJECTORIES")
print("="*90)
print(outdir)
