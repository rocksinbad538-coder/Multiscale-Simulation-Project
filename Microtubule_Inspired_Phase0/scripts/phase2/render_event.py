#!/usr/bin/env python3

from pathlib import Path

from md_analysis.visualization import MDRenderer

ROOT = Path(__file__).resolve().parents[2]

xyz = (
    ROOT /
    "runs/phase2/day045_md_analysis/event_trajectories/event_104k.xyz"
)

out = (
    ROOT /
    "runs/phase2/day045_md_analysis/render_event104"
)

renderer = MDRenderer()

n = renderer.render_frames(
    xyz,
    out,
)

print("="*90)
print("DAY047 / PHASE2-A42")
print("EVENT RENDER")
print("="*90)
print(f"Frames rendered : {n}")
print(out)
