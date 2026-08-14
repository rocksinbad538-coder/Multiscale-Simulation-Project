#!/usr/bin/env python3

from pathlib import Path

from md_analysis.visualization import MDRenderer

ROOT = Path(__file__).resolve().parents[2]

EVENTS = (
    ROOT /
    "runs/phase2/day045_md_analysis/event_trajectories"
)

OUT = (
    ROOT /
    "runs/phase2/day045_md_analysis/event_renders"
)

OUT.mkdir(exist_ok=True)

renderer = MDRenderer()

for xyz in sorted(EVENTS.glob("*.xyz")):

    name = xyz.stem

    output = OUT / name

    n = renderer.render_frames(
        xyz,
        output,
    )

    print(f"{name:15s} {n:4d} frames")
