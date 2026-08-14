#!/usr/bin/env python3

from pathlib import Path
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]

renders = (
    ROOT /
    "runs/phase2/day045_md_analysis/event_renders"
)

for folder in sorted(renders.iterdir()):

    if not folder.is_dir():
        continue

    frames = sorted(folder.glob("*.png"))

    if not frames:
        continue

    images = [Image.open(f) for f in frames]

    outfile = folder / "animation.gif"

    images[0].save(
        outfile,
        save_all=True,
        append_images=images[1:],
        duration=200,
        loop=0,
    )

    print(outfile)
