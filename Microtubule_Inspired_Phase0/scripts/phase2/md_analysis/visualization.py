from __future__ import annotations

import warnings
warnings.filterwarnings(
    "ignore",
    message=".*OVITO.*PyPI"
)

from pathlib import Path

from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer


class MDRenderer:

    def __init__(self):

        self.renderer = TachyonRenderer()

        self.viewport = Viewport(
            type=Viewport.Type.Perspective
        )

    def render_frames(
        self,
        xyz_file,
        output_directory,
        size=(1600,1200),
    ):

        xyz_file = Path(xyz_file)
        output_directory = Path(output_directory)

        output_directory.mkdir(
            parents=True,
            exist_ok=True
        )

        pipeline = import_file(
            str(xyz_file),
            multiple_frames=True,
        )

        pipeline.add_to_scene()

        nframes = pipeline.source.num_frames

        for frame in range(nframes):

            pipeline.compute(frame)

            outfile = (
                output_directory /
                f"{frame:04d}.png"
            )

            self.viewport.render_image(
                filename=str(outfile),
                frame=frame,
                size=size,
                renderer=self.renderer,
            )

        pipeline.remove_from_scene()

        return nframes
