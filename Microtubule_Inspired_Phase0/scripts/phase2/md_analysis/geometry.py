#!/usr/bin/env python3

from __future__ import annotations

import math


def centroid(atoms):
    """
    Geometric centroid.
    """

    n = len(atoms)

    return (
        sum(a["x"] for a in atoms) / n,
        sum(a["y"] for a in atoms) / n,
        sum(a["z"] for a in atoms) / n,
    )


def translate_to_centroid(atoms):
    """
    Returns coordinates translated to the centroid.
    """

    cx, cy, cz = centroid(atoms)

    coords = []

    for a in atoms:

        coords.append(
            (
                a["x"] - cx,
                a["y"] - cy,
                a["z"] - cz,
            )
        )

    return coords


def radius_of_gyration(atoms):

    coords = translate_to_centroid(atoms)

    s = 0.0

    for x, y, z in coords:

        s += x*x + y*y + z*z

    return math.sqrt(s / len(coords))


def bounding_box(atoms):

    xs = [a["x"] for a in atoms]
    ys = [a["y"] for a in atoms]
    zs = [a["z"] for a in atoms]

    return {

        "xmin": min(xs),
        "xmax": max(xs),

        "ymin": min(ys),
        "ymax": max(ys),

        "zmin": min(zs),
        "zmax": max(zs),

    }
