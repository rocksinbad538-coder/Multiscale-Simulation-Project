#!/usr/bin/env python3

from __future__ import annotations

import math
import numpy as np

from .geometry import centroid


def centered_coordinates(atoms):

    cx, cy, cz = centroid(atoms)

    xyz = np.array(
        [
            [a["x"], a["y"], a["z"]]
            for a in atoms
        ],
        dtype=float,
    )

    xyz[:,0] -= cx
    xyz[:,1] -= cy
    xyz[:,2] -= cz

    return xyz


def covariance_matrix(reference, target):

    P = centered_coordinates(reference)
    Q = centered_coordinates(target)

    return P.T @ Q


def kabsch_rotation(reference, target):

    H = covariance_matrix(reference, target)

    U, S, VT = np.linalg.svd(H)

    R = VT.T @ U.T

    if np.linalg.det(R) < 0:

        VT[-1,:] *= -1

        R = VT.T @ U.T

    return R


def aligned_rmsd(reference, target):

    P = centered_coordinates(reference)

    Q = centered_coordinates(target)

    R = kabsch_rotation(reference, target)

    Q = Q @ R.T

    diff = P - Q

    return math.sqrt(
        np.mean(
            np.sum(diff*diff, axis=1)
        )
    )
