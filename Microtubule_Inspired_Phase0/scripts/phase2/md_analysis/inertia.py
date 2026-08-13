#!/usr/bin/env python3

from __future__ import annotations

import numpy as np

from .geometry import centroid


def inertia_tensor(atoms):

    cx, cy, cz = centroid(atoms)

    I = np.zeros((3,3))

    for a in atoms:

        x = a["x"] - cx
        y = a["y"] - cy
        z = a["z"] - cz

        I[0,0] += y*y + z*z
        I[1,1] += x*x + z*z
        I[2,2] += x*x + y*y

        I[0,1] -= x*y
        I[1,0] -= x*y

        I[0,2] -= x*z
        I[2,0] -= x*z

        I[1,2] -= y*z
        I[2,1] -= y*z

    return I


def principal_moments(atoms):

    I = inertia_tensor(atoms)

    eigvals, eigvecs = np.linalg.eigh(I)

    order = np.argsort(eigvals)

    eigvals = eigvals[order]
    eigvecs = eigvecs[:,order]

    return eigvals, eigvecs


def asphericity(l):

    l1,l2,l3 = l

    return l3 - 0.5*(l1+l2)


def acylindricity(l):

    l1,l2,l3 = l

    return l2-l1


def relative_shape_anisotropy(l):

    l1,l2,l3 = l

    s = l1+l2+l3

    if s==0:
        return 0.0

    num = (
        (l1-l2)**2 +
        (l2-l3)**2 +
        (l3-l1)**2
    )

    return 1.5*num/(s*s)
