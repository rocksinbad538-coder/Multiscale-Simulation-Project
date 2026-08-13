#!/usr/bin/env python3

from __future__ import annotations


def read_lammps_dump(path):

    frames = []

    with open(path) as f:

        while True:

            line = f.readline()

            if not line:
                break

            if not line.startswith("ITEM: TIMESTEP"):
                continue

            timestep = int(f.readline())

            header = f.readline()

            if not header.startswith("ITEM: NUMBER"):
                raise RuntimeError("Invalid LAMMPS dump.")

            natoms = int(f.readline())

            header = f.readline()

            if not header.startswith("ITEM: BOX"):
                raise RuntimeError("Invalid LAMMPS dump.")

            # Skip three box-bound lines
            f.readline()
            f.readline()
            f.readline()

            header = f.readline()

            if not header.startswith("ITEM: ATOMS"):
                raise RuntimeError("Invalid LAMMPS dump.")

            columns = header.split()[2:]

            col = {name:i for i,name in enumerate(columns)}

            atoms = []

            for _ in range(natoms):

                s = f.readline().split()

                atom = {}

                for name in columns:

                    value = s[col[name]]

                    if name in ("id","mol","type","ix","iy","iz"):
                        atom[name] = int(value)
                    else:
                        atom[name] = float(value)

                atoms.append(atom)

            frames.append({

                "timestep": timestep,

                "natoms": natoms,

                "atoms": atoms,

            })

    return frames
