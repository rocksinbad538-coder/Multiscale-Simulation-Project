#!/usr/bin/env python3
"""
Extract selected frames from a LAMMPS dump trajectory and write multi-frame XYZ files.

Expected dump columns:
id type mol q x y z

Outputs:
- all atoms
- scaffold only: types 1,2
- water only: types 3,4
- water oxygens only: type 3
"""

from __future__ import annotations

import argparse
from pathlib import Path


TYPE_TO_ELEMENT = {
    1: "B",
    2: "N",
    3: "O",
    4: "H",
}


def read_dump_frames(dump_path: Path):
    frames = []

    with dump_path.open("r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        if lines[i].strip() != "ITEM: TIMESTEP":
            i += 1
            continue

        timestep = int(lines[i + 1].strip())
        i += 2

        if lines[i].strip() != "ITEM: NUMBER OF ATOMS":
            raise ValueError("Unexpected dump format: missing NUMBER OF ATOMS.")
        n_atoms = int(lines[i + 1].strip())
        i += 2

        if not lines[i].startswith("ITEM: BOX BOUNDS"):
            raise ValueError("Unexpected dump format: missing BOX BOUNDS.")
        i += 4

        if not lines[i].startswith("ITEM: ATOMS"):
            raise ValueError("Unexpected dump format: missing ATOMS section.")

        columns = lines[i].strip().split()[2:]
        i += 1

        atoms = []
        for _ in range(n_atoms):
            parts = lines[i].split()
            record = dict(zip(columns, parts))
            atoms.append(record)
            i += 1

        frames.append((timestep, atoms))

    if not frames:
        raise ValueError(f"No frames found in {dump_path}")

    return frames


def select_frames(frames, max_frames: int):
    if len(frames) <= max_frames:
        return frames

    idxs = [
        round(i * (len(frames) - 1) / (max_frames - 1))
        for i in range(max_frames)
    ]
    return [frames[i] for i in idxs]


def filter_atoms(atoms, mode: str):
    if mode == "all":
        return atoms
    if mode == "scaffold":
        return [a for a in atoms if int(a["type"]) in (1, 2)]
    if mode == "water":
        return [a for a in atoms if int(a["type"]) in (3, 4)]
    if mode == "water_oxygen":
        return [a for a in atoms if int(a["type"]) == 3]
    if mode == "all_no_h":
        return [a for a in atoms if int(a["type"]) in (1, 2, 3)]
    raise ValueError(f"Unknown mode: {mode}")


def write_multiframe_xyz(frames, output_path: Path, mode: str):
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w") as f:
        for timestep, atoms in frames:
            selected = filter_atoms(atoms, mode)
            selected = sorted(selected, key=lambda a: int(a["id"]))

            f.write(f"{len(selected)}\n")
            f.write(f"{mode}; timestep {timestep}\n")

            for a in selected:
                atom_type = int(a["type"])
                element = TYPE_TO_ELEMENT.get(atom_type, "X")
                x = float(a["x"])
                y = float(a["y"])
                z = float(a["z"])
                f.write(f"{element:2s} {x:15.6f} {y:15.6f} {z:15.6f}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dump", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--max-frames", type=int, default=50)
    args = parser.parse_args()

    frames = read_dump_frames(args.dump)
    selected_frames = select_frames(frames, args.max_frames)

    modes = ["all", "all_no_h", "scaffold", "water", "water_oxygen"]

    for mode in modes:
        out = args.output_prefix.with_name(args.output_prefix.name + f"_{mode}.xyz")
        write_multiframe_xyz(selected_frames, out, mode)
        print(f"Wrote {out}")

    print(f"Input frames: {len(frames)}")
    print(f"Written frames: {len(selected_frames)}")


if __name__ == "__main__":
    main()
