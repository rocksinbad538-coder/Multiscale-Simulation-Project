#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd


OUTPUT_DIR = Path("results/representative_frames")
GEOMETRY_DIR = Path("systems/hybrid/bn_like_chromophore_scaffold_water/representative_frames")

OXYGEN_TYPE = 3
HYDROGEN_TYPE = 4
POSITIVE_SITE_TYPE = 5
NEGATIVE_SITE_TYPE = 6
CONTACT_CUTOFF_A = 6.0


@dataclass
class Frame:
    timestep: int
    columns: list[str]
    atom_lines: list[str]
    ids: np.ndarray
    mols: np.ndarray
    types: np.ndarray
    xyz: np.ndarray
    box_header: str
    box_lines: list[str]


@dataclass
class Case:
    label: str
    trajectories: list[Path]
    min_step: int
    max_step: int


CASES = [
    Case(
        label="fieldfree_representative",
        trajectories=[
            Path(
                "systems/hybrid/bn_like_chromophore_scaffold_water/outputs/"
                "nvt300_contained_bn_like_chromophore_12dipoles_carved_"
                "extend_5k_to_20k.lammpstrj"
            ),
        ],
        min_step=6000,
        max_step=20000,
    ),
    Case(
        label="fieldZ_early",
        trajectories=[
            Path(
                "systems/hybrid/bn_like_chromophore_scaffold_water/outputs/"
                "nvt300_fieldZ_contained_bn_like_chromophore_12dipoles_"
                "carved_10k.lammpstrj"
            ),
        ],
        min_step=0,
        max_step=5000,
    ),
    Case(
        label="fieldZ_intermediate",
        trajectories=[
            Path(
                "systems/hybrid/bn_like_chromophore_scaffold_water/outputs/"
                "nvt300_fieldZ_contained_bn_like_chromophore_12dipoles_"
                "carved_10k.lammpstrj"
            ),
            Path(
                "systems/hybrid/bn_like_chromophore_scaffold_water/outputs/"
                "nvt300_fieldZ_contained_bn_like_chromophore_12dipoles_"
                "carved_extend_10k_to_20k.lammpstrj"
            ),
        ],
        min_step=8000,
        max_step=12000,
    ),
    Case(
        label="fieldZ_late",
        trajectories=[
            Path(
                "systems/hybrid/bn_like_chromophore_scaffold_water/outputs/"
                "nvt300_fieldZ_contained_bn_like_chromophore_12dipoles_"
                "carved_10k.lammpstrj"
            ),
            Path(
                "systems/hybrid/bn_like_chromophore_scaffold_water/outputs/"
                "nvt300_fieldZ_contained_bn_like_chromophore_12dipoles_"
                "carved_extend_10k_to_20k.lammpstrj"
            ),
        ],
        min_step=15000,
        max_step=20000,
    ),
    Case(
        label="fieldZ_late_replica01",
        trajectories=[
            Path(
                "systems/hybrid/bn_like_chromophore_scaffold_water/outputs/"
                "nvt300_fieldZ_contained_bn_like_chromophore_12dipoles_"
                "carved_replica01_10k.lammpstrj"
            ),
            Path(
                "systems/hybrid/bn_like_chromophore_scaffold_water/outputs/"
                "nvt300_fieldZ_contained_bn_like_chromophore_12dipoles_"
                "carved_replica01_extend_10k_to_20k.lammpstrj"
            ),
        ],
        min_step=15000,
        max_step=20000,
    ),
]


def iter_lammpstrj(path: Path) -> Iterator[Frame]:
    with path.open("r", encoding="utf-8") as handle:
        while True:
            line = handle.readline()
            if not line:
                return

            if not line.startswith("ITEM: TIMESTEP"):
                continue

            timestep = int(handle.readline().strip())

            number_header = handle.readline()
            if not number_header.startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError("Expected NUMBER OF ATOMS")

            natoms = int(handle.readline().strip())

            box_header = handle.readline().rstrip("\n")
            if not box_header.startswith("ITEM: BOX BOUNDS"):
                raise ValueError("Expected BOX BOUNDS")

            box_lines = [
                handle.readline().rstrip("\n")
                for _ in range(3)
            ]

            atom_header = handle.readline().strip()
            if not atom_header.startswith("ITEM: ATOMS"):
                raise ValueError("Expected ATOMS section")

            columns = atom_header.split()[2:]
            required = {"id", "mol", "type", "x", "y", "z"}
            missing = required.difference(columns)

            if missing:
                raise ValueError(
                    f"{path} is missing columns: {sorted(missing)}"
                )

            col = {name: columns.index(name) for name in required}

            atom_lines: list[str] = []
            ids = np.empty(natoms, dtype=np.int64)
            mols = np.empty(natoms, dtype=np.int64)
            types = np.empty(natoms, dtype=np.int32)
            xyz = np.empty((natoms, 3), dtype=np.float64)

            for index in range(natoms):
                raw = handle.readline().rstrip("\n")
                fields = raw.split()

                atom_lines.append(raw)
                ids[index] = int(fields[col["id"]])
                mols[index] = int(fields[col["mol"]])
                types[index] = int(fields[col["type"]])
                xyz[index, 0] = float(fields[col["x"]])
                xyz[index, 1] = float(fields[col["y"]])
                xyz[index, 2] = float(fields[col["z"]])

            yield Frame(
                timestep=timestep,
                columns=columns,
                atom_lines=atom_lines,
                ids=ids,
                mols=mols,
                types=types,
                xyz=xyz,
                box_header=box_header,
                box_lines=box_lines,
            )


def load_unique_frames(paths: list[Path]) -> dict[int, Frame]:
    frames: dict[int, Frame] = {}

    for path in paths:
        if not path.exists():
            raise FileNotFoundError(path)

        print(f"Reading {path}")

        for frame in iter_lammpstrj(path):
            if frame.timestep in frames:
                print(
                    f"Skipping duplicate timestep {frame.timestep} "
                    f"from {path.name}"
                )
                continue

            frames[frame.timestep] = frame

    return frames


def water_dipoles(frame: Frame) -> tuple[np.ndarray, np.ndarray]:
    oxygen_mask = frame.types == OXYGEN_TYPE
    hydrogen_mask = frame.types == HYDROGEN_TYPE

    oxygen_mols = frame.mols[oxygen_mask]
    oxygen_xyz = frame.xyz[oxygen_mask]

    hydrogen_mols = frame.mols[hydrogen_mask]
    hydrogen_xyz = frame.xyz[hydrogen_mask]

    oxygen_order = np.argsort(oxygen_mols)
    oxygen_mols = oxygen_mols[oxygen_order]
    oxygen_xyz = oxygen_xyz[oxygen_order]

    hydrogen_order = np.argsort(hydrogen_mols)
    hydrogen_mols = hydrogen_mols[hydrogen_order]
    hydrogen_xyz = hydrogen_xyz[hydrogen_order]

    unique_h_mols, starts, counts = np.unique(
        hydrogen_mols,
        return_index=True,
        return_counts=True,
    )

    if not np.all(counts == 2):
        raise ValueError("Expected two hydrogen atoms per water molecule")

    if not np.array_equal(oxygen_mols, unique_h_mols):
        raise ValueError("Oxygen and hydrogen molecule sets do not match")

    h1 = hydrogen_xyz[starts]
    h2 = hydrogen_xyz[starts + 1]

    vectors = 0.5 * (h1 + h2) - oxygen_xyz
    vectors /= np.linalg.norm(vectors, axis=1)[:, None]

    return oxygen_xyz, vectors


def define_pairs(reference: Frame) -> list[tuple[int, int]]:
    positive_ids = reference.ids[
        reference.types == POSITIVE_SITE_TYPE
    ]
    negative_ids = reference.ids[
        reference.types == NEGATIVE_SITE_TYPE
    ]

    positive_xyz = reference.xyz[
        reference.types == POSITIVE_SITE_TYPE
    ]
    negative_xyz = reference.xyz[
        reference.types == NEGATIVE_SITE_TYPE
    ]

    if len(positive_ids) != 12 or len(negative_ids) != 12:
        raise ValueError(
            f"Expected 12+12 sites, found "
            f"{len(positive_ids)}+{len(negative_ids)}"
        )

    distance_matrix = np.linalg.norm(
        positive_xyz[:, None, :] - negative_xyz[None, :, :],
        axis=2,
    )

    unused_positive = set(range(12))
    unused_negative = set(range(12))
    raw_pairs: list[tuple[int, int, np.ndarray]] = []

    while unused_positive:
        best: tuple[float, int, int] | None = None

        for ipos in unused_positive:
            for ineg in unused_negative:
                candidate = (
                    float(distance_matrix[ipos, ineg]),
                    ipos,
                    ineg,
                )

                if best is None or candidate < best:
                    best = candidate

        if best is None:
            raise RuntimeError("Could not pair chromophore sites")

        _, ipos, ineg = best
        center = 0.5 * (
            positive_xyz[ipos] + negative_xyz[ineg]
        )

        raw_pairs.append(
            (
                int(positive_ids[ipos]),
                int(negative_ids[ineg]),
                center,
            )
        )

        unused_positive.remove(ipos)
        unused_negative.remove(ineg)

    raw_pairs.sort(
        key=lambda item: (
            item[2][2],
            np.arctan2(item[2][1], item[2][0]),
        )
    )

    return [
        (positive_id, negative_id)
        for positive_id, negative_id, _ in raw_pairs
    ]


def frame_metrics(
    frame: Frame,
    pairs: list[tuple[int, int]],
) -> dict[str, float | int]:
    oxygen_xyz, dipoles = water_dipoles(frame)

    mean_cos_z = float(dipoles[:, 2].mean())
    mean_S_z = float(
        (0.5 * (3.0 * dipoles[:, 2] ** 2 - 1.0)).mean()
    )

    id_to_xyz = {
        int(atom_id): frame.xyz[index]
        for index, atom_id in enumerate(frame.ids)
    }

    pair_counts: list[int] = []

    for positive_id, negative_id in pairs:
        positive_xyz = id_to_xyz[positive_id]
        negative_xyz = id_to_xyz[negative_id]

        distance_positive = np.linalg.norm(
            oxygen_xyz - positive_xyz[None, :],
            axis=1,
        )
        distance_negative = np.linalg.norm(
            oxygen_xyz - negative_xyz[None, :],
            axis=1,
        )

        minimum_distance = np.minimum(
            distance_positive,
            distance_negative,
        )

        pair_counts.append(
            int((minimum_distance <= CONTACT_CUTOFF_A).sum())
        )

    result: dict[str, float | int] = {
        "timestep": frame.timestep,
        "mean_cos_theta_z": mean_cos_z,
        "mean_S_z": mean_S_z,
        "total_contact_count_6A": int(sum(pair_counts)),
    }

    for pair_index, count in enumerate(pair_counts, start=1):
        result[f"pair_{pair_index:02d}_contact_count"] = count

    return result


def robust_standardize(matrix: np.ndarray) -> np.ndarray:
    median = np.median(matrix, axis=0)
    q25 = np.percentile(matrix, 25, axis=0)
    q75 = np.percentile(matrix, 75, axis=0)
    scale = q75 - q25

    standard_deviation = matrix.std(axis=0)
    scale = np.where(scale > 0.0, scale, standard_deviation)
    scale = np.where(scale > 0.0, scale, 1.0)

    return (matrix - median) / scale


def select_medoid(metrics: pd.DataFrame) -> tuple[int, pd.DataFrame]:
    feature_columns = [
        "mean_cos_theta_z",
        "mean_S_z",
        "total_contact_count_6A",
    ] + [
        f"pair_{index:02d}_contact_count"
        for index in range(1, 13)
    ]

    feature_matrix = metrics[feature_columns].to_numpy(
        dtype=np.float64
    )

    standardized = robust_standardize(feature_matrix)
    target = np.median(standardized, axis=0)

    distance = np.sqrt(
        np.sum((standardized - target[None, :]) ** 2, axis=1)
    )

    ranked = metrics.copy()
    ranked["representative_distance"] = distance
    ranked = ranked.sort_values(
        ["representative_distance", "timestep"],
        ignore_index=True,
    )

    selected_timestep = int(ranked.iloc[0]["timestep"])

    return selected_timestep, ranked


def write_lammpstrj(frame: Frame, output: Path) -> None:
    with output.open("w", encoding="utf-8") as handle:
        handle.write("ITEM: TIMESTEP\n")
        handle.write(f"{frame.timestep}\n")
        handle.write("ITEM: NUMBER OF ATOMS\n")
        handle.write(f"{len(frame.atom_lines)}\n")
        handle.write(frame.box_header + "\n")

        for line in frame.box_lines:
            handle.write(line + "\n")

        handle.write("ITEM: ATOMS " + " ".join(frame.columns) + "\n")

        for line in frame.atom_lines:
            handle.write(line + "\n")


def write_xyz(frame: Frame, output: Path) -> None:
    element_map = {
        1: "B",
        2: "N",
        3: "O",
        4: "H",
        5: "X",
        6: "X",
    }

    order = np.argsort(frame.ids)

    with output.open("w", encoding="utf-8") as handle:
        handle.write(f"{len(frame.ids)}\n")
        handle.write(
            f"Representative hybrid frame; timestep={frame.timestep}\n"
        )

        for index in order:
            atom_type = int(frame.types[index])
            element = element_map.get(atom_type, "X")
            x, y, z = frame.xyz[index]

            handle.write(
                f"{element:2s} "
                f"{x:16.8f} {y:16.8f} {z:16.8f}\n"
            )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    GEOMETRY_DIR.mkdir(parents=True, exist_ok=True)

    selection_rows: list[dict[str, float | int | str]] = []

    for case in CASES:
        print("\n" + "=" * 80)
        print(
            f"Selecting {case.label}: "
            f"steps {case.min_step}–{case.max_step}"
        )

        frames = load_unique_frames(case.trajectories)

        selected_steps = [
            step
            for step in sorted(frames)
            if case.min_step <= step <= case.max_step
        ]

        if len(selected_steps) < 2:
            raise ValueError(
                f"{case.label} contains fewer than two frames"
            )

        pairs = define_pairs(frames[selected_steps[0]])

        metrics_rows = [
            frame_metrics(frames[step], pairs)
            for step in selected_steps
        ]

        metrics = pd.DataFrame(metrics_rows)

        selected_timestep, ranking = select_medoid(metrics)
        selected_frame = frames[selected_timestep]

        metrics_csv = (
            OUTPUT_DIR / f"{case.label}_candidate_metrics.csv"
        )
        ranking_csv = (
            OUTPUT_DIR / f"{case.label}_candidate_ranking.csv"
        )

        metrics.to_csv(metrics_csv, index=False)
        ranking.to_csv(ranking_csv, index=False)

        output_lammpstrj = (
            GEOMETRY_DIR
            / f"{case.label}_step_{selected_timestep}.lammpstrj"
        )
        output_xyz = (
            GEOMETRY_DIR
            / f"{case.label}_step_{selected_timestep}.xyz"
        )

        write_lammpstrj(selected_frame, output_lammpstrj)
        write_xyz(selected_frame, output_xyz)

        selected_metrics = metrics.loc[
            metrics["timestep"].eq(selected_timestep)
        ].iloc[0]

        selection_rows.append(
            {
                "selection_label": case.label,
                "window_min_step": case.min_step,
                "window_max_step": case.max_step,
                "n_candidate_frames": len(selected_steps),
                "selected_timestep": selected_timestep,
                "mean_cos_theta_z": selected_metrics[
                    "mean_cos_theta_z"
                ],
                "mean_S_z": selected_metrics["mean_S_z"],
                "total_contact_count_6A": selected_metrics[
                    "total_contact_count_6A"
                ],
                "representative_distance": ranking.iloc[0][
                    "representative_distance"
                ],
                "lammpstrj": str(output_lammpstrj),
                "xyz": str(output_xyz),
            }
        )

        print("\nTop candidates:")
        print(
            ranking[
                [
                    "timestep",
                    "mean_cos_theta_z",
                    "mean_S_z",
                    "total_contact_count_6A",
                    "representative_distance",
                ]
            ].head(5).to_string(index=False)
        )

        print(f"\nSelected timestep: {selected_timestep}")
        print(f"Wrote {output_lammpstrj}")
        print(f"Wrote {output_xyz}")

    selection = pd.DataFrame(selection_rows)

    output_summary = (
        OUTPUT_DIR / "day009_representative_frame_selection.csv"
    )
    selection.to_csv(output_summary, index=False)

    print("\n" + "=" * 80)
    print("Representative-frame selection complete.")
    print(selection.to_string(index=False))
    print(f"\nWrote {output_summary}")


if __name__ == "__main__":
    main()
