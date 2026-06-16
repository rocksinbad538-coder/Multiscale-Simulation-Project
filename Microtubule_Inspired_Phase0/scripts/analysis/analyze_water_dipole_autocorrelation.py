#!/usr/bin/env python3

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass
class Frame:
    timestep: int
    mols: np.ndarray
    types: np.ndarray
    xyz: np.ndarray


def iter_lammpstrj(path: Path) -> Iterator[Frame]:
    with path.open("r", encoding="utf-8") as handle:
        while True:
            line = handle.readline()
            if not line:
                return

            if not line.startswith("ITEM: TIMESTEP"):
                continue

            timestep = int(handle.readline().strip())

            if not handle.readline().startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError("Expected NUMBER OF ATOMS")
            natoms = int(handle.readline().strip())

            if not handle.readline().startswith("ITEM: BOX BOUNDS"):
                raise ValueError("Expected BOX BOUNDS")

            for _ in range(3):
                handle.readline()

            columns = handle.readline().strip().split()[2:]
            required = {"mol", "type", "x", "y", "z"}
            missing = required.difference(columns)

            if missing:
                raise ValueError(
                    f"Missing required trajectory columns: {sorted(missing)}"
                )

            col = {name: columns.index(name) for name in required}

            mols = np.empty(natoms, dtype=np.int64)
            types = np.empty(natoms, dtype=np.int32)
            xyz = np.empty((natoms, 3), dtype=np.float64)

            for index in range(natoms):
                fields = handle.readline().split()
                mols[index] = int(fields[col["mol"]])
                types[index] = int(fields[col["type"]])
                xyz[index, 0] = float(fields[col["x"]])
                xyz[index, 1] = float(fields[col["y"]])
                xyz[index, 2] = float(fields[col["z"]])

            yield Frame(
                timestep=timestep,
                mols=mols,
                types=types,
                xyz=xyz,
            )


def water_dipole_vectors(
    frame: Frame,
    oxygen_type: int,
    hydrogen_type: int,
) -> tuple[np.ndarray, np.ndarray]:
    oxygen_mask = frame.types == oxygen_type
    hydrogen_mask = frame.types == hydrogen_type

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
        bad = unique_h_mols[counts != 2]
        raise ValueError(
            f"Expected exactly two hydrogens per water molecule. "
            f"Bad molecule IDs include: {bad[:10]}"
        )

    if not np.array_equal(oxygen_mols, unique_h_mols):
        raise ValueError(
            "Oxygen and hydrogen molecule-ID sets do not match"
        )

    h1 = hydrogen_xyz[starts]
    h2 = hydrogen_xyz[starts + 1]
    midpoint = 0.5 * (h1 + h2)

    vectors = midpoint - oxygen_xyz
    norms = np.linalg.norm(vectors, axis=1)

    if np.any(norms <= 0.0):
        raise ValueError("Zero-length dipole vector detected")

    unit_vectors = vectors / norms[:, None]

    return oxygen_mols, unit_vectors


def validate_uniform_timesteps(timesteps: np.ndarray) -> int:
    if len(timesteps) < 2:
        raise ValueError("At least two frames are required")

    differences = np.diff(timesteps)
    unique = np.unique(differences)

    if len(unique) != 1:
        raise ValueError(
            f"Nonuniform saved-step spacing detected: {unique.tolist()}"
        )

    return int(unique[0])


def autocorrelations(
    vectors: np.ndarray,
) -> pd.DataFrame:
    n_frames = vectors.shape[0]

    instantaneous_mean = vectors.mean(axis=1, keepdims=True)
    fluctuations = vectors - instantaneous_mean

    rows: list[dict[str, float | int]] = []

    fluctuation_zero = np.mean(
        np.sum(fluctuations * fluctuations, axis=2)
    )

    if fluctuation_zero <= 0.0:
        raise ValueError("Connected-correlation normalization is zero")

    for lag in range(n_frames):
        left = vectors[: n_frames - lag]
        right = vectors[lag:]

        dots = np.einsum("tmi,tmi->tm", left, right)

        c1 = float(dots.mean())
        c2 = float((0.5 * (3.0 * dots**2 - 1.0)).mean())

        fluct_left = fluctuations[: n_frames - lag]
        fluct_right = fluctuations[lag:]

        connected_numerator = np.mean(
            np.einsum("tmi,tmi->tm", fluct_left, fluct_right)
        )
        connected_c1 = float(
            connected_numerator / fluctuation_zero
        )

        rows.append(
            {
                "lag_index": lag,
                "n_time_origins": n_frames - lag,
                "C1_raw": c1,
                "C2_raw": c2,
                "C1_connected": connected_c1,
            }
        )

    return pd.DataFrame(rows)


def positive_window_integral(
    time_ps: np.ndarray,
    values: np.ndarray,
) -> tuple[float, float]:
    if len(values) < 2:
        return 0.0, 0.0

    stop = len(values)

    nonpositive = np.where(values[1:] <= 0.0)[0]
    if len(nonpositive):
        stop = int(nonpositive[0] + 2)

    selected_time = time_ps[:stop]
    selected_values = values[:stop]

    integral = float(np.trapezoid(selected_values, selected_time))
    max_time = float(selected_time[-1])

    return integral, max_time


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--trajectory",
        required=True,
        nargs="+",
        type=Path,
    )
    parser.add_argument("--label", required=True)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/dipole_autocorrelation"),
    )
    parser.add_argument("--oxygen-type", type=int, default=3)
    parser.add_argument("--hydrogen-type", type=int, default=4)
    parser.add_argument(
        "--timestep-fs",
        type=float,
        default=0.05,
    )
    parser.add_argument(
        "--min-step",
        type=int,
        default=None,
        help="Optional inclusive minimum timestep",
    )
    parser.add_argument(
        "--max-step",
        type=int,
        default=None,
        help="Optional inclusive maximum timestep",
    )

    args = parser.parse_args()

    for path in args.trajectory:
        if not path.exists():
            raise FileNotFoundError(path)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    unique_frames: dict[int, Frame] = {}

    for path in args.trajectory:
        print(f"Reading {path}")

        for frame in iter_lammpstrj(path):
            if frame.timestep in unique_frames:
                print(
                    f"Skipping duplicate timestep {frame.timestep} "
                    f"from {path.name}"
                )
                continue

            unique_frames[frame.timestep] = frame

    selected_steps = sorted(unique_frames)

    if args.min_step is not None:
        selected_steps = [
            step for step in selected_steps
            if step >= args.min_step
        ]

    if args.max_step is not None:
        selected_steps = [
            step for step in selected_steps
            if step <= args.max_step
        ]

    if len(selected_steps) < 2:
        raise ValueError(
            "The selected timestep window contains fewer than two frames"
        )

    timesteps = np.asarray(selected_steps, dtype=np.int64)
    saved_step_spacing = validate_uniform_timesteps(timesteps)

    molecule_ids_reference: np.ndarray | None = None
    dipole_frames: list[np.ndarray] = []

    for frame_number, timestep in enumerate(timesteps, start=1):
        frame = unique_frames[int(timestep)]

        molecule_ids, vectors = water_dipole_vectors(
            frame,
            oxygen_type=args.oxygen_type,
            hydrogen_type=args.hydrogen_type,
        )

        if molecule_ids_reference is None:
            molecule_ids_reference = molecule_ids
        elif not np.array_equal(
            molecule_ids_reference,
            molecule_ids,
        ):
            raise ValueError(
                "Water molecule IDs changed between frames"
            )

        dipole_frames.append(vectors)

        print(
            f"[{args.label}] frame {frame_number}/{len(timesteps)}, "
            f"step {timestep}"
        )

    vectors = np.stack(dipole_frames, axis=0)
    correlation = autocorrelations(vectors)

    frame_spacing_fs = (
        saved_step_spacing * args.timestep_fs
    )

    correlation["lag_steps"] = (
        correlation["lag_index"] * saved_step_spacing
    )
    correlation["lag_fs"] = (
        correlation["lag_index"] * frame_spacing_fs
    )
    correlation["lag_ps"] = correlation["lag_fs"] / 1000.0

    correlation = correlation[
        [
            "lag_index",
            "lag_steps",
            "lag_fs",
            "lag_ps",
            "n_time_origins",
            "C1_raw",
            "C2_raw",
            "C1_connected",
        ]
    ]

    integral_c1, integral_window_c1 = positive_window_integral(
        correlation["lag_ps"].to_numpy(),
        correlation["C1_raw"].to_numpy(),
    )

    integral_c2, integral_window_c2 = positive_window_integral(
        correlation["lag_ps"].to_numpy(),
        correlation["C2_raw"].to_numpy(),
    )

    integral_connected, integral_window_connected = (
        positive_window_integral(
            correlation["lag_ps"].to_numpy(),
            correlation["C1_connected"].to_numpy(),
        )
    )

    collective_mean = vectors.mean(axis=1)
    collective_magnitude = np.linalg.norm(
        collective_mean,
        axis=1,
    )

    collective_z = collective_mean[:, 2]

    summary = pd.DataFrame(
        [
            {
                "label": args.label,
                "n_frames": vectors.shape[0],
                "n_water": vectors.shape[1],
                "first_timestep": int(timesteps[0]),
                "last_timestep": int(timesteps[-1]),
                "saved_step_spacing": saved_step_spacing,
                "frame_spacing_fs": frame_spacing_fs,
                "total_sampled_window_fs": (
                    timesteps[-1] - timesteps[0]
                ) * args.timestep_fs,
                "initial_C1_connected": float(
                    correlation.iloc[0]["C1_connected"]
                ),
                "final_C1_raw": float(
                    correlation.iloc[-1]["C1_raw"]
                ),
                "final_C2_raw": float(
                    correlation.iloc[-1]["C2_raw"]
                ),
                "final_C1_connected": float(
                    correlation.iloc[-1]["C1_connected"]
                ),
                "positive_window_integral_C1_ps": integral_c1,
                "positive_window_limit_C1_ps": integral_window_c1,
                "positive_window_integral_C2_ps": integral_c2,
                "positive_window_limit_C2_ps": integral_window_c2,
                "positive_window_integral_connected_C1_ps": (
                    integral_connected
                ),
                "positive_window_limit_connected_C1_ps": (
                    integral_window_connected
                ),
                "mean_collective_dipole_magnitude": float(
                    collective_magnitude.mean()
                ),
                "initial_collective_z": float(
                    collective_z[0]
                ),
                "final_collective_z": float(
                    collective_z[-1]
                ),
            }
        ]
    )

    records_csv = (
        args.output_dir / f"{args.label}_dipole_autocorrelation.csv"
    )
    summary_csv = (
        args.output_dir
        / f"{args.label}_dipole_autocorrelation_summary.csv"
    )

    correlation.to_csv(records_csv, index=False)
    summary.to_csv(summary_csv, index=False)

    fig, axis = plt.subplots(figsize=(9, 6))

    axis.plot(
        correlation["lag_ps"],
        correlation["C1_raw"],
        marker="o",
        label="C1 raw",
    )
    axis.plot(
        correlation["lag_ps"],
        correlation["C2_raw"],
        marker="o",
        label="C2 raw",
    )
    axis.plot(
        correlation["lag_ps"],
        correlation["C1_connected"],
        marker="o",
        label="C1 connected",
    )

    axis.axhline(0.0, linewidth=1)
    axis.set_xlabel("Lag time, ps")
    axis.set_ylabel("Orientational autocorrelation")
    axis.set_title(
        f"Water-dipole orientational autocorrelation\n{args.label}"
    )
    axis.legend()
    fig.tight_layout()

    figure = (
        args.output_dir
        / f"{args.label}_dipole_autocorrelation.png"
    )
    fig.savefig(figure, dpi=200)
    plt.close(fig)

    print("\nDipole autocorrelation analysis complete.")
    print(f"Records: {records_csv}")
    print(f"Summary: {summary_csv}")
    print(f"Figure: {figure}")
    print("\nSummary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
