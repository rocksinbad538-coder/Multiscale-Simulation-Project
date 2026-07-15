#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("runs/phase1A/day016_md_bath_extraction")

DEFAULT_INPUT = (
    ROOT
    / "state_identity_analysis"
    / "low_state_identity_tracking.csv"
)

TRAJ_DIR = ROOT / "tracked_site_energy_trajectory"
HAM_DIR = ROOT / "tracked_hamiltonian_diagonals"

CHROMOPHORES = ["PYR2", "PYR3", "PYR4", "PYR5"]
FRAME_SPACING_PS = 5.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build site-energy trajectories and diagonal Hamiltonians "
            "using electronic-state identity rather than adiabatic root number."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="State-identity tracking CSV.",
    )
    parser.add_argument(
        "--minimum-character-weight",
        type=float,
        default=0.70,
        help=(
            "Minimum accepted weight of the target HOMO-LUMO "
            "configuration."
        ),
    )
    return parser.parse_args()


def clean_frame_products(directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)

    for pattern in [
        "Htracked_frame*.csv",
        "Htracked_frame*.npy",
        "Htracked_centered_frame*.csv",
        "Htracked_centered_frame*.npy",
    ]:
        for path in directory.glob(pattern):
            if path.is_file():
                path.unlink()


def validate(df: pd.DataFrame, minimum_weight: float) -> pd.DataFrame:
    required = [
        "frame",
        "cluster",
        "tracked_root",
        "tracked_energy_eV",
        "tracked_fosc",
        "tracked_HOMO_LUMO_weight",
        "alternate_root",
        "alternate_energy_eV",
        "state_separation_meV",
    ]

    missing = [column for column in required if column not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns: {missing}")

    work = df.copy()

    work["frame"] = pd.to_numeric(
        work["frame"], errors="raise"
    ).astype(int)

    work["tracked_root"] = pd.to_numeric(
        work["tracked_root"], errors="raise"
    ).astype(int)

    numeric = [
        "tracked_energy_eV",
        "tracked_fosc",
        "tracked_HOMO_LUMO_weight",
        "alternate_energy_eV",
        "state_separation_meV",
    ]

    for column in numeric:
        work[column] = pd.to_numeric(
            work[column], errors="coerce"
        )

    if work[numeric].isna().any().any():
        raise SystemExit("Missing numeric state-tracking values.")

    duplicated = work.duplicated(
        ["frame", "cluster"],
        keep=False,
    )

    if duplicated.any():
        raise SystemExit(
            "Duplicate frame/cluster rows:\n"
            + work.loc[
                duplicated,
                ["frame", "cluster"],
            ].to_string(index=False)
        )

    unexpected = sorted(
        set(work["cluster"]) - set(CHROMOPHORES)
    )
    if unexpected:
        raise SystemExit(f"Unexpected clusters: {unexpected}")

    expected = set(CHROMOPHORES)
    incomplete = []

    for frame, group in work.groupby("frame"):
        observed = set(group["cluster"])

        if observed != expected or len(group) != 4:
            incomplete.append(
                {
                    "frame": int(frame),
                    "observed": sorted(observed),
                    "missing": sorted(expected - observed),
                }
            )

    if incomplete:
        raise SystemExit(
            "Incomplete frames:\n"
            + json.dumps(incomplete, indent=2)
        )

    weak = work.loc[
        work["tracked_HOMO_LUMO_weight"] < minimum_weight,
        [
            "frame",
            "cluster",
            "tracked_root",
            "tracked_HOMO_LUMO_weight",
        ],
    ]

    if len(weak):
        raise SystemExit(
            "Tracked-state character below threshold:\n"
            + weak.to_string(index=False)
        )

    invalid_roots = work.loc[
        ~work["tracked_root"].isin([1, 2]),
        ["frame", "cluster", "tracked_root"],
    ]

    if len(invalid_roots):
        raise SystemExit(
            "Unexpected tracked roots:\n"
            + invalid_roots.to_string(index=False)
        )

    return work.sort_values(
        ["frame", "cluster"]
    ).reset_index(drop=True)


def pivot_quantity(
    df: pd.DataFrame,
    quantity: str,
) -> pd.DataFrame:
    table = (
        df.pivot(
            index="frame",
            columns="cluster",
            values=quantity,
        )
        .reset_index()
    )

    table = table[
        ["frame"] + CHROMOPHORES
    ].sort_values("frame")

    if table[CHROMOPHORES].isna().any().any():
        raise SystemExit(
            f"NaN detected after pivoting {quantity}."
        )

    return table


def main() -> None:
    args = parse_args()

    if not args.input.is_file():
        raise SystemExit(f"Input not found: {args.input}")

    raw = pd.read_csv(args.input)
    df = validate(
        raw,
        minimum_weight=args.minimum_character_weight,
    )

    TRAJ_DIR.mkdir(parents=True, exist_ok=True)
    clean_frame_products(HAM_DIR)

    energy = pivot_quantity(df, "tracked_energy_eV")
    roots = pivot_quantity(df, "tracked_root")
    oscillator = pivot_quantity(df, "tracked_fosc")
    character = pivot_quantity(
        df,
        "tracked_HOMO_LUMO_weight",
    )
    separation = pivot_quantity(
        df,
        "state_separation_meV",
    )

    centered = energy.copy()
    centered[CHROMOPHORES] = energy[
        CHROMOPHORES
    ].sub(
        energy[CHROMOPHORES].mean(axis=1),
        axis=0,
    )

    energy.to_csv(
        TRAJ_DIR / "tracked_site_energy_trajectory.csv",
        index=False,
    )
    centered.to_csv(
        TRAJ_DIR
        / "tracked_site_energy_trajectory_centered.csv",
        index=False,
    )
    roots.to_csv(
        TRAJ_DIR / "tracked_root_trajectory.csv",
        index=False,
    )
    oscillator.to_csv(
        TRAJ_DIR
        / "tracked_oscillator_strength_trajectory.csv",
        index=False,
    )
    character.to_csv(
        TRAJ_DIR
        / "tracked_character_weight_trajectory.csv",
        index=False,
    )
    separation.to_csv(
        TRAJ_DIR
        / "low_state_separation_trajectory_meV.csv",
        index=False,
    )

    np.save(
        TRAJ_DIR / "tracked_site_energy_trajectory_eV.npy",
        energy[CHROMOPHORES].to_numpy(dtype=float),
    )
    np.save(
        TRAJ_DIR
        / "tracked_site_energy_trajectory_centered_eV.npy",
        centered[CHROMOPHORES].to_numpy(dtype=float),
    )

    summary_rows = []

    for _, row in energy.iterrows():
        frame = int(row["frame"])
        diagonal = row[CHROMOPHORES].to_numpy(dtype=float)

        matrix = np.diag(diagonal)

        centered_row = centered.loc[
            centered["frame"] == frame,
            CHROMOPHORES,
        ].iloc[0].to_numpy(dtype=float)

        centered_matrix = np.diag(centered_row)

        np.save(
            HAM_DIR / f"Htracked_frame{frame:03d}.npy",
            matrix,
        )
        pd.DataFrame(
            matrix,
            index=CHROMOPHORES,
            columns=CHROMOPHORES,
        ).to_csv(
            HAM_DIR / f"Htracked_frame{frame:03d}.csv"
        )

        np.save(
            HAM_DIR
            / f"Htracked_centered_frame{frame:03d}.npy",
            centered_matrix,
        )
        pd.DataFrame(
            centered_matrix,
            index=CHROMOPHORES,
            columns=CHROMOPHORES,
        ).to_csv(
            HAM_DIR
            / f"Htracked_centered_frame{frame:03d}.csv"
        )

        reference_mean = float(
            np.mean(diagonal[:3])
        )

        summary_rows.append(
            {
                "frame": frame,
                "time_ps": frame * FRAME_SPACING_PS,
                "trace_eV": float(np.trace(matrix)),
                "mean_diag_eV": float(diagonal.mean()),
                "diagonal_spread_eV": float(
                    diagonal.max() - diagonal.min()
                ),
                "PYR5_offset_vs_PYR2_PYR4_eV": float(
                    reference_mean - diagonal[3]
                ),
                "PYR5_offset_vs_PYR2_PYR4_meV": float(
                    1000.0
                    * (reference_mean - diagonal[3])
                ),
            }
        )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(
        HAM_DIR / "tracked_hamiltonian_summary.csv",
        index=False,
    )

    root_counts = (
        df.groupby(["cluster", "tracked_root"])
        .size()
        .rename("n_jobs")
        .reset_index()
    )

    cluster_stats = (
        df.groupby("cluster")
        .agg(
            n_frames=("frame", "count"),
            mean_energy_eV=("tracked_energy_eV", "mean"),
            std_energy_eV=("tracked_energy_eV", "std"),
            min_energy_eV=("tracked_energy_eV", "min"),
            max_energy_eV=("tracked_energy_eV", "max"),
            mean_fosc=("tracked_fosc", "mean"),
            mean_character_weight=(
                "tracked_HOMO_LUMO_weight",
                "mean",
            ),
            minimum_character_weight=(
                "tracked_HOMO_LUMO_weight",
                "min",
            ),
            mean_state_separation_meV=(
                "state_separation_meV",
                "mean",
            ),
        )
        .reset_index()
    )

    cluster_stats.to_csv(
        TRAJ_DIR / "tracked_state_statistics.csv",
        index=False,
    )

    mean_offset_mev = float(
        summary["PYR5_offset_vs_PYR2_PYR4_meV"].mean()
    )

    metadata = {
        "production_day": "018",
        "description": (
            "Site-energy trajectory obtained by tracking the "
            "low-lying state with maximum 52a->53a configuration weight."
        ),
        "state_tracking_rule": (
            "Among S1 and S2, select the root with the largest "
            "HOMO-LUMO 52a->53a CI weight."
        ),
        "minimum_accepted_character_weight":
            args.minimum_character_weight,
        "units": "eV",
        "chromophores": CHROMOPHORES,
        "frames": [
            int(frame) for frame in energy["frame"]
        ],
        "physical_scope": (
            "Tracked bright local excited-state energies under "
            "electrostatic water embedding and frozen-solute conditions."
        ),
        "important_limitation": (
            "The PYR5 frozen geometry is structurally distinct from "
            "PYR2-PYR4. The site offset therefore includes static "
            "geometry/site effects and solvent electrostatics."
        ),
        "off_diagonal_couplings": (
            "Not included in these diagonal Hamiltonians."
        ),
    }

    (TRAJ_DIR / "metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )
    (HAM_DIR / "metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )

    report = TRAJ_DIR / "TRACKED_STATE_PRODUCTS_DAY018.md"

    with report.open("w", encoding="utf-8") as handle:
        handle.write("# Day018 tracked-state products\n\n")
        handle.write(f"- Calculations represented: {len(df)}\n")
        handle.write(f"- Complete frames: {len(energy)}\n")
        handle.write(
            "- State definition: largest 52a->53a configuration "
            "weight among S1 and S2.\n"
        )
        handle.write(
            f"- Minimum accepted character weight: "
            f"{args.minimum_character_weight:.3f}\n"
        )
        handle.write(
            f"- Mean PYR5 offset relative to PYR2-PYR4: "
            f"{mean_offset_mev:.1f} meV\n\n"
        )
        handle.write("## Tracked-root counts\n\n")
        handle.write(root_counts.to_string(index=False))
        handle.write("\n\n## Cluster statistics\n\n")
        handle.write(cluster_stats.to_string(index=False))
        handle.write("\n\n## Tracked energy trajectory\n\n")
        handle.write(energy.to_string(index=False))
        handle.write("\n")

    print("Tracked-state products generated successfully.")
    print(f"Validated calculations: {len(df)}")
    print(f"Complete frames: {len(energy)}")
    print(f"Mean PYR5 offset: {mean_offset_mev:.3f} meV")
    print("\nTracked-root counts:")
    print(root_counts.to_string(index=False))
    print(f"\nTrajectory directory: {TRAJ_DIR}")
    print(f"Hamiltonian directory: {HAM_DIR}")


if __name__ == "__main__":
    main()
