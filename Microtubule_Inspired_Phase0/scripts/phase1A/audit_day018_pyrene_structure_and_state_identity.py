#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("runs/phase1A/day016_md_bath_extraction")
ORCA_ROOT = ROOT / "orca_embedding_pilot_inputs"
XYZ_ROOT = ROOT / "local_qm_clusters"
OUTDIR = ROOT / "state_identity_analysis"

CHROMOPHORES = ["PYR2", "PYR3", "PYR4", "PYR5"]
N_PYRENE_ATOMS = 26
N_CARBON_ATOMS = 16

STATE_RE = re.compile(
    r"STATE\s+(\d+):\s+E=\s+"
    r"([-+0-9.Ee]+)\s+au\s+"
    r"([-+0-9.Ee]+)\s+eV"
)

TRANSITION_RE = re.compile(
    r"(\d+)a\s*->\s*(\d+)a\s*:\s*"
    r"([-+0-9.Ee]+)\s*"
    r"\(c=\s*([-+0-9.Ee]+)\)"
)

ABSORPTION_RE = re.compile(
    r"0-1A\s+->\s+(\d+)-1A\s+"
    r"([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)"
)

JOB_RE = re.compile(
    r"frame(\d+)_([A-Za-z0-9]+)_embedding"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Audit structural equivalence of frozen pyrenes and track "
            "the low-lying TDDFT state by orbital-transition character."
        )
    )
    parser.add_argument(
        "--homo",
        type=int,
        default=52,
        help="ORCA orbital index assigned to the HOMO.",
    )
    parser.add_argument(
        "--lumo",
        type=int,
        default=53,
        help="ORCA orbital index assigned to the LUMO.",
    )
    parser.add_argument(
        "--candidate-roots",
        type=int,
        default=2,
        help="Number of lowest roots considered for state tracking.",
    )
    return parser.parse_args()


def read_pyrene_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    if not path.is_file():
        raise FileNotFoundError(path)

    lines = path.read_text(errors="ignore").splitlines()
    if len(lines) < N_PYRENE_ATOMS + 2:
        raise ValueError(f"XYZ is too short: {path}")

    symbols = []
    coordinates = []

    for line in lines[2 : 2 + N_PYRENE_ATOMS]:
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"Malformed XYZ line in {path}: {line}")

        symbols.append(fields[0])
        coordinates.append(
            [float(fields[1]), float(fields[2]), float(fields[3])]
        )

    expected = ["C"] * 16 + ["H"] * 10
    if symbols != expected:
        raise ValueError(
            f"Unexpected pyrene atom order in {path}: {symbols}"
        )

    return symbols, np.asarray(coordinates, dtype=float)


def distance_matrix(xyz: np.ndarray) -> np.ndarray:
    delta = xyz[:, None, :] - xyz[None, :, :]
    return np.linalg.norm(delta, axis=2)


def kabsch_rmsd(mobile: np.ndarray, reference: np.ndarray) -> float:
    mobile_centered = mobile - mobile.mean(axis=0)
    reference_centered = reference - reference.mean(axis=0)

    covariance = mobile_centered.T @ reference_centered
    u, _, vt = np.linalg.svd(covariance)

    rotation = u @ vt

    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1.0
        rotation = u @ vt

    aligned = mobile_centered @ rotation

    return float(
        np.sqrt(
            np.mean(
                np.sum(
                    (aligned - reference_centered) ** 2,
                    axis=1,
                )
            )
        )
    )


def planarity_rmsd(carbon_xyz: np.ndarray) -> float:
    centered = carbon_xyz - carbon_xyz.mean(axis=0)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    normal = vt[-1]
    distances = centered @ normal
    return float(np.sqrt(np.mean(distances**2)))


def parse_orca_states(path: Path) -> dict:
    text = path.read_text(errors="ignore")
    lines = text.splitlines()

    states: dict[int, dict] = {}
    current_state = None

    for line in lines:
        state_match = STATE_RE.search(line)
        if state_match:
            root = int(state_match.group(1))
            current_state = root
            states[root] = {
                "energy_au": float(state_match.group(2)),
                "energy_eV": float(state_match.group(3)),
                "transitions": [],
                "fosc": np.nan,
            }
            continue

        transition_match = TRANSITION_RE.search(line)
        if transition_match and current_state is not None:
            occupied, virtual, weight, coefficient = (
                transition_match.groups()
            )
            states[current_state]["transitions"].append(
                {
                    "occupied": int(occupied),
                    "virtual": int(virtual),
                    "weight": float(weight),
                    "coefficient": float(coefficient),
                }
            )

    in_electric_absorption = False

    for line in lines:
        if (
            "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC "
            "DIPOLE MOMENTS" in line
        ):
            in_electric_absorption = True
            continue

        if (
            in_electric_absorption
            and "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY"
            in line
        ):
            break

        if in_electric_absorption:
            match = ABSORPTION_RE.search(line)
            if match:
                root = int(match.group(1))
                if root in states:
                    states[root]["fosc"] = float(match.group(5))

    return states


def transition_weight(
    state: dict,
    occupied: int,
    virtual: int,
) -> float:
    return float(
        sum(
            transition["weight"]
            for transition in state["transitions"]
            if transition["occupied"] == occupied
            and transition["virtual"] == virtual
        )
    )


def dominant_transition(state: dict) -> tuple[str, float]:
    transitions = state["transitions"]
    if not transitions:
        return "", np.nan

    dominant = max(
        transitions,
        key=lambda item: item["weight"],
    )

    label = (
        f"{dominant['occupied']}a->"
        f"{dominant['virtual']}a"
    )
    return label, float(dominant["weight"])


def build_structure_audit(
    jobs: list[tuple[int, str, Path]],
) -> pd.DataFrame:
    structures = {}

    for frame, cluster, _ in jobs:
        xyz_path = (
            XYZ_ROOT
            / f"frame{frame:03d}_{cluster}_water5A.xyz"
        )
        _, xyz = read_pyrene_xyz(xyz_path)
        structures[(frame, cluster)] = xyz

    reference = structures[(0, "PYR2")]
    reference_distances = distance_matrix(reference)

    first_by_cluster = {
        cluster: structures[
            min(
                key
                for key in structures
                if key[1] == cluster
            )
        ]
        for cluster in CHROMOPHORES
    }

    rows = []

    for (frame, cluster), xyz in sorted(structures.items()):
        cluster_reference = first_by_cluster[cluster]

        direct_frame_rmsd = float(
            np.sqrt(np.mean(np.sum((xyz - cluster_reference) ** 2, axis=1)))
        )

        distance_difference = (
            distance_matrix(xyz) - reference_distances
        )

        rows.append(
            {
                "frame": frame,
                "cluster": cluster,
                "direct_rmsd_vs_same_cluster_first_frame_A":
                    direct_frame_rmsd,
                "kabsch_rmsd_vs_frame000_PYR2_all_atoms_A":
                    kabsch_rmsd(xyz, reference),
                "kabsch_rmsd_vs_frame000_PYR2_heavy_atoms_A":
                    kabsch_rmsd(
                        xyz[:N_CARBON_ATOMS],
                        reference[:N_CARBON_ATOMS],
                    ),
                "distance_matrix_rms_vs_frame000_PYR2_A":
                    float(np.sqrt(np.mean(distance_difference**2))),
                "distance_matrix_max_abs_vs_frame000_PYR2_A":
                    float(np.max(np.abs(distance_difference))),
                "carbon_planarity_rmsd_A":
                    planarity_rmsd(xyz[:N_CARBON_ATOMS]),
            }
        )

    return pd.DataFrame(rows)


def build_state_audit(
    jobs: list[tuple[int, str, Path]],
    homo: int,
    lumo: int,
    candidate_roots: int,
) -> pd.DataFrame:
    rows = []

    for frame, cluster, out_path in jobs:
        states = parse_orca_states(out_path)

        missing_roots = [
            root
            for root in range(1, candidate_roots + 1)
            if root not in states
        ]
        if missing_roots:
            raise RuntimeError(
                f"Missing roots {missing_roots} in {out_path}"
            )

        candidates = []

        for root in range(1, candidate_roots + 1):
            state = states[root]
            hl_weight = transition_weight(
                state,
                homo,
                lumo,
            )
            dominant_label, dominant_weight = (
                dominant_transition(state)
            )

            candidates.append(
                {
                    "root": root,
                    "energy_eV": state["energy_eV"],
                    "fosc": state["fosc"],
                    "homo_lumo_weight": hl_weight,
                    "dominant_transition": dominant_label,
                    "dominant_weight": dominant_weight,
                }
            )

        tracked = max(
            candidates,
            key=lambda item: item["homo_lumo_weight"],
        )

        alternate = min(
            candidates,
            key=lambda item: item["homo_lumo_weight"],
        )

        row = {
            "frame": frame,
            "time_ps": frame * 5.0,
            "cluster": cluster,
            "tracked_root": tracked["root"],
            "tracked_energy_eV": tracked["energy_eV"],
            "tracked_fosc": tracked["fosc"],
            "tracked_HOMO_LUMO_weight":
                tracked["homo_lumo_weight"],
            "tracked_dominant_transition":
                tracked["dominant_transition"],
            "tracked_dominant_weight":
                tracked["dominant_weight"],
            "alternate_root": alternate["root"],
            "alternate_energy_eV": alternate["energy_eV"],
            "alternate_fosc": alternate["fosc"],
            "alternate_HOMO_LUMO_weight":
                alternate["homo_lumo_weight"],
            "state_separation_eV": abs(
                candidates[1]["energy_eV"]
                - candidates[0]["energy_eV"]
            ),
            "state_separation_meV": 1000.0
            * abs(
                candidates[1]["energy_eV"]
                - candidates[0]["energy_eV"]
            ),
            "HOMO_LUMO_weight_contrast":
                tracked["homo_lumo_weight"]
                - alternate["homo_lumo_weight"],
        }

        for candidate in candidates:
            root = candidate["root"]
            row[f"S{root}_energy_eV"] = candidate["energy_eV"]
            row[f"S{root}_fosc"] = candidate["fosc"]
            row[f"S{root}_HOMO_LUMO_weight"] = (
                candidate["homo_lumo_weight"]
            )
            row[f"S{root}_dominant_transition"] = (
                candidate["dominant_transition"]
            )
            row[f"S{root}_dominant_weight"] = (
                candidate["dominant_weight"]
            )

        rows.append(row)

    return pd.DataFrame(rows).sort_values(
        ["frame", "cluster"]
    )


def main() -> None:
    args = parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)

    jobs = []

    for out_path in sorted(
        ORCA_ROOT.glob(
            "frame*_PYR*_embedding/"
            "frame*_PYR*_embedding.out"
        )
    ):
        match = JOB_RE.fullmatch(out_path.parent.name)
        if not match:
            continue

        jobs.append(
            (
                int(match.group(1)),
                match.group(2),
                out_path,
            )
        )

    if not jobs:
        raise SystemExit("No ORCA outputs were found.")

    structure = build_structure_audit(jobs)
    states = build_state_audit(
        jobs,
        homo=args.homo,
        lumo=args.lumo,
        candidate_roots=args.candidate_roots,
    )

    structure.to_csv(
        OUTDIR / "pyrene_structure_equivalence.csv",
        index=False,
    )

    states.to_csv(
        OUTDIR / "low_state_identity_tracking.csv",
        index=False,
    )

    root_counts = (
        states.groupby(["cluster", "tracked_root"])
        .size()
        .rename("n_jobs")
        .reset_index()
    )

    tracked_stats = (
        states.groupby("cluster")
        .agg(
            n_jobs=("frame", "count"),
            mean_tracked_energy_eV=(
                "tracked_energy_eV",
                "mean",
            ),
            std_tracked_energy_eV=(
                "tracked_energy_eV",
                "std",
            ),
            min_tracked_energy_eV=(
                "tracked_energy_eV",
                "min",
            ),
            max_tracked_energy_eV=(
                "tracked_energy_eV",
                "max",
            ),
            mean_tracked_fosc=("tracked_fosc", "mean"),
            mean_HOMO_LUMO_weight=(
                "tracked_HOMO_LUMO_weight",
                "mean",
            ),
            min_HOMO_LUMO_weight=(
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

    tracked_stats.to_csv(
        OUTDIR / "tracked_state_statistics.csv",
        index=False,
    )

    report = OUTDIR / "STATE_IDENTITY_AUDIT_DAY018.md"

    with report.open("w", encoding="utf-8") as handle:
        handle.write(
            "# Day018 pyrene structure and state-identity audit\n\n"
        )
        handle.write(f"- Jobs audited: {len(states)}\n")
        handle.write(
            f"- Target orbital transition: "
            f"{args.homo}a -> {args.lumo}a\n"
        )
        handle.write(
            f"- Candidate roots: S1–S{args.candidate_roots}\n\n"
        )

        handle.write("## Tracked-root counts\n\n")
        handle.write(root_counts.to_string(index=False))
        handle.write("\n\n## Tracked-state statistics\n\n")
        handle.write(tracked_stats.to_string(index=False))

        handle.write("\n\n## Structural equivalence summary\n\n")
        handle.write(
            structure.groupby("cluster")
            .agg(
                maximum_direct_frame_RMSD_A=(
                    "direct_rmsd_vs_same_cluster_first_frame_A",
                    "max",
                ),
                Kabsch_RMSD_all_atoms_vs_PYR2_A=(
                    "kabsch_rmsd_vs_frame000_PYR2_all_atoms_A",
                    "mean",
                ),
                Kabsch_RMSD_heavy_atoms_vs_PYR2_A=(
                    "kabsch_rmsd_vs_frame000_PYR2_heavy_atoms_A",
                    "mean",
                ),
                maximum_distance_matrix_difference_A=(
                    "distance_matrix_max_abs_vs_frame000_PYR2_A",
                    "max",
                ),
                mean_planarity_RMSD_A=(
                    "carbon_planarity_rmsd_A",
                    "mean",
                ),
            )
            .to_string()
        )

        handle.write(
            "\n\n## Interpretation constraint\n\n"
            "Root number alone is not used as state identity. "
            "The present tracking follows the low-lying state with "
            "the largest HOMO-to-LUMO configuration weight. "
            "Definitive diabatic assignment still requires orbital "
            "or NTO inspection for representative cases.\n"
        )

    print("State-identity audit completed.")
    print(f"Jobs audited: {len(states)}")
    print("\nTracked-root counts:")
    print(root_counts.to_string(index=False))
    print("\nTracked-state statistics:")
    print(tracked_stats.to_string(index=False))
    print("\nMaximum structural deviations:")
    print(
        structure[
            [
                "direct_rmsd_vs_same_cluster_first_frame_A",
                "kabsch_rmsd_vs_frame000_PYR2_heavy_atoms_A",
                "distance_matrix_max_abs_vs_frame000_PYR2_A",
                "carbon_planarity_rmsd_A",
            ]
        ].max().to_string()
    )
    print(f"\nWrote: {OUTDIR}")


if __name__ == "__main__":
    main()
