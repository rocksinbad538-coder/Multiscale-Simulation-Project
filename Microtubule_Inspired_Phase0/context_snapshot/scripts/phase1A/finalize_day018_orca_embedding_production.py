#!/usr/bin/env python3

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction"
)

ORCA_ROOT = ROOT / "orca_embedding_pilot_inputs"

ANALYSIS_CSV = (
    ROOT
    / "orca_embedding_analysis"
    / "embedding_pilot_summary.csv"
)

ADIABATIC_TRAJECTORY = (
    ROOT
    / "site_energy_trajectory"
    / "site_energy_trajectory.csv"
)

STATE_TRACKING_CSV = (
    ROOT
    / "state_identity_analysis"
    / "low_state_identity_tracking.csv"
)

TRACKED_TRAJECTORY = (
    ROOT
    / "tracked_site_energy_trajectory"
    / "tracked_site_energy_trajectory.csv"
)

TRACKED_ROOTS = (
    ROOT
    / "tracked_site_energy_trajectory"
    / "tracked_root_trajectory.csv"
)

SUMMARY_DIR = ROOT / "day018_production_summary"
SUMMARY_REPORT = (
    SUMMARY_DIR / "PRODUCTION_FINALIZATION_DAY018.md"
)

CHROMOPHORES = ["PYR2", "PYR3", "PYR4", "PYR5"]
EXPECTED_FRAMES = list(range(21))

ERROR_PATTERNS = [
    re.compile(r"ORCA finished by error termination", re.I),
    re.compile(r"error termination", re.I),
    re.compile(r"aborting the run", re.I),
    re.compile(r"SCF NOT CONVERGED", re.I),
    re.compile(r"SCF failed", re.I),
    re.compile(r"segmentation fault", re.I),
]


def expected_output(frame: int, chromophore: str) -> Path:
    job = f"frame{frame:03d}_{chromophore}_embedding"
    return ORCA_ROOT / job / f"{job}.out"


def preflight_orca_outputs() -> list[dict]:
    rows = []
    missing = []
    failed = []

    for frame in EXPECTED_FRAMES:
        for chromophore in CHROMOPHORES:
            path = expected_output(frame, chromophore)

            if not path.is_file():
                missing.append(str(path))
                continue

            text = path.read_text(errors="ignore")

            normal = "ORCA TERMINATED NORMALLY" in text
            scf = "SCF CONVERGED" in text
            tddft = (
                "ORCA-CIS/TD-DFT FINISHED WITHOUT ERROR"
                in text
            )
            error = any(
                pattern.search(text)
                for pattern in ERROR_PATTERNS
            )

            row = {
                "frame": frame,
                "cluster": chromophore,
                "output": str(path.relative_to(PROJECT_ROOT)),
                "normal_termination": normal,
                "scf_converged": scf,
                "tddft_finished": tddft,
                "explicit_error": error,
            }
            rows.append(row)

            if not (normal and scf and tddft and not error):
                failed.append(row)

    if missing:
        raise SystemExit(
            "Production is incomplete. Missing outputs:\n"
            + "\n".join(missing)
        )

    if failed:
        raise SystemExit(
            "One or more ORCA outputs failed preflight QC:\n"
            + pd.DataFrame(failed).to_string(index=False)
        )

    if len(rows) != 84:
        raise SystemExit(
            f"Expected 84 outputs, found {len(rows)}."
        )

    return rows


def run_script(script_name: str, *arguments: str) -> None:
    command = [
        sys.executable,
        str(PROJECT_ROOT / "scripts/phase1A" / script_name),
        *arguments,
    ]

    print("\nRUN:", " ".join(command))

    subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        check=True,
    )


def require_table(
    path: Path,
    expected_rows: int,
    expected_frames: int,
) -> pd.DataFrame:
    if not path.is_file():
        raise SystemExit(f"Expected product not found: {path}")

    table = pd.read_csv(path)

    if len(table) != expected_rows:
        raise SystemExit(
            f"{path}: expected {expected_rows} rows, "
            f"found {len(table)}."
        )

    if "frame" not in table.columns:
        raise SystemExit(f"{path}: missing frame column.")

    frames = sorted(
        pd.to_numeric(
            table["frame"],
            errors="raise",
        ).astype(int).unique()
    )

    if len(frames) != expected_frames:
        raise SystemExit(
            f"{path}: expected {expected_frames} unique frames, "
            f"found {len(frames)}."
        )

    if frames != EXPECTED_FRAMES:
        raise SystemExit(
            f"{path}: unexpected frame set: {frames}"
        )

    return table


def main() -> None:
    print("Day018 ORCA production finalization")
    print("Expected calculations: 84")
    print("Expected complete frames: 21")

    preflight = preflight_orca_outputs()

    print(
        f"Preflight passed: {len(preflight)}/84 "
        "ORCA outputs complete."
    )

    run_script(
        "analyze_day018_orca_embedding_outputs.py",
        "--quiet",
    )

    run_script(
        "build_day018_site_energy_products.py",
    )

    run_script(
        "audit_day018_pyrene_structure_and_state_identity.py",
    )

    run_script(
        "build_day018_tracked_state_products.py",
    )

    analysis = require_table(
        ANALYSIS_CSV,
        expected_rows=84,
        expected_frames=21,
    )

    adiabatic = require_table(
        ADIABATIC_TRAJECTORY,
        expected_rows=21,
        expected_frames=21,
    )

    tracking = require_table(
        STATE_TRACKING_CSV,
        expected_rows=84,
        expected_frames=21,
    )

    tracked = require_table(
        TRACKED_TRAJECTORY,
        expected_rows=21,
        expected_frames=21,
    )

    roots = require_table(
        TRACKED_ROOTS,
        expected_rows=21,
        expected_frames=21,
    )

    count_mismatch = analysis.loc[
        analysis["n_point_charges_file"]
        != analysis["n_point_charges_orca"]
    ]

    nonneutral = analysis.loc[
        analysis["point_charge_total"].abs() > 1.0e-8
    ]

    if len(count_mismatch):
        raise SystemExit(
            "Point-charge count mismatches detected:\n"
            + count_mismatch[
                [
                    "frame",
                    "cluster",
                    "n_point_charges_file",
                    "n_point_charges_orca",
                ]
            ].to_string(index=False)
        )

    if len(nonneutral):
        raise SystemExit(
            "Non-neutral embedding files detected:\n"
            + nonneutral[
                [
                    "frame",
                    "cluster",
                    "point_charge_total",
                ]
            ].to_string(index=False)
        )

    expected_root_pattern = {
        "PYR2": 2,
        "PYR3": 2,
        "PYR4": 2,
        "PYR5": 1,
    }

    root_pattern_changes = []

    for cluster, expected_root in expected_root_pattern.items():
        observed = sorted(
            tracking.loc[
                tracking["cluster"] == cluster,
                "tracked_root",
            ]
            .astype(int)
            .unique()
            .tolist()
        )

        if observed != [expected_root]:
            root_pattern_changes.append(
                {
                    "cluster": cluster,
                    "expected_root": expected_root,
                    "observed_roots": observed,
                }
            )

    energy_stats = (
        tracking.groupby("cluster")
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
        )
        .reset_index()
    )

    reference_mean = tracked[
        ["PYR2", "PYR3", "PYR4"]
    ].mean(axis=1)

    pyr5_offset_mev = (
        1000.0 * (reference_mean - tracked["PYR5"])
    )

    absolute_csv = list(
        (
            ROOT / "tracked_hamiltonian_diagonals"
        ).glob("Htracked_frame*.csv")
    )

    absolute_npy = list(
        (
            ROOT / "tracked_hamiltonian_diagonals"
        ).glob("Htracked_frame*.npy")
    )

    centered_csv = list(
        (
            ROOT / "tracked_hamiltonian_diagonals"
        ).glob("Htracked_centered_frame*.csv")
    )

    centered_npy = list(
        (
            ROOT / "tracked_hamiltonian_diagonals"
        ).glob("Htracked_centered_frame*.npy")
    )

    matrix_counts = {
        "absolute_csv": len(absolute_csv),
        "absolute_npy": len(absolute_npy),
        "centered_csv": len(centered_csv),
        "centered_npy": len(centered_npy),
    }

    if any(count != 21 for count in matrix_counts.values()):
        raise SystemExit(
            f"Unexpected tracked Hamiltonian counts: "
            f"{matrix_counts}"
        )

    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)

    with SUMMARY_REPORT.open("w", encoding="utf-8") as handle:
        handle.write(
            "# Day018 embedded TDDFT production finalization\n\n"
        )

        handle.write("## Production status\n\n")
        handle.write("- ORCA calculations: 84/84\n")
        handle.write("- Complete MD frames: 21/21\n")
        handle.write("- Chromophores per frame: 4\n")
        handle.write("- Normal ORCA terminations: 84/84\n")
        handle.write("- SCF converged: 84/84\n")
        handle.write("- TDDFT/TDA completed: 84/84\n")
        handle.write("- Point-charge count mismatches: 0\n")
        handle.write("- Non-neutral embedding files: 0\n\n")

        handle.write("## Physical scope\n\n")
        handle.write(
            "The source MD trajectory has frozen h-BN and pyrene "
            "coordinates. The time-dependent variation therefore "
            "represents solvent-induced electrostatic fluctuations "
            "for fixed solute geometries.\n\n"
        )

        handle.write("## State tracking\n\n")
        handle.write(
            "The local bright state is selected from S1 and S2 "
            "using the largest 52a->53a HOMO-LUMO configuration "
            "weight rather than root number alone.\n\n"
        )

        handle.write(
            tracking.groupby(
                ["cluster", "tracked_root"]
            )
            .size()
            .rename("n_jobs")
            .reset_index()
            .to_string(index=False)
        )
        handle.write("\n\n")

        if root_pattern_changes:
            handle.write(
                "Root-pattern changes relative to the first "
                "12-frame dataset were detected:\n\n"
            )
            handle.write(
                pd.DataFrame(
                    root_pattern_changes
                ).to_string(index=False)
            )
            handle.write("\n\n")
        else:
            handle.write(
                "The previously observed root-ordering pattern "
                "remains unchanged across all 21 frames.\n\n"
            )

        handle.write("## Tracked-state statistics\n\n")
        handle.write(energy_stats.to_string(index=False))
        handle.write("\n\n")

        handle.write("## PYR5 site offset\n\n")
        handle.write(
            f"- Mean offset relative to PYR2-PYR4: "
            f"{pyr5_offset_mev.mean():.3f} meV\n"
        )
        handle.write(
            f"- Standard deviation: "
            f"{pyr5_offset_mev.std(ddof=1):.3f} meV\n"
        )
        handle.write(
            f"- Minimum: {pyr5_offset_mev.min():.3f} meV\n"
        )
        handle.write(
            f"- Maximum: {pyr5_offset_mev.max():.3f} meV\n\n"
        )

        handle.write("## Generated Hamiltonians\n\n")
        handle.write(
            f"- Absolute CSV matrices: "
            f"{matrix_counts['absolute_csv']}\n"
        )
        handle.write(
            f"- Absolute NPY matrices: "
            f"{matrix_counts['absolute_npy']}\n"
        )
        handle.write(
            f"- Centered CSV matrices: "
            f"{matrix_counts['centered_csv']}\n"
        )
        handle.write(
            f"- Centered NPY matrices: "
            f"{matrix_counts['centered_npy']}\n\n"
        )

        handle.write("## Interpretation constraint\n\n")
        handle.write(
            "The tracked diagonal energies are electronically "
            "consistent, but the PYR5 offset cannot yet be assigned "
            "exclusively to solvent electrostatics because its frozen "
            "geometry differs measurably from PYR2-PYR4. Vacuum "
            "reference calculations and representative NTO analysis "
            "remain necessary.\n"
        )

    print("\nDay018 finalization completed successfully.")
    print("ORCA calculations: 84/84")
    print("Complete frames: 21/21")
    print(
        f"Mean tracked PYR5 offset: "
        f"{pyr5_offset_mev.mean():.3f} meV"
    )
    print(f"Wrote: {SUMMARY_REPORT}")


if __name__ == "__main__":
    main()
