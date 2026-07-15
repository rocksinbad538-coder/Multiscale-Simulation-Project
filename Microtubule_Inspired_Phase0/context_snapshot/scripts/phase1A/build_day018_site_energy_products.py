#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("runs/phase1A/day016_md_bath_extraction")
DEFAULT_INPUT = ROOT / "orca_embedding_analysis" / "embedding_pilot_summary.csv"
SITE_DIR = ROOT / "site_energy_trajectory"
HAMILTONIAN_DIR = ROOT / "hamiltonian_diagonals"
QC_DIR = ROOT / "site_energy_qc"

CHROMOPHORES = ["PYR2", "PYR3", "PYR4", "PYR5"]
FRAME_SPACING_PS = 5.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate ORCA embedding results and regenerate the site-energy "
            "trajectory, diagonal Hamiltonians, and temporal QC products."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Parsed ORCA embedding summary CSV.",
    )
    parser.add_argument(
        "--step-review-mev",
        type=float,
        default=40.0,
        help=(
            "Descriptive review threshold for consecutive observed-site "
            "energy changes. This is not a physical cutoff."
        ),
    )
    parser.add_argument(
        "--gap-review-mev",
        type=float,
        default=60.0,
        help=(
            "Descriptive review threshold for small S2-S1 gaps. "
            "This is not proof of root switching."
        ),
    )
    return parser.parse_args()


def coerce_bool(series: pd.Series, name: str) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        if series.isna().any():
            raise SystemExit(f"Missing Boolean values in column: {name}")
        return series.astype(bool)

    mapping = {
        "true": True,
        "false": False,
        "1": True,
        "0": False,
    }
    converted = (
        series.astype(str)
        .str.strip()
        .str.lower()
        .map(mapping)
    )
    if converted.isna().any():
        bad = series[converted.isna()].unique().tolist()
        raise SystemExit(f"Invalid Boolean values in {name}: {bad}")
    return converted.astype(bool)


def clean_products(directory: Path, patterns: list[str]) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for pattern in patterns:
        for path in directory.glob(pattern):
            if path.is_file():
                path.unlink()


def validate_input(df: pd.DataFrame) -> pd.DataFrame:
    required = [
        "frame",
        "cluster",
        "terminated_normally",
        "scf_converged",
        "tddft_finished",
        "has_error_flag",
        "n_point_charges_orca",
        "S1_eV",
        "S2_eV",
        "S3_eV",
        "f1",
        "dipole_D",
        "total_runtime_min",
    ]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns: {missing}")

    work = df.copy()

    for column in [
        "terminated_normally",
        "scf_converged",
        "tddft_finished",
        "has_error_flag",
    ]:
        work[column] = coerce_bool(work[column], column)

    work["frame"] = pd.to_numeric(work["frame"], errors="raise").astype(int)

    numeric_columns = [
        "n_point_charges_orca",
        "S1_eV",
        "S2_eV",
        "S3_eV",
        "f1",
        "dipole_D",
        "total_runtime_min",
    ]
    for column in numeric_columns:
        work[column] = pd.to_numeric(work[column], errors="coerce")

    if work[numeric_columns].isna().any().any():
        bad_rows = work.loc[
            work[numeric_columns].isna().any(axis=1),
            ["frame", "cluster"] + numeric_columns,
        ]
        raise SystemExit(
            "Missing or non-numeric production values:\n"
            + bad_rows.to_string(index=False)
        )

    unexpected = sorted(set(work["cluster"]) - set(CHROMOPHORES))
    if unexpected:
        raise SystemExit(f"Unexpected chromophores: {unexpected}")

    duplicate_mask = work.duplicated(["frame", "cluster"], keep=False)
    if duplicate_mask.any():
        raise SystemExit(
            "Duplicate frame/cluster rows:\n"
            + work.loc[
                duplicate_mask,
                ["frame", "cluster"],
            ].to_string(index=False)
        )

    failed = work.loc[
        ~(
            work["terminated_normally"]
            & work["scf_converged"]
            & work["tddft_finished"]
            & ~work["has_error_flag"]
        )
    ]
    if len(failed):
        raise SystemExit(
            "One or more ORCA jobs failed production validation:\n"
            + failed[
                [
                    "frame",
                    "cluster",
                    "terminated_normally",
                    "scf_converged",
                    "tddft_finished",
                    "has_error_flag",
                ]
            ].to_string(index=False)
        )

    expected = set(CHROMOPHORES)
    incomplete = []

    for frame, group in work.groupby("frame"):
        observed = set(group["cluster"])
        if observed != expected or len(group) != len(CHROMOPHORES):
            incomplete.append(
                {
                    "frame": int(frame),
                    "observed": sorted(observed),
                    "missing": sorted(expected - observed),
                    "unexpected": sorted(observed - expected),
                    "rows": len(group),
                }
            )

    if incomplete:
        raise SystemExit(
            "Incomplete frame coverage:\n"
            + json.dumps(incomplete, indent=2)
        )

    if (work["S2_eV"] <= work["S1_eV"]).any():
        bad = work.loc[
            work["S2_eV"] <= work["S1_eV"],
            ["frame", "cluster", "S1_eV", "S2_eV"],
        ]
        raise SystemExit(
            "Detected S2 <= S1:\n"
            + bad.to_string(index=False)
        )

    if (work["n_point_charges_orca"] <= 0).any():
        raise SystemExit("Non-positive point-charge count detected.")

    non_tip4p = work.loc[
        work["n_point_charges_orca"].astype(int) % 4 != 0,
        ["frame", "cluster", "n_point_charges_orca"],
    ]
    if len(non_tip4p):
        raise SystemExit(
            "Point-charge count is not divisible by four for TIP4P/2005:\n"
            + non_tip4p.to_string(index=False)
        )

    return work.sort_values(["frame", "cluster"]).reset_index(drop=True)


def build_site_energy_products(df: pd.DataFrame) -> pd.DataFrame:
    SITE_DIR.mkdir(parents=True, exist_ok=True)

    site = (
        df.pivot(index="frame", columns="cluster", values="S1_eV")
        .reset_index()
    )
    site = site[["frame"] + CHROMOPHORES].sort_values("frame")

    if site[CHROMOPHORES].isna().any().any():
        raise SystemExit("NaN values detected after site-energy pivot.")

    site.to_csv(SITE_DIR / "site_energy_trajectory.csv", index=False)

    centered = site.copy()
    centered[CHROMOPHORES] = site[CHROMOPHORES].sub(
        site[CHROMOPHORES].mean(axis=1),
        axis=0,
    )
    centered.to_csv(
        SITE_DIR / "site_energy_trajectory_centered.csv",
        index=False,
    )

    global_centered = site.copy()
    global_mean = float(site[CHROMOPHORES].to_numpy().mean())
    global_centered[CHROMOPHORES] = (
        site[CHROMOPHORES] - global_mean
    )
    global_centered.to_csv(
        SITE_DIR / "site_energy_trajectory_global_centered.csv",
        index=False,
    )

    np.save(
        SITE_DIR / "site_energy_trajectory_eV.npy",
        site[CHROMOPHORES].to_numpy(dtype=float),
    )
    np.save(
        SITE_DIR / "site_energy_trajectory_centered_eV.npy",
        centered[CHROMOPHORES].to_numpy(dtype=float),
    )

    centered_dir = SITE_DIR / "hamiltonian_diagonal_frames"
    clean_products(
        centered_dir,
        ["H_diagonal_centered_frame*.csv"],
    )

    for _, row in centered.iterrows():
        frame = int(row["frame"])
        matrix = np.diag(
            row[CHROMOPHORES].to_numpy(dtype=float)
        )
        pd.DataFrame(
            matrix,
            index=CHROMOPHORES,
            columns=CHROMOPHORES,
        ).to_csv(
            centered_dir
            / f"H_diagonal_centered_frame{frame:03d}.csv"
        )

    metadata = {
        "production_day": "018",
        "source": str(DEFAULT_INPUT),
        "quantity": "TDA-TDDFT S1 site energies",
        "units": "eV",
        "chromophore_order": CHROMOPHORES,
        "frames": [int(value) for value in site["frame"]],
        "frame_spacing_ps": FRAME_SPACING_PS,
        "embedding": (
            "TIP4P/2005 water represented as electrostatic point charges "
            "within the production 5 Å local-shell construction"
        ),
        "solute_condition": (
            "h-BN scaffold and pyrene chromophores frozen in x, y, and z "
            "during the source MD trajectory"
        ),
        "physical_scope": (
            "solvent-induced site-energy fluctuations under "
            "frozen-solute conditions"
        ),
        "hamiltonian_note": (
            "Only diagonal site-energy terms are generated. "
            "Off-diagonal excitonic couplings J_ij are not included."
        ),
        "centering": {
            "site_energy_trajectory_centered.csv": (
                "Per-frame common energy removed. This subtracts an "
                "identity contribution from each diagonal Hamiltonian."
            ),
            "site_energy_trajectory_global_centered.csv": (
                "Single global mean over all frames and chromophores removed."
            ),
        },
    }
    (SITE_DIR / "metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )

    summary_path = SITE_DIR / "SITE_ENERGY_TRAJECTORY_DAY018.md"
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("# Day018 site-energy trajectory\n\n")
        handle.write(f"- Complete frames: {len(site)}\n")
        handle.write(
            f"- Embedded TDDFT calculations represented: "
            f"{len(site) * len(CHROMOPHORES)}\n"
        )
        handle.write(
            "- Physical scope: solvent-induced site-energy fluctuations "
            "under frozen-solute conditions.\n"
        )
        handle.write(
            "- Hamiltonian scope: diagonal terms only; "
            "off-diagonal couplings are not yet included.\n\n"
        )
        handle.write("## Absolute site energies\n\n")
        handle.write(site.to_string(index=False))
        handle.write("\n\n## Per-frame centered site energies\n\n")
        handle.write(centered.to_string(index=False))
        handle.write("\n\n## Per-chromophore statistics\n\n")
        handle.write(site[CHROMOPHORES].describe().T.to_string())
        handle.write("\n")

    return site


def export_absolute_hamiltonians(site: pd.DataFrame) -> pd.DataFrame:
    clean_products(
        HAMILTONIAN_DIR,
        ["Hdiag_frame*.csv", "Hdiag_frame*.npy"],
    )

    summary_rows = []

    for _, row in site.iterrows():
        frame = int(row["frame"])
        diagonal = row[CHROMOPHORES].to_numpy(dtype=float)
        matrix = np.diag(diagonal)

        np.save(
            HAMILTONIAN_DIR / f"Hdiag_frame{frame:03d}.npy",
            matrix,
        )
        pd.DataFrame(
            matrix,
            index=CHROMOPHORES,
            columns=CHROMOPHORES,
        ).to_csv(
            HAMILTONIAN_DIR / f"Hdiag_frame{frame:03d}.csv"
        )

        summary_rows.append(
            {
                "frame": frame,
                "time_ps": frame * FRAME_SPACING_PS,
                "trace_eV": float(np.trace(matrix)),
                "mean_diag_eV": float(diagonal.mean()),
                "min_diag_eV": float(diagonal.min()),
                "max_diag_eV": float(diagonal.max()),
                "diagonal_spread_eV": float(
                    diagonal.max() - diagonal.min()
                ),
            }
        )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(
        HAMILTONIAN_DIR / "hamiltonian_summary.csv",
        index=False,
    )

    metadata = {
        "production_day": "018",
        "description": (
            "Absolute diagonal excitonic Hamiltonian components extracted "
            "from embedded TDA-TDDFT calculations."
        ),
        "units": "eV",
        "method": "wB97X-D3/def2-SVP TDA-TDDFT electrostatic embedding",
        "hamiltonian_type": "diagonal_only",
        "chromophores": CHROMOPHORES,
        "physical_scope": (
            "solvent-induced diagonal disorder under frozen-solute conditions"
        ),
        "off_diagonal_terms": (
            "Not included. J_ij must be parameterized separately."
        ),
        "source": str(
            SITE_DIR / "site_energy_trajectory.csv"
        ),
    }
    (HAMILTONIAN_DIR / "metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )

    return summary


def build_qc_products(
    df: pd.DataFrame,
    site: pd.DataFrame,
    step_review_mev: float,
    gap_review_mev: float,
) -> None:
    QC_DIR.mkdir(parents=True, exist_ok=True)

    row_qc = df.copy()
    row_qc["time_ps"] = (
        row_qc["frame"] * FRAME_SPACING_PS
    )
    row_qc["S2_minus_S1_eV"] = (
        row_qc["S2_eV"] - row_qc["S1_eV"]
    )
    row_qc["S2_minus_S1_meV"] = (
        1000.0 * row_qc["S2_minus_S1_eV"]
    )
    row_qc["n_water_molecules_embedding"] = (
        row_qc["n_point_charges_orca"] / 4.0
    )

    grouped = row_qc.groupby("cluster", sort=False)
    row_qc["previous_frame"] = grouped["frame"].shift(1)
    row_qc["previous_S1_eV"] = grouped["S1_eV"].shift(1)
    row_qc["delta_frame"] = (
        row_qc["frame"] - row_qc["previous_frame"]
    )
    row_qc["delta_time_ps"] = (
        row_qc["delta_frame"] * FRAME_SPACING_PS
    )
    row_qc["delta_S1_eV"] = (
        row_qc["S1_eV"] - row_qc["previous_S1_eV"]
    )
    row_qc["delta_S1_meV"] = (
        1000.0 * row_qc["delta_S1_eV"]
    )
    row_qc["abs_delta_S1_meV"] = (
        row_qc["delta_S1_meV"].abs()
    )
    row_qc["delta_S1_rate_meV_per_ps"] = (
        row_qc["delta_S1_meV"] / row_qc["delta_time_ps"]
    )

    row_qc["step_review_flag"] = (
        row_qc["abs_delta_S1_meV"] >= step_review_mev
    ).fillna(False)
    row_qc["small_gap_review_flag"] = (
        row_qc["S2_minus_S1_meV"] <= gap_review_mev
    )
    row_qc["state_identity_review_flag"] = (
        row_qc["step_review_flag"]
        | row_qc["small_gap_review_flag"]
    )

    row_qc.to_csv(
        QC_DIR / "site_energy_row_level_qc.csv",
        index=False,
    )

    reference_mean = site[["PYR2", "PYR3", "PYR4"]].mean(axis=1)
    frame_metrics = pd.DataFrame(
        {
            "frame": site["frame"].astype(int),
            "time_ps": site["frame"] * FRAME_SPACING_PS,
            "mean_PYR2_PYR4_eV": reference_mean,
            "PYR5_eV": site["PYR5"],
            "PYR5_redshift_vs_PYR2_PYR4_eV": (
                reference_mean - site["PYR5"]
            ),
            "PYR5_redshift_vs_PYR2_PYR4_meV": (
                1000.0 * (reference_mean - site["PYR5"])
            ),
            "frame_min_eV": site[CHROMOPHORES].min(axis=1),
            "frame_max_eV": site[CHROMOPHORES].max(axis=1),
            "diagonal_spread_eV": (
                site[CHROMOPHORES].max(axis=1)
                - site[CHROMOPHORES].min(axis=1)
            ),
            "diagonal_spread_meV": (
                1000.0
                * (
                    site[CHROMOPHORES].max(axis=1)
                    - site[CHROMOPHORES].min(axis=1)
                )
            ),
        }
    )
    frame_metrics.to_csv(
        QC_DIR / "site_energy_frame_level_qc.csv",
        index=False,
    )

    cluster_stats = (
        row_qc.groupby("cluster")
        .agg(
            n_frames=("frame", "count"),
            mean_S1_eV=("S1_eV", "mean"),
            std_S1_eV=("S1_eV", "std"),
            min_S1_eV=("S1_eV", "min"),
            max_S1_eV=("S1_eV", "max"),
            mean_f1=("f1", "mean"),
            std_f1=("f1", "std"),
            min_S2_minus_S1_meV=("S2_minus_S1_meV", "min"),
            mean_S2_minus_S1_meV=("S2_minus_S1_meV", "mean"),
            mean_point_charges=("n_point_charges_orca", "mean"),
            min_point_charges=("n_point_charges_orca", "min"),
            max_point_charges=("n_point_charges_orca", "max"),
            maximum_observed_step_meV=("abs_delta_S1_meV", "max"),
        )
        .reset_index()
    )
    cluster_stats.to_csv(
        QC_DIR / "site_energy_cluster_statistics.csv",
        index=False,
    )

    correlation_rows = []
    for cluster, group in row_qc.groupby("cluster"):
        correlation_rows.append(
            {
                "cluster": cluster,
                "n": len(group),
                "pearson_S1_vs_point_charges": (
                    group["S1_eV"].corr(
                        group["n_point_charges_orca"]
                    )
                ),
                "pearson_S1_vs_f1": (
                    group["S1_eV"].corr(group["f1"])
                ),
                "pearson_S1_vs_S2_minus_S1": (
                    group["S1_eV"].corr(
                        group["S2_minus_S1_eV"]
                    )
                ),
            }
        )

    correlations = pd.DataFrame(correlation_rows)
    correlations.to_csv(
        QC_DIR / "site_energy_descriptive_correlations.csv",
        index=False,
    )

    review_cases = row_qc.loc[
        row_qc["state_identity_review_flag"],
        [
            "frame",
            "time_ps",
            "cluster",
            "S1_eV",
            "S2_eV",
            "S2_minus_S1_meV",
            "f1",
            "n_point_charges_orca",
            "delta_time_ps",
            "delta_S1_meV",
            "step_review_flag",
            "small_gap_review_flag",
        ],
    ].copy()
    review_cases.to_csv(
        QC_DIR / "state_identity_review_cases.csv",
        index=False,
    )

    top_steps = (
        row_qc.dropna(subset=["abs_delta_S1_meV"])
        .sort_values("abs_delta_S1_meV", ascending=False)
        .head(12)
    )
    top_steps.to_csv(
        QC_DIR / "largest_observed_site_energy_steps.csv",
        index=False,
    )

    report = QC_DIR / "SITE_ENERGY_QC_DAY018.md"
    with report.open("w", encoding="utf-8") as handle:
        handle.write("# Day018 site-energy production QC\n\n")
        handle.write(f"- Complete frames: {site.shape[0]}\n")
        handle.write(f"- Parsed calculations: {len(row_qc)}\n")
        handle.write(
            f"- Review threshold for observed steps: "
            f"{step_review_mev:.1f} meV\n"
        )
        handle.write(
            f"- Review threshold for S2-S1 gaps: "
            f"{gap_review_mev:.1f} meV\n"
        )
        handle.write(
            "- These thresholds are screening criteria only; "
            "they are not physical definitions of root switching.\n"
        )
        handle.write(
            "- State identity cannot be established from excitation "
            "energy and oscillator strength alone. Transition character "
            "or NTO analysis remains necessary for definitive assignment.\n"
        )
        handle.write(
            "- The source MD trajectory contains frozen h-BN and pyrene "
            "coordinates; the resulting variation is interpreted as "
            "solvent-induced diagonal disorder.\n\n"
        )

        handle.write("## Cluster statistics\n\n")
        handle.write(cluster_stats.to_string(index=False))
        handle.write("\n\n## Frame-level diagonal disorder\n\n")
        handle.write(frame_metrics.to_string(index=False))
        handle.write("\n\n## Descriptive correlations\n\n")
        handle.write(correlations.to_string(index=False))
        handle.write("\n\n## State-identity review cases\n\n")
        if len(review_cases):
            handle.write(review_cases.to_string(index=False))
        else:
            handle.write("No rows crossed the declared review thresholds.")
        handle.write("\n\n## Largest observed changes\n\n")
        handle.write(
            top_steps[
                [
                    "frame",
                    "cluster",
                    "delta_time_ps",
                    "delta_S1_meV",
                    "S2_minus_S1_meV",
                    "f1",
                    "n_point_charges_orca",
                ]
            ].to_string(index=False)
        )
        handle.write("\n")


def main() -> None:
    args = parse_args()

    if not args.input.is_file():
        raise SystemExit(f"Input summary not found: {args.input}")

    raw = pd.read_csv(args.input)
    validated = validate_input(raw)

    site = build_site_energy_products(validated)
    hamiltonian_summary = export_absolute_hamiltonians(site)
    build_qc_products(
        validated,
        site,
        step_review_mev=args.step_review_mev,
        gap_review_mev=args.gap_review_mev,
    )

    print("Day018 site-energy products generated successfully.")
    print(f"Validated ORCA calculations: {len(validated)}")
    print(f"Complete frames: {len(site)}")
    print(
        "Frames: "
        + ", ".join(f"{int(frame):03d}" for frame in site["frame"])
    )
    print(
        f"Absolute Hamiltonians exported: "
        f"{len(hamiltonian_summary)}"
    )
    print(f"Site-energy directory: {SITE_DIR}")
    print(f"Hamiltonian directory: {HAMILTONIAN_DIR}")
    print(f"QC directory: {QC_DIR}")


if __name__ == "__main__":
    main()
