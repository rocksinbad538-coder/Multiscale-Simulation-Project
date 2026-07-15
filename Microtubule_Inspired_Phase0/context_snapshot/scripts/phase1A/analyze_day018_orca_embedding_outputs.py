#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


DEFAULT_ROOT = Path(
    "runs/phase1A/day016_md_bath_extraction/"
    "orca_embedding_pilot_inputs"
)
DEFAULT_OUTDIR = Path(
    "runs/phase1A/day016_md_bath_extraction/"
    "orca_embedding_analysis"
)

STATE_RE = re.compile(
    r"STATE\s+(\d+):\s+E=\s+"
    r"([-+0-9.Ee]+)\s+au\s+"
    r"([-+0-9.Ee]+)\s+eV\s+"
    r"([-+0-9.Ee]+)\s+cm\*\*-1"
)

ABS_RE = re.compile(
    r"0-1A\s+->\s+(\d+)-1A\s+"
    r"([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)\s+"
    r"([-+0-9.Ee]+)"
)

ERROR_PATTERNS = [
    re.compile(r"ORCA finished by error termination", re.I),
    re.compile(r"error termination", re.I),
    re.compile(r"aborting the run", re.I),
    re.compile(r"SCF NOT CONVERGED", re.I),
    re.compile(r"SCF failed", re.I),
    re.compile(r"segmentation fault", re.I),
    re.compile(r"Please increase MaxIter", re.I),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parse completed ORCA electrostatic-embedding TDDFT outputs "
            "and generate a validated production summary."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help="Directory containing frame*_PYR*_embedding job directories.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help="Directory for CSV and Markdown analysis products.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress the full parsed-job table.",
    )
    return parser.parse_args()


def parse_frame_cluster(name: str) -> tuple[int | None, str | None]:
    match = re.fullmatch(r"frame(\d+)_([^_]+)_embedding", name)
    if not match:
        return None, None
    return int(match.group(1)), match.group(2)


def parse_pc(pc_path: Path) -> tuple[int | None, float | None]:
    if not pc_path.is_file():
        return None, None

    lines = [
        line.strip()
        for line in pc_path.read_text(errors="ignore").splitlines()
        if line.strip()
    ]
    if not lines:
        return None, None

    try:
        declared = int(float(lines[0].split()[0]))
    except (ValueError, IndexError):
        declared = None

    charge_sum = 0.0
    parsed = 0

    for line in lines[1:]:
        fields = line.split()
        if len(fields) != 4:
            continue
        try:
            charge_sum += float(fields[0])
            parsed += 1
        except ValueError:
            continue

    return declared if declared is not None else parsed, charge_sum


def parse_output(out_path: Path) -> dict:
    frame, cluster = parse_frame_cluster(out_path.parent.name)
    if frame is None or cluster is None:
        raise ValueError(f"Cannot parse frame/cluster from {out_path.parent.name}")

    text = out_path.read_text(errors="ignore")

    terminated_normally = "ORCA TERMINATED NORMALLY" in text
    tddft_finished = "ORCA-CIS/TD-DFT FINISHED WITHOUT ERROR" in text
    scf_converged = "SCF CONVERGED" in text
    has_error_flag = any(pattern.search(text) for pattern in ERROR_PATTERNS)

    final_energy = None
    for line in text.splitlines():
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                final_energy = float(line.split()[-1])
            except (ValueError, IndexError):
                pass

    total_runtime_min = None
    runtime_match = re.search(
        r"TOTAL RUN TIME:\s+"
        r"(\d+)\s+days\s+"
        r"(\d+)\s+hours\s+"
        r"(\d+)\s+minutes\s+"
        r"([0-9.]+)\s+seconds",
        text,
    )
    if runtime_match:
        days, hours, minutes, seconds = runtime_match.groups()
        total_runtime_min = (
            int(days) * 1440
            + int(hours) * 60
            + int(minutes)
            + float(seconds) / 60.0
        )

    pc_path = out_path.parent / f"{out_path.stem}.pc"
    n_pc_file, pc_total_charge = parse_pc(pc_path)

    pc_reads = re.findall(r"ok \((\d+) point charges\)", text)
    n_pc_orca = int(pc_reads[-1]) if pc_reads else None

    row = {
        "frame": frame,
        "cluster": cluster,
        "job": out_path.parent.name,
        "terminated_normally": terminated_normally,
        "scf_converged": scf_converged,
        "tddft_finished": tddft_finished,
        "has_error_flag": has_error_flag,
        "final_single_point_energy_Eh": final_energy,
        "n_point_charges_file": n_pc_file,
        "n_point_charges_orca": n_pc_orca,
        "point_charge_total": pc_total_charge,
        "total_runtime_min": total_runtime_min,
        "out_file": str(out_path),
    }

    for state_id, energy_au, energy_ev, energy_cm in STATE_RE.findall(text)[:10]:
        state = int(state_id)
        row[f"S{state}_au"] = float(energy_au)
        row[f"S{state}_eV"] = float(energy_ev)
        row[f"S{state}_cm-1"] = float(energy_cm)

    in_absorption_block = False
    for line in text.splitlines():
        if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in line:
            in_absorption_block = True
            continue

        if (
            in_absorption_block
            and "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY" in line
        ):
            break

        if in_absorption_block:
            match = ABS_RE.search(line)
            if match:
                state, energy_ev, energy_cm, wavelength_nm, oscillator_strength = (
                    match.groups()
                )
                state_i = int(state)
                row[f"f{state_i}"] = float(oscillator_strength)
                row[f"lambda{state_i}_nm"] = float(wavelength_nm)

    dipole_match = re.search(
        r"Magnitude \(Debye\)\s+:\s+([-+0-9.Ee]+)",
        text,
    )
    row["dipole_D"] = (
        float(dipole_match.group(1)) if dipole_match else None
    )

    return row


def main() -> None:
    args = parse_args()

    if not args.root.is_dir():
        raise SystemExit(f"Input root does not exist: {args.root}")

    outputs = sorted(
        args.root.glob(
            "frame*_PYR*_embedding/frame*_PYR*_embedding.out"
        )
    )
    if not outputs:
        raise SystemExit(f"No ORCA outputs found under: {args.root}")

    rows = [parse_output(path) for path in outputs]
    df = pd.DataFrame(rows).sort_values(["frame", "cluster"]).reset_index(drop=True)

    duplicate_mask = df.duplicated(["frame", "cluster"], keep=False)
    if duplicate_mask.any():
        raise SystemExit(
            "Duplicate frame/cluster outputs detected:\n"
            + df.loc[duplicate_mask, ["frame", "cluster", "out_file"]]
            .to_string(index=False)
        )

    args.outdir.mkdir(parents=True, exist_ok=True)

    csv_path = args.outdir / "embedding_pilot_summary.csv"
    df.to_csv(csv_path, index=False)

    success_mask = (
        df["terminated_normally"]
        & df["scf_converged"]
        & df["tddft_finished"]
        & ~df["has_error_flag"]
    )
    successful = int(success_mask.sum())

    audit_path = args.outdir / "EMBEDDING_PRODUCTION_AUDIT_DAY018.md"
    with audit_path.open("w", encoding="utf-8") as handle:
        handle.write("# Day018 ORCA embedding production audit\n\n")
        handle.write(f"- Jobs parsed: {len(df)}\n")
        handle.write(
            f"- Fully successful embedded TDDFT jobs: "
            f"{successful}/{len(df)}\n"
        )
        handle.write(
            f"- Normal ORCA terminations: "
            f"{int(df['terminated_normally'].sum())}/{len(df)}\n"
        )
        handle.write(
            f"- SCF converged: "
            f"{int(df['scf_converged'].sum())}/{len(df)}\n"
        )
        handle.write(
            f"- TDDFT/TDA completed: "
            f"{int(df['tddft_finished'].sum())}/{len(df)}\n"
        )
        handle.write(
            f"- Explicit error signatures: "
            f"{int(df['has_error_flag'].sum())}\n"
        )

        if len(df):
            handle.write(
                f"- Point charges read by ORCA: "
                f"{int(df['n_point_charges_orca'].min())}–"
                f"{int(df['n_point_charges_orca'].max())}\n"
            )
            handle.write(
                f"- S1 range: "
                f"{df['S1_eV'].min():.3f}–"
                f"{df['S1_eV'].max():.3f} eV\n"
            )
            handle.write(
                f"- Mean S1 across all sites and frames: "
                f"{df['S1_eV'].mean():.3f} eV\n"
            )
            handle.write(
                f"- Standard deviation across all sites and frames: "
                f"{df['S1_eV'].std(ddof=1):.3f} eV\n"
            )

        handle.write("\n## Parsed jobs\n\n")
        columns = [
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
        handle.write(df[columns].to_string(index=False))
        handle.write("\n")

    display_columns = [
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

    if not args.quiet:
        print(df[display_columns].to_string(index=False))

    print(f"Jobs parsed: {len(df)}")
    print(f"Fully successful: {successful}/{len(df)}")
    print(f"Explicit error signatures: {int(df['has_error_flag'].sum())}")
    print(f"Wrote: {csv_path}")
    print(f"Wrote: {audit_path}")

    if successful != len(df):
        failed = df.loc[~success_mask, display_columns]
        raise SystemExit(
            "One or more parsed calculations failed production QC:\n"
            + failed.to_string(index=False)
        )


if __name__ == "__main__":
    main()
