#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import time
from datetime import datetime
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
      "day019_vacuum_reference_inputs"
)

CHROMOPHORES = ["PYR2", "PYR3", "PYR4", "PYR5"]

SUMMARY_CSV = ROOT / "VACUUM_REFERENCE_RUN_SUMMARY_DAY019.csv"
SUMMARY_MD = ROOT / "VACUUM_REFERENCE_RUN_SUMMARY_DAY019.md"

NORMAL_MARKER = "ORCA TERMINATED NORMALLY"
SCF_MARKER = "SCF CONVERGED"
TDDFT_MARKER = "ORCA-CIS/TD-DFT FINISHED WITHOUT ERROR"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the four Day019 vacuum-reference ORCA calculations "
            "sequentially and validate their outputs."
        )
    )
    parser.add_argument(
        "--orca",
        default=os.environ.get("ORCA_EXE", "orca"),
        help=(
            "ORCA executable or path. Defaults to ORCA_EXE or 'orca'."
        ),
    )
    parser.add_argument(
        "--rerun-incomplete",
        action="store_true",
        help=(
            "Archive an incomplete pre-existing output and rerun that job."
        ),
    )
    parser.add_argument(
        "--only",
        choices=CHROMOPHORES,
        help="Run only one chromophore.",
    )
    return parser.parse_args()


def resolve_executable(value: str) -> str:
    candidate = Path(value).expanduser()

    if candidate.is_file():
        return str(candidate.resolve())

    located = shutil.which(value)
    if located:
        return located

    raise SystemExit(
        f"ORCA executable not found: {value}\n"
        "Set ORCA_EXE or use --orca /absolute/path/to/orca."
    )


def job_paths(chromophore: str) -> tuple[str, Path, Path, Path]:
    job = f"frame000_{chromophore}_vacuum_reference"
    job_dir = ROOT / job
    inp = job_dir / f"{job}.inp"
    out = job_dir / f"{job}.out"
    return job, job_dir, inp, out


def inspect_output(path: Path) -> dict:
    if not path.is_file():
        return {
            "exists": False,
            "normal_termination": False,
            "scf_converged": False,
            "tddft_finished": False,
            "ok": False,
        }

    text = path.read_text(errors="ignore")

    normal = NORMAL_MARKER in text
    scf = SCF_MARKER in text
    tddft = TDDFT_MARKER in text

    return {
        "exists": True,
        "normal_termination": normal,
        "scf_converged": scf,
        "tddft_finished": tddft,
        "ok": normal and scf and tddft,
    }


def archive_incomplete_output(path: Path) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    archived = path.with_name(f"{path.stem}.incomplete_{stamp}{path.suffix}")
    path.rename(archived)
    return archived


def write_reports(rows: list[dict]) -> None:
    fieldnames = [
        "chromophore",
        "job",
        "status",
        "returncode",
        "runtime_min",
        "output",
        "normal_termination",
        "scf_converged",
        "tddft_finished",
        "ok",
    ]

    with SUMMARY_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    with SUMMARY_MD.open("w", encoding="utf-8") as handle:
        handle.write("# Day019 vacuum-reference ORCA run summary\n\n")
        handle.write("| Chromophore | Status | Runtime (min) | Normal | SCF | TDDFT/TDA | OK |\n")
        handle.write("|---|---|---:|---:|---:|---:|---:|\n")

        for row in rows:
            handle.write(
                f"| {row['chromophore']} "
                f"| {row['status']} "
                f"| {row['runtime_min']:.2f} "
                f"| {row['normal_termination']} "
                f"| {row['scf_converged']} "
                f"| {row['tddft_finished']} "
                f"| {row['ok']} |\n"
            )

        handle.write("\n")
        handle.write(f"- Jobs represented: {len(rows)}\n")
        handle.write(
            f"- Fully successful: "
            f"{sum(bool(row['ok']) for row in rows)}/{len(rows)}\n"
        )


def main() -> None:
    args = parse_args()
    orca_exe = resolve_executable(args.orca)

    selected = [args.only] if args.only else CHROMOPHORES
    rows: list[dict] = []

    print("Day019 vacuum-reference ORCA runner")
    print(f"ORCA executable: {orca_exe}")
    print(f"Jobs selected: {len(selected)}")

    for chromophore in selected:
        job, job_dir, inp, out = job_paths(chromophore)

        if not job_dir.is_dir():
            raise SystemExit(f"Missing job directory: {job_dir}")

        if not inp.is_file():
            raise SystemExit(f"Missing ORCA input: {inp}")

        if "%pointcharges" in inp.read_text(errors="ignore").lower():
            raise SystemExit(
                f"Vacuum input unexpectedly contains %pointcharges: {inp}"
            )

        previous = inspect_output(out)

        if previous["ok"]:
            print(f"SKIP {job}: already complete")
            rows.append(
                {
                    "chromophore": chromophore,
                    "job": job,
                    "status": "already_complete",
                    "returncode": 0,
                    "runtime_min": 0.0,
                    "output": str(out.relative_to(PROJECT_ROOT)),
                    **{
                        key: previous[key]
                        for key in [
                            "normal_termination",
                            "scf_converged",
                            "tddft_finished",
                            "ok",
                        ]
                    },
                }
            )
            continue

        if previous["exists"]:
            if not args.rerun_incomplete:
                raise SystemExit(
                    f"Incomplete output already exists: {out}\n"
                    "Inspect it, or rerun with --rerun-incomplete."
                )

            archived = archive_incomplete_output(out)
            print(
                f"ARCHIVE {job}: "
                f"{archived.name}"
            )

        print(f"RUN {job}")

        start = time.perf_counter()

        with out.open("w", encoding="utf-8") as handle:
            completed = subprocess.run(
                [orca_exe, inp.name],
                cwd=job_dir,
                stdout=handle,
                stderr=subprocess.STDOUT,
                check=False,
            )

        runtime_min = (time.perf_counter() - start) / 60.0
        result = inspect_output(out)

        print(
            f"DONE {job} "
            f"returncode={completed.returncode} "
            f"ok={result['ok']} "
            f"runtime_min={runtime_min:.2f}"
        )

        row = {
            "chromophore": chromophore,
            "job": job,
            "status": "completed" if result["ok"] else "failed",
            "returncode": completed.returncode,
            "runtime_min": runtime_min,
            "output": str(out.relative_to(PROJECT_ROOT)),
            **{
                key: result[key]
                for key in [
                    "normal_termination",
                    "scf_converged",
                    "tddft_finished",
                    "ok",
                ]
            },
        }
        rows.append(row)

        write_reports(rows)

        if completed.returncode != 0 or not result["ok"]:
            raise SystemExit(
                f"ORCA vacuum-reference job failed: {job}\n"
                f"Inspect: {out}"
            )

    write_reports(rows)

    successful = sum(bool(row["ok"]) for row in rows)

    print("\nDay019 vacuum-reference run completed.")
    print(f"Fully successful: {successful}/{len(rows)}")
    print(f"Wrote: {SUMMARY_CSV.relative_to(PROJECT_ROOT)}")
    print(f"Wrote: {SUMMARY_MD.relative_to(PROJECT_ROOT)}")

    if successful != len(rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
