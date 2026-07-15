#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import pty
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/day019_nto_inputs"
)

MANIFEST = ROOT / "NTO_INPUT_MANIFEST_DAY019.csv"

SUMMARY_CSV = ROOT / "NTO_RUN_SUMMARY_DAY019.csv"
SUMMARY_MD = ROOT / "NTO_RUN_SUMMARY_DAY019.md"

NORMAL_MARKER = "ORCA TERMINATED NORMALLY"
SCF_MARKER = "SCF CONVERGED"
TDDFT_MARKER = "ORCA-CIS/TD-DFT FINISHED WITHOUT ERROR"

NTO_OUTPUT_PATTERNS = [
    re.compile(r"NATURAL\s+TRANSITION\s+ORBITAL", re.I),
    re.compile(r"\bNTO\b", re.I),
]


def log(message: str = "") -> None:
    print(message, flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run Day019 NTO ORCA jobs sequentially while streaming ORCA "
            "output live to the terminal and saving the same output in "
            "each job-specific .out file."
        )
    )
    parser.add_argument(
        "--orca",
        default=os.environ.get("ORCA_EXE", "orca"),
        help="ORCA executable or absolute path.",
    )
    parser.add_argument(
        "--only",
        help=(
            "Run exactly one target job from "
            "NTO_INPUT_MANIFEST_DAY019.csv."
        ),
    )
    parser.add_argument(
        "--rerun-incomplete",
        action="store_true",
        help="Archive and rerun an incomplete pre-existing output.",
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


def load_manifest() -> list[dict[str, str]]:
    if not MANIFEST.is_file():
        raise SystemExit(f"Missing manifest: {MANIFEST}")

    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if len(rows) != 8:
        raise SystemExit(
            f"Expected 8 NTO manifest rows, found {len(rows)}."
        )

    jobs = [row["target_job"] for row in rows]

    if len(set(jobs)) != len(jobs):
        raise SystemExit("Duplicate target_job values in NTO manifest.")

    return rows


def validate_input(row: dict[str, str]) -> tuple[Path, Path]:
    job = row["target_job"]
    job_dir = ROOT / job
    inp = job_dir / f"{job}.inp"

    if not job_dir.is_dir():
        raise RuntimeError(f"Missing job directory: {job_dir}")

    if not inp.is_file():
        raise RuntimeError(f"Missing ORCA input: {inp}")

    text = inp.read_text(errors="ignore")

    required = {
        "DoNTO": r"(?mi)^\s*DoNTO\s+true\s*$",
        "NTOStates": r"(?mi)^\s*NTOStates\s+1\s*,\s*2\s*$",
        "NTOThresh": r"(?mi)^\s*NTOThresh\s+1e-4\s*$",
    }

    for label, pattern in required.items():
        matches = re.findall(pattern, text)
        if len(matches) != 1:
            raise RuntimeError(
                f"{job}: expected exactly one {label} directive."
            )

    pointcharge_matches = re.findall(
        r'(?mi)^\s*%pointcharges\s+"([^"]+)"\s*$',
        text,
    )

    calculation_type = row["calculation_type"]

    if calculation_type == "embedded":
        if len(pointcharge_matches) != 1:
            raise RuntimeError(
                f"{job}: embedded input must contain one %pointcharges."
            )

        pc_path = job_dir / pointcharge_matches[0]
        if not pc_path.is_file():
            raise RuntimeError(
                f"{job}: missing point-charge file {pc_path.name}."
            )

    elif calculation_type == "vacuum_reference":
        if pointcharge_matches:
            raise RuntimeError(
                f"{job}: vacuum input must not contain %pointcharges."
            )

    else:
        raise RuntimeError(
            f"{job}: unsupported calculation_type={calculation_type}"
        )

    return job_dir, inp


def inspect_output(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {
            "exists": False,
            "normal_termination": False,
            "scf_converged": False,
            "tddft_finished": False,
            "nto_text_detected": False,
            "s1_nto_section": False,
            "s2_nto_section": False,
            "ok_core": False,
        }

    text = path.read_text(errors="ignore")

    normal = NORMAL_MARKER in text
    scf = SCF_MARKER in text
    tddft = TDDFT_MARKER in text
    nto_text = any(pattern.search(text) for pattern in NTO_OUTPUT_PATTERNS)
    s1_section = "NATURAL TRANSITION ORBITALS FOR STATE    1" in text
    s2_section = "NATURAL TRANSITION ORBITALS FOR STATE    2" in text

    return {
        "exists": True,
        "normal_termination": normal,
        "scf_converged": scf,
        "tddft_finished": tddft,
        "nto_text_detected": nto_text,
        "s1_nto_section": s1_section,
        "s2_nto_section": s2_section,
        "ok_core": (
            normal
            and scf
            and tddft
            and nto_text
            and s1_section
            and s2_section
        ),
    }


def nto_file_counts(job_dir: Path) -> tuple[int, int]:
    n_s1 = len(list(job_dir.glob("*.s1.nto")))
    n_s2 = len(list(job_dir.glob("*.s2.nto")))
    return n_s1, n_s2


def archive_incomplete_output(path: Path) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    archived = path.with_name(
        f"{path.stem}.incomplete_{stamp}{path.suffix}"
    )
    path.rename(archived)
    return archived


def stream_orca_live(
    orca_exe: str,
    inp: Path,
    job_dir: Path,
    out: Path,
) -> int:
    """
    Run ORCA attached to a pseudo-terminal so output is emitted live.
    Every byte is written simultaneously to:
      1. the current terminal/stdout;
      2. the job-specific ORCA .out file.

    An outer shell `tee` can therefore also capture the complete batch log.
    """
    master_fd, slave_fd = pty.openpty()

    try:
        process = subprocess.Popen(
            [orca_exe, inp.name],
            cwd=job_dir,
            stdin=subprocess.DEVNULL,
            stdout=slave_fd,
            stderr=slave_fd,
            close_fds=True,
        )
    finally:
        os.close(slave_fd)

    with out.open("wb") as output_handle:
        while True:
            try:
                chunk = os.read(master_fd, 8192)
            except OSError:
                break

            if not chunk:
                break

            output_handle.write(chunk)
            output_handle.flush()

            sys.stdout.buffer.write(chunk)
            sys.stdout.buffer.flush()

    os.close(master_fd)
    return process.wait()


def write_reports(rows: list[dict[str, object]]) -> None:
    if not rows:
        return

    fieldnames = [
        "calculation_type",
        "frame",
        "cluster",
        "job",
        "status",
        "returncode",
        "runtime_min",
        "normal_termination",
        "scf_converged",
        "tddft_finished",
        "nto_text_detected",
        "s1_nto_section",
        "s2_nto_section",
        "n_s1_nto_files",
        "n_s2_nto_files",
        "ok_core",
        "output",
    ]

    with SUMMARY_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    with SUMMARY_MD.open("w", encoding="utf-8") as handle:
        handle.write("# Day019 NTO ORCA run summary\n\n")
        handle.write(
            "| Type | Frame | Site | Job | Status | Runtime (min) | "
            "Normal | SCF | TDDFT | S1 NTO | S2 NTO | "
            "S1 files | S2 files | OK |\n"
        )
        handle.write(
            "|---|---:|---|---|---|---:|---:|---:|---:|---:|---:|"
            "---:|---:|---:|\n"
        )

        for row in rows:
            handle.write(
                f"| {row['calculation_type']} "
                f"| {row['frame']} "
                f"| {row['cluster']} "
                f"| `{row['job']}` "
                f"| {row['status']} "
                f"| {float(row['runtime_min']):.2f} "
                f"| {row['normal_termination']} "
                f"| {row['scf_converged']} "
                f"| {row['tddft_finished']} "
                f"| {row['s1_nto_section']} "
                f"| {row['s2_nto_section']} "
                f"| {row['n_s1_nto_files']} "
                f"| {row['n_s2_nto_files']} "
                f"| {row['ok_core']} |\n"
            )


def main() -> None:
    args = parse_args()
    orca_exe = resolve_executable(args.orca)
    manifest_rows = load_manifest()

    if args.only:
        selected = [
            row
            for row in manifest_rows
            if row["target_job"] == args.only
        ]

        if len(selected) != 1:
            available = "\n".join(
                f"  {row['target_job']}"
                for row in manifest_rows
            )
            raise SystemExit(
                f"Unknown --only job: {args.only}\n"
                f"Available jobs:\n{available}"
            )
    else:
        selected = manifest_rows

    log("Day019 NTO ORCA runner with live terminal streaming")
    log(f"ORCA executable: {orca_exe}")
    log(f"Jobs selected: {len(selected)}")

    summary_rows: list[dict[str, object]] = []

    for job_number, row in enumerate(selected, start=1):
        job = row["target_job"]
        job_dir, inp = validate_input(row)
        out = job_dir / f"{job}.out"

        previous = inspect_output(out)
        n_s1_existing, n_s2_existing = nto_file_counts(job_dir)

        if (
            previous["ok_core"]
            and n_s1_existing == 1
            and n_s2_existing == 1
        ):
            log(
                f"\n[{job_number}/{len(selected)}] SKIP {job}: "
                "already complete with S1 and S2 NTO files"
            )

            summary_rows.append(
                {
                    "calculation_type": row["calculation_type"],
                    "frame": row["frame"],
                    "cluster": row["cluster"],
                    "job": job,
                    "status": "already_complete",
                    "returncode": 0,
                    "runtime_min": 0.0,
                    "normal_termination": previous["normal_termination"],
                    "scf_converged": previous["scf_converged"],
                    "tddft_finished": previous["tddft_finished"],
                    "nto_text_detected": previous["nto_text_detected"],
                    "s1_nto_section": previous["s1_nto_section"],
                    "s2_nto_section": previous["s2_nto_section"],
                    "n_s1_nto_files": n_s1_existing,
                    "n_s2_nto_files": n_s2_existing,
                    "ok_core": previous["ok_core"],
                    "output": str(out.relative_to(PROJECT_ROOT)),
                }
            )
            write_reports(summary_rows)
            continue

        if previous["exists"]:
            if not args.rerun_incomplete:
                raise SystemExit(
                    f"Incomplete output already exists: {out}\n"
                    "Inspect it or rerun with --rerun-incomplete."
                )

            archived = archive_incomplete_output(out)
            log(f"ARCHIVE {job}: {archived.name}")

        log("")
        log("=" * 88)
        log(
            f"[{job_number}/{len(selected)}] RUN {job}"
        )
        log(
            f"Output is streaming live and being saved to: "
            f"{out.relative_to(PROJECT_ROOT)}"
        )
        log("=" * 88)

        start = time.perf_counter()

        returncode = stream_orca_live(
            orca_exe=orca_exe,
            inp=inp,
            job_dir=job_dir,
            out=out,
        )

        runtime_min = (time.perf_counter() - start) / 60.0
        result = inspect_output(out)
        n_s1, n_s2 = nto_file_counts(job_dir)

        fully_ok = (
            returncode == 0
            and bool(result["ok_core"])
            and n_s1 == 1
            and n_s2 == 1
        )

        log("")
        log(
            f"DONE {job} "
            f"returncode={returncode} "
            f"ok={fully_ok} "
            f"S1_NTO_files={n_s1} "
            f"S2_NTO_files={n_s2} "
            f"runtime_min={runtime_min:.2f}"
        )

        summary_rows.append(
            {
                "calculation_type": row["calculation_type"],
                "frame": row["frame"],
                "cluster": row["cluster"],
                "job": job,
                "status": "completed" if fully_ok else "failed",
                "returncode": returncode,
                "runtime_min": runtime_min,
                "normal_termination": result["normal_termination"],
                "scf_converged": result["scf_converged"],
                "tddft_finished": result["tddft_finished"],
                "nto_text_detected": result["nto_text_detected"],
                "s1_nto_section": result["s1_nto_section"],
                "s2_nto_section": result["s2_nto_section"],
                "n_s1_nto_files": n_s1,
                "n_s2_nto_files": n_s2,
                "ok_core": fully_ok,
                "output": str(out.relative_to(PROJECT_ROOT)),
            }
        )

        write_reports(summary_rows)

        if not fully_ok:
            raise SystemExit(
                f"NTO ORCA job failed validation: {job}\n"
                f"Inspect: {out}"
            )

    write_reports(summary_rows)

    successful = sum(bool(row["ok_core"]) for row in summary_rows)

    log("")
    log("Day019 NTO runner completed.")
    log(
        f"Fully successful jobs: "
        f"{successful}/{len(summary_rows)}"
    )
    log(f"Wrote: {SUMMARY_CSV.relative_to(PROJECT_ROOT)}")
    log(f"Wrote: {SUMMARY_MD.relative_to(PROJECT_ROOT)}")

    if successful != len(summary_rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
