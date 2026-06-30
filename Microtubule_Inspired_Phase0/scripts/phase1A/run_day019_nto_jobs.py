#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run Day019 NTO ORCA jobs sequentially. "
            "Use --only for the required pilot before the full batch."
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
            "ok_core": False,
        }

    text = path.read_text(errors="ignore")

    normal = NORMAL_MARKER in text
    scf = SCF_MARKER in text
    tddft = TDDFT_MARKER in text
    nto_text = any(pattern.search(text) for pattern in NTO_OUTPUT_PATTERNS)

    return {
        "exists": True,
        "normal_termination": normal,
        "scf_converged": scf,
        "tddft_finished": tddft,
        "nto_text_detected": nto_text,
        "ok_core": normal and scf and tddft,
    }


def list_nto_artifacts(job_dir: Path, inp: Path, out: Path) -> list[Path]:
    excluded = {
        inp.resolve(),
        out.resolve(),
    }

    artifacts: list[Path] = []

    for path in sorted(job_dir.iterdir()):
        if not path.is_file():
            continue

        if path.resolve() in excluded:
            continue

        lower = path.name.lower()

        if (
            "nto" in lower
            or "state" in lower
            or lower.endswith(".gbw")
            or lower.endswith(".cis")
            or lower.endswith(".densities")
            or lower.endswith(".densitiesinfo")
        ):
            artifacts.append(path)

    return artifacts


def archive_incomplete_output(path: Path) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    archived = path.with_name(
        f"{path.stem}.incomplete_{stamp}{path.suffix}"
    )
    path.rename(archived)
    return archived


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
        "n_nto_candidate_artifacts",
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
            "Normal | SCF | TDDFT | NTO text | Candidate artifacts | Core OK |\n"
        )
        handle.write(
            "|---|---:|---|---|---|---:|---:|---:|---:|---:|---:|---:|\n"
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
                f"| {row['nto_text_detected']} "
                f"| {row['n_nto_candidate_artifacts']} "
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

    print("Day019 NTO ORCA runner")
    print(f"ORCA executable: {orca_exe}")
    print(f"Jobs selected: {len(selected)}")

    summary_rows: list[dict[str, object]] = []

    for row in selected:
        job = row["target_job"]
        job_dir, inp = validate_input(row)
        out = job_dir / f"{job}.out"

        previous = inspect_output(out)

        if previous["ok_core"]:
            artifacts = list_nto_artifacts(job_dir, inp, out)
            print(
                f"SKIP {job}: already core-complete; "
                f"candidate_artifacts={len(artifacts)}"
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
                    "n_nto_candidate_artifacts": len(artifacts),
                    "ok_core": previous["ok_core"],
                    "output": str(out.relative_to(PROJECT_ROOT)),
                }
            )
            continue

        if previous["exists"]:
            if not args.rerun_incomplete:
                raise SystemExit(
                    f"Incomplete output already exists: {out}\n"
                    "Inspect it or rerun with --rerun-incomplete."
                )

            archived = archive_incomplete_output(out)
            print(f"ARCHIVE {job}: {archived.name}")

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
        artifacts = list_nto_artifacts(job_dir, inp, out)

        print(
            f"DONE {job} "
            f"returncode={completed.returncode} "
            f"core_ok={result['ok_core']} "
            f"nto_text={result['nto_text_detected']} "
            f"candidate_artifacts={len(artifacts)} "
            f"runtime_min={runtime_min:.2f}"
        )

        summary_rows.append(
            {
                "calculation_type": row["calculation_type"],
                "frame": row["frame"],
                "cluster": row["cluster"],
                "job": job,
                "status": (
                    "completed"
                    if result["ok_core"]
                    else "failed"
                ),
                "returncode": completed.returncode,
                "runtime_min": runtime_min,
                "normal_termination": result["normal_termination"],
                "scf_converged": result["scf_converged"],
                "tddft_finished": result["tddft_finished"],
                "nto_text_detected": result["nto_text_detected"],
                "n_nto_candidate_artifacts": len(artifacts),
                "ok_core": result["ok_core"],
                "output": str(out.relative_to(PROJECT_ROOT)),
            }
        )

        write_reports(summary_rows)

        if completed.returncode != 0 or not result["ok_core"]:
            raise SystemExit(
                f"NTO ORCA job failed: {job}\nInspect: {out}"
            )

    write_reports(summary_rows)

    print("\nDay019 NTO runner completed.")
    print(
        f"Core-successful jobs: "
        f"{sum(bool(row['ok_core']) for row in summary_rows)}/"
        f"{len(summary_rows)}"
    )
    print(f"Wrote: {SUMMARY_CSV.relative_to(PROJECT_ROOT)}")
    print(f"Wrote: {SUMMARY_MD.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
