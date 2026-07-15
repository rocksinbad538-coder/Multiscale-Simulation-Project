#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
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

NTO_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/day019_nto_inputs"
)

ANALYSIS_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/day019_nto_analysis"
)

STATE_METRICS = ANALYSIS_ROOT / "nto_state_metrics.csv"

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction/"
    "day019_nto_cube_generation"
)

SELECTION_MANIFEST = (
    OUTPUT_ROOT / "CUBE_SELECTION_MANIFEST_DAY019.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT / "CUBE_GENERATION_SUMMARY_DAY019.csv"
)

AUDIT_MD = (
    OUTPUT_ROOT / "CUBE_GENERATION_AUDIT_DAY019.md"
)

PAIR_RE = re.compile(
    r"^\s*(\d+)([ab])\s*->\s*(\d+)([ab])\s*$"
)

MIN_VALID_CUBE_BYTES = 1000
GRID_POINTS = 80


def log(message: str = "") -> None:
    print(message, flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate the 48 selected Day019 NTO cube files with ORCA "
            "orca_plot. All orca_plot output is streamed live to the "
            "terminal and can be captured simultaneously with shell tee."
        )
    )
    parser.add_argument(
        "--orca-plot",
        default=os.environ.get("ORCA_PLOT_EXE", "orca_plot"),
        help="orca_plot executable or absolute path.",
    )
    parser.add_argument(
        "--rerun-invalid",
        action="store_true",
        help=(
            "Archive and regenerate pre-existing cube files that fail "
            "validation."
        ),
    )
    parser.add_argument(
        "--only-job",
        help="Generate cubes only for one target job.",
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
        f"orca_plot executable not found: {value}\n"
        "Set ORCA_PLOT_EXE or use --orca-plot /absolute/path/orca_plot."
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def parse_pair(value: str) -> list[tuple[int, str]]:
    match = PAIR_RE.match(value)

    if not match:
        raise RuntimeError(f"Invalid NTO pair label: {value!r}")

    hole_orbital = int(match.group(1))
    hole_spin = match.group(2)
    particle_orbital = int(match.group(3))
    particle_spin = match.group(4)

    return [
        (hole_orbital, hole_spin),
        (particle_orbital, particle_spin),
    ]


def read_state_metrics() -> list[dict[str, str]]:
    if not STATE_METRICS.is_file():
        raise SystemExit(
            f"Missing NTO state metrics: {STATE_METRICS}"
        )

    with STATE_METRICS.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle))

    if len(rows) != 16:
        raise SystemExit(
            f"Expected 16 NTO state rows, found {len(rows)}."
        )

    return rows


def build_selection(
    rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    selection: list[dict[str, object]] = []

    for row in rows:
        is_tracked = (
            row["is_tracked_root"].strip().lower() == "true"
        )

        pair_labels = [row["dominant_pair"]]

        if not is_tracked:
            pair_labels.append(row["second_pair"])

        seen_orbitals: set[tuple[int, str]] = set()

        for pair_rank, pair_label in enumerate(
            pair_labels,
            start=1,
        ):
            for role, (orbital, spin) in zip(
                ("hole", "particle"),
                parse_pair(pair_label),
            ):
                key = (orbital, spin)

                if key in seen_orbitals:
                    raise RuntimeError(
                        f"Duplicate selected orbital {orbital}{spin} "
                        f"within {row['job']} S{row['root']}"
                    )

                seen_orbitals.add(key)

                nto_file = (
                    PROJECT_ROOT / row["nto_file"]
                ).resolve()

                expected_cube = nto_file.with_name(
                    nto_file.name.removesuffix(".nto")
                    + f".mo{orbital}{spin}.cube"
                )

                selection.append(
                    {
                        "calculation_type": row[
                            "calculation_type"
                        ],
                        "frame": int(row["frame"]),
                        "cluster": row["cluster"],
                        "job": row["job"],
                        "root": int(row["root"]),
                        "is_tracked_root": is_tracked,
                        "pair_rank": pair_rank,
                        "pair_label": pair_label,
                        "orbital_role": role,
                        "orbital_index": orbital,
                        "orbital_spin": spin,
                        "nto_file": str(
                            nto_file.relative_to(PROJECT_ROOT)
                        ),
                        "cube_file": str(
                            expected_cube.relative_to(PROJECT_ROOT)
                        ),
                    }
                )

    if len(selection) != 48:
        raise RuntimeError(
            f"Expected 48 cube selections, found {len(selection)}."
        )

    by_job: dict[str, int] = {}

    for row in selection:
        job = str(row["job"])
        by_job[job] = by_job.get(job, 0) + 1

    if sorted(by_job.values()) != [6] * 8:
        raise RuntimeError(
            "Expected six cube selections for each of eight jobs; "
            f"found {by_job}"
        )

    return selection


def validate_nto_file(path: Path) -> None:
    if not path.is_file():
        raise RuntimeError(f"Missing NTO file: {path}")

    if path.stat().st_size < MIN_VALID_CUBE_BYTES:
        raise RuntimeError(
            f"NTO file is unexpectedly small: {path}"
        )


def validate_cube(path: Path) -> tuple[bool, str]:
    if not path.is_file():
        return False, "missing"

    size = path.stat().st_size

    if size < MIN_VALID_CUBE_BYTES:
        return False, f"too_small:{size}"

    try:
        with path.open(
            "r",
            encoding="utf-8",
            errors="replace",
        ) as handle:
            first_six = [
                handle.readline()
                for _ in range(6)
            ]
    except OSError as exc:
        return False, f"read_error:{exc}"

    if any(line == "" for line in first_six):
        return False, "truncated_header"

    return True, "valid"


def archive_invalid_cube(path: Path) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    archived = path.with_name(
        f"{path.name}.invalid_{stamp}"
    )
    path.rename(archived)
    return archived


def orca_plot_commands(
    orbital: int,
    spin: str,
) -> str:
    operator = 0 if spin == "a" else 1

    commands = [
        "2",
        str(orbital),
        "3",
        str(operator),
        "4",
        str(GRID_POINTS),
        "5",
        "7",
        "8",
        "0",
        "11",
        "12",
        "",
    ]

    return "\n".join(commands)


def run_orca_plot_live(
    executable: str,
    nto_file: Path,
    orbital: int,
    spin: str,
    transcript: Path,
) -> int:
    master_fd, slave_fd = pty.openpty()

    process = subprocess.Popen(
        [executable, nto_file.name, "-i"],
        cwd=nto_file.parent,
        stdin=slave_fd,
        stdout=slave_fd,
        stderr=slave_fd,
        close_fds=True,
    )

    os.close(slave_fd)

    command_bytes = orca_plot_commands(
        orbital=orbital,
        spin=spin,
    ).encode("utf-8")

    os.write(master_fd, command_bytes)

    with transcript.open("wb") as log_handle:
        while True:
            try:
                chunk = os.read(master_fd, 8192)
            except OSError:
                break

            if not chunk:
                break

            log_handle.write(chunk)
            log_handle.flush()

            sys.stdout.buffer.write(chunk)
            sys.stdout.buffer.flush()

    os.close(master_fd)
    return process.wait()


def write_selection_manifest(
    selection: list[dict[str, object]],
) -> None:
    fieldnames = [
        "calculation_type",
        "frame",
        "cluster",
        "job",
        "root",
        "is_tracked_root",
        "pair_rank",
        "pair_label",
        "orbital_role",
        "orbital_index",
        "orbital_spin",
        "nto_file",
        "cube_file",
    ]

    with SELECTION_MANIFEST.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(selection)


def write_summary(
    rows: list[dict[str, object]],
) -> None:
    fieldnames = [
        "calculation_type",
        "frame",
        "cluster",
        "job",
        "root",
        "is_tracked_root",
        "pair_rank",
        "pair_label",
        "orbital_role",
        "orbital_index",
        "orbital_spin",
        "status",
        "returncode",
        "runtime_s",
        "cube_size_bytes",
        "cube_sha256",
        "validation",
        "nto_file",
        "cube_file",
        "transcript",
    ]

    with SUMMARY_CSV.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(rows)


def write_audit(
    rows: list[dict[str, object]],
) -> None:
    successful = sum(
        row["validation"] == "valid"
        for row in rows
    )

    skipped = sum(
        row["status"] == "already_valid"
        for row in rows
    )

    generated = sum(
        row["status"] == "generated"
        for row in rows
    )

    by_job: dict[str, int] = {}

    for row in rows:
        if row["validation"] == "valid":
            job = str(row["job"])
            by_job[job] = by_job.get(job, 0) + 1

    with AUDIT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day019 selected NTO cube-generation audit\n\n"
        )

        handle.write("## Scope\n\n")
        handle.write(
            "- Tracked root: dominant NTO pair only "
            "(hole and particle; two cubes per job).\n"
        )
        handle.write(
            "- Alternate root: two leading NTO pairs "
            "(four cubes per job).\n"
        )
        handle.write(
            "- Expected total: six cubes per job Ã eight jobs "
            "= 48 cubes.\n\n"
        )

        handle.write("## Validation\n\n")
        handle.write(
            f"- Valid selected cubes: {successful}/48\n"
        )
        handle.write(
            f"- Generated in this invocation: {generated}\n"
        )
        handle.write(
            f"- Reused as already valid: {skipped}\n"
        )
        handle.write(
            f"- Jobs with six valid cubes: "
            f"{sum(value == 6 for value in by_job.values())}/8\n"
        )
        handle.write(
            f"- Grid points per direction: {GRID_POINTS}\n\n"
        )

        handle.write("## Per-job cube counts\n\n")
        handle.write("| Job | Valid cubes |\n")
        handle.write("|---|---:|\n")

        for job in sorted(by_job):
            handle.write(
                f"| `{job}` | {by_job[job]} |\n"
            )

        handle.write("\n## Acceptance\n\n")

        if (
            successful == 48
            and len(by_job) == 8
            and all(value == 6 for value in by_job.values())
        ):
            handle.write(
                "**Day019 selected NTO cube generation: PASS.**\n"
            )
        else:
            handle.write(
                "**Day019 selected NTO cube generation: FAIL.**\n"
            )


def main() -> None:
    args = parse_args()
    executable = resolve_executable(args.orca_plot)

    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    transcript_root = OUTPUT_ROOT / "transcripts"
    transcript_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    metrics_rows = read_state_metrics()
    selection = build_selection(metrics_rows)

    if args.only_job:
        selection = [
            row
            for row in selection
            if row["job"] == args.only_job
        ]

        if not selection:
            raise SystemExit(
                f"Unknown --only-job: {args.only_job}"
            )

    write_selection_manifest(selection)

    log("Day019 selected NTO cube generator")
    log(f"orca_plot executable: {executable}")
    log(f"Cube selections: {len(selection)}")
    log(
        "Output is streamed live to the terminal; use shell tee "
        "to capture the complete batch log."
    )

    summary_rows: list[dict[str, object]] = []

    for index, row in enumerate(selection, start=1):
        nto_file = (
            PROJECT_ROOT / str(row["nto_file"])
        ).resolve()

        cube_file = (
            PROJECT_ROOT / str(row["cube_file"])
        ).resolve()

        validate_nto_file(nto_file)

        valid, validation = validate_cube(cube_file)

        transcript = (
            transcript_root
            / (
                f"{row['job']}.s{row['root']}."
                f"mo{row['orbital_index']}{row['orbital_spin']}.log"
            )
        )

        if valid:
            log(
                f"[{index}/{len(selection)}] SKIP "
                f"{cube_file.name}: already valid"
            )

            summary_rows.append(
                {
                    **row,
                    "status": "already_valid",
                    "returncode": 0,
                    "runtime_s": 0.0,
                    "cube_size_bytes": cube_file.stat().st_size,
                    "cube_sha256": sha256_file(cube_file),
                    "validation": "valid",
                    "transcript": (
                        str(transcript.relative_to(PROJECT_ROOT))
                        if transcript.exists()
                        else ""
                    ),
                }
            )

            write_summary(summary_rows)
            continue

        if cube_file.exists():
            if not args.rerun_invalid:
                raise SystemExit(
                    f"Invalid pre-existing cube: {cube_file} "
                    f"({validation}). Use --rerun-invalid."
                )

            archived = archive_invalid_cube(cube_file)
            log(
                f"[{index}/{len(selection)}] ARCHIVE "
                f"{archived.name}"
            )

        log("")
        log("=" * 88)
        log(
            f"[{index}/{len(selection)}] GENERATE "
            f"{row['job']} S{row['root']} "
            f"{row['orbital_role']} "
            f"MO {row['orbital_index']}{row['orbital_spin']}"
        )
        log(
            f"NTO: {nto_file.relative_to(PROJECT_ROOT)}"
        )
        log(
            f"Cube: {cube_file.relative_to(PROJECT_ROOT)}"
        )
        log("=" * 88)

        start = time.perf_counter()

        returncode = run_orca_plot_live(
            executable=executable,
            nto_file=nto_file,
            orbital=int(row["orbital_index"]),
            spin=str(row["orbital_spin"]),
            transcript=transcript,
        )

        runtime_s = time.perf_counter() - start
        valid, validation = validate_cube(cube_file)

        log("")
        log(
            f"DONE {cube_file.name} "
            f"returncode={returncode} "
            f"valid={valid} "
            f"runtime_s={runtime_s:.2f}"
        )

        summary_rows.append(
            {
                **row,
                "status": "generated" if valid else "failed",
                "returncode": returncode,
                "runtime_s": runtime_s,
                "cube_size_bytes": (
                    cube_file.stat().st_size
                    if cube_file.exists()
                    else 0
                ),
                "cube_sha256": (
                    sha256_file(cube_file)
                    if valid
                    else ""
                ),
                "validation": (
                    "valid"
                    if valid
                    else validation
                ),
                "transcript": str(
                    transcript.relative_to(PROJECT_ROOT)
                ),
            }
        )

        write_summary(summary_rows)
        write_audit(summary_rows)

        if returncode != 0 or not valid:
            raise SystemExit(
                f"Cube generation failed: {cube_file}\n"
                f"returncode={returncode}, validation={validation}\n"
                f"Transcript: {transcript}"
            )

    write_summary(summary_rows)
    write_audit(summary_rows)

    valid_count = sum(
        row["validation"] == "valid"
        for row in summary_rows
    )

    log("")
    log("Day019 selected NTO cube generation completed.")
    log(
        f"Valid selected cubes: "
        f"{valid_count}/{len(summary_rows)}"
    )
    log(
        f"Wrote: {SELECTION_MANIFEST.relative_to(PROJECT_ROOT)}"
    )
    log(
        f"Wrote: {SUMMARY_CSV.relative_to(PROJECT_ROOT)}"
    )
    log(
        f"Wrote: {AUDIT_MD.relative_to(PROJECT_ROOT)}"
    )

    expected = 48 if not args.only_job else len(selection)

    if valid_count != expected:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
