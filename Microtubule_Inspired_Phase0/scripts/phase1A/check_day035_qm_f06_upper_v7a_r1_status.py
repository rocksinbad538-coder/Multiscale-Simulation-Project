#!/usr/bin/env python3
"""
Reproducible status check for the QM_F06_UPPER_V7A_R1 ORCA optimization.

This script is observational only:
- it does not modify the ORCA execution;
- it does not terminate processes;
- it does not promote coordinates;
- it does not execute downstream post-QM gates.

Exit codes
----------
0 : optimization converged and ORCA terminated normally
1 : optimization still running or status incomplete
2 : ORCA terminated without verified geometry convergence
3 : explicit error/failure marker detected
4 : required execution artifacts are missing or invalid
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]

EXECUTIONS_DIR = (
    REPO_ROOT
    / "runs"
    / "phase1A"
    / "day035_qm_f06_upper_v7a_r1_executions"
)

LATEST_POINTER = EXECUTIONS_DIR / "LATEST_V7A_R1_EXECUTION.txt"

REPORT_NAME = "QM_F06_UPPER_V7A_R1_LIVE_STATUS.json"

CYCLE_RE = re.compile(
    r"GEOMETRY OPTIMIZATION CYCLE\s+(\d+)",
    re.IGNORECASE,
)

CONVERGENCE_HEADER = "Geometry convergence"

SUCCESS_MARKERS = {
    "optimization_converged": "THE OPTIMIZATION HAS CONVERGED",
    "orca_terminated_normally": "ORCA TERMINATED NORMALLY",
}

FAILURE_PATTERNS = {
    "maximum_cycles_reached": re.compile(
        r"Maximum number of optimization cycles",
        re.IGNORECASE,
    ),
    "orca_error": re.compile(
        r"\bORCA\b.*\bERROR\b|\bERROR\b.*\bORCA\b",
        re.IGNORECASE,
    ),
    "abort": re.compile(r"\bABORT(?:ED)?\b", re.IGNORECASE),
    "fatal": re.compile(r"\bFATAL\b", re.IGNORECASE),
    "segmentation_fault": re.compile(
        r"segmentation fault|segfault",
        re.IGNORECASE,
    ),
    "scf_maxiter": re.compile(
        r"Please increase MaxIter|SCF.*not converged",
        re.IGNORECASE,
    ),
}

CONVERGENCE_ROW_RE = re.compile(
    r"^\s*"
    r"(Energy change|RMS gradient|MAX gradient|RMS step|MAX step)"
    r"\s+"
    r"([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)"
    r"\s+"
    r"([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)"
    r"\s+"
    r"(YES|NO)"
    r"\s*$",
    re.IGNORECASE,
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def resolve_execution_dir() -> Path:
    if not LATEST_POINTER.is_file():
        raise FileNotFoundError(f"Missing latest-execution pointer: {LATEST_POINTER}")

    raw = LATEST_POINTER.read_text(encoding="utf-8").strip()
    if not raw:
        raise ValueError(f"Empty latest-execution pointer: {LATEST_POINTER}")

    candidate = Path(raw)
    if not candidate.is_absolute():
        candidate = REPO_ROOT / candidate

    candidate = candidate.resolve()

    if not candidate.is_dir():
        raise NotADirectoryError(f"Execution directory not found: {candidate}")

    return candidate


def read_pid(path: Path) -> int | None:
    if not path.is_file():
        return None

    raw = path.read_text(encoding="utf-8").strip()
    if not raw.isdigit():
        return None

    return int(raw)


def process_is_active(pid: int | None) -> bool:
    if pid is None:
        return False

    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    else:
        return True


def process_snapshot(pid: int | None) -> dict[str, Any] | None:
    if pid is None:
        return None

    command = [
        "ps",
        "-p",
        str(pid),
        "-o",
        "pid=,ppid=,etime=,%cpu=,%mem=,state=,command=",
    ]

    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
    )

    output = completed.stdout.strip()
    if completed.returncode != 0 or not output:
        return None

    return {
        "raw": output,
    }


def complete_xyz_frames(path: Path) -> dict[str, Any]:
    result: dict[str, Any] = {
        "path": str(path),
        "exists": path.is_file(),
        "atoms_per_frame": None,
        "total_lines": None,
        "complete_frames": None,
        "remainder_lines": None,
        "valid": False,
    }

    if not path.is_file():
        return result

    lines = read_text(path).splitlines()
    result["total_lines"] = len(lines)

    if not lines:
        return result

    try:
        atoms = int(lines[0].strip())
    except ValueError:
        return result

    frame_lines = atoms + 2
    if atoms <= 0 or frame_lines <= 2:
        return result

    result["atoms_per_frame"] = atoms
    result["complete_frames"] = len(lines) // frame_lines
    result["remainder_lines"] = len(lines) % frame_lines
    result["valid"] = (len(lines) % frame_lines) == 0

    return result


def extract_cycles(text: str) -> list[int]:
    return [int(match.group(1)) for match in CYCLE_RE.finditer(text)]


def extract_last_convergence_table(text: str) -> dict[str, Any] | None:
    lines = text.splitlines()

    header_indices = [
        index
        for index, line in enumerate(lines)
        if CONVERGENCE_HEADER.lower() in line.lower()
    ]

    if not header_indices:
        return None

    start = header_indices[-1]
    window = lines[start : start + 15]

    criteria: dict[str, Any] = {}

    for line in window:
        match = CONVERGENCE_ROW_RE.match(line)
        if not match:
            continue

        label, value, tolerance, converged = match.groups()
        key = label.lower().replace(" ", "_")

        criteria[key] = {
            "value": float(value),
            "tolerance": float(tolerance),
            "converged": converged.upper() == "YES",
        }

    return {
        "start_line_1_based": start + 1,
        "criteria": criteria,
        "all_criteria_converged": (
            len(criteria) == 5
            and all(item["converged"] for item in criteria.values())
        ),
        "raw_lines": window,
    }


def detect_failure_markers(text: str) -> dict[str, bool]:
    return {
        name: bool(pattern.search(text))
        for name, pattern in FAILURE_PATTERNS.items()
    }


def classify_status(
    *,
    process_active: bool,
    optimization_converged: bool,
    terminated_normally: bool,
    failures: dict[str, bool],
) -> tuple[str, str, int]:
    explicit_failure = any(failures.values())

    if explicit_failure:
        return (
            "FAILED",
            "Diagnose the explicit ORCA failure before any downstream gate.",
            3,
        )

    if optimization_converged and terminated_normally:
        return (
            "CONVERGED_AND_TERMINATED",
            "The post-QM audit is eligible to run.",
            0,
        )

    if process_active:
        return (
            "RUNNING",
            "Leave ORCA running and repeat this status check later.",
            1,
        )

    if terminated_normally and not optimization_converged:
        return (
            "TERMINATED_WITHOUT_VERIFIED_OPTIMIZATION_CONVERGENCE",
            "Do not run post-QM gates; inspect the optimization outcome.",
            2,
        )

    return (
        "INCOMPLETE_OR_UNVERIFIED",
        "Do not run post-QM gates; inspect process and output state.",
        1,
    )


def file_metadata(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {
            "path": str(path),
            "exists": False,
        }

    stat_result = path.stat()

    return {
        "path": str(path),
        "exists": True,
        "size_bytes": stat_result.st_size,
        "mtime_utc": datetime.fromtimestamp(
            stat_result.st_mtime,
            timezone.utc,
        ).isoformat(),
    }


def main() -> int:
    try:
        execution_dir = resolve_execution_dir()
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 4

    output_path = execution_dir / "v7a_r1.out"
    trajectory_path = execution_dir / "v7a_r1_trj.xyz"
    pid_path = execution_dir / "v7a_r1.orca_pid"
    report_path = execution_dir / REPORT_NAME

    if not output_path.is_file():
        print(f"ERROR: Missing ORCA output: {output_path}", file=sys.stderr)
        return 4

    output_text = read_text(output_path)
    pid = read_pid(pid_path)
    active = process_is_active(pid)

    cycles = extract_cycles(output_text)
    last_cycle = cycles[-1] if cycles else None

    optimization_converged = (
        SUCCESS_MARKERS["optimization_converged"] in output_text
    )
    terminated_normally = (
        SUCCESS_MARKERS["orca_terminated_normally"] in output_text
    )

    failures = detect_failure_markers(output_text)
    convergence = extract_last_convergence_table(output_text)
    trajectory = complete_xyz_frames(trajectory_path)

    status, allowed_action, exit_code = classify_status(
        process_active=active,
        optimization_converged=optimization_converged,
        terminated_normally=terminated_normally,
        failures=failures,
    )

    report: dict[str, Any] = {
        "schema_version": 1,
        "generated_at_utc": utc_now(),
        "system_id": "QM_F06_UPPER_V7A_R1",
        "execution_directory": str(execution_dir),
        "status": status,
        "allowed_action": allowed_action,
        "process": {
            "pid_file": str(pid_path),
            "pid": pid,
            "active": active,
            "snapshot": process_snapshot(pid),
        },
        "markers": {
            "optimization_converged": optimization_converged,
            "orca_terminated_normally": terminated_normally,
            **failures,
        },
        "optimization": {
            "cycles_detected": len(cycles),
            "last_cycle": last_cycle,
            "last_convergence_table": convergence,
        },
        "trajectory": trajectory,
        "files": {
            "output": file_metadata(output_path),
            "trajectory": file_metadata(trajectory_path),
            "gbw": file_metadata(execution_dir / "v7a_r1.gbw"),
            "final_xyz": file_metadata(execution_dir / "v7a_r1.xyz"),
        },
    }

    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print("=" * 72)
    print("QM_F06_UPPER_V7A_R1 ORCA STATUS")
    print("=" * 72)
    print(f"execution_dir             : {execution_dir}")
    print(f"status                    : {status}")
    print(f"orca_pid                  : {pid}")
    print(f"process_active            : {active}")
    print(f"last_geometry_cycle       : {last_cycle}")
    print(
        "complete_trajectory_frames:"
        f" {trajectory.get('complete_frames')}"
    )
    print(f"optimization_converged    : {optimization_converged}")
    print(f"orca_terminated_normally  : {terminated_normally}")

    active_failures = [
        name for name, detected in failures.items() if detected
    ]
    print(
        "failure_markers           : "
        + (", ".join(active_failures) if active_failures else "none")
    )

    if convergence is not None:
        print("last convergence criteria :")
        for name, values in convergence["criteria"].items():
            print(
                f"  {name:16s} "
                f"value={values['value']:.10g} "
                f"tolerance={values['tolerance']:.10g} "
                f"converged={values['converged']}"
            )
    else:
        print("last convergence criteria : unavailable")

    print(f"allowed_action            : {allowed_action}")
    print(f"report                    : {report_path}")
    print("=" * 72)

    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
