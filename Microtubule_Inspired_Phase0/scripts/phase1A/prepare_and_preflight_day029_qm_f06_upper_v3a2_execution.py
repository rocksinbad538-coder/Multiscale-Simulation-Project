#!/usr/bin/env python3
"""
Prepare a timestamped ORCA execution directory for QM_F06 UPPER V3-A2.

Requirements:
- the final V3-A2 pre-QM report must authorize ORCA execution;
- the source input and geometry hashes must match the audited hashes;
- charge, multiplicity, resources and fresh-SCF state are rechecked;
- no ORCA calculation is executed by this script.

A timestamped execution directory, manifest and launch command are generated.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

WORKFLOW = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V3/orca_v3a2_workflow"
)

SOURCE_INPUT = WORKFLOW / "v3a2.inp"
SOURCE_XYZ = WORKFLOW / "v3a2_start.xyz"
SOURCE_MAP = WORKFLOW / "v3a2_atom_role_constraint_map.csv"

AUDIT_REPORT = (
    WORKFLOW
    / "pre_qm_audit"
    / "QM_F06_UPPER_V3A2_PRE_QM_AUDIT.json"
)

EXECUTION_PARENT = ROOT / (
    "runs/phase1A/"
    "day029_qm_f06_upper_v3a2_executions"
)

EXPECTED_DECISION = (
    "QM_F06_UPPER_V3A2_PRE_QM_GATE_PASS_"
    "ORCA_EXECUTION_AUTHORIZED"
)

EXPECTED_ATOMS = 30
EXPECTED_CHARGE = 0
EXPECTED_MULTIPLICITY = 1
EXPECTED_NPROCS = 4
EXPECTED_MAXCORE = 2500
EXPECTED_FIXED = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
    10, 11, 12, 13, 15, 16, 18, 19,
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def extract_xyz_header(text: str) -> tuple[int, int]:
    match = re.search(
        r"(?im)^\s*\*\s+xyz\s+(-?\d+)\s+(\d+)\s*$",
        text,
    )

    if match is None:
        raise RuntimeError(
            "Could not find ORCA '* xyz charge multiplicity' header."
        )

    return int(match.group(1)), int(match.group(2))


def extract_integer(
    text: str,
    pattern: str,
    label: str,
) -> int:
    match = re.search(pattern, text)

    if match is None:
        raise RuntimeError(f"Could not find {label}.")

    return int(match.group(1))


def extract_constraints(text: str) -> set[int]:
    return {
        int(value)
        for value in re.findall(
            r"\{\s*C\s+(\d+)\s+C\s*\}",
            text,
        )
    }


def count_xyz_atoms_in_input(text: str) -> int:
    block_match = re.search(
        r"(?ms)^\s*\*\s+xyz\s+-?\d+\s+\d+\s*$"
        r"(.*?)"
        r"^\s*\*\s*$",
        text,
    )

    if block_match is None:
        raise RuntimeError(
            "Could not identify the ORCA XYZ coordinate block."
        )

    atom_lines = re.findall(
        r"(?m)^\s*(?:B|N|H)\s+"
        r"[-+]?(?:\d+\.\d+|\d+)(?:[Ee][-+]?\d+)?\s+"
        r"[-+]?(?:\d+\.\d+|\d+)(?:[Ee][-+]?\d+)?\s+"
        r"[-+]?(?:\d+\.\d+|\d+)(?:[Ee][-+]?\d+)?\s*$",
        block_match.group(1),
    )

    return len(atom_lines)


def command_version(command: str) -> str | None:
    executable = shutil.which(command)

    if executable is None:
        return None

    try:
        result = subprocess.run(
            [executable],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=20,
            check=False,
        )
    except Exception:
        return executable

    output = result.stdout.strip()

    return output[:1000] if output else executable


def main() -> None:
    for path in (
        SOURCE_INPUT,
        SOURCE_XYZ,
        SOURCE_MAP,
        AUDIT_REPORT,
    ):
        require_file(path)

    audit = json.loads(
        AUDIT_REPORT.read_text(encoding="utf-8")
    )

    input_text = SOURCE_INPUT.read_text(
        encoding="utf-8",
        errors="replace",
    )

    checks = {}

    checks["audit_decision"] = (
        audit.get("decision") == EXPECTED_DECISION
    )

    checks["audit_overall_pass"] = (
        audit.get("overall_pass") is True
    )

    checks["audit_authorization"] = (
        audit.get("authorization", {}).get(
            "orca_execution_authorized"
        )
        is True
    )

    current_hashes = {
        "input": sha256(SOURCE_INPUT),
        "xyz": sha256(SOURCE_XYZ),
        "constraint_map": sha256(SOURCE_MAP),
    }

    audited_hashes = audit["files_sha256"]

    checks["input_hash_matches_audit"] = (
        current_hashes["input"]
        == audited_hashes["input"]
    )

    checks["xyz_hash_matches_audit"] = (
        current_hashes["xyz"]
        == audited_hashes["xyz"]
    )

    checks["constraint_map_hash_matches_audit"] = (
        current_hashes["constraint_map"]
        == audited_hashes["constraint_map"]
    )

    charge, multiplicity = extract_xyz_header(input_text)

    checks["charge"] = (
        charge == EXPECTED_CHARGE
    )

    checks["multiplicity"] = (
        multiplicity == EXPECTED_MULTIPLICITY
    )

    nprocs = extract_integer(
        input_text,
        r"(?im)^\s*nprocs\s+(\d+)\s*$",
        "nprocs",
    )

    maxcore = extract_integer(
        input_text,
        r"(?im)^\s*%maxcore\s+(\d+)\s*$",
        "%maxcore",
    )

    checks["nprocs"] = (
        nprocs == EXPECTED_NPROCS
    )

    checks["maxcore"] = (
        maxcore == EXPECTED_MAXCORE
    )

    fixed_indices = extract_constraints(input_text)

    checks["fixed_indices"] = (
        fixed_indices == EXPECTED_FIXED
    )

    coordinate_count = count_xyz_atoms_in_input(input_text)

    checks["coordinate_count"] = (
        coordinate_count == EXPECTED_ATOMS
    )

    forbidden_patterns = {
        "moread": r"(?i)\bmoread\b",
        "moinp": r"(?i)\bmoinp\b",
        "gbw_reference": r"(?i)\.gbw\b",
    }

    forbidden_hits = {
        name: bool(re.search(pattern, input_text))
        for name, pattern in forbidden_patterns.items()
    }

    checks["fresh_scf"] = not any(
        forbidden_hits.values()
    )

    orca_path = shutil.which("orca")
    checks["orca_available"] = (
        orca_path is not None
    )

    overall_pass = all(checks.values())

    print("=" * 78)
    print("QM_F06 UPPER V3-A2 EXECUTION PREFLIGHT")
    print("=" * 78)

    for name, passed in checks.items():
        print(
            f"{name:40s}: "
            f"{'PASS' if passed else 'FAIL'}"
        )

    print()
    print("Charge / multiplicity:", charge, multiplicity)
    print("nprocs:", nprocs)
    print("maxcore:", maxcore)
    print("Coordinate count:", coordinate_count)
    print("Fixed indices:", sorted(fixed_indices))
    print("Forbidden reuse hits:", forbidden_hits)
    print("ORCA executable:", orca_path)
    print()

    if not overall_pass:
        print(
            "Decision: "
            "QM_F06_UPPER_V3A2_EXECUTION_PREFLIGHT_FAIL"
        )
        print("ORCA execution prepared: False")
        raise SystemExit(1)

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%SZ")

    execution_dir = (
        EXECUTION_PARENT
        / f"v3a2_{timestamp}"
    )

    if execution_dir.exists():
        raise RuntimeError(
            f"Execution directory already exists: {execution_dir}"
        )

    execution_dir.mkdir(parents=True)

    execution_input = execution_dir / "v3a2.inp"
    execution_xyz = execution_dir / "v3a2_start.xyz"
    execution_map = (
        execution_dir
        / "v3a2_atom_role_constraint_map.csv"
    )
    execution_audit = (
        execution_dir
        / "QM_F06_UPPER_V3A2_PRE_QM_AUDIT.json"
    )

    shutil.copy2(SOURCE_INPUT, execution_input)
    shutil.copy2(SOURCE_XYZ, execution_xyz)
    shutil.copy2(SOURCE_MAP, execution_map)
    shutil.copy2(AUDIT_REPORT, execution_audit)

    copied_hashes = {
        "input": sha256(execution_input),
        "xyz": sha256(execution_xyz),
        "constraint_map": sha256(execution_map),
        "audit_report": sha256(execution_audit),
    }

    post_copy_checks = {
        "input": (
            copied_hashes["input"]
            == current_hashes["input"]
        ),
        "xyz": (
            copied_hashes["xyz"]
            == current_hashes["xyz"]
        ),
        "constraint_map": (
            copied_hashes["constraint_map"]
            == current_hashes["constraint_map"]
        ),
    }

    if not all(post_copy_checks.values()):
        raise RuntimeError(
            "Hash mismatch after copying execution files."
        )

    launch_script = execution_dir / "launch_v3a2.sh"

    launch_script.write_text(
        """#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

if pgrep -af "[o]rca.*v3a2.inp" >/dev/null 2>&1; then
    echo "Another v3a2 ORCA process appears to be active."
    exit 1
fi

ORCA_BIN="$(command -v orca)"

if [ -z "${ORCA_BIN:-}" ]; then
    echo "ORCA executable not found."
    exit 1
fi

echo "ORCA executable: $ORCA_BIN"
echo "Execution directory: $PWD"
echo "Start UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"

"$ORCA_BIN" v3a2.inp > v3a2.out 2> v3a2.stderr

status=$?

echo "$status" > v3a2.exit_status

echo "End UTC: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "Exit status: $status"

exit "$status"
""",
        encoding="utf-8",
    )

    launch_script.chmod(0o755)

    manifest = {
        "decision": (
            "QM_F06_UPPER_V3A2_EXECUTION_PREPARED"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "execution_directory": str(
            execution_dir.relative_to(ROOT)
        ),
        "preflight_overall_pass": True,
        "checks": checks,
        "post_copy_checks": post_copy_checks,
        "orca_executable": orca_path,
        "orca_probe": command_version("orca"),
        "environment": {
            "PATH": os.environ.get("PATH"),
            "OMP_NUM_THREADS": os.environ.get(
                "OMP_NUM_THREADS"
            ),
            "MKL_NUM_THREADS": os.environ.get(
                "MKL_NUM_THREADS"
            ),
        },
        "input_parameters": {
            "charge": charge,
            "multiplicity": multiplicity,
            "nprocs": nprocs,
            "maxcore_mb_per_process": maxcore,
            "coordinate_count": coordinate_count,
            "fixed_indices": sorted(fixed_indices),
            "fresh_scf": True,
        },
        "files_sha256": copied_hashes,
        "launch_script": str(
            launch_script.relative_to(ROOT)
        ),
        "authorization": {
            "orca_execution_authorized": True,
            "geometry_reference_accepted": False,
            "electronic_reference_accepted": False,
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    manifest_path = (
        execution_dir
        / "QM_F06_UPPER_V3A2_EXECUTION_MANIFEST.json"
    )

    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )

    pointer_path = (
        EXECUTION_PARENT
        / "LATEST_V3A2_EXECUTION.txt"
    )

    pointer_path.write_text(
        str(execution_dir.relative_to(ROOT)) + "\n",
        encoding="utf-8",
    )

    print(
        "Decision: "
        "QM_F06_UPPER_V3A2_EXECUTION_PREFLIGHT_PASS"
    )
    print("Execution directory:", execution_dir)
    print("Manifest:", manifest_path)
    print("Launch script:", launch_script)
    print()
    print("ORCA execution prepared: True")
    print("ORCA execution started: False")


if __name__ == "__main__":
    main()
