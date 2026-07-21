#!/usr/bin/env python3
"""
Prepare an immutable execution directory for QM_F06 UPPER V4.

The preflight:
- requires the accepted V4 input-preparation report;
- creates a timestamped execution directory;
- copies the audited input, XYZ and constraint map;
- verifies hashes after copying;
- reparses all critical ORCA parameters;
- confirms a fresh SCF and the expected constraint indices;
- verifies the ORCA executable;
- creates a supervised launch script;
- does not start ORCA.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import stat
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

SOURCE_DIR = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_orca_input"
)

SOURCE_INPUT = SOURCE_DIR / "v4.inp"
SOURCE_XYZ = SOURCE_DIR / "v4_start.xyz"

SOURCE_MAP = (
    SOURCE_DIR
    / "QM_F06_UPPER_V4_constraint_map.csv"
)

SOURCE_REPORT = (
    SOURCE_DIR
    / "QM_F06_UPPER_V4_ORCA_INPUT_PREPARATION.json"
)

EXECUTION_PARENT = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_executions"
)

LATEST_POINTER = (
    EXECUTION_PARENT
    / "LATEST_V4_EXECUTION.txt"
)

ORCA_EXECUTABLE = Path(
    "/Users/alejandro/projects/"
    "orca_6_1_1_macosx_intel_openmpi411/orca"
)

EXPECTED_DECISION = (
    "QM_F06_UPPER_V4_ORCA_INPUT_PREPARED_"
    "EXECUTION_PREFLIGHT_REQUIRED"
)

EXPECTED_KEYWORD_LINE = (
    "! PBE0 D4 def2-TZVP def2/J RIJCOSX "
    "TightSCF Opt DefGrid3"
)

EXPECTED_ATOM_COUNT = 46
EXPECTED_CHARGE = 0
EXPECTED_MULTIPLICITY = 1
EXPECTED_NPROCS = 4
EXPECTED_MAXCORE = 2500
EXPECTED_SCF_MAXITER = 500
EXPECTED_GEOM_MAXITER = 250

EXPECTED_FIXED_INDICES = (
    list(range(14)) + [16]
)


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty file: {path}"
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {
        "true",
        "1",
        "yes",
    }


def parse_input(path: Path) -> dict:
    require_file(path)

    text = path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    non_comment_text = "\n".join(
        line
        for line in text.splitlines()
        if not line.lstrip().startswith("#")
    )

    keyword_lines = [
        line.strip()
        for line in text.splitlines()
        if line.lstrip().startswith("!")
    ]

    nprocs_match = re.search(
        r"\bnprocs\s+(\d+)",
        text,
        re.IGNORECASE,
    )

    maxcore_match = re.search(
        r"%maxcore\s+(\d+)",
        text,
        re.IGNORECASE,
    )

    scf_maxiter_match = re.search(
        r"%scf\b.*?\bMaxIter\s+(\d+).*?\bend\b",
        text,
        re.IGNORECASE | re.DOTALL,
    )

    geom_maxiter_match = re.search(
        r"%geom\b.*?\bMaxIter\s+(\d+)",
        text,
        re.IGNORECASE | re.DOTALL,
    )

    xyz_header = re.search(
        r"^\s*\*\s*xyz\s+([-+]?\d+)\s+(\d+)\s*$",
        text,
        re.IGNORECASE | re.MULTILINE,
    )

    if not all((
        nprocs_match,
        maxcore_match,
        scf_maxiter_match,
        geom_maxiter_match,
        xyz_header,
    )):
        raise RuntimeError(
            "Could not parse all required ORCA fields."
        )

    fixed_indices = sorted({
        int(value)
        for value in re.findall(
            r"\{\s*C\s+(\d+)\s+C\s*\}",
            text,
            re.IGNORECASE,
        )
    })

    xyz_block = re.search(
        r"^\s*\*\s*xyz\s+[-+]?\d+\s+\d+\s*$"
        r"(.*?)"
        r"^\s*\*\s*$",
        text,
        re.IGNORECASE
        | re.MULTILINE
        | re.DOTALL,
    )

    if xyz_block is None:
        raise RuntimeError(
            "Could not locate inline XYZ block."
        )

    coordinate_lines = [
        line
        for line in xyz_block.group(1).splitlines()
        if line.strip()
        and not line.lstrip().startswith("#")
    ]

    reuse_hits = {
        "moread": bool(re.search(
            r"\bmoread\b",
            non_comment_text,
            re.IGNORECASE,
        )),
        "moinp": bool(re.search(
            r"%moinp\b",
            non_comment_text,
            re.IGNORECASE,
        )),
        "gbw_reference": bool(re.search(
            r"\S+\.gbw",
            non_comment_text,
            re.IGNORECASE,
        )),
    }

    return {
        "text": text,
        "keyword_lines": keyword_lines,
        "nprocs": int(nprocs_match.group(1)),
        "maxcore": int(maxcore_match.group(1)),
        "scf_maxiter": int(
            scf_maxiter_match.group(1)
        ),
        "geom_maxiter": int(
            geom_maxiter_match.group(1)
        ),
        "charge": int(xyz_header.group(1)),
        "multiplicity": int(
            xyz_header.group(2)
        ),
        "fixed_indices": fixed_indices,
        "coordinate_count": len(
            coordinate_lines
        ),
        "reuse_hits": reuse_hits,
    }


def main() -> None:
    for path in (
        SOURCE_INPUT,
        SOURCE_XYZ,
        SOURCE_MAP,
        SOURCE_REPORT,
    ):
        require_file(path)

    report = json.loads(
        SOURCE_REPORT.read_text(
            encoding="utf-8"
        )
    )

    if report["decision"] != EXPECTED_DECISION:
        raise RuntimeError(
            "Unexpected input-preparation decision: "
            f"{report['decision']}"
        )

    if not report["overall_pass"]:
        raise RuntimeError(
            "Input-preparation report did not pass."
        )

    if not report["authorization"][
        "execution_preflight_authorized"
    ]:
        raise RuntimeError(
            "Execution preflight is not authorized."
        )

    source_hashes = {
        "input": sha256(SOURCE_INPUT),
        "xyz": sha256(SOURCE_XYZ),
        "constraint_map": sha256(SOURCE_MAP),
        "input_preparation_report": sha256(
            SOURCE_REPORT
        ),
    }

    report_hashes = report["files_sha256"]

    source_hash_gate = (
        source_hashes["input"]
        == report_hashes["prepared_input"]
        and source_hashes["xyz"]
        == report_hashes["copied_start_xyz"]
        and source_hashes["constraint_map"]
        == report_hashes["copied_constraint_map"]
    )

    if not source_hash_gate:
        raise RuntimeError(
            "One or more prepared source artifacts "
            "changed after input preparation."
        )

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%SZ")

    execution_dir = (
        EXECUTION_PARENT
        / f"v4_{timestamp}"
    )

    if execution_dir.exists():
        raise RuntimeError(
            f"Execution directory already exists: "
            f"{execution_dir}"
        )

    execution_dir.mkdir(
        parents=True,
        exist_ok=False,
    )

    copied_input = execution_dir / "v4.inp"
    copied_xyz = execution_dir / "v4_start.xyz"

    copied_map = (
        execution_dir
        / "QM_F06_UPPER_V4_constraint_map.csv"
    )

    copied_report = (
        execution_dir
        / "QM_F06_UPPER_V4_ORCA_INPUT_PREPARATION.json"
    )

    shutil.copy2(SOURCE_INPUT, copied_input)
    shutil.copy2(SOURCE_XYZ, copied_xyz)
    shutil.copy2(SOURCE_MAP, copied_map)
    shutil.copy2(
        SOURCE_REPORT,
        copied_report,
    )

    copied_hashes = {
        "input": sha256(copied_input),
        "xyz": sha256(copied_xyz),
        "constraint_map": sha256(copied_map),
        "input_preparation_report": sha256(
            copied_report
        ),
    }

    post_copy_checks = {
        key: copied_hashes[key]
        == source_hashes[key]
        for key in source_hashes
    }

    parsed = parse_input(copied_input)
    map_rows = read_csv(copied_map)

    map_fixed_indices = sorted(
        int(row["index_0based"])
        for row in map_rows
        if parse_bool(row["v4_fixed"])
    )

    map_mobile_indices = sorted(
        int(row["index_0based"])
        for row in map_rows
        if parse_bool(row["v4_mobile"])
    )

    orca_available = (
        ORCA_EXECUTABLE.is_file()
        and os.access(
            ORCA_EXECUTABLE,
            os.X_OK,
        )
    )

    gates = {
        "input_preparation_decision": (
            report["decision"]
            == EXPECTED_DECISION
        ),
        "input_preparation_overall_pass": (
            report["overall_pass"]
        ),
        "input_preparation_authorization": (
            report["authorization"][
                "execution_preflight_authorized"
            ]
        ),
        "source_hashes_match_report": (
            source_hash_gate
        ),
        "post_copy_input_hash": (
            post_copy_checks["input"]
        ),
        "post_copy_xyz_hash": (
            post_copy_checks["xyz"]
        ),
        "post_copy_constraint_map_hash": (
            post_copy_checks["constraint_map"]
        ),
        "post_copy_report_hash": (
            post_copy_checks[
                "input_preparation_report"
            ]
        ),
        "keyword_line": (
            parsed["keyword_lines"]
            == [EXPECTED_KEYWORD_LINE]
        ),
        "charge": (
            parsed["charge"]
            == EXPECTED_CHARGE
        ),
        "multiplicity": (
            parsed["multiplicity"]
            == EXPECTED_MULTIPLICITY
        ),
        "nprocs": (
            parsed["nprocs"]
            == EXPECTED_NPROCS
        ),
        "maxcore": (
            parsed["maxcore"]
            == EXPECTED_MAXCORE
        ),
        "scf_maxiter": (
            parsed["scf_maxiter"]
            == EXPECTED_SCF_MAXITER
        ),
        "geom_maxiter": (
            parsed["geom_maxiter"]
            == EXPECTED_GEOM_MAXITER
        ),
        "coordinate_count": (
            parsed["coordinate_count"]
            == EXPECTED_ATOM_COUNT
        ),
        "input_fixed_indices": (
            parsed["fixed_indices"]
            == EXPECTED_FIXED_INDICES
        ),
        "map_fixed_indices": (
            map_fixed_indices
            == EXPECTED_FIXED_INDICES
        ),
        "fixed_mobile_partition": (
            sorted(
                map_fixed_indices
                + map_mobile_indices
            )
            == list(range(EXPECTED_ATOM_COUNT))
            and not (
                set(map_fixed_indices)
                & set(map_mobile_indices)
            )
        ),
        "fresh_scf": (
            not any(
                parsed["reuse_hits"].values()
            )
        ),
        "orca_available": orca_available,
    }

    overall_pass = all(gates.values())

    if not overall_pass:
        failed = [
            key
            for key, value in gates.items()
            if not value
        ]

        raise RuntimeError(
            "Execution preflight failed: "
            f"{failed}"
        )

    launch_script = (
        execution_dir
        / "run_v4_supervised.sh"
    )

    launch_text = f"""\
#!/usr/bin/env bash
set -u
set -o pipefail

cd "$(dirname "$0")"

ORCA="{ORCA_EXECUTABLE}"

echo "$$" > v4.supervisor_pid

START_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf '%s\\n' "$START_UTC" > v4.start_utc

"$ORCA" v4.inp > v4.out 2> v4.stderr &
ORCA_PID="$!"

echo "$ORCA_PID" > v4.orca_pid

wait "$ORCA_PID"
STATUS="$?"

END_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

printf '%s\\n' "$STATUS" > v4.exit_status
printf '%s\\n' "$END_UTC" > v4.end_utc

exit "$STATUS"
"""

    launch_script.write_text(
        launch_text,
        encoding="utf-8",
    )

    current_mode = launch_script.stat().st_mode

    launch_script.chmod(
        current_mode
        | stat.S_IXUSR
        | stat.S_IXGRP
        | stat.S_IXOTH
    )

    manifest = {
        "decision": (
            "QM_F06_UPPER_V4_EXECUTION_"
            "PREFLIGHT_PASS"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "execution_directory": str(
            execution_dir.relative_to(ROOT)
        ),
        "source_directory": str(
            SOURCE_DIR.relative_to(ROOT)
        ),
        "orca_executable": str(
            ORCA_EXECUTABLE
        ),
        "gates": gates,
        "overall_pass": overall_pass,
        "post_copy_checks": post_copy_checks,
        "input_parameters": {
            "keyword_line": (
                parsed["keyword_lines"][0]
            ),
            "charge": parsed["charge"],
            "multiplicity": (
                parsed["multiplicity"]
            ),
            "nprocs": parsed["nprocs"],
            "maxcore_mb_per_process": (
                parsed["maxcore"]
            ),
            "nominal_total_maxcore_mb": (
                parsed["nprocs"]
                * parsed["maxcore"]
            ),
            "scf_maxiter": (
                parsed["scf_maxiter"]
            ),
            "geometry_maxiter": (
                parsed["geom_maxiter"]
            ),
            "coordinate_count": (
                parsed["coordinate_count"]
            ),
            "fixed_indices": (
                parsed["fixed_indices"]
            ),
            "mobile_indices": (
                map_mobile_indices
            ),
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
        / "QM_F06_UPPER_V4_EXECUTION_MANIFEST.json"
    )

    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )

    EXECUTION_PARENT.mkdir(
        parents=True,
        exist_ok=True,
    )

    LATEST_POINTER.write_text(
        str(execution_dir.relative_to(ROOT))
        + "\n",
        encoding="utf-8",
    )

    print("=" * 78)
    print("QM_F06 UPPER V4 EXECUTION PREFLIGHT")
    print("=" * 78)

    for gate, passed in gates.items():
        print(
            f"{gate:42s}: "
            f"{'PASS' if passed else 'FAIL'}"
        )

    print()
    print(
        "Charge / multiplicity:",
        parsed["charge"],
        parsed["multiplicity"],
    )
    print("nprocs:", parsed["nprocs"])
    print("maxcore:", parsed["maxcore"])
    print(
        "Coordinate count:",
        parsed["coordinate_count"],
    )
    print(
        "Fixed indices:",
        parsed["fixed_indices"],
    )
    print(
        "Forbidden reuse hits:",
        parsed["reuse_hits"],
    )
    print(
        "ORCA executable:",
        ORCA_EXECUTABLE,
    )

    print()
    print(
        "Decision: "
        "QM_F06_UPPER_V4_EXECUTION_PREFLIGHT_PASS"
    )
    print(
        "Execution directory:",
        execution_dir,
    )
    print("Manifest:", manifest_path)
    print("Launch script:", launch_script)
    print()
    print("ORCA execution prepared: True")
    print("ORCA execution started: False")


if __name__ == "__main__":
    main()
