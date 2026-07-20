#!/usr/bin/env python3
"""
Prepare a clean 4-process restart of QM_F06 UPPER Boundary V3-A.

The restart:
- uses the final complete trajectory frame from the failed 6-process run;
- preserves the original atom ordering and constraints;
- starts a fresh SCF rather than reusing the potentially inconsistent GBW;
- does not execute ORCA.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

EXECUTION_ROOT = ROOT / (
    "runs/phase1A/"
    "day028_qm_f06_upper_boundary_v3a_executions"
)

WORKFLOW_DIR = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V3/orca_v3_workflow"
)

SOURCE_INPUT = WORKFLOW_DIR / "v3a_boundary_relax.inp"
STATE_PATH = WORKFLOW_DIR / "v3_workflow_state.json"

RESTART_DIR = WORKFLOW_DIR / "restart4"
RESTART_XYZ = RESTART_DIR / "v3a_restart4_start.xyz"
RESTART_INPUT = RESTART_DIR / "v3a_restart4.inp"
RESTART_SUMMARY = RESTART_DIR / "v3a_restart4_summary.json"


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def read_last_xyz_frame(
    path: Path,
) -> tuple[str, list[str]]:
    require_file(path)

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    frames: list[tuple[str, list[str]]] = []
    index = 0

    while index < len(lines):
        stripped = lines[index].strip()

        if not stripped:
            index += 1
            continue

        try:
            atom_count = int(stripped)
        except ValueError:
            index += 1
            continue

        if index + 2 + atom_count > len(lines):
            break

        comment = lines[index + 1].strip()
        atom_lines = lines[index + 2:index + 2 + atom_count]

        valid = all(
            len(line.split()) >= 4
            for line in atom_lines
        )

        if valid:
            frames.append((comment, atom_lines))

        index += atom_count + 2

    if not frames:
        raise RuntimeError("No complete XYZ frames found.")

    return frames[-1]


def main() -> None:
    RESTART_DIR.mkdir(parents=True, exist_ok=True)

    require_file(SOURCE_INPUT)
    require_file(STATE_PATH)

    executions = sorted(
        path
        for path in EXECUTION_ROOT.glob("v3a_*")
        if path.is_dir()
    )

    if not executions:
        raise RuntimeError("No V3-A execution directory found.")

    failed_execution = executions[-1]
    trajectory = (
        failed_execution / "v3a_boundary_relax_trj.xyz"
    )

    comment, atom_lines = read_last_xyz_frame(trajectory)

    if len(atom_lines) != 30:
        raise RuntimeError(
            f"Expected 30 atoms; found {len(atom_lines)}"
        )

    RESTART_XYZ.write_text(
        "\n".join(
            [
                "30",
                (
                    "QM_F06_UPPER_V3A_RESTART4_START; "
                    f"source comment: {comment}"
                ),
                *atom_lines,
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    source_text = SOURCE_INPUT.read_text(
        encoding="utf-8"
    )

    xyz_match = re.search(
        r"(?ms)^(\s*\*\s+xyz\s+)(-?\d+)\s+(\d+)\s*$"
        r"(.*?)"
        r"^\s*\*\s*$",
        source_text,
    )

    if not xyz_match:
        raise RuntimeError(
            "XYZ block not found in source input."
        )

    coordinate_block = "\n".join(atom_lines)

    restart_text = (
        source_text[:xyz_match.start()]
        + "* xyz 0 1\n"
        + coordinate_block
        + "\n*\n"
        + source_text[xyz_match.end():]
    )

    restart_text = re.sub(
        r"(?im)^\s*nprocs\s+\d+\s*$",
        "  nprocs 4",
        restart_text,
        count=1,
    )

    restart_text = re.sub(
        r"(?im)^%maxcore\s+\d+\s*$",
        "%maxcore 2500",
        restart_text,
        count=1,
    )

    restart_text = restart_text.replace(
        "# QM_F06 UPPER Boundary V3-A constrained boundary relaxation",
        (
            "# QM_F06 UPPER Boundary V3-A clean restart "
            "from final complete 6-process trajectory frame"
        ),
    )

    RESTART_INPUT.write_text(
        restart_text,
        encoding="utf-8",
    )

    state = json.loads(
        STATE_PATH.read_text(encoding="utf-8")
    )

    summary = {
        "decision": (
            "QM_F06_UPPER_V3A_RESTART4_PREPARED_"
            "EXECUTION_NOT_AUTHORIZED"
        ),
        "failed_execution_directory": str(
            failed_execution.relative_to(ROOT)
        ),
        "source_trajectory": str(
            trajectory.relative_to(ROOT)
        ),
        "source_frame_comment": comment,
        "source_frame_atom_count": len(atom_lines),
        "restart_xyz": str(
            RESTART_XYZ.relative_to(ROOT)
        ),
        "restart_input": str(
            RESTART_INPUT.relative_to(ROOT)
        ),
        "restart_input_sha256": sha256(RESTART_INPUT),
        "nprocs": 4,
        "maxcore_mb_per_process": 2500,
        "nominal_total_maxcore_mb": 10000,
        "fresh_scf_guess": True,
        "gbw_reused": False,
        "constraints_preserved": True,
        "restart_execution_authorized": False,
    }

    RESTART_SUMMARY.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    state["v3a_6proc_failed"] = True
    state["v3a_6proc_failure_module"] = "LEANSCF_MPI"
    state["v3a_6proc_final_geometry_accepted"] = False
    state["v3a_restart4_prepared"] = True
    state["v3a_restart4_input"] = str(
        RESTART_INPUT.relative_to(ROOT)
    )
    state["v3a_restart4_input_sha256"] = sha256(
        RESTART_INPUT
    )
    state["v3a_restart4_executed"] = False
    state["v3a_restart4_validation_pass"] = False
    state["qm_execution_authorized"] = False

    STATE_PATH.write_text(
        json.dumps(state, indent=2) + "\n",
        encoding="utf-8",
    )

    print("QM_F06 UPPER V3-A restart4 prepared.")
    print("Failed execution:", failed_execution)
    print("Source frame:", comment)
    print("Atoms:", len(atom_lines))
    print("Restart XYZ:", RESTART_XYZ)
    print("Restart input:", RESTART_INPUT)
    print("nprocs: 4")
    print("MaxCore: 2500 MB/process")
    print("GBW reused: False")
    print("Input SHA256:", sha256(RESTART_INPUT))
    print("Execution authorized: False")


if __name__ == "__main__":
    main()
