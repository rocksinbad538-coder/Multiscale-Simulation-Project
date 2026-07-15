#!/usr/bin/env python3
"""
Guarded execution of QM_F06 ORCA Stage 1.

The script refuses to execute unless --execute is supplied explicitly.

For each selected fragment it:
- validates ORCA availability;
- validates the Stage-1 input;
- creates a timestamped execution directory;
- copies the exact input used;
- records hashes and environment information;
- runs ORCA;
- captures stdout/stderr;
- checks normal termination, SCF convergence and geometry convergence;
- locates the optimized XYZ geometry;
- updates workflow_state.json only after successful validation.

By default, this script performs preflight only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

WORKFLOW_DIR = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/"
    "QM_F06/orca_workflow"
)

EXECUTION_ROOT = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/"
    "QM_F06/orca_executions"
)

ORCA_DEFAULT = Path(
    "/Users/alejandro/projects/"
    "orca_6_1_1_macosx_intel_openmpi411/orca"
)

VALID_FRAGMENTS = {
    "LOWER": "QM_F06_LOWER_CAPPED_REPAIRED",
    "UPPER": "QM_F06_UPPER_CAPPED_REPAIRED",
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty required file: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def parse_input(path: Path) -> dict[str, Any]:
    require_file(path)
    text = path.read_text(encoding="utf-8")

    state_match = re.search(
        r"(?m)^\s*\*\s+xyz\s+(-?\d+)\s+(\d+)\s*$",
        text,
    )

    if not state_match:
        raise RuntimeError(
            f"Charge/multiplicity block not found in {path}"
        )

    xyz_match = re.search(
        r"(?ms)^\s*\*\s+xyz\s+(-?\d+)\s+(\d+)\s*$"
        r"(.*?)"
        r"^\s*\*\s*$",
        text,
    )

    if not xyz_match:
        raise RuntimeError(f"XYZ block not found in {path}")

    atom_lines = [
        line.strip()
        for line in xyz_match.group(3).splitlines()
        if line.strip()
    ]

    constraints = re.findall(
        r"(?m)^\s*\{\s*C\s+(\d+)\s+C\s*\}\s*$",
        text,
    )

    return {
        "charge": int(state_match.group(1)),
        "multiplicity": int(state_match.group(2)),
        "atom_count": len(atom_lines),
        "constraint_indices": [int(value) for value in constraints],
        "contains_opt": bool(
            re.search(r"(?i)(^|\s)Opt(\s|$)", text)
        ),
        "contains_tightscf": "TightSCF" in text,
    }


def find_optimized_xyz(
    execution_dir: Path,
) -> Path | None:
    preferred_names = (
        "stage1.xyz",
        "stage1_opt.xyz",
        "stage1_trj.xyz",
    )

    for name in preferred_names:
        candidate = execution_dir / name
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate

    candidates = sorted(
        path
        for path in execution_dir.glob("*.xyz")
        if path.is_file() and path.stat().st_size > 0
    )

    return candidates[0] if candidates else None


def validate_output(
    output_path: Path,
    execution_dir: Path,
) -> dict[str, Any]:
    require_file(output_path)
    text = output_path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    normal_termination = (
        "ORCA TERMINATED NORMALLY" in text
    )

    scf_converged = (
        "SCF CONVERGED AFTER" in text
        or "SCF CONVERGED" in text
    )

    geometry_converged = (
        "THE OPTIMIZATION HAS CONVERGED" in text
        or "OPTIMIZATION CONVERGED" in text
    )

    error_markers = [
        marker
        for marker in (
            "ORCA finished by error termination",
            "ORCA TERMINATED ABNORMALLY",
            "SCF NOT CONVERGED",
            "The optimization did not converge",
        )
        if marker in text
    ]

    optimized_xyz = find_optimized_xyz(execution_dir)

    return {
        "normal_termination": normal_termination,
        "scf_converged": scf_converged,
        "geometry_converged": geometry_converged,
        "error_markers": error_markers,
        "optimized_xyz": (
            str(optimized_xyz)
            if optimized_xyz is not None
            else None
        ),
        "stage1_gate_pass": all(
            (
                normal_termination,
                scf_converged,
                geometry_converged,
                not error_markers,
                optimized_xyz is not None,
            )
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fragment",
        choices=("LOWER", "UPPER", "BOTH"),
        default="BOTH",
    )
    parser.add_argument(
        "--orca",
        type=Path,
        default=ORCA_DEFAULT,
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Explicitly authorize and execute ORCA Stage 1.",
    )

    args = parser.parse_args()

    selected = (
        list(VALID_FRAGMENTS.items())
        if args.fragment == "BOTH"
        else [
            (
                args.fragment,
                VALID_FRAGMENTS[args.fragment],
            )
        ]
    )

    orca_path = args.orca.expanduser().resolve()

    if not orca_path.is_file():
        detected = shutil.which("orca")

        if detected:
            orca_path = Path(detected).resolve()
        else:
            raise RuntimeError(
                f"ORCA executable not found: {orca_path}"
            )

    preflight_rows: list[dict[str, Any]] = []

    for end, fragment in selected:
        fragment_dir = WORKFLOW_DIR / fragment
        input_path = fragment_dir / "stage1.inp"
        state_path = fragment_dir / "workflow_state.json"

        require_file(input_path)
        require_file(state_path)

        parsed = parse_input(input_path)

        state = json.loads(
            state_path.read_text(encoding="utf-8")
        )

        preflight_pass = all(
            (
                parsed["atom_count"] == 22,
                parsed["charge"] == 0,
                parsed["multiplicity"] == 1,
                len(parsed["constraint_indices"]) == 11,
                len(set(parsed["constraint_indices"])) == 11,
                all(
                    0 <= index < 22
                    for index in parsed["constraint_indices"]
                ),
                parsed["contains_opt"],
                parsed["contains_tightscf"],
                state["stage1_input_prepared"] is True,
                state["stage1_executed"] is False,
            )
        )

        preflight_rows.append(
            {
                "end": end,
                "fragment": fragment,
                "input": str(input_path.relative_to(ROOT)),
                "input_sha256": sha256(input_path),
                "atom_count": parsed["atom_count"],
                "charge": parsed["charge"],
                "multiplicity": parsed["multiplicity"],
                "constraint_count": len(
                    parsed["constraint_indices"]
                ),
                "preflight_pass": preflight_pass,
            }
        )

    all_preflight_pass = all(
        row["preflight_pass"]
        for row in preflight_rows
    )

    print("=" * 72)
    print("QM_F06 ORCA STAGE-1 PREFLIGHT")
    print("=" * 72)
    print("ORCA:", orca_path)
    print("ORCA SHA256:", sha256(orca_path))
    print("Execution requested:", args.execute)
    print()

    for row in preflight_rows:
        print(row)

    print()
    print(
        "Preflight decision:",
        (
            "PASS"
            if all_preflight_pass
            else "FAIL"
        ),
    )

    if not all_preflight_pass:
        raise RuntimeError(
            "Stage-1 preflight failed. Execution blocked."
        )

    if not args.execute:
        print()
        print(
            "DRY PREFLIGHT COMPLETE — ORCA WAS NOT EXECUTED."
        )
        print(
            "Use --execute only after explicit authorization."
        )
        return

    timestamp = datetime.now().strftime("%Y%m%dT%H%M%S")

    for end, fragment in selected:
        fragment_dir = WORKFLOW_DIR / fragment
        source_input = fragment_dir / "stage1.inp"
        state_path = fragment_dir / "workflow_state.json"

        execution_dir = (
            EXECUTION_ROOT
            / fragment
            / f"stage1_{timestamp}"
        )
        execution_dir.mkdir(
            parents=True,
            exist_ok=False,
        )

        input_path = execution_dir / "stage1.inp"
        shutil.copy2(source_input, input_path)

        output_path = execution_dir / "stage1.out"
        stderr_path = execution_dir / "stage1.stderr"

        metadata = {
            "fragment": fragment,
            "end": end,
            "stage": 1,
            "timestamp": timestamp,
            "orca_path": str(orca_path),
            "orca_sha256": sha256(orca_path),
            "input_sha256": sha256(input_path),
            "python": sys.version,
            "platform": platform.platform(),
            "cwd": str(execution_dir),
            "environment_path": os.environ.get("PATH", ""),
            "explicit_execute_flag": True,
        }

        (
            execution_dir / "execution_metadata.json"
        ).write_text(
            json.dumps(metadata, indent=2) + "\n",
            encoding="utf-8",
        )

        print()
        print("=" * 72)
        print("EXECUTING:", fragment, "STAGE 1")
        print("Directory:", execution_dir)
        print("=" * 72)

        with (
            output_path.open("w", encoding="utf-8") as stdout,
            stderr_path.open("w", encoding="utf-8") as stderr,
        ):
            completed = subprocess.run(
                [str(orca_path), input_path.name],
                cwd=execution_dir,
                stdout=stdout,
                stderr=stderr,
                check=False,
            )

        validation = validate_output(
            output_path,
            execution_dir,
        )
        validation["return_code"] = completed.returncode

        (
            execution_dir / "stage1_validation.json"
        ).write_text(
            json.dumps(validation, indent=2) + "\n",
            encoding="utf-8",
        )

        print("Return code:", completed.returncode)
        print(
            "Normal termination:",
            validation["normal_termination"],
        )
        print(
            "SCF converged:",
            validation["scf_converged"],
        )
        print(
            "Geometry converged:",
            validation["geometry_converged"],
        )
        print(
            "Optimized XYZ:",
            validation["optimized_xyz"],
        )
        print(
            "Stage-1 gate:",
            (
                "PASS"
                if validation["stage1_gate_pass"]
                else "FAIL"
            ),
        )

        if (
            completed.returncode != 0
            or not validation["stage1_gate_pass"]
        ):
            raise RuntimeError(
                f"{fragment}: Stage-1 execution failed validation. "
                f"Inspect {execution_dir}"
            )

        state = json.loads(
            state_path.read_text(encoding="utf-8")
        )

        state["stage1_executed"] = True
        state["stage1_execution_directory"] = str(
            execution_dir.relative_to(ROOT)
        )
        state["stage1_optimized_xyz"] = str(
            Path(validation["optimized_xyz"]).relative_to(ROOT)
        )
        state["stage1_validation_pass"] = True
        state["qm_execution_authorized"] = False

        state_path.write_text(
            json.dumps(state, indent=2) + "\n",
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
