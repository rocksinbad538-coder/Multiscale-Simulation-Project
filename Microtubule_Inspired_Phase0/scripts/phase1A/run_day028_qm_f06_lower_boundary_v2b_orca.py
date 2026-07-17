#!/usr/bin/env python3
"""
Guarded ORCA execution for QM_F06 LOWER Boundary V2-B.

V2-B starts from the converged V2-A geometry and relaxes:
- all 15 hydrogen atoms;
- the three restored real heavy atoms.

The ten validated bridge/scaffold heavy atoms remain fixed.

The script:
- validates the input and workflow state;
- performs a dry preflight unless --execute is supplied;
- runs ORCA in a timestamped directory;
- validates normal termination, SCF convergence, geometry convergence
  and final XYZ generation;
- updates v2_workflow_state.json only after a fully validated run.

No ESP/RESP or force-field parameterization is executed.
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
    "runs/phase1A/day027_qm_f06_lower_boundary_redesign/"
    "QM_F06_LOWER_BOUNDARY_V2/orca_v2_workflow"
)

INPUT_PATH = WORKFLOW_DIR / "v2b_partial_relax.inp"
STATE_PATH = WORKFLOW_DIR / "v2_workflow_state.json"

EXECUTION_ROOT = ROOT / (
    "runs/phase1A/day028_qm_f06_lower_boundary_v2b_executions"
)

DEFAULT_ORCA = Path(
    "/Users/alejandro/projects/"
    "orca_6_1_1_macosx_intel_openmpi411/orca"
)


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def parse_input(path: Path) -> dict[str, Any]:
    require_file(path)

    text = path.read_text(encoding="utf-8")

    match = re.search(
        r"(?ms)^\s*\*\s+xyz\s+(-?\d+)\s+(\d+)\s*$"
        r"(.*?)"
        r"^\s*\*\s*$",
        text,
    )

    if not match:
        raise RuntimeError("XYZ block not found.")

    atom_lines = [
        line.split()
        for line in match.group(3).splitlines()
        if line.strip()
    ]

    constraints = [
        int(value)
        for value in re.findall(
            r"(?m)^\s*\{\s*C\s+(\d+)\s+C\s*\}\s*$",
            text,
        )
    ]

    checks = {
        "atom_count_28": len(atom_lines) == 28,
        "charge_zero": int(match.group(1)) == 0,
        "multiplicity_one": int(match.group(2)) == 1,
        "constraint_count_10": len(constraints) == 10,
        "constraint_indices_unique": (
            len(constraints) == len(set(constraints))
        ),
        "constraint_indices_in_range": all(
            0 <= index < 28
            for index in constraints
        ),
        "contains_pbe0": "PBE0" in text,
        "contains_d4": "D4" in text,
        "contains_def2_tzvp": "def2-TZVP" in text,
        "contains_def2_j": "def2/J" in text,
        "contains_rijcosx": "RIJCOSX" in text,
        "contains_tightscf": "TightSCF" in text,
        "contains_opt": bool(
            re.search(r"(?i)(^|\s)Opt(\s|$)", text)
        ),
        "contains_defgrid3": "DefGrid3" in text,
        "no_obsolete_grid": (
            "Grid5" not in text
            and "FinalGrid6" not in text
        ),
    }

    return {
        "atom_count": len(atom_lines),
        "charge": int(match.group(1)),
        "multiplicity": int(match.group(2)),
        "constraints": constraints,
        "checks": checks,
        "gate_pass": all(checks.values()),
    }


def locate_final_xyz(
    execution_dir: Path,
) -> Path | None:
    preferred = (
        execution_dir / "v2b_partial_relax.xyz",
        execution_dir / "v2b.xyz",
    )

    for path in preferred:
        if path.is_file() and path.stat().st_size > 0:
            return path

    candidates = sorted(
        path
        for path in execution_dir.glob("*.xyz")
        if path.is_file()
        and path.stat().st_size > 0
        and "_trj" not in path.name
    )

    return candidates[0] if candidates else None


def validate_output(
    output_path: Path,
    execution_dir: Path,
    return_code: int,
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
            "INPUT ERROR",
        )
        if marker in text
    ]

    final_xyz = locate_final_xyz(execution_dir)

    gate_pass = all(
        (
            return_code == 0,
            normal_termination,
            scf_converged,
            geometry_converged,
            not error_markers,
            final_xyz is not None,
        )
    )

    return {
        "return_code": return_code,
        "normal_termination": normal_termination,
        "scf_converged": scf_converged,
        "geometry_converged": geometry_converged,
        "error_markers": error_markers,
        "optimized_xyz": (
            str(final_xyz)
            if final_xyz is not None
            else None
        ),
        "v2b_gate_pass": gate_pass,
    }


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--orca",
        type=Path,
        default=DEFAULT_ORCA,
    )

    parser.add_argument(
        "--execute",
        action="store_true",
    )

    args = parser.parse_args()

    require_file(INPUT_PATH)
    require_file(STATE_PATH)

    orca_path = args.orca.expanduser().resolve()

    if not orca_path.is_file():
        detected = shutil.which("orca")

        if detected:
            orca_path = Path(detected).resolve()
        else:
            raise RuntimeError(
                f"ORCA executable not found: {orca_path}"
            )

    parsed = parse_input(INPUT_PATH)

    state = json.loads(
        STATE_PATH.read_text(encoding="utf-8")
    )

    state_checks = {
        "v2a_executed": (
            state.get("v2a_executed") is True
        ),
        "v2a_validation_pass": (
            state.get("v2a_validation_pass") is True
        ),
        "v2a_geometry_promoted": (
            state.get("v2a_geometry_promoted") is True
        ),
        "v2b_input_generated": (
            state.get("v2b_input_generated") is True
        ),
        "v2b_fixed_count_10": (
            state.get("v2b_fixed_atom_count") == 10
        ),
        "v2b_mobile_count_18": (
            state.get("v2b_mobile_atom_count") == 18
        ),
        "v2b_not_executed": (
            state.get("v2b_executed") is False
        ),
    }

    constraint_map_match = (
        parsed["constraints"]
        == state.get("v2b_fixed_indices")
    )

    state_checks["constraint_map_match"] = (
        constraint_map_match
    )

    preflight_pass = (
        parsed["gate_pass"]
        and all(state_checks.values())
    )

    print("=" * 78)
    print("QM_F06 LOWER BOUNDARY V2-B ORCA PREFLIGHT")
    print("=" * 78)
    print("ORCA:", orca_path)
    print("ORCA SHA256:", sha256(orca_path))
    print("Input:", INPUT_PATH)
    print("Input SHA256:", sha256(INPUT_PATH))
    print("Execution requested:", args.execute)
    print()

    print("Input checks:")
    for key, value in parsed["checks"].items():
        print(f"  {key:36s}: {value}")

    print()
    print("Workflow-state checks:")
    for key, value in state_checks.items():
        print(f"  {key:36s}: {value}")

    print()
    print(
        "Preflight decision:",
        "PASS" if preflight_pass else "FAIL",
    )

    if not preflight_pass:
        raise RuntimeError(
            "V2-B preflight failed. Execution blocked."
        )

    if not args.execute:
        print()
        print(
            "DRY PREFLIGHT COMPLETE — ORCA WAS NOT EXECUTED."
        )
        return

    timestamp = datetime.now().strftime("%Y%m%dT%H%M%S")

    execution_dir = (
        EXECUTION_ROOT / f"v2b_{timestamp}"
    )
    execution_dir.mkdir(
        parents=True,
        exist_ok=False,
    )

    execution_input = (
        execution_dir / "v2b_partial_relax.inp"
    )
    shutil.copy2(INPUT_PATH, execution_input)

    output_path = (
        execution_dir / "v2b_partial_relax.out"
    )
    stderr_path = (
        execution_dir / "v2b_partial_relax.stderr"
    )

    metadata = {
        "fragment": "QM_F06_LOWER_BOUNDARY_V2_REPAIRED",
        "stage": "V2-B",
        "timestamp": timestamp,
        "orca_path": str(orca_path),
        "orca_sha256": sha256(orca_path),
        "input_sha256": sha256(execution_input),
        "python": sys.version,
        "platform": platform.platform(),
        "working_directory": str(execution_dir),
        "environment_path": os.environ.get("PATH", ""),
        "explicit_execute_flag": True,
        "fixed_atom_count": 10,
        "mobile_atom_count": 18,
    }

    (
        execution_dir / "execution_metadata.json"
    ).write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )

    print()
    print("=" * 78)
    print("EXECUTING QM_F06 LOWER BOUNDARY V2-B")
    print("=" * 78)
    print("Execution directory:", execution_dir)

    with (
        output_path.open("w", encoding="utf-8") as stdout,
        stderr_path.open("w", encoding="utf-8") as stderr,
    ):
        completed = subprocess.run(
            [
                str(orca_path),
                execution_input.name,
            ],
            cwd=execution_dir,
            stdout=stdout,
            stderr=stderr,
            check=False,
        )

    validation = validate_output(
        output_path,
        execution_dir,
        completed.returncode,
    )

    (
        execution_dir / "v2b_validation.json"
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
        "V2-B gate:",
        (
            "PASS"
            if validation["v2b_gate_pass"]
            else "FAIL"
        ),
    )

    if not validation["v2b_gate_pass"]:
        raise RuntimeError(
            "V2-B execution failed validation. "
            f"Inspect {execution_dir}"
        )

    state = json.loads(
        STATE_PATH.read_text(encoding="utf-8")
    )

    optimized_xyz = Path(
        validation["optimized_xyz"]
    )

    state["v2b_executed"] = True
    state["v2b_validation_pass"] = True
    state["v2b_execution_directory"] = str(
        execution_dir.relative_to(ROOT)
    )
    state["v2b_optimized_xyz"] = str(
        optimized_xyz.relative_to(ROOT)
    )
    state["qm_execution_authorized"] = False

    STATE_PATH.write_text(
        json.dumps(state, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
