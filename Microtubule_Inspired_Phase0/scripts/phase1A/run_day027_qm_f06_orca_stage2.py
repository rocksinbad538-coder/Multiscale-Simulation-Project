#!/usr/bin/env python3
"""
Guarded ORCA Stage-2 execution for QM_F06 LOWER/UPPER fragments.

Stage 2 is executed only when --execute is explicitly supplied.

Preflight requirements:
- Stage-1 execution and validation passed.
- Stage-1 geometry was promoted.
- Stage-2 input exists and was generated.
- 22 atoms, charge 0, multiplicity 1.
- Exactly four constraints: 0, 2, 11, 12.
- PBE0-D4/def2-TZVP, Opt, TightSCF and DefGrid3 present.
- Stage 2 has not already been executed.

After execution, the script validates:
- return code;
- normal ORCA termination;
- SCF convergence;
- geometry convergence;
- existence of the final optimized XYZ.

No Stage-3 input is generated automatically.
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
    "runs/phase1A/day027_qm_f06_orca_executions"
)

ORCA_DEFAULT = Path(
    "/Users/alejandro/projects/"
    "orca_6_1_1_macosx_intel_openmpi411/orca"
)

FRAGMENTS = {
    "LOWER": "QM_F06_LOWER_CAPPED_REPAIRED",
    "UPPER": "QM_F06_UPPER_CAPPED_REPAIRED",
}

EXPECTED_CONSTRAINTS = [0, 2, 11, 12]


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

    xyz_match = re.search(
        r"(?ms)^\s*\*\s+xyz\s+(-?\d+)\s+(\d+)\s*$"
        r"(.*?)"
        r"^\s*\*\s*$",
        text,
    )

    if not xyz_match:
        raise RuntimeError(f"XYZ block not found in {path}")

    atoms = [
        line.split()
        for line in xyz_match.group(3).splitlines()
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
        "atom_count_22": len(atoms) == 22,
        "charge_zero": int(xyz_match.group(1)) == 0,
        "multiplicity_one": int(xyz_match.group(2)) == 1,
        "constraints_expected": constraints == EXPECTED_CONSTRAINTS,
        "promoted_from_stage1": (
            "Coordinates promoted from Stage 1" in text
        ),
        "contains_pbe0": "PBE0" in text,
        "contains_d4": "D4" in text,
        "contains_def2_tzvp": "def2-TZVP" in text,
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
        "atom_count": len(atoms),
        "charge": int(xyz_match.group(1)),
        "multiplicity": int(xyz_match.group(2)),
        "constraints": constraints,
        "checks": checks,
        "gate_pass": all(checks.values()),
    }


def locate_optimized_xyz(
    execution_dir: Path,
) -> Path | None:
    preferred = (
        execution_dir / "stage2.xyz",
        execution_dir / "stage2_opt.xyz",
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

    optimized_xyz = locate_optimized_xyz(execution_dir)

    gate_pass = all(
        (
            return_code == 0,
            normal_termination,
            scf_converged,
            geometry_converged,
            not error_markers,
            optimized_xyz is not None,
        )
    )

    return {
        "return_code": return_code,
        "normal_termination": normal_termination,
        "scf_converged": scf_converged,
        "geometry_converged": geometry_converged,
        "error_markers": error_markers,
        "optimized_xyz": (
            str(optimized_xyz)
            if optimized_xyz is not None
            else None
        ),
        "stage2_gate_pass": gate_pass,
    }


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fragment",
        choices=("LOWER", "UPPER"),
        required=True,
    )

    parser.add_argument(
        "--orca",
        type=Path,
        default=ORCA_DEFAULT,
    )

    parser.add_argument(
        "--execute",
        action="store_true",
    )

    args = parser.parse_args()

    fragment = FRAGMENTS[args.fragment]
    fragment_dir = WORKFLOW_DIR / fragment
    input_path = fragment_dir / "stage2.inp"
    state_path = fragment_dir / "workflow_state.json"

    require_file(input_path)
    require_file(state_path)

    orca_path = args.orca.expanduser().resolve()

    if not orca_path.is_file():
        detected = shutil.which("orca")

        if detected:
            orca_path = Path(detected).resolve()
        else:
            raise RuntimeError(
                f"ORCA executable not found: {orca_path}"
            )

    parsed = parse_input(input_path)

    state = json.loads(
        state_path.read_text(encoding="utf-8")
    )

    state_checks = {
        "stage1_executed": (
            state.get("stage1_executed") is True
        ),
        "stage1_validation_pass": (
            state.get("stage1_validation_pass") is True
        ),
        "stage1_geometry_promoted": (
            state.get("stage1_geometry_promoted") is True
        ),
        "stage2_input_generated": (
            state.get("stage2_input_generated") is True
        ),
        "stage2_not_executed": (
            state.get("stage2_executed") is False
        ),
    }

    preflight_pass = (
        parsed["gate_pass"]
        and all(state_checks.values())
    )

    print("=" * 76)
    print("QM_F06 ORCA STAGE-2 PREFLIGHT")
    print("=" * 76)
    print("Fragment:", fragment)
    print("ORCA:", orca_path)
    print("ORCA SHA256:", sha256(orca_path))
    print("Input:", input_path)
    print("Input SHA256:", sha256(input_path))
    print("Execution requested:", args.execute)
    print()

    print("Input checks:")
    for key, value in parsed["checks"].items():
        print(f"  {key:34s}: {value}")

    print()
    print("Workflow-state checks:")
    for key, value in state_checks.items():
        print(f"  {key:34s}: {value}")

    print()
    print(
        "Preflight decision:",
        "PASS" if preflight_pass else "FAIL",
    )

    if not preflight_pass:
        raise RuntimeError(
            "Stage-2 preflight failed. Execution blocked."
        )

    if not args.execute:
        print()
        print(
            "DRY PREFLIGHT COMPLETE — ORCA WAS NOT EXECUTED."
        )
        return

    timestamp = datetime.now().strftime("%Y%m%dT%H%M%S")

    execution_dir = (
        EXECUTION_ROOT
        / fragment
        / f"stage2_{timestamp}"
    )
    execution_dir.mkdir(
        parents=True,
        exist_ok=False,
    )

    execution_input = execution_dir / "stage2.inp"
    shutil.copy2(input_path, execution_input)

    output_path = execution_dir / "stage2.out"
    stderr_path = execution_dir / "stage2.stderr"

    metadata = {
        "fragment": fragment,
        "stage": 2,
        "timestamp": timestamp,
        "orca_path": str(orca_path),
        "orca_sha256": sha256(orca_path),
        "input_sha256": sha256(execution_input),
        "python": sys.version,
        "platform": platform.platform(),
        "working_directory": str(execution_dir),
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
    print("=" * 76)
    print("EXECUTING QM_F06 STAGE 2")
    print("=" * 76)
    print("Fragment:", fragment)
    print("Execution directory:", execution_dir)

    with (
        output_path.open("w", encoding="utf-8") as stdout,
        stderr_path.open("w", encoding="utf-8") as stderr,
    ):
        completed = subprocess.run(
            [str(orca_path), execution_input.name],
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

    validation_path = (
        execution_dir / "stage2_validation.json"
    )
    validation_path.write_text(
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
        "Stage-2 gate:",
        (
            "PASS"
            if validation["stage2_gate_pass"]
            else "FAIL"
        ),
    )

    if not validation["stage2_gate_pass"]:
        raise RuntimeError(
            f"{fragment}: Stage-2 execution failed validation. "
            f"Inspect {execution_dir}"
        )

    state = json.loads(
        state_path.read_text(encoding="utf-8")
    )

    optimized_xyz = Path(
        validation["optimized_xyz"]
    )

    state["stage2_executed"] = True
    state["stage2_execution_directory"] = str(
        execution_dir.relative_to(ROOT)
    )
    state["stage2_optimized_xyz"] = str(
        optimized_xyz.relative_to(ROOT)
    )
    state["stage2_validation_pass"] = True
    state["qm_execution_authorized"] = False

    state_path.write_text(
        json.dumps(state, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
