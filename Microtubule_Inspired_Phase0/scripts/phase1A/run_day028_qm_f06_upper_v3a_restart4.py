#!/usr/bin/env python3
"""
Guarded runner for the clean four-process QM_F06 UPPER V3-A restart.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

WORKFLOW_DIR = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V3/orca_v3_workflow/restart4"
)

INPUT_PATH = WORKFLOW_DIR / "v3a_restart4.inp"
SUMMARY_PATH = WORKFLOW_DIR / "v3a_restart4_summary.json"

STATE_PATH = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V3/orca_v3_workflow/"
    "v3_workflow_state.json"
)

EXECUTION_ROOT = ROOT / (
    "runs/phase1A/"
    "day028_qm_f06_upper_boundary_v3a_restart4_executions"
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

    for path in (
        INPUT_PATH,
        SUMMARY_PATH,
        STATE_PATH,
    ):
        require_file(path)

    orca = args.orca.expanduser().resolve()

    if not orca.is_file():
        raise RuntimeError(f"ORCA not found: {orca}")

    text = INPUT_PATH.read_text(encoding="utf-8")
    summary = json.loads(
        SUMMARY_PATH.read_text(encoding="utf-8")
    )
    state = json.loads(
        STATE_PATH.read_text(encoding="utf-8")
    )

    checks = {
        "nprocs_4": bool(
            re.search(r"(?i)\bnprocs\s+4\b", text)
        ),
        "maxcore_2500": "%maxcore 2500" in text,
        "contains_opt": bool(
            re.search(r"(?i)(^|\s)Opt(\s|$)", text)
        ),
        "contains_pbe0_d4": "PBE0 D4" in text,
        "contains_def2_tzvp": "def2-TZVP" in text,
        "contains_rijcosx": "RIJCOSX" in text,
        "fresh_scf_guess": (
            summary["fresh_scf_guess"] is True
        ),
        "gbw_not_reused": (
            summary["gbw_reused"] is False
        ),
        "restart_prepared": (
            state["v3a_restart4_prepared"] is True
        ),
        "restart_not_executed": (
            state["v3a_restart4_executed"] is False
        ),
        "input_hash_matches": (
            summary["restart_input_sha256"]
            == sha256(INPUT_PATH)
            == state["v3a_restart4_input_sha256"]
        ),
    }

    preflight_pass = all(checks.values())

    print("=" * 78)
    print("QM_F06 UPPER V3-A RESTART4 PREFLIGHT")
    print("=" * 78)
    print("ORCA:", orca)
    print("Input:", INPUT_PATH)
    print("Input SHA256:", sha256(INPUT_PATH))
    print("Execution requested:", args.execute)
    print()

    for key, value in checks.items():
        print(f"{key:38s}: {value}")

    print()
    print(
        "Preflight decision:",
        "PASS" if preflight_pass else "FAIL",
    )

    if not preflight_pass:
        raise RuntimeError("Restart4 preflight failed.")

    if not args.execute:
        print()
        print(
            "DRY PREFLIGHT COMPLETE — ORCA WAS NOT EXECUTED."
        )
        return

    timestamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    execution_dir = (
        EXECUTION_ROOT / f"restart4_{timestamp}"
    )
    execution_dir.mkdir(parents=True, exist_ok=False)

    execution_input = execution_dir / "restart4.inp"
    output_path = execution_dir / "restart4.out"
    stderr_path = execution_dir / "restart4.stderr"

    shutil.copy2(INPUT_PATH, execution_input)

    with (
        output_path.open("w", encoding="utf-8") as stdout,
        stderr_path.open("w", encoding="utf-8") as stderr,
    ):
        completed = subprocess.run(
            [str(orca), execution_input.name],
            cwd=execution_dir,
            stdout=stdout,
            stderr=stderr,
            check=False,
        )

    output_text = output_path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    energies = re.findall(
        r"FINAL SINGLE POINT ENERGY\s+"
        r"(-?\d+\.\d+)",
        output_text,
    )

    optimized_xyz = execution_dir / "restart4.xyz"

    validation = {
        "return_code": completed.returncode,
        "normal_termination": (
            "ORCA TERMINATED NORMALLY" in output_text
        ),
        "scf_converged": (
            "SCF CONVERGED AFTER" in output_text
        ),
        "geometry_converged": (
            "THE OPTIMIZATION HAS CONVERGED" in output_text
        ),
        "final_energy_hartree": (
            float(energies[-1])
            if energies
            else None
        ),
        "optimized_xyz": (
            str(optimized_xyz)
            if optimized_xyz.is_file()
            else None
        ),
    }

    validation["execution_gate_pass"] = all(
        (
            validation["return_code"] == 0,
            validation["normal_termination"],
            validation["scf_converged"],
            validation["geometry_converged"],
            validation["final_energy_hartree"]
            is not None,
            optimized_xyz.is_file(),
        )
    )

    (
        execution_dir / "restart4_validation.json"
    ).write_text(
        json.dumps(validation, indent=2) + "\n",
        encoding="utf-8",
    )

    print()
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
        "Execution gate:",
        (
            "PASS"
            if validation["execution_gate_pass"]
            else "FAIL"
        ),
    )

    if not validation["execution_gate_pass"]:
        raise RuntimeError(
            f"Restart4 failed. Inspect {execution_dir}"
        )

    state["v3a_restart4_executed"] = True
    state["v3a_restart4_validation_pass"] = True
    state["v3a_restart4_execution_directory"] = str(
        execution_dir.relative_to(ROOT)
    )
    state["v3a_restart4_optimized_xyz"] = str(
        optimized_xyz.relative_to(ROOT)
    )
    state["v3a_restart4_final_energy_hartree"] = (
        validation["final_energy_hartree"]
    )
    state["qm_execution_authorized"] = False

    STATE_PATH.write_text(
        json.dumps(state, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
