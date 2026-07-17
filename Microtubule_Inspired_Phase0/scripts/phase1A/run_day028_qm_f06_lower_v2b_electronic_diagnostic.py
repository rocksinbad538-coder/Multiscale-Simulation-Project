#!/usr/bin/env python3
"""
Guarded ORCA execution for the QM_F06 LOWER V2-B electronic diagnostic.

The calculation is a single point on the accepted V2-B geometry and requests:
- Mayer bond orders;
- Hirshfeld charges;
- MBIS charges;
- CHELPG charges.

No geometry optimization, RESP fitting, force-field fitting or parameter
adoption is performed.
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

DIAGNOSTIC_DIR = ROOT / (
    "runs/phase1A/"
    "day028_qm_f06_lower_boundary_v2b_electronic_diagnostic"
)

INPUT_PATH = DIAGNOSTIC_DIR / (
    "QM_F06_LOWER_BOUNDARY_V2B_ELECTRONIC_DIAGNOSTIC.inp"
)

SUMMARY_PATH = DIAGNOSTIC_DIR / (
    "QM_F06_LOWER_BOUNDARY_V2B_"
    "electronic_diagnostic_summary.json"
)

STATE_PATH = ROOT / (
    "runs/phase1A/day027_qm_f06_lower_boundary_redesign/"
    "QM_F06_LOWER_BOUNDARY_V2/orca_v2_workflow/"
    "v2_workflow_state.json"
)

STRUCTURAL_SUMMARY_PATH = ROOT / (
    "runs/phase1A/"
    "day028_qm_f06_lower_boundary_v2b_postprocessing/"
    "QM_F06_LOWER_BOUNDARY_V2B_validation_summary.json"
)

CONTACT_SUMMARY_PATH = ROOT / (
    "runs/phase1A/"
    "day028_qm_f06_lower_boundary_v2b_postprocessing/"
    "residual_real_contact_audit/"
    "residual_real_contact_audit_summary.json"
)

EXECUTION_ROOT = DIAGNOSTIC_DIR / "orca_executions"

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

    checks = {
        "atom_count_28": len(atom_lines) == 28,
        "charge_zero": int(match.group(1)) == 0,
        "multiplicity_one": int(match.group(2)) == 1,
        "contains_pbe0": "PBE0" in text,
        "contains_d4": "D4" in text,
        "contains_def2_tzvp": "def2-TZVP" in text,
        "contains_def2_j": "def2/J" in text,
        "contains_rijcosx": "RIJCOSX" in text,
        "contains_tightscf": "TightSCF" in text,
        "contains_defgrid3": "DefGrid3" in text,
        "contains_mayer": "MAYER" in text,
        "contains_hirshfeld": "HIRSHFELD" in text,
        "contains_mbis": "MBIS" in text,
        "contains_chelpg": "CHELPG" in text,
        "mayer_threshold_present": (
            "MAYER_BONDORDERTHRESH 0.01" in text
        ),
        "no_geometry_optimization": not bool(
            re.search(r"(?i)(^|\s)Opt(\s|$)", text)
        ),
        "no_frequency_calculation": not bool(
            re.search(r"(?i)(^|\s)Freq(\s|$)", text)
        ),
    }

    return {
        "atom_count": len(atom_lines),
        "charge": int(match.group(1)),
        "multiplicity": int(match.group(2)),
        "checks": checks,
        "gate_pass": all(checks.values()),
    }


def validate_output(
    output_path: Path,
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

    final_energy_matches = re.findall(
        r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)",
        text,
    )

    final_energy = (
        float(final_energy_matches[-1])
        if final_energy_matches
        else None
    )

    error_markers = [
        marker
        for marker in (
            "ORCA finished by error termination",
            "ORCA TERMINATED ABNORMALLY",
            "SCF NOT CONVERGED",
            "INPUT ERROR",
            "UNKNOWN KEYWORD",
            "ERROR IN INPUT",
        )
        if marker in text
    ]

    section_checks = {
        "mayer_output_detected": any(
            marker in text
            for marker in (
                "MAYER BOND ORDERS",
                "Mayer bond orders",
                "MAYER POPULATION ANALYSIS",
            )
        ),
        "hirshfeld_output_detected": any(
            marker in text
            for marker in (
                "HIRSHFELD ANALYSIS",
                "HIRSHFELD ATOMIC CHARGES",
                "Hirshfeld charges",
            )
        ),
        "mbis_output_detected": any(
            marker in text
            for marker in (
                "MBIS ANALYSIS",
                "MBIS ATOMIC CHARGES",
                "Minimal Basis Iterative Stockholder",
            )
        ),
        "chelpg_output_detected": any(
            marker in text
            for marker in (
                "CHELPG Charges",
                "CHELPG CHARGES",
                "CHELPG",
            )
        ),
    }

    execution_gate_pass = all(
        (
            return_code == 0,
            normal_termination,
            scf_converged,
            final_energy is not None,
            not error_markers,
        )
    )

    complete_diagnostic_output = (
        execution_gate_pass
        and all(section_checks.values())
    )

    return {
        "return_code": return_code,
        "normal_termination": normal_termination,
        "scf_converged": scf_converged,
        "final_single_point_energy_hartree": final_energy,
        "error_markers": error_markers,
        **section_checks,
        "execution_gate_pass": execution_gate_pass,
        "complete_diagnostic_output": complete_diagnostic_output,
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

    for path in (
        INPUT_PATH,
        SUMMARY_PATH,
        STATE_PATH,
        STRUCTURAL_SUMMARY_PATH,
        CONTACT_SUMMARY_PATH,
    ):
        require_file(path)

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

    preparation_summary = json.loads(
        SUMMARY_PATH.read_text(encoding="utf-8")
    )

    workflow_state = json.loads(
        STATE_PATH.read_text(encoding="utf-8")
    )

    structural_summary = json.loads(
        STRUCTURAL_SUMMARY_PATH.read_text(encoding="utf-8")
    )

    contact_summary = json.loads(
        CONTACT_SUMMARY_PATH.read_text(encoding="utf-8")
    )

    state_checks = {
        "v2b_executed": (
            workflow_state.get("v2b_executed") is True
        ),
        "v2b_validation_pass": (
            workflow_state.get("v2b_validation_pass") is True
        ),
        "structural_gate_pass": (
            structural_summary.get("v2b_structural_gate_pass")
            is True
        ),
        "geometry_acceptance_retained": (
            contact_summary.get("geometry_acceptance_retained")
            is True
        ),
        "electronic_protocol_definition_authorized": (
            contact_summary.get(
                "electronic_property_protocol_definition_authorized"
            )
            is True
        ),
        "preparation_execution_flag_false": (
            preparation_summary.get(
                "electronic_diagnostic_execution_authorized"
            )
            is False
        ),
    }

    preflight_pass = (
        parsed["gate_pass"]
        and all(state_checks.values())
    )

    print("=" * 80)
    print("QM_F06 LOWER V2-B ELECTRONIC DIAGNOSTIC PREFLIGHT")
    print("=" * 80)
    print("ORCA:", orca_path)
    print("ORCA SHA256:", sha256(orca_path))
    print("Input:", INPUT_PATH)
    print("Input SHA256:", sha256(INPUT_PATH))
    print("Execution requested:", args.execute)
    print()

    print("Input checks:")
    for key, value in parsed["checks"].items():
        print(f"  {key:44s}: {value}")

    print()
    print("Scientific-state checks:")
    for key, value in state_checks.items():
        print(f"  {key:44s}: {value}")

    print()
    print(
        "Preflight decision:",
        "PASS" if preflight_pass else "FAIL",
    )

    if not preflight_pass:
        raise RuntimeError(
            "Electronic diagnostic preflight failed. "
            "Execution blocked."
        )

    if not args.execute:
        print()
        print(
            "DRY PREFLIGHT COMPLETE — ORCA WAS NOT EXECUTED."
        )
        return

    timestamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    execution_dir = EXECUTION_ROOT / (
        f"electronic_diagnostic_{timestamp}"
    )
    execution_dir.mkdir(
        parents=True,
        exist_ok=False,
    )

    execution_input = execution_dir / "diagnostic.inp"
    output_path = execution_dir / "diagnostic.out"
    stderr_path = execution_dir / "diagnostic.stderr"

    shutil.copy2(INPUT_PATH, execution_input)

    metadata = {
        "fragment": "QM_F06_LOWER_BOUNDARY_V2B",
        "stage": "ELECTRONIC_DIAGNOSTIC_SINGLE_POINT",
        "timestamp": timestamp,
        "orca_path": str(orca_path),
        "orca_sha256": sha256(orca_path),
        "input_sha256": sha256(execution_input),
        "python": sys.version,
        "platform": platform.platform(),
        "working_directory": str(execution_dir),
        "environment_path": os.environ.get("PATH", ""),
        "explicit_execute_flag": True,
        "geometry_optimization_requested": False,
        "resp_fitting_requested": False,
        "force_field_fitting_requested": False,
    }

    (
        execution_dir / "execution_metadata.json"
    ).write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )

    print()
    print("=" * 80)
    print("EXECUTING ELECTRONIC DIAGNOSTIC SINGLE POINT")
    print("=" * 80)
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
        completed.returncode,
    )

    (
        execution_dir / "diagnostic_validation.json"
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
        "Final energy:",
        validation["final_single_point_energy_hartree"],
    )
    print(
        "Execution gate:",
        (
            "PASS"
            if validation["execution_gate_pass"]
            else "FAIL"
        ),
    )
    print(
        "Complete diagnostic output:",
        validation["complete_diagnostic_output"],
    )

    if not validation["execution_gate_pass"]:
        raise RuntimeError(
            "Electronic diagnostic execution failed. "
            f"Inspect {execution_dir}"
        )

    preparation_summary[
        "electronic_diagnostic_executed"
    ] = True

    preparation_summary[
        "electronic_diagnostic_execution_gate_pass"
    ] = True

    preparation_summary[
        "complete_diagnostic_output_detected"
    ] = validation["complete_diagnostic_output"]

    preparation_summary[
        "execution_directory"
    ] = str(execution_dir.relative_to(ROOT))

    preparation_summary[
        "output_file"
    ] = str(output_path.relative_to(ROOT))

    preparation_summary[
        "final_single_point_energy_hartree"
    ] = validation["final_single_point_energy_hartree"]

    preparation_summary[
        "esp_resp_parameter_adoption_authorized"
    ] = False

    preparation_summary[
        "force_field_parameter_adoption_authorized"
    ] = False

    SUMMARY_PATH.write_text(
        json.dumps(preparation_summary, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
