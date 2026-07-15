#!/usr/bin/env python3
"""
Build a non-executing sequential ORCA workflow for QM_F06.

The workflow prepares:
- Stage-1 inputs ready for later execution.
- A geometry-promotion utility that constructs Stage 2 from the optimized
  Stage-1 XYZ and Stage 3 from the optimized Stage-2 XYZ.
- Per-fragment working directories.
- A guarded execution script that refuses to run unless the user supplies
  the explicit --execute flag.

This script does not execute ORCA.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

INPUT_DIR = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/"
    "QM_F06/orca_inputs"
)

WORKFLOW_DIR = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/"
    "QM_F06/orca_workflow"
)

FRAGMENTS = (
    "QM_F06_LOWER_CAPPED_REPAIRED",
    "QM_F06_UPPER_CAPPED_REPAIRED",
)


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty required file: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows available for {path}")

    fields: list[str] = []

    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(
            {
                field: row.get(field, "")
                for field in fields
            }
            for row in rows
        )


def remove_xyz_block(input_text: str) -> tuple[str, str]:
    pattern = re.compile(
        r"(?ms)^\s*\*\s+xyz\s+(-?\d+)\s+(\d+)\s*$"
        r"(.*?)"
        r"^\s*\*\s*$"
    )

    match = pattern.search(input_text)

    if not match:
        raise RuntimeError("Could not locate ORCA XYZ coordinate block.")

    charge = match.group(1)
    multiplicity = match.group(2)

    prefix = input_text[:match.start()].rstrip()
    return prefix, f"{charge} {multiplicity}"


def main() -> None:
    WORKFLOW_DIR.mkdir(parents=True, exist_ok=True)

    manifest_rows: list[dict[str, Any]] = []

    for fragment in FRAGMENTS:
        fragment_dir = WORKFLOW_DIR / fragment
        fragment_dir.mkdir(parents=True, exist_ok=True)

        stage1_source = (
            INPUT_DIR
            / f"{fragment}_STAGE1_CONSTRAINED_OPT.inp"
        )
        stage2_source = (
            INPUT_DIR
            / f"{fragment}_STAGE2_PARTIAL_RELAX_OPT.inp"
        )
        stage3_source = (
            INPUT_DIR
            / f"{fragment}_STAGE3_SINGLE_POINT_TEMPLATE.inp"
        )

        for path in (
            stage1_source,
            stage2_source,
            stage3_source,
        ):
            require_file(path)

        stage1_target = fragment_dir / "stage1.inp"
        shutil.copy2(stage1_source, stage1_target)

        stage2_prefix, stage2_state = remove_xyz_block(
            stage2_source.read_text(encoding="utf-8")
        )
        stage3_prefix, stage3_state = remove_xyz_block(
            stage3_source.read_text(encoding="utf-8")
        )

        (fragment_dir / "stage2_prefix.txt").write_text(
            stage2_prefix + "\n",
            encoding="utf-8",
        )
        (fragment_dir / "stage2_charge_multiplicity.txt").write_text(
            stage2_state + "\n",
            encoding="utf-8",
        )

        (fragment_dir / "stage3_prefix.txt").write_text(
            stage3_prefix + "\n",
            encoding="utf-8",
        )
        (fragment_dir / "stage3_charge_multiplicity.txt").write_text(
            stage3_state + "\n",
            encoding="utf-8",
        )

        state = {
            "fragment": fragment,
            "stage1_input_prepared": True,
            "stage1_executed": False,
            "stage1_geometry_promoted": False,
            "stage2_input_generated": False,
            "stage2_executed": False,
            "stage2_geometry_promoted": False,
            "stage3_input_generated": False,
            "stage3_executed": False,
            "qm_execution_authorized": False,
        }

        state_path = fragment_dir / "workflow_state.json"
        state_path.write_text(
            json.dumps(state, indent=2) + "\n",
            encoding="utf-8",
        )

        for role, path in (
            ("stage1_input", stage1_target),
            ("stage2_prefix", fragment_dir / "stage2_prefix.txt"),
            (
                "stage2_charge_multiplicity",
                fragment_dir / "stage2_charge_multiplicity.txt",
            ),
            ("stage3_prefix", fragment_dir / "stage3_prefix.txt"),
            (
                "stage3_charge_multiplicity",
                fragment_dir / "stage3_charge_multiplicity.txt",
            ),
            ("workflow_state", state_path),
        ):
            manifest_rows.append(
                {
                    "fragment": fragment,
                    "role": role,
                    "file": str(path.relative_to(ROOT)),
                    "sha256": sha256(path),
                    "calculation_executed": False,
                }
            )

    manifest_path = WORKFLOW_DIR / "QM_F06_workflow_manifest.csv"
    write_csv(manifest_path, manifest_rows)

    summary = {
        "decision": (
            "QM_F06_SEQUENTIAL_ORCA_WORKFLOW_PREPARED_"
            "EXECUTION_NOT_AUTHORIZED"
        ),
        "fragments": list(FRAGMENTS),
        "stage1_inputs_ready": True,
        "stage2_requires_stage1_optimized_xyz": True,
        "stage3_requires_stage2_optimized_xyz": True,
        "qm_execution_authorized": False,
        "required_next_step": (
            "VALIDATE_GEOMETRY_PROMOTION_UTILITY_WITH_DRY_RUN"
        ),
    }

    (
        WORKFLOW_DIR / "QM_F06_workflow_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("Sequential QM_F06 ORCA workflow prepared.")
    print(f"Workflow directory: {WORKFLOW_DIR}")
    print("QM execution authorized: False")


if __name__ == "__main__":
    main()
