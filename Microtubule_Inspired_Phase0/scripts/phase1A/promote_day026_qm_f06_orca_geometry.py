#!/usr/bin/env python3
"""
Generate the next ORCA input from an optimized XYZ geometry.

Examples
--------
Stage 1 -> Stage 2:
python scripts/phase1A/promote_day026_qm_f06_orca_geometry.py \
  --fragment QM_F06_LOWER_CAPPED_REPAIRED \
  --from-stage 1 \
  --xyz path/to/stage1_optimized.xyz

Stage 2 -> Stage 3:
python scripts/phase1A/promote_day026_qm_f06_orca_geometry.py \
  --fragment QM_F06_LOWER_CAPPED_REPAIRED \
  --from-stage 2 \
  --xyz path/to/stage2_optimized.xyz

This utility does not execute ORCA.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

WORKFLOW_DIR = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/"
    "QM_F06/orca_workflow"
)

VALID_FRAGMENTS = {
    "QM_F06_LOWER_CAPPED_REPAIRED",
    "QM_F06_UPPER_CAPPED_REPAIRED",
}


def read_xyz(path: Path) -> list[str]:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty XYZ file: {path}")

    lines = path.read_text(encoding="utf-8").splitlines()

    if len(lines) < 3:
        raise RuntimeError(f"Invalid XYZ file: {path}")

    try:
        expected_atoms = int(lines[0].strip())
    except ValueError as exc:
        raise RuntimeError(
            f"Invalid XYZ atom count in {path}"
        ) from exc

    coordinate_lines = [
        line.strip()
        for line in lines[2:]
        if line.strip()
    ]

    if len(coordinate_lines) != expected_atoms:
        raise RuntimeError(
            f"XYZ atom-count mismatch in {path}: "
            f"header={expected_atoms}, coordinates={len(coordinate_lines)}"
        )

    if expected_atoms != 22:
        raise RuntimeError(
            f"Expected 22 atoms; found {expected_atoms}"
        )

    parsed_lines: list[str] = []

    for line in coordinate_lines:
        parts = line.split()

        if len(parts) != 4:
            raise RuntimeError(
                f"Invalid coordinate line in {path}: {line}"
            )

        element = parts[0]
        x, y, z = map(float, parts[1:])

        parsed_lines.append(
            f"{element:<2s} {x: .10f} {y: .10f} {z: .10f}"
        )

    return parsed_lines


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fragment",
        required=True,
        choices=sorted(VALID_FRAGMENTS),
    )
    parser.add_argument(
        "--from-stage",
        required=True,
        type=int,
        choices=(1, 2),
    )
    parser.add_argument(
        "--xyz",
        required=True,
        type=Path,
    )

    args = parser.parse_args()

    fragment_dir = WORKFLOW_DIR / args.fragment

    if not fragment_dir.is_dir():
        raise RuntimeError(
            f"Workflow fragment directory not found: {fragment_dir}"
        )

    next_stage = args.from_stage + 1

    prefix_path = (
        fragment_dir / f"stage{next_stage}_prefix.txt"
    )
    state_path = (
        fragment_dir
        / f"stage{next_stage}_charge_multiplicity.txt"
    )

    if not prefix_path.is_file() or not state_path.is_file():
        raise RuntimeError(
            f"Stage-{next_stage} template components are missing."
        )

    coordinates = read_xyz(args.xyz)
    charge_multiplicity = state_path.read_text(
        encoding="utf-8"
    ).strip()

    output_path = fragment_dir / f"stage{next_stage}.inp"

    output_path.write_text(
        "\n".join(
            [
                prefix_path.read_text(
                    encoding="utf-8"
                ).rstrip(),
                "",
                (
                    f"# Coordinates promoted from Stage "
                    f"{args.from_stage}: {args.xyz}"
                ),
                f"* xyz {charge_multiplicity}",
                *coordinates,
                "*",
                "",
            ]
        ),
        encoding="utf-8",
    )

    workflow_state_path = fragment_dir / "workflow_state.json"
    state = json.loads(
        workflow_state_path.read_text(encoding="utf-8")
    )

    state[f"stage{args.from_stage}_geometry_promoted"] = True
    state[f"stage{next_stage}_input_generated"] = True
    state["qm_execution_authorized"] = False

    workflow_state_path.write_text(
        json.dumps(state, indent=2) + "\n",
        encoding="utf-8",
    )

    print(
        f"Stage-{next_stage} input generated from "
        f"Stage-{args.from_stage} XYZ."
    )
    print(f"Output: {output_path}")
    print("QM execution authorized: False")


if __name__ == "__main__":
    main()
