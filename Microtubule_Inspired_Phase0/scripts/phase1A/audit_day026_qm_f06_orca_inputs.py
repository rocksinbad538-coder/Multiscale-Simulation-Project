#!/usr/bin/env python3
"""
Technical audit of prepared QM_F06 ORCA inputs.

Checks:
- input files exist and are non-empty;
- 22 atoms are present;
- charge/multiplicity are 0/1;
- constraint indices are valid, unique and 0-based;
- Stage 1 fixes exactly artificial caps plus peripheral heavy atoms;
- Stage 2 fixes only peripheral heavy atoms;
- bridge and attachment atoms remain mobile;
- Stage 1 and Stage 2 coordinate blocks match the repaired atom manifests;
- no ORCA calculation is executed.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

INPUT_DIR = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/"
    "QM_F06/orca_inputs"
)

F06_DIR = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/QM_F06"
)

ENDS = ("LOWER", "UPPER")
EXPECTED_ATOMS = 22


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(f"No rows found in {path}")

    return rows


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

        for row in rows:
            writer.writerow(
                {
                    field: row.get(field, "")
                    for field in fields
                }
            )


def parse_orca_input(
    path: Path,
) -> tuple[
    int,
    int,
    list[tuple[str, float, float, float]],
    list[int],
]:
    require_file(path)
    lines = path.read_text(encoding="utf-8").splitlines()

    xyz_start = None
    charge = None
    multiplicity = None

    for index, line in enumerate(lines):
        match = re.match(
            r"^\s*\*\s+xyz\s+(-?\d+)\s+(\d+)\s*$",
            line,
            flags=re.IGNORECASE,
        )

        if match:
            xyz_start = index + 1
            charge = int(match.group(1))
            multiplicity = int(match.group(2))
            break

    if xyz_start is None or charge is None or multiplicity is None:
        raise RuntimeError(
            f"Could not locate '* xyz charge multiplicity' in {path}"
        )

    atoms: list[tuple[str, float, float, float]] = []

    for line in lines[xyz_start:]:
        if line.strip() == "*":
            break

        parts = line.split()

        if len(parts) != 4:
            raise RuntimeError(
                f"Invalid XYZ coordinate line in {path}: {line}"
            )

        atoms.append(
            (
                parts[0],
                float(parts[1]),
                float(parts[2]),
                float(parts[3]),
            )
        )

    constraints = []

    for line in lines:
        match = re.match(
            r"^\s*\{\s*C\s+(\d+)\s+C\s*\}\s*$",
            line,
            flags=re.IGNORECASE,
        )

        if match:
            constraints.append(int(match.group(1)))

    return charge, multiplicity, atoms, constraints


def main() -> None:
    role_map = read_csv(
        INPUT_DIR
        / "QM_F06_orca_atom_role_and_constraint_map.csv"
    )

    audit_rows: list[dict[str, Any]] = []
    report_sections: list[str] = []
    all_pass = True

    for end in ENDS:
        fragment = f"QM_F06_{end}_CAPPED_REPAIRED"

        roles = [
            row
            for row in role_map
            if row["fragment"] == fragment
        ]

        if len(roles) != EXPECTED_ATOMS:
            raise RuntimeError(
                f"{fragment}: expected {EXPECTED_ATOMS} role rows; "
                f"found {len(roles)}"
            )

        repaired_atoms = read_csv(
            F06_DIR / f"{fragment}_atoms.csv"
        )

        for stage, fixed_field in (
            ("STAGE1_CONSTRAINED_OPT", "stage1_fixed"),
            ("STAGE2_PARTIAL_RELAX_OPT", "stage2_fixed"),
        ):
            input_path = (
                INPUT_DIR / f"{fragment}_{stage}.inp"
            )

            (
                charge,
                multiplicity,
                input_atoms,
                constraint_indices,
            ) = parse_orca_input(input_path)

            expected_fixed = sorted(
                int(row["orca_atom_index_0based"])
                for row in roles
                if row[fixed_field].lower() == "true"
            )

            unique_constraints = sorted(set(constraint_indices))

            index_range_pass = all(
                0 <= index < EXPECTED_ATOMS
                for index in constraint_indices
            )

            unique_pass = (
                len(unique_constraints)
                == len(constraint_indices)
            )

            constraint_match_pass = (
                unique_constraints == expected_fixed
            )

            atom_count_pass = (
                len(input_atoms) == EXPECTED_ATOMS
                == len(repaired_atoms)
            )

            charge_multiplicity_pass = (
                charge == 0
                and multiplicity == 1
            )

            coordinate_match_pass = True
            maximum_coordinate_error = 0.0

            if atom_count_pass:
                for input_atom, source_atom in zip(
                    input_atoms,
                    repaired_atoms,
                    strict=True,
                ):
                    element, x, y, z = input_atom

                    if element != source_atom["element"]:
                        coordinate_match_pass = False
                        break

                    source_xyz = (
                        float(source_atom["x_angstrom"]),
                        float(source_atom["y_angstrom"]),
                        float(source_atom["z_angstrom"]),
                    )

                    error = max(
                        abs(x - source_xyz[0]),
                        abs(y - source_xyz[1]),
                        abs(z - source_xyz[2]),
                    )

                    maximum_coordinate_error = max(
                        maximum_coordinate_error,
                        error,
                    )

                coordinate_match_pass = (
                    coordinate_match_pass
                    and maximum_coordinate_error <= 1.0e-8
                )
            else:
                coordinate_match_pass = False

            mobile_bridge_or_attachment_failures = []

            for row in roles:
                must_remain_mobile = (
                    row["is_bridge_atom"].lower() == "true"
                    or row["is_attachment_center"].lower() == "true"
                )

                if (
                    must_remain_mobile
                    and int(row["orca_atom_index_0based"])
                    in unique_constraints
                ):
                    mobile_bridge_or_attachment_failures.append(
                        row["atom_id"]
                    )

            scientific_constraint_pass = (
                not mobile_bridge_or_attachment_failures
            )

            stage_pass = all(
                (
                    atom_count_pass,
                    charge_multiplicity_pass,
                    index_range_pass,
                    unique_pass,
                    constraint_match_pass,
                    coordinate_match_pass,
                    scientific_constraint_pass,
                )
            )

            all_pass = all_pass and stage_pass

            audit_rows.append(
                {
                    "fragment": fragment,
                    "stage": stage,
                    "input_file": str(input_path.relative_to(ROOT)),
                    "atom_count": len(input_atoms),
                    "charge": charge,
                    "multiplicity": multiplicity,
                    "constraint_count": len(constraint_indices),
                    "expected_constraint_count": len(expected_fixed),
                    "constraint_indices": "|".join(
                        str(index)
                        for index in unique_constraints
                    ),
                    "expected_constraint_indices": "|".join(
                        str(index)
                        for index in expected_fixed
                    ),
                    "atom_count_gate_pass": atom_count_pass,
                    "charge_multiplicity_gate_pass": (
                        charge_multiplicity_pass
                    ),
                    "constraint_index_range_gate_pass": (
                        index_range_pass
                    ),
                    "constraint_uniqueness_gate_pass": unique_pass,
                    "constraint_map_gate_pass": (
                        constraint_match_pass
                    ),
                    "coordinate_match_gate_pass": (
                        coordinate_match_pass
                    ),
                    "maximum_coordinate_error_angstrom": (
                        f"{maximum_coordinate_error:.12e}"
                    ),
                    "bridge_attachment_mobile_gate_pass": (
                        scientific_constraint_pass
                    ),
                    "incorrectly_fixed_bridge_or_attachment_atoms": (
                        "|".join(
                            mobile_bridge_or_attachment_failures
                        )
                    ),
                    "technical_input_gate_pass": stage_pass,
                }
            )

        report_sections.extend(
            [
                f"## {fragment}",
                "",
                (
                    "- Stage-1 expected fixed atoms: "
                    f"**{sum(row['stage1_fixed'].lower() == 'true' for row in roles)}**"
                ),
                (
                    "- Stage-2 expected fixed atoms: "
                    f"**{sum(row['stage2_fixed'].lower() == 'true' for row in roles)}**"
                ),
                (
                    "- Bridge/attachment atoms incorrectly fixed: "
                    "**0 expected**"
                ),
                "",
            ]
        )

    audit_path = (
        INPUT_DIR / "QM_F06_orca_input_technical_audit.csv"
    )
    write_csv(audit_path, audit_rows)

    decision = (
        "QM_F06_ORCA_INPUTS_PASS_TECHNICAL_AUDIT"
        if all_pass
        else "QM_F06_ORCA_INPUTS_FAIL_TECHNICAL_AUDIT"
    )

    report_path = (
        INPUT_DIR / "QM_F06_ORCA_INPUT_TECHNICAL_AUDIT.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 ORCA Input Technical Audit — Day026",
                "",
                f"## Decision: **{decision}**",
                "",
                (
                    "Stage-1 and Stage-2 input files were checked "
                    "against the repaired atom manifests and the "
                    "atom-role/constraint map."
                ),
                "",
                *report_sections,
                "## Workflow dependency",
                "",
                (
                    "The current Stage-2 inputs contain the repaired "
                    "initial coordinates. They are templates only. Before "
                    "execution, their coordinate blocks must be replaced "
                    "with the optimized Stage-1 geometries."
                ),
                "",
                (
                    "Likewise, Stage 3 must be generated from the optimized "
                    "Stage-2 geometry rather than from the current initial "
                    "coordinates."
                ),
                "",
                "## Authorization state",
                "",
                "- Static input audit: "
                f"**{'PASSED' if all_pass else 'FAILED'}**",
                "- Sequential workflow preparation: **PENDING**",
                "- QM execution: **NOT AUTHORIZED**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = {
        "decision": decision,
        "static_input_audit_pass": all_pass,
        "stage2_currently_uses_stage1_output": False,
        "stage3_currently_uses_stage2_output": False,
        "qm_execution_authorized": False,
        "required_next_step": (
            "BUILD_SEQUENTIAL_ORCA_WORKFLOW"
            if all_pass
            else "CORRECT_ORCA_INPUT_GENERATION"
        ),
    }

    (
        INPUT_DIR
        / "QM_F06_orca_input_technical_audit_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("QM_F06 ORCA input technical audit completed.")
    print(f"Decision: {decision}")
    print("QM execution authorized: False")


if __name__ == "__main__":
    main()
