#!/usr/bin/env python3
"""
Prepare ORCA input files for the repaired QM_F06 bridge fragments.

Prepared workflow
-----------------
Stage 1:
    Constrained optimization.
    - Bridge atoms: mobile.
    - Bridge attachment centers: mobile.
    - Existing bridge/passivant H atoms: mobile.
    - Peripheral scaffold atoms: fixed.
    - Artificial boundary caps: fixed.

Stage 2:
    Partially relaxed optimization.
    - Bridge atoms and attachment centers: mobile.
    - All H atoms, including artificial caps: mobile.
    - Peripheral scaffold heavy atoms: fixed.

Stage 3:
    Single-point calculation on the Stage-2 geometry.
    This input is prepared but must be updated to read the optimized
    Stage-2 XYZ geometry after Stage 2 finishes.

No calculation is executed by this script.

Electronic-structure level
--------------------------
PBE0-D4/def2-TZVP with RIJCOSX, TightSCF and a dense integration grid.

The method is a defensible initial geometry/reference level for B/N/H
clusters. Final parameterization and charge protocols remain separate
scientific decisions.
"""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

F06_DIR = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/QM_F06"
)

OUTPUT_DIR = F06_DIR / "orca_inputs"

ENDS = ("LOWER", "UPPER")

CHARGE = 0
MULTIPLICITY = 1

METHOD_LINE_OPT = (
    "! PBE0 D4 def2-TZVP def2/J RIJCOSX TightSCF "
    "Opt DefGrid3"
)

METHOD_LINE_SP = (
    "! PBE0 D4 def2-TZVP def2/J RIJCOSX TightSCF "
    "DefGrid3"
)


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty required file: {path}")


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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def xyz_lines(atoms: list[dict[str, str]]) -> list[str]:
    lines: list[str] = []

    for atom in atoms:
        lines.append(
            f"{atom['element']:<2s} "
            f"{float(atom['x_angstrom']): .10f} "
            f"{float(atom['y_angstrom']): .10f} "
            f"{float(atom['z_angstrom']): .10f}"
        )

    return lines


def constraint_block(indices: list[int]) -> str:
    if not indices:
        return ""

    rows = [
        "%geom",
        "  Constraints",
    ]

    for index in indices:
        rows.append(f"    {{ C {index} C }}")

    rows.extend(
        [
            "  end",
            "end",
        ]
    )

    return "\n".join(rows)


def write_orca_input(
    path: Path,
    method_line: str,
    atoms: list[dict[str, str]],
    fixed_indices: list[int],
    title: str,
) -> None:
    sections = [
        method_line,
        "",
        "%pal",
        "  nprocs 4",
        "end",
        "",
        "%maxcore 3000",
        "",
        "%scf",
        "  MaxIter 500",
        "end",
        "",
    ]

    constraints = constraint_block(fixed_indices)

    if constraints:
        sections.extend(
            [
                constraints,
                "",
            ]
        )

    sections.extend(
        [
            f"# {title}",
            f"* xyz {CHARGE} {MULTIPLICITY}",
            *xyz_lines(atoms),
            "*",
            "",
        ]
    )

    path.write_text(
        "\n".join(sections),
        encoding="utf-8",
    )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    extraction_manifest_rows: list[dict[str, Any]] = []
    atom_role_rows: list[dict[str, Any]] = []
    report_sections: list[str] = []

    orca_executable = shutil.which("orca")

    for end in ENDS:
        source_label = f"QM_F06_{end}_CAPPED_REPAIRED"

        repaired_atoms_path = (
            F06_DIR / f"{source_label}_atoms.csv"
        )

        original_fragment_path = (
            F06_DIR / f"QM_F06_{end}_atoms.csv"
        )

        repaired_atoms = read_csv(repaired_atoms_path)
        original_fragment = read_csv(original_fragment_path)

        original_flags = {
            row["node_id"]: {
                "is_bridge_atom": (
                    row["is_bridge_atom"].strip().lower()
                    == "true"
                ),
                "is_attachment_center": (
                    row["is_attachment_center"].strip().lower()
                    == "true"
                ),
            }
            for row in original_fragment
        }

        stage1_fixed: list[int] = []
        stage2_fixed: list[int] = []

        for index, atom in enumerate(repaired_atoms):
            atom_id = atom["atom_id"]
            element = atom["element"]
            role = atom["atom_role"]

            flags = original_flags.get(
                atom_id,
                {
                    "is_bridge_atom": False,
                    "is_attachment_center": False,
                },
            )

            is_bridge = flags["is_bridge_atom"]
            is_attachment = flags["is_attachment_center"]
            is_artificial_cap = (
                role == "ARTIFICIAL_BOUNDARY_CAP"
            )
            is_hydrogen = element == "H"

            # Stage 1:
            # Fix artificial caps and peripheral scaffold atoms.
            stage1_mobile = (
                is_bridge
                or is_attachment
                or (
                    is_hydrogen
                    and not is_artificial_cap
                )
            )

            if not stage1_mobile:
                stage1_fixed.append(index)

            # Stage 2:
            # All H atoms mobile; peripheral scaffold heavy atoms fixed.
            stage2_mobile = (
                is_bridge
                or is_attachment
                or is_hydrogen
            )

            if not stage2_mobile:
                stage2_fixed.append(index)

            atom_role_rows.append(
                {
                    "fragment": source_label,
                    "orca_atom_index_0based": index,
                    "orca_atom_index_1based": index + 1,
                    "atom_id": atom_id,
                    "element": element,
                    "atom_role": role,
                    "is_bridge_atom": is_bridge,
                    "is_attachment_center": is_attachment,
                    "stage1_fixed": not stage1_mobile,
                    "stage2_fixed": not stage2_mobile,
                }
            )

        stage1_path = (
            OUTPUT_DIR
            / f"{source_label}_STAGE1_CONSTRAINED_OPT.inp"
        )

        stage2_path = (
            OUTPUT_DIR
            / f"{source_label}_STAGE2_PARTIAL_RELAX_OPT.inp"
        )

        stage3_path = (
            OUTPUT_DIR
            / f"{source_label}_STAGE3_SINGLE_POINT_TEMPLATE.inp"
        )

        write_orca_input(
            stage1_path,
            METHOD_LINE_OPT,
            repaired_atoms,
            stage1_fixed,
            (
                f"{source_label}: constrained optimization; "
                "caps and peripheral scaffold fixed"
            ),
        )

        write_orca_input(
            stage2_path,
            METHOD_LINE_OPT,
            repaired_atoms,
            stage2_fixed,
            (
                f"{source_label}: partial relaxation; "
                "all H atoms mobile; peripheral heavy atoms fixed"
            ),
        )

        write_orca_input(
            stage3_path,
            METHOD_LINE_SP,
            repaired_atoms,
            [],
            (
                f"{source_label}: single-point template. "
                "Replace coordinates with the optimized Stage-2 geometry."
            ),
        )

        for stage, path in (
            ("STAGE1_CONSTRAINED_OPT", stage1_path),
            ("STAGE2_PARTIAL_RELAX_OPT", stage2_path),
            ("STAGE3_SINGLE_POINT_TEMPLATE", stage3_path),
        ):
            extraction_manifest_rows.append(
                {
                    "fragment": source_label,
                    "stage": stage,
                    "file": str(path.relative_to(ROOT)),
                    "sha256": sha256(path),
                    "calculation_executed": False,
                }
            )

        report_sections.extend(
            [
                f"## {source_label}",
                "",
                f"- Atoms: **{len(repaired_atoms)}**",
                (
                    "- Stage-1 fixed atoms: "
                    f"**{len(stage1_fixed)}**"
                ),
                (
                    "- Stage-1 mobile atoms: "
                    f"**{len(repaired_atoms) - len(stage1_fixed)}**"
                ),
                (
                    "- Stage-2 fixed atoms: "
                    f"**{len(stage2_fixed)}**"
                ),
                (
                    "- Stage-2 mobile atoms: "
                    f"**{len(repaired_atoms) - len(stage2_fixed)}**"
                ),
                "- Charge/multiplicity: **0 / 1**",
                "",
            ]
        )

    write_csv(
        OUTPUT_DIR / "QM_F06_orca_input_manifest.csv",
        extraction_manifest_rows,
    )

    write_csv(
        OUTPUT_DIR / "QM_F06_orca_atom_role_and_constraint_map.csv",
        atom_role_rows,
    )

    report_path = (
        OUTPUT_DIR / "QM_F06_ORCA_INPUT_PREPARATION_REPORT.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 ORCA Input Preparation — Day026",
                "",
                "## Decision",
                "",
                (
                    "**QM_F06_ORCA_INPUTS_PREPARED_EXECUTION_NOT_AUTHORIZED**"
                ),
                "",
                "## Electronic-structure model",
                "",
                "- Functional: **PBE0**",
                "- Dispersion: **D4**",
                "- Basis: **def2-TZVP**",
                "- Approximation: **RIJCOSX / def2-J**",
                "- SCF: **TightSCF, maximum 500 iterations**",
                "- Charge: **0**",
                "- Multiplicity: **1**",
                "",
                "## Optimization hierarchy",
                "",
                (
                    "1. Stage 1 keeps artificial caps and peripheral "
                    "scaffold atoms fixed."
                ),
                (
                    "2. Stage 2 releases all hydrogen atoms while retaining "
                    "peripheral heavy-atom constraints."
                ),
                (
                    "3. Stage 3 is a single-point template and must receive "
                    "the optimized Stage-2 coordinates."
                ),
                "",
                *report_sections,
                "## Executable detection",
                "",
                (
                    f"- ORCA detected in current shell: "
                    f"**{'YES' if orca_executable else 'NO'}**"
                ),
                (
                    f"- Detected path: "
                    f"`{orca_executable or 'NOT_FOUND'} `"
                ),
                "",
                "## Authorization state",
                "",
                "- Input preparation: **COMPLETED**",
                "- QM execution: **NOT AUTHORIZED**",
                "- Charge fitting: **NOT YET DEFINED**",
                "- Force-field fitting: **NOT YET EXECUTED**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = {
        "decision": (
            "QM_F06_ORCA_INPUTS_PREPARED_EXECUTION_NOT_AUTHORIZED"
        ),
        "orca_detected": bool(orca_executable),
        "orca_path": orca_executable,
        "method": "PBE0-D4/def2-TZVP",
        "charge": CHARGE,
        "multiplicity": MULTIPLICITY,
        "input_preparation_completed": True,
        "qm_execution_authorized": False,
        "required_next_step": (
            "TECHNICAL_REVIEW_OF_INPUTS_AND_ORCA_AVAILABILITY"
        ),
    }

    (
        OUTPUT_DIR / "QM_F06_orca_input_preparation_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("QM_F06 ORCA input preparation completed.")
    print(f"Output directory: {OUTPUT_DIR}")
    print("ORCA detected:", bool(orca_executable))
    print("ORCA path:", orca_executable or "NOT_FOUND")
    print("QM execution authorized: False")


if __name__ == "__main__":
    main()
