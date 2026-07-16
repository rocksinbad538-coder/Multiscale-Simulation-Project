#!/usr/bin/env python3
"""
Prepare the first ORCA optimization for QM_F06 LOWER Boundary V2.

V2-A objective
--------------
Relax only:
- the four restored real R2 boundary atoms;
- the three newly introduced V2 artificial caps.

Keep fixed:
- all 21 retained atoms from the converged V1 Stage-2 fragment.

This isolates the Stage2/Day024 boundary seam and prevents the already
validated B-N-B-N bridge geometry from being reoptimized prematurely.

No QM calculation is executed.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

V2_DIR = ROOT / (
    "runs/phase1A/day027_qm_f06_lower_boundary_redesign/"
    "QM_F06_LOWER_BOUNDARY_V2"
)

ATOMS_PATH = V2_DIR / (
    "cap02_repair/"
    "QM_F06_LOWER_BOUNDARY_V2_REPAIRED_atoms.csv"
)

XYZ_PATH = V2_DIR / (
    "cap02_repair/"
    "QM_F06_LOWER_BOUNDARY_V2_REPAIRED.xyz"
)

PRE_QM_SUMMARY = V2_DIR / (
    "repaired_pre_qm_audit/"
    "QM_F06_LOWER_BOUNDARY_V2_pre_qm_summary.json"
)

OUTPUT_DIR = V2_DIR / "orca_v2_workflow"

ORCA_PATH = Path(
    "/Users/alejandro/projects/"
    "orca_6_1_1_macosx_intel_openmpi411/orca"
)

MOBILE_ROLES = {
    "REAL_R2_BOUNDARY_EXPANSION_ATOM",
    "ARTIFICIAL_BOUNDARY_CAP_V2",
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def read_atoms(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if len(rows) != 28:
        raise RuntimeError(
            f"Expected 28 atoms; found {len(rows)}"
        )

    return rows


def read_xyz(path: Path) -> list[tuple[str, float, float, float]]:
    require_file(path)

    lines = path.read_text(encoding="utf-8").splitlines()
    expected = int(lines[0].strip())

    rows = [
        (
            parts[0],
            float(parts[1]),
            float(parts[2]),
            float(parts[3]),
        )
        for line in lines[2:]
        if line.strip()
        for parts in [line.split()]
    ]

    if len(rows) != expected:
        raise RuntimeError(
            f"XYZ count mismatch: {expected} vs {len(rows)}"
        )

    return rows


def write_csv(
    path: Path,
    rows: list[dict[str, Any]],
) -> None:
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


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    require_file(PRE_QM_SUMMARY)

    summary = json.loads(
        PRE_QM_SUMMARY.read_text(encoding="utf-8")
    )

    if summary.get("pre_qm_gate_pass") is not True:
        raise RuntimeError(
            "Repaired V2 pre-QM gate did not pass."
        )

    atoms = read_atoms(ATOMS_PATH)
    xyz_rows = read_xyz(XYZ_PATH)

    if len(atoms) != len(xyz_rows):
        raise RuntimeError(
            "Atom manifest and XYZ counts differ."
        )

    fixed_indices: list[int] = []
    mobile_indices: list[int] = []
    role_rows: list[dict[str, Any]] = []

    for index, (atom, xyz) in enumerate(
        zip(atoms, xyz_rows, strict=True)
    ):
        if atom["element"] != xyz[0]:
            raise RuntimeError(
                f"Element mismatch at index {index}: "
                f"{atom['element']} vs {xyz[0]}"
            )

        mobile = atom["atom_role"] in MOBILE_ROLES
        fixed = not mobile

        if fixed:
            fixed_indices.append(index)
        else:
            mobile_indices.append(index)

        role_rows.append(
            {
                "index_0based": index,
                "atom_id": atom["atom_id"],
                "element": atom["element"],
                "atom_role": atom["atom_role"],
                "coordinate_source": atom["coordinate_source"],
                "v2a_fixed": fixed,
                "v2a_mobile": mobile,
            }
        )

    if len(fixed_indices) != 21:
        raise RuntimeError(
            f"Expected 21 fixed atoms; found {len(fixed_indices)}"
        )

    if len(mobile_indices) != 7:
        raise RuntimeError(
            f"Expected 7 mobile atoms; found {len(mobile_indices)}"
        )

    constraint_lines = [
        f"    {{ C {index} C }}"
        for index in fixed_indices
    ]

    coordinate_lines = [
        (
            f"{element:<2s} "
            f"{x: .10f} "
            f"{y: .10f} "
            f"{z: .10f}"
        )
        for element, x, y, z in xyz_rows
    ]

    input_text = "\n".join(
        [
            (
                "! PBE0 D4 def2-TZVP def2/J "
                "RIJCOSX TightSCF Opt DefGrid3"
            ),
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
            "%geom",
            "  MaxIter 200",
            "  Constraints",
            *constraint_lines,
            "  end",
            "end",
            "",
            (
                "# QM_F06 LOWER Boundary V2-A: "
                "relax restored boundary and new caps only"
            ),
            "* xyz 0 1",
            *coordinate_lines,
            "*",
            "",
        ]
    )

    input_path = OUTPUT_DIR / "v2a_boundary_relax.inp"
    input_path.write_text(
        input_text,
        encoding="utf-8",
    )

    write_csv(
        OUTPUT_DIR / "v2a_atom_role_constraint_map.csv",
        role_rows,
    )

    workflow_state = {
        "fragment": "QM_F06_LOWER_BOUNDARY_V2_REPAIRED",
        "formula": "B6N7H15",
        "charge": 0,
        "multiplicity": 1,
        "pre_qm_gate_pass": True,
        "v2a_input_prepared": True,
        "v2a_fixed_atom_count": len(fixed_indices),
        "v2a_mobile_atom_count": len(mobile_indices),
        "v2a_fixed_indices": fixed_indices,
        "v2a_mobile_indices": mobile_indices,
        "v2a_executed": False,
        "v2a_geometry_promoted": False,
        "v2b_input_generated": False,
        "qm_execution_authorized": False,
    }

    (
        OUTPUT_DIR / "v2_workflow_state.json"
    ).write_text(
        json.dumps(workflow_state, indent=2) + "\n",
        encoding="utf-8",
    )

    report_path = (
        OUTPUT_DIR
        / "QM_F06_LOWER_BOUNDARY_V2_ORCA_PREPARATION.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 LOWER Boundary V2 ORCA Preparation — Day027",
                "",
                "## Electronic-structure model",
                "",
                "- Functional: **PBE0**",
                "- Dispersion: **D4**",
                "- Basis: **def2-TZVP**",
                "- RI approximation: **RIJCOSX / def2-J**",
                "- Grid: **DefGrid3**",
                "- SCF: **TightSCF**",
                "- Charge/multiplicity: **0 / 1**",
                "",
                "## V2-A optimization",
                "",
                "- Total atoms: **28**",
                f"- Fixed atoms: **{len(fixed_indices)}**",
                f"- Mobile atoms: **{len(mobile_indices)}**",
                (
                    "- Mobile region: four restored real R2 atoms "
                    "and three newly added V2 caps"
                ),
                (
                    "- Fixed region: 21 atoms retained from the "
                    "converged V1 Stage-2 fragment"
                ),
                "",
                "## Authorization state",
                "",
                "- Input preparation: **COMPLETED**",
                "- QM execution: **NOT AUTHORIZED**",
                "",
                "## Required next step",
                "",
                "Perform static input audit and dry preflight.",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    manifest = {
        "input_file": str(input_path.relative_to(ROOT)),
        "input_sha256": sha256(input_path),
        "atom_manifest": str(ATOMS_PATH.relative_to(ROOT)),
        "atom_manifest_sha256": sha256(ATOMS_PATH),
        "coordinate_file": str(XYZ_PATH.relative_to(ROOT)),
        "coordinate_file_sha256": sha256(XYZ_PATH),
        "orca_detected": ORCA_PATH.is_file(),
        "orca_path": str(ORCA_PATH),
        "qm_execution_authorized": False,
    }

    (
        OUTPUT_DIR / "v2a_input_manifest.json"
    ).write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )

    print("QM_F06 LOWER Boundary V2 ORCA input prepared.")
    print("Input:", input_path)
    print("Fixed atoms:", len(fixed_indices))
    print("Mobile atoms:", len(mobile_indices))
    print("ORCA detected:", ORCA_PATH.is_file())
    print("QM execution authorized: False")


if __name__ == "__main__":
    main()
