#!/usr/bin/env python3
"""
Prepare QM_F06 LOWER Boundary V2-B optimization.

Input geometry:
- converged V2-A optimized XYZ.

V2-B mobile atoms:
- every hydrogen atom;
- restored real heavy atoms:
  A:LOWER:11:-3
  A:LOWER:13:-3
  A:LOWER:14:-2

V2-B fixed atoms:
- the ten remaining heavy atoms forming the previously validated
  bridge/scaffold core.

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

WORKFLOW_DIR = V2_DIR / "orca_v2_workflow"
STATE_PATH = WORKFLOW_DIR / "v2_workflow_state.json"

ATOM_MANIFEST = V2_DIR / (
    "cap02_repair/"
    "QM_F06_LOWER_BOUNDARY_V2_REPAIRED_atoms.csv"
)

V2A_VALIDATION = V2_DIR / (
    "v2a_postprocessing/"
    "QM_F06_LOWER_BOUNDARY_V2A_validation_summary.json"
)

RESTORED_HEAVY_ATOMS = {
    "A:LOWER:11:-3",
    "A:LOWER:13:-3",
    "A:LOWER:14:-2",
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


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"No records in {path}")
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
            f"XYZ count mismatch: header={expected}, rows={len(rows)}"
        )

    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []

    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(
            {field: row.get(field, "") for field in fields}
            for row in rows
        )


def main() -> None:
    require_file(STATE_PATH)
    require_file(V2A_VALIDATION)

    state = json.loads(STATE_PATH.read_text(encoding="utf-8"))
    validation = json.loads(
        V2A_VALIDATION.read_text(encoding="utf-8")
    )

    required_state = {
        "v2a_executed": state.get("v2a_executed") is True,
        "v2a_validation_pass": (
            state.get("v2a_validation_pass") is True
        ),
        "v2a_optimized_xyz": bool(state.get("v2a_optimized_xyz")),
    }

    if not all(required_state.values()):
        raise RuntimeError(
            f"Incomplete V2-A workflow state: {required_state}"
        )

    # V2-A global steric gate failed only because one peripheral,
    # noncovalent artificial-cap contact remains.
    scientific_checks = {
        "fixed_atom_motion_failures_zero": (
            validation["fixed_atom_motion_failures"] == 0
        ),
        "graph_connected": validation["graph_connected"] is True,
        "valence_failures_zero": (
            validation["valence_failures"] == 0
        ),
        "bond_failures_zero": (
            validation["bond_failures"] == 0
        ),
        "cap_bond_failures_zero": (
            validation["cap_bond_failures"] == 0
        ),
        "bridge_cap_hard_contacts_zero": (
            validation["bridge_cap_hard_contacts"] == 0
        ),
        "unintended_covalent_contacts_zero": (
            validation["unintended_covalent_contacts"] == 0
        ),
        "exactly_one_peripheral_hard_cap_contact": (
            validation["hard_cap_contacts"] == 1
        ),
    }

    if not all(scientific_checks.values()):
        raise RuntimeError(
            "V2-A does not satisfy the chemically reclassified "
            f"V2-B readiness conditions: {scientific_checks}"
        )

    atoms = read_csv(ATOM_MANIFEST)

    if len(atoms) != 28:
        raise RuntimeError(f"Expected 28 atoms; found {len(atoms)}")

    optimized_xyz_path = ROOT / state["v2a_optimized_xyz"]
    xyz_rows = read_xyz(optimized_xyz_path)

    if len(xyz_rows) != len(atoms):
        raise RuntimeError("Manifest/XYZ atom-count mismatch.")

    fixed_indices: list[int] = []
    mobile_indices: list[int] = []
    map_rows: list[dict[str, Any]] = []

    for index, (atom, xyz) in enumerate(
        zip(atoms, xyz_rows, strict=True)
    ):
        if atom["element"] != xyz[0]:
            raise RuntimeError(
                f"Element-order mismatch at index {index}: "
                f"{atom['atom_id']}"
            )

        atom_id = atom["atom_id"]
        element = atom["element"]

        mobile = (
            element == "H"
            or atom_id in RESTORED_HEAVY_ATOMS
        )
        fixed = not mobile

        if fixed:
            fixed_indices.append(index)
        else:
            mobile_indices.append(index)

        map_rows.append(
            {
                "index_0based": index,
                "atom_id": atom_id,
                "element": element,
                "atom_role": atom["atom_role"],
                "v2b_fixed": fixed,
                "v2b_mobile": mobile,
                "mobility_basis": (
                    "ALL_HYDROGENS"
                    if element == "H"
                    else (
                        "RESTORED_REAL_HEAVY_ATOM"
                        if atom_id in RESTORED_HEAVY_ATOMS
                        else "VALIDATED_CORE_HEAVY_ATOM"
                    )
                ),
            }
        )

    if len(fixed_indices) != 10:
        raise RuntimeError(
            f"Expected 10 fixed heavy atoms; found {len(fixed_indices)}"
        )

    if len(mobile_indices) != 18:
        raise RuntimeError(
            f"Expected 18 mobile atoms; found {len(mobile_indices)}"
        )

    constraints = [
        f"    {{ C {index} C }}"
        for index in fixed_indices
    ]

    coordinates = [
        f"{element:<2s} {x: .10f} {y: .10f} {z: .10f}"
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
            *constraints,
            "  end",
            "end",
            "",
            (
                "# QM_F06 LOWER Boundary V2-B; "
                "all H and restored heavy atoms mobile; "
                "validated core heavy atoms fixed"
            ),
            "* xyz 0 1",
            *coordinates,
            "*",
            "",
        ]
    )

    input_path = WORKFLOW_DIR / "v2b_partial_relax.inp"
    input_path.write_text(input_text, encoding="utf-8")

    map_path = WORKFLOW_DIR / "v2b_atom_role_constraint_map.csv"
    write_csv(map_path, map_rows)

    state["v2a_geometry_promoted"] = True
    state["v2b_input_generated"] = True
    state["v2b_input_file"] = str(input_path.relative_to(ROOT))
    state["v2b_input_sha256"] = sha256(input_path)
    state["v2b_fixed_atom_count"] = len(fixed_indices)
    state["v2b_mobile_atom_count"] = len(mobile_indices)
    state["v2b_fixed_indices"] = fixed_indices
    state["v2b_mobile_indices"] = mobile_indices
    state["v2b_executed"] = False
    state["v2b_validation_pass"] = False
    state["qm_execution_authorized"] = False

    STATE_PATH.write_text(
        json.dumps(state, indent=2) + "\n",
        encoding="utf-8",
    )

    report_path = (
        WORKFLOW_DIR
        / "QM_F06_LOWER_BOUNDARY_V2B_PREPARATION_REPORT.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 LOWER Boundary V2-B Preparation — Day028",
                "",
                "## Starting geometry",
                "",
                (
                    f"- V2-A optimized XYZ: "
                    f"`{optimized_xyz_path.relative_to(ROOT)}`"
                ),
                "",
                "## Mobility model",
                "",
                "- Total atoms: **28**",
                "- Fixed atoms: **10**",
                "- Mobile atoms: **18**",
                (
                    "- Fixed: validated bridge/scaffold heavy-atom core"
                ),
                (
                    "- Mobile: all 15 H atoms and three restored "
                    "real heavy atoms"
                ),
                "",
                "## Scientific purpose",
                "",
                (
                    "Resolve the single residual peripheral contact "
                    "`HCAP:LOWER:03 — A:LOWER:13:-3` by allowing "
                    "the complete hydrogen boundary and restored "
                    "real coordination region to relax jointly."
                ),
                "",
                "## Authorization state",
                "",
                "- V2-B input preparation: **COMPLETED**",
                "- V2-B execution: **NOT AUTHORIZED**",
                "- ESP/RESP calculation: **NOT AUTHORIZED**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print("QM_F06 LOWER Boundary V2-B input prepared.")
    print("Input:", input_path)
    print("Fixed atoms:", len(fixed_indices))
    print("Mobile atoms:", len(mobile_indices))
    print("Input SHA256:", sha256(input_path))
    print("QM execution authorized: False")


if __name__ == "__main__":
    main()
