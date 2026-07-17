#!/usr/bin/env python3
"""
Prepare the constrained QM_F06 UPPER Boundary V3-A ORCA optimization.

Mobility model
--------------
Fixed:
- retained atoms from the original repaired UPPER fragment;
- original bridge/scaffold core;
- original retained artificial caps.

Mobile:
- all real R2 atoms restored during Boundary V2 and V3;
- all artificial caps introduced during Boundary V2 and V3.

No QM calculation is executed.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

V3_DIR = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V3"
)

ATOMS_PATH = V3_DIR / "QM_F06_UPPER_BOUNDARY_V3_atoms.csv"
XYZ_PATH = V3_DIR / "QM_F06_UPPER_BOUNDARY_V3.xyz"

PRE_QM_SUMMARY = V3_DIR / (
    "pre_qm_audit/"
    "QM_F06_UPPER_BOUNDARY_V3_pre_qm_summary.json"
)

WORKFLOW_DIR = V3_DIR / "orca_v3_workflow"

INPUT_PATH = WORKFLOW_DIR / "v3a_boundary_relax.inp"
MAP_PATH = WORKFLOW_DIR / "v3a_atom_role_constraint_map.csv"
STATE_PATH = WORKFLOW_DIR / "v3_workflow_state.json"
REPORT_PATH = WORKFLOW_DIR / (
    "QM_F06_UPPER_BOUNDARY_V3A_PREPARATION.md"
)

MOBILE_ROLE_PREFIXES = (
    "REAL_R2_BOUNDARY_EXPANSION",
    "REAL_R2_BOUNDARY_V3_EXPANSION",
    "ARTIFICIAL_BOUNDARY_CAP_V2",
    "ARTIFICIAL_BOUNDARY_CAP_V3",
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


def read_xyz(
    path: Path,
) -> list[tuple[str, float, float, float]]:
    require_file(path)

    lines = path.read_text(encoding="utf-8").splitlines()
    expected = int(lines[0].strip())

    rows = [
        (
            fields[0],
            float(fields[1]),
            float(fields[2]),
            float(fields[3]),
        )
        for line in lines[2:]
        if line.strip()
        for fields in [line.split()]
    ]

    if len(rows) != expected:
        raise RuntimeError(
            f"XYZ atom-count mismatch: {expected} vs {len(rows)}"
        )

    return rows


def is_v3a_mobile(role: str) -> bool:
    return role.startswith(MOBILE_ROLE_PREFIXES)


def main() -> None:
    WORKFLOW_DIR.mkdir(parents=True, exist_ok=True)

    for path in (
        ATOMS_PATH,
        XYZ_PATH,
        PRE_QM_SUMMARY,
    ):
        require_file(path)

    pre_qm = json.loads(
        PRE_QM_SUMMARY.read_text(encoding="utf-8")
    )

    if pre_qm.get("pre_qm_gate_pass") is not True:
        raise RuntimeError(
            "UPPER Boundary V3 has not passed the pre-QM gate."
        )

    with ATOMS_PATH.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        atoms = list(csv.DictReader(handle))

    xyz_rows = read_xyz(XYZ_PATH)

    if len(atoms) != 30 or len(xyz_rows) != 30:
        raise RuntimeError("Expected 30 UPPER V3 atoms.")

    map_rows = []
    fixed_indices = []
    mobile_indices = []

    for index, (atom, xyz) in enumerate(
        zip(atoms, xyz_rows, strict=True)
    ):
        if atom["element"] != xyz[0]:
            raise RuntimeError(
                f"Element mismatch at index {index}: "
                f"{atom['atom_id']}"
            )

        mobile = is_v3a_mobile(atom["atom_role"])
        fixed = not mobile

        if fixed:
            fixed_indices.append(index)
        else:
            mobile_indices.append(index)

        map_rows.append(
            {
                "index_0based": index,
                "atom_id": atom["atom_id"],
                "element": atom["element"],
                "atom_role": atom["atom_role"],
                "node_type": atom["node_type"],
                "artificial_cap": atom["artificial_cap"],
                "v3a_fixed": fixed,
                "v3a_mobile": mobile,
                "mobility_basis": (
                    "RESTORED_BOUNDARY_OR_NEW_CAP"
                    if mobile
                    else "RETAINED_VALIDATED_CORE"
                ),
            }
        )

    if len(fixed_indices) != 20:
        raise RuntimeError(
            f"Expected 20 fixed atoms; found {len(fixed_indices)}"
        )

    if len(mobile_indices) != 10:
        raise RuntimeError(
            f"Expected 10 mobile atoms; found {len(mobile_indices)}"
        )

    constraint_lines = [
        f"    {{ C {index} C }}"
        for index in fixed_indices
    ]

    coordinate_lines = [
        f"{element:<2s} {x: .10f} {y: .10f} {z: .10f}"
        for element, x, y, z in xyz_rows
    ]

    input_text = "\n".join(
        [
            (
                "! PBE0 D4 def2-TZVP def2/J RIJCOSX "
                "TightSCF Opt DefGrid3"
            ),
            "",
            "%pal",
            "  nprocs 8",
            "end",
            "",
            "%maxcore 3000",
            "",
            "%scf",
            "  MaxIter 500",
            "end",
            "",
            "%geom",
            "  MaxIter 250",
            "  Constraints",
            *constraint_lines,
            "  end",
            "end",
            "",
            (
                "# QM_F06 UPPER Boundary V3-A constrained "
                "boundary relaxation"
            ),
            (
                "# Fixed: retained validated core; "
                "mobile: restored V2/V3 region and V2/V3 caps"
            ),
            "* xyz 0 1",
            *coordinate_lines,
            "*",
            "",
        ]
    )

    INPUT_PATH.write_text(
        input_text,
        encoding="utf-8",
    )

    with MAP_PATH.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(map_rows[0]),
        )
        writer.writeheader()
        writer.writerows(map_rows)

    state = {
        "fragment": "QM_F06_UPPER_BOUNDARY_V3",
        "formula": "B8N7H15",
        "charge": 0,
        "multiplicity": 1,
        "pre_qm_gate_pass": True,
        "v3a_input_prepared": True,
        "v3a_input": str(INPUT_PATH.relative_to(ROOT)),
        "v3a_input_sha256": sha256(INPUT_PATH),
        "v3a_nprocs": 8,
        "v3a_maxcore_mb_per_process": 3000,
        "v3a_fixed_atom_count": len(fixed_indices),
        "v3a_mobile_atom_count": len(mobile_indices),
        "v3a_fixed_indices": fixed_indices,
        "v3a_mobile_indices": mobile_indices,
        "v3a_executed": False,
        "v3a_validation_pass": False,
        "v3b_input_prepared": False,
        "v3b_executed": False,
        "electronic_diagnostic_executed": False,
        "qm_execution_authorized": False,
        "esp_resp_execution_authorized": False,
        "force_field_parameter_adoption_authorized": False,
    }

    STATE_PATH.write_text(
        json.dumps(state, indent=2) + "\n",
        encoding="utf-8",
    )

    REPORT_PATH.write_text(
        "\n".join(
            [
                "# QM_F06 UPPER Boundary V3-A Preparation — Day028",
                "",
                "## Electronic-structure model",
                "",
                "- Method: **PBE0-D4/def2-TZVP**",
                "- Approximation: **RIJCOSX/def2-J**",
                "- Charge/multiplicity: **0 / 1**",
                "- ORCA processes: **8**",
                "- MaxCore: **3000 MB per process**",
                "",
                "## Mobility model",
                "",
                "- Total atoms: **30**",
                "- Fixed atoms: **20**",
                "- Mobile atoms: **10**",
                (
                    "- Mobile region: restored V2/V3 real atoms "
                    "and V2/V3 artificial caps"
                ),
                "",
                "## Purpose",
                "",
                (
                    "Relax the independently reconstructed UPPER "
                    "boundary while retaining the previously validated "
                    "bridge/scaffold core."
                ),
                "",
                "## Authorization state",
                "",
                "- Input preparation: **COMPLETED**",
                "- Static audit: **PENDING**",
                "- ORCA execution: **NOT AUTHORIZED**",
                "- ESP/RESP execution: **NOT AUTHORIZED**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print("QM_F06 UPPER Boundary V3-A input prepared.")
    print("Input:", INPUT_PATH)
    print("Fixed atoms:", len(fixed_indices))
    print("Mobile atoms:", len(mobile_indices))
    print("ORCA nprocs: 8")
    print("Input SHA256:", sha256(INPUT_PATH))
    print("QM execution authorized: False")


if __name__ == "__main__":
    main()
