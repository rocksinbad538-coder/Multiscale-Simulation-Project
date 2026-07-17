#!/usr/bin/env python3
"""
Prepare a non-optimizing electronic diagnostic single point for the
accepted QM_F06 LOWER Boundary V2-B geometry.

Purpose
-------
Evaluate the residual real 1-5 contact:

    BR4:LOWER:00:3 ... H4:LOWER:0017:0

Requested analyses:
- Mayer bond orders and atomic valences;
- Hirshfeld charges;
- MBIS charges;
- CHELPG charges.

These analyses are diagnostic only. No ESP/RESP or force-field charge
assignment is authorized by this script.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

V2_DIR = ROOT / (
    "runs/phase1A/day027_qm_f06_lower_boundary_redesign/"
    "QM_F06_LOWER_BOUNDARY_V2"
)

STATE_PATH = V2_DIR / (
    "orca_v2_workflow/v2_workflow_state.json"
)

ATOM_MANIFEST = V2_DIR / (
    "cap02_repair/"
    "QM_F06_LOWER_BOUNDARY_V2_REPAIRED_atoms.csv"
)

STRUCTURAL_SUMMARY = ROOT / (
    "runs/phase1A/"
    "day028_qm_f06_lower_boundary_v2b_postprocessing/"
    "QM_F06_LOWER_BOUNDARY_V2B_validation_summary.json"
)

CONTACT_SUMMARY = ROOT / (
    "runs/phase1A/"
    "day028_qm_f06_lower_boundary_v2b_postprocessing/"
    "residual_real_contact_audit/"
    "residual_real_contact_audit_summary.json"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day028_qm_f06_lower_boundary_v2b_electronic_diagnostic"
)

INPUT_PATH = OUTPUT_DIR / (
    "QM_F06_LOWER_BOUNDARY_V2B_ELECTRONIC_DIAGNOSTIC.inp"
)

TARGET_B = "BR4:LOWER:00:3"
TARGET_H = "H4:LOWER:0017:0"


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


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


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for path in (
        STATE_PATH,
        ATOM_MANIFEST,
        STRUCTURAL_SUMMARY,
        CONTACT_SUMMARY,
    ):
        require_file(path)

    state = json.loads(
        STATE_PATH.read_text(encoding="utf-8")
    )

    structural = json.loads(
        STRUCTURAL_SUMMARY.read_text(encoding="utf-8")
    )

    contact = json.loads(
        CONTACT_SUMMARY.read_text(encoding="utf-8")
    )

    readiness_checks = {
        "v2b_executed": state.get("v2b_executed") is True,
        "v2b_validation_pass": (
            state.get("v2b_validation_pass") is True
        ),
        "v2b_xyz_present": bool(
            state.get("v2b_optimized_xyz")
        ),
        "structural_gate_pass": (
            structural.get("v2b_structural_gate_pass") is True
        ),
        "geometry_acceptance_retained": (
            contact.get("geometry_acceptance_retained") is True
        ),
        "electronic_protocol_definition_authorized": (
            contact.get(
                "electronic_property_protocol_definition_authorized"
            )
            is True
        ),
        "target_pair_matches": (
            contact.get("bridge_atom") == TARGET_B
            and contact.get("hydrogen_atom") == TARGET_H
        ),
        "no_unintended_covalent_contact": (
            contact.get(
                "possible_unintended_covalent_contact"
            )
            is False
        ),
    }

    if not all(readiness_checks.values()):
        raise RuntimeError(
            "Electronic diagnostic readiness failed: "
            f"{readiness_checks}"
        )

    optimized_xyz_path = ROOT / state["v2b_optimized_xyz"]
    xyz_rows = read_xyz(optimized_xyz_path)

    with ATOM_MANIFEST.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        atoms = list(csv.DictReader(handle))

    if len(atoms) != 28 or len(xyz_rows) != 28:
        raise RuntimeError("Expected 28 atoms.")

    atom_ids = [row["atom_id"] for row in atoms]

    for index, (atom, xyz) in enumerate(
        zip(atoms, xyz_rows, strict=True)
    ):
        if atom["element"] != xyz[0]:
            raise RuntimeError(
                f"Element mismatch at index {index}: "
                f"{atom['atom_id']}"
            )

    target_b_index = atom_ids.index(TARGET_B)
    target_h_index = atom_ids.index(TARGET_H)

    coordinate_lines = [
        f"{element:<2s} {x: .10f} {y: .10f} {z: .10f}"
        for element, x, y, z in xyz_rows
    ]

    input_text = "\n".join(
        [
            (
                "! PBE0 D4 def2-TZVP def2/J RIJCOSX "
                "TightSCF DefGrid3 MAYER HIRSHFELD MBIS CHELPG"
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
            "%method",
            "  MAYER_BONDORDERTHRESH 0.01",
            "end",
            "",
            "%chelpg",
            "  GRID 0.3",
            "  RMAX 2.8",
            "end",
            "",
            "%output",
            "  Print[ P_Mayer ] 1",
            "  Print[ P_Hirshfeld ] 1",
            "  Print[ P_MBIS ] 1",
            "end",
            "",
            (
                "# Diagnostic electronic single point on accepted "
                "QM_F06 LOWER V2-B geometry"
            ),
            (
                f"# Target contact: atom {target_b_index} {TARGET_B} "
                f"... atom {target_h_index} {TARGET_H}"
            ),
            "# No geometry optimization; no RESP fitting",
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

    atom_map_path = OUTPUT_DIR / (
        "QM_F06_LOWER_BOUNDARY_V2B_atom_index_map.csv"
    )

    with atom_map_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        fields = [
            "index_0based",
            "atom_id",
            "element",
            "atom_role",
            "node_type",
            "target_contact_atom",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
        )
        writer.writeheader()

        for index, atom in enumerate(atoms):
            writer.writerow(
                {
                    "index_0based": index,
                    "atom_id": atom["atom_id"],
                    "element": atom["element"],
                    "atom_role": atom["atom_role"],
                    "node_type": atom["node_type"],
                    "target_contact_atom": (
                        atom["atom_id"] in {TARGET_B, TARGET_H}
                    ),
                }
            )

    summary = {
        "decision": (
            "QM_F06_LOWER_V2B_ELECTRONIC_DIAGNOSTIC_INPUT_PREPARED"
        ),
        "geometry_source": str(
            optimized_xyz_path.relative_to(ROOT)
        ),
        "atom_count": len(atoms),
        "charge": 0,
        "multiplicity": 1,
        "target_b_atom_id": TARGET_B,
        "target_b_index_0based": target_b_index,
        "target_h_atom_id": TARGET_H,
        "target_h_index_0based": target_h_index,
        "method": "PBE0-D4/def2-TZVP",
        "requested_analyses": [
            "MAYER",
            "HIRSHFELD",
            "MBIS",
            "CHELPG",
        ],
        "input_file": str(INPUT_PATH.relative_to(ROOT)),
        "input_sha256": sha256(INPUT_PATH),
        "geometry_optimization_requested": False,
        "electronic_diagnostic_execution_authorized": False,
        "esp_resp_parameter_adoption_authorized": False,
        "force_field_parameter_adoption_authorized": False,
    }

    (
        OUTPUT_DIR
        / "QM_F06_LOWER_BOUNDARY_V2B_electronic_diagnostic_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    report_path = OUTPUT_DIR / (
        "QM_F06_LOWER_BOUNDARY_V2B_"
        "ELECTRONIC_DIAGNOSTIC_PREPARATION.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 LOWER Boundary V2-B Electronic Diagnostic",
                "",
                "## Geometry",
                "",
                (
                    f"- Source: "
                    f"`{optimized_xyz_path.relative_to(ROOT)}`"
                ),
                "- Atoms: **28**",
                "- Charge/multiplicity: **0 / 1**",
                "- Geometry optimization: **NO**",
                "",
                "## Target contact",
                "",
                f"- B atom: `{TARGET_B}` — index `{target_b_index}`",
                f"- H atom: `{TARGET_H}` — index `{target_h_index}`",
                "",
                "## Electronic analyses",
                "",
                "- Mayer bond orders and valences",
                "- Hirshfeld charges",
                "- MBIS charges",
                "- CHELPG charges",
                "",
                "## Interpretation boundary",
                "",
                (
                    "The requested population analyses are diagnostic. "
                    "They do not authorize RESP, charge adoption or "
                    "force-field parameter assignment."
                ),
                "",
                "## Authorization state",
                "",
                "- Input preparation: **COMPLETED**",
                "- Diagnostic single-point execution: **NOT AUTHORIZED**",
                "- RESP/ESP charge adoption: **NOT AUTHORIZED**",
                "- Force-field parameter adoption: **NOT AUTHORIZED**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print("Electronic diagnostic input prepared.")
    print("Input:", INPUT_PATH)
    print("Target B index:", target_b_index)
    print("Target H index:", target_h_index)
    print("Input SHA256:", sha256(INPUT_PATH))
    print("QM execution authorized: False")


if __name__ == "__main__":
    main()
