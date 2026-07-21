#!/usr/bin/env python3
"""
Finalize the rigid mapping of canonical refined R2 coordinates into the
QM_F06 UPPER V3-A2 frame using only fixed shared real atoms.

Mobile V3-A2 atoms are explicitly excluded because they underwent local
relaxation and do not represent the original canonical R2 geometry.

The script:
- selects fixed, shared, non-artificial atoms;
- computes a proper Kabsch transformation;
- validates the fixed-anchor residuals;
- transforms the complete canonical coordinate table;
- exports the 12 selected V4 restoration atoms;
- verifies restored B-N and real B-H distances;
- authorizes coordinate use for V4 construction only;
- does not construct V4 and does not authorize ORCA.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

SOURCE_COORDINATES = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "28_r2_inner_h_reflected_direction_refinement/"
    "r2_selected_four_atom_refined_full_coordinates.csv"
)

V3_ROOT = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V3/orca_v3a2_workflow"
)

V3_XYZ = V3_ROOT / "v3a2_start.xyz"
V3_MAP = V3_ROOT / "v3a2_atom_role_constraint_map.csv"

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day030_qm_f06_upper_v4_fixed_anchor"
)

REPORT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_FIXED_ANCHOR_REPORT.json"
)

ANCHOR_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_fixed_anchor_residuals.csv"
)

TRANSFORMED_SELECTED_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_selected_transformed_coordinates.csv"
)

RESTORED_BONDS_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_fixed_anchor_restored_bonds.csv"
)

DEFECTIVE_CAPS = {
    "HCAP:UPPER:01",
    "HCAP:UPPER:04",
}

SELECTED_V4_ATOMS = {
    "P:1640",
    "P:1581",
    "P:1583",
    "S:1739",
    "P:1639",
    "H4:UPPER:0203:0",
    "P:1580",
    "P:1582",
    "P:1638",
    "P:1642",
    "S:1738",
    "P:1523",
}

RESTORED_BN_EDGES = [
    ("P:1641", "P:1640"),
    ("P:1641", "S:1739"),
    ("P:1640", "P:1581"),
    ("P:1640", "P:1583"),
    ("S:1739", "P:1639"),
    ("P:1581", "P:1580"),
    ("P:1581", "P:1638"),
    ("P:1583", "P:1582"),
    ("P:1583", "P:1642"),
    ("P:1639", "P:1638"),
    ("P:1639", "S:1738"),
    ("P:1580", "P:1523"),
    ("P:1582", "P:1523"),
]

REAL_BH_EDGE = (
    "S:1739",
    "H4:UPPER:0203:0",
)

RMSD_TOLERANCE_A = 5.0e-5
MAX_RESIDUAL_TOLERANCE_A = 5.0e-5

BN_MIN_A = 1.35
BN_MAX_A = 1.60
BH_MIN_A = 1.10
BH_MAX_A = 1.30


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def read_xyz(path: Path) -> list[dict]:
    require_file(path)

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    count = int(lines[0].strip())

    atoms = []

    for index, line in enumerate(lines[2:2 + count]):
        fields = line.split()

        atoms.append({
            "index": index,
            "element": fields[0],
            "xyz_A": np.array(
                [float(value) for value in fields[1:4]],
                dtype=float,
            ),
        })

    if len(atoms) != count:
        raise RuntimeError(
            f"Incomplete XYZ: expected {count}, found {len(atoms)}"
        )

    return atoms


def kabsch(
    source_matrix: np.ndarray,
    target_matrix: np.ndarray,
):
    source_centroid = source_matrix.mean(axis=0)
    target_centroid = target_matrix.mean(axis=0)

    source_centered = source_matrix - source_centroid
    target_centered = target_matrix - target_centroid

    covariance = source_centered.T @ target_centered

    u_matrix, singular_values, vt_matrix = np.linalg.svd(
        covariance
    )

    rotation = vt_matrix.T @ u_matrix.T
    reflection_corrected = False

    if np.linalg.det(rotation) < 0.0:
        vt_matrix[-1, :] *= -1.0
        rotation = vt_matrix.T @ u_matrix.T
        reflection_corrected = True

    translation = (
        target_centroid
        - rotation @ source_centroid
    )

    return (
        rotation,
        translation,
        singular_values,
        reflection_corrected,
    )


def json_default(value):
    if isinstance(value, np.bool_):
        return bool(value)

    if isinstance(value, np.integer):
        return int(value)

    if isinstance(value, np.floating):
        return float(value)

    if isinstance(value, np.ndarray):
        return value.tolist()

    raise TypeError(
        f"Unsupported JSON type: {type(value).__name__}"
    )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    source_rows = read_csv(SOURCE_COORDINATES)
    map_rows = read_csv(V3_MAP)
    xyz_atoms = read_xyz(V3_XYZ)

    if len(map_rows) != len(xyz_atoms):
        raise RuntimeError(
            "V3 map and XYZ atom counts differ."
        )

    source = {}

    for row in source_rows:
        atom_id = row["node_id"]

        if atom_id in source:
            raise RuntimeError(
                f"Duplicate source atom ID: {atom_id}"
            )

        source[atom_id] = {
            "element": row["element"],
            "node_type": row["node_type"],
            "coordinate_source": row["coordinate_source"],
            "xyz_A": np.array(
                [
                    float(row["x_nm"]) * 10.0,
                    float(row["y_nm"]) * 10.0,
                    float(row["z_nm"]) * 10.0,
                ],
                dtype=float,
            ),
        }

    missing_selected = SELECTED_V4_ATOMS - set(source)

    if missing_selected:
        raise RuntimeError(
            "Selected V4 atoms missing from source: "
            f"{sorted(missing_selected)}"
        )

    v3 = {}

    for atom, row in zip(
        xyz_atoms,
        map_rows,
        strict=True,
    ):
        atom_id = row["atom_id"]

        v3[atom_id] = {
            "element": atom["element"],
            "xyz_A": atom["xyz_A"],
            "fixed": (
                row["v3a2_fixed"].strip().lower()
                == "true"
            ),
            "artificial_cap": (
                row["artificial_cap"].strip().lower()
                == "true"
            ),
        }

    anchor_ids = sorted(
        atom_id
        for atom_id, record in v3.items()
        if atom_id in source
        and atom_id not in DEFECTIVE_CAPS
        and record["fixed"]
        and not record["artificial_cap"]
    )

    if len(anchor_ids) != 15:
        raise RuntimeError(
            "Expected 15 fixed shared real anchor atoms; "
            f"found {len(anchor_ids)}: {anchor_ids}"
        )

    source_matrix = np.vstack([
        source[atom_id]["xyz_A"]
        for atom_id in anchor_ids
    ])

    target_matrix = np.vstack([
        v3[atom_id]["xyz_A"]
        for atom_id in anchor_ids
    ])

    (
        rotation,
        translation,
        singular_values,
        reflection_corrected,
    ) = kabsch(source_matrix, target_matrix)

    def transform(xyz: np.ndarray) -> np.ndarray:
        return rotation @ xyz + translation

    transformed_anchor = np.vstack([
        transform(source[atom_id]["xyz_A"])
        for atom_id in anchor_ids
    ])

    residuals = np.linalg.norm(
        transformed_anchor - target_matrix,
        axis=1,
    )

    rmsd_A = float(math.sqrt(
        np.mean(residuals ** 2)
    ))

    maximum_residual_A = float(residuals.max())

    proper_rotation = (
        abs(float(np.linalg.det(rotation)) - 1.0)
        <= 1.0e-10
    )

    anchor_pass = (
        rmsd_A <= RMSD_TOLERANCE_A
        and maximum_residual_A
        <= MAX_RESIDUAL_TOLERANCE_A
        and proper_rotation
    )

    anchor_records = []

    for index, atom_id in enumerate(anchor_ids):
        transformed = transformed_anchor[index]
        target = target_matrix[index]

        anchor_records.append({
            "atom_id": atom_id,
            "element": source[atom_id]["element"],
            "transformed_x_A": transformed[0],
            "transformed_y_A": transformed[1],
            "transformed_z_A": transformed[2],
            "v3_x_A": target[0],
            "v3_y_A": target[1],
            "v3_z_A": target[2],
            "residual_A": residuals[index],
            "within_tolerance": (
                residuals[index]
                <= MAX_RESIDUAL_TOLERANCE_A
            ),
        })

    with ANCHOR_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(anchor_records[0]),
        )
        writer.writeheader()
        writer.writerows(anchor_records)

    transformed_all = {
        atom_id: transform(record["xyz_A"])
        for atom_id, record in source.items()
    }

    selected_rows = []

    for atom_id in sorted(SELECTED_V4_ATOMS):
        source_record = source[atom_id]
        transformed = transformed_all[atom_id]

        selected_rows.append({
            "atom_id": atom_id,
            "element": source_record["element"],
            "node_type": source_record["node_type"],
            "coordinate_source": (
                source_record["coordinate_source"]
            ),
            "transformed_x_A": transformed[0],
            "transformed_y_A": transformed[1],
            "transformed_z_A": transformed[2],
            "selected_for_v4_restoration": True,
        })

    with TRANSFORMED_SELECTED_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(selected_rows[0]),
        )
        writer.writeheader()
        writer.writerows(selected_rows)

    bond_records = []

    for first, second in RESTORED_BN_EDGES:
        value = float(np.linalg.norm(
            transformed_all[first]
            - transformed_all[second]
        ))

        bond_records.append({
            "first_atom": first,
            "second_atom": second,
            "bond_class": "B-N",
            "distance_A": value,
            "minimum_A": BN_MIN_A,
            "maximum_A": BN_MAX_A,
            "pass": BN_MIN_A <= value <= BN_MAX_A,
        })

    first, second = REAL_BH_EDGE

    value = float(np.linalg.norm(
        transformed_all[first]
        - transformed_all[second]
    ))

    bond_records.append({
        "first_atom": first,
        "second_atom": second,
        "bond_class": "REAL_R2_B-H",
        "distance_A": value,
        "minimum_A": BH_MIN_A,
        "maximum_A": BH_MAX_A,
        "pass": BH_MIN_A <= value <= BH_MAX_A,
    })

    with RESTORED_BONDS_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(bond_records[0]),
        )
        writer.writeheader()
        writer.writerows(bond_records)

    restored_bonds_pass = all(
        row["pass"]
        for row in bond_records
    )

    overall_pass = (
        anchor_pass
        and restored_bonds_pass
        and len(selected_rows) == 12
    )

    decision = (
        "QM_F06_UPPER_V4_FIXED_ANCHOR_PASS_"
        "TRANSFORMED_RESTORATION_COORDINATES_AUTHORIZED"
        if overall_pass
        else
        "QM_F06_UPPER_V4_FIXED_ANCHOR_FAIL_"
        "CONSTRUCTION_BLOCKED"
    )

    report = {
        "decision": decision,
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "anchor_selection": (
            "FIXED_SHARED_REAL"
        ),
        "anchor_atom_count": len(anchor_ids),
        "anchor_atom_ids": anchor_ids,
        "rmsd_A": rmsd_A,
        "maximum_residual_A": maximum_residual_A,
        "rmsd_tolerance_A": RMSD_TOLERANCE_A,
        "maximum_residual_tolerance_A": (
            MAX_RESIDUAL_TOLERANCE_A
        ),
        "rotation_matrix": rotation,
        "rotation_determinant": np.linalg.det(rotation),
        "translation_vector_A": translation,
        "singular_values": singular_values,
        "reflection_corrected": reflection_corrected,
        "proper_rotation": proper_rotation,
        "anchor_gate_pass": anchor_pass,
        "selected_v4_atom_count": len(selected_rows),
        "restored_bond_gate_pass": restored_bonds_pass,
        "overall_pass": overall_pass,
        "files": {
            "anchor_residuals_csv": str(
                ANCHOR_CSV.relative_to(ROOT)
            ),
            "selected_transformed_coordinates_csv": str(
                TRANSFORMED_SELECTED_CSV.relative_to(ROOT)
            ),
            "restored_bonds_csv": str(
                RESTORED_BONDS_CSV.relative_to(ROOT)
            ),
        },
        "files_sha256": {
            "source_coordinates": sha256(
                SOURCE_COORDINATES
            ),
            "v3_xyz": sha256(V3_XYZ),
            "v3_map": sha256(V3_MAP),
            "anchor_residuals_csv": sha256(
                ANCHOR_CSV
            ),
            "selected_transformed_coordinates_csv": sha256(
                TRANSFORMED_SELECTED_CSV
            ),
            "restored_bonds_csv": sha256(
                RESTORED_BONDS_CSV
            ),
        },
        "authorization": {
            "fixed_anchor_selected": overall_pass,
            "transformed_restoration_coordinates_authorized": (
                overall_pass
            ),
            "v4_geometry_construction_authorized": (
                overall_pass
            ),
            "orca_execution_authorized": False,
            "geometric_reference_accepted": False,
            "electronic_reference_accepted": False,
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    REPORT_PATH.write_text(
        json.dumps(
            report,
            indent=2,
            default=json_default,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 78)
    print("QM_F06 UPPER V4 FIXED-ANCHOR FINALIZATION")
    print("=" * 78)
    print("Anchor:", "FIXED_SHARED_REAL")
    print("Anchor atoms:", len(anchor_ids))
    print("RMSD A:", rmsd_A)
    print("Maximum residual A:", maximum_residual_A)
    print("Rotation determinant:", np.linalg.det(rotation))
    print(
        "Anchor gate:",
        "PASS" if anchor_pass else "FAIL",
    )
    print(
        "Restored bonds:",
        "PASS" if restored_bonds_pass else "FAIL",
    )
    print("Selected transformed atoms:", len(selected_rows))
    print()
    print("Decision:", decision)
    print("Report:", REPORT_PATH)
    print("Anchor residuals:", ANCHOR_CSV)
    print(
        "Selected transformed coordinates:",
        TRANSFORMED_SELECTED_CSV,
    )
    print("Restored bonds:", RESTORED_BONDS_CSV)
    print(
        "V4 construction authorized:",
        overall_pass,
    )
    print("ORCA execution authorized: False")


if __name__ == "__main__":
    main()
