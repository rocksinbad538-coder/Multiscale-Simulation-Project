#!/usr/bin/env python3
"""
Rigidly align the canonical refined R2 geometry to the QM_F06 UPPER
V3-A2 coordinate frame using all shared real atoms.

The audit:
- computes the unaligned displacement statistics;
- performs a proper Kabsch rotation plus translation;
- reports RMSD and maximum residual after alignment;
- applies the same transformation to all mandatory V4 atoms;
- validates restored B-N and real B-H distances after transformation;
- writes transformed V4 source coordinates for the constructor;
- does not construct the final capped V4 geometry or authorize ORCA.
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

V3_WORKFLOW = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V3/orca_v3a2_workflow"
)

V3_XYZ = V3_WORKFLOW / "v3a2_start.xyz"

V3_MAP = (
    V3_WORKFLOW
    / "v3a2_atom_role_constraint_map.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day030_qm_f06_upper_v4_rigid_alignment"
)

REPORT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_RIGID_ALIGNMENT_AUDIT.json"
)

SHARED_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_shared_atom_aligned_residuals.csv"
)

TRANSFORMED_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_transformed_source_coordinates.csv"
)

BONDS_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_transformed_restored_bonds.csv"
)

DEFECTIVE_CAPS = {
    "HCAP:UPPER:01",
    "HCAP:UPPER:04",
}

MANDATORY_V4_ATOMS = {
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

RMSD_TOLERANCE_A = 1.0e-6
MAX_RESIDUAL_TOLERANCE_A = 5.0e-6

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

    atom_count = int(lines[0].strip())

    if len(lines) < atom_count + 2:
        raise RuntimeError(f"Incomplete XYZ: {path}")

    atoms = []

    for index, line in enumerate(
        lines[2:2 + atom_count]
    ):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Malformed XYZ line {index}: {line}"
            )

        atoms.append({
            "index": index,
            "element": fields[0],
            "xyz_A": np.array(
                [float(value) for value in fields[1:4]],
                dtype=float,
            ),
        })

    return atoms


def euclidean(
    first: np.ndarray,
    second: np.ndarray,
) -> float:
    return float(np.linalg.norm(first - second))


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    source_rows = read_csv(SOURCE_COORDINATES)
    map_rows = read_csv(V3_MAP)
    v3_atoms = read_xyz(V3_XYZ)

    if len(map_rows) != len(v3_atoms):
        raise RuntimeError(
            "V3 map and XYZ counts differ."
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
            "coordinate_source": row[
                "coordinate_source"
            ],
            "xyz_A": np.array(
                [
                    float(row["x_nm"]) * 10.0,
                    float(row["y_nm"]) * 10.0,
                    float(row["z_nm"]) * 10.0,
                ],
                dtype=float,
            ),
        }

    missing = MANDATORY_V4_ATOMS - set(source)

    if missing:
        raise RuntimeError(
            f"Mandatory V4 atoms missing: {sorted(missing)}"
        )

    v3 = {}

    for atom, row in zip(
        v3_atoms,
        map_rows,
        strict=True,
    ):
        atom_id = row["atom_id"]

        v3[atom_id] = {
            "element": atom["element"],
            "xyz_A": atom["xyz_A"],
            "artificial_cap": (
                row["artificial_cap"]
                .strip()
                .lower()
                == "true"
            ),
        }

    retained_v3_ids = set(v3) - DEFECTIVE_CAPS

    shared_ids = sorted(
        atom_id
        for atom_id in retained_v3_ids
        if atom_id in source
        and not v3[atom_id]["artificial_cap"]
    )

    if len(shared_ids) < 3:
        raise RuntimeError(
            "At least three shared real atoms are required."
        )

    source_matrix = np.vstack([
        source[atom_id]["xyz_A"]
        for atom_id in shared_ids
    ])

    target_matrix = np.vstack([
        v3[atom_id]["xyz_A"]
        for atom_id in shared_ids
    ])

    source_centroid = source_matrix.mean(axis=0)
    target_centroid = target_matrix.mean(axis=0)

    source_centered = (
        source_matrix - source_centroid
    )

    target_centered = (
        target_matrix - target_centroid
    )

    covariance = (
        source_centered.T @ target_centered
    )

    u_matrix, singular_values, vt_matrix = (
        np.linalg.svd(covariance)
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

    def transform(xyz: np.ndarray) -> np.ndarray:
        return rotation @ xyz + translation

    unaligned_distances = np.linalg.norm(
        source_matrix - target_matrix,
        axis=1,
    )

    transformed_shared = np.vstack([
        transform(source[atom_id]["xyz_A"])
        for atom_id in shared_ids
    ])

    residual_vectors = (
        transformed_shared - target_matrix
    )

    residual_distances = np.linalg.norm(
        residual_vectors,
        axis=1,
    )

    rmsd_A = float(math.sqrt(
        np.mean(residual_distances ** 2)
    ))

    maximum_residual_A = float(
        residual_distances.max()
    )

    mean_unaligned_A = float(
        unaligned_distances.mean()
    )

    maximum_unaligned_A = float(
        unaligned_distances.max()
    )

    shared_records = []

    for index, atom_id in enumerate(shared_ids):
        transformed = transformed_shared[index]
        target = target_matrix[index]
        residual = residual_distances[index]

        shared_records.append({
            "atom_id": atom_id,
            "element": source[atom_id]["element"],
            "transformed_x_A": transformed[0],
            "transformed_y_A": transformed[1],
            "transformed_z_A": transformed[2],
            "v3_x_A": target[0],
            "v3_y_A": target[1],
            "v3_z_A": target[2],
            "residual_A": residual,
            "within_tolerance": (
                residual
                <= MAX_RESIDUAL_TOLERANCE_A
            ),
        })

    with SHARED_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(shared_records[0]),
        )
        writer.writeheader()
        writer.writerows(shared_records)

    transformed_source = {}

    transformed_rows = []

    for atom_id, record in sorted(source.items()):
        xyz = transform(record["xyz_A"])

        transformed_source[atom_id] = xyz

        transformed_rows.append({
            "atom_id": atom_id,
            "element": record["element"],
            "node_type": record["node_type"],
            "source_coordinate_source": (
                record["coordinate_source"]
            ),
            "source_x_A": record["xyz_A"][0],
            "source_y_A": record["xyz_A"][1],
            "source_z_A": record["xyz_A"][2],
            "transformed_x_A": xyz[0],
            "transformed_y_A": xyz[1],
            "transformed_z_A": xyz[2],
            "selected_for_v4_restoration": (
                atom_id in MANDATORY_V4_ATOMS
            ),
        })

    with TRANSFORMED_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(transformed_rows[0]),
        )
        writer.writeheader()
        writer.writerows(transformed_rows)

    bond_records = []

    for first, second in RESTORED_BN_EDGES:
        value = euclidean(
            transformed_source[first],
            transformed_source[second],
        )

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

    value = euclidean(
        transformed_source[first],
        transformed_source[second],
    )

    bond_records.append({
        "first_atom": first,
        "second_atom": second,
        "bond_class": "REAL_R2_B-H",
        "distance_A": value,
        "minimum_A": BH_MIN_A,
        "maximum_A": BH_MAX_A,
        "pass": BH_MIN_A <= value <= BH_MAX_A,
    })

    with BONDS_CSV.open(
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

    proper_rotation = (
        abs(np.linalg.det(rotation) - 1.0)
        <= 1.0e-10
    )

    alignment_pass = (
        rmsd_A <= RMSD_TOLERANCE_A
        and maximum_residual_A
        <= MAX_RESIDUAL_TOLERANCE_A
        and proper_rotation
        and all(
            row["within_tolerance"]
            for row in shared_records
        )
    )

    restored_bonds_pass = all(
        row["pass"]
        for row in bond_records
    )

    overall_pass = (
        alignment_pass
        and restored_bonds_pass
    )

    decision = (
        "QM_F06_UPPER_V4_RIGID_ALIGNMENT_PASS_"
        "TRANSFORMED_CONSTRUCTION_COORDINATES_AUTHORIZED"
        if overall_pass
        else
        "QM_F06_UPPER_V4_RIGID_ALIGNMENT_FAIL_"
        "CONSTRUCTION_BLOCKED"
    )

    report = {
        "decision": decision,
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "source_coordinates": str(
            SOURCE_COORDINATES.relative_to(ROOT)
        ),
        "target_v3_xyz": str(
            V3_XYZ.relative_to(ROOT)
        ),
        "shared_atom_count": len(shared_ids),
        "shared_atom_ids": shared_ids,
        "unaligned": {
            "mean_distance_A": mean_unaligned_A,
            "maximum_distance_A": (
                maximum_unaligned_A
            ),
        },
        "aligned": {
            "rmsd_A": rmsd_A,
            "maximum_residual_A": (
                maximum_residual_A
            ),
            "rmsd_tolerance_A": (
                RMSD_TOLERANCE_A
            ),
            "maximum_residual_tolerance_A": (
                MAX_RESIDUAL_TOLERANCE_A
            ),
        },
        "rotation_matrix": rotation.tolist(),
        "rotation_determinant": float(
            np.linalg.det(rotation)
        ),
        "translation_vector_A": (
            translation.tolist()
        ),
        "singular_values": (
            singular_values.tolist()
        ),
        "reflection_corrected": (
            reflection_corrected
        ),
        "proper_rotation": proper_rotation,
        "alignment_pass": alignment_pass,
        "restored_bond_gate_pass": (
            restored_bonds_pass
        ),
        "overall_pass": overall_pass,
        "files": {
            "shared_residuals_csv": str(
                SHARED_CSV.relative_to(ROOT)
            ),
            "transformed_source_coordinates_csv": str(
                TRANSFORMED_CSV.relative_to(ROOT)
            ),
            "transformed_restored_bonds_csv": str(
                BONDS_CSV.relative_to(ROOT)
            ),
        },
        "files_sha256": {
            "source_coordinates": sha256(
                SOURCE_COORDINATES
            ),
            "v3_xyz": sha256(V3_XYZ),
            "v3_map": sha256(V3_MAP),
            "shared_residuals_csv": sha256(
                SHARED_CSV
            ),
            "transformed_source_coordinates_csv": sha256(
                TRANSFORMED_CSV
            ),
            "transformed_restored_bonds_csv": sha256(
                BONDS_CSV
            ),
        },
        "authorization": {
            "transformed_v4_source_coordinates_authorized": (
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

    def json_compatible(value):
        """
        Convert NumPy scalar values to native Python objects
        for deterministic JSON serialization.
        """
        if isinstance(value, np.bool_):
            return bool(value)

        if isinstance(value, np.integer):
            return int(value)

        if isinstance(value, np.floating):
            return float(value)

        if isinstance(value, np.ndarray):
            return value.tolist()

        raise TypeError(
            "Unsupported JSON value type: "
            f"{type(value).__name__}"
        )

    REPORT_PATH.write_text(
        json.dumps(
            report,
            indent=2,
            default=json_compatible,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 78)
    print("QM_F06 UPPER V4 RIGID-ALIGNMENT AUDIT")
    print("=" * 78)
    print("Shared atoms:", len(shared_ids))
    print(
        "Unaligned mean distance A:",
        mean_unaligned_A,
    )
    print(
        "Unaligned maximum distance A:",
        maximum_unaligned_A,
    )
    print("Aligned RMSD A:", rmsd_A)
    print(
        "Aligned maximum residual A:",
        maximum_residual_A,
    )
    print(
        "Rotation determinant:",
        np.linalg.det(rotation),
    )
    print(
        "Reflection corrected:",
        reflection_corrected,
    )
    print(
        "Alignment gate:",
        "PASS" if alignment_pass else "FAIL",
    )
    print(
        "Restored bonds:",
        "PASS"
        if restored_bonds_pass
        else "FAIL",
    )
    print()
    print("Rotation matrix:")
    print(rotation)
    print()
    print("Translation vector A:")
    print(translation)
    print()
    print("Decision:", decision)
    print("Report:", REPORT_PATH)
    print("Shared residuals:", SHARED_CSV)
    print("Transformed source:", TRANSFORMED_CSV)
    print("Restored bonds:", BONDS_CSV)
    print(
        "V4 construction authorized:",
        overall_pass,
    )
    print("ORCA execution authorized: False")


if __name__ == "__main__":
    main()
