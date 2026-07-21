#!/usr/bin/env python3
"""
Compare rigid-alignment anchor subsets for QM_F06 UPPER V4.

The purpose is to determine whether:
- the fixed retained core remains exactly compatible with canonical R2;
- mobile/shared atoms account for the global non-rigid discrepancy;
- a chemically local fixed anchor can place the restored branch.

No construction or ORCA execution is authorized.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

SOURCE = ROOT / (
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
    "day030_qm_f06_upper_v4_alignment_subsets"
)

REPORT = OUTPUT_DIR / (
    "QM_F06_UPPER_V4_ALIGNMENT_SUBSET_AUDIT.json"
)

CSV_PATH = OUTPUT_DIR / (
    "QM_F06_UPPER_V4_alignment_subset_metrics.csv"
)

DEFECTIVE_CAPS = {
    "HCAP:UPPER:01",
    "HCAP:UPPER:04",
}

LOCAL_REQUIRED_IDS = {
    "P:1641",
    "S:1710",
}


def read_csv(path):
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def read_xyz(path):
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    count = int(lines[0])
    atoms = []

    for index, line in enumerate(lines[2:2 + count]):
        fields = line.split()
        atoms.append({
            "index": index,
            "element": fields[0],
            "xyz": np.array(
                [float(value) for value in fields[1:4]],
                dtype=float,
            ),
        })

    return atoms


def kabsch(source_matrix, target_matrix):
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

    if np.linalg.det(rotation) < 0:
        vt_matrix[-1, :] *= -1
        rotation = vt_matrix.T @ u_matrix.T
        reflection_corrected = True

    translation = target_centroid - rotation @ source_centroid

    transformed = (
        source_matrix @ rotation.T
        + translation
    )

    residuals = np.linalg.norm(
        transformed - target_matrix,
        axis=1,
    )

    return {
        "rotation": rotation,
        "translation": translation,
        "singular_values": singular_values,
        "reflection_corrected": reflection_corrected,
        "residuals": residuals,
        "rmsd_A": float(math.sqrt(
            np.mean(residuals ** 2)
        )),
        "maximum_residual_A": float(
            residuals.max()
        ),
        "determinant": float(
            np.linalg.det(rotation)
        ),
    }


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    source_rows = read_csv(SOURCE)
    map_rows = read_csv(V3_MAP)
    xyz_atoms = read_xyz(V3_XYZ)

    source = {
        row["node_id"]: {
            "element": row["element"],
            "xyz": np.array(
                [
                    float(row["x_nm"]) * 10,
                    float(row["y_nm"]) * 10,
                    float(row["z_nm"]) * 10,
                ],
                dtype=float,
            ),
        }
        for row in source_rows
    }

    v3 = {}

    for atom, row in zip(
        xyz_atoms,
        map_rows,
        strict=True,
    ):
        v3[row["atom_id"]] = {
            "element": atom["element"],
            "xyz": atom["xyz"],
            "fixed": (
                row["v3a2_fixed"].strip().lower()
                == "true"
            ),
            "mobile": (
                row["v3a2_mobile"].strip().lower()
                == "true"
            ),
            "atom_role": row["atom_role"],
            "node_type": row["node_type"],
            "mobility_basis": (
                row["v3a2_mobility_basis"]
            ),
            "artificial_cap": (
                row["artificial_cap"].strip().lower()
                == "true"
            ),
        }

    shared_real = sorted(
        atom_id
        for atom_id in v3
        if atom_id in source
        and atom_id not in DEFECTIVE_CAPS
        and not v3[atom_id]["artificial_cap"]
    )

    subsets = {
        "ALL_SHARED_REAL": shared_real,
        "FIXED_SHARED_REAL": [
            atom_id for atom_id in shared_real
            if v3[atom_id]["fixed"]
        ],
        "MOBILE_SHARED_REAL": [
            atom_id for atom_id in shared_real
            if v3[atom_id]["mobile"]
        ],
        "FIXED_PARENT_HBN": [
            atom_id for atom_id in shared_real
            if (
                v3[atom_id]["fixed"]
                and v3[atom_id]["node_type"]
                == "PARENT_HBN"
            )
        ],
        "FIXED_RETAINED_VALIDATED_CORE": [
            atom_id for atom_id in shared_real
            if (
                v3[atom_id]["fixed"]
                and v3[atom_id]["mobility_basis"]
                == "RETAINED_VALIDATED_CORE"
            )
        ],
    }

    records = []
    detailed = {}

    for label, atom_ids in subsets.items():
        if len(atom_ids) < 3:
            records.append({
                "subset": label,
                "atom_count": len(atom_ids),
                "rmsd_A": "",
                "maximum_residual_A": "",
                "rotation_determinant": "",
                "status": "INSUFFICIENT_ATOMS",
            })
            continue

        source_matrix = np.vstack([
            source[atom_id]["xyz"]
            for atom_id in atom_ids
        ])

        target_matrix = np.vstack([
            v3[atom_id]["xyz"]
            for atom_id in atom_ids
        ])

        result = kabsch(
            source_matrix,
            target_matrix,
        )

        records.append({
            "subset": label,
            "atom_count": len(atom_ids),
            "rmsd_A": result["rmsd_A"],
            "maximum_residual_A": (
                result["maximum_residual_A"]
            ),
            "rotation_determinant": (
                result["determinant"]
            ),
            "status": "EVALUATED",
        })

        detailed[label] = {
            "atom_ids": atom_ids,
            "rmsd_A": result["rmsd_A"],
            "maximum_residual_A": (
                result["maximum_residual_A"]
            ),
            "rotation_matrix": (
                result["rotation"].tolist()
            ),
            "translation_vector_A": (
                result["translation"].tolist()
            ),
            "rotation_determinant": (
                result["determinant"]
            ),
            "reflection_corrected": bool(
                result["reflection_corrected"]
            ),
        }

    with CSV_PATH.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(records[0]),
        )
        writer.writeheader()
        writer.writerows(records)

    viable = [
        row for row in records
        if row["status"] == "EVALUATED"
    ]

    best = min(
        viable,
        key=lambda row: float(row["rmsd_A"]),
    )

    report = {
        "decision": (
            "QM_F06_UPPER_V4_ALIGNMENT_SUBSETS_AUDITED_"
            "ANCHOR_SELECTION_PENDING"
        ),
        "subsets": detailed,
        "best_subset": best,
        "authorization": {
            "alignment_anchor_selected": False,
            "v4_geometry_construction_authorized": False,
            "orca_execution_authorized": False,
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    REPORT.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 90)
    print("QM_F06 UPPER V4 ALIGNMENT-SUBSET AUDIT")
    print("=" * 90)

    for row in records:
        print(
            f"{row['subset']:34s} "
            f"n={row['atom_count']:2d} "
            f"RMSD={str(row['rmsd_A']):18s} "
            f"MAX={str(row['maximum_residual_A']):18s} "
            f"{row['status']}"
        )

    print()
    print("Best subset:", best["subset"])
    print("Best RMSD A:", best["rmsd_A"])
    print(
        "Best maximum residual A:",
        best["maximum_residual_A"],
    )
    print()
    print(
        "Decision:",
        report["decision"],
    )
    print("Report:", REPORT)
    print("CSV:", CSV_PATH)
    print("V4 construction authorized: False")


if __name__ == "__main__":
    main()
