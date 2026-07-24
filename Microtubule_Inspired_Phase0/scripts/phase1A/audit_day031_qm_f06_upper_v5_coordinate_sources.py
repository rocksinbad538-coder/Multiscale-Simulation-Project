#!/usr/bin/env python3
"""
Audit candidate coordinate sources for QM_F06 UPPER V5.

For each candidate source:
- verifies presence of BR4:UPPER:14:1, P:1637 and S:1737;
- finds atoms shared with the V4 fixed core;
- performs a Kabsch rigid alignment using only fixed shared real atoms;
- rejects reflections;
- measures anchor RMSD and maximum residual;
- transforms the three selected V5 atoms into the V4 coordinate frame;
- checks the newly internalized canonical bonds:
    S:1738 -- BR4:UPPER:14:1
    S:1738 -- P:1637
    P:1637 -- S:1737
- ranks viable coordinate sources;
- authorizes transformed-coordinate extraction only when all gates pass.

It does not construct V5 or authorize ORCA.
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

MODEL_SELECTION = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5_model_selection/"
    "QM_F06_UPPER_V5_MODEL_SELECTION.json"
)

V4_XYZ = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_construction/"
    "QM_F06_UPPER_V4_start.xyz"
)

V4_CONSTRAINT_MAP = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_constraint_design/"
    "QM_F06_UPPER_V4_constraint_map.csv"
)

CANDIDATE_SOURCES = [
    ROOT / (
        "runs/phase1A/day024_chemical_end_rim_design/"
        "08_r2_partial_attachment_annulus_static_coordinate_embedding/"
        "r2_partial_attachment_static_coordinates.csv"
    ),
    ROOT / (
        "runs/phase1A/day024_chemical_end_rim_design/"
        "19_r2_selected_four_atom_heavy_coordinate_embedding/"
        "r2_selected_four_atom_heavy_coordinates.csv"
    ),
    ROOT / (
        "runs/phase1A/day024_chemical_end_rim_design/"
        "13_r2_trimer_bridge_conformer_and_h_refinement/"
        "r2_trimer_bridge_refined_coordinates.csv"
    ),
    ROOT / (
        "runs/phase1A/day024_chemical_end_rim_design/"
        "23_r2_four_atom_hydrogen_coordinate_embedding/"
        "r2_selected_four_atom_full_coordinates.csv"
    ),
    ROOT / (
        "runs/phase1A/day024_chemical_end_rim_design/"
        "12_r2_alternating_bn_trimer_bridge_static_coordinate_embedding/"
        "r2_trimer_bridge_static_coordinates.csv"
    ),
    ROOT / (
        "runs/phase1A/day024_chemical_end_rim_design/"
        "18_r2_selected_four_atom_exact_conformer_replay/"
        "r2_selected_four_atom_exact_bridge_coordinates.csv"
    ),
    ROOT / (
        "runs/phase1A/day024_chemical_end_rim_design/"
        "28_r2_inner_h_reflected_direction_refinement/"
        "r2_selected_four_atom_refined_full_coordinates.csv"
    ),
]

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5_coordinate_sources"
)

REPORT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5_COORDINATE_SOURCE_AUDIT.json"
)

SOURCE_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5_coordinate_source_ranking.csv"
)

TRANSFORMED_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5_selected_transformed_coordinates.csv"
)

ANCHOR_RESIDUAL_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5_selected_anchor_residuals.csv"
)

TARGET_ATOMS = {
    "BR4:UPPER:14:1",
    "P:1637",
    "S:1737",
}

REQUIRED_BONDS = [
    ("S:1738", "BR4:UPPER:14:1"),
    ("S:1738", "P:1637"),
    ("P:1637", "S:1737"),
]

MINIMUM_ANCHOR_ATOMS = 6
ANCHOR_RMSD_MAX_A = 5.0e-4
ANCHOR_MAX_RESIDUAL_A = 1.0e-3
BN_MIN_A = 1.20
BN_MAX_A = 1.85


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


def read_csv(path: Path):
    require_file(path)

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames or [], list(reader)


def read_xyz(path: Path):
    require_file(path)

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
            "xyz_A": np.array(
                [float(value) for value in fields[1:4]],
                dtype=float,
            ),
        })

    if len(atoms) != count:
        raise RuntimeError(
            f"Incomplete XYZ: {path}"
        )

    return atoms


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {
        "true",
        "1",
        "yes",
        "pass",
    }


def resolve_id_column(headers):
    for candidate in (
        "node_id",
        "atom_id",
        "bridge_node",
        "id",
    ):
        if candidate in headers:
            return candidate

    raise RuntimeError(
        f"No atom-ID column found: {headers}"
    )


def resolve_element_column(headers):
    for candidate in (
        "element",
        "atom_element",
    ):
        if candidate in headers:
            return candidate

    raise RuntimeError(
        f"No element column found: {headers}"
    )


def resolve_coordinate_columns(headers):
    candidates = [
        ("x_nm", "y_nm", "z_nm", 10.0),
        ("x_angstrom", "y_angstrom", "z_angstrom", 1.0),
        ("x_A", "y_A", "z_A", 1.0),
        ("x", "y", "z", 1.0),
    ]

    for x_key, y_key, z_key, scale in candidates:
        if all(
            key in headers
            for key in (x_key, y_key, z_key)
        ):
            return x_key, y_key, z_key, scale

    raise RuntimeError(
        f"No coordinate columns found: {headers}"
    )


def kabsch(source, target):
    source_centroid = source.mean(axis=0)
    target_centroid = target.mean(axis=0)

    source_centered = source - source_centroid
    target_centered = target - target_centroid

    covariance = source_centered.T @ target_centered

    u_matrix, _, vt_matrix = np.linalg.svd(
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
        - source_centroid @ rotation.T
    )

    transformed = (
        source @ rotation.T
        + translation
    )

    return (
        rotation,
        translation,
        transformed,
        reflection_corrected,
    )


def distance(first, second):
    return float(np.linalg.norm(first - second))


def main() -> None:
    for path in (
        MODEL_SELECTION,
        V4_XYZ,
        V4_CONSTRAINT_MAP,
    ):
        require_file(path)

    selection = json.loads(
        MODEL_SELECTION.read_text(
            encoding="utf-8"
        )
    )

    if not selection["authorization"][
        "coordinate_source_validation_authorized"
    ]:
        raise RuntimeError(
            "Coordinate-source validation is not authorized."
        )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    v4_atoms = read_xyz(V4_XYZ)
    _, constraint_rows = read_csv(
        V4_CONSTRAINT_MAP
    )

    if len(v4_atoms) != len(constraint_rows):
        raise RuntimeError(
            "V4 XYZ / constraint-map atom-count mismatch."
        )

    v4_by_id = {}
    fixed_ids = set()

    for atom, row in zip(
        v4_atoms,
        constraint_rows,
        strict=True,
    ):
        index = int(row["index_0based"])

        if index != atom["index"]:
            raise RuntimeError(
                f"V4 index mismatch at {index}"
            )

        atom_id = row["atom_id"]

        v4_by_id[atom_id] = {
            "element": atom["element"],
            "xyz_A": atom["xyz_A"],
        }

        if parse_bool(row["v4_fixed"]):
            fixed_ids.add(atom_id)

    source_records = []
    source_details = {}

    for source_path in CANDIDATE_SOURCES:
        base_record = {
            "source": str(
                source_path.relative_to(ROOT)
            ),
            "file_exists": source_path.is_file(),
            "target_count": 0,
            "all_targets_present": False,
            "shared_fixed_count": 0,
            "anchor_rmsd_A": None,
            "anchor_max_residual_A": None,
            "rotation_determinant": None,
            "reflection_corrected": None,
            "required_bonds_pass": False,
            "maximum_required_bond_A": None,
            "minimum_required_bond_A": None,
            "overall_pass": False,
            "failure_reason": "",
        }

        if not source_path.is_file():
            base_record["failure_reason"] = (
                "SOURCE_MISSING"
            )
            source_records.append(base_record)
            continue

        try:
            headers, rows = read_csv(source_path)

            id_key = resolve_id_column(headers)
            element_key = resolve_element_column(
                headers
            )

            (
                x_key,
                y_key,
                z_key,
                scale,
            ) = resolve_coordinate_columns(headers)

            source_by_id = {}

            for row in rows:
                atom_id = row[id_key].strip()

                if not atom_id:
                    continue

                values = [
                    row[x_key].strip(),
                    row[y_key].strip(),
                    row[z_key].strip(),
                ]

                if any(not value for value in values):
                    continue

                source_by_id[atom_id] = {
                    "element": row[element_key],
                    "xyz_A": np.array(
                        [
                            float(values[0]) * scale,
                            float(values[1]) * scale,
                            float(values[2]) * scale,
                        ],
                        dtype=float,
                    ),
                }

            present_targets = (
                TARGET_ATOMS
                & set(source_by_id)
            )

            base_record["target_count"] = len(
                present_targets
            )

            base_record["all_targets_present"] = (
                present_targets == TARGET_ATOMS
            )

            if present_targets != TARGET_ATOMS:
                base_record["failure_reason"] = (
                    "MISSING_TARGETS:"
                    + "|".join(
                        sorted(
                            TARGET_ATOMS
                            - present_targets
                        )
                    )
                )
                source_records.append(base_record)
                continue

            shared_fixed = sorted(
                fixed_ids
                & set(source_by_id)
            )

            base_record["shared_fixed_count"] = len(
                shared_fixed
            )

            if len(shared_fixed) < MINIMUM_ANCHOR_ATOMS:
                base_record["failure_reason"] = (
                    "INSUFFICIENT_FIXED_ANCHORS"
                )
                source_records.append(base_record)
                continue

            source_anchor = np.array([
                source_by_id[atom_id]["xyz_A"]
                for atom_id in shared_fixed
            ])

            target_anchor = np.array([
                v4_by_id[atom_id]["xyz_A"]
                for atom_id in shared_fixed
            ])

            (
                rotation,
                translation,
                transformed_anchor,
                reflection_corrected,
            ) = kabsch(
                source_anchor,
                target_anchor,
            )

            residuals = np.linalg.norm(
                transformed_anchor
                - target_anchor,
                axis=1,
            )

            rmsd = float(math.sqrt(
                float(np.mean(residuals ** 2))
            ))

            maximum_residual = float(
                np.max(residuals)
            )

            determinant = float(
                np.linalg.det(rotation)
            )

            base_record["anchor_rmsd_A"] = rmsd
            base_record[
                "anchor_max_residual_A"
            ] = maximum_residual
            base_record[
                "rotation_determinant"
            ] = determinant
            base_record[
                "reflection_corrected"
            ] = reflection_corrected

            transformed_coordinates = {}

            for atom_id, record in source_by_id.items():
                transformed_coordinates[atom_id] = (
                    record["xyz_A"] @ rotation.T
                    + translation
                )

            combined = {
                atom_id: record["xyz_A"]
                for atom_id, record in v4_by_id.items()
            }

            for atom_id in TARGET_ATOMS:
                combined[atom_id] = (
                    transformed_coordinates[atom_id]
                )

            bond_values = []

            for first, second in REQUIRED_BONDS:
                if first not in combined or second not in combined:
                    raise RuntimeError(
                        f"Missing required bond atom: "
                        f"{first} -- {second}"
                    )

                bond_values.append({
                    "first": first,
                    "second": second,
                    "distance_A": distance(
                        combined[first],
                        combined[second],
                    ),
                })

            bond_pass = all(
                BN_MIN_A
                <= record["distance_A"]
                <= BN_MAX_A
                for record in bond_values
            )

            base_record[
                "required_bonds_pass"
            ] = bond_pass

            base_record[
                "maximum_required_bond_A"
            ] = max(
                record["distance_A"]
                for record in bond_values
            )

            base_record[
                "minimum_required_bond_A"
            ] = min(
                record["distance_A"]
                for record in bond_values
            )

            anchor_pass = (
                rmsd <= ANCHOR_RMSD_MAX_A
                and maximum_residual
                <= ANCHOR_MAX_RESIDUAL_A
                and abs(determinant - 1.0)
                <= 1.0e-8
            )

            overall_pass = (
                anchor_pass
                and bond_pass
            )

            base_record["overall_pass"] = (
                overall_pass
            )

            if not anchor_pass:
                base_record["failure_reason"] = (
                    "ANCHOR_ALIGNMENT_FAIL"
                )
            elif not bond_pass:
                base_record["failure_reason"] = (
                    "RESTORED_BOND_FAIL"
                )

            source_details[
                str(source_path.relative_to(ROOT))
            ] = {
                "shared_fixed_ids": shared_fixed,
                "rotation_matrix": (
                    rotation.tolist()
                ),
                "translation_A": (
                    translation.tolist()
                ),
                "anchor_residuals": [
                    {
                        "atom_id": atom_id,
                        "residual_A": float(value),
                    }
                    for atom_id, value in zip(
                        shared_fixed,
                        residuals,
                        strict=True,
                    )
                ],
                "transformed_targets": {
                    atom_id: {
                        "element": source_by_id[
                            atom_id
                        ]["element"],
                        "xyz_A": transformed_coordinates[
                            atom_id
                        ].tolist(),
                    }
                    for atom_id in sorted(
                        TARGET_ATOMS
                    )
                },
                "required_bonds": bond_values,
            }

        except Exception as exc:
            base_record["failure_reason"] = (
                f"EXCEPTION:{type(exc).__name__}:{exc}"
            )

        source_records.append(base_record)

    source_records.sort(
        key=lambda row: (
            not row["overall_pass"],
            (
                row["anchor_rmsd_A"]
                if row["anchor_rmsd_A"] is not None
                else float("inf")
            ),
            (
                row["anchor_max_residual_A"]
                if row["anchor_max_residual_A"]
                is not None
                else float("inf")
            ),
            row["source"],
        )
    )

    for rank, row in enumerate(
        source_records,
        start=1,
    ):
        row["rank"] = rank

    passing = [
        row for row in source_records
        if row["overall_pass"]
    ]

    selected = passing[0] if passing else None

    if selected is None:
        decision = (
            "QM_F06_UPPER_V5_COORDINATE_SOURCE_FAIL_"
            "CONSTRUCTION_BLOCKED"
        )
        coordinate_authorized = False
    else:
        decision = (
            "QM_F06_UPPER_V5_COORDINATE_SOURCE_PASS_"
            "TRANSFORMED_COORDINATES_AUTHORIZED"
        )
        coordinate_authorized = True

    with SOURCE_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        fieldnames = [
            "rank",
            "source",
            "file_exists",
            "target_count",
            "all_targets_present",
            "shared_fixed_count",
            "anchor_rmsd_A",
            "anchor_max_residual_A",
            "rotation_determinant",
            "reflection_corrected",
            "required_bonds_pass",
            "maximum_required_bond_A",
            "minimum_required_bond_A",
            "overall_pass",
            "failure_reason",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(source_records)

    if selected is not None:
        selected_source = selected["source"]
        details = source_details[selected_source]

        with TRANSFORMED_CSV.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "atom_id",
                    "element",
                    "x_A",
                    "y_A",
                    "z_A",
                    "coordinate_source",
                ],
            )
            writer.writeheader()

            for atom_id in sorted(TARGET_ATOMS):
                record = details[
                    "transformed_targets"
                ][atom_id]

                writer.writerow({
                    "atom_id": atom_id,
                    "element": record["element"],
                    "x_A": record["xyz_A"][0],
                    "y_A": record["xyz_A"][1],
                    "z_A": record["xyz_A"][2],
                    "coordinate_source": selected_source,
                })

        with ANCHOR_RESIDUAL_CSV.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "atom_id",
                    "residual_A",
                ],
            )
            writer.writeheader()
            writer.writerows(
                details["anchor_residuals"]
            )

    files = {
        "source_ranking": SOURCE_CSV,
    }

    if selected is not None:
        files[
            "selected_transformed_coordinates"
        ] = TRANSFORMED_CSV
        files[
            "selected_anchor_residuals"
        ] = ANCHOR_RESIDUAL_CSV

    report = {
        "decision": decision,
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "target_atoms": sorted(TARGET_ATOMS),
        "required_bonds": [
            list(pair)
            for pair in REQUIRED_BONDS
        ],
        "thresholds": {
            "minimum_anchor_atoms": (
                MINIMUM_ANCHOR_ATOMS
            ),
            "anchor_rmsd_max_A": (
                ANCHOR_RMSD_MAX_A
            ),
            "anchor_max_residual_A": (
                ANCHOR_MAX_RESIDUAL_A
            ),
            "BN_min_A": BN_MIN_A,
            "BN_max_A": BN_MAX_A,
        },
        "candidate_sources": source_records,
        "selected_source": (
            selected["source"]
            if selected is not None
            else None
        ),
        "selected_source_details": (
            source_details[selected["source"]]
            if selected is not None
            else None
        ),
        "files": {
            key: str(path.relative_to(ROOT))
            for key, path in files.items()
        },
        "files_sha256": {
            key: sha256(path)
            for key, path in files.items()
        },
        "authorization": {
            "coordinate_source_selected": (
                coordinate_authorized
            ),
            "transformed_coordinates_authorized": (
                coordinate_authorized
            ),
            "v5_geometry_construction_authorized": (
                coordinate_authorized
            ),
            "orca_input_preparation_authorized": False,
            "orca_execution_authorized": False,
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    REPORT_PATH.write_text(
        json.dumps(
            report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 100)
    print("QM_F06 UPPER V5 COORDINATE-SOURCE AUDIT")
    print("=" * 100)

    for row in source_records:
        print(
            f"rank={row['rank']:2d} "
            f"targets={row['target_count']}/3 "
            f"anchors={row['shared_fixed_count']:2d} "
            f"RMSD={str(row['anchor_rmsd_A']):18s} "
            f"MAX={str(row['anchor_max_residual_A']):18s} "
            f"bonds={row['required_bonds_pass']} "
            f"pass={row['overall_pass']} "
            f"{row['source']}"
        )

        if row["failure_reason"]:
            print(
                "   reason:",
                row["failure_reason"],
            )

    print()
    print(
        "Selected source:",
        (
            selected["source"]
            if selected is not None
            else None
        ),
    )
    print("Decision:", decision)
    print("Report:", REPORT_PATH)
    print("Ranking:", SOURCE_CSV)

    if selected is not None:
        print(
            "Transformed coordinates:",
            TRANSFORMED_CSV,
        )
        print(
            "Anchor residuals:",
            ANCHOR_RESIDUAL_CSV,
        )

    print()
    print(
        "V5 construction authorized:",
        coordinate_authorized,
    )
    print("ORCA execution authorized: False")


if __name__ == "__main__":
    main()
