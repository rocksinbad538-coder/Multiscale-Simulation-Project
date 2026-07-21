#!/usr/bin/env python3
"""
Construct the selected QM_F06 UPPER V4 starting geometry.

Inputs:
- validated V3-A2 starting geometry and atom map;
- canonical refined R2 coordinates;
- accepted fixed-core rigid transformation.

Construction:
- removes only HCAP:UPPER:01 and HCAP:UPPER:04;
- retains the other 28 V3-A2 atoms in their original order;
- restores 12 canonical R2 atoms;
- introduces seven artificial caps along canonical cut-bond vectors;
- writes deterministic geometry and provenance artifacts.

This script does not assign V4 optimization constraints and does not
authorize ORCA execution.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

V3_ROOT = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V3/orca_v3a2_workflow"
)

V3_XYZ = V3_ROOT / "v3a2_start.xyz"
V3_MAP = V3_ROOT / "v3a2_atom_role_constraint_map.csv"

SOURCE_COORDINATES = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "28_r2_inner_h_reflected_direction_refinement/"
    "r2_selected_four_atom_refined_full_coordinates.csv"
)

ANCHOR_REPORT = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_fixed_anchor/"
    "QM_F06_UPPER_V4_FIXED_ANCHOR_REPORT.json"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_construction"
)

OUTPUT_XYZ = OUTPUT_DIR / "QM_F06_UPPER_V4_start.xyz"

OUTPUT_MAP = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_atom_role_provenance_map.csv"
)

OUTPUT_CAPS = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_artificial_caps.csv"
)

OUTPUT_REPORT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_CONSTRUCTION_REPORT.json"
)

REMOVED_V3_CAPS = {
    # Defective geminal caps formerly attached to P:1641.
    "HCAP:UPPER:01",
    "HCAP:UPPER:04",

    # Superseded because the real canonical neighbor P:1642
    # is explicitly restored in V4.
    "HCAP:UPPER:02",
}

RESTORED_ATOM_ORDER = [
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
]

# cap_id, center, omitted canonical neighbor, bond length in Å
CAP_SPECS = [
    (
        "HCAPV4:UPPER:P1580:P1521",
        "P:1580",
        "P:1521",
        1.19,
    ),
    (
        "HCAPV4:UPPER:P1582:P1525",
        "P:1582",
        "P:1525",
        1.19,
    ),
    (
        "HCAPV4:UPPER:P1638:P1579",
        "P:1638",
        "P:1579",
        1.19,
    ),
    (
        "HCAPV4:UPPER:P1642:P1585",
        "P:1642",
        "P:1585",
        1.19,
    ),
    (
        "HCAPV4:UPPER:S1738:BR4_14_1",
        "S:1738",
        "BR4:UPPER:14:1",
        1.19,
    ),
    (
        "HCAPV4:UPPER:S1738:P1637",
        "S:1738",
        "P:1637",
        1.19,
    ),
    (
        "HCAPV4:UPPER:P1523:P1522",
        "P:1523",
        "P:1522",
        1.01,
    ),
]

EXPECTED_COMPOSITION = Counter({
    "B": 15,
    "N": 11,
    "H": 20,
})

EXPECTED_ATOM_COUNT = 46


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

    declared = int(lines[0].strip())

    atoms = []

    for index, line in enumerate(lines[2:2 + declared]):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Malformed XYZ line {index}: {line}"
            )

        atoms.append({
            "index": index,
            "element": fields[0],
            "xyz_A": tuple(
                float(value)
                for value in fields[1:4]
            ),
        })

    if len(atoms) != declared:
        raise RuntimeError(
            f"Incomplete XYZ: expected {declared}, "
            f"found {len(atoms)}"
        )

    return atoms


def matvec(matrix, vector):
    return tuple(
        sum(
            float(matrix[row][column])
            * float(vector[column])
            for column in range(3)
        )
        for row in range(3)
    )


def add(first, second):
    return tuple(
        a + b
        for a, b in zip(first, second)
    )


def subtract(first, second):
    return tuple(
        a - b
        for a, b in zip(first, second)
    )


def norm(vector):
    return math.sqrt(sum(value * value for value in vector))


def scale(vector, factor):
    return tuple(value * factor for value in vector)


def distance(first, second):
    return norm(subtract(first, second))


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    v3_atoms = read_xyz(V3_XYZ)
    v3_map_rows = read_csv(V3_MAP)
    source_rows = read_csv(SOURCE_COORDINATES)

    require_file(ANCHOR_REPORT)

    anchor_report = json.loads(
        ANCHOR_REPORT.read_text(encoding="utf-8")
    )

    if not anchor_report["overall_pass"]:
        raise RuntimeError(
            "Fixed-anchor report did not pass."
        )

    if not anchor_report["authorization"][
        "v4_geometry_construction_authorized"
    ]:
        raise RuntimeError(
            "V4 construction is not authorized by anchor report."
        )

    if len(v3_atoms) != len(v3_map_rows):
        raise RuntimeError(
            "V3 XYZ and map row counts differ."
        )

    rotation = anchor_report["rotation_matrix"]
    translation = tuple(
        float(value)
        for value in anchor_report[
            "translation_vector_A"
        ]
    )

    source = {}

    for row in source_rows:
        atom_id = row["node_id"]

        if atom_id in source:
            raise RuntimeError(
                f"Duplicate canonical atom ID: {atom_id}"
            )

        source_xyz_A = (
            float(row["x_nm"]) * 10.0,
            float(row["y_nm"]) * 10.0,
            float(row["z_nm"]) * 10.0,
        )

        transformed_xyz_A = add(
            matvec(rotation, source_xyz_A),
            translation,
        )

        source[atom_id] = {
            "atom_id": atom_id,
            "element": row["element"],
            "node_type": row["node_type"],
            "coordinate_source": row[
                "coordinate_source"
            ],
            "xyz_A": transformed_xyz_A,
        }

    required_source_ids = set(RESTORED_ATOM_ORDER)

    for _, center, outside, _ in CAP_SPECS:
        required_source_ids.add(center)
        required_source_ids.add(outside)

    missing_source_ids = (
        required_source_ids - set(source)
    )

    if missing_source_ids:
        raise RuntimeError(
            "Canonical coordinate source is missing: "
            f"{sorted(missing_source_ids)}"
        )

    atoms = []
    provenance = []

    # Retain V3 atoms in deterministic original order.
    for atom, row in zip(
        v3_atoms,
        v3_map_rows,
        strict=True,
    ):
        atom_id = row["atom_id"]

        if atom_id in REMOVED_V3_CAPS:
            continue

        atoms.append({
            "atom_id": atom_id,
            "element": atom["element"],
            "xyz_A": atom["xyz_A"],
        })

        provenance.append({
            "atom_id": atom_id,
            "element": atom["element"],
            "atom_role": row["atom_role"],
            "node_type": row["node_type"],
            "coordinate_source": (
                "RETAINED_VALIDATED_V3A2_START"
            ),
            "source_center": "",
            "source_outside_atom": "",
            "source_cut_edge": "",
            "artificial_cap": row["artificial_cap"],
            "v4_constraint_status": "UNASSIGNED",
            "v4_mobility_basis": "PENDING_PRE_QM_DESIGN",
        })

    retained_ids = {
        atom["atom_id"]
        for atom in atoms
    }

    if retained_ids & set(RESTORED_ATOM_ORDER):
        overlap = retained_ids & set(RESTORED_ATOM_ORDER)

        raise RuntimeError(
            "Restored atoms already exist in retained V3: "
            f"{sorted(overlap)}"
        )

    # Append restored real atoms in explicit deterministic order.
    for atom_id in RESTORED_ATOM_ORDER:
        record = source[atom_id]

        atoms.append({
            "atom_id": atom_id,
            "element": record["element"],
            "xyz_A": record["xyz_A"],
        })

        provenance.append({
            "atom_id": atom_id,
            "element": record["element"],
            "atom_role": (
                "RESTORED_CANONICAL_R2_ATOM"
            ),
            "node_type": record["node_type"],
            "coordinate_source": (
                "FIXED_ANCHOR_TRANSFORMED_"
                + record["coordinate_source"]
            ),
            "source_center": "",
            "source_outside_atom": "",
            "source_cut_edge": "",
            "artificial_cap": "False",
            "v4_constraint_status": "UNASSIGNED",
            "v4_mobility_basis": "PENDING_PRE_QM_DESIGN",
        })

    cap_records = []

    # Append artificial caps along omitted heavy-bond directions.
    for cap_id, center_id, outside_id, bond_length_A in CAP_SPECS:
        center = source[center_id]["xyz_A"]
        outside = source[outside_id]["xyz_A"]

        cut_vector = subtract(outside, center)
        cut_distance_A = norm(cut_vector)

        if cut_distance_A <= 0.0:
            raise RuntimeError(
                f"Zero cut vector: {center_id} -- {outside_id}"
            )

        unit_vector = scale(
            cut_vector,
            1.0 / cut_distance_A,
        )

        cap_xyz_A = add(
            center,
            scale(unit_vector, bond_length_A),
        )

        realized_distance_A = distance(
            center,
            cap_xyz_A,
        )

        atoms.append({
            "atom_id": cap_id,
            "element": "H",
            "xyz_A": cap_xyz_A,
        })

        provenance.append({
            "atom_id": cap_id,
            "element": "H",
            "atom_role": "ARTIFICIAL_BOUNDARY_CAP",
            "node_type": "QM_BOUNDARY_CAP_H",
            "coordinate_source": (
                "TRANSFORMED_CANONICAL_CUT_BOND_VECTOR"
            ),
            "source_center": center_id,
            "source_outside_atom": outside_id,
            "source_cut_edge": (
                f"{center_id}--{outside_id}"
            ),
            "artificial_cap": "True",
            "v4_constraint_status": "UNASSIGNED",
            "v4_mobility_basis": "PENDING_PRE_QM_DESIGN",
        })

        cap_records.append({
            "cap_id": cap_id,
            "center_atom": center_id,
            "center_element": source[
                center_id
            ]["element"],
            "omitted_outside_atom": outside_id,
            "omitted_outside_element": source[
                outside_id
            ]["element"],
            "canonical_cut_distance_A": (
                cut_distance_A
            ),
            "target_cap_bond_length_A": (
                bond_length_A
            ),
            "realized_cap_bond_length_A": (
                realized_distance_A
            ),
            "cap_x_A": cap_xyz_A[0],
            "cap_y_A": cap_xyz_A[1],
            "cap_z_A": cap_xyz_A[2],
        })

    atom_ids = [
        atom["atom_id"]
        for atom in atoms
    ]

    if len(atom_ids) != len(set(atom_ids)):
        duplicates = sorted({
            atom_id
            for atom_id in atom_ids
            if atom_ids.count(atom_id) > 1
        })

        raise RuntimeError(
            f"Duplicate final atom IDs: {duplicates}"
        )

    composition = Counter(
        atom["element"]
        for atom in atoms
    )

    atom_count_gate = (
        len(atoms) == EXPECTED_ATOM_COUNT
    )

    composition_gate = (
        composition == EXPECTED_COMPOSITION
    )

    removed_cap_gate = (
        not (REMOVED_V3_CAPS & set(atom_ids))
    )

    restored_atom_gate = (
        set(RESTORED_ATOM_ORDER)
        .issubset(atom_ids)
    )

    new_cap_gate = (
        len(cap_records) == 7
        and all(
            abs(
                row["realized_cap_bond_length_A"]
                - row["target_cap_bond_length_A"]
            )
            <= 1.0e-10
            for row in cap_records
        )
    )

    basic_construction_pass = all((
        atom_count_gate,
        composition_gate,
        removed_cap_gate,
        restored_atom_gate,
        new_cap_gate,
    ))

    with OUTPUT_XYZ.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(f"{len(atoms)}\n")
        handle.write(
            "QM_F06 UPPER V4 selected P:1523 expansion; "
            "pre-QM geometry; ORCA not authorized\n"
        )

        for atom in atoms:
            x_value, y_value, z_value = atom["xyz_A"]

            handle.write(
                f"{atom['element']:2s} "
                f"{x_value: .12f} "
                f"{y_value: .12f} "
                f"{z_value: .12f}\n"
            )

    map_fields = [
        "index_0based",
        "atom_id",
        "element",
        "atom_role",
        "node_type",
        "coordinate_source",
        "source_center",
        "source_outside_atom",
        "source_cut_edge",
        "artificial_cap",
        "v4_constraint_status",
        "v4_mobility_basis",
    ]

    with OUTPUT_MAP.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=map_fields,
        )
        writer.writeheader()

        for index, row in enumerate(provenance):
            writer.writerow({
                "index_0based": index,
                **row,
            })

    with OUTPUT_CAPS.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(cap_records[0]),
        )
        writer.writeheader()
        writer.writerows(cap_records)

    report = {
        "decision": (
            "QM_F06_UPPER_V4_START_GEOMETRY_CONSTRUCTED_"
            "PRE_QM_AUDIT_REQUIRED"
            if basic_construction_pass
            else
            "QM_F06_UPPER_V4_CONSTRUCTION_FAIL"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "selected_strategy": (
            "P1523_SELECTIVE_EXPANSION"
        ),
        "retained_v3_atom_count": 27,
        "removed_v3_caps": sorted(
            REMOVED_V3_CAPS
        ),
        "restored_real_atom_count": len(
            RESTORED_ATOM_ORDER
        ),
        "restored_real_atoms": (
            RESTORED_ATOM_ORDER
        ),
        "new_artificial_cap_count": len(
            cap_records
        ),
        "final_atom_count": len(atoms),
        "composition": dict(
            sorted(composition.items())
        ),
        "expected_composition": dict(
            sorted(EXPECTED_COMPOSITION.items())
        ),
        "gates": {
            "atom_count": atom_count_gate,
            "composition": composition_gate,
            "defective_caps_removed": removed_cap_gate,
            "restored_atoms_present": restored_atom_gate,
            "new_caps_constructed": new_cap_gate,
        },
        "basic_construction_pass": (
            basic_construction_pass
        ),
        "files": {
            "xyz": str(
                OUTPUT_XYZ.relative_to(ROOT)
            ),
            "atom_map": str(
                OUTPUT_MAP.relative_to(ROOT)
            ),
            "cap_inventory": str(
                OUTPUT_CAPS.relative_to(ROOT)
            ),
        },
        "authorization": {
            "v4_start_geometry_constructed": (
                basic_construction_pass
            ),
            "pre_qm_structural_audit_authorized": (
                basic_construction_pass
            ),
            "orca_input_preparation_authorized": False,
            "orca_execution_authorized": False,
            "geometric_reference_accepted": False,
            "electronic_reference_accepted": False,
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    OUTPUT_REPORT.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 78)
    print("QM_F06 UPPER V4 CONSTRUCTION")
    print("=" * 78)
    print("Retained V3 atoms:", 27)
    print(
        "Restored real atoms:",
        len(RESTORED_ATOM_ORDER),
    )
    print("New artificial caps:", len(cap_records))
    print("Final atom count:", len(atoms))
    print("Composition:", dict(composition))
    print()
    print("Atom-count gate:", atom_count_gate)
    print("Composition gate:", composition_gate)
    print(
        "Defective-cap removal gate:",
        removed_cap_gate,
    )
    print(
        "Restored-atom gate:",
        restored_atom_gate,
    )
    print("New-cap gate:", new_cap_gate)
    print()
    print("Decision:", report["decision"])
    print("XYZ:", OUTPUT_XYZ)
    print("Map:", OUTPUT_MAP)
    print("Caps:", OUTPUT_CAPS)
    print("Report:", OUTPUT_REPORT)
    print(
        "Pre-QM audit authorized:",
        basic_construction_pass,
    )
    print("ORCA execution authorized: False")


if __name__ == "__main__":
    main()
