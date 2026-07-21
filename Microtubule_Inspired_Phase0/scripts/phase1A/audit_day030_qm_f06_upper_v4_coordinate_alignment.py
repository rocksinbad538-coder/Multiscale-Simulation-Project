#!/usr/bin/env python3
"""
Audit coordinate compatibility between the retained QM_F06 UPPER V3
fragment and the canonical refined R2 coordinate source selected for V4.

The audit:
- reads the V3-A2 starting geometry and atom map;
- reads the final refined R2 coordinates in nm;
- compares all shared real atoms;
- validates the restored V4 B-N and B-H distances;
- verifies that no rigid transformation or unit mismatch exists;
- writes reproducible JSON and CSV artifacts;
- does not construct V4 or authorize ORCA.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path


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
    "day030_qm_f06_upper_v4_coordinate_alignment"
)

REPORT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_COORDINATE_ALIGNMENT_AUDIT.json"
)

SHARED_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_shared_atom_coordinate_deltas.csv"
)

BONDS_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_restored_bond_lengths.csv"
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

SHARED_ATOM_TOLERANCE_A = 1.0e-6
BN_MIN_A = 1.35
BN_MAX_A = 1.60
BH_MIN_A = 1.10
BH_MAX_A = 1.30


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty file: {path}"
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def read_xyz(path: Path):
    require_file(path)

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    count = int(lines[0].strip())

    if len(lines) < count + 2:
        raise RuntimeError(
            f"Incomplete XYZ: {path}"
        )

    atoms = []

    for index, line in enumerate(
        lines[2:2 + count]
    ):
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

    return atoms


def read_csv(path: Path):
    require_file(path)

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def distance(
    xyz_a: tuple[float, float, float],
    xyz_b: tuple[float, float, float],
) -> float:
    return math.sqrt(sum(
        (a - b) ** 2
        for a, b in zip(xyz_a, xyz_b)
    ))


def main() -> None:
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    source_rows = read_csv(
        SOURCE_COORDINATES
    )
    v3_map_rows = read_csv(V3_MAP)
    v3_atoms = read_xyz(V3_XYZ)

    if len(v3_atoms) != len(v3_map_rows):
        raise RuntimeError(
            "V3 XYZ and map row counts differ."
        )

    source = {}

    for row in source_rows:
        atom_id = row["node_id"]

        if atom_id in source:
            raise RuntimeError(
                f"Duplicate source ID: {atom_id}"
            )

        source[atom_id] = {
            "element": row["element"],
            "node_type": row["node_type"],
            "xyz_A": (
                float(row["x_nm"]) * 10.0,
                float(row["y_nm"]) * 10.0,
                float(row["z_nm"]) * 10.0,
            ),
            "coordinate_source": row[
                "coordinate_source"
            ],
        }

    missing_v4 = (
        MANDATORY_V4_ATOMS - set(source)
    )

    if missing_v4:
        raise RuntimeError(
            "Mandatory V4 atoms missing from source: "
            f"{sorted(missing_v4)}"
        )

    v3_by_id = {}

    for atom, row in zip(
        v3_atoms,
        v3_map_rows,
        strict=True,
    ):
        atom_id = row["atom_id"]

        v3_by_id[atom_id] = {
            "element": atom["element"],
            "xyz_A": atom["xyz_A"],
            "artificial_cap": (
                row["artificial_cap"]
                .strip()
                .lower()
                == "true"
            ),
        }

    retained_v3_ids = (
        set(v3_by_id) - DEFECTIVE_CAPS
    )

    shared_real_ids = sorted(
        atom_id
        for atom_id in retained_v3_ids
        if atom_id in source
        and not v3_by_id[
            atom_id
        ]["artificial_cap"]
    )

    if not shared_real_ids:
        raise RuntimeError(
            "No shared real atoms between V3 and source."
        )

    shared_records = []

    for atom_id in shared_real_ids:
        v3_record = v3_by_id[atom_id]
        source_record = source[atom_id]

        delta = distance(
            v3_record["xyz_A"],
            source_record["xyz_A"],
        )

        element_match = (
            v3_record["element"]
            == source_record["element"]
        )

        shared_records.append({
            "atom_id": atom_id,
            "element_v3": v3_record["element"],
            "element_source": (
                source_record["element"]
            ),
            "element_match": element_match,
            "v3_x_A": v3_record["xyz_A"][0],
            "v3_y_A": v3_record["xyz_A"][1],
            "v3_z_A": v3_record["xyz_A"][2],
            "source_x_A": (
                source_record["xyz_A"][0]
            ),
            "source_y_A": (
                source_record["xyz_A"][1]
            ),
            "source_z_A": (
                source_record["xyz_A"][2]
            ),
            "coordinate_delta_A": delta,
            "within_tolerance": (
                element_match
                and delta
                <= SHARED_ATOM_TOLERANCE_A
            ),
        })

    with SHARED_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        fields = list(shared_records[0])
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
        )
        writer.writeheader()
        writer.writerows(shared_records)

    max_shared_delta = max(
        row["coordinate_delta_A"]
        for row in shared_records
    )

    shared_alignment_pass = all(
        row["within_tolerance"]
        for row in shared_records
    )

    coordinate_pool = {
        atom_id: record["xyz_A"]
        for atom_id, record in source.items()
    }

    # For P:1641 and any retained atoms, use the same
    # source coordinate because alignment must be exact.
    bond_records = []

    for first, second in RESTORED_BN_EDGES:
        if first not in coordinate_pool:
            raise RuntimeError(
                f"Missing edge atom in source: {first}"
            )

        if second not in coordinate_pool:
            raise RuntimeError(
                f"Missing edge atom in source: {second}"
            )

        value = distance(
            coordinate_pool[first],
            coordinate_pool[second],
        )

        bond_records.append({
            "first_atom": first,
            "second_atom": second,
            "bond_class": "B-N",
            "distance_A": value,
            "minimum_A": BN_MIN_A,
            "maximum_A": BN_MAX_A,
            "pass": (
                BN_MIN_A <= value <= BN_MAX_A
            ),
        })

    first, second = REAL_BH_EDGE

    bh_distance = distance(
        coordinate_pool[first],
        coordinate_pool[second],
    )

    bond_records.append({
        "first_atom": first,
        "second_atom": second,
        "bond_class": "REAL_R2_B-H",
        "distance_A": bh_distance,
        "minimum_A": BH_MIN_A,
        "maximum_A": BH_MAX_A,
        "pass": (
            BH_MIN_A <= bh_distance <= BH_MAX_A
        ),
    })

    with BONDS_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        fields = list(bond_records[0])
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
        )
        writer.writeheader()
        writer.writerows(bond_records)

    restored_bond_gate = all(
        row["pass"]
        for row in bond_records
    )

    source_coordinate_labels = sorted({
        source[atom_id]["coordinate_source"]
        for atom_id in MANDATORY_V4_ATOMS
    })

    overall_pass = (
        shared_alignment_pass
        and restored_bond_gate
        and len(MANDATORY_V4_ATOMS) == 12
    )

    decision = (
        "QM_F06_UPPER_V4_COORDINATE_ALIGNMENT_PASS_"
        "CONSTRUCTION_INPUTS_COMPATIBLE"
        if overall_pass
        else
        "QM_F06_UPPER_V4_COORDINATE_ALIGNMENT_FAIL_"
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
        "source_units": "nm",
        "conversion_to_angstrom": 10.0,
        "mandatory_v4_atom_count": len(
            MANDATORY_V4_ATOMS
        ),
        "shared_real_atom_count": len(
            shared_records
        ),
        "shared_alignment_tolerance_A": (
            SHARED_ATOM_TOLERANCE_A
        ),
        "maximum_shared_coordinate_delta_A": (
            max_shared_delta
        ),
        "shared_alignment_pass": (
            shared_alignment_pass
        ),
        "restored_bond_gate_pass": (
            restored_bond_gate
        ),
        "source_coordinate_labels": (
            source_coordinate_labels
        ),
        "overall_pass": overall_pass,
        "files": {
            "shared_coordinate_csv": str(
                SHARED_CSV.relative_to(ROOT)
            ),
            "restored_bonds_csv": str(
                BONDS_CSV.relative_to(ROOT)
            ),
        },
        "files_sha256": {
            "source_coordinates": sha256(
                SOURCE_COORDINATES
            ),
            "v3_xyz": sha256(V3_XYZ),
            "v3_map": sha256(V3_MAP),
            "shared_coordinate_csv": sha256(
                SHARED_CSV
            ),
            "restored_bonds_csv": sha256(
                BONDS_CSV
            ),
        },
        "authorization": {
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
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 78)
    print("QM_F06 UPPER V4 COORDINATE-ALIGNMENT AUDIT")
    print("=" * 78)
    print(
        "Source:",
        SOURCE_COORDINATES,
    )
    print(
        "Mandatory V4 atoms:",
        len(MANDATORY_V4_ATOMS),
    )
    print(
        "Shared real atoms checked:",
        len(shared_records),
    )
    print(
        "Maximum shared coordinate delta A:",
        max_shared_delta,
    )
    print(
        "Shared alignment:",
        "PASS"
        if shared_alignment_pass
        else "FAIL",
    )
    print(
        "Restored bond gate:",
        "PASS"
        if restored_bond_gate
        else "FAIL",
    )

    print()
    print("RESTORED BONDS")

    for row in bond_records:
        print(
            f"{row['first_atom']:22s} -- "
            f"{row['second_atom']:22s} "
            f"{row['distance_A']:.8f} Å "
            f"{row['bond_class']:12s} "
            f"{'PASS' if row['pass'] else 'FAIL'}"
        )

    print()
    print("Decision:", decision)
    print("Report:", REPORT_PATH)
    print("Shared CSV:", SHARED_CSV)
    print("Bonds CSV:", BONDS_CSV)
    print(
        "V4 construction authorized:",
        overall_pass,
    )
    print("ORCA execution authorized: False")


if __name__ == "__main__":
    main()
