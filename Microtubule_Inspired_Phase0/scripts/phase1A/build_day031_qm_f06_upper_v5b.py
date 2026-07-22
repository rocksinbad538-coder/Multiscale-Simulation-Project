#!/usr/bin/env python3
"""
Construct the QM_F06 UPPER V5-B starting geometry.

Starting point:
- validated V4 starting geometry and provenance map.

Changes:
- remove two superseded S:1738 artificial caps;
- add BR4:UPPER:14:1, P:1637 and S:1737;
- add three new caps along canonical external cut vectors;
- rebuild H4:UPPER:0203:0 using the missing trigonal direction
  around S:1739, away from P:1641 and P:1639;
- annotate A:UPPER:14:4 for later release in V5-B constraints.

This script authorizes only the pre-QM structural audit.
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

V4_DIR = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_construction"
)

V4_XYZ = V4_DIR / "QM_F06_UPPER_V4_start.xyz"

V4_MAP = (
    V4_DIR
    / "QM_F06_UPPER_V4_atom_role_provenance_map.csv"
)

STRATEGY_REPORT = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_strategy/"
    "QM_F06_UPPER_V5B_STRATEGY.json"
)

COORDINATE_REPORT = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5_coordinate_sources/"
    "QM_F06_UPPER_V5_COORDINATE_SOURCE_AUDIT.json"
)

TRANSFORMED_TARGETS = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5_coordinate_sources/"
    "QM_F06_UPPER_V5_selected_transformed_coordinates.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_construction"
)

OUTPUT_XYZ = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5B_start.xyz"
)

OUTPUT_MAP = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5B_atom_role_provenance_map.csv"
)

OUTPUT_CAPS = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5B_new_artificial_caps.csv"
)

OUTPUT_REPORT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5B_CONSTRUCTION_REPORT.json"
)

REMOVED_CAPS = {
    "HCAPV4:UPPER:S1738:BR4_14_1",
    "HCAPV4:UPPER:S1738:P1637",
}

ADDED_REAL_ORDER = [
    "BR4:UPPER:14:1",
    "P:1637",
    "S:1737",
]

ADDED_CANONICAL_H_ORDER = [
    "H4:UPPER:0170:0",
    "H4:UPPER:0202:0",
]

REBUILT_H = "H4:UPPER:0203:0"

NEW_CAP_SPECS = [
    (
        "HCAPV5B:UPPER:BR4_14_1:BR4_14_2",
        "BR4:UPPER:14:1",
        "BR4:UPPER:14:2",
        1.01,
    ),
    (
        "HCAPV5B:UPPER:P1637:P1636",
        "P:1637",
        "P:1636",
        1.01,
    ),
    (
        "HCAPV5B:UPPER:S1737:P1635",
        "S:1737",
        "P:1635",
        1.19,
    ),
]

EXPECTED_ATOM_COUNT = 52

EXPECTED_COMPOSITION = Counter({
    "B": 16,
    "N": 13,
    "H": 23,
})


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

    for index, line in enumerate(
        lines[2:2 + declared]
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

    if len(atoms) != declared:
        raise RuntimeError(
            f"Incomplete XYZ: expected {declared}, "
            f"found {len(atoms)}"
        )

    return atoms


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


def scale(vector, factor):
    return tuple(
        value * factor
        for value in vector
    )


def norm(vector):
    return math.sqrt(sum(
        value * value
        for value in vector
    ))


def normalize(vector):
    value = norm(vector)

    if value <= 0.0:
        raise RuntimeError(
            "Cannot normalize a zero-length vector."
        )

    return scale(vector, 1.0 / value)


def distance(first, second):
    return norm(subtract(first, second))


def matvec(matrix, vector):
    return tuple(
        sum(
            float(matrix[row][column])
            * float(vector[column])
            for column in range(3)
        )
        for row in range(3)
    )


def main() -> None:
    for path in (
        V4_XYZ,
        V4_MAP,
        STRATEGY_REPORT,
        COORDINATE_REPORT,
        TRANSFORMED_TARGETS,
    ):
        require_file(path)

    strategy = json.loads(
        STRATEGY_REPORT.read_text(
            encoding="utf-8"
        )
    )

    if not strategy["authorization"][
        "v5b_geometry_construction_authorized"
    ]:
        raise RuntimeError(
            "V5-B geometry construction is not authorized."
        )

    coordinate_report = json.loads(
        COORDINATE_REPORT.read_text(
            encoding="utf-8"
        )
    )

    if not coordinate_report["authorization"][
        "transformed_coordinates_authorized"
    ]:
        raise RuntimeError(
            "Transformed V5 coordinates are not authorized."
        )

    selected_source_relative = Path(
        coordinate_report["selected_source"]
    )

    selected_source = ROOT / selected_source_relative

    require_file(selected_source)

    source_rows = read_csv(selected_source)

    source_headers = set(source_rows[0])

    required_source_columns = {
        "node_id",
        "element",
        "node_type",
        "x_nm",
        "y_nm",
        "z_nm",
    }

    if not required_source_columns.issubset(
        source_headers
    ):
        raise RuntimeError(
            "Selected coordinate source has an "
            "unexpected schema."
        )

    rotation = coordinate_report[
        "selected_source_details"
    ]["rotation_matrix"]

    translation = tuple(
        coordinate_report[
            "selected_source_details"
        ]["translation_A"]
    )

    canonical_source = {}

    for row in source_rows:
        atom_id = row["node_id"].strip()

        if not atom_id:
            continue

        values = [
            row["x_nm"].strip(),
            row["y_nm"].strip(),
            row["z_nm"].strip(),
        ]

        if any(not value for value in values):
            continue

        original_A = tuple(
            float(value) * 10.0
            for value in values
        )

        transformed_A = add(
            matvec(rotation, original_A),
            translation,
        )

        canonical_source[atom_id] = {
            "element": row["element"],
            "node_type": row["node_type"],
            "xyz_A": transformed_A,
        }

    required_source_ids = (
        set(ADDED_REAL_ORDER)
        | set(ADDED_CANONICAL_H_ORDER)
        | {
            center
            for _, center, _, _ in NEW_CAP_SPECS
        }
        | {
            outside
            for _, _, outside, _ in NEW_CAP_SPECS
        }
    )

    missing_source = (
        required_source_ids
        - set(canonical_source)
    )

    if missing_source:
        raise RuntimeError(
            "Missing canonical source atoms: "
            f"{sorted(missing_source)}"
        )

    transformed_rows = read_csv(
        TRANSFORMED_TARGETS
    )

    target_coordinates = {
        row["atom_id"]: (
            float(row["x_A"]),
            float(row["y_A"]),
            float(row["z_A"]),
        )
        for row in transformed_rows
    }

    if set(target_coordinates) != set(
        ADDED_REAL_ORDER
    ):
        raise RuntimeError(
            "Transformed target inventory mismatch."
        )

    v4_atoms = read_xyz(V4_XYZ)
    v4_map_rows = read_csv(V4_MAP)

    if len(v4_atoms) != len(v4_map_rows):
        raise RuntimeError(
            "V4 XYZ and provenance map row counts differ."
        )

    atoms = []
    provenance = []

    removed_found = set()
    rebuilt_h_found = False

    # Retain V4 atoms in deterministic order, excluding
    # the two superseded S:1738 caps.
    for atom, row in zip(
        v4_atoms,
        v4_map_rows,
        strict=True,
    ):
        atom_id = row["atom_id"]

        if atom_id in REMOVED_CAPS:
            removed_found.add(atom_id)
            continue

        xyz_A = atom["xyz_A"]
        coordinate_source = (
            "RETAINED_QM_F06_UPPER_V4_START"
        )

        if atom_id == REBUILT_H:
            rebuilt_h_found = True
            coordinate_source = (
                "V5B_REBUILT_TRIGONAL_MISSING_DIRECTION"
            )

        atoms.append({
            "atom_id": atom_id,
            "element": atom["element"],
            "xyz_A": xyz_A,
        })

        provenance.append({
            "atom_id": atom_id,
            "element": atom["element"],
            "atom_role": row["atom_role"],
            "node_type": row["node_type"],
            "coordinate_source": coordinate_source,
            "source_center": row["source_center"],
            "source_outside_atom": (
                row["source_outside_atom"]
            ),
            "source_cut_edge": row["source_cut_edge"],
            "artificial_cap": row["artificial_cap"],
            "v5b_constraint_status": "UNASSIGNED",
            "v5b_mobility_basis": (
                "PENDING_PRE_QM_DESIGN"
            ),
        })

    if removed_found != REMOVED_CAPS:
        raise RuntimeError(
            "Not all superseded caps were found: "
            f"{sorted(REMOVED_CAPS - removed_found)}"
        )

    if not rebuilt_h_found:
        raise RuntimeError(
            f"Rebuilt H not found in V4: {REBUILT_H}"
        )

    id_to_position = {
        atom["atom_id"]: index
        for index, atom in enumerate(atoms)
    }

    for required_id in (
        "S:1739",
        "P:1641",
        "P:1639",
        REBUILT_H,
    ):
        if required_id not in id_to_position:
            raise RuntimeError(
                f"Missing V5-B local atom: {required_id}"
            )

    # Rebuild the S:1739-H direction as the missing
    # trigonal direction opposite to the two heavy neighbors.
    center_xyz = atoms[
        id_to_position["S:1739"]
    ]["xyz_A"]

    p1641_xyz = atoms[
        id_to_position["P:1641"]
    ]["xyz_A"]

    p1639_xyz = atoms[
        id_to_position["P:1639"]
    ]["xyz_A"]

    unit_to_p1641 = normalize(
        subtract(p1641_xyz, center_xyz)
    )

    unit_to_p1639 = normalize(
        subtract(p1639_xyz, center_xyz)
    )

    missing_direction = normalize(
        scale(
            add(
                unit_to_p1641,
                unit_to_p1639,
            ),
            -1.0,
        )
    )

    rebuilt_h_xyz = add(
        center_xyz,
        scale(missing_direction, 1.19),
    )

    atoms[
        id_to_position[REBUILT_H]
    ]["xyz_A"] = rebuilt_h_xyz

    rebuilt_metrics = {
        "S1739_H_A": distance(
            center_xyz,
            rebuilt_h_xyz,
        ),
        "P1641_H_A": distance(
            p1641_xyz,
            rebuilt_h_xyz,
        ),
        "P1639_H_A": distance(
            p1639_xyz,
            rebuilt_h_xyz,
        ),
    }

    retained_ids = {
        atom["atom_id"]
        for atom in atoms
    }

    overlap = (
        retained_ids
        & set(ADDED_REAL_ORDER)
    )

    if overlap:
        raise RuntimeError(
            "Selected V5 atoms already exist in V4: "
            f"{sorted(overlap)}"
        )

    # Append selected V5 real atoms.
    for atom_id in ADDED_REAL_ORDER:
        record = canonical_source[atom_id]

        xyz_A = target_coordinates[atom_id]

        atoms.append({
            "atom_id": atom_id,
            "element": record["element"],
            "xyz_A": xyz_A,
        })

        provenance.append({
            "atom_id": atom_id,
            "element": record["element"],
            "atom_role": (
                "RESTORED_CANONICAL_R2_ATOM"
            ),
            "node_type": record["node_type"],
            "coordinate_source": (
                "V5_FIXED_ANCHOR_TRANSFORMED_"
                + str(selected_source_relative)
            ),
            "source_center": "",
            "source_outside_atom": "",
            "source_cut_edge": "",
            "artificial_cap": "False",
            "v5b_constraint_status": "UNASSIGNED",
            "v5b_mobility_basis": (
                "PENDING_PRE_QM_DESIGN"
            ),
        })

    # Append the two canonical passivant H atoms required
    # by the selected R2 chemical graph.
    canonical_h_centers = {
        "H4:UPPER:0170:0": "BR4:UPPER:14:1",
        "H4:UPPER:0202:0": "S:1737",
    }

    for atom_id in ADDED_CANONICAL_H_ORDER:
        record = canonical_source[atom_id]

        if record["element"] != "H":
            raise RuntimeError(
                f"Canonical passivant is not H: {atom_id}"
            )

        atoms.append({
            "atom_id": atom_id,
            "element": "H",
            "xyz_A": record["xyz_A"],
        })

        provenance.append({
            "atom_id": atom_id,
            "element": "H",
            "atom_role": "RESTORED_CANONICAL_PASSIVANT_H",
            "node_type": record["node_type"],
            "coordinate_source": (
                "V5_FIXED_ANCHOR_TRANSFORMED_"
                + str(selected_source_relative)
            ),
            "source_center": canonical_h_centers[atom_id],
            "source_outside_atom": "",
            "source_cut_edge": "",
            "artificial_cap": "False",
            "v5b_constraint_status": "UNASSIGNED",
            "v5b_mobility_basis": (
                "PENDING_PRE_QM_DESIGN"
            ),
        })

    # Append caps for the three new external cuts.
    cap_records = []

    for (
        cap_id,
        center_id,
        outside_id,
        bond_length_A,
    ) in NEW_CAP_SPECS:
        center = canonical_source[
            center_id
        ]["xyz_A"]

        # Use the exact target coordinate for the newly
        # added center, keeping both sources consistent.
        if center_id in target_coordinates:
            center = target_coordinates[center_id]

        outside = canonical_source[
            outside_id
        ]["xyz_A"]

        direction = normalize(
            subtract(outside, center)
        )

        cap_xyz_A = add(
            center,
            scale(direction, bond_length_A),
        )

        realized_A = distance(
            center,
            cap_xyz_A,
        )

        center_element = canonical_source[
            center_id
        ]["element"]

        outside_element = canonical_source[
            outside_id
        ]["element"]

        atoms.append({
            "atom_id": cap_id,
            "element": "H",
            "xyz_A": cap_xyz_A,
        })

        provenance.append({
            "atom_id": cap_id,
            "element": "H",
            "atom_role": (
                "ARTIFICIAL_BOUNDARY_CAP"
            ),
            "node_type": "QM_BOUNDARY_CAP_H",
            "coordinate_source": (
                "V5_TRANSFORMED_CANONICAL_"
                "CUT_BOND_VECTOR"
            ),
            "source_center": center_id,
            "source_outside_atom": outside_id,
            "source_cut_edge": (
                f"{center_id}--{outside_id}"
            ),
            "artificial_cap": "True",
            "v5b_constraint_status": "UNASSIGNED",
            "v5b_mobility_basis": (
                "PENDING_PRE_QM_DESIGN"
            ),
        })

        cap_records.append({
            "cap_id": cap_id,
            "center_atom": center_id,
            "center_element": center_element,
            "omitted_outside_atom": outside_id,
            "omitted_outside_element": (
                outside_element
            ),
            "target_cap_bond_length_A": (
                bond_length_A
            ),
            "realized_cap_bond_length_A": (
                realized_A
            ),
            "cap_x_A": cap_xyz_A[0],
            "cap_y_A": cap_xyz_A[1],
            "cap_z_A": cap_xyz_A[2],
        })

    for index, record in enumerate(provenance):
        record["index_0based"] = index

    atom_ids = [
        atom["atom_id"]
        for atom in atoms
    ]

    if len(atom_ids) != len(set(atom_ids)):
        duplicates = [
            atom_id
            for atom_id, count in Counter(
                atom_ids
            ).items()
            if count > 1
        ]

        raise RuntimeError(
            f"Duplicate atom IDs: {duplicates}"
        )

    composition = Counter(
        atom["element"]
        for atom in atoms
    )

    gates = {
        "strategy_authorized": True,
        "removed_superseded_caps": (
            removed_found == REMOVED_CAPS
        ),
        "selected_real_atoms_added": (
            set(ADDED_REAL_ORDER)
            .issubset(atom_ids)
        ),
        "new_caps_constructed": (
            len(cap_records) == 3
            and all(
                abs(
                    row[
                        "realized_cap_bond_length_A"
                    ]
                    - row[
                        "target_cap_bond_length_A"
                    ]
                )
                <= 1.0e-10
                for row in cap_records
            )
        ),
        "rebuilt_H_S1739_bond": (
            abs(
                rebuilt_metrics["S1739_H_A"]
                - 1.19
            )
            <= 1.0e-10
        ),
        "rebuilt_H_away_from_P1641": (
            rebuilt_metrics["P1641_H_A"]
            >= 2.0
        ),
        "rebuilt_H_away_from_P1639": (
            rebuilt_metrics["P1639_H_A"]
            >= 2.0
        ),
        "atom_count": (
            len(atoms) == EXPECTED_ATOM_COUNT
        ),
        "composition": (
            composition == EXPECTED_COMPOSITION
        ),
    }

    overall_pass = all(gates.values())

    decision = (
        "QM_F06_UPPER_V5B_START_GEOMETRY_"
        "CONSTRUCTED_PRE_QM_AUDIT_AUTHORIZED"
        if overall_pass
        else
        "QM_F06_UPPER_V5B_CONSTRUCTION_FAIL"
    )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    with OUTPUT_XYZ.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(f"{len(atoms)}\n")
        handle.write(
            "QM_F06 UPPER V5-B selective expansion, "
            "E2915 mobility repair and rebuilt "
            "S1739-H direction\n"
        )

        for atom in atoms:
            x_value, y_value, z_value = (
                atom["xyz_A"]
            )

            handle.write(
                f"{atom['element']:2s} "
                f"{x_value: .12f} "
                f"{y_value: .12f} "
                f"{z_value: .12f}\n"
            )

    map_fieldnames = [
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
        "v5b_constraint_status",
        "v5b_mobility_basis",
    ]

    with OUTPUT_MAP.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=map_fieldnames,
        )
        writer.writeheader()
        writer.writerows(provenance)

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
        "decision": decision,
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "model": "QM_F06_UPPER_V5B",
        "removed_caps": sorted(REMOVED_CAPS),
        "added_real_atoms": ADDED_REAL_ORDER,
        "added_canonical_passivant_H": (
            ADDED_CANONICAL_H_ORDER
        ),
        "new_cap_count": len(cap_records),
        "released_atom_for_constraint_design": (
            "A:UPPER:14:4"
        ),
        "rebuilt_passivant": {
            "atom_id": REBUILT_H,
            "center": "S:1739",
            "heavy_neighbors_used": [
                "P:1641",
                "P:1639",
            ],
            "method": (
                "NEGATIVE_NORMALIZED_SUM_OF_"
                "HEAVY_NEIGHBOR_UNIT_VECTORS"
            ),
            "metrics_A": rebuilt_metrics,
        },
        "final_atom_count": len(atoms),
        "composition": dict(
            sorted(composition.items())
        ),
        "expected_atom_count": (
            EXPECTED_ATOM_COUNT
        ),
        "expected_composition": dict(
            sorted(EXPECTED_COMPOSITION.items())
        ),
        "gates": gates,
        "overall_pass": overall_pass,
        "files": {
            "xyz": str(
                OUTPUT_XYZ.relative_to(ROOT)
            ),
            "atom_map": str(
                OUTPUT_MAP.relative_to(ROOT)
            ),
            "new_caps": str(
                OUTPUT_CAPS.relative_to(ROOT)
            ),
        },
        "files_sha256": {
            "xyz": sha256(OUTPUT_XYZ),
            "atom_map": sha256(OUTPUT_MAP),
            "new_caps": sha256(OUTPUT_CAPS),
            "strategy_report": sha256(
                STRATEGY_REPORT
            ),
            "coordinate_report": sha256(
                COORDINATE_REPORT
            ),
        },
        "authorization": {
            "v5b_start_geometry_constructed": (
                overall_pass
            ),
            "pre_qm_structural_audit_authorized": (
                overall_pass
            ),
            "constraint_design_authorized": False,
            "orca_input_preparation_authorized": False,
            "orca_execution_authorized": False,
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    OUTPUT_REPORT.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 90)
    print("QM_F06 UPPER V5-B CONSTRUCTION")
    print("=" * 90)

    for gate, passed in gates.items():
        print(
            f"{gate:38s}: "
            f"{'PASS' if passed else 'FAIL'}"
        )

    print()
    print("Removed caps:", sorted(REMOVED_CAPS))
    print("Added real atoms:", ADDED_REAL_ORDER)
    print(
        "Added canonical passivant H:",
        ADDED_CANONICAL_H_ORDER,
    )
    print("New artificial caps:", len(cap_records))
    print("Final atom count:", len(atoms))
    print("Composition:", dict(composition))
    print()
    print(
        "Rebuilt S1739-H A:",
        rebuilt_metrics["S1739_H_A"],
    )
    print(
        "Rebuilt P1641-H A:",
        rebuilt_metrics["P1641_H_A"],
    )
    print(
        "Rebuilt P1639-H A:",
        rebuilt_metrics["P1639_H_A"],
    )
    print()
    print("Decision:", decision)
    print("XYZ:", OUTPUT_XYZ)
    print("Map:", OUTPUT_MAP)
    print("Caps:", OUTPUT_CAPS)
    print("Report:", OUTPUT_REPORT)
    print()
    print(
        "Pre-QM audit authorized:",
        overall_pass,
    )
    print("ORCA execution authorized: False")


if __name__ == "__main__":
    main()
