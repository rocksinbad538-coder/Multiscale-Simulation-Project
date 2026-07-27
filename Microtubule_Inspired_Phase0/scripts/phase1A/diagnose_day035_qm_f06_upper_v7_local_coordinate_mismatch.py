#!/usr/bin/env python3

from pathlib import Path
import csv
import math
import json


ROOT = Path.cwd()

CANONICAL = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "23_r2_four_atom_hydrogen_coordinate_embedding/"
    "r2_selected_four_atom_full_coordinates.csv"
)

TRANSFORMED = ROOT / (
    "runs/phase1A/day035_qm_f06_upper_v7_coordinate_preflight/"
    "QM_F06_UPPER_V7_transformed_local_coordinates.csv"
)

V6B_START = ROOT / (
    "runs/phase1A/day033_qm_f06_upper_v6b_pre_qm_audit/"
    "QM_F06_UPPER_V6B_start.xyz"
)

V6B_FINAL = ROOT / (
    "runs/phase1A/day034_qm_f06_upper_v6b_post_qm/"
    "QM_F06_UPPER_V6B_FINAL.xyz"
)

V6B_MAP = ROOT / (
    "runs/phase1A/day033_qm_f06_upper_v6b_orca_input/"
    "QM_F06_UPPER_V6B_constraint_map.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day035_qm_f06_upper_v7_local_coordinate_mismatch"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

OUTPUT_JSON = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7_LOCAL_COORDINATE_MISMATCH.json"
)

TARGETS = [
    "A:UPPER:13:1",
    "A:UPPER:14:0",
    "A:UPPER:14:2",
    "A:UPPER:13:3",
    "A:UPPER:11:1",
    "A:UPPER:13:-1",
]

PAIRS = [
    (
        "A:UPPER:13:1",
        "A:UPPER:14:0",
    ),
    (
        "A:UPPER:13:1",
        "A:UPPER:14:2",
    ),
    (
        "A:UPPER:14:2",
        "A:UPPER:13:3",
    ),
    (
        "A:UPPER:13:1",
        "A:UPPER:11:1",
    ),
    (
        "A:UPPER:14:0",
        "A:UPPER:13:-1",
    ),
]


def distance(first, second):
    return math.sqrt(
        sum(
            (a - b) ** 2
            for a, b in zip(first, second)
        )
    )


def read_csv_coordinates(path):
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        return {}

    headers = set(rows[0])

    id_field = next(
        (
            field
            for field in (
                "atom_id",
                "node_id",
            )
            if field in headers
        ),
        None,
    )

    if id_field is None:
        raise RuntimeError(
            f"No atom-ID field in {path}"
        )

    if {
        "transformed_x_A",
        "transformed_y_A",
        "transformed_z_A",
    } <= headers:
        fields = (
            "transformed_x_A",
            "transformed_y_A",
            "transformed_z_A",
        )
        scale = 1.0

    elif {
        "x_A",
        "y_A",
        "z_A",
    } <= headers:
        fields = (
            "x_A",
            "y_A",
            "z_A",
        )
        scale = 1.0

    elif {
        "x_nm",
        "y_nm",
        "z_nm",
    } <= headers:
        fields = (
            "x_nm",
            "y_nm",
            "z_nm",
        )
        scale = 10.0

    else:
        raise RuntimeError(
            f"No recognized coordinates in {path}"
        )

    coordinates = {}

    for row in rows:
        atom_id = row[id_field]

        if not all(
            row.get(field, "").strip()
            for field in fields
        ):
            continue

        coordinates[atom_id] = tuple(
            scale * float(row[field])
            for field in fields
        )

    return coordinates


def read_xyz(path):
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    count = int(lines[0].strip())

    atoms = []

    for line in lines[2:2 + count]:
        fields = line.split()

        atoms.append({
            "element": fields[0],
            "xyz_A": tuple(
                map(float, fields[1:4])
            ),
        })

    if len(atoms) != count:
        raise RuntimeError(
            f"Incomplete XYZ: {path}"
        )

    return atoms


def read_v6b_map(path):
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        rows = list(csv.DictReader(handle))

    index_field = next(
        (
            field
            for field in (
                "v6b_index_0based",
                "v6a_index_0based",
                "index_0based",
            )
            if rows and field in rows[0]
        ),
        None,
    )

    if index_field is None:
        raise RuntimeError(
            "No usable index field in V6-B map"
        )

    retained_field = next(
        (
            field
            for field in (
                "v6b_retained",
                "v6a_retained",
            )
            if rows and field in rows[0]
        ),
        None,
    )

    if retained_field is not None:
        rows = [
            row
            for row in rows
            if row[retained_field]
            .strip()
            .lower()
            == "true"
        ]

    rows.sort(
        key=lambda row: int(
            row[index_field]
        )
    )

    return rows


def map_xyz_to_ids(xyz_path, map_path):
    atoms = read_xyz(xyz_path)
    rows = read_v6b_map(map_path)

    if len(atoms) != len(rows):
        raise RuntimeError(
            f"XYZ/map count mismatch for {xyz_path}: "
            f"{len(atoms)} vs {len(rows)}"
        )

    coordinates = {}

    for xyz_index, (
        atom,
        row,
    ) in enumerate(zip(atoms, rows)):
        if atom["element"] != row["element"]:
            raise RuntimeError(
                "Element mismatch at sequential "
                f"XYZ index {xyz_index}: "
                f"{atom['element']} vs "
                f"{row['element']} "
                f"for {row['atom_id']}"
            )

        coordinates[row["atom_id"]] = (
            atom["xyz_A"]
        )

    return coordinates


sources = {
    "canonical_source": (
        read_csv_coordinates(CANONICAL)
    ),
    "transformed_v7_local": (
        read_csv_coordinates(TRANSFORMED)
    ),
    "v6b_start": (
        map_xyz_to_ids(
            V6B_START,
            V6B_MAP,
        )
    ),
    "v6b_final_rejected": (
        map_xyz_to_ids(
            V6B_FINAL,
            V6B_MAP,
        )
    ),
}

report = {
    "decision": (
        "QM_F06_UPPER_V7_LOCAL_COORDINATE_"
        "MISMATCH_DIAGNOSED_REPAIR_DESIGN_REQUIRED"
    ),
    "sources": {},
    "orca_authorized": False,
    "RESP_authorized": False,
    "MD_authorized": False,
}

print("=" * 118)
print("QM_F06 UPPER V7 LOCAL COORDINATE-MISMATCH DIAGNOSTIC")
print("=" * 118)

for source_name, coordinates in sources.items():
    print()
    print("SOURCE:", source_name)
    print("Available targets:")

    available = [
        atom_id
        for atom_id in TARGETS
        if atom_id in coordinates
    ]

    for atom_id in available:
        xyz = coordinates[atom_id]

        print(
            f"  {atom_id:28s} "
            f"({xyz[0]: .6f}, "
            f"{xyz[1]: .6f}, "
            f"{xyz[2]: .6f})"
        )

    print()
    print("Pair distances:")

    pair_records = []

    for first, second in PAIRS:
        if (
            first not in coordinates
            or second not in coordinates
        ):
            print(
                f"  {first:28s} -- "
                f"{second:28s} unavailable"
            )
            continue

        value = distance(
            coordinates[first],
            coordinates[second],
        )

        pair_records.append({
            "first_atom": first,
            "second_atom": second,
            "distance_A": value,
        })

        print(
            f"  {first:28s} -- "
            f"{second:28s} "
            f"{value:9.6f} Å"
        )

    report["sources"][source_name] = {
        "available_targets": available,
        "pair_distances": pair_records,
    }

OUTPUT_JSON.write_text(
    json.dumps(
        report,
        indent=2,
    ) + "\n",
    encoding="utf-8",
)

print()
print("Decision:", report["decision"])
print("Report:", OUTPUT_JSON)
print()
print("ORCA authorized: False")
print("RESP authorized: False")
