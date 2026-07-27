#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import csv
import json
import math


ROOT = Path.cwd()

CANONICAL_COORDINATES = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "23_r2_four_atom_hydrogen_coordinate_embedding/"
    "r2_selected_four_atom_full_coordinates.csv"
)

V6B_XYZ = ROOT / (
    "runs/phase1A/day033_qm_f06_upper_v6b_pre_qm_audit/"
    "QM_F06_UPPER_V6B_start.xyz"
)

V6B_MAP = ROOT / (
    "runs/phase1A/day033_qm_f06_upper_v6b_orca_input/"
    "QM_F06_UPPER_V6B_constraint_map.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day035_qm_f06_upper_v7a_axial_embedding"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

RANKING_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_axial_embedding_ranking.csv"
)

BEST_COORDINATES_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_best_local_coordinates.csv"
)

BEST_XYZ = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_start.xyz"
)

REPORT_JSON = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_AXIAL_EMBEDDING.json"
)

AXIS_FIRST = "A:UPPER:14:2"
AXIS_SECOND = "A:UPPER:13:3"

NEW_REAL_ATOMS = [
    "A:UPPER:13:1",
    "A:UPPER:14:0",
]

OMITTED_OUTSIDE_ATOMS = {
    "A:UPPER:13:1": "A:UPPER:11:1",
    "A:UPPER:14:0": "A:UPPER:13:-1",
}

OLD_CAP_TO_REMOVE = "HCAPV2:UPPER:03"

NEW_CAPS = {
    "HCAPV7:UPPER:A13_1:A11_1": {
        "owner": "A:UPPER:13:1",
        "outside": "A:UPPER:11:1",
        "element": "H",
        "distance_A": 1.01,
    },
    "HCAPV7:UPPER:A14_0:A13_M1": {
        "owner": "A:UPPER:14:0",
        "outside": "A:UPPER:13:-1",
        "element": "H",
        "distance_A": 1.19,
    },
}

INTENDED_NEW_BONDS = {
    tuple(sorted((
        "A:UPPER:13:1",
        "A:UPPER:14:2",
    ))),
    tuple(sorted((
        "A:UPPER:13:1",
        "A:UPPER:14:0",
    ))),
    tuple(sorted((
        "A:UPPER:13:1",
        "HCAPV7:UPPER:A13_1:A11_1",
    ))),
    tuple(sorted((
        "A:UPPER:14:0",
        "HCAPV7:UPPER:A14_0:A13_M1",
    ))),
}

ANGLE_SAMPLES = 1440

MINIMUM_HEAVY_NONBOND_DISTANCE_A = 1.90
MINIMUM_H_HEAVY_NONBOND_DISTANCE_A = 1.45
MINIMUM_H_H_DISTANCE_A = 1.20


def add(a, b):
    return tuple(
        x + y
        for x, y in zip(a, b)
    )


def subtract(a, b):
    return tuple(
        x - y
        for x, y in zip(a, b)
    )


def scale(vector, factor):
    return tuple(
        factor * value
        for value in vector
    )


def dot(a, b):
    return sum(
        x * y
        for x, y in zip(a, b)
    )


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(vector):
    return math.sqrt(dot(vector, vector))


def normalize(vector):
    value = norm(vector)

    if value <= 1.0e-14:
        raise RuntimeError(
            "Cannot normalize a zero-length vector"
        )

    return scale(vector, 1.0 / value)


def distance(a, b):
    return norm(subtract(a, b))


def matrix_vector(matrix, vector):
    return tuple(
        sum(
            matrix[i][j] * vector[j]
            for j in range(3)
        )
        for i in range(3)
    )


def matrix_multiply(first, second):
    return tuple(
        tuple(
            sum(
                first[i][k] * second[k][j]
                for k in range(3)
            )
            for j in range(3)
        )
        for i in range(3)
    )


def identity_matrix():
    return (
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    )


def rotation_axis_angle(axis, angle):
    x, y, z = normalize(axis)

    cosine = math.cos(angle)
    sine = math.sin(angle)
    complement = 1.0 - cosine

    return (
        (
            cosine + x * x * complement,
            x * y * complement - z * sine,
            x * z * complement + y * sine,
        ),
        (
            y * x * complement + z * sine,
            cosine + y * y * complement,
            y * z * complement - x * sine,
        ),
        (
            z * x * complement - y * sine,
            z * y * complement + x * sine,
            cosine + z * z * complement,
        ),
    )


def align_vector_rotation(source, target):
    source_unit = normalize(source)
    target_unit = normalize(target)

    cosine = max(
        -1.0,
        min(
            1.0,
            dot(source_unit, target_unit),
        ),
    )

    if cosine >= 1.0 - 1.0e-12:
        return identity_matrix()

    if cosine <= -1.0 + 1.0e-12:
        trial = (1.0, 0.0, 0.0)

        if abs(dot(source_unit, trial)) > 0.9:
            trial = (0.0, 1.0, 0.0)

        axis = normalize(
            cross(source_unit, trial)
        )

        return rotation_axis_angle(
            axis,
            math.pi,
        )

    axis = normalize(
        cross(source_unit, target_unit)
    )

    angle = math.acos(cosine)

    return rotation_axis_angle(
        axis,
        angle,
    )


def read_canonical_coordinates(path):
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        rows = list(csv.DictReader(handle))

    coordinates = {}
    elements = {}

    for row in rows:
        atom_id = row["node_id"]

        coordinates[atom_id] = (
            10.0 * float(row["x_nm"]),
            10.0 * float(row["y_nm"]),
            10.0 * float(row["z_nm"]),
        )

        elements[atom_id] = row["element"]

    return coordinates, elements


def read_xyz(path):
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    atom_count = int(lines[0].strip())

    atoms = []

    for line in lines[2:2 + atom_count]:
        fields = line.split()

        atoms.append({
            "element": fields[0],
            "xyz_A": tuple(
                map(float, fields[1:4])
            ),
        })

    if len(atoms) != atom_count:
        raise RuntimeError(
            f"Incomplete XYZ file: {path}"
        )

    return atoms


def read_map_rows(path):
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(
            f"Empty map: {path}"
        )

    index_field = next(
        (
            field
            for field in (
                "v6b_index_0based",
                "v6a_index_0based",
                "index_0based",
            )
            if field in rows[0]
        ),
        None,
    )

    if index_field is None:
        raise RuntimeError(
            "No usable V6-B index field"
        )

    retained_field = next(
        (
            field
            for field in (
                "v6b_retained",
                "v6a_retained",
            )
            if field in rows[0]
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


def cap_coordinate(
    owner_coordinate,
    outside_coordinate,
    bond_length,
):
    direction = normalize(
        subtract(
            outside_coordinate,
            owner_coordinate,
        )
    )

    return add(
        owner_coordinate,
        scale(direction, bond_length),
    )


def nonbond_threshold(
    first_element,
    second_element,
):
    elements = {
        first_element,
        second_element,
    }

    if elements == {"H"}:
        return MINIMUM_H_H_DISTANCE_A

    if "H" in elements:
        return MINIMUM_H_HEAVY_NONBOND_DISTANCE_A

    return MINIMUM_HEAVY_NONBOND_DISTANCE_A


canonical_xyz, canonical_elements = (
    read_canonical_coordinates(
        CANONICAL_COORDINATES
    )
)

v6b_atoms = read_xyz(V6B_XYZ)
map_rows = read_map_rows(V6B_MAP)

if len(v6b_atoms) != len(map_rows):
    raise RuntimeError(
        "V6-B XYZ/map atom-count mismatch: "
        f"{len(v6b_atoms)} vs {len(map_rows)}"
    )

retained = {}

for xyz_index, (
    atom,
    row,
) in enumerate(zip(v6b_atoms, map_rows)):
    if atom["element"] != row["element"]:
        raise RuntimeError(
            "V6-B map/XYZ element mismatch at "
            f"sequential index {xyz_index}: "
            f"{atom['element']} vs "
            f"{row['element']} "
            f"for {row['atom_id']}"
        )

    retained[row["atom_id"]] = {
        "element": atom["element"],
        "xyz_A": atom["xyz_A"],
    }

for atom_id in (
    AXIS_FIRST,
    AXIS_SECOND,
):
    if atom_id not in retained:
        raise RuntimeError(
            f"Missing V6-B axis atom: {atom_id}"
        )

required_canonical = {
    AXIS_FIRST,
    AXIS_SECOND,
    *NEW_REAL_ATOMS,
    *OMITTED_OUTSIDE_ATOMS.values(),
}

missing_canonical = sorted(
    required_canonical - set(canonical_xyz)
)

if missing_canonical:
    raise RuntimeError(
        "Missing canonical coordinates: "
        + "|".join(missing_canonical)
    )

source_origin = canonical_xyz[AXIS_FIRST]
target_origin = retained[AXIS_FIRST]["xyz_A"]

source_axis_vector = subtract(
    canonical_xyz[AXIS_SECOND],
    canonical_xyz[AXIS_FIRST],
)

target_axis_vector = subtract(
    retained[AXIS_SECOND]["xyz_A"],
    retained[AXIS_FIRST]["xyz_A"],
)

initial_rotation = align_vector_rotation(
    source_axis_vector,
    target_axis_vector,
)

target_axis_unit = normalize(
    target_axis_vector
)


def base_transformed(atom_id):
    relative = subtract(
        canonical_xyz[atom_id],
        source_origin,
    )

    return add(
        target_origin,
        matrix_vector(
            initial_rotation,
            relative,
        ),
    )


base_coordinates = {
    atom_id: base_transformed(atom_id)
    for atom_id in required_canonical
}

records = []

for sample_index in range(ANGLE_SAMPLES):
    angle_radians = (
        2.0
        * math.pi
        * sample_index
        / ANGLE_SAMPLES
    )

    axial_rotation = rotation_axis_angle(
        target_axis_unit,
        angle_radians,
    )

    local_coordinates = {}

    for atom_id in (
        *NEW_REAL_ATOMS,
        *OMITTED_OUTSIDE_ATOMS.values(),
    ):
        relative = subtract(
            base_coordinates[atom_id],
            target_origin,
        )

        local_coordinates[atom_id] = add(
            target_origin,
            matrix_vector(
                axial_rotation,
                relative,
            ),
        )

    candidate_atoms = {}

    for atom_id in NEW_REAL_ATOMS:
        candidate_atoms[atom_id] = {
            "element": canonical_elements[atom_id],
            "xyz_A": local_coordinates[atom_id],
        }

    for cap_id, specification in NEW_CAPS.items():
        owner = specification["owner"]
        outside = specification["outside"]

        candidate_atoms[cap_id] = {
            "element": specification["element"],
            "xyz_A": cap_coordinate(
                local_coordinates[owner],
                local_coordinates[outside],
                specification["distance_A"],
            ),
        }

    combined = {
        atom_id: record
        for atom_id, record in retained.items()
        if atom_id != OLD_CAP_TO_REMOVE
    }

    combined.update(candidate_atoms)

    intended_distances = {
        "A13_1_A14_2_A": distance(
            combined["A:UPPER:13:1"]["xyz_A"],
            combined["A:UPPER:14:2"]["xyz_A"],
        ),
        "A13_1_A14_0_A": distance(
            combined["A:UPPER:13:1"]["xyz_A"],
            combined["A:UPPER:14:0"]["xyz_A"],
        ),
        "A13_1_cap_A": distance(
            combined["A:UPPER:13:1"]["xyz_A"],
            combined[
                "HCAPV7:UPPER:A13_1:A11_1"
            ]["xyz_A"],
        ),
        "A14_0_cap_A": distance(
            combined["A:UPPER:14:0"]["xyz_A"],
            combined[
                "HCAPV7:UPPER:A14_0:A13_M1"
            ]["xyz_A"],
        ),
    }

    intended_geometry_pass = (
        1.25
        <= intended_distances[
            "A13_1_A14_2_A"
        ]
        <= 1.85
        and 1.25
        <= intended_distances[
            "A13_1_A14_0_A"
        ]
        <= 1.85
        and 0.90
        <= intended_distances[
            "A13_1_cap_A"
        ]
        <= 1.20
        and 0.95
        <= intended_distances[
            "A14_0_cap_A"
        ]
        <= 1.35
    )

    minimum_clearance = math.inf
    limiting_pair = ""
    clearance_failures = []

    new_ids = set(candidate_atoms)

    all_ids = list(combined)

    for i, first in enumerate(all_ids):
        for second in all_ids[i + 1:]:
            if (
                first not in new_ids
                and second not in new_ids
            ):
                continue

            pair = tuple(sorted((
                first,
                second,
            )))

            if pair in INTENDED_NEW_BONDS:
                continue

            first_record = combined[first]
            second_record = combined[second]

            value = distance(
                first_record["xyz_A"],
                second_record["xyz_A"],
            )

            threshold = nonbond_threshold(
                first_record["element"],
                second_record["element"],
            )

            clearance = value - threshold

            if clearance < minimum_clearance:
                minimum_clearance = clearance
                limiting_pair = (
                    f"{pair[0]}--{pair[1]}"
                )

            if clearance < 0.0:
                clearance_failures.append({
                    "first_atom": pair[0],
                    "second_atom": pair[1],
                    "distance_A": value,
                    "required_minimum_A": threshold,
                    "clearance_A": clearance,
                })

    valid = (
        intended_geometry_pass
        and not clearance_failures
    )

    record = {
        "sample_index": sample_index,
        "angle_degrees": math.degrees(
            angle_radians
        ),
        **intended_distances,
        "minimum_nonnominal_clearance_A": (
            minimum_clearance
        ),
        "limiting_nonnominal_pair": (
            limiting_pair
        ),
        "clearance_failure_count": len(
            clearance_failures
        ),
        "valid": valid,
    }

    for atom_id, atom_record in candidate_atoms.items():
        safe_id = (
            atom_id
            .replace(":", "_")
            .replace("-", "M")
        )

        x, y, z = atom_record["xyz_A"]

        record[f"{safe_id}_x_A"] = x
        record[f"{safe_id}_y_A"] = y
        record[f"{safe_id}_z_A"] = z

    records.append(record)

valid_records = [
    record
    for record in records
    if record["valid"]
]

valid_records.sort(
    key=lambda record: (
        -record[
            "minimum_nonnominal_clearance_A"
        ],
        record["angle_degrees"],
    )
)

for rank, record in enumerate(
    valid_records,
    start=1,
):
    record["rank"] = rank

all_ranked = (
    valid_records
    + [
        record
        for record in records
        if not record["valid"]
    ]
)

fieldnames = []

for record in all_ranked:
    for key in record:
        if key not in fieldnames:
            fieldnames.append(key)

with RANKING_CSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
    )

    writer.writeheader()
    writer.writerows(all_ranked)

best = (
    valid_records[0]
    if valid_records
    else None
)

if best is not None:
    best_atoms = {
        atom_id: {
            "element": record["element"],
            "xyz_A": record["xyz_A"],
        }
        for atom_id, record in retained.items()
        if atom_id != OLD_CAP_TO_REMOVE
    }

    coordinate_rows = []

    for atom_id in NEW_REAL_ATOMS:
        safe_id = (
            atom_id
            .replace(":", "_")
            .replace("-", "M")
        )

        xyz = (
            float(best[f"{safe_id}_x_A"]),
            float(best[f"{safe_id}_y_A"]),
            float(best[f"{safe_id}_z_A"]),
        )

        best_atoms[atom_id] = {
            "element": canonical_elements[atom_id],
            "xyz_A": xyz,
        }

        coordinate_rows.append({
            "atom_id": atom_id,
            "element": canonical_elements[atom_id],
            "x_A": xyz[0],
            "y_A": xyz[1],
            "z_A": xyz[2],
            "coordinate_role": "NEW_REAL_V7_ATOM",
        })

    for cap_id, specification in NEW_CAPS.items():
        safe_id = (
            cap_id
            .replace(":", "_")
            .replace("-", "M")
        )

        xyz = (
            float(best[f"{safe_id}_x_A"]),
            float(best[f"{safe_id}_y_A"]),
            float(best[f"{safe_id}_z_A"]),
        )

        best_atoms[cap_id] = {
            "element": "H",
            "xyz_A": xyz,
        }

        coordinate_rows.append({
            "atom_id": cap_id,
            "element": "H",
            "x_A": xyz[0],
            "y_A": xyz[1],
            "z_A": xyz[2],
            "coordinate_role": "NEW_V7_BOUNDARY_CAP",
        })

    with BEST_COORDINATES_CSV.open(
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
                "coordinate_role",
            ],
        )

        writer.writeheader()
        writer.writerows(coordinate_rows)

    final_order = [
        row["atom_id"]
        for row in map_rows
        if row["atom_id"] != OLD_CAP_TO_REMOVE
    ]

    final_order.extend(NEW_REAL_ATOMS)
    final_order.extend(NEW_CAPS)

    if len(final_order) != len(set(final_order)):
        raise RuntimeError(
            "Duplicate atom IDs in V7-A order"
        )

    if set(final_order) != set(best_atoms):
        raise RuntimeError(
            "V7-A atom-order/set mismatch"
        )

    composition = Counter(
        best_atoms[atom_id]["element"]
        for atom_id in final_order
    )

    xyz_lines = [
        str(len(final_order)),
        (
            "QM_F06_UPPER_V7A axial canonical "
            "embedding; pre-QM candidate only"
        ),
    ]

    for atom_id in final_order:
        atom = best_atoms[atom_id]
        x, y, z = atom["xyz_A"]

        xyz_lines.append(
            f"{atom['element']:2s} "
            f"{x: .12f} "
            f"{y: .12f} "
            f"{z: .12f}"
        )

    BEST_XYZ.write_text(
        "\n".join(xyz_lines) + "\n",
        encoding="utf-8",
    )

else:
    composition = {}

report = {
    "decision": (
        "QM_F06_UPPER_V7A_AXIAL_EMBEDDING_"
        "CANDIDATE_FOUND_FORMAL_CONSTRUCTION_REQUIRED"
        if best is not None
        else
        "QM_F06_UPPER_V7A_AXIAL_EMBEDDING_"
        "NO_VALID_ORIENTATION"
    ),
    "search": {
        "axis_first": AXIS_FIRST,
        "axis_second": AXIS_SECOND,
        "angle_samples": ANGLE_SAMPLES,
        "valid_candidate_count": len(
            valid_records
        ),
        "embedding_rule": (
            "ALIGN_CANONICAL_A14_2_A13_3_VECTOR_"
            "TO_V6B_START_THEN_ROTATE_LOCAL_PATCH_"
            "AROUND_THE_ALIGNED_AXIS"
        ),
    },
    "model": {
        "old_cap_removed": OLD_CAP_TO_REMOVE,
        "new_real_atoms": NEW_REAL_ATOMS,
        "new_caps": list(NEW_CAPS),
        "expected_atom_count": 51,
        "expected_composition": {
            "B": 17,
            "N": 14,
            "H": 20,
        },
        "realized_composition": dict(
            sorted(composition.items())
        ) if best is not None else None,
    },
    "best_candidate": best,
    "formal_construction_authorized": (
        best is not None
    ),
    "orca_authorized": False,
    "RESP_authorized": False,
    "force_field_adoption_authorized": False,
    "MD_authorized": False,
}

REPORT_JSON.write_text(
    json.dumps(
        report,
        indent=2,
    ) + "\n",
    encoding="utf-8",
)

print("=" * 112)
print("QM_F06 UPPER V7-A AXIAL CANONICAL EMBEDDING SEARCH")
print("=" * 112)
print("Angles evaluated:", ANGLE_SAMPLES)
print("Valid orientations:", len(valid_records))

if best is not None:
    print()
    print("Best orientation:")
    print(
        "  rank:",
        best["rank"],
    )
    print(
        "  angle degrees:",
        best["angle_degrees"],
    )
    print(
        "  A13:1--A14:2 A:",
        best["A13_1_A14_2_A"],
    )
    print(
        "  A13:1--A14:0 A:",
        best["A13_1_A14_0_A"],
    )
    print(
        "  minimum nonnominal clearance A:",
        best[
            "minimum_nonnominal_clearance_A"
        ],
    )
    print(
        "  limiting pair:",
        best[
            "limiting_nonnominal_pair"
        ],
    )
    print(
        "  composition:",
        dict(sorted(composition.items())),
    )
    print()
    print("Best XYZ:", BEST_XYZ)
    print(
        "Best local coordinates:",
        BEST_COORDINATES_CSV,
    )

print()
print("Decision:", report["decision"])
print("Ranking:", RANKING_CSV)
print("Report:", REPORT_JSON)
print()
print(
    "Formal V7 construction authorized:",
    best is not None,
)
print("ORCA authorized: False")
print("RESP authorized: False")

if best is None:
    raise SystemExit(1)
