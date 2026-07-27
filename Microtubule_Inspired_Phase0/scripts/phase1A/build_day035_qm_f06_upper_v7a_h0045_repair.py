#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import csv
import hashlib
import json
import math

ROOT = Path.cwd()

SOURCE_DIR = (
    ROOT
    / "runs/phase1A/"
    "day035_qm_f06_upper_v7a_axial_embedding"
)

SOURCE_XYZ = (
    SOURCE_DIR
    / "QM_F06_UPPER_V7A_start.xyz"
)

SOURCE_LOCAL = (
    SOURCE_DIR
    / "QM_F06_UPPER_V7A_best_local_coordinates.csv"
)

SOURCE_REPORT = (
    SOURCE_DIR
    / "QM_F06_UPPER_V7A_AXIAL_EMBEDDING.json"
)

OUTPUT_DIR = (
    ROOT
    / "runs/phase1A/"
    "day035_qm_f06_upper_v7a_h0045_repair"
)

OUTPUT_XYZ = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_H0045_REPAIRED_start.xyz"
)

OUTPUT_LOCAL = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_H0045_REPAIRED_local_coordinates.csv"
)

OUTPUT_REPORT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_H0045_REPAIR.json"
)

OWNER = "A:UPPER:14:0"
INTERNAL_NEIGHBOR = "A:UPPER:13:1"
CUT_CAP = "HCAPV7:UPPER:A14_0:A13_M1"
TARGET_H = "H4:UPPER:0045:0"

TARGET_BH_A = 1.19
MINIMUM_NONOWNER_HEAVY_DISTANCE_A = 1.60
MINIMUM_HH_DISTANCE_A = 1.35


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def vector(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> tuple[float, float, float]:
    return tuple(
        a - b
        for a, b in zip(first, second)
    )


def add(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> tuple[float, float, float]:
    return tuple(
        a + b
        for a, b in zip(first, second)
    )


def scale(
    value: tuple[float, float, float],
    factor: float,
) -> tuple[float, float, float]:
    return tuple(
        factor * component
        for component in value
    )


def norm(
    value: tuple[float, float, float],
) -> float:
    return math.sqrt(
        sum(component * component for component in value)
    )


def unit(
    value: tuple[float, float, float],
) -> tuple[float, float, float]:
    length = norm(value)

    if length <= 1.0e-12:
        raise RuntimeError(
            "Cannot normalize a near-zero vector."
        )

    return scale(value, 1.0 / length)


def distance(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> float:
    return norm(vector(first, second))


def angle(
    first: tuple[float, float, float],
    center: tuple[float, float, float],
    second: tuple[float, float, float],
) -> float:
    u = vector(first, center)
    v = vector(second, center)

    cosine = sum(
        a * b
        for a, b in zip(u, v)
    ) / (norm(u) * norm(v))

    cosine = max(-1.0, min(1.0, cosine))

    return math.degrees(
        math.acos(cosine)
    )


def read_xyz(path: Path) -> list[dict]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    if len(lines) < 2:
        raise RuntimeError(
            f"Malformed XYZ: {path}"
        )

    count = int(lines[0].strip())
    atom_lines = lines[2:2 + count]

    if len(atom_lines) != count:
        raise RuntimeError(
            "Incomplete XYZ: "
            f"declared={count}, "
            f"parsed={len(atom_lines)}"
        )

    atoms = []

    for index, line in enumerate(atom_lines):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Malformed XYZ line: {line}"
            )

        atoms.append({
            "index_0based": index,
            "element": fields[0],
            "xyz_A": tuple(
                float(value)
                for value in fields[1:4]
            ),
        })

    return atoms


def coordinate_fields(
    fieldnames: list[str],
) -> tuple[str, str, str]:
    candidates = (
        ("x_A", "y_A", "z_A"),
        ("x", "y", "z"),
        (
            "transformed_x_A",
            "transformed_y_A",
            "transformed_z_A",
        ),
    )

    fields = set(fieldnames)

    for candidate in candidates:
        if set(candidate) <= fields:
            return candidate

    raise RuntimeError(
        "Could not identify coordinate columns in "
        f"{fieldnames}"
    )


def main() -> None:
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    for path in (
        SOURCE_XYZ,
        SOURCE_LOCAL,
        SOURCE_REPORT,
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(
                f"Missing source artifact: {path}"
            )

    atoms = read_xyz(SOURCE_XYZ)

    with SOURCE_LOCAL.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(
            reader.fieldnames or []
        )
        local_rows = list(reader)

    if not local_rows:
        raise RuntimeError(
            "Local-coordinate table is empty."
        )

    if "atom_id" not in fieldnames:
        raise RuntimeError(
            "Expected atom_id column in local coordinates."
        )

    x_field, y_field, z_field = (
        coordinate_fields(fieldnames)
    )

    row_by_id = {
        row["atom_id"]: row
        for row in local_rows
    }

    required = {
        OWNER,
        INTERNAL_NEIGHBOR,
        CUT_CAP,
    }

    missing = sorted(
        required - set(row_by_id)
    )

    if missing:
        raise RuntimeError(
            "Missing required local records: "
            + "|".join(missing)
        )

    if TARGET_H in row_by_id:
        raise RuntimeError(
            f"{TARGET_H} already exists in local coordinates."
        )

    def xyz_from_row(
        row: dict[str, str],
    ) -> tuple[float, float, float]:
        return (
            float(row[x_field]),
            float(row[y_field]),
            float(row[z_field]),
        )

    owner_xyz = xyz_from_row(
        row_by_id[OWNER]
    )

    internal_xyz = xyz_from_row(
        row_by_id[INTERNAL_NEIGHBOR]
    )

    cap_xyz = xyz_from_row(
        row_by_id[CUT_CAP]
    )

    direction_internal = unit(
        vector(internal_xyz, owner_xyz)
    )

    direction_cut = unit(
        vector(cap_xyz, owner_xyz)
    )

    direction_sum = add(
        direction_internal,
        direction_cut,
    )

    if norm(direction_sum) <= 1.0e-8:
        raise RuntimeError(
            "The two known owner directions are "
            "nearly antiparallel; trigonal completion "
            "is undefined."
        )

    hydrogen_direction = unit(
        scale(direction_sum, -1.0)
    )

    hydrogen_xyz = add(
        owner_xyz,
        scale(
            hydrogen_direction,
            TARGET_BH_A,
        ),
    )

    owner_distance = distance(
        owner_xyz,
        hydrogen_xyz,
    )

    angle_internal_cut = angle(
        internal_xyz,
        owner_xyz,
        cap_xyz,
    )

    angle_internal_h = angle(
        internal_xyz,
        owner_xyz,
        hydrogen_xyz,
    )

    angle_cut_h = angle(
        cap_xyz,
        owner_xyz,
        hydrogen_xyz,
    )

    heavy_distances = []

    for atom in atoms:
        if atom["element"] == "H":
            continue

        value = distance(
            hydrogen_xyz,
            atom["xyz_A"],
        )

        heavy_distances.append({
            "index_0based": atom[
                "index_0based"
            ],
            "element": atom["element"],
            "distance_A": value,
        })

    heavy_distances.sort(
        key=lambda row: row["distance_A"]
    )

    hydrogen_distances = []

    for atom in atoms:
        if atom["element"] != "H":
            continue

        value = distance(
            hydrogen_xyz,
            atom["xyz_A"],
        )

        hydrogen_distances.append({
            "index_0based": atom[
                "index_0based"
            ],
            "distance_A": value,
        })

    hydrogen_distances.sort(
        key=lambda row: row["distance_A"]
    )

    nearest_heavy = heavy_distances[0]
    nonowner_heavy = [
        row
        for row in heavy_distances
        if row["distance_A"] > 1.35
    ]

    if not nonowner_heavy:
        raise RuntimeError(
            "Could not identify a nonowner heavy atom."
        )

    nearest_nonowner_heavy = (
        nonowner_heavy[0]
    )

    nearest_existing_hydrogen = (
        hydrogen_distances[0]
        if hydrogen_distances
        else None
    )

    owner_geometry_pass = (
        abs(owner_distance - TARGET_BH_A)
        <= 1.0e-8
    )

    owner_is_nearest_heavy = (
        nearest_heavy["distance_A"]
        <= 1.35
    )

    nonowner_heavy_clearance_pass = (
        nearest_nonowner_heavy[
            "distance_A"
        ]
        >= MINIMUM_NONOWNER_HEAVY_DISTANCE_A
    )

    hydrogen_clearance_pass = (
        nearest_existing_hydrogen is None
        or nearest_existing_hydrogen[
            "distance_A"
        ]
        >= MINIMUM_HH_DISTANCE_A
    )

    trigonal_angle_pass = (
        95.0 <= angle_internal_h <= 145.0
        and 95.0 <= angle_cut_h <= 145.0
    )

    repaired_atoms = list(atoms)

    repaired_atoms.append({
        "index_0based": len(repaired_atoms),
        "element": "H",
        "xyz_A": hydrogen_xyz,
    })

    composition = Counter(
        atom["element"]
        for atom in repaired_atoms
    )

    composition_pass = (
        composition
        == Counter({
            "B": 17,
            "N": 14,
            "H": 21,
        })
    )

    atom_count_pass = (
        len(repaired_atoms) == 52
    )

    with OUTPUT_XYZ.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            f"{len(repaired_atoms)}\n"
        )
        handle.write(
            "QM_F06 UPPER V7-A with canonical "
            "H4:UPPER:0045:0 trigonal completion\n"
        )

        for atom in repaired_atoms:
            x, y, z = atom["xyz_A"]

            handle.write(
                f"{atom['element']:2s} "
                f"{x: .12f} "
                f"{y: .12f} "
                f"{z: .12f}\n"
            )

    new_fieldnames = list(fieldnames)

    optional_fields = (
        "element",
        "node_type",
        "coordinate_source",
    )

    for field in optional_fields:
        if field not in new_fieldnames:
            new_fieldnames.append(field)

    hydrogen_row = {
        field: ""
        for field in new_fieldnames
    }

    hydrogen_row["atom_id"] = TARGET_H
    hydrogen_row["element"] = "H"
    hydrogen_row[
        "node_type"
    ] = "ANNULUS_OUTER_PASSIVANT_H"

    hydrogen_row[
        "coordinate_source"
    ] = (
        "V7A_TRIGONAL_COMPLETION_FROM_"
        "A13_1_AND_A13_M1_CUT_DIRECTION"
    )

    hydrogen_row[x_field] = (
        f"{hydrogen_xyz[0]:.15f}"
    )
    hydrogen_row[y_field] = (
        f"{hydrogen_xyz[1]:.15f}"
    )
    hydrogen_row[z_field] = (
        f"{hydrogen_xyz[2]:.15f}"
    )

    with OUTPUT_LOCAL.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=new_fieldnames,
        )

        writer.writeheader()

        for row in local_rows:
            writer.writerow({
                field: row.get(field, "")
                for field in new_fieldnames
            })

        writer.writerow(hydrogen_row)

    gates = {
        "source_atom_count_51": (
            len(atoms) == 51
        ),
        "source_composition_B17_N14_H20": (
            Counter(
                atom["element"]
                for atom in atoms
            )
            == Counter({
                "B": 17,
                "N": 14,
                "H": 20,
            })
        ),
        "repaired_atom_count_52": (
            atom_count_pass
        ),
        "repaired_composition_B17_N14_H21": (
            composition_pass
        ),
        "H0045_owner_geometry": (
            owner_geometry_pass
        ),
        "H0045_owner_is_nearest_heavy": (
            owner_is_nearest_heavy
        ),
        "H0045_nonnowner_heavy_clearance": (
            nonowner_heavy_clearance_pass
        ),
        "H0045_hydrogen_clearance": (
            hydrogen_clearance_pass
        ),
        "H0045_trigonal_angles": (
            trigonal_angle_pass
        ),
    }

    passed = all(gates.values())

    decision = (
        "QM_F06_UPPER_V7A_H0045_REPAIR_"
        "GATE_PASS_FORMAL_CONSTRUCTION_AUTHORIZED"
        if passed
        else
        "QM_F06_UPPER_V7A_H0045_REPAIR_"
        "GATE_FAIL_REVIEW_REQUIRED"
    )

    report = {
        "model": "QM_F06_UPPER_V7A",
        "decision": decision,
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "repair": {
            "added_atom": TARGET_H,
            "owner_atom": OWNER,
            "construction_basis": (
                "THIRD_TRIGONAL_DIRECTION_FROM_"
                "INTERNAL_BN_AND_BOUNDARY_CUT_DIRECTIONS"
            ),
            "xyz_A": hydrogen_xyz,
            "owner_distance_A": owner_distance,
            "angle_internal_cut_deg": (
                angle_internal_cut
            ),
            "angle_internal_H_deg": (
                angle_internal_h
            ),
            "angle_cut_H_deg": angle_cut_h,
            "nearest_heavy_record": (
                nearest_heavy
            ),
            "nearest_nonowner_heavy_record": (
                nearest_nonowner_heavy
            ),
            "nearest_existing_hydrogen_record": (
                nearest_existing_hydrogen
            ),
        },
        "atom_count": len(repaired_atoms),
        "composition": dict(
            sorted(composition.items())
        ),
        "gates": gates,
        "files": {
            "source_xyz": str(
                SOURCE_XYZ.relative_to(ROOT)
            ),
            "source_local_coordinates": str(
                SOURCE_LOCAL.relative_to(ROOT)
            ),
            "repaired_xyz": str(
                OUTPUT_XYZ.relative_to(ROOT)
            ),
            "repaired_local_coordinates": str(
                OUTPUT_LOCAL.relative_to(ROOT)
            ),
        },
        "sha256": {
            "source_xyz": sha256(SOURCE_XYZ),
            "source_local_coordinates": sha256(
                SOURCE_LOCAL
            ),
            "repaired_xyz": sha256(OUTPUT_XYZ),
            "repaired_local_coordinates": sha256(
                OUTPUT_LOCAL
            ),
        },
        "formal_construction_authorized": (
            passed
        ),
        "orca_authorized": False,
        "RESP_authorized": False,
        "force_field_adoption_authorized": False,
        "MD_authorized": False,
    }

    OUTPUT_REPORT.write_text(
        json.dumps(
            report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 104)
    print("QM_F06 UPPER V7-A H0045 REPAIR")
    print("=" * 104)

    for name, value in gates.items():
        print(
            f"{name:48s}: "
            f"{'PASS' if value else 'FAIL'}"
        )

    print()
    print("Added atom:", TARGET_H)
    print("Owner:", OWNER)
    print(
        "Hydrogen xyz A:",
        tuple(
            round(value, 12)
            for value in hydrogen_xyz
        ),
    )
    print(
        "Owner B-H distance A:",
        owner_distance,
    )
    print(
        "A13:1--A14:0--H angle deg:",
        angle_internal_h,
    )
    print(
        "CutCap--A14:0--H angle deg:",
        angle_cut_h,
    )
    print(
        "Nearest nonowner heavy distance A:",
        nearest_nonowner_heavy[
            "distance_A"
        ],
    )
    print(
        "Nearest existing H distance A:",
        (
            nearest_existing_hydrogen[
                "distance_A"
            ]
            if nearest_existing_hydrogen
            else None
        ),
    )
    print(
        "Composition:",
        dict(sorted(composition.items())),
    )

    print()
    print("Decision:", decision)
    print("XYZ:", OUTPUT_XYZ)
    print("Local coordinates:", OUTPUT_LOCAL)
    print("Report:", OUTPUT_REPORT)
    print()
    print(
        "Formal construction authorized:",
        passed,
    )
    print("ORCA authorized: False")
    print("RESP authorized: False")
    print("MD authorized: False")

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
