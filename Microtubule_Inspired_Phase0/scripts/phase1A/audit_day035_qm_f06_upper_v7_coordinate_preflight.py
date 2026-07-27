#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import csv
import json
import math

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

V6B_XYZ = ROOT / (
    "runs/phase1A/"
    "day033_qm_f06_upper_v6b_pre_qm_audit/"
    "QM_F06_UPPER_V6B_start.xyz"
)

V6B_MAP = ROOT / (
    "runs/phase1A/"
    "day033_qm_f06_upper_v6b_orca_input/"
    "QM_F06_UPPER_V6B_constraint_map.csv"
)

CANONICAL_COORDINATES = ROOT / (
    "runs/phase1A/"
    "day024_chemical_end_rim_design/"
    "23_r2_four_atom_hydrogen_coordinate_embedding/"
    "r2_selected_four_atom_full_coordinates.csv"
)

CANONICAL_NODES = ROOT / (
    "runs/phase1A/"
    "day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

CANONICAL_EDGES = ROOT / (
    "runs/phase1A/"
    "day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_edges.csv"
)

BOUNDARY_REPORT = ROOT / (
    "runs/phase1A/"
    "day035_qm_f06_upper_v7_boundary_provenance/"
    "QM_F06_UPPER_V7_BOUNDARY_PROVENANCE.json"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day035_qm_f06_upper_v7_coordinate_preflight"
)

OUTPUT_COORDINATES = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7_transformed_local_coordinates.csv"
)

OUTPUT_CAPS = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7_proposed_caps.csv"
)

OUTPUT_REPORT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7_COORDINATE_PREFLIGHT.json"
)


ADDED_ATOMS = {
    "A:UPPER:13:1",
    "A:UPPER:14:0",
}

OUTSIDE_ATOMS = {
    "A:UPPER:11:1",
    "A:UPPER:13:-1",
}

LOCAL_ATOMS = ADDED_ATOMS | OUTSIDE_ATOMS

INTERNAL_REQUIRED_EDGES = {
    tuple(sorted((
        "A:UPPER:13:1",
        "A:UPPER:14:2",
    ))),
    tuple(sorted((
        "A:UPPER:13:1",
        "A:UPPER:14:0",
    ))),
}

NEW_CUTS = {
    (
        "A:UPPER:13:1",
        "A:UPPER:11:1",
        1.01,
        "HCAPV7:UPPER:A13_1:A11_1",
    ),
    (
        "A:UPPER:14:0",
        "A:UPPER:13:-1",
        1.19,
        "HCAPV7:UPPER:A14_0:A13_M1",
    ),
}

MAXIMUM_ALIGNMENT_RMSD_A = 5.0e-3
MAXIMUM_ALIGNMENT_RESIDUAL_A = 2.0e-2
MINIMUM_BN_DISTANCE_A = 1.25
MAXIMUM_BN_DISTANCE_A = 1.85
MINIMUM_CAP_HEAVY_CLEARANCE_A = 1.55


def read_csv(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty file: {path}"
        )

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def read_xyz(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty XYZ: {path}"
        )

    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    count = int(lines[0].strip())

    atoms = []

    for index, line in enumerate(
        lines[2:2 + count]
    ):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Malformed XYZ line {index + 3}"
            )

        atoms.append({
            "index_0based": index,
            "element": fields[0],
            "xyz_A": np.array(
                [
                    float(fields[1]),
                    float(fields[2]),
                    float(fields[3]),
                ],
                dtype=float,
            ),
        })

    if len(atoms) != count:
        raise RuntimeError(
            "XYZ atom-count mismatch"
        )

    return atoms


def find_field(headers, candidates):
    for candidate in candidates:
        if candidate in headers:
            return candidate

    raise RuntimeError(
        "Could not identify required field. "
        f"Candidates={candidates}; "
        f"headers={sorted(headers)}"
    )


def truthy(value):
    return str(value).strip().lower() in {
        "true",
        "1",
        "yes",
    }


def distance(first, second):
    return float(
        np.linalg.norm(first - second)
    )


def kabsch(source, target):
    source_center = source.mean(axis=0)
    target_center = target.mean(axis=0)

    source_centered = source - source_center
    target_centered = target - target_center

    covariance = (
        source_centered.T
        @ target_centered
    )

    u, _, vt = np.linalg.svd(covariance)

    rotation = vt.T @ u.T

    reflection_corrected = False

    if np.linalg.det(rotation) < 0.0:
        vt[-1, :] *= -1.0
        rotation = vt.T @ u.T
        reflection_corrected = True

    translation = (
        target_center
        - source_center @ rotation.T
    )

    return (
        rotation,
        translation,
        reflection_corrected,
    )


def transform(xyz, rotation, translation):
    return xyz @ rotation.T + translation


def nearest_heavy_clearance(
    point,
    heavy_coordinates,
    excluded_atoms,
):
    records = []

    for atom_id, xyz in heavy_coordinates.items():
        if atom_id in excluded_atoms:
            continue

        records.append((
            distance(point, xyz),
            atom_id,
        ))

    if not records:
        raise RuntimeError(
            "No nonowner heavy atoms available"
        )

    records.sort()

    return records[0], records[1:5]


def main():
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    boundary_report = json.loads(
        BOUNDARY_REPORT.read_text(
            encoding="utf-8"
        )
    )

    if not boundary_report.get(
        "construction_authorized",
        False,
    ):
        raise RuntimeError(
            "Boundary-provenance gate did not "
            "authorize V7 construction"
        )

    map_rows = read_csv(V6B_MAP)
    node_rows = read_csv(CANONICAL_NODES)
    edge_rows = read_csv(CANONICAL_EDGES)
    coordinate_rows = read_csv(
        CANONICAL_COORDINATES
    )
    v6b_atoms = read_xyz(V6B_XYZ)

    map_headers = set(map_rows[0])

    index_field = find_field(
        map_headers,
        (
            "index_0based",
            "v6b_index_0based",
            "v6a_index_0based",
        ),
    )

    fixed_field = find_field(
        map_headers,
        (
            "v6b_fixed",
            "fixed",
            "v5b_fixed",
        ),
    )

    map_rows.sort(
        key=lambda row: int(row[index_field])
    )

    if len(map_rows) != len(v6b_atoms):
        raise RuntimeError(
            "V6-B map/XYZ atom-count mismatch"
        )

    historical_indices = [
        int(row[index_field])
        for row in map_rows
    ]

    historical_index_is_contiguous = (
        historical_indices
        == list(range(len(map_rows)))
    )

    # The V6-B map can retain historical indices with gaps
    # after atom deletion. The current XYZ is necessarily
    # indexed sequentially from 0 to N-1. After sorting the
    # map by the historical index, map each retained row to
    # its current sequential XYZ position.
    id_by_xyz_index = {
        xyz_index: row["atom_id"]
        for xyz_index, row in enumerate(map_rows)
    }

    element_by_id = {
        row["atom_id"]: row["element"]
        for row in map_rows
    }

    element_order_failures = []

    for xyz_index, atom in enumerate(v6b_atoms):
        row = map_rows[xyz_index]

        if atom["element"] != row["element"]:
            element_order_failures.append({
                "xyz_index_0based": xyz_index,
                "historical_index": int(
                    row[index_field]
                ),
                "atom_id": row["atom_id"],
                "map_element": row["element"],
                "xyz_element": atom["element"],
            })

    if element_order_failures:
        raise RuntimeError(
            "V6-B map/XYZ element-order mismatch: "
            + json.dumps(
                element_order_failures[:10]
            )
        )

    v6b_coordinates = {
        id_by_xyz_index[index]: atom["xyz_A"]
        for index, atom in enumerate(v6b_atoms)
    }

    canonical_element = {
        row["node_id"]: row["element"]
        for row in node_rows
    }

    canonical_type = {
        row["node_id"]: row["node_type"]
        for row in node_rows
    }

    canonical_coordinates = {}

    for row in coordinate_rows:
        canonical_coordinates[row["node_id"]] = (
            np.array(
                [
                    float(row["x_nm"]),
                    float(row["y_nm"]),
                    float(row["z_nm"]),
                ],
                dtype=float,
            )
            * 10.0
        )

    graph = defaultdict(set)
    canonical_heavy_edges = set()

    for row in edge_rows:
        first = row["source_node"]
        second = row["target_node"]

        graph[first].add(second)
        graph[second].add(first)

        if (
            row["source_element"] != "H"
            and row["target_element"] != "H"
        ):
            canonical_heavy_edges.add(
                tuple(sorted((first, second)))
            )

    missing_local_coordinates = sorted(
        atom_id
        for atom_id in LOCAL_ATOMS
        if atom_id not in canonical_coordinates
    )

    if missing_local_coordinates:
        raise RuntimeError(
            "Missing canonical local coordinates: "
            + "|".join(missing_local_coordinates)
        )

    fixed_heavy_ids = []

    for row in map_rows:
        atom_id = row["atom_id"]

        if not truthy(row[fixed_field]):
            continue

        if row["element"] == "H":
            continue

        if atom_id not in canonical_coordinates:
            continue

        fixed_heavy_ids.append(atom_id)

    if len(fixed_heavy_ids) < 3:
        raise RuntimeError(
            "Insufficient fixed heavy anchors for "
            "rigid alignment"
        )

    source = np.array([
        canonical_coordinates[atom_id]
        for atom_id in fixed_heavy_ids
    ])

    target = np.array([
        v6b_coordinates[atom_id]
        for atom_id in fixed_heavy_ids
    ])

    (
        rotation,
        translation,
        reflection_corrected,
    ) = kabsch(source, target)

    fitted = np.array([
        transform(
            canonical_coordinates[atom_id],
            rotation,
            translation,
        )
        for atom_id in fixed_heavy_ids
    ])

    residuals = np.linalg.norm(
        fitted - target,
        axis=1,
    )

    alignment_rmsd = float(
        math.sqrt(
            float(np.mean(residuals ** 2))
        )
    )

    alignment_maximum = float(
        np.max(residuals)
    )

    transformed_local = {
        atom_id: transform(
            canonical_coordinates[atom_id],
            rotation,
            translation,
        )
        for atom_id in LOCAL_ATOMS
    }

    combined_coordinates = dict(
        v6b_coordinates
    )

    combined_coordinates.update({
        atom_id: transformed_local[atom_id]
        for atom_id in ADDED_ATOMS
    })

    combined_elements = dict(
        element_by_id
    )

    for atom_id in ADDED_ATOMS:
        combined_elements[atom_id] = (
            canonical_element[atom_id]
        )

    internal_records = []

    for first, second in sorted(
        INTERNAL_REQUIRED_EDGES
    ):
        if (
            first not in combined_coordinates
            or second not in combined_coordinates
        ):
            raise RuntimeError(
                "Missing coordinate for required "
                f"internal edge: {first}--{second}"
            )

        value = distance(
            combined_coordinates[first],
            combined_coordinates[second],
        )

        passed = (
            MINIMUM_BN_DISTANCE_A
            <= value
            <= MAXIMUM_BN_DISTANCE_A
        )

        internal_records.append({
            "first_atom": first,
            "second_atom": second,
            "distance_A": value,
            "pass": passed,
        })

    heavy_coordinates = {
        atom_id: xyz
        for atom_id, xyz
        in combined_coordinates.items()
        if combined_elements.get(atom_id) != "H"
    }

    cap_records = []

    for (
        owner,
        outside,
        bond_length,
        cap_id,
    ) in sorted(NEW_CUTS):
        owner_xyz = transformed_local[owner]
        outside_xyz = transformed_local[outside]

        direction = outside_xyz - owner_xyz
        direction_norm = np.linalg.norm(direction)

        if direction_norm <= 1.0e-12:
            raise RuntimeError(
                f"Zero cut vector: {owner}--{outside}"
            )

        unit = direction / direction_norm

        cap_xyz = (
            owner_xyz
            + bond_length * unit
        )

        nearest, next_nearest = (
            nearest_heavy_clearance(
                cap_xyz,
                heavy_coordinates,
                excluded_atoms={owner},
            )
        )

        nearest_distance, nearest_atom = nearest

        cap_pass = (
            nearest_distance
            >= MINIMUM_CAP_HEAVY_CLEARANCE_A
        )

        cap_records.append({
            "cap_id": cap_id,
            "owner_atom": owner,
            "owner_element": (
                canonical_element[owner]
            ),
            "outside_atom": outside,
            "outside_element": (
                canonical_element[outside]
            ),
            "target_owner_distance_A": (
                bond_length
            ),
            "source_cut_distance_A": (
                distance(owner_xyz, outside_xyz)
            ),
            "cap_x_A": float(cap_xyz[0]),
            "cap_y_A": float(cap_xyz[1]),
            "cap_z_A": float(cap_xyz[2]),
            "nearest_nonowner_heavy_atom": (
                nearest_atom
            ),
            "nearest_nonowner_heavy_distance_A": (
                nearest_distance
            ),
            "minimum_required_clearance_A": (
                MINIMUM_CAP_HEAVY_CLEARANCE_A
            ),
            "pass": cap_pass,
            "next_nearest_heavy_atoms": "|".join(
                f"{atom_id}:{value:.6f}"
                for value, atom_id
                in next_nearest
            ),
        })

    coordinate_records = []

    for atom_id in sorted(LOCAL_ATOMS):
        xyz = transformed_local[atom_id]

        coordinate_records.append({
            "atom_id": atom_id,
            "element": canonical_element[atom_id],
            "node_type": canonical_type[atom_id],
            "v7_role": (
                "ADDED_REAL_ATOM"
                if atom_id in ADDED_ATOMS
                else "OUTSIDE_CUT_REFERENCE"
            ),
            "x_A": float(xyz[0]),
            "y_A": float(xyz[1]),
            "z_A": float(xyz[2]),
        })

    with OUTPUT_COORDINATES.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                coordinate_records[0]
            ),
        )
        writer.writeheader()
        writer.writerows(
            coordinate_records
        )

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

    local_neighbor_inventory = {}

    for atom_id in sorted(ADDED_ATOMS):
        local_neighbor_inventory[atom_id] = [
            {
                "atom_id": neighbor,
                "element": canonical_element.get(
                    neighbor,
                    "",
                ),
                "inside_v7_real_model": (
                    neighbor
                    in combined_coordinates
                    and combined_elements.get(
                        neighbor
                    ) != "H"
                ),
            }
            for neighbor in sorted(
                graph.get(atom_id, set())
            )
            if canonical_element.get(
                neighbor
            ) != "H"
        ]

    gates = {
        "boundary_provenance_authorized": True,
        "fixed_heavy_anchor_count": (
            len(fixed_heavy_ids) >= 3
        ),
        "alignment_rmsd": (
            alignment_rmsd
            <= MAXIMUM_ALIGNMENT_RMSD_A
        ),
        "alignment_maximum_residual": (
            alignment_maximum
            <= MAXIMUM_ALIGNMENT_RESIDUAL_A
        ),
        "all_local_coordinates_available": (
            len(missing_local_coordinates) == 0
        ),
        "required_internal_edges_canonical": (
            INTERNAL_REQUIRED_EDGES
            <= canonical_heavy_edges
        ),
        "required_internal_BN_geometry": all(
            record["pass"]
            for record in internal_records
        ),
        "new_cap_geometry": all(
            record["pass"]
            for record in cap_records
        ),
    }

    passed = all(gates.values())

    report = {
        "model": "QM_F06_UPPER_V7",
        "decision": (
            "QM_F06_UPPER_V7_COORDINATE_"
            "PREFLIGHT_PASS_FORMAL_"
            "CONSTRUCTION_AUTHORIZED"
            if passed
            else
            "QM_F06_UPPER_V7_COORDINATE_"
            "PREFLIGHT_FAIL_REVIEW_REQUIRED"
        ),
        "index_mapping": {
            "source_index_field": index_field,
            "historical_indices": historical_indices,
            "historical_index_is_contiguous": (
                historical_index_is_contiguous
            ),
            "current_xyz_indices": list(
                range(len(map_rows))
            ),
            "mapping_rule": (
                "MAP_ROWS_SORTED_BY_HISTORICAL_INDEX_"
                "THEN_ASSIGNED_TO_SEQUENTIAL_XYZ_INDEX"
            ),
            "element_order_verified": True,
        },
        "alignment": {
            "anchor_count": len(
                fixed_heavy_ids
            ),
            "anchor_atom_ids": (
                fixed_heavy_ids
            ),
            "rmsd_A": alignment_rmsd,
            "maximum_residual_A": (
                alignment_maximum
            ),
            "reflection_corrected": (
                reflection_corrected
            ),
            "rotation_determinant": float(
                np.linalg.det(rotation)
            ),
        },
        "added_atoms": sorted(ADDED_ATOMS),
        "outside_cut_reference_atoms": sorted(
            OUTSIDE_ATOMS
        ),
        "local_neighbor_inventory": (
            local_neighbor_inventory
        ),
        "internal_edge_records": (
            internal_records
        ),
        "cap_records": cap_records,
        "gates": gates,
        "formal_construction_authorized": passed,
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

    print("=" * 112)
    print("QM_F06 UPPER V7 COORDINATE AND CAP PREFLIGHT")
    print("=" * 112)

    for name, value in gates.items():
        print(
            f"{name:48s}: "
            f"{'PASS' if value else 'FAIL'}"
        )

    print()
    print(
        "Fixed heavy anchors:",
        len(fixed_heavy_ids),
    )
    print(
        "Alignment RMSD A:",
        alignment_rmsd,
    )
    print(
        "Alignment maximum residual A:",
        alignment_maximum,
    )
    print(
        "Rotation determinant:",
        float(np.linalg.det(rotation)),
    )
    print(
        "Reflection corrected:",
        reflection_corrected,
    )

    print()
    print("Canonical heavy-neighbor inventory:")

    for atom_id in sorted(
        local_neighbor_inventory
    ):
        print()
        print(
            f"  {atom_id} "
            f"({canonical_element[atom_id]}):"
        )

        for record in (
            local_neighbor_inventory[atom_id]
        ):
            print(
                f"    {record['atom_id']:28s} "
                f"{record['element']:2s} "
                f"inside_V7="
                f"{record['inside_v7_real_model']}"
            )

    print()
    print("Required internal B-N bonds:")

    for record in internal_records:
        print(
            f"  {record['first_atom']:28s} -- "
            f"{record['second_atom']:28s} "
            f"{record['distance_A']:.6f} Å | "
            f"{'PASS' if record['pass'] else 'FAIL'}"
        )

    print()
    print("Proposed caps:")

    for record in cap_records:
        print()
        print(
            f"  {record['cap_id']}"
        )
        print(
            f"    owner: "
            f"{record['owner_atom']} "
            f"({record['owner_element']})"
        )
        print(
            f"    omitted outside atom: "
            f"{record['outside_atom']} "
            f"({record['outside_element']})"
        )
        print(
            f"    owner distance: "
            f"{record['target_owner_distance_A']:.6f} Å"
        )
        print(
            f"    nearest nonowner heavy atom: "
            f"{record['nearest_nonowner_heavy_atom']}"
        )
        print(
            f"    clearance: "
            f"{record['nearest_nonowner_heavy_distance_A']:.6f} Å"
        )
        print(
            f"    gate: "
            f"{'PASS' if record['pass'] else 'FAIL'}"
        )

    print()
    print("Decision:", report["decision"])
    print(
        "Coordinates:",
        OUTPUT_COORDINATES,
    )
    print("Caps:", OUTPUT_CAPS)
    print("Report:", OUTPUT_REPORT)
    print()
    print(
        "Formal V7 construction authorized:",
        passed,
    )
    print("ORCA authorized: False")
    print("RESP authorized: False")

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
