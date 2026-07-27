#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import csv
import json
import math


ROOT = Path.cwd()

V6B_XYZ = (
    ROOT
    / "runs/phase1A/day033_qm_f06_upper_v6b_orca_input"
    / "v6b_start.xyz"
)

V6B_MAP = (
    ROOT
    / "runs/phase1A/day033_qm_f06_upper_v6b_orca_input"
    / "QM_F06_UPPER_V6B_constraint_map.csv"
)

V6A_EDGES = (
    ROOT
    / "runs/phase1A/day033_qm_f06_upper_v6a_topology_closure"
    / "QM_F06_UPPER_V6A_nominal_edges.csv"
)

V7_REPAIRED_XYZ = (
    ROOT
    / "runs/phase1A/day035_qm_f06_upper_v7a_h0045_repair"
    / "QM_F06_UPPER_V7A_H0045_REPAIRED_start.xyz"
)

V7_LOCAL = (
    ROOT
    / "runs/phase1A/day035_qm_f06_upper_v7a_h0045_repair"
    / "QM_F06_UPPER_V7A_H0045_REPAIRED_local_coordinates.csv"
)

V7_BOUNDARY = (
    ROOT
    / "runs/phase1A/day035_qm_f06_upper_v7_boundary_provenance"
    / "QM_F06_UPPER_V7_boundary_cut_provenance.csv"
)

V7_H0045_REPORT = (
    ROOT
    / "runs/phase1A/day035_qm_f06_upper_v7a_h0045_repair"
    / "QM_F06_UPPER_V7A_H0045_REPAIR.json"
)

OUTPUT_DIR = (
    ROOT
    / "runs/phase1A/day035_qm_f06_upper_v7a_formal_construction"
)

OUTPUT_XYZ = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_start.xyz"
)

OUTPUT_MAP = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_atom_role_provenance_map.csv"
)

OUTPUT_EDGES = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_nominal_edges.csv"
)

OUTPUT_CAPS = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_boundary_caps.csv"
)

OUTPUT_REPORT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_CONSTRUCTION_REPORT.json"
)


OBSOLETE_ATOM = "HCAPV2:UPPER:03"

NEW_ATOM_ORDER = [
    "A:UPPER:13:1",
    "A:UPPER:14:0",
    "HCAPV7:UPPER:A13_1:A11_1",
    "HCAPV7:UPPER:A14_0:A13_M1",
    "H4:UPPER:0045:0",
]

EXPECTED_COMPOSITION = {
    "B": 17,
    "N": 14,
    "H": 21,
}

NEW_EDGES = [
    (
        "A:UPPER:13:1",
        "A:UPPER:14:2",
        "CANONICAL_V7_INTERNAL_BN",
    ),
    (
        "A:UPPER:13:1",
        "A:UPPER:14:0",
        "CANONICAL_V7_INTERNAL_BN",
    ),
    (
        "A:UPPER:13:1",
        "HCAPV7:UPPER:A13_1:A11_1",
        "V7_ARTIFICIAL_BOUNDARY_CAP",
    ),
    (
        "A:UPPER:14:0",
        "HCAPV7:UPPER:A14_0:A13_M1",
        "V7_ARTIFICIAL_BOUNDARY_CAP",
    ),
    (
        "A:UPPER:14:0",
        "H4:UPPER:0045:0",
        "CANONICAL_ANNULUS_PASSIVANT",
    ),
]


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty source file: {path}"
        )


def file_sha256(path: Path) -> str:
    digest = sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def read_xyz(
    path: Path,
) -> tuple[list[dict], str]:
    lines = path.read_text(
        encoding="utf-8",
        errors="strict",
    ).splitlines()

    if len(lines) < 2:
        raise RuntimeError(
            f"Incomplete XYZ: {path}"
        )

    count = int(lines[0].strip())
    atom_lines = lines[2:2 + count]

    if len(atom_lines) != count:
        raise RuntimeError(
            f"XYZ atom-count mismatch: {path}"
        )

    atoms = []

    for index, line in enumerate(atom_lines):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Malformed XYZ line {index}: {line}"
            )

        atoms.append({
            "index_0based": index,
            "element": fields[0],
            "x_A": float(fields[1]),
            "y_A": float(fields[2]),
            "z_A": float(fields[3]),
        })

    return atoms, lines[1]


def xyz_distance(
    first: dict,
    second: dict,
) -> float:
    return math.sqrt(
        (first["x_A"] - second["x_A"]) ** 2
        + (first["y_A"] - second["y_A"]) ** 2
        + (first["z_A"] - second["z_A"]) ** 2
    )


def canonical_pair(
    first: str,
    second: str,
) -> tuple[str, str]:
    return tuple(sorted((first, second)))


def identify_edge_columns(
    rows: list[dict[str, str]],
) -> tuple[str, str]:
    if not rows:
        raise RuntimeError(
            "Nominal edge table is empty."
        )

    headers = set(rows[0])

    candidates = [
        ("first_atom", "second_atom"),
        ("source_atom", "target_atom"),
        ("source_node", "target_node"),
    ]

    for first_field, second_field in candidates:
        if {
            first_field,
            second_field,
        } <= headers:
            return first_field, second_field

    raise RuntimeError(
        "Could not identify edge columns: "
        f"{sorted(headers)}"
    )


def main() -> None:
    for path in (
        V6B_XYZ,
        V6B_MAP,
        V6A_EDGES,
        V7_REPAIRED_XYZ,
        V7_LOCAL,
        V7_BOUNDARY,
        V7_H0045_REPORT,
    ):
        require_file(path)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    v6b_atoms, _ = read_xyz(V6B_XYZ)
    v7_atoms, v7_comment = read_xyz(
        V7_REPAIRED_XYZ
    )

    v6b_map_rows = read_csv(V6B_MAP)
    local_rows = read_csv(V7_LOCAL)
    source_edge_rows = read_csv(V6A_EDGES)

    if len(v6b_atoms) != 48:
        raise RuntimeError(
            f"Expected 48 V6-B atoms; "
            f"found {len(v6b_atoms)}"
        )

    if len(v6b_map_rows) != 48:
        raise RuntimeError(
            f"Expected 48 V6-B map rows; "
            f"found {len(v6b_map_rows)}"
        )

    if len(v7_atoms) != 52:
        raise RuntimeError(
            f"Expected 52 V7-A atoms; "
            f"found {len(v7_atoms)}"
        )

    if len(local_rows) != 5:
        raise RuntimeError(
            f"Expected five repaired local records; "
            f"found {len(local_rows)}"
        )

    v6b_map_rows.sort(
        key=lambda row: int(
            row["v6b_index_0based"]
        )
    )

    for expected_index, row in enumerate(
        v6b_map_rows
    ):
        actual_index = int(
            row["v6b_index_0based"]
        )

        if actual_index != expected_index:
            raise RuntimeError(
                "V6-B map index sequence is not "
                f"contiguous at {expected_index}: "
                f"{actual_index}"
            )

        if (
            v6b_atoms[expected_index]["element"]
            != row["element"]
        ):
            raise RuntimeError(
                "V6-B map/XYZ element mismatch at "
                f"index {expected_index}: "
                f"{row['atom_id']}"
            )

    obsolete_rows = [
        row
        for row in v6b_map_rows
        if row["atom_id"] == OBSOLETE_ATOM
    ]

    if len(obsolete_rows) != 1:
        raise RuntimeError(
            "Expected exactly one obsolete HCAPV2 "
            f"record; found {len(obsolete_rows)}"
        )

    obsolete_index = int(
        obsolete_rows[0]["v6b_index_0based"]
    )

    retained_v6b_rows = [
        row
        for row in v6b_map_rows
        if row["atom_id"] != OBSOLETE_ATOM
    ]

    retained_v6b_atoms = [
        atom
        for index, atom in enumerate(v6b_atoms)
        if index != obsolete_index
    ]

    if len(retained_v6b_rows) != 47:
        raise RuntimeError(
            "Expected 47 retained V6-B records."
        )

    if len(retained_v6b_atoms) != 47:
        raise RuntimeError(
            "Expected 47 retained V6-B coordinates."
        )

    local_by_id = {
        row["atom_id"]: row
        for row in local_rows
    }

    if set(local_by_id) != set(NEW_ATOM_ORDER):
        raise RuntimeError(
            "Unexpected repaired local identity set. "
            f"Observed: {sorted(local_by_id)}"
        )

    expected_atom_ids = [
        row["atom_id"]
        for row in retained_v6b_rows
    ] + NEW_ATOM_ORDER

    if len(expected_atom_ids) != 52:
        raise RuntimeError(
            "Formal identity list is not 52 atoms."
        )

    # The provisional V7 XYZ must preserve all retained V6-B
    # coordinates in order, followed by the five repaired atoms.
    maximum_retained_difference = 0.0

    for index in range(47):
        source_atom = retained_v6b_atoms[index]
        candidate_atom = v7_atoms[index]

        if (
            source_atom["element"]
            != candidate_atom["element"]
        ):
            raise RuntimeError(
                "Retained V6-B/V7 element mismatch "
                f"at sequential index {index}"
            )

        difference = xyz_distance(
            source_atom,
            candidate_atom,
        )

        maximum_retained_difference = max(
            maximum_retained_difference,
            difference,
        )

    if maximum_retained_difference > 5.0e-7:
        raise RuntimeError(
            "Retained V6-B coordinates were changed "
            "during provisional V7 construction. "
            f"Maximum difference: "
            f"{maximum_retained_difference}"
        )

    maximum_local_difference = 0.0

    for offset, atom_id in enumerate(
        NEW_ATOM_ORDER
    ):
        xyz_index = 47 + offset
        xyz_atom = v7_atoms[xyz_index]
        local_row = local_by_id[atom_id]

        if xyz_atom["element"] != local_row["element"]:
            raise RuntimeError(
                f"Local element mismatch for {atom_id}"
            )

        local_atom = {
            "x_A": float(local_row["x_A"]),
            "y_A": float(local_row["y_A"]),
            "z_A": float(local_row["z_A"]),
        }

        difference = xyz_distance(
            xyz_atom,
            local_atom,
        )

        maximum_local_difference = max(
            maximum_local_difference,
            difference,
        )

    if maximum_local_difference > 5.0e-10:
        raise RuntimeError(
            "Repaired local coordinates do not match "
            "the final five XYZ records. "
            f"Maximum difference: "
            f"{maximum_local_difference}"
        )

    composition = Counter(
        atom["element"]
        for atom in v7_atoms
    )

    if dict(composition) != EXPECTED_COMPOSITION:
        raise RuntimeError(
            "Unexpected formal V7-A composition: "
            f"{dict(composition)}"
        )

    # Build formal provenance map.
    formal_map_rows = []

    inherited_fixed_count = 0

    for new_index, row in enumerate(
        retained_v6b_rows
    ):
        source_index = int(
            row["v6b_index_0based"]
        )

        fixed = (
            row["v6b_fixed"]
            .strip()
            .lower()
            == "true"
        )

        inherited_fixed_count += int(fixed)

        formal_map_rows.append({
            "v7a_index_0based": new_index,
            "atom_id": row["atom_id"],
            "element": row["element"],
            "atom_role": row["atom_role"],
            "node_type": row["node_type"],
            "coordinate_source": (
                "RETAINED_QM_F06_UPPER_V6B_START"
            ),
            "source_v6b_index_0based": source_index,
            "source_center": row.get(
                "source_center",
                "",
            ),
            "source_outside_atom": row.get(
                "source_outside_atom",
                "",
            ),
            "source_cut_edge": row.get(
                "source_cut_edge",
                "",
            ),
            "artificial_cap": row.get(
                "artificial_cap",
                "False",
            ),
            "v7a_topology_action": (
                "RETAINED_FROM_V6B"
            ),
            "v7a_fixed": str(fixed),
            "v7a_mobile": str(not fixed),
            "v7a_constraint_basis": (
                "INHERITED_V6B_FIXED_CORE"
                if fixed
                else "INHERITED_V6B_MOBILE_REGION"
            ),
        })

    new_metadata = {
        "A:UPPER:13:1": {
            "atom_role": "REAL_V7_BOUNDARY_EXPANSION_ATOM",
            "node_type": "ANNULUS_INTERIOR",
            "coordinate_source": (
                "V7A_AXIAL_CANONICAL_EMBEDDING"
            ),
            "source_center": "",
            "source_outside_atom": "",
            "source_cut_edge": "",
            "artificial_cap": "False",
            "v7a_topology_action": (
                "ADDED_REAL_CANONICAL_ATOM"
            ),
        },
        "A:UPPER:14:0": {
            "atom_role": "REAL_V7_BOUNDARY_EXPANSION_ATOM",
            "node_type": "ANNULUS_OUTER_BOUNDARY",
            "coordinate_source": (
                "V7A_AXIAL_CANONICAL_EMBEDDING"
            ),
            "source_center": "",
            "source_outside_atom": "",
            "source_cut_edge": "",
            "artificial_cap": "False",
            "v7a_topology_action": (
                "ADDED_REAL_CANONICAL_ATOM"
            ),
        },
        "HCAPV7:UPPER:A13_1:A11_1": {
            "atom_role": "ARTIFICIAL_BOUNDARY_CAP_V7",
            "node_type": "QM_BOUNDARY_CAP_H",
            "coordinate_source": (
                "V7A_REPAIRED_LOCAL_COORDINATES"
            ),
            "source_center": "A:UPPER:13:1",
            "source_outside_atom": "A:UPPER:11:1",
            "source_cut_edge": (
                "A:UPPER:13:1--A:UPPER:11:1"
            ),
            "artificial_cap": "True",
            "v7a_topology_action": (
                "ADDED_V7_BOUNDARY_CAP"
            ),
        },
        "HCAPV7:UPPER:A14_0:A13_M1": {
            "atom_role": "ARTIFICIAL_BOUNDARY_CAP_V7",
            "node_type": "QM_BOUNDARY_CAP_H",
            "coordinate_source": (
                "V7A_REPAIRED_LOCAL_COORDINATES"
            ),
            "source_center": "A:UPPER:14:0",
            "source_outside_atom": "A:UPPER:13:-1",
            "source_cut_edge": (
                "A:UPPER:14:0--A:UPPER:13:-1"
            ),
            "artificial_cap": "True",
            "v7a_topology_action": (
                "ADDED_V7_BOUNDARY_CAP"
            ),
        },
        "H4:UPPER:0045:0": {
            "atom_role": "CANONICAL_R2_HYDROGEN_ADDED",
            "node_type": "ANNULUS_OUTER_PASSIVANT_H",
            "coordinate_source": (
                "V7A_TRIGONAL_COMPLETION_"
                "FROM_A13_1_AND_A13_M1_CUT_DIRECTION"
            ),
            "source_center": "A:UPPER:14:0",
            "source_outside_atom": "",
            "source_cut_edge": "",
            "artificial_cap": "False",
            "v7a_topology_action": (
                "RESTORED_CANONICAL_PASSIVANT"
            ),
        },
    }

    for offset, atom_id in enumerate(
        NEW_ATOM_ORDER
    ):
        index = 47 + offset
        local_row = local_by_id[atom_id]
        metadata = new_metadata[atom_id]

        formal_map_rows.append({
            "v7a_index_0based": index,
            "atom_id": atom_id,
            "element": local_row["element"],
            "atom_role": metadata["atom_role"],
            "node_type": metadata["node_type"],
            "coordinate_source": metadata[
                "coordinate_source"
            ],
            "source_v6b_index_0based": "",
            "source_center": metadata[
                "source_center"
            ],
            "source_outside_atom": metadata[
                "source_outside_atom"
            ],
            "source_cut_edge": metadata[
                "source_cut_edge"
            ],
            "artificial_cap": metadata[
                "artificial_cap"
            ],
            "v7a_topology_action": metadata[
                "v7a_topology_action"
            ],
            "v7a_fixed": "False",
            "v7a_mobile": "True",
            "v7a_constraint_basis": (
                "NEW_OR_REPAIRED_V7_MOBILE_REGION"
            ),
        })

    if len(formal_map_rows) != 52:
        raise RuntimeError(
            "Formal V7-A map is not 52 rows."
        )

    for index, row in enumerate(
        formal_map_rows
    ):
        if int(
            row["v7a_index_0based"]
        ) != index:
            raise RuntimeError(
                "Formal map order failure."
            )

        if row["atom_id"] != expected_atom_ids[index]:
            raise RuntimeError(
                "Formal identity order failure at "
                f"index {index}"
            )

        if row["element"] != v7_atoms[index]["element"]:
            raise RuntimeError(
                "Formal map/XYZ element failure at "
                f"index {index}"
            )

    # Build formal nominal edge inventory.
    first_field, second_field = (
        identify_edge_columns(source_edge_rows)
    )

    inherited_edges = set()

    for row in source_edge_rows:
        first = row[first_field]
        second = row[second_field]

        if OBSOLETE_ATOM in (first, second):
            continue

        if (
            first not in expected_atom_ids
            or second not in expected_atom_ids
        ):
            continue

        inherited_edges.add(
            canonical_pair(first, second)
        )

    new_edge_types = {}

    for first, second, edge_type in NEW_EDGES:
        pair = canonical_pair(first, second)
        inherited_edges.add(pair)
        new_edge_types[pair] = edge_type

    formal_edges = sorted(inherited_edges)

    if len(formal_edges) != 57:
        raise RuntimeError(
            "Expected 57 formal V7-A nominal edges; "
            f"found {len(formal_edges)}"
        )

    degree = defaultdict(int)

    for first, second in formal_edges:
        degree[first] += 1
        degree[second] += 1

    expected_degree = {
        "B": 3,
        "N": 3,
        "H": 1,
    }

    degree_failures = []

    element_by_id = {
        row["atom_id"]: row["element"]
        for row in formal_map_rows
    }

    for atom_id in expected_atom_ids:
        observed = degree[atom_id]
        expected = expected_degree[
            element_by_id[atom_id]
        ]

        if observed != expected:
            degree_failures.append({
                "atom_id": atom_id,
                "element": element_by_id[atom_id],
                "observed_degree": observed,
                "expected_degree": expected,
            })

    if degree_failures:
        raise RuntimeError(
            "Formal V7-A degree failures: "
            + json.dumps(
                degree_failures,
                indent=2,
            )
        )

    # Verify connected graph.
    adjacency = defaultdict(set)

    for first, second in formal_edges:
        adjacency[first].add(second)
        adjacency[second].add(first)

    visited = set()
    stack = [expected_atom_ids[0]]

    while stack:
        current = stack.pop()

        if current in visited:
            continue

        visited.add(current)

        stack.extend(
            adjacency[current] - visited
        )

    if visited != set(expected_atom_ids):
        raise RuntimeError(
            "Formal V7-A graph is disconnected."
        )

    # Write formal XYZ without altering coordinates.
    xyz_lines = [
        "52",
        (
            "QM_F06 UPPER V7-A formal construction; "
            "HCAPV2:UPPER:03 removed; canonical "
            "A13:1/A14:0 expansion and H0045 restored"
        ),
    ]

    for atom in v7_atoms:
        xyz_lines.append(
            f"{atom['element']:2s} "
            f"{atom['x_A']: .12f} "
            f"{atom['y_A']: .12f} "
            f"{atom['z_A']: .12f}"
        )

    OUTPUT_XYZ.write_text(
        "\n".join(xyz_lines) + "\n",
        encoding="utf-8",
    )

    map_fieldnames = list(
        formal_map_rows[0].keys()
    )

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
        writer.writerows(formal_map_rows)

    coordinate_by_id = {
        expected_atom_ids[index]: atom
        for index, atom in enumerate(v7_atoms)
    }

    edge_records = []

    for edge_index, (first, second) in enumerate(
        formal_edges,
        start=1,
    ):
        pair = canonical_pair(first, second)

        if pair in new_edge_types:
            provenance = "ADDED_BY_V7A"
            edge_type = new_edge_types[pair]
        else:
            provenance = "INHERITED_FROM_V6A"
            edge_type = "INHERITED_NOMINAL_EDGE"

        edge_records.append({
            "edge_index": edge_index,
            "first_atom": first,
            "first_element": element_by_id[first],
            "second_atom": second,
            "second_element": element_by_id[second],
            "edge_type": edge_type,
            "provenance": provenance,
            "distance_A": xyz_distance(
                coordinate_by_id[first],
                coordinate_by_id[second],
            ),
        })

    with OUTPUT_EDGES.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                edge_records[0].keys()
            ),
        )
        writer.writeheader()
        writer.writerows(edge_records)

    cap_records = []

    for cap_id in (
        "HCAPV7:UPPER:A13_1:A11_1",
        "HCAPV7:UPPER:A14_0:A13_M1",
    ):
        metadata = new_metadata[cap_id]
        owner = metadata["source_center"]

        cap_records.append({
            "cap_id": cap_id,
            "owner_atom": owner,
            "owner_element": element_by_id[owner],
            "outside_atom": metadata[
                "source_outside_atom"
            ],
            "source_cut_edge": metadata[
                "source_cut_edge"
            ],
            "owner_distance_A": xyz_distance(
                coordinate_by_id[cap_id],
                coordinate_by_id[owner],
            ),
            "cap_x_A": coordinate_by_id[
                cap_id
            ]["x_A"],
            "cap_y_A": coordinate_by_id[
                cap_id
            ]["y_A"],
            "cap_z_A": coordinate_by_id[
                cap_id
            ]["z_A"],
            "coordinate_source": (
                "V7A_H0045_REPAIRED_LOCAL_COORDINATES"
            ),
        })

    with OUTPUT_CAPS.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                cap_records[0].keys()
            ),
        )
        writer.writeheader()
        writer.writerows(cap_records)

    report = {
        "model": "QM_F06_UPPER_V7A",
        "decision": (
            "QM_F06_UPPER_V7A_FORMALLY_CONSTRUCTED_"
            "GLOBAL_PRE_QM_AUDIT_REQUIRED"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "atom_count": len(v7_atoms),
        "composition": dict(
            sorted(composition.items())
        ),
        "nominal_edge_count": len(
            formal_edges
        ),
        "connected_component_count": 1,
        "degree_failure_count": len(
            degree_failures
        ),
        "topology_changes_from_v6b": {
            "removed_atom": OBSOLETE_ATOM,
            "added_atoms": NEW_ATOM_ORDER,
            "added_edges": [
                {
                    "first_atom": first,
                    "second_atom": second,
                    "edge_type": edge_type,
                }
                for first, second, edge_type
                in NEW_EDGES
            ],
        },
        "coordinate_validation": {
            "maximum_retained_v6b_difference_A": (
                maximum_retained_difference
            ),
            "maximum_repaired_local_difference_A": (
                maximum_local_difference
            ),
            "repaired_local_coordinates_used": True,
            "obsolete_pre_repair_cap_coordinates_used": False,
        },
        "constraint_design": {
            "inherited_fixed_atom_count": (
                inherited_fixed_count
            ),
            "new_atoms_mobile": True,
            "formal_constraint_design_authorized": False,
        },
        "sources": {
            "v6b_xyz": str(
                V6B_XYZ.relative_to(ROOT)
            ),
            "v6b_map": str(
                V6B_MAP.relative_to(ROOT)
            ),
            "v6a_edges": str(
                V6A_EDGES.relative_to(ROOT)
            ),
            "v7_repaired_xyz": str(
                V7_REPAIRED_XYZ.relative_to(ROOT)
            ),
            "v7_repaired_local_coordinates": str(
                V7_LOCAL.relative_to(ROOT)
            ),
            "v7_boundary_provenance": str(
                V7_BOUNDARY.relative_to(ROOT)
            ),
            "v7_h0045_repair_report": str(
                V7_H0045_REPORT.relative_to(ROOT)
            ),
        },
        "sha256": {
            "formal_xyz": file_sha256(
                OUTPUT_XYZ
            ),
            "formal_map": file_sha256(
                OUTPUT_MAP
            ),
            "formal_edges": file_sha256(
                OUTPUT_EDGES
            ),
            "formal_caps": file_sha256(
                OUTPUT_CAPS
            ),
        },
        "global_pre_qm_audit_authorized": True,
        "ORCA_input_design_authorized": False,
        "ORCA_execution_authorized": False,
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
    print("QM_F06 UPPER V7-A FORMAL CONSTRUCTION")
    print("=" * 104)
    print("Atoms:", len(v7_atoms))
    print(
        "Composition:",
        dict(sorted(composition.items())),
    )
    print(
        "Nominal edges:",
        len(formal_edges),
    )
    print(
        "Connected components:",
        1,
    )
    print(
        "Degree failures:",
        len(degree_failures),
    )
    print(
        "Inherited fixed atoms:",
        inherited_fixed_count,
    )
    print(
        "Maximum retained-coordinate difference A:",
        maximum_retained_difference,
    )
    print(
        "Maximum repaired-local difference A:",
        maximum_local_difference,
    )

    print()
    print("Removed atom:")
    print("  ", OBSOLETE_ATOM)

    print()
    print("Added atoms:")

    for atom_id in NEW_ATOM_ORDER:
        print("  ", atom_id)

    print()
    print("Added nominal edges:")

    for first, second, edge_type in NEW_EDGES:
        value = xyz_distance(
            coordinate_by_id[first],
            coordinate_by_id[second],
        )

        print(
            f"  {first:38s} -- "
            f"{second:38s} "
            f"{value:.6f} Å | {edge_type}"
        )

    print()
    print("Decision:", report["decision"])
    print("XYZ:", OUTPUT_XYZ)
    print("Map:", OUTPUT_MAP)
    print("Edges:", OUTPUT_EDGES)
    print("Caps:", OUTPUT_CAPS)
    print("Report:", OUTPUT_REPORT)
    print()
    print("Global pre-QM audit authorized: True")
    print("ORCA input design authorized: False")
    print("ORCA execution authorized: False")
    print("RESP authorized: False")
    print("MD authorized: False")


if __name__ == "__main__":
    main()
