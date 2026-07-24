#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

SOURCE_XYZ = ROOT / (
    "runs/phase1A/"
    "day032_qm_f06_upper_v5c_construction/"
    "QM_F06_UPPER_V5C_start.xyz"
)

SOURCE_MAP = ROOT / (
    "runs/phase1A/"
    "day032_qm_f06_upper_v5c_construction/"
    "QM_F06_UPPER_V5C_atom_role_provenance_map.csv"
)

SOURCE_VALENCE = ROOT / (
    "runs/phase1A/"
    "day031_qm_f06_upper_v5b_pre_qm_audit/"
    "QM_F06_UPPER_V5B_valence.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day033_qm_f06_upper_v6a_topology_closure"
)

OUTPUT_XYZ = OUTPUT_DIR / (
    "QM_F06_UPPER_V6A_TOPOLOGY_CLOSURE_start.xyz"
)

OUTPUT_MAP = OUTPUT_DIR / (
    "QM_F06_UPPER_V6A_atom_role_provenance_map.csv"
)

OUTPUT_EDGES = OUTPUT_DIR / (
    "QM_F06_UPPER_V6A_nominal_edges.csv"
)

OUTPUT_REPORT = OUTPUT_DIR / (
    "QM_F06_UPPER_V6A_TOPOLOGY_CLOSURE_REPORT.json"
)


REMOVED_HYDROGENS = {
    "H4:UPPER:0116:0",
    "H4:UPPER:0203:0",
    "H4:UPPER:0117:0",
    "H4:UPPER:0170:0",
}

REMOVED_BONDS = {
    tuple(sorted((
        "BR4:UPPER:00:3",
        "H4:UPPER:0116:0",
    ))),
    tuple(sorted((
        "S:1739",
        "H4:UPPER:0203:0",
    ))),
    tuple(sorted((
        "BR4:UPPER:00:4",
        "H4:UPPER:0117:0",
    ))),
    tuple(sorted((
        "BR4:UPPER:14:1",
        "H4:UPPER:0170:0",
    ))),
}

ADDED_CLOSURE_BONDS = {
    tuple(sorted((
        "BR4:UPPER:00:3",
        "S:1739",
    ))),
    tuple(sorted((
        "BR4:UPPER:00:4",
        "BR4:UPPER:14:1",
    ))),
}


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

    atom_count = int(lines[0].strip())

    coordinate_lines = lines[2:2 + atom_count]

    if len(coordinate_lines) != atom_count:
        raise RuntimeError(
            "Incomplete XYZ coordinate block."
        )

    atoms = []

    for index, line in enumerate(coordinate_lines):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Malformed XYZ line {index + 3}: {line}"
            )

        atoms.append({
            "index_0based": index,
            "element": fields[0],
            "xyz_A": (
                float(fields[1]),
                float(fields[2]),
                float(fields[3]),
            ),
        })

    return atoms


def canonical_pair(
    first: str,
    second: str,
) -> tuple[str, str]:
    return tuple(sorted((first, second)))


def main() -> None:
    for path in (
        SOURCE_XYZ,
        SOURCE_MAP,
        SOURCE_VALENCE,
    ):
        require_file(path)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    atoms = read_xyz(SOURCE_XYZ)
    mapping = read_csv(SOURCE_MAP)
    valence_rows = read_csv(SOURCE_VALENCE)

    if len(atoms) != len(mapping):
        raise RuntimeError(
            "XYZ/map atom-count mismatch."
        )

    source_records = []

    for atom, row in zip(atoms, mapping):
        if atom["element"] != row["element"]:
            raise RuntimeError(
                "XYZ/map element mismatch at index "
                f"{atom['index_0based']}."
            )

        source_records.append({
            **row,
            "xyz_A": atom["xyz_A"],
        })

    source_ids = {
        row["atom_id"]
        for row in source_records
    }

    missing_removed = (
        REMOVED_HYDROGENS - source_ids
    )

    if missing_removed:
        raise RuntimeError(
            "Hydrogens requested for removal are missing: "
            + "|".join(sorted(missing_removed))
        )

    retained = [
        row
        for row in source_records
        if row["atom_id"] not in REMOVED_HYDROGENS
    ]

    old_to_new_index = {}

    for new_index, row in enumerate(retained):
        old_to_new_index[
            int(row["index_0based"])
        ] = new_index

    edges = set()

    for row in valence_rows:
        center = row["atom_id"]

        if center not in source_ids:
            continue

        neighbors = [
            value
            for value in row["neighbors"].split("|")
            if value
        ]

        for neighbor in neighbors:
            if neighbor not in source_ids:
                continue

            edges.add(
                canonical_pair(center, neighbor)
            )

    original_edges = set(edges)

    edges -= REMOVED_BONDS

    edges = {
        pair
        for pair in edges
        if not (
            pair[0] in REMOVED_HYDROGENS
            or pair[1] in REMOVED_HYDROGENS
        )
    }

    edges |= ADDED_CLOSURE_BONDS

    retained_ids = {
        row["atom_id"]
        for row in retained
    }

    invalid_edges = [
        pair
        for pair in edges
        if (
            pair[0] not in retained_ids
            or pair[1] not in retained_ids
        )
    ]

    if invalid_edges:
        raise RuntimeError(
            "Edges reference removed or absent atoms: "
            + repr(invalid_edges)
        )

    adjacency = {
        atom_id: set()
        for atom_id in retained_ids
    }

    for first, second in edges:
        adjacency[first].add(second)
        adjacency[second].add(first)

    expected_degree = {
        "B": 3,
        "N": 3,
        "H": 1,
    }

    degree_failures = []

    for row in retained:
        atom_id = row["atom_id"]
        element = row["element"]

        degree = len(adjacency[atom_id])
        expected = expected_degree[element]

        if degree != expected:
            degree_failures.append({
                "atom_id": atom_id,
                "element": element,
                "degree": degree,
                "expected_degree": expected,
                "neighbors": sorted(
                    adjacency[atom_id]
                ),
            })

    if degree_failures:
        raise RuntimeError(
            "Topology-closure degree failures: "
            + json.dumps(
                degree_failures,
                indent=2,
            )
        )

    composition = Counter(
        row["element"]
        for row in retained
    )

    expected_composition = Counter({
        "B": 16,
        "N": 13,
        "H": 19,
    })

    if composition != expected_composition:
        raise RuntimeError(
            "Unexpected V6-A composition: "
            f"{dict(composition)}"
        )

    with OUTPUT_XYZ.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(f"{len(retained)}\n")
        handle.write(
            "QM_F06 UPPER V6-A topology-closure "
            "candidate; ORCA not yet authorized\n"
        )

        for row in retained:
            x_value, y_value, z_value = row["xyz_A"]

            handle.write(
                f"{row['element']:2s} "
                f"{x_value: .12f} "
                f"{y_value: .12f} "
                f"{z_value: .12f}\n"
            )

    map_fieldnames = list(mapping[0].keys()) + [
        "v6a_index_0based",
        "v6a_retained",
        "v6a_topology_action",
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

        new_index_by_id = {
            row["atom_id"]: index
            for index, row in enumerate(retained)
        }

        for row in mapping:
            atom_id = row["atom_id"]
            retained_flag = (
                atom_id not in REMOVED_HYDROGENS
            )

            record = dict(row)
            record["v6a_index_0based"] = (
                new_index_by_id[atom_id]
                if retained_flag
                else ""
            )
            record["v6a_retained"] = retained_flag

            if atom_id in REMOVED_HYDROGENS:
                record["v6a_topology_action"] = (
                    "REMOVED_PASSIVANT_FOR_BN_CLOSURE"
                )
            elif any(
                atom_id in pair
                for pair in ADDED_CLOSURE_BONDS
            ):
                record["v6a_topology_action"] = (
                    "RETAINED_BN_CLOSURE_CENTER"
                )
            else:
                record["v6a_topology_action"] = (
                    "RETAINED_UNCHANGED"
                )

            writer.writerow(record)

    edge_fieldnames = [
        "first_atom",
        "second_atom",
        "first_element",
        "second_element",
        "edge_source",
        "new_in_v6a",
    ]

    element_by_id = {
        row["atom_id"]: row["element"]
        for row in retained
    }

    with OUTPUT_EDGES.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=edge_fieldnames,
        )

        writer.writeheader()

        for first, second in sorted(edges):
            is_new = (
                canonical_pair(first, second)
                in ADDED_CLOSURE_BONDS
            )

            writer.writerow({
                "first_atom": first,
                "second_atom": second,
                "first_element": element_by_id[first],
                "second_element": element_by_id[second],
                "edge_source": (
                    "V6A_TOPOLOGY_CLOSURE"
                    if is_new
                    else "V5C_RETAINED_NOMINAL"
                ),
                "new_in_v6a": is_new,
            })

    report = {
        "decision": (
            "QM_F06_UPPER_V6A_TOPOLOGY_CLOSURE_"
            "CONSTRUCTED_STRUCTURAL_AUDIT_REQUIRED"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "source_model": "QM_F06_UPPER_V5C",
        "atom_count": len(retained),
        "composition": dict(
            sorted(composition.items())
        ),
        "removed_hydrogens": sorted(
            REMOVED_HYDROGENS
        ),
        "removed_bonds": [
            list(pair)
            for pair in sorted(REMOVED_BONDS)
        ],
        "added_closure_bonds": [
            list(pair)
            for pair in sorted(
                ADDED_CLOSURE_BONDS
            )
        ],
        "original_nominal_edge_count": len(
            original_edges
        ),
        "v6a_nominal_edge_count": len(edges),
        "degree_failure_count": len(
            degree_failures
        ),
        "files": {
            "start_xyz": str(
                OUTPUT_XYZ.relative_to(ROOT)
            ),
            "provenance_map": str(
                OUTPUT_MAP.relative_to(ROOT)
            ),
            "nominal_edges": str(
                OUTPUT_EDGES.relative_to(ROOT)
            ),
        },
        "sha256": {
            "start_xyz": sha256(OUTPUT_XYZ),
            "provenance_map": sha256(
                OUTPUT_MAP
            ),
            "nominal_edges": sha256(
                OUTPUT_EDGES
            ),
        },
        "structural_audit_authorized": True,
        "orca_input_design_authorized": False,
        "orca_execution_authorized": False,
        "RESP_authorized": False,
        "force_field_adoption_authorized": False,
        "MD_authorized": False,
    }

    OUTPUT_REPORT.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 96)
    print("QM_F06 UPPER V6-A TOPOLOGY CLOSURE")
    print("=" * 96)
    print("Source atoms:", len(source_records))
    print("Retained atoms:", len(retained))
    print(
        "Composition:",
        dict(sorted(composition.items())),
    )
    print(
        "Original nominal edges:",
        len(original_edges),
    )
    print("V6-A nominal edges:", len(edges))
    print(
        "Degree failures:",
        len(degree_failures),
    )

    print()
    print("Removed hydrogens:")

    for atom_id in sorted(
        REMOVED_HYDROGENS
    ):
        print("  ", atom_id)

    print()
    print("Added B-N closure bonds:")

    for first, second in sorted(
        ADDED_CLOSURE_BONDS
    ):
        print(f"  {first} -- {second}")

    print()
    print("Decision:", report["decision"])
    print("XYZ:", OUTPUT_XYZ)
    print("Edges:", OUTPUT_EDGES)
    print("Map:", OUTPUT_MAP)
    print("Report:", OUTPUT_REPORT)
    print()
    print("Structural audit authorized: True")
    print("ORCA execution authorized: False")
    print("RESP authorized: False")
    print("MD authorized: False")


if __name__ == "__main__":
    main()
