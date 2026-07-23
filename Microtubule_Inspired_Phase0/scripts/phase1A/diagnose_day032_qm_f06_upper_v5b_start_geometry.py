#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

START_XYZ = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_construction/"
    "QM_F06_UPPER_V5B_start.xyz"
)

MAP_CSV = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_construction/"
    "QM_F06_UPPER_V5B_atom_role_provenance_map.csv"
)

VALENCE_CSV = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_pre_qm_audit/"
    "QM_F06_UPPER_V5B_valence.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day032_qm_f06_upper_v5b_start_geometry_diagnostic"
)

OUTPUT_CSV = OUTPUT_DIR / (
    "QM_F06_UPPER_V5B_start_geometric_reconnectivity.csv"
)

OUTPUT_JSON = OUTPUT_DIR / (
    "QM_F06_UPPER_V5B_START_GEOMETRY_DIAGNOSTIC.json"
)


def read_csv(path: Path):
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def read_xyz(path: Path):
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    count = int(lines[0])
    atoms = []

    for line in lines[2:2 + count]:
        fields = line.split()

        atoms.append({
            "element": fields[0],
            "xyz": tuple(map(float, fields[1:4])),
        })

    if len(atoms) != count:
        raise RuntimeError(
            f"Incomplete XYZ: expected {count}, "
            f"found {len(atoms)}"
        )

    return atoms


def distance(first, second):
    return math.sqrt(
        sum(
            (a - b) ** 2
            for a, b in zip(first, second)
        )
    )


def geometric_bond(
    first_element: str,
    second_element: str,
    value: float,
) -> bool:
    pair = frozenset((
        first_element,
        second_element,
    ))

    if pair == frozenset(("B", "N")):
        return 1.25 <= value <= 1.85

    if pair == frozenset(("B", "H")):
        return 0.95 <= value <= 1.35

    if pair == frozenset(("N", "H")):
        return 0.85 <= value <= 1.25

    return False


def canonical_pair(first: str, second: str):
    return tuple(sorted((first, second)))


def main():
    mapping = read_csv(MAP_CSV)
    valence = read_csv(VALENCE_CSV)
    atoms = read_xyz(START_XYZ)

    if len(mapping) != len(atoms):
        raise RuntimeError(
            "Map/XYZ atom-count mismatch."
        )

    id_by_index = {
        int(row["index_0based"]): row["atom_id"]
        for row in mapping
    }

    index_by_id = {
        row["atom_id"]: int(row["index_0based"])
        for row in mapping
    }

    nominal_neighbors = {
        row["atom_id"]: set(
            value
            for value in row["neighbors"].split("|")
            if value
        )
        for row in valence
    }

    nominal_pairs = set()

    for atom_id, neighbors in nominal_neighbors.items():
        for neighbor in neighbors:
            nominal_pairs.add(
                canonical_pair(atom_id, neighbor)
            )

    geometric_pairs = set()
    geometric_distances = {}

    for first_index in range(len(atoms)):
        for second_index in range(
            first_index + 1,
            len(atoms),
        ):
            first = atoms[first_index]
            second = atoms[second_index]

            value = distance(
                first["xyz"],
                second["xyz"],
            )

            if geometric_bond(
                first["element"],
                second["element"],
                value,
            ):
                pair = canonical_pair(
                    id_by_index[first_index],
                    id_by_index[second_index],
                )

                geometric_pairs.add(pair)
                geometric_distances[pair] = value

    gained_pairs = sorted(
        geometric_pairs - nominal_pairs
    )

    lost_pairs = sorted(
        nominal_pairs - geometric_pairs
    )

    records = []

    for first, second in gained_pairs:
        first_index = index_by_id[first]
        second_index = index_by_id[second]

        records.append({
            "classification": "NONNOMINAL_GEOMETRIC_BOND",
            "first_atom": first,
            "first_element": atoms[first_index]["element"],
            "second_atom": second,
            "second_element": atoms[second_index]["element"],
            "distance_A": geometric_distances[
                canonical_pair(first, second)
            ],
        })

    for first, second in lost_pairs:
        first_index = index_by_id[first]
        second_index = index_by_id[second]

        records.append({
            "classification": "NOMINAL_BOND_OUTSIDE_RANGE",
            "first_atom": first,
            "first_element": atoms[first_index]["element"],
            "second_atom": second,
            "second_element": atoms[second_index]["element"],
            "distance_A": distance(
                atoms[first_index]["xyz"],
                atoms[second_index]["xyz"],
            ),
        })

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    with OUTPUT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "classification",
                "first_atom",
                "first_element",
                "second_atom",
                "second_element",
                "distance_A",
            ],
        )
        writer.writeheader()
        writer.writerows(records)

    affected_atoms = sorted({
        atom_id
        for pair in gained_pairs + lost_pairs
        for atom_id in pair
    })

    report = {
        "decision": (
            "QM_F06_UPPER_V5B_START_GEOMETRY_"
            "CONFLICTS_DETECTED"
            if records
            else
            "QM_F06_UPPER_V5B_START_GEOMETRY_"
            "NOMINAL_GRAPH_CONSISTENT"
        ),
        "atom_count": len(atoms),
        "nonnominal_geometric_bond_count": len(
            gained_pairs
        ),
        "nominal_bond_outside_range_count": len(
            lost_pairs
        ),
        "affected_atom_count": len(affected_atoms),
        "affected_atoms": affected_atoms,
        "orca_authorized": False,
        "RESP_authorized": False,
    }

    OUTPUT_JSON.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 96)
    print("QM_F06 UPPER V5-B START-GEOMETRY DIAGNOSTIC")
    print("=" * 96)
    print(
        "Nonnominal geometric bonds:",
        len(gained_pairs),
    )
    print(
        "Nominal bonds outside range:",
        len(lost_pairs),
    )

    print()
    print("Nonnominal geometric bonds:")

    for first, second in gained_pairs:
        pair = canonical_pair(first, second)

        print(
            f"  {first:28s} -- {second:28s} "
            f"{geometric_distances[pair]:.6f} Å"
        )

    print()
    print("Nominal bonds outside geometric range:")

    for first, second in lost_pairs:
        value = distance(
            atoms[index_by_id[first]]["xyz"],
            atoms[index_by_id[second]]["xyz"],
        )

        print(
            f"  {first:28s} -- {second:28s} "
            f"{value:.6f} Å"
        )

    print()
    print("Decision:", report["decision"])
    print("CSV:", OUTPUT_CSV)
    print("Report:", OUTPUT_JSON)


if __name__ == "__main__":
    main()
