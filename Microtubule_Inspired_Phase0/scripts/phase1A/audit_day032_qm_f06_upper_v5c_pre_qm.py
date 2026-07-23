#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

CONSTRUCTION_DIR = ROOT / (
    "runs/phase1A/day032_qm_f06_upper_v5c_construction"
)

XYZ_PATH = (
    CONSTRUCTION_DIR
    / "QM_F06_UPPER_V5C_start.xyz"
)

MAP_PATH = (
    CONSTRUCTION_DIR
    / "QM_F06_UPPER_V5C_atom_role_provenance_map.csv"
)

CONSTRUCTION_REPORT = (
    CONSTRUCTION_DIR
    / "QM_F06_UPPER_V5C_CONSTRUCTION_REPORT.json"
)

VALENCE_PATH = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_pre_qm_audit/"
    "QM_F06_UPPER_V5B_valence.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day032_qm_f06_upper_v5c_pre_qm_audit"
)

BOND_CSV = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_bond_audit.csv"
)

RECONNECTIVITY_CSV = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_geometric_reconnectivity.csv"
)

CONTACT_CSV = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_contacts.csv"
)

REPORT_PATH = OUTPUT_DIR / (
    "QM_F06_UPPER_V5C_PRE_QM_AUDIT.json"
)

EXPECTED_COUNT = 52

EXPECTED_COMPOSITION = Counter({
    "B": 16,
    "N": 13,
    "H": 23,
})

MODIFIED_ATOMS = {
    "S:1739",
    "BR4:UPPER:14:1",
    "H4:UPPER:0203:0",
}

MINIMUM_LOCAL_MARGIN_A = 0.04

NONNOMINAL_PAIRS = {
    tuple(sorted((
        "S:1739",
        "BR4:UPPER:00:3",
    ))),
    tuple(sorted((
        "BR4:UPPER:14:1",
        "BR4:UPPER:00:4",
    ))),
}

MINIMUM_NONNOMINAL_CLEARANCE_A = 1.95


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty file: {path}"
        )


def read_csv(path: Path):
    require_file(path)

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def read_xyz(path: Path):
    require_file(path)

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


def canonical_pair(first, second):
    return tuple(sorted((first, second)))


def bond_limits(first_element, second_element):
    pair = frozenset((
        first_element,
        second_element,
    ))

    if pair == frozenset(("B", "N")):
        return 1.25, 1.85, "B-N"

    if pair == frozenset(("B", "H")):
        return 0.95, 1.35, "B-H"

    if pair == frozenset(("N", "H")):
        return 0.85, 1.25, "N-H"

    return None


def main() -> None:
    mapping = read_csv(MAP_PATH)
    valence = read_csv(VALENCE_PATH)
    atoms = read_xyz(XYZ_PATH)

    construction = json.loads(
        CONSTRUCTION_REPORT.read_text(
            encoding="utf-8"
        )
    )

    index_by_id = {
        row["atom_id"]: int(row["index_0based"])
        for row in mapping
    }

    id_by_index = {
        int(row["index_0based"]): row["atom_id"]
        for row in mapping
    }

    element_by_id = {
        row["atom_id"]: row["element"]
        for row in mapping
    }

    composition = Counter(
        atom["element"]
        for atom in atoms
    )

    atom_count_gate = (
        len(atoms) == EXPECTED_COUNT
        and len(mapping) == EXPECTED_COUNT
    )

    composition_gate = (
        composition == EXPECTED_COMPOSITION
    )

    nominal_pairs = set()

    expected_degree = {}

    for row in valence:
        atom_id = row["atom_id"]

        expected_degree[atom_id] = int(
            row["expected_degree"]
        )

        for neighbor in row["neighbors"].split("|"):
            if neighbor:
                nominal_pairs.add(
                    canonical_pair(
                        atom_id,
                        neighbor,
                    )
                )

    bond_records = []
    bond_failures = []
    minimum_local_margin = float("inf")

    for first, second in sorted(nominal_pairs):
        limits = bond_limits(
            element_by_id[first],
            element_by_id[second],
        )

        if limits is None:
            raise RuntimeError(
                f"Unsupported nominal bond: "
                f"{first} -- {second}"
            )

        lower, upper, bond_class = limits

        value = distance(
            atoms[index_by_id[first]]["xyz"],
            atoms[index_by_id[second]]["xyz"],
        )

        margin = min(
            value - lower,
            upper - value,
        )

        local = (
            first in MODIFIED_ATOMS
            or second in MODIFIED_ATOMS
        )

        passed = (
            lower <= value <= upper
        )

        if local:
            minimum_local_margin = min(
                minimum_local_margin,
                margin,
            )

        record = {
            "first_atom": first,
            "second_atom": second,
            "bond_class": bond_class,
            "distance_A": value,
            "minimum_A": lower,
            "maximum_A": upper,
            "margin_A": margin,
            "modified_region_bond": local,
            "pass": passed,
        }

        bond_records.append(record)

        if not passed:
            bond_failures.append(record)

    nominal_bond_gate = (
        len(bond_failures) == 0
    )

    local_margin_gate = (
        minimum_local_margin
        >= MINIMUM_LOCAL_MARGIN_A
    )

    geometric_pairs = set()
    geometric_distances = {}

    for first_index in range(len(atoms)):
        for second_index in range(
            first_index + 1,
            len(atoms),
        ):
            first_id = id_by_index[first_index]
            second_id = id_by_index[second_index]

            limits = bond_limits(
                atoms[first_index]["element"],
                atoms[second_index]["element"],
            )

            if limits is None:
                continue

            lower, upper, _ = limits

            value = distance(
                atoms[first_index]["xyz"],
                atoms[second_index]["xyz"],
            )

            if lower <= value <= upper:
                pair = canonical_pair(
                    first_id,
                    second_id,
                )

                geometric_pairs.add(pair)
                geometric_distances[pair] = value

    gained = sorted(
        geometric_pairs - nominal_pairs
    )

    lost = sorted(
        nominal_pairs - geometric_pairs
    )

    reconnectivity_records = []

    for first, second in gained:
        reconnectivity_records.append({
            "classification": (
                "NONNOMINAL_GEOMETRIC_BOND"
            ),
            "first_atom": first,
            "second_atom": second,
            "distance_A": (
                geometric_distances[
                    canonical_pair(first, second)
                ]
            ),
        })

    for first, second in lost:
        reconnectivity_records.append({
            "classification": (
                "NOMINAL_BOND_OUTSIDE_RANGE"
            ),
            "first_atom": first,
            "second_atom": second,
            "distance_A": distance(
                atoms[index_by_id[first]]["xyz"],
                atoms[index_by_id[second]]["xyz"],
            ),
        })

    reconnectivity_gate = (
        len(gained) == 0
        and len(lost) == 0
    )

    adjacency = defaultdict(set)

    for first, second in geometric_pairs:
        adjacency[first].add(second)
        adjacency[second].add(first)

    overcoordinated = []

    for atom_id in index_by_id:
        degree = len(adjacency[atom_id])

        maximum = expected_degree[atom_id]

        if degree > maximum:
            overcoordinated.append({
                "atom_id": atom_id,
                "element": element_by_id[atom_id],
                "geometric_degree": degree,
                "maximum_degree": maximum,
                "neighbors": sorted(
                    adjacency[atom_id]
                ),
            })

    overcoordination_gate = (
        len(overcoordinated) == 0
    )

    start = next(iter(index_by_id))
    visited = {start}
    queue = deque([start])

    while queue:
        current = queue.popleft()

        for neighbor in adjacency[current]:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)

    connected_gate = (
        len(visited) == len(index_by_id)
    )

    nonnominal_distances = {
        f"{first}--{second}": distance(
            atoms[index_by_id[first]]["xyz"],
            atoms[index_by_id[second]]["xyz"],
        )
        for first, second in NONNOMINAL_PAIRS
    }

    minimum_nonnominal_clearance = min(
        nonnominal_distances.values()
    )

    nonnominal_clearance_gate = (
        minimum_nonnominal_clearance
        >= MINIMUM_NONNOMINAL_CLEARANCE_A
    )

    hard_contacts = []

    for first_index in range(len(atoms)):
        for second_index in range(
            first_index + 1,
            len(atoms),
        ):
            first = atoms[first_index]
            second = atoms[second_index]

            pair = canonical_pair(
                id_by_index[first_index],
                id_by_index[second_index],
            )

            if pair in nominal_pairs:
                continue

            value = distance(
                first["xyz"],
                second["xyz"],
            )

            if (
                first["element"] == "H"
                and second["element"] == "H"
            ):
                threshold = 1.20
            elif (
                first["element"] == "H"
                or second["element"] == "H"
            ):
                threshold = 0.85
            else:
                threshold = 1.20

            if value < threshold:
                hard_contacts.append({
                    "first_atom": (
                        id_by_index[first_index]
                    ),
                    "second_atom": (
                        id_by_index[second_index]
                    ),
                    "distance_A": value,
                    "threshold_A": threshold,
                })

    contact_gate = (
        len(hard_contacts) == 0
    )

    gates = {
        "construction_authorized": (
            construction[
                "pre_qm_audit_authorized"
            ]
        ),
        "atom_count": atom_count_gate,
        "composition": composition_gate,
        "all_nominal_bonds_in_range": (
            nominal_bond_gate
        ),
        "modified_region_bond_margin": (
            local_margin_gate
        ),
        "no_geometric_reconnectivity": (
            reconnectivity_gate
        ),
        "no_overcoordinated_atoms": (
            overcoordination_gate
        ),
        "single_connected_component": (
            connected_gate
        ),
        "nonnominal_pair_clearance": (
            nonnominal_clearance_gate
        ),
        "no_hard_contacts": contact_gate,
    }

    passed = all(gates.values())

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    with BOND_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(bond_records[0]),
        )
        writer.writeheader()
        writer.writerows(bond_records)

    reconnectivity_fieldnames = [
        "classification",
        "first_atom",
        "second_atom",
        "distance_A",
    ]

    with RECONNECTIVITY_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=reconnectivity_fieldnames,
        )
        writer.writeheader()
        writer.writerows(
            reconnectivity_records
        )

    contact_fieldnames = [
        "first_atom",
        "second_atom",
        "distance_A",
        "threshold_A",
    ]

    with CONTACT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=contact_fieldnames,
        )
        writer.writeheader()
        writer.writerows(hard_contacts)

    decision = (
        "QM_F06_UPPER_V5C_PRE_QM_GATE_PASS_"
        "ORCA_INPUT_DESIGN_AUTHORIZED"
        if passed
        else
        "QM_F06_UPPER_V5C_PRE_QM_GATE_FAIL_"
        "CONSTRUCTION_REVIEW_REQUIRED"
    )

    report = {
        "decision": decision,
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "gates": gates,
        "minimum_modified_region_bond_margin_A": (
            minimum_local_margin
        ),
        "minimum_nonnominal_clearance_A": (
            minimum_nonnominal_clearance
        ),
        "nonnominal_pair_distances_A": (
            nonnominal_distances
        ),
        "bond_failure_count": len(
            bond_failures
        ),
        "gained_geometric_bond_count": len(
            gained
        ),
        "lost_nominal_bond_count": len(
            lost
        ),
        "overcoordinated_atom_count": len(
            overcoordinated
        ),
        "overcoordinated_atoms": (
            overcoordinated
        ),
        "hard_contact_count": len(
            hard_contacts
        ),
        "authorization": {
            "orca_input_design_authorized": passed,
            "orca_execution_authorized": False,
            "RESP_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    REPORT_PATH.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 92)
    print("QM_F06 UPPER V5-C FORMAL PRE-QM AUDIT")
    print("=" * 92)

    for name, value in gates.items():
        print(
            f"{name:42s}: "
            f"{'PASS' if value else 'FAIL'}"
        )

    print()
    print(
        "Minimum modified-region bond margin A:",
        minimum_local_margin,
    )
    print(
        "Minimum nonnominal clearance A:",
        minimum_nonnominal_clearance,
    )
    print("Bond failures:", len(bond_failures))
    print("Gained geometric bonds:", len(gained))
    print("Lost nominal bonds:", len(lost))
    print(
        "Overcoordinated atoms:",
        len(overcoordinated),
    )
    print("Hard contacts:", len(hard_contacts))

    print()
    print("Decision:", decision)
    print("Report:", REPORT_PATH)
    print()
    print(
        "ORCA input design authorized:",
        passed,
    )
    print("ORCA execution authorized: False")
    print("RESP authorized: False")
    print("MD authorized: False")


if __name__ == "__main__":
    main()
