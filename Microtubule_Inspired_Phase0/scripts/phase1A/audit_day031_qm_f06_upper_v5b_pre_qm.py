#!/usr/bin/env python3
"""
Pre-QM structural audit for QM_F06 UPPER V5-B.

Connectivity is reconstructed from the validated V4 valence graph,
retaining only atoms that survive in V5-B, then adding:

- S:1738 -- BR4:UPPER:14:1
- S:1738 -- P:1637
- P:1637 -- S:1737
- the three new artificial cap bonds

The audit validates:
- construction authorization;
- atom count, composition and map consistency;
- one connected component;
- nominal B/N/H valence;
- critical restored and repaired geometry;
- new-cap geometry and ownership;
- absence of non-topological hard contacts.

ORCA execution remains blocked.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

CONSTRUCTION_DIR = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_construction"
)

XYZ_PATH = (
    CONSTRUCTION_DIR
    / "QM_F06_UPPER_V5B_start.xyz"
)

MAP_PATH = (
    CONSTRUCTION_DIR
    / "QM_F06_UPPER_V5B_atom_role_provenance_map.csv"
)

NEW_CAPS_PATH = (
    CONSTRUCTION_DIR
    / "QM_F06_UPPER_V5B_new_artificial_caps.csv"
)

CONSTRUCTION_REPORT = (
    CONSTRUCTION_DIR
    / "QM_F06_UPPER_V5B_CONSTRUCTION_REPORT.json"
)

V4_VALENCE = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_pre_qm_audit/"
    "QM_F06_UPPER_V4_valence.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5b_pre_qm_audit"
)

REPORT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5B_PRE_QM_AUDIT.json"
)

VALENCE_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5B_valence.csv"
)

CAP_AUDIT_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5B_new_cap_audit.csv"
)

CONTACTS_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5B_contacts.csv"
)

CRITICAL_GEOMETRY_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5B_critical_geometry.csv"
)

EXPECTED_ATOM_COUNT = 52

EXPECTED_COMPOSITION = Counter({
    "B": 16,
    "N": 13,
    "H": 23,
})

NEW_CANONICAL_H_BONDS = {
    tuple(sorted((
        "BR4:UPPER:14:1",
        "H4:UPPER:0170:0",
    ))),
    tuple(sorted((
        "S:1737",
        "H4:UPPER:0202:0",
    ))),
}

NEW_REAL_BONDS = {
    tuple(sorted((
        "S:1738",
        "BR4:UPPER:14:1",
    ))),
    tuple(sorted((
        "S:1738",
        "P:1637",
    ))),
    tuple(sorted((
        "P:1637",
        "S:1737",
    ))),
}

CRITICAL_BONDS = {
    tuple(sorted((
        "A:UPPER:13:3",
        "A:UPPER:14:4",
    ))): (1.25, 1.85),

    tuple(sorted((
        "S:1739",
        "H4:UPPER:0203:0",
    ))): (0.95, 1.35),

    tuple(sorted((
        "S:1738",
        "BR4:UPPER:14:1",
    ))): (1.25, 1.85),

    tuple(sorted((
        "S:1738",
        "P:1637",
    ))): (1.25, 1.85),

    tuple(sorted((
        "P:1637",
        "S:1737",
    ))): (1.25, 1.85),
}

BN_MIN_A = 1.20
BN_MAX_A = 1.90

BH_MIN_A = 0.95
BH_MAX_A = 1.35

NH_MIN_A = 0.85
NH_MAX_A = 1.25

HH_HARD_CONTACT_A = 1.20
HX_HARD_CONTACT_A = 0.85
HEAVY_HEAVY_HARD_CONTACT_A = 1.20

NEW_CAP_OWNER_MARGIN_A = 0.10

REBUILT_H_MIN_DISTANCE_TO_OTHER_N_A = 2.00


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


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {
        "true",
        "1",
        "yes",
        "pass",
    }


def canonical_pair(first: str, second: str):
    return tuple(sorted((first, second)))


def distance(first, second) -> float:
    return math.sqrt(sum(
        (a - b) ** 2
        for a, b in zip(first, second)
    ))


def main() -> None:
    for path in (
        XYZ_PATH,
        MAP_PATH,
        NEW_CAPS_PATH,
        CONSTRUCTION_REPORT,
        V4_VALENCE,
    ):
        require_file(path)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    construction = json.loads(
        CONSTRUCTION_REPORT.read_text(
            encoding="utf-8"
        )
    )

    if not construction["overall_pass"]:
        raise RuntimeError(
            "V5-B construction report did not pass."
        )

    if not construction["authorization"][
        "pre_qm_structural_audit_authorized"
    ]:
        raise RuntimeError(
            "V5-B pre-QM audit is not authorized."
        )

    atoms = read_xyz(XYZ_PATH)
    map_rows = read_csv(MAP_PATH)
    new_cap_rows = read_csv(NEW_CAPS_PATH)
    v4_valence_rows = read_csv(V4_VALENCE)

    if len(atoms) != len(map_rows):
        raise RuntimeError(
            "XYZ and provenance-map row counts differ."
        )

    mapped = {}

    for atom, row in zip(
        atoms,
        map_rows,
        strict=True,
    ):
        index = int(row["index_0based"])

        if index != atom["index"]:
            raise RuntimeError(
                f"Index mismatch at row {index}"
            )

        atom_id = row["atom_id"]

        if atom_id in mapped:
            raise RuntimeError(
                f"Duplicate atom ID: {atom_id}"
            )

        if row["element"] != atom["element"]:
            raise RuntimeError(
                f"Element mismatch for {atom_id}"
            )

        mapped[atom_id] = {
            **atom,
            **row,
        }

    composition = Counter(
        atom["element"]
        for atom in atoms
    )

    bonded_pairs = set()

    # Import every validated V4 bond whose endpoints survive.
    for row in v4_valence_rows:
        first = row["atom_id"]

        if first not in mapped:
            continue

        neighbors = [
            value.strip()
            for value in row["neighbors"].split("|")
            if value.strip()
        ]

        for second in neighbors:
            if second in mapped:
                bonded_pairs.add(
                    canonical_pair(first, second)
                )

    # Add the three canonical real-atom bonds restored in V5-B.
    bonded_pairs.update(NEW_REAL_BONDS)

    # Add the two canonical passivant-H bonds required
    # by the selected R2 chemical graph.
    bonded_pairs.update(NEW_CANONICAL_H_BONDS)

    # Add the three new artificial-cap bonds.
    for row in new_cap_rows:
        cap_id = row["cap_id"]
        center = row["center_atom"]

        if cap_id not in mapped:
            raise RuntimeError(
                f"New cap absent from V5-B map: {cap_id}"
            )

        if center not in mapped:
            raise RuntimeError(
                f"New cap center absent: {center}"
            )

        bonded_pairs.add(
            canonical_pair(cap_id, center)
        )

    adjacency = defaultdict(set)

    for first, second in bonded_pairs:
        adjacency[first].add(second)
        adjacency[second].add(first)

    atom_ids = set(mapped)

    # Connected-component gate.
    visited = set()
    queue = deque([next(iter(atom_ids))])

    while queue:
        current = queue.popleft()

        if current in visited:
            continue

        visited.add(current)

        for neighbor in adjacency[current]:
            if neighbor not in visited:
                queue.append(neighbor)

    connected_component_gate = (
        visited == atom_ids
    )

    valence_records = []
    valence_failures = []

    for atom_id, record in sorted(
        mapped.items(),
        key=lambda item: item[1]["index"],
    ):
        element = record["element"]
        degree = len(adjacency[atom_id])

        expected_degree = {
            "B": 3,
            "N": 3,
            "H": 1,
        }[element]

        passed = degree == expected_degree

        valence_record = {
            "index_0based": record["index"],
            "atom_id": atom_id,
            "element": element,
            "degree": degree,
            "expected_degree": expected_degree,
            "neighbors": "|".join(
                sorted(adjacency[atom_id])
            ),
            "pass": passed,
        }

        valence_records.append(valence_record)

        if not passed:
            valence_failures.append(
                valence_record
            )

    with VALENCE_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(valence_records[0]),
        )
        writer.writeheader()
        writer.writerows(valence_records)

    valence_gate = (
        len(valence_failures) == 0
    )

    # Validate critical geometry.
    critical_records = []
    critical_failures = []

    for pair, limits in sorted(
        CRITICAL_BONDS.items()
    ):
        first, second = pair

        if first not in mapped or second not in mapped:
            raise RuntimeError(
                f"Missing critical pair: {first} -- {second}"
            )

        value = distance(
            mapped[first]["xyz_A"],
            mapped[second]["xyz_A"],
        )

        minimum_A, maximum_A = limits
        passed = minimum_A <= value <= maximum_A

        critical_record = {
            "first_atom": first,
            "second_atom": second,
            "distance_A": value,
            "minimum_A": minimum_A,
            "maximum_A": maximum_A,
            "pass": passed,
        }

        critical_records.append(
            critical_record
        )

        if not passed:
            critical_failures.append(
                critical_record
            )

    rebuilt_h = "H4:UPPER:0203:0"

    for nitrogen in (
        "P:1641",
        "P:1639",
    ):
        value = distance(
            mapped[rebuilt_h]["xyz_A"],
            mapped[nitrogen]["xyz_A"],
        )

        passed = (
            value
            >= REBUILT_H_MIN_DISTANCE_TO_OTHER_N_A
        )

        record = {
            "first_atom": rebuilt_h,
            "second_atom": nitrogen,
            "distance_A": value,
            "minimum_A": (
                REBUILT_H_MIN_DISTANCE_TO_OTHER_N_A
            ),
            "maximum_A": "",
            "pass": passed,
        }

        critical_records.append(record)

        if not passed:
            critical_failures.append(record)

    with CRITICAL_GEOMETRY_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "first_atom",
                "second_atom",
                "distance_A",
                "minimum_A",
                "maximum_A",
                "pass",
            ],
        )
        writer.writeheader()
        writer.writerows(critical_records)

    critical_geometry_gate = (
        len(critical_failures) == 0
    )

    # New-cap geometry and ownership.
    cap_audit_records = []
    cap_failures = []

    for row in new_cap_rows:
        cap_id = row["cap_id"]
        center_id = row["center_atom"]

        cap = mapped[cap_id]
        center = mapped[center_id]

        owner_distance_A = distance(
            cap["xyz_A"],
            center["xyz_A"],
        )

        heavy_candidates = sorted(
            (
                distance(
                    cap["xyz_A"],
                    record["xyz_A"],
                ),
                atom_id,
                record["element"],
            )
            for atom_id, record in mapped.items()
            if record["element"] in {"B", "N"}
        )

        nearest_distance_A, nearest_id, _ = (
            heavy_candidates[0]
        )

        second_distance_A = (
            heavy_candidates[1][0]
        )

        ownership_pass = (
            nearest_id == center_id
            and (
                second_distance_A
                - nearest_distance_A
            ) >= NEW_CAP_OWNER_MARGIN_A
        )

        center_element = center["element"]

        if center_element == "B":
            bond_range_pass = (
                BH_MIN_A
                <= owner_distance_A
                <= BH_MAX_A
            )
        elif center_element == "N":
            bond_range_pass = (
                NH_MIN_A
                <= owner_distance_A
                <= NH_MAX_A
            )
        else:
            bond_range_pass = False

        passed = (
            ownership_pass
            and bond_range_pass
        )

        record = {
            "cap_id": cap_id,
            "center_atom": center_id,
            "center_element": center_element,
            "owner_distance_A": owner_distance_A,
            "nearest_heavy_atom": nearest_id,
            "nearest_heavy_distance_A": nearest_distance_A,
            "second_heavy_distance_A": second_distance_A,
            "nearest_center_margin_A": (
                second_distance_A
                - nearest_distance_A
            ),
            "ownership_pass": ownership_pass,
            "bond_range_pass": bond_range_pass,
            "pass": passed,
        }

        cap_audit_records.append(record)

        if not passed:
            cap_failures.append(record)

    with CAP_AUDIT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(cap_audit_records[0]),
        )
        writer.writeheader()
        writer.writerows(cap_audit_records)

    cap_gate = (
        len(cap_failures) == 0
    )

    # Graph distances for topological contact exclusions.
    graph_distances = {}

    for origin in atom_ids:
        distances = {origin: 0}
        local_queue = deque([origin])

        while local_queue:
            current = local_queue.popleft()

            for neighbor in adjacency[current]:
                if neighbor not in distances:
                    distances[neighbor] = (
                        distances[current] + 1
                    )
                    local_queue.append(neighbor)

        for target, separation in distances.items():
            graph_distances[
                canonical_pair(origin, target)
            ] = separation

    contact_records = []
    hard_contacts = []

    ordered = sorted(
        mapped.items(),
        key=lambda item: item[1]["index"],
    )

    for position, (first_id, first) in enumerate(
        ordered
    ):
        for second_id, second in ordered[
            position + 1:
        ]:
            value = distance(
                first["xyz_A"],
                second["xyz_A"],
            )

            separation = graph_distances.get(
                canonical_pair(first_id, second_id)
            )

            if separation in {1, 2}:
                classification = (
                    "TOPOLOGICAL_EXCLUSION"
                )
                is_hard = False
            else:
                if (
                    first["element"] == "H"
                    and second["element"] == "H"
                ):
                    threshold = HH_HARD_CONTACT_A
                elif (
                    first["element"] == "H"
                    or second["element"] == "H"
                ):
                    threshold = HX_HARD_CONTACT_A
                else:
                    threshold = (
                        HEAVY_HEAVY_HARD_CONTACT_A
                    )

                is_hard = value < threshold

                classification = (
                    "HARD_CONTACT"
                    if is_hard
                    else "NONBONDED_OK"
                )

            record = {
                "first_atom": first_id,
                "first_element": first["element"],
                "second_atom": second_id,
                "second_element": second["element"],
                "distance_A": value,
                "graph_separation": separation,
                "classification": classification,
                "hard_contact": is_hard,
            }

            contact_records.append(record)

            if is_hard:
                hard_contacts.append(record)

    with CONTACTS_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(contact_records[0]),
        )
        writer.writeheader()
        writer.writerows(contact_records)

    hard_contact_gate = (
        len(hard_contacts) == 0
    )

    gates = {
        "construction_report": True,
        "atom_count": (
            len(atoms) == EXPECTED_ATOM_COUNT
        ),
        "composition": (
            composition == EXPECTED_COMPOSITION
        ),
        "map_consistency": (
            len(mapped) == EXPECTED_ATOM_COUNT
        ),
        "single_connected_component": (
            connected_component_gate
        ),
        "nominal_valence": valence_gate,
        "critical_initial_geometry": (
            critical_geometry_gate
        ),
        "new_artificial_cap_geometry": (
            cap_gate
        ),
        "no_unresolved_hard_contacts": (
            hard_contact_gate
        ),
    }

    overall_pass = all(gates.values())

    decision = (
        "QM_F06_UPPER_V5B_PRE_QM_STRUCTURAL_GATE_PASS_"
        "CONSTRAINT_DESIGN_AUTHORIZED"
        if overall_pass
        else
        "QM_F06_UPPER_V5B_PRE_QM_STRUCTURAL_GATE_FAIL_"
        "CONSTRUCTION_REVIEW_REQUIRED"
    )

    report = {
        "decision": decision,
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "atom_count": len(atoms),
        "composition": dict(
            sorted(composition.items())
        ),
        "gates": gates,
        "valence_failure_count": len(
            valence_failures
        ),
        "valence_failures": valence_failures,
        "critical_geometry_failure_count": len(
            critical_failures
        ),
        "critical_geometry_failures": (
            critical_failures
        ),
        "new_cap_failure_count": len(
            cap_failures
        ),
        "new_cap_failures": cap_failures,
        "hard_contact_count": len(
            hard_contacts
        ),
        "hard_contacts": hard_contacts,
        "files": {
            "valence": str(
                VALENCE_CSV.relative_to(ROOT)
            ),
            "critical_geometry": str(
                CRITICAL_GEOMETRY_CSV.relative_to(ROOT)
            ),
            "new_cap_audit": str(
                CAP_AUDIT_CSV.relative_to(ROOT)
            ),
            "contacts": str(
                CONTACTS_CSV.relative_to(ROOT)
            ),
        },
        "files_sha256": {
            "xyz": sha256(XYZ_PATH),
            "map": sha256(MAP_PATH),
            "new_caps": sha256(
                NEW_CAPS_PATH
            ),
            "construction_report": sha256(
                CONSTRUCTION_REPORT
            ),
            "v4_valence": sha256(
                V4_VALENCE
            ),
            "valence": sha256(VALENCE_CSV),
            "critical_geometry": sha256(
                CRITICAL_GEOMETRY_CSV
            ),
            "new_cap_audit": sha256(
                CAP_AUDIT_CSV
            ),
            "contacts": sha256(
                CONTACTS_CSV
            ),
        },
        "authorization": {
            "constraint_design_authorized": (
                overall_pass
            ),
            "orca_input_preparation_authorized": False,
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

    print("=" * 82)
    print("QM_F06 UPPER V5-B PRE-QM STRUCTURAL AUDIT")
    print("=" * 82)

    for gate, passed in gates.items():
        print(
            f"{gate:40s}: "
            f"{'PASS' if passed else 'FAIL'}"
        )

    print()
    print("Atom count:", len(atoms))
    print("Composition:", dict(composition))
    print(
        "Valence failures:",
        len(valence_failures),
    )
    print(
        "Critical-geometry failures:",
        len(critical_failures),
    )
    print(
        "New-cap failures:",
        len(cap_failures),
    )
    print(
        "Hard contacts:",
        len(hard_contacts),
    )

    print()
    print("Decision:", decision)
    print("Report:", REPORT_PATH)
    print("Valence:", VALENCE_CSV)
    print(
        "Critical geometry:",
        CRITICAL_GEOMETRY_CSV,
    )
    print("New caps:", CAP_AUDIT_CSV)
    print("Contacts:", CONTACTS_CSV)
    print()
    print(
        "Constraint design authorized:",
        overall_pass,
    )
    print("ORCA execution authorized: False")


if __name__ == "__main__":
    main()
