#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter, deque
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

INPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day033_qm_f06_upper_v6a_topology_closure"
)

XYZ_PATH = INPUT_DIR / (
    "QM_F06_UPPER_V6A_TOPOLOGY_CLOSURE_start.xyz"
)

MAP_PATH = INPUT_DIR / (
    "QM_F06_UPPER_V6A_atom_role_provenance_map.csv"
)

EDGE_PATH = INPUT_DIR / (
    "QM_F06_UPPER_V6A_nominal_edges.csv"
)

CONSTRUCTION_REPORT = INPUT_DIR / (
    "QM_F06_UPPER_V6A_TOPOLOGY_CLOSURE_REPORT.json"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day033_qm_f06_upper_v6a_pre_qm_audit"
)

BOND_CSV = OUTPUT_DIR / (
    "QM_F06_UPPER_V6A_bond_audit.csv"
)

RECONNECTIVITY_CSV = OUTPUT_DIR / (
    "QM_F06_UPPER_V6A_geometric_reconnectivity.csv"
)

CONTACT_CSV = OUTPUT_DIR / (
    "QM_F06_UPPER_V6A_contacts.csv"
)

REPORT_PATH = OUTPUT_DIR / (
    "QM_F06_UPPER_V6A_PRE_QM_AUDIT.json"
)


EXPECTED_ATOM_COUNT = 48

EXPECTED_COMPOSITION = Counter({
    "B": 16,
    "N": 13,
    "H": 19,
})

EXPECTED_DEGREE = {
    "B": 3,
    "N": 3,
    "H": 1,
}

BN_MIN_A = 1.25
BN_MAX_A = 1.90

BH_MIN_A = 0.90
BH_MAX_A = 1.35

NH_MIN_A = 0.80
NH_MAX_A = 1.25

STRICT_NONNOMINAL_BN_A = 1.90

HH_HARD_CONTACT_A = 0.70
HX_HARD_CONTACT_A = 0.75
HEAVY_HEAVY_HARD_CONTACT_A = 1.10

NEW_CLOSURE_BONDS = {
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

    count = int(lines[0].strip())
    coordinate_lines = lines[2:2 + count]

    if len(coordinate_lines) != count:
        raise RuntimeError(
            "Incomplete XYZ file."
        )

    atoms = []

    for index, line in enumerate(
        coordinate_lines
    ):
        fields = line.split()

        atoms.append({
            "index": index,
            "element": fields[0],
            "xyz_A": tuple(
                map(float, fields[1:4])
            ),
        })

    return atoms


def canonical_pair(
    first: str,
    second: str,
) -> tuple[str, str]:
    return tuple(sorted((first, second)))


def distance(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> float:
    return math.sqrt(
        sum(
            (a - b) ** 2
            for a, b in zip(first, second)
        )
    )


def nominal_range(
    first_element: str,
    second_element: str,
) -> tuple[float, float] | None:
    pair = {
        first_element,
        second_element,
    }

    if pair == {"B", "N"}:
        return BN_MIN_A, BN_MAX_A

    if pair == {"B", "H"}:
        return BH_MIN_A, BH_MAX_A

    if pair == {"N", "H"}:
        return NH_MIN_A, NH_MAX_A

    return None


def geometric_bond(
    first_element: str,
    second_element: str,
    value: float,
) -> bool:
    limits = nominal_range(
        first_element,
        second_element,
    )

    if limits is None:
        return False

    minimum, maximum = limits

    return minimum <= value <= maximum


def main() -> None:
    for path in (
        XYZ_PATH,
        MAP_PATH,
        EDGE_PATH,
        CONSTRUCTION_REPORT,
    ):
        require_file(path)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    atoms = read_xyz(XYZ_PATH)
    map_rows = read_csv(MAP_PATH)
    edge_rows = read_csv(EDGE_PATH)

    retained_map_rows = [
        row
        for row in map_rows
        if row["v6a_retained"].strip().lower()
        == "true"
    ]

    retained_map_rows.sort(
        key=lambda row: int(
            row["v6a_index_0based"]
        )
    )

    if len(atoms) != len(retained_map_rows):
        raise RuntimeError(
            "XYZ/map retained-atom mismatch."
        )

    mapped = {}

    for atom, row in zip(
        atoms,
        retained_map_rows,
    ):
        if atom["element"] != row["element"]:
            raise RuntimeError(
                "Element mismatch at V6-A index "
                f"{atom['index']}."
            )

        mapped[row["atom_id"]] = {
            "index": atom["index"],
            "element": atom["element"],
            "xyz_A": atom["xyz_A"],
        }

    atom_ids = set(mapped)

    nominal_edges = set()

    for row in edge_rows:
        pair = canonical_pair(
            row["first_atom"],
            row["second_atom"],
        )

        if (
            pair[0] not in atom_ids
            or pair[1] not in atom_ids
        ):
            raise RuntimeError(
                f"Nominal edge references absent atom: {pair}"
            )

        nominal_edges.add(pair)

    composition = Counter(
        atom["element"]
        for atom in mapped.values()
    )

    atom_count_gate = (
        len(mapped) == EXPECTED_ATOM_COUNT
    )

    composition_gate = (
        composition == EXPECTED_COMPOSITION
    )

    adjacency = {
        atom_id: set()
        for atom_id in atom_ids
    }

    for first, second in nominal_edges:
        adjacency[first].add(second)
        adjacency[second].add(first)

    degree_failures = []

    for atom_id in sorted(atom_ids):
        element = mapped[atom_id]["element"]
        degree = len(adjacency[atom_id])
        expected = EXPECTED_DEGREE[element]

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

    nominal_degree_gate = (
        len(degree_failures) == 0
    )

    start_atom = sorted(atom_ids)[0]

    visited = {start_atom}
    queue = deque([start_atom])

    while queue:
        current = queue.popleft()

        for neighbor in adjacency[current]:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)

    connected_gate = (
        len(visited) == len(atom_ids)
    )

    bond_records = []
    bond_failures = []

    for first, second in sorted(
        nominal_edges
    ):
        first_atom = mapped[first]
        second_atom = mapped[second]

        value = distance(
            first_atom["xyz_A"],
            second_atom["xyz_A"],
        )

        limits = nominal_range(
            first_atom["element"],
            second_atom["element"],
        )

        if limits is None:
            minimum = None
            maximum = None
            passed = False
            bond_class = (
                first_atom["element"]
                + "-"
                + second_atom["element"]
            )
        else:
            minimum, maximum = limits
            passed = minimum <= value <= maximum
            bond_class = "-".join(sorted((
                first_atom["element"],
                second_atom["element"],
            )))

        record = {
            "first_atom": first,
            "first_element": first_atom["element"],
            "second_atom": second,
            "second_element": second_atom["element"],
            "bond_class": bond_class,
            "distance_A": value,
            "minimum_A": minimum,
            "maximum_A": maximum,
            "new_closure_bond": (
                canonical_pair(first, second)
                in NEW_CLOSURE_BONDS
            ),
            "pass": passed,
        }

        bond_records.append(record)

        if not passed:
            bond_failures.append(record)

    with BOND_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                bond_records[0].keys()
            ),
        )

        writer.writeheader()
        writer.writerows(bond_records)

    bond_gate = (
        len(bond_failures) == 0
    )

    ordered_ids = sorted(
        atom_ids,
        key=lambda atom_id: mapped[
            atom_id
        ]["index"],
    )

    geometric_edges = set()
    geometric_distances = {}

    contact_records = []
    hard_contacts = []

    for first_position, first in enumerate(
        ordered_ids
    ):
        first_atom = mapped[first]

        for second in ordered_ids[
            first_position + 1:
        ]:
            second_atom = mapped[second]

            value = distance(
                first_atom["xyz_A"],
                second_atom["xyz_A"],
            )

            pair = canonical_pair(
                first,
                second,
            )

            is_geometric_bond = geometric_bond(
                first_atom["element"],
                second_atom["element"],
                value,
            )

            if is_geometric_bond:
                geometric_edges.add(pair)
                geometric_distances[pair] = value

            if (
                first_atom["element"] == "H"
                and second_atom["element"] == "H"
            ):
                hard_threshold = (
                    HH_HARD_CONTACT_A
                )
            elif (
                first_atom["element"] == "H"
                or second_atom["element"] == "H"
            ):
                hard_threshold = (
                    HX_HARD_CONTACT_A
                )
            else:
                hard_threshold = (
                    HEAVY_HEAVY_HARD_CONTACT_A
                )

            hard_contact = (
                value < hard_threshold
            )

            contact_record = {
                "first_atom": first,
                "first_element": first_atom["element"],
                "second_atom": second,
                "second_element": second_atom["element"],
                "distance_A": value,
                "nominal_edge": pair in nominal_edges,
                "geometric_edge": is_geometric_bond,
                "hard_contact_threshold_A": (
                    hard_threshold
                ),
                "hard_contact": hard_contact,
            }

            contact_records.append(
                contact_record
            )

            if hard_contact:
                hard_contacts.append(
                    contact_record
                )

    with CONTACT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                contact_records[0].keys()
            ),
        )

        writer.writeheader()
        writer.writerows(contact_records)

    gained_edges = sorted(
        geometric_edges - nominal_edges
    )

    lost_edges = sorted(
        nominal_edges - geometric_edges
    )

    geometric_adjacency = {
        atom_id: set()
        for atom_id in atom_ids
    }

    for first, second in geometric_edges:
        geometric_adjacency[first].add(
            second
        )
        geometric_adjacency[second].add(
            first
        )

    overcoordinated = []

    reconnectivity_records = []

    for atom_id in ordered_ids:
        element = mapped[atom_id]["element"]

        nominal_neighbors = sorted(
            adjacency[atom_id]
        )

        geometric_neighbors = sorted(
            geometric_adjacency[atom_id]
        )

        lost_neighbors = sorted(
            set(nominal_neighbors)
            - set(geometric_neighbors)
        )

        gained_neighbors = sorted(
            set(geometric_neighbors)
            - set(nominal_neighbors)
        )

        geometric_degree = len(
            geometric_neighbors
        )

        expected_degree = (
            EXPECTED_DEGREE[element]
        )

        if geometric_degree > expected_degree:
            overcoordinated.append({
                "atom_id": atom_id,
                "element": element,
                "geometric_degree": geometric_degree,
                "expected_degree": expected_degree,
                "geometric_neighbors": (
                    geometric_neighbors
                ),
            })

        reconnectivity_records.append({
            "atom_id": atom_id,
            "element": element,
            "nominal_degree": len(
                nominal_neighbors
            ),
            "geometric_degree": (
                geometric_degree
            ),
            "nominal_neighbors": "|".join(
                nominal_neighbors
            ),
            "geometric_neighbors": "|".join(
                geometric_neighbors
            ),
            "lost_nominal_neighbors": "|".join(
                lost_neighbors
            ),
            "gained_geometric_neighbors": "|".join(
                gained_neighbors
            ),
            "reconnectivity_detected": bool(
                lost_neighbors
                or gained_neighbors
            ),
        })

    with RECONNECTIVITY_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                reconnectivity_records[0].keys()
            ),
        )

        writer.writeheader()
        writer.writerows(
            reconnectivity_records
        )

    reconnectivity_gate = (
        len(gained_edges) == 0
        and len(lost_edges) == 0
    )

    overcoordination_gate = (
        len(overcoordinated) == 0
    )

    contact_gate = (
        len(hard_contacts) == 0
    )

    closure_distances = {}

    closure_gate = True

    for pair in sorted(
        NEW_CLOSURE_BONDS
    ):
        value = distance(
            mapped[pair[0]]["xyz_A"],
            mapped[pair[1]]["xyz_A"],
        )

        closure_distances[
            "--".join(pair)
        ] = value

        closure_gate = (
            closure_gate
            and BN_MIN_A <= value <= BN_MAX_A
        )

    gates = {
        "construction_report_present": True,
        "atom_count": atom_count_gate,
        "composition": composition_gate,
        "nominal_degree_exact": nominal_degree_gate,
        "single_connected_component": connected_gate,
        "all_nominal_bonds_in_range": bond_gate,
        "closure_BN_geometry": closure_gate,
        "no_geometric_reconnectivity": reconnectivity_gate,
        "no_overcoordinated_atoms": overcoordination_gate,
        "no_hard_contacts": contact_gate,
    }

    passed = all(gates.values())

    decision = (
        "QM_F06_UPPER_V6A_PRE_QM_GATE_PASS_"
        "ORCA_INPUT_DESIGN_AUTHORIZED"
        if passed
        else
        "QM_F06_UPPER_V6A_PRE_QM_GATE_FAIL_"
        "TOPOLOGY_OR_GEOMETRY_REVIEW_REQUIRED"
    )

    files_for_hash = {
        "start_xyz": XYZ_PATH,
        "provenance_map": MAP_PATH,
        "nominal_edges": EDGE_PATH,
        "construction_report": (
            CONSTRUCTION_REPORT
        ),
        "bond_audit": BOND_CSV,
        "reconnectivity_audit": (
            RECONNECTIVITY_CSV
        ),
        "contact_audit": CONTACT_CSV,
    }

    report = {
        "decision": decision,
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "atom_count": len(mapped),
        "composition": dict(
            sorted(composition.items())
        ),
        "nominal_edge_count": len(
            nominal_edges
        ),
        "geometric_edge_count": len(
            geometric_edges
        ),
        "gates": gates,
        "pass": passed,
        "degree_failure_count": len(
            degree_failures
        ),
        "bond_failure_count": len(
            bond_failures
        ),
        "gained_geometric_edge_count": len(
            gained_edges
        ),
        "gained_geometric_edges": [
            {
                "first_atom": pair[0],
                "second_atom": pair[1],
                "distance_A": (
                    geometric_distances[pair]
                ),
            }
            for pair in gained_edges
        ],
        "lost_nominal_edge_count": len(
            lost_edges
        ),
        "lost_nominal_edges": [
            list(pair)
            for pair in lost_edges
        ],
        "overcoordinated_atom_count": len(
            overcoordinated
        ),
        "overcoordinated_atoms": (
            overcoordinated
        ),
        "hard_contact_count": len(
            hard_contacts
        ),
        "closure_BN_distances_A": (
            closure_distances
        ),
        "files": {
            key: str(path.relative_to(ROOT))
            for key, path in files_for_hash.items()
        },
        "sha256": {
            key: sha256(path)
            for key, path in files_for_hash.items()
        },
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

    print("=" * 104)
    print("QM_F06 UPPER V6-A FORMAL PRE-QM AUDIT")
    print("=" * 104)

    for name, value in gates.items():
        print(
            f"{name:44s}: "
            f"{'PASS' if value else 'FAIL'}"
        )

    print()
    print("Atoms:", len(mapped))
    print(
        "Composition:",
        dict(sorted(composition.items())),
    )
    print(
        "Nominal edges:",
        len(nominal_edges),
    )
    print(
        "Geometric edges:",
        len(geometric_edges),
    )
    print(
        "Degree failures:",
        len(degree_failures),
    )
    print(
        "Bond failures:",
        len(bond_failures),
    )
    print(
        "Gained geometric edges:",
        len(gained_edges),
    )
    print(
        "Lost nominal edges:",
        len(lost_edges),
    )
    print(
        "Overcoordinated atoms:",
        len(overcoordinated),
    )
    print(
        "Hard contacts:",
        len(hard_contacts),
    )

    print()
    print("New closure B-N distances:")

    for pair, value in sorted(
        closure_distances.items()
    ):
        print(
            f"  {pair:66s} "
            f"{value:.6f} Å"
        )

    if gained_edges:
        print()
        print("Nonnominal geometric edges:")

        for pair in gained_edges:
            print(
                f"  {pair[0]:28s} -- "
                f"{pair[1]:28s} "
                f"{geometric_distances[pair]:.6f} Å"
            )

    if lost_edges:
        print()
        print("Nominal edges outside range:")

        for first, second in lost_edges:
            value = distance(
                mapped[first]["xyz_A"],
                mapped[second]["xyz_A"],
            )

            print(
                f"  {first:28s} -- "
                f"{second:28s} "
                f"{value:.6f} Å"
            )

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
