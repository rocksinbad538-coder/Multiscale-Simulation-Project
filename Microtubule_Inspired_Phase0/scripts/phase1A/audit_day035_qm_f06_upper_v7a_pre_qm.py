#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import csv
import json
import math


ROOT = Path.cwd()

SOURCE_DIR = (
    ROOT
    / "runs/phase1A/day035_qm_f06_upper_v7a_formal_construction"
)

XYZ_PATH = (
    SOURCE_DIR
    / "QM_F06_UPPER_V7A_start.xyz"
)

MAP_PATH = (
    SOURCE_DIR
    / "QM_F06_UPPER_V7A_atom_role_provenance_map.csv"
)

EDGES_PATH = (
    SOURCE_DIR
    / "QM_F06_UPPER_V7A_nominal_edges.csv"
)

CAPS_PATH = (
    SOURCE_DIR
    / "QM_F06_UPPER_V7A_boundary_caps.csv"
)

CONSTRUCTION_REPORT = (
    SOURCE_DIR
    / "QM_F06_UPPER_V7A_CONSTRUCTION_REPORT.json"
)

OUTPUT_DIR = (
    ROOT
    / "runs/phase1A/day035_qm_f06_upper_v7a_pre_qm_audit"
)

REPORT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_PRE_QM_AUDIT.json"
)

BOND_AUDIT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_bond_audit.csv"
)

CONTACT_AUDIT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_contact_audit.csv"
)

VALENCE_AUDIT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_valence_audit.csv"
)


EXPECTED_COMPOSITION = {
    "B": 17,
    "N": 14,
    "H": 21,
}

EXPECTED_ATOM_COUNT = 52
EXPECTED_EDGE_COUNT = 57

MODIFIED_ATOMS = {
    "A:UPPER:13:1",
    "A:UPPER:14:0",
    "HCAPV7:UPPER:A13_1:A11_1",
    "HCAPV7:UPPER:A14_0:A13_M1",
    "H4:UPPER:0045:0",
    "A:UPPER:14:2",
}

NEW_ATOMS = {
    "A:UPPER:13:1",
    "A:UPPER:14:0",
    "HCAPV7:UPPER:A13_1:A11_1",
    "HCAPV7:UPPER:A14_0:A13_M1",
    "H4:UPPER:0045:0",
}

OBSOLETE_ATOM = "HCAPV2:UPPER:03"

REQUIRED_NEW_EDGES = {
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
    tuple(sorted((
        "A:UPPER:14:0",
        "H4:UPPER:0045:0",
    ))),
}

FORBIDDEN_FAILURE_PAIRS = {
    tuple(sorted((
        "A:UPPER:14:2",
        "P:1641",
    ))),
    tuple(sorted((
        "S:1710",
        "HCAPV2:UPPER:03",
    ))),
    tuple(sorted((
        "A:UPPER:14:2",
        "HCAPV2:UPPER:03",
    ))),
}

EXPECTED_DEGREE = {
    "B": 3,
    "N": 3,
    "H": 1,
}

# Nominal-bond acceptance windows.
BOND_WINDOWS = {
    ("B", "N"): (1.35, 1.75),
    ("B", "H"): (0.90, 1.35),
    ("H", "N"): (0.80, 1.25),
}

# Distances at or below these values are interpreted as
# geometric bonds when evaluating reconnectivity.
GEOMETRIC_BOND_MAXIMUM = {
    ("B", "N"): 1.90,
    ("B", "H"): 1.35,
    ("H", "N"): 1.25,
}

# Strengthened pre-QM near-contact margins.
MINIMUM_MODIFIED_REGION_MARGIN = {
    ("B", "N"): 0.08,
    ("B", "H"): 0.10,
    ("H", "N"): 0.10,
}

HARD_CONTACT_MINIMUM = {
    ("B", "B"): 1.20,
    ("B", "N"): 1.20,
    ("N", "N"): 1.20,
    ("B", "H"): 0.75,
    ("H", "N"): 0.70,
    ("H", "H"): 0.65,
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty source file: {path}"
        )


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def read_xyz(path: Path) -> list[dict]:
    lines = path.read_text(
        encoding="utf-8",
        errors="strict",
    ).splitlines()

    if len(lines) < 2:
        raise RuntimeError(
            f"Incomplete XYZ: {path}"
        )

    declared = int(lines[0].strip())
    atom_lines = lines[2:2 + declared]

    if len(atom_lines) != declared:
        raise RuntimeError(
            "XYZ atom-count mismatch: "
            f"declared={declared}, "
            f"parsed={len(atom_lines)}"
        )

    atoms = []

    for index, line in enumerate(atom_lines):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Malformed XYZ record {index}: {line}"
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


def canonical_elements(
    first: str,
    second: str,
) -> tuple[str, str]:
    return tuple(sorted((first, second)))


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


def geometric_bond(
    first_element: str,
    second_element: str,
    value: float,
) -> bool:
    element_pair = canonical_elements(
        first_element,
        second_element,
    )

    maximum = GEOMETRIC_BOND_MAXIMUM.get(
        element_pair
    )

    if maximum is None:
        return False

    return value <= maximum


def main() -> None:
    for path in (
        XYZ_PATH,
        MAP_PATH,
        EDGES_PATH,
        CAPS_PATH,
        CONSTRUCTION_REPORT,
    ):
        require_file(path)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    atoms = read_xyz(XYZ_PATH)
    map_rows = read_csv(MAP_PATH)
    edge_rows = read_csv(EDGES_PATH)
    cap_rows = read_csv(CAPS_PATH)

    construction = json.loads(
        CONSTRUCTION_REPORT.read_text(
            encoding="utf-8",
        )
    )

    map_rows.sort(
        key=lambda row: int(
            row["v7a_index_0based"]
        )
    )

    identity_order_failures = []

    if len(map_rows) == len(atoms):
        for index, (atom, row) in enumerate(
            zip(atoms, map_rows)
        ):
            if (
                int(row["v7a_index_0based"])
                != index
                or atom["element"]
                != row["element"]
            ):
                identity_order_failures.append({
                    "index_0based": index,
                    "atom_id": row["atom_id"],
                    "map_element": row["element"],
                    "xyz_element": atom["element"],
                })

    id_by_index = {
        index: row["atom_id"]
        for index, row in enumerate(map_rows)
    }

    index_by_id = {
        atom_id: index
        for index, atom_id in id_by_index.items()
    }

    element_by_id = {
        row["atom_id"]: row["element"]
        for row in map_rows
    }

    xyz_by_id = {
        id_by_index[index]: atom["xyz_A"]
        for index, atom in enumerate(atoms)
    }

    composition = Counter(
        atom["element"]
        for atom in atoms
    )

    nominal_edges = set()

    for row in edge_rows:
        nominal_edges.add(
            canonical_pair(
                row["first_atom"],
                row["second_atom"],
            )
        )

    adjacency = defaultdict(set)

    for first, second in nominal_edges:
        adjacency[first].add(second)
        adjacency[second].add(first)

    # Nominal-bond audit.
    bond_records = []
    bond_failures = []

    for first, second in sorted(nominal_edges):
        first_element = element_by_id[first]
        second_element = element_by_id[second]

        element_pair = canonical_elements(
            first_element,
            second_element,
        )

        value = distance(
            xyz_by_id[first],
            xyz_by_id[second],
        )

        window = BOND_WINDOWS.get(element_pair)

        if window is None:
            passed = False
            minimum = None
            maximum = None
            margin = None
        else:
            minimum, maximum = window
            passed = (
                minimum <= value <= maximum
            )
            margin = min(
                value - minimum,
                maximum - value,
            )

        record = {
            "first_atom": first,
            "first_element": first_element,
            "second_atom": second,
            "second_element": second_element,
            "bond_class": "-".join(element_pair),
            "distance_A": value,
            "minimum_A": minimum,
            "maximum_A": maximum,
            "margin_A": margin,
            "modified_region_edge": (
                first in MODIFIED_ATOMS
                or second in MODIFIED_ATOMS
            ),
            "pass": passed,
        }

        bond_records.append(record)

        if not passed:
            bond_failures.append(record)

    with BOND_AUDIT_PATH.open(
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

    # Geometric connectivity and global contacts.
    geometric_edges = set()
    geometric_distances = {}
    contact_records = []
    hard_contacts = []
    modified_near_contact_failures = []

    atom_ids = [
        row["atom_id"]
        for row in map_rows
    ]

    for first_index in range(len(atom_ids)):
        for second_index in range(
            first_index + 1,
            len(atom_ids),
        ):
            first = atom_ids[first_index]
            second = atom_ids[second_index]

            first_element = element_by_id[first]
            second_element = element_by_id[second]

            pair = canonical_pair(first, second)
            element_pair = canonical_elements(
                first_element,
                second_element,
            )

            value = distance(
                xyz_by_id[first],
                xyz_by_id[second],
            )

            nominal = pair in nominal_edges

            if geometric_bond(
                first_element,
                second_element,
                value,
            ):
                geometric_edges.add(pair)
                geometric_distances[pair] = value

            hard_minimum = HARD_CONTACT_MINIMUM.get(
                element_pair
            )

            hard_contact = (
                hard_minimum is not None
                and value < hard_minimum
            )

            modified_pair = (
                first in MODIFIED_ATOMS
                or second in MODIFIED_ATOMS
            )

            geometric_maximum = (
                GEOMETRIC_BOND_MAXIMUM.get(
                    element_pair
                )
            )

            if (
                not nominal
                and geometric_maximum is not None
            ):
                nonnominal_margin = (
                    value - geometric_maximum
                )
            else:
                nonnominal_margin = None

            required_margin = (
                MINIMUM_MODIFIED_REGION_MARGIN.get(
                    element_pair
                )
            )

            modified_near_contact_failure = (
                not nominal
                and modified_pair
                and required_margin is not None
                and nonnominal_margin is not None
                and nonnominal_margin
                < required_margin
            )

            record = {
                "first_atom": first,
                "first_element": first_element,
                "second_atom": second,
                "second_element": second_element,
                "pair_class": "-".join(element_pair),
                "distance_A": value,
                "nominal_edge": nominal,
                "geometric_edge": (
                    pair in geometric_edges
                ),
                "modified_region_pair": modified_pair,
                "geometric_bond_maximum_A": (
                    geometric_maximum
                ),
                "nonnominal_margin_A": (
                    nonnominal_margin
                ),
                "minimum_modified_region_margin_A": (
                    required_margin
                ),
                "modified_region_near_contact_failure": (
                    modified_near_contact_failure
                ),
                "hard_contact_minimum_A": (
                    hard_minimum
                ),
                "hard_contact": hard_contact,
            }

            contact_records.append(record)

            if hard_contact:
                hard_contacts.append(record)

            if modified_near_contact_failure:
                modified_near_contact_failures.append(
                    record
                )

    with CONTACT_AUDIT_PATH.open(
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

    # Exact nominal valence audit.
    valence_records = []
    degree_failures = []

    for atom_id in atom_ids:
        element = element_by_id[atom_id]
        observed = len(adjacency[atom_id])
        expected = EXPECTED_DEGREE[element]

        record = {
            "index_0based": index_by_id[atom_id],
            "atom_id": atom_id,
            "element": element,
            "degree": observed,
            "expected_degree": expected,
            "neighbors": "|".join(
                sorted(adjacency[atom_id])
            ),
            "pass": observed == expected,
        }

        valence_records.append(record)

        if observed != expected:
            degree_failures.append(record)

    with VALENCE_AUDIT_PATH.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                valence_records[0].keys()
            ),
        )
        writer.writeheader()
        writer.writerows(valence_records)

    # Geometric overcoordination.
    geometric_adjacency = defaultdict(set)

    for first, second in geometric_edges:
        geometric_adjacency[first].add(second)
        geometric_adjacency[second].add(first)

    overcoordinated = []

    for atom_id in atom_ids:
        observed = len(
            geometric_adjacency[atom_id]
        )
        expected = EXPECTED_DEGREE[
            element_by_id[atom_id]
        ]

        if observed > expected:
            overcoordinated.append({
                "atom_id": atom_id,
                "element": element_by_id[atom_id],
                "geometric_degree": observed,
                "expected_maximum": expected,
                "neighbors": sorted(
                    geometric_adjacency[atom_id]
                ),
            })

    # Graph connected components from nominal topology.
    visited = set()
    component_count = 0

    for atom_id in atom_ids:
        if atom_id in visited:
            continue

        component_count += 1
        stack = [atom_id]

        while stack:
            current = stack.pop()

            if current in visited:
                continue

            visited.add(current)
            stack.extend(
                adjacency[current] - visited
            )

    # Modified-region minimum margins.
    modified_bond_records = [
        row
        for row in bond_records
        if row["modified_region_edge"]
    ]

    minimum_modified_bond_margin = min(
        float(row["margin_A"])
        for row in modified_bond_records
        if row["margin_A"] is not None
    )

    eligible_modified_contacts = [
        row
        for row in contact_records
        if (
            row["modified_region_pair"]
            and not row["nominal_edge"]
            and row["nonnominal_margin_A"]
            is not None
        )
    ]

    minimum_modified_contact = min(
        eligible_modified_contacts,
        key=lambda row: float(
            row["nonnominal_margin_A"]
        ),
    )

    # New-cap and canonical-H geometry.
    required_cap_distances = {
        "HCAPV7:UPPER:A13_1:A11_1": 1.01,
        "HCAPV7:UPPER:A14_0:A13_M1": 1.19,
    }

    cap_failures = []
    cap_geometry = []

    for row in cap_rows:
        cap_id = row["cap_id"]
        owner = row["owner_atom"]

        value = distance(
            xyz_by_id[cap_id],
            xyz_by_id[owner],
        )

        target = required_cap_distances[
            cap_id
        ]

        passed = abs(value - target) <= 0.03

        record = {
            "cap_id": cap_id,
            "owner_atom": owner,
            "distance_A": value,
            "target_A": target,
            "absolute_deviation_A": abs(
                value - target
            ),
            "pass": passed,
        }

        cap_geometry.append(record)

        if not passed:
            cap_failures.append(record)

    h0045_distance = distance(
        xyz_by_id["H4:UPPER:0045:0"],
        xyz_by_id["A:UPPER:14:0"],
    )

    h0045_pass = (
        abs(h0045_distance - 1.19)
        <= 0.03
    )

    required_edges_present = (
        REQUIRED_NEW_EDGES <= nominal_edges
    )

    obsolete_atom_absent = (
        OBSOLETE_ATOM not in index_by_id
    )

    forbidden_pairs_absent = all(
        pair not in geometric_edges
        for pair in FORBIDDEN_FAILURE_PAIRS
    )

    specific_v6b_defect_eliminated = (
        obsolete_atom_absent
        and required_edges_present
        and forbidden_pairs_absent
        and canonical_pair(
            "A:UPPER:13:1",
            "A:UPPER:14:2",
        ) in nominal_edges
    )

    construction_authorized = (
        construction.get("decision")
        == (
            "QM_F06_UPPER_V7A_FORMALLY_CONSTRUCTED_"
            "GLOBAL_PRE_QM_AUDIT_REQUIRED"
        )
    )

    gates = {
        "construction_decision": (
            construction_authorized
        ),
        "atom_identity_and_order": (
            len(identity_order_failures) == 0
            and len(map_rows) == len(atoms)
        ),
        "atom_count": (
            len(atoms) == EXPECTED_ATOM_COUNT
        ),
        "composition": (
            dict(composition)
            == EXPECTED_COMPOSITION
        ),
        "nominal_edge_count": (
            len(nominal_edges)
            == EXPECTED_EDGE_COUNT
        ),
        "required_V7A_edges_present": (
            required_edges_present
        ),
        "obsolete_HCAPV2_absent": (
            obsolete_atom_absent
        ),
        "nominal_degree_exact": (
            len(degree_failures) == 0
        ),
        "single_connected_component": (
            component_count == 1
        ),
        "all_nominal_bonds_in_range": (
            len(bond_failures) == 0
        ),
        "no_geometric_reconnectivity": (
            len(gained_edges) == 0
            and len(lost_edges) == 0
        ),
        "no_geometric_overcoordination": (
            len(overcoordinated) == 0
        ),
        "no_hard_contacts": (
            len(hard_contacts) == 0
        ),
        "modified_region_near_contact_gate": (
            len(
                modified_near_contact_failures
            ) == 0
        ),
        "new_boundary_cap_geometry": (
            len(cap_failures) == 0
        ),
        "canonical_H0045_geometry": (
            h0045_pass
        ),
        "V6B_failure_mode_eliminated": (
            specific_v6b_defect_eliminated
        ),
    }

    passed = all(gates.values())

    decision = (
        "QM_F06_UPPER_V7A_PRE_QM_GATE_PASS_"
        "ORCA_INPUT_DESIGN_AUTHORIZED"
        if passed
        else
        "QM_F06_UPPER_V7A_PRE_QM_GATE_FAIL_"
        "STRUCTURAL_REVIEW_REQUIRED"
    )

    report = {
        "model": "QM_F06_UPPER_V7A",
        "decision": decision,
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "gates": gates,
        "summary": {
            "atom_count": len(atoms),
            "composition": dict(
                sorted(composition.items())
            ),
            "nominal_edge_count": len(
                nominal_edges
            ),
            "geometric_edge_count": len(
                geometric_edges
            ),
            "connected_component_count": (
                component_count
            ),
            "identity_order_failure_count": (
                len(identity_order_failures)
            ),
            "degree_failure_count": (
                len(degree_failures)
            ),
            "bond_failure_count": (
                len(bond_failures)
            ),
            "gained_geometric_edge_count": (
                len(gained_edges)
            ),
            "lost_nominal_edge_count": (
                len(lost_edges)
            ),
            "overcoordinated_atom_count": (
                len(overcoordinated)
            ),
            "hard_contact_count": (
                len(hard_contacts)
            ),
            "modified_region_near_contact_failure_count": (
                len(
                    modified_near_contact_failures
                )
            ),
            "minimum_modified_region_bond_margin_A": (
                minimum_modified_bond_margin
            ),
            "minimum_modified_region_nonnominal_margin_A": (
                minimum_modified_contact[
                    "nonnominal_margin_A"
                ]
            ),
            "limiting_modified_region_nonnominal_pair": (
                minimum_modified_contact[
                    "first_atom"
                ]
                + "--"
                + minimum_modified_contact[
                    "second_atom"
                ]
            ),
            "H0045_owner_distance_A": (
                h0045_distance
            ),
        },
        "gained_geometric_edges": [
            {
                "first_atom": first,
                "second_atom": second,
                "distance_A": (
                    geometric_distances[
                        (first, second)
                    ]
                ),
            }
            for first, second in gained_edges
        ],
        "lost_nominal_edges": [
            {
                "first_atom": first,
                "second_atom": second,
                "distance_A": distance(
                    xyz_by_id[first],
                    xyz_by_id[second],
                ),
            }
            for first, second in lost_edges
        ],
        "overcoordinated_atoms": (
            overcoordinated
        ),
        "modified_region_near_contact_failures": (
            modified_near_contact_failures
        ),
        "cap_geometry": cap_geometry,
        "V6B_failure_mode": {
            "obsolete_cap_absent": (
                obsolete_atom_absent
            ),
            "A14_2_A13_1_nominal_edge_present": (
                canonical_pair(
                    "A:UPPER:14:2",
                    "A:UPPER:13:1",
                )
                in nominal_edges
            ),
            "forbidden_failure_pairs_absent": (
                forbidden_pairs_absent
            ),
            "eliminated": (
                specific_v6b_defect_eliminated
            ),
        },
        "authorizations": {
            "ORCA_input_design_authorized": passed,
            "ORCA_execution_authorized": False,
            "RESP_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    REPORT_PATH.write_text(
        json.dumps(
            report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 112)
    print("QM_F06 UPPER V7-A FORMAL GLOBAL PRE-QM AUDIT")
    print("=" * 112)

    for name, value in gates.items():
        print(
            f"{name:52s}: "
            f"{'PASS' if value else 'FAIL'}"
        )

    print()
    print("Atoms:", len(atoms))
    print(
        "Composition:",
        dict(sorted(composition.items())),
    )
    print("Nominal edges:", len(nominal_edges))
    print(
        "Geometric edges:",
        len(geometric_edges),
    )
    print(
        "Connected components:",
        component_count,
    )
    print("Degree failures:", len(degree_failures))
    print("Bond failures:", len(bond_failures))
    print("Gained edges:", len(gained_edges))
    print("Lost edges:", len(lost_edges))
    print(
        "Overcoordinated atoms:",
        len(overcoordinated),
    )
    print("Hard contacts:", len(hard_contacts))
    print(
        "Modified-region near-contact failures:",
        len(modified_near_contact_failures),
    )

    print()
    print(
        "Minimum modified-region bond margin A:",
        minimum_modified_bond_margin,
    )
    print(
        "Minimum modified-region nonnominal margin A:",
        minimum_modified_contact[
            "nonnominal_margin_A"
        ],
    )
    print(
        "Limiting modified-region nonnominal pair:",
        minimum_modified_contact[
            "first_atom"
        ]
        + "--"
        + minimum_modified_contact[
            "second_atom"
        ],
    )

    print()
    print("New cap geometry:")

    for record in cap_geometry:
        print(
            f"  {record['cap_id']:42s} "
            f"owner={record['owner_atom']:24s} "
            f"distance={record['distance_A']:.6f} Å | "
            f"{'PASS' if record['pass'] else 'FAIL'}"
        )

    print(
        "  H4:UPPER:0045:0"
        f"{'':25s} "
        "owner=A:UPPER:14:0"
        f"{'':8s} "
        f"distance={h0045_distance:.6f} Å | "
        f"{'PASS' if h0045_pass else 'FAIL'}"
    )

    if gained_edges:
        print()
        print("Nonnominal geometric edges:")

        for first, second in gained_edges:
            print(
                f"  {first:32s} -- "
                f"{second:32s} "
                f"{geometric_distances[(first, second)]:.6f} Å"
            )

    if lost_edges:
        print()
        print("Nominal edges outside geometric range:")

        for first, second in lost_edges:
            print(
                f"  {first:32s} -- "
                f"{second:32s} "
                f"{distance(xyz_by_id[first], xyz_by_id[second]):.6f} Å"
            )

    if modified_near_contact_failures:
        print()
        print("Modified-region near-contact failures:")

        for record in (
            modified_near_contact_failures
        ):
            print(
                f"  {record['first_atom']:32s} -- "
                f"{record['second_atom']:32s} "
                f"distance={float(record['distance_A']):.6f} Å | "
                f"margin={float(record['nonnominal_margin_A']):.6f} Å"
            )

    print()
    print("Decision:", decision)
    print("Report:", REPORT_PATH)
    print("Bond audit:", BOND_AUDIT_PATH)
    print("Contact audit:", CONTACT_AUDIT_PATH)
    print("Valence audit:", VALENCE_AUDIT_PATH)

    print()
    print(
        "ORCA input design authorized:",
        passed,
    )
    print("ORCA execution authorized: False")
    print("RESP authorized: False")
    print("Force-field adoption authorized: False")
    print("MD authorized: False")

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
