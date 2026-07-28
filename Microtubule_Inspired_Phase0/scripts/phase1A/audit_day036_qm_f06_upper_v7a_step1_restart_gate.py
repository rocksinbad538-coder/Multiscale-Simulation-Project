#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
import csv
import json
import math


ROOT = Path(__file__).resolve().parents[2]

EXEC_PARENT = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_executions"
)

LATEST_FILE = (
    EXEC_PARENT
    / "LATEST_V7A_EXECUTION.txt"
)

PREQM_DIR = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_pre_qm_audit"
)

OUTPUT_DIR = (
    ROOT
    / "runs/phase1A/"
    / "day036_qm_f06_upper_v7a_step1_restart_gate"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

OUTPUT_REPORT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_STEP1_RESTART_GATE.json"
)

OUTPUT_BONDS = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_STEP1_bond_audit.csv"
)

OUTPUT_CONTACTS = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_STEP1_contact_audit.csv"
)

OUTPUT_VALENCE = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_STEP1_valence_audit.csv"
)

OUTPUT_DISPLACEMENTS = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_STEP1_displacements.csv"
)


FIXED_ATOM_RESTART_TOLERANCE_A = 5.0e-5

REQUIRED_V7A_EDGES = {
    frozenset(
        (
            "A:UPPER:13:1",
            "A:UPPER:14:2",
        )
    ),
    frozenset(
        (
            "A:UPPER:13:1",
            "A:UPPER:14:0",
        )
    ),
}

V6B_FORBIDDEN_RECONNECTIVITY = {
    frozenset(
        (
            "A:UPPER:14:2",
            "P:1641",
        )
    ),
}

REQUIRED_V7A_CAP_BONDS = {
    frozenset(
        (
            "A:UPPER:13:1",
            "HCAPV7:UPPER:A13_1:A11_1",
        )
    ),
    frozenset(
        (
            "A:UPPER:14:0",
            "HCAPV7:UPPER:A14_0:A13_M1",
        )
    ),
    frozenset(
        (
            "A:UPPER:14:0",
            "H4:UPPER:0045:0",
        )
    ),
}


def read_xyz(path: Path) -> list[dict]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    atom_count = int(lines[0].strip())

    atoms = []

    for index, line in enumerate(
        lines[2:2 + atom_count]
    ):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Incomplete XYZ record at index {index}"
            )

        atoms.append({
            "index_0based": index,
            "element": fields[0],
            "xyz_A": tuple(
                float(value)
                for value in fields[1:4]
            ),
        })

    if len(atoms) != atom_count:
        raise RuntimeError(
            f"Incomplete XYZ: {path}"
        )

    return atoms


def read_csv(path: Path) -> list[dict]:
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path,
    records: list[dict],
) -> None:
    if not records:
        path.write_text(
            "",
            encoding="utf-8",
        )
        return

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(records[0]),
        )

        writer.writeheader()
        writer.writerows(records)


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


def parse_bool(value: str) -> bool:
    return (
        value.strip().lower()
        == "true"
    )


def parse_optional_float(
    value: str,
) -> float | None:
    stripped = value.strip()

    if not stripped:
        return None

    return float(stripped)


def main() -> None:
    execution_relative = (
        LATEST_FILE
        .read_text(encoding="utf-8")
        .strip()
    )

    execution_dir = (
        ROOT
        / execution_relative
    )

    start_xyz_path = (
        execution_dir
        / "v7a_start.xyz"
    )

    current_xyz_path = (
        execution_dir
        / "v7a.xyz"
    )

    map_path = (
        execution_dir
        / "QM_F06_UPPER_V7A_constraint_map.csv"
    )

    preqm_bond_path = (
        PREQM_DIR
        / "QM_F06_UPPER_V7A_bond_audit.csv"
    )

    preqm_contact_path = (
        PREQM_DIR
        / "QM_F06_UPPER_V7A_contact_audit.csv"
    )

    preqm_valence_path = (
        PREQM_DIR
        / "QM_F06_UPPER_V7A_valence_audit.csv"
    )

    required_paths = (
        start_xyz_path,
        current_xyz_path,
        map_path,
        preqm_bond_path,
        preqm_contact_path,
        preqm_valence_path,
    )

    for path in required_paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(
                f"Missing required input: {path}"
            )

    start_atoms = read_xyz(start_xyz_path)
    current_atoms = read_xyz(current_xyz_path)

    map_rows = read_csv(map_path)

    map_rows.sort(
        key=lambda row: int(
            row["v7a_index_0based"]
        )
    )

    if not (
        len(start_atoms)
        == len(current_atoms)
        == len(map_rows)
        == 52
    ):
        raise RuntimeError(
            "Atom-count mismatch"
        )

    atom_ids = []
    elements = {}
    start_coordinates = {}
    current_coordinates = {}
    fixed_status = {}

    identity_failures = []

    for index, row in enumerate(map_rows):
        atom_id = row["atom_id"]
        element = row["element"]

        atom_ids.append(atom_id)
        elements[atom_id] = element
        fixed_status[atom_id] = parse_bool(
            row["v7a_fixed"]
        )

        if not (
            start_atoms[index]["element"]
            == current_atoms[index]["element"]
            == element
        ):
            identity_failures.append({
                "index_0based": index,
                "atom_id": atom_id,
                "map_element": element,
                "start_element": (
                    start_atoms[index]["element"]
                ),
                "current_element": (
                    current_atoms[index]["element"]
                ),
            })

        start_coordinates[atom_id] = (
            start_atoms[index]["xyz_A"]
        )

        current_coordinates[atom_id] = (
            current_atoms[index]["xyz_A"]
        )

    composition = Counter(elements.values())

    displacement_records = []

    for index, atom_id in enumerate(atom_ids):
        shift = distance(
            start_coordinates[atom_id],
            current_coordinates[atom_id],
        )

        fixed = fixed_status[atom_id]

        displacement_records.append({
            "index_0based": index,
            "atom_id": atom_id,
            "element": elements[atom_id],
            "v7a_fixed": fixed,
            "displacement_A": shift,
            "fixed_tolerance_A": (
                FIXED_ATOM_RESTART_TOLERANCE_A
                if fixed
                else ""
            ),
            "fixed_displacement_pass": (
                shift
                <= FIXED_ATOM_RESTART_TOLERANCE_A
                if fixed
                else True
            ),
        })

    maximum_displacement = max(
        float(row["displacement_A"])
        for row in displacement_records
    )

    maximum_fixed_displacement = max(
        float(row["displacement_A"])
        for row in displacement_records
        if row["v7a_fixed"]
    )

    preqm_bonds = read_csv(
        preqm_bond_path
    )

    nominal_edges = {
        frozenset(
            (
                row["first_atom"],
                row["second_atom"],
            )
        )
        for row in preqm_bonds
    }

    bond_records = []
    bond_failures = []

    for row in preqm_bonds:
        first = row["first_atom"]
        second = row["second_atom"]

        value = distance(
            current_coordinates[first],
            current_coordinates[second],
        )

        minimum = float(
            row["minimum_A"]
        )

        maximum = float(
            row["maximum_A"]
        )

        margin = min(
            value - minimum,
            maximum - value,
        )

        passed = (
            minimum
            <= value
            <= maximum
        )

        record = {
            "first_atom": first,
            "first_element": row[
                "first_element"
            ],
            "second_atom": second,
            "second_element": row[
                "second_element"
            ],
            "bond_class": row["bond_class"],
            "distance_A": value,
            "minimum_A": minimum,
            "maximum_A": maximum,
            "margin_A": margin,
            "modified_region_edge": row[
                "modified_region_edge"
            ],
            "pass": passed,
        }

        bond_records.append(record)

        if not passed:
            bond_failures.append(record)

    preqm_contacts = read_csv(
        preqm_contact_path
    )

    contact_records = []
    geometric_edges = set()
    gained_edges = []
    hard_contacts = []
    modified_near_contact_failures = []

    for row in preqm_contacts:
        first = row["first_atom"]
        second = row["second_atom"]
        pair = frozenset((first, second))

        value = distance(
            current_coordinates[first],
            current_coordinates[second],
        )

        nominal = (
            pair in nominal_edges
        )

        geometric_maximum = (
            parse_optional_float(
                row[
                    "geometric_bond_maximum_A"
                ]
            )
        )

        geometric = (
            geometric_maximum is not None
            and value <= geometric_maximum
        )

        if geometric:
            geometric_edges.add(pair)

        gained = (
            geometric
            and not nominal
        )

        nonnominal_margin = (
            value - geometric_maximum
            if (
                geometric_maximum is not None
                and not nominal
            )
            else None
        )

        modified_region_pair = parse_bool(
            row["modified_region_pair"]
        )

        required_margin = (
            parse_optional_float(
                row[
                    "minimum_modified_region_margin_A"
                ]
            )
        )

        near_contact_failure = (
            modified_region_pair
            and not nominal
            and geometric_maximum is not None
            and required_margin is not None
            and nonnominal_margin is not None
            and nonnominal_margin
            < required_margin
        )

        hard_minimum = float(
            row["hard_contact_minimum_A"]
        )

        hard_contact = (
            value < hard_minimum
        )

        record = {
            "first_atom": first,
            "first_element": row[
                "first_element"
            ],
            "second_atom": second,
            "second_element": row[
                "second_element"
            ],
            "pair_class": row["pair_class"],
            "distance_A": value,
            "nominal_edge": nominal,
            "geometric_edge": geometric,
            "gained_geometric_edge": gained,
            "modified_region_pair": (
                modified_region_pair
            ),
            "geometric_bond_maximum_A": (
                geometric_maximum
                if geometric_maximum is not None
                else ""
            ),
            "nonnominal_margin_A": (
                nonnominal_margin
                if nonnominal_margin is not None
                else ""
            ),
            "minimum_modified_region_margin_A": (
                required_margin
                if required_margin is not None
                else ""
            ),
            "modified_region_near_contact_failure": (
                near_contact_failure
            ),
            "hard_contact_minimum_A": hard_minimum,
            "hard_contact": hard_contact,
        }

        contact_records.append(record)

        if gained:
            gained_edges.append(record)

        if hard_contact:
            hard_contacts.append(record)

        if near_contact_failure:
            modified_near_contact_failures.append(
                record
            )

    lost_edges = sorted(
        nominal_edges - geometric_edges,
        key=lambda pair: sorted(pair),
    )

    adjacency = defaultdict(set)

    for pair in geometric_edges:
        first, second = tuple(pair)
        adjacency[first].add(second)
        adjacency[second].add(first)

    preqm_valence = read_csv(
        preqm_valence_path
    )

    valence_records = []
    degree_failures = []
    overcoordinated_atoms = []

    expected_degree = {
        row["atom_id"]: int(
            row["expected_degree"]
        )
        for row in preqm_valence
    }

    for index, atom_id in enumerate(atom_ids):
        degree = len(
            adjacency[atom_id]
        )

        expected = expected_degree[
            atom_id
        ]

        passed = (
            degree == expected
        )

        record = {
            "index_0based": index,
            "atom_id": atom_id,
            "element": elements[atom_id],
            "degree": degree,
            "expected_degree": expected,
            "neighbors": "|".join(
                sorted(adjacency[atom_id])
            ),
            "pass": passed,
        }

        valence_records.append(record)

        if not passed:
            degree_failures.append(record)

        if degree > expected:
            overcoordinated_atoms.append(
                record
            )

    visited = set()
    components = []

    for atom_id in atom_ids:
        if atom_id in visited:
            continue

        stack = [atom_id]
        component = set()

        while stack:
            current = stack.pop()

            if current in component:
                continue

            component.add(current)
            visited.add(current)

            stack.extend(
                adjacency[current]
                - component
            )

        components.append(component)

    required_edges_present = (
        REQUIRED_V7A_EDGES
        <= geometric_edges
    )

    required_cap_bonds_present = (
        REQUIRED_V7A_CAP_BONDS
        <= geometric_edges
    )

    forbidden_reconnectivity_absent = (
        not any(
            pair in geometric_edges
            for pair
            in V6B_FORBIDDEN_RECONNECTIVITY
        )
    )

    obsolete_cap_absent = (
        "HCAPV2:UPPER:03"
        not in atom_ids
    )

    fixed_gate = all(
        row["fixed_displacement_pass"]
        for row in displacement_records
        if row["v7a_fixed"]
    )

    gates = {
        "atom_identity_and_order": (
            len(identity_failures) == 0
        ),
        "atom_count_52": (
            len(atom_ids) == 52
        ),
        "composition_B17_N14_H21": (
            composition
            == Counter({
                "B": 17,
                "N": 14,
                "H": 21,
            })
        ),
        "fixed_atoms_preserved": fixed_gate,
        "nominal_edge_count_57": (
            len(nominal_edges) == 57
        ),
        "all_nominal_bonds_in_range": (
            len(bond_failures) == 0
        ),
        "no_lost_nominal_edges": (
            len(lost_edges) == 0
        ),
        "no_geometric_reconnectivity": (
            len(gained_edges) == 0
        ),
        "nominal_degree_exact": (
            len(degree_failures) == 0
        ),
        "no_geometric_overcoordination": (
            len(overcoordinated_atoms) == 0
        ),
        "single_connected_component": (
            len(components) == 1
        ),
        "no_hard_contacts": (
            len(hard_contacts) == 0
        ),
        "modified_region_near_contact_gate": (
            len(
                modified_near_contact_failures
            )
            == 0
        ),
        "required_V7A_edges_present": (
            required_edges_present
        ),
        "required_V7A_caps_present": (
            required_cap_bonds_present
        ),
        "obsolete_HCAPV2_absent": (
            obsolete_cap_absent
        ),
        "V6B_forbidden_reconnectivity_absent": (
            forbidden_reconnectivity_absent
        ),
    }

    passed = all(gates.values())

    minimum_bond_margin = min(
        float(row["margin_A"])
        for row in bond_records
    )

    nonnominal_margins = [
        float(row["nonnominal_margin_A"])
        for row in contact_records
        if (
            row["modified_region_pair"]
            and not row["nominal_edge"]
            and row["nonnominal_margin_A"] != ""
        )
    ]

    minimum_modified_nonnominal_margin = (
        min(nonnominal_margins)
        if nonnominal_margins
        else None
    )

    decision = (
        "QM_F06_UPPER_V7A_STEP1_"
        "STRUCTURAL_GATE_PASS_"
        "SCF_RESTART_DESIGN_AUTHORIZED"
        if passed
        else
        "QM_F06_UPPER_V7A_STEP1_"
        "STRUCTURAL_GATE_FAIL_"
        "GEOMETRY_OR_TOPOLOGY_REVIEW_REQUIRED"
    )

    report = {
        "model": "QM_F06_UPPER_V7A",
        "audit_stage": (
            "POST_FIRST_ACCEPTED_GEOMETRY_STEP_"
            "AFTER_LEANSCF_ERROR_TERMINATION"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "execution_directory": str(
            execution_dir
        ),
        "audited_geometry": str(
            current_xyz_path
        ),
        "decision": decision,
        "gates": gates,
        "summary": {
            "atom_count": len(atom_ids),
            "composition": dict(
                sorted(composition.items())
            ),
            "nominal_edge_count": (
                len(nominal_edges)
            ),
            "geometric_edge_count": (
                len(geometric_edges)
            ),
            "connected_component_count": (
                len(components)
            ),
            "maximum_atom_displacement_A": (
                maximum_displacement
            ),
            "maximum_fixed_atom_displacement_A": (
                maximum_fixed_displacement
            ),
            "fixed_atom_restart_tolerance_A": (
                FIXED_ATOM_RESTART_TOLERANCE_A
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
            "degree_failure_count": (
                len(degree_failures)
            ),
            "overcoordinated_atom_count": (
                len(overcoordinated_atoms)
            ),
            "hard_contact_count": (
                len(hard_contacts)
            ),
            "modified_region_near_contact_failure_count": (
                len(
                    modified_near_contact_failures
                )
            ),
            "minimum_bond_margin_A": (
                minimum_bond_margin
            ),
            "minimum_modified_region_nonnominal_margin_A": (
                minimum_modified_nonnominal_margin
            ),
        },
        "bond_failures": bond_failures,
        "gained_geometric_edges": gained_edges,
        "lost_nominal_edges": [
            sorted(pair)
            for pair in lost_edges
        ],
        "degree_failures": degree_failures,
        "overcoordinated_atoms": (
            overcoordinated_atoms
        ),
        "hard_contacts": hard_contacts,
        "modified_region_near_contact_failures": (
            modified_near_contact_failures
        ),
        "authorizations": {
            "SCF_restart_input_design_authorized": (
                passed
            ),
            "ORCA_restart_execution_authorized": (
                False
            ),
            "RESP_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    write_csv(
        OUTPUT_BONDS,
        bond_records,
    )

    write_csv(
        OUTPUT_CONTACTS,
        contact_records,
    )

    write_csv(
        OUTPUT_VALENCE,
        valence_records,
    )

    write_csv(
        OUTPUT_DISPLACEMENTS,
        displacement_records,
    )

    OUTPUT_REPORT.write_text(
        json.dumps(
            report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 112)
    print(
        "QM_F06 UPPER V7-A STEP-1 "
        "STRUCTURAL RESTART GATE"
    )
    print("=" * 112)

    for name, value in gates.items():
        print(
            f"{name:58s}: "
            f"{'PASS' if value else 'FAIL'}"
        )

    print()
    print("Atoms:", len(atom_ids))
    print(
        "Composition:",
        dict(sorted(composition.items())),
    )
    print(
        "Maximum atom displacement A:",
        maximum_displacement,
    )
    print(
        "Maximum fixed displacement A:",
        maximum_fixed_displacement,
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
        "Degree failures:",
        len(degree_failures),
    )
    print(
        "Overcoordinated atoms:",
        len(overcoordinated_atoms),
    )
    print(
        "Hard contacts:",
        len(hard_contacts),
    )
    print(
        "Modified-region near-contact failures:",
        len(modified_near_contact_failures),
    )
    print(
        "Minimum bond margin A:",
        minimum_bond_margin,
    )
    print(
        "Minimum modified-region "
        "nonnominal margin A:",
        minimum_modified_nonnominal_margin,
    )

    if bond_failures:
        print()
        print("Bond failures:")
        for row in bond_failures:
            print(
                f"  {row['first_atom']} -- "
                f"{row['second_atom']} | "
                f"{float(row['distance_A']):.6f} Å"
            )

    if gained_edges:
        print()
        print("Gained geometric edges:")
        for row in gained_edges:
            print(
                f"  {row['first_atom']} -- "
                f"{row['second_atom']} | "
                f"{float(row['distance_A']):.6f} Å"
            )

    if lost_edges:
        print()
        print("Lost nominal edges:")
        for pair in lost_edges:
            first, second = sorted(pair)
            print(
                f"  {first} -- {second}"
            )

    print()
    print("Decision:", decision)
    print("Report:", OUTPUT_REPORT)
    print()
    print(
        "SCF restart input design authorized:",
        passed,
    )
    print("ORCA restart execution authorized: False")
    print("RESP authorized: False")
    print("MD authorized: False")

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
