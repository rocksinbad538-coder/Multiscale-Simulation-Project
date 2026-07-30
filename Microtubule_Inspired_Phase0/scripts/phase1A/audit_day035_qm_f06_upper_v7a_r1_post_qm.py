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
    / "day035_qm_f06_upper_v7a_r1_executions"
)

LATEST_FILE = (
    EXEC_PARENT
    / "LATEST_V7A_R1_EXECUTION.txt"
)

READINESS_REPORT = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_r1_post_qm/"
    / "QM_F06_UPPER_V7A_R1_POST_QM_READINESS.json"
)

OUTPUT_DIR = (
    ROOT
    / "runs/phase1A/"
    / "day035_qm_f06_upper_v7a_r1_post_qm"
)

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

OUTPUT_REPORT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_POST_QM_AUDIT.json"
)

OUTPUT_FINAL_XYZ = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_FINAL.xyz"
)

OUTPUT_BONDS = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_post_qm_bond_audit.csv"
)

OUTPUT_CONTACTS = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_post_qm_contact_audit.csv"
)

OUTPUT_VALENCE = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_post_qm_valence_audit.csv"
)

OUTPUT_DISPLACEMENTS = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_post_qm_displacements.csv"
)

OUTPUT_TRAJECTORY = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_post_qm_trajectory_audit.csv"
)


EXPECTED_READINESS_DECISION = (
    "QM_F06_UPPER_V7A_R1_"
    "POST_QM_READINESS_GATE_PASS_"
    "STRUCTURAL_AUDIT_AUTHORIZED"
)

EXPECTED_ATOM_COUNT = 52
EXPECTED_EDGE_COUNT = 57

EXPECTED_COMPOSITION = {
    "B": 17,
    "N": 14,
    "H": 21,
}

EXPECTED_DEGREE = {
    "B": 3,
    "N": 3,
    "H": 1,
}

FIXED_ATOM_TOLERANCE_A = 5.0e-5

BOND_WINDOWS = {
    ("B", "N"): (1.35, 1.75),
    ("B", "H"): (0.90, 1.35),
    ("H", "N"): (0.80, 1.25),
}

GEOMETRIC_BOND_MAXIMUM = {
    ("B", "N"): 1.90,
    ("B", "H"): 1.35,
    ("H", "N"): 1.25,
}

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

MODIFIED_ATOMS = {
    "A:UPPER:13:1",
    "A:UPPER:14:0",
    "A:UPPER:14:2",
    "HCAPV7:UPPER:A13_1:A11_1",
    "HCAPV7:UPPER:A14_0:A13_M1",
    "H4:UPPER:0045:0",
}

REQUIRED_V7A_EDGES = {
    frozenset((
        "A:UPPER:13:1",
        "A:UPPER:14:2",
    )),
    frozenset((
        "A:UPPER:13:1",
        "A:UPPER:14:0",
    )),
    frozenset((
        "A:UPPER:13:1",
        "HCAPV7:UPPER:A13_1:A11_1",
    )),
    frozenset((
        "A:UPPER:14:0",
        "HCAPV7:UPPER:A14_0:A13_M1",
    )),
    frozenset((
        "A:UPPER:14:0",
        "H4:UPPER:0045:0",
    )),
}

FORBIDDEN_V6B_EDGES = {
    frozenset((
        "A:UPPER:14:2",
        "P:1641",
    )),
}

OBSOLETE_ATOM = "HCAPV2:UPPER:03"


def parse_bool(value: str) -> bool:
    return value.strip().lower() == "true"


def pair_class(
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


def read_csv(path: Path) -> list[dict]:
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def write_csv(
    path: Path,
    rows: list[dict],
) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
        )
        writer.writeheader()
        writer.writerows(rows)


def read_xyz(path: Path) -> list[dict]:
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
                f"Incomplete XYZ record {index}: {path}"
            )

        atoms.append({
            "index_0based": index,
            "element": fields[0],
            "xyz_A": tuple(
                float(value)
                for value in fields[1:4]
            ),
        })

    if len(atoms) != count:
        raise RuntimeError(
            f"Incomplete XYZ: {path}"
        )

    return atoms


def read_trajectory(path: Path) -> list[dict]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    frames = []
    cursor = 0

    while cursor < len(lines):
        try:
            atom_count = int(
                lines[cursor].strip()
            )
        except (ValueError, IndexError):
            break

        end = cursor + atom_count + 2

        if end > len(lines):
            break

        atoms = []

        for index, line in enumerate(
            lines[cursor + 2:end]
        ):
            fields = line.split()

            if len(fields) < 4:
                atoms = []
                break

            atoms.append({
                "index_0based": index,
                "element": fields[0],
                "xyz_A": tuple(
                    float(value)
                    for value in fields[1:4]
                ),
            })

        if len(atoms) != atom_count:
            break

        frames.append({
            "frame_index_0based": len(frames),
            "comment": lines[cursor + 1],
            "atoms": atoms,
        })

        cursor = end

    return frames


def extract_energy(comment: str) -> float | None:
    fields = comment.split()

    for index, token in enumerate(fields):
        if token == "E" and index + 1 < len(fields):
            try:
                return float(fields[index + 1])
            except ValueError:
                return None

    return None


def audit_geometry(
    coordinates: dict[
        str,
        tuple[float, float, float]
    ],
    atom_ids: list[str],
    elements: dict[str, str],
    nominal_edges: set[frozenset[str]],
    edge_rows: list[dict],
) -> dict:
    geometric_edges: set[frozenset[str]] = set()
    contact_rows = []
    hard_contacts = []
    modified_margin_failures = []

    for first_index, first in enumerate(atom_ids):
        for second in atom_ids[first_index + 1:]:
            pair = frozenset((first, second))

            element_pair = pair_class(
                elements[first],
                elements[second],
            )

            value = distance(
                coordinates[first],
                coordinates[second],
            )

            nominal = pair in nominal_edges

            geometric_maximum = (
                GEOMETRIC_BOND_MAXIMUM.get(
                    element_pair
                )
            )

            geometric = (
                geometric_maximum is not None
                and value <= geometric_maximum
            )

            if geometric:
                geometric_edges.add(pair)

            hard_minimum = (
                HARD_CONTACT_MINIMUM.get(
                    element_pair
                )
            )

            hard_contact = (
                hard_minimum is not None
                and value < hard_minimum
            )

            modified_pair = (
                first in MODIFIED_ATOMS
                or second in MODIFIED_ATOMS
            )

            required_margin = (
                MINIMUM_MODIFIED_REGION_MARGIN.get(
                    element_pair
                )
            )

            nonnominal_margin = (
                value - geometric_maximum
                if (
                    not nominal
                    and geometric_maximum is not None
                )
                else None
            )

            modified_margin_failure = (
                modified_pair
                and not nominal
                and required_margin is not None
                and nonnominal_margin is not None
                and nonnominal_margin
                < required_margin
            )

            record = {
                "first_atom": first,
                "first_element": elements[first],
                "second_atom": second,
                "second_element": elements[second],
                "pair_class": "-".join(
                    element_pair
                ),
                "distance_A": value,
                "nominal_edge": nominal,
                "geometric_edge": geometric,
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
                "modified_region_pair": (
                    modified_pair
                ),
                "minimum_modified_region_margin_A": (
                    required_margin
                    if required_margin is not None
                    else ""
                ),
                "modified_region_near_contact_failure": (
                    modified_margin_failure
                ),
                "hard_contact_minimum_A": (
                    hard_minimum
                    if hard_minimum is not None
                    else ""
                ),
                "hard_contact": hard_contact,
            }

            contact_rows.append(record)

            if hard_contact:
                hard_contacts.append(record)

            if modified_margin_failure:
                modified_margin_failures.append(
                    record
                )

    bond_rows = []
    bond_failures = []

    for row in edge_rows:
        first = row["first_atom"]
        second = row["second_atom"]

        element_pair = pair_class(
            row["first_element"],
            row["second_element"],
        )

        minimum, maximum = (
            BOND_WINDOWS[element_pair]
        )

        value = distance(
            coordinates[first],
            coordinates[second],
        )

        margin = min(
            value - minimum,
            maximum - value,
        )

        passed = (
            minimum <= value <= maximum
        )

        record = {
            "first_atom": first,
            "first_element": (
                row["first_element"]
            ),
            "second_atom": second,
            "second_element": (
                row["second_element"]
            ),
            "edge_type": row["edge_type"],
            "provenance": row["provenance"],
            "distance_A": value,
            "minimum_A": minimum,
            "maximum_A": maximum,
            "margin_A": margin,
            "pass": passed,
        }

        bond_rows.append(record)

        if not passed:
            bond_failures.append(record)

    lost_edges = (
        nominal_edges - geometric_edges
    )

    gained_edges = (
        geometric_edges - nominal_edges
    )

    adjacency = defaultdict(set)

    for pair in geometric_edges:
        first, second = tuple(pair)
        adjacency[first].add(second)
        adjacency[second].add(first)

    valence_rows = []
    degree_failures = []
    overcoordinated_atoms = []

    for index, atom_id in enumerate(atom_ids):
        degree = len(adjacency[atom_id])
        expected = EXPECTED_DEGREE[
            elements[atom_id]
        ]

        passed = degree == expected

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

        valence_rows.append(record)

        if not passed:
            degree_failures.append(record)

        if degree > expected:
            overcoordinated_atoms.append(
                record
            )

    components = []
    visited = set()

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
                adjacency[current] - component
            )

        components.append(component)

    nonnominal_margins = [
        float(row["nonnominal_margin_A"])
        for row in contact_rows
        if row["nonnominal_margin_A"] != ""
    ]

    return {
        "geometric_edges": geometric_edges,
        "bond_rows": bond_rows,
        "contact_rows": contact_rows,
        "valence_rows": valence_rows,
        "bond_failures": bond_failures,
        "lost_edges": lost_edges,
        "gained_edges": gained_edges,
        "degree_failures": degree_failures,
        "overcoordinated_atoms": (
            overcoordinated_atoms
        ),
        "components": components,
        "hard_contacts": hard_contacts,
        "modified_margin_failures": (
            modified_margin_failures
        ),
        "minimum_bond_margin_A": min(
            float(row["margin_A"])
            for row in bond_rows
        ),
        "minimum_nonnominal_margin_A": (
            min(nonnominal_margins)
            if nonnominal_margins
            else None
        ),
    }


def main() -> None:
    if not READINESS_REPORT.is_file():
        raise RuntimeError(
            f"Missing readiness report: "
            f"{READINESS_REPORT}"
        )

    readiness = json.loads(
        READINESS_REPORT.read_text(
            encoding="utf-8"
        )
    )

    readiness_decision = readiness.get(
        "decision"
    )

    readiness_authorized = (
        readiness.get(
            "authorizations",
            {},
        ).get(
            "post_QM_structural_audit_authorized"
        )
        is True
    )

    if (
        readiness_decision
        != EXPECTED_READINESS_DECISION
        or not readiness_authorized
    ):
        raise RuntimeError(
            "Post-QM readiness gate has not passed: "
            f"decision={readiness_decision}; "
            "post_QM_structural_audit_authorized="
            f"{readiness_authorized}"
        )

    if not LATEST_FILE.is_file():
        raise RuntimeError(
            f"Missing execution pointer: {LATEST_FILE}"
        )

    execution_relative = (
        LATEST_FILE
        .read_text(encoding="utf-8")
        .strip()
    )

    if not execution_relative:
        raise RuntimeError(
            f"Empty execution pointer: {LATEST_FILE}"
        )

    execution_dir = (
        ROOT
        / execution_relative
    )

    readiness_execution_directory = (
        readiness.get(
            "execution_directory"
        )
    )

    if (
        readiness_execution_directory
        != str(execution_dir)
    ):
        raise RuntimeError(
            "Readiness report and current execution pointer "
            "refer to different executions: "
            f"readiness={readiness_execution_directory}; "
            f"current={execution_dir}"
        )

    start_path = (
        execution_dir
        / "v7a_r1_start.xyz"
    )

    final_path = (
        execution_dir
        / "v7a_r1.xyz"
    )

    trajectory_path = (
        execution_dir
        / "v7a_r1_trj.xyz"
    )

    map_path = (
        execution_dir
        / "QM_F06_UPPER_V7A_R1_constraint_map.csv"
    )

    edges_path = (
        execution_dir
        / "QM_F06_UPPER_V7A_R1_nominal_edges.csv"
    )

    for path in (
        start_path,
        final_path,
        trajectory_path,
        map_path,
        edges_path,
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(
                f"Missing post-QM input: {path}"
            )

    map_rows = read_csv(map_path)

    map_rows.sort(
        key=lambda row: int(
            row["v7a_index_0based"]
        )
    )

    edge_rows = read_csv(edges_path)

    start_atoms = read_xyz(start_path)
    final_atoms = read_xyz(final_path)
    frames = read_trajectory(
        trajectory_path
    )

    atom_ids = [
        row["atom_id"]
        for row in map_rows
    ]

    elements = {
        row["atom_id"]: row["element"]
        for row in map_rows
    }

    fixed_status = {
        row["atom_id"]: parse_bool(
            row["v7a_fixed"]
        )
        for row in map_rows
    }

    identity_failures = []

    for index, row in enumerate(map_rows):
        if not (
            start_atoms[index]["element"]
            == final_atoms[index]["element"]
            == row["element"]
        ):
            identity_failures.append({
                "index_0based": index,
                "atom_id": row["atom_id"],
                "map_element": row["element"],
                "start_element": (
                    start_atoms[index]["element"]
                ),
                "final_element": (
                    final_atoms[index]["element"]
                ),
            })

    start_coordinates = {
        atom_ids[index]:
        start_atoms[index]["xyz_A"]
        for index in range(len(atom_ids))
    }

    final_coordinates = {
        atom_ids[index]:
        final_atoms[index]["xyz_A"]
        for index in range(len(atom_ids))
    }

    displacement_rows = []

    for index, atom_id in enumerate(atom_ids):
        shift = distance(
            start_coordinates[atom_id],
            final_coordinates[atom_id],
        )

        fixed = fixed_status[atom_id]

        displacement_rows.append({
            "index_0based": index,
            "atom_id": atom_id,
            "element": elements[atom_id],
            "v7a_fixed": fixed,
            "displacement_A": shift,
            "fixed_tolerance_A": (
                FIXED_ATOM_TOLERANCE_A
                if fixed
                else ""
            ),
            "fixed_displacement_pass": (
                shift <= FIXED_ATOM_TOLERANCE_A
                if fixed
                else True
            ),
        })

    nominal_edges = {
        frozenset((
            row["first_atom"],
            row["second_atom"],
        ))
        for row in edge_rows
    }

    final_audit = audit_geometry(
        final_coordinates,
        atom_ids,
        elements,
        nominal_edges,
        edge_rows,
    )

    trajectory_rows = []
    trajectory_failures = []

    for frame in frames:
        atoms = frame["atoms"]

        frame_identity = (
            len(atoms) == len(map_rows)
            and all(
                atoms[index]["element"]
                == map_rows[index]["element"]
                for index in range(len(map_rows))
            )
        )

        coordinates = {
            atom_ids[index]:
            atoms[index]["xyz_A"]
            for index in range(len(atom_ids))
        }

        audit = audit_geometry(
            coordinates,
            atom_ids,
            elements,
            nominal_edges,
            edge_rows,
        )

        frame_pass = all((
            frame_identity,
            len(audit["bond_failures"]) == 0,
            len(audit["lost_edges"]) == 0,
            len(audit["gained_edges"]) == 0,
            len(audit["degree_failures"]) == 0,
            len(
                audit["overcoordinated_atoms"]
            ) == 0,
            len(audit["components"]) == 1,
            len(audit["hard_contacts"]) == 0,
            len(
                audit[
                    "modified_margin_failures"
                ]
            ) == 0,
            REQUIRED_V7A_EDGES
            <= audit["geometric_edges"],
            not any(
                pair
                in audit["geometric_edges"]
                for pair
                in FORBIDDEN_V6B_EDGES
            ),
        ))

        trajectory_rows.append({
            "frame_index_0based": (
                frame["frame_index_0based"]
            ),
            "energy_Eh": (
                extract_energy(frame["comment"])
            ),
            "atom_count": len(atoms),
            "identity_pass": frame_identity,
            "geometric_edge_count": len(
                audit["geometric_edges"]
            ),
            "bond_failure_count": len(
                audit["bond_failures"]
            ),
            "lost_edge_count": len(
                audit["lost_edges"]
            ),
            "gained_edge_count": len(
                audit["gained_edges"]
            ),
            "degree_failure_count": len(
                audit["degree_failures"]
            ),
            "overcoordinated_atom_count": len(
                audit["overcoordinated_atoms"]
            ),
            "connected_component_count": len(
                audit["components"]
            ),
            "hard_contact_count": len(
                audit["hard_contacts"]
            ),
            "modified_near_contact_failure_count": len(
                audit[
                    "modified_margin_failures"
                ]
            ),
            "minimum_bond_margin_A": (
                audit["minimum_bond_margin_A"]
            ),
            "minimum_nonnominal_margin_A": (
                audit[
                    "minimum_nonnominal_margin_A"
                ]
            ),
            "frame_pass": frame_pass,
        })

        if not frame_pass:
            trajectory_failures.append(
                frame["frame_index_0based"]
            )

    composition = Counter(
        atom["element"]
        for atom in final_atoms
    )

    maximum_fixed_displacement = max(
        float(row["displacement_A"])
        for row in displacement_rows
        if row["v7a_fixed"]
    )

    maximum_total_displacement = max(
        float(row["displacement_A"])
        for row in displacement_rows
    )

    fixed_atoms_preserved = all(
        row["fixed_displacement_pass"]
        for row in displacement_rows
        if row["v7a_fixed"]
    )

    gates = {
        "post_QM_readiness_gate": True,
        "atom_identity_and_order": (
            len(identity_failures) == 0
        ),
        "atom_count_52": (
            len(final_atoms)
            == EXPECTED_ATOM_COUNT
        ),
        "composition_B17_N14_H21": (
            dict(composition)
            == EXPECTED_COMPOSITION
        ),
        "obsolete_HCAPV2_absent": (
            OBSOLETE_ATOM not in atom_ids
        ),
        "nominal_edge_count_57": (
            len(nominal_edges)
            == EXPECTED_EDGE_COUNT
        ),
        "fixed_atoms_preserved": (
            fixed_atoms_preserved
        ),
        "all_nominal_bonds_in_range": (
            len(
                final_audit[
                    "bond_failures"
                ]
            )
            == 0
        ),
        "no_lost_nominal_edges": (
            len(final_audit["lost_edges"])
            == 0
        ),
        "no_geometric_reconnectivity": (
            len(final_audit["gained_edges"])
            == 0
        ),
        "nominal_degree_exact": (
            len(
                final_audit[
                    "degree_failures"
                ]
            )
            == 0
        ),
        "no_geometric_overcoordination": (
            len(
                final_audit[
                    "overcoordinated_atoms"
                ]
            )
            == 0
        ),
        "single_connected_component": (
            len(final_audit["components"])
            == 1
        ),
        "no_hard_contacts": (
            len(
                final_audit[
                    "hard_contacts"
                ]
            )
            == 0
        ),
        "modified_region_near_contact_gate": (
            len(
                final_audit[
                    "modified_margin_failures"
                ]
            )
            == 0
        ),
        "required_V7A_edges_present": (
            REQUIRED_V7A_EDGES
            <= final_audit[
                "geometric_edges"
            ]
        ),
        "V6B_failure_mode_absent": (
            not any(
                pair
                in final_audit[
                    "geometric_edges"
                ]
                for pair
                in FORBIDDEN_V6B_EDGES
            )
        ),
        "trajectory_available_and_parseable": (
            len(frames) >= 1
        ),
        "all_trajectory_frames_pass": (
            len(trajectory_failures) == 0
        ),
    }

    passed = all(gates.values())

    decision = (
        "QM_F06_UPPER_V7A_R1_"
        "POST_QM_GATE_PASS_"
        "RESP_INPUT_PREPARATION_AUTHORIZED"
        if passed
        else
        "QM_F06_UPPER_V7A_R1_"
        "POST_QM_GATE_FAIL_"
        "STRUCTURAL_REVIEW_REQUIRED"
    )

    write_csv(
        OUTPUT_BONDS,
        final_audit["bond_rows"],
    )

    write_csv(
        OUTPUT_CONTACTS,
        final_audit["contact_rows"],
    )

    write_csv(
        OUTPUT_VALENCE,
        final_audit["valence_rows"],
    )

    write_csv(
        OUTPUT_DISPLACEMENTS,
        displacement_rows,
    )

    write_csv(
        OUTPUT_TRAJECTORY,
        trajectory_rows,
    )

    final_lines = [
        str(len(final_atoms)),
        (
            "QM_F06 UPPER V7-A R1 final "
            "post-QM audited geometry"
        ),
    ]

    for atom in final_atoms:
        x, y, z = atom["xyz_A"]

        final_lines.append(
            f"{atom['element']:2s} "
            f"{x: .12f} "
            f"{y: .12f} "
            f"{z: .12f}"
        )

    OUTPUT_FINAL_XYZ.write_text(
        "\n".join(final_lines) + "\n",
        encoding="utf-8",
    )

    report = {
        "model": "QM_F06_UPPER_V7A_R1",
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "execution_directory": str(
            execution_dir
        ),
        "upstream_readiness": {
            "report": str(
                READINESS_REPORT
            ),
            "decision": (
                readiness_decision
            ),
            "post_QM_structural_audit_authorized": (
                readiness_authorized
            ),
            "execution_directory": (
                readiness_execution_directory
            ),
            "execution_identity_match": (
                readiness_execution_directory
                == str(execution_dir)
            ),
        },
        "decision": decision,
        "gates": gates,
        "summary": {
            "atom_count": len(final_atoms),
            "composition": dict(
                sorted(composition.items())
            ),
            "nominal_edge_count": len(
                nominal_edges
            ),
            "geometric_edge_count": len(
                final_audit[
                    "geometric_edges"
                ]
            ),
            "maximum_total_displacement_A": (
                maximum_total_displacement
            ),
            "maximum_fixed_displacement_A": (
                maximum_fixed_displacement
            ),
            "fixed_tolerance_A": (
                FIXED_ATOM_TOLERANCE_A
            ),
            "bond_failure_count": len(
                final_audit["bond_failures"]
            ),
            "lost_edge_count": len(
                final_audit["lost_edges"]
            ),
            "gained_edge_count": len(
                final_audit["gained_edges"]
            ),
            "degree_failure_count": len(
                final_audit[
                    "degree_failures"
                ]
            ),
            "overcoordinated_atom_count": len(
                final_audit[
                    "overcoordinated_atoms"
                ]
            ),
            "connected_component_count": len(
                final_audit["components"]
            ),
            "hard_contact_count": len(
                final_audit["hard_contacts"]
            ),
            "modified_near_contact_failure_count": len(
                final_audit[
                    "modified_margin_failures"
                ]
            ),
            "minimum_bond_margin_A": (
                final_audit[
                    "minimum_bond_margin_A"
                ]
            ),
            "minimum_nonnominal_margin_A": (
                final_audit[
                    "minimum_nonnominal_margin_A"
                ]
            ),
            "trajectory_frame_count": (
                len(frames)
            ),
            "trajectory_failure_frames": (
                trajectory_failures
            ),
        },
        "failures": {
            "identity_failures": (
                identity_failures
            ),
            "bond_failures": (
                final_audit["bond_failures"]
            ),
            "lost_edges": [
                sorted(pair)
                for pair
                in final_audit["lost_edges"]
            ],
            "gained_edges": [
                sorted(pair)
                for pair
                in final_audit["gained_edges"]
            ],
            "degree_failures": (
                final_audit[
                    "degree_failures"
                ]
            ),
            "overcoordinated_atoms": (
                final_audit[
                    "overcoordinated_atoms"
                ]
            ),
            "hard_contacts": (
                final_audit["hard_contacts"]
            ),
            "modified_near_contact_failures": (
                final_audit[
                    "modified_margin_failures"
                ]
            ),
            "trajectory_failure_frames": (
                trajectory_failures
            ),
        },
        "outputs": {
            "final_XYZ": str(
                OUTPUT_FINAL_XYZ
            ),
            "bond_audit": str(
                OUTPUT_BONDS
            ),
            "contact_audit": str(
                OUTPUT_CONTACTS
            ),
            "valence_audit": str(
                OUTPUT_VALENCE
            ),
            "displacement_audit": str(
                OUTPUT_DISPLACEMENTS
            ),
            "trajectory_audit": str(
                OUTPUT_TRAJECTORY
            ),
        },
        "authorizations": {
            "RESP_input_preparation_authorized": (
                passed
            ),
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": (
                False
            ),
            "MD_authorized": False,
        },
    }

    OUTPUT_REPORT.write_text(
        json.dumps(
            report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 112)
    print(
        "QM_F06 UPPER V7-A R1 "
        "FINAL POST-QM STRUCTURAL AUDIT"
    )
    print("=" * 112)

    for name, value in gates.items():
        print(
            f"{name:60s}: "
            f"{'PASS' if value else 'FAIL'}"
        )

    print()
    print("Atoms:", len(final_atoms))
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
        len(
            final_audit[
                "geometric_edges"
            ]
        ),
    )
    print(
        "Maximum total displacement A:",
        maximum_total_displacement,
    )
    print(
        "Maximum fixed displacement A:",
        maximum_fixed_displacement,
    )
    print(
        "Bond failures:",
        len(final_audit["bond_failures"]),
    )
    print(
        "Lost edges:",
        len(final_audit["lost_edges"]),
    )
    print(
        "Gained edges:",
        len(final_audit["gained_edges"]),
    )
    print(
        "Overcoordinated atoms:",
        len(
            final_audit[
                "overcoordinated_atoms"
            ]
        ),
    )
    print(
        "Hard contacts:",
        len(final_audit["hard_contacts"]),
    )
    print(
        "Modified near-contact failures:",
        len(
            final_audit[
                "modified_margin_failures"
            ]
        ),
    )
    print(
        "Trajectory frames:",
        len(frames),
    )
    print(
        "Trajectory failure frames:",
        trajectory_failures,
    )

    print()
    print("Decision:", decision)
    print("Final XYZ:", OUTPUT_FINAL_XYZ)
    print("Report:", OUTPUT_REPORT)
    print()
    print(
        "RESP input preparation authorized:",
        passed,
    )
    print("RESP execution authorized: False")
    print(
        "Force-field adoption authorized: False"
    )
    print("MD authorized: False")

    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
