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

V6A_DIR = ROOT / (
    "runs/phase1A/"
    "day033_qm_f06_upper_v6a_topology_closure"
)

V6A_XYZ = V6A_DIR / (
    "QM_F06_UPPER_V6A_TOPOLOGY_CLOSURE_start.xyz"
)

V6A_MAP = V6A_DIR / (
    "QM_F06_UPPER_V6A_atom_role_provenance_map.csv"
)

V6A_EDGES = V6A_DIR / (
    "QM_F06_UPPER_V6A_nominal_edges.csv"
)

V6A_REPORT = V6A_DIR / (
    "QM_F06_UPPER_V6A_TOPOLOGY_CLOSURE_REPORT.json"
)

SEARCH_DIR = ROOT / (
    "runs/phase1A/"
    "day033_qm_f06_upper_v6b_local_geometry_search"
)

CANDIDATE_XYZ = SEARCH_DIR / (
    "QM_F06_UPPER_V6B_LOCAL_GEOMETRY_BEST.xyz"
)

SEARCH_RANKING = SEARCH_DIR / (
    "QM_F06_UPPER_V6B_local_geometry_ranking.csv"
)

SEARCH_REPORT = SEARCH_DIR / (
    "QM_F06_UPPER_V6B_LOCAL_GEOMETRY_SEARCH.json"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day033_qm_f06_upper_v6b_pre_qm_audit"
)

FORMAL_XYZ = OUTPUT_DIR / (
    "QM_F06_UPPER_V6B_start.xyz"
)

ATOM_AUDIT_CSV = OUTPUT_DIR / (
    "QM_F06_UPPER_V6B_atom_displacements.csv"
)

BOND_AUDIT_CSV = OUTPUT_DIR / (
    "QM_F06_UPPER_V6B_bond_audit.csv"
)

CONTACT_AUDIT_CSV = OUTPUT_DIR / (
    "QM_F06_UPPER_V6B_contacts.csv"
)

REPORT_PATH = OUTPUT_DIR / (
    "QM_F06_UPPER_V6B_PRE_QM_AUDIT.json"
)


EXPECTED_ATOM_COUNT = 48

EXPECTED_COMPOSITION = Counter({
    "B": 16,
    "N": 13,
    "H": 19,
})

AUTHORIZED_MOVED_ATOMS = {
    "S:1739",
    "BR4:UPPER:14:1",
    "HCAPV5B:UPPER:BR4_14_1:BR4_14_2",
}

REMOVED_HYDROGENS = {
    "H4:UPPER:0116:0",
    "H4:UPPER:0117:0",
    "H4:UPPER:0170:0",
    "H4:UPPER:0203:0",
}

CLOSURE_PAIRS = {
    tuple(sorted((
        "BR4:UPPER:00:3",
        "S:1739",
    ))),
    tuple(sorted((
        "BR4:UPPER:00:4",
        "BR4:UPPER:14:1",
    ))),
}

CAP_OWNER = "BR4:UPPER:14:1"

CAP_ID = (
    "HCAPV5B:UPPER:BR4_14_1:BR4_14_2"
)

UNCHANGED_ATOM_TOLERANCE_A = 5.0e-7
RIGID_CAP_VECTOR_TOLERANCE_A = 5.0e-7
MAXIMUM_AUTHORIZED_SHIFT_A = 0.95

BN_MIN_A = 1.25
BN_MAX_A = 1.90

BH_MIN_A = 0.90
BH_MAX_A = 1.35

NH_MIN_A = 0.80
NH_MAX_A = 1.25

CLOSURE_MIN_A = 1.40
CLOSURE_MAX_A = 1.70

MINIMUM_CLOSURE_MARGIN_A = 0.02
MINIMUM_MODIFIED_BOND_MARGIN_A = 0.08
MINIMUM_NONNOMINAL_CLEARANCE_A = 0.06

HH_HARD_CONTACT_A = 0.70
HX_HARD_CONTACT_A = 0.75
HEAVY_HEAVY_HARD_CONTACT_A = 1.10


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
            f"Incomplete XYZ: {path}"
        )

    atoms = []

    for index, line in enumerate(
        coordinate_lines
    ):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Malformed XYZ line: {line}"
            )

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


def subtract(a, b):
    return tuple(
        x - y
        for x, y in zip(a, b)
    )


def norm(vector) -> float:
    return math.sqrt(
        sum(value * value for value in vector)
    )


def distance(a, b) -> float:
    return norm(subtract(a, b))


def vector_difference(a, b) -> float:
    return norm(subtract(a, b))


def bond_limits(
    first_element: str,
    second_element: str,
):
    pair = {
        first_element,
        second_element,
    }

    if pair == {"B", "N"}:
        return BN_MIN_A, BN_MAX_A, "B-N"

    if pair == {"B", "H"}:
        return BH_MIN_A, BH_MAX_A, "B-H"

    if pair == {"N", "H"}:
        return NH_MIN_A, NH_MAX_A, "N-H"

    return None


def hard_contact_threshold(
    first_element: str,
    second_element: str,
) -> float:
    if (
        first_element == "H"
        and second_element == "H"
    ):
        return HH_HARD_CONTACT_A

    if (
        first_element == "H"
        or second_element == "H"
    ):
        return HX_HARD_CONTACT_A

    return HEAVY_HEAVY_HARD_CONTACT_A


def main() -> None:
    required_files = (
        V6A_XYZ,
        V6A_MAP,
        V6A_EDGES,
        V6A_REPORT,
        CANDIDATE_XYZ,
        SEARCH_RANKING,
        SEARCH_REPORT,
    )

    for path in required_files:
        require_file(path)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    v6a_report = json.loads(
        V6A_REPORT.read_text(
            encoding="utf-8"
        )
    )

    search_report = json.loads(
        SEARCH_REPORT.read_text(
            encoding="utf-8"
        )
    )

    ranking_rows = read_csv(
        SEARCH_RANKING
    )

    if not ranking_rows:
        raise RuntimeError(
            "V6-B ranking is empty."
        )

    selected_row = ranking_rows[0]

    selected_rank_gate = (
        int(selected_row["rank"]) == 1
    )

    search_decision_gate = (
        search_report.get("decision")
        == (
            "QM_F06_UPPER_V6B_LOCAL_GEOMETRY_"
            "CANDIDATE_FOUND_FORMAL_AUDIT_REQUIRED"
        )
    )

    v6a_decision_gate = (
        v6a_report.get("decision")
        == (
            "QM_F06_UPPER_V6A_TOPOLOGY_CLOSURE_"
            "CONSTRUCTED_STRUCTURAL_AUDIT_REQUIRED"
        )
    )

    map_rows = read_csv(V6A_MAP)

    retained_rows = [
        row
        for row in map_rows
        if row["v6a_retained"].strip().lower()
        == "true"
    ]

    retained_rows.sort(
        key=lambda row: int(
            row["v6a_index_0based"]
        )
    )

    atom_ids = [
        row["atom_id"]
        for row in retained_rows
    ]

    id_set = set(atom_ids)

    removed_atoms_absent_gate = (
        REMOVED_HYDROGENS.isdisjoint(id_set)
    )

    source_atoms = read_xyz(V6A_XYZ)
    candidate_atoms = read_xyz(
        CANDIDATE_XYZ
    )

    if (
        len(source_atoms) != len(retained_rows)
        or len(candidate_atoms)
        != len(retained_rows)
    ):
        raise RuntimeError(
            "XYZ/map atom-count mismatch."
        )

    source = {}
    candidate = {}
    element_by_id = {}

    identity_gate = True

    for row, source_atom, candidate_atom in zip(
        retained_rows,
        source_atoms,
        candidate_atoms,
    ):
        atom_id = row["atom_id"]
        expected_element = row["element"]

        if (
            source_atom["element"]
            != expected_element
            or candidate_atom["element"]
            != expected_element
        ):
            identity_gate = False

        element_by_id[atom_id] = (
            expected_element
        )

        source[atom_id] = (
            source_atom["xyz_A"]
        )

        candidate[atom_id] = (
            candidate_atom["xyz_A"]
        )

    atom_count_gate = (
        len(candidate) == EXPECTED_ATOM_COUNT
    )

    composition = Counter(
        element_by_id.values()
    )

    composition_gate = (
        composition == EXPECTED_COMPOSITION
    )

    atom_records = []
    unchanged_failures = []
    authorized_shift_failures = []

    for atom_id in atom_ids:
        value = distance(
            source[atom_id],
            candidate[atom_id],
        )

        authorized = (
            atom_id in AUTHORIZED_MOVED_ATOMS
        )

        unchanged_pass = (
            authorized
            or value
            <= UNCHANGED_ATOM_TOLERANCE_A
        )

        authorized_shift_pass = (
            not authorized
            or value
            <= MAXIMUM_AUTHORIZED_SHIFT_A
        )

        record = {
            "index_0based": (
                atom_ids.index(atom_id)
            ),
            "atom_id": atom_id,
            "element": element_by_id[atom_id],
            "authorized_to_move": authorized,
            "displacement_A": value,
            "unchanged_atom_pass": (
                unchanged_pass
            ),
            "authorized_shift_pass": (
                authorized_shift_pass
            ),
        }

        atom_records.append(record)

        if not unchanged_pass:
            unchanged_failures.append(record)

        if not authorized_shift_pass:
            authorized_shift_failures.append(
                record
            )

    with ATOM_AUDIT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(atom_records[0]),
        )

        writer.writeheader()
        writer.writerows(atom_records)

    authorized_motion_only_gate = (
        len(unchanged_failures) == 0
    )

    authorized_shift_gate = (
        len(authorized_shift_failures) == 0
    )

    source_cap_vector = subtract(
        source[CAP_ID],
        source[CAP_OWNER],
    )

    candidate_cap_vector = subtract(
        candidate[CAP_ID],
        candidate[CAP_OWNER],
    )

    cap_vector_difference_A = (
        vector_difference(
            source_cap_vector,
            candidate_cap_vector,
        )
    )

    rigid_cap_translation_gate = (
        cap_vector_difference_A
        <= RIGID_CAP_VECTOR_TOLERANCE_A
    )

    edge_rows = read_csv(V6A_EDGES)

    nominal_edges = set()
    adjacency = {
        atom_id: set()
        for atom_id in atom_ids
    }

    for row in edge_rows:
        first = row["first_atom"]
        second = row["second_atom"]

        if (
            first not in id_set
            or second not in id_set
        ):
            raise RuntimeError(
                "Nominal edge references an "
                "unretained atom: "
                f"{first}--{second}"
            )

        pair = canonical_pair(
            first,
            second,
        )

        nominal_edges.add(pair)
        adjacency[first].add(second)
        adjacency[second].add(first)

    closure_inventory_gate = (
        CLOSURE_PAIRS.issubset(
            nominal_edges
        )
    )

    degree_records = []
    degree_failures = []

    for atom_id in atom_ids:
        degree = len(adjacency[atom_id])

        expected_degree = (
            1
            if element_by_id[atom_id] == "H"
            else 3
        )

        passed = (
            degree == expected_degree
        )

        record = {
            "atom_id": atom_id,
            "element": element_by_id[
                atom_id
            ],
            "degree": degree,
            "expected_degree": (
                expected_degree
            ),
            "neighbors": "|".join(
                sorted(adjacency[atom_id])
            ),
            "pass": passed,
        }

        degree_records.append(record)

        if not passed:
            degree_failures.append(record)

    degree_gate = (
        len(degree_failures) == 0
    )

    visited = set()
    queue = deque([atom_ids[0]])

    while queue:
        current = queue.popleft()

        if current in visited:
            continue

        visited.add(current)

        for neighbor in adjacency[current]:
            if neighbor not in visited:
                queue.append(neighbor)

    connected_gate = (
        len(visited) == len(atom_ids)
    )

    bond_records = []
    bond_failures = []

    minimum_modified_margin = float(
        "inf"
    )

    closure_records = []
    closure_failures = []

    for first, second in sorted(
        nominal_edges
    ):
        first_element = element_by_id[first]
        second_element = element_by_id[second]

        limits = bond_limits(
            first_element,
            second_element,
        )

        if limits is None:
            raise RuntimeError(
                "Unsupported nominal bond class: "
                f"{first_element}-{second_element}"
            )

        minimum, maximum, bond_class = (
            limits
        )

        value = distance(
            candidate[first],
            candidate[second],
        )

        margin = min(
            value - minimum,
            maximum - value,
        )

        passed = (
            minimum <= value <= maximum
        )

        modified_region = (
            first in AUTHORIZED_MOVED_ATOMS
            or second in AUTHORIZED_MOVED_ATOMS
        )

        if modified_region:
            minimum_modified_margin = min(
                minimum_modified_margin,
                margin,
            )

        record = {
            "first_atom": first,
            "first_element": first_element,
            "second_atom": second,
            "second_element": second_element,
            "bond_class": bond_class,
            "distance_A": value,
            "minimum_A": minimum,
            "maximum_A": maximum,
            "margin_A": margin,
            "modified_region": (
                modified_region
            ),
            "pass": passed,
        }

        bond_records.append(record)

        if not passed:
            bond_failures.append(record)

        pair = canonical_pair(
            first,
            second,
        )

        if pair in CLOSURE_PAIRS:
            closure_margin = min(
                value - CLOSURE_MIN_A,
                CLOSURE_MAX_A - value,
            )

            closure_pass = (
                CLOSURE_MIN_A
                <= value
                <= CLOSURE_MAX_A
                and closure_margin
                >= MINIMUM_CLOSURE_MARGIN_A
            )

            closure_record = {
                "first_atom": first,
                "second_atom": second,
                "distance_A": value,
                "minimum_A": CLOSURE_MIN_A,
                "maximum_A": CLOSURE_MAX_A,
                "margin_A": closure_margin,
                "pass": closure_pass,
            }

            closure_records.append(
                closure_record
            )

            if not closure_pass:
                closure_failures.append(
                    closure_record
                )

    with BOND_AUDIT_CSV.open(
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

    bond_gate = (
        len(bond_failures) == 0
    )

    modified_margin_gate = (
        minimum_modified_margin
        >= MINIMUM_MODIFIED_BOND_MARGIN_A
    )

    closure_geometry_gate = (
        len(closure_records) == 2
        and len(closure_failures) == 0
    )

    geometric_edges = set()
    geometric_distances = {}

    contact_records = []
    hard_contacts = []

    minimum_nonnominal_clearance = (
        float("inf")
    )

    limiting_nonnominal_pair = None

    for first_position, first in enumerate(
        atom_ids
    ):
        for second in atom_ids[
            first_position + 1:
        ]:
            pair = canonical_pair(
                first,
                second,
            )

            first_element = element_by_id[first]
            second_element = element_by_id[
                second
            ]

            value = distance(
                candidate[first],
                candidate[second],
            )

            limits = bond_limits(
                first_element,
                second_element,
            )

            geometrically_bonded = False

            if limits is not None:
                minimum, maximum, _ = limits

                geometrically_bonded = (
                    minimum <= value <= maximum
                )

                if geometrically_bonded:
                    geometric_edges.add(pair)
                    geometric_distances[
                        pair
                    ] = value

                # The quantitative clearance gate is local to
                # V6-B-modified geometry. Static inherited pairs
                # remain protected independently by the global
                # reconnectivity, overcoordination, and hard-contact
                # gates below.
                modified_nonnominal_pair = (
                    pair not in nominal_edges
                    and (
                        first in AUTHORIZED_MOVED_ATOMS
                        or second in AUTHORIZED_MOVED_ATOMS
                    )
                )

                if modified_nonnominal_pair:
                    clearance = value - maximum

                    if (
                        clearance
                        < minimum_nonnominal_clearance
                    ):
                        minimum_nonnominal_clearance = (
                            clearance
                        )
                        limiting_nonnominal_pair = (
                            f"{first}--{second}"
                        )

            hard_threshold = (
                hard_contact_threshold(
                    first_element,
                    second_element,
                )
            )

            hard_contact = (
                pair not in nominal_edges
                and value < hard_threshold
            )

            record = {
                "first_atom": first,
                "first_element": (
                    first_element
                ),
                "second_atom": second,
                "second_element": (
                    second_element
                ),
                "distance_A": value,
                "nominal_edge": (
                    pair in nominal_edges
                ),
                "geometrically_bonded": (
                    geometrically_bonded
                ),
                "hard_contact_threshold_A": (
                    hard_threshold
                ),
                "hard_contact": hard_contact,
            }

            contact_records.append(record)

            if hard_contact:
                hard_contacts.append(record)

    with CONTACT_AUDIT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                contact_records[0]
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

    reconnectivity_gate = (
        len(gained_edges) == 0
        and len(lost_edges) == 0
    )

    nonnominal_clearance_gate = (
        minimum_nonnominal_clearance
        >= MINIMUM_NONNOMINAL_CLEARANCE_A
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

    for atom_id in atom_ids:
        maximum_degree = (
            1
            if element_by_id[atom_id] == "H"
            else 3
        )

        degree = len(
            geometric_adjacency[atom_id]
        )

        if degree > maximum_degree:
            overcoordinated.append({
                "atom_id": atom_id,
                "element": element_by_id[
                    atom_id
                ],
                "degree": degree,
                "maximum_degree": (
                    maximum_degree
                ),
                "neighbors": "|".join(
                    sorted(
                        geometric_adjacency[
                            atom_id
                        ]
                    )
                ),
            })

    overcoordination_gate = (
        len(overcoordinated) == 0
    )

    hard_contact_gate = (
        len(hard_contacts) == 0
    )

    gates = {
        "v6a_construction_decision": (
            v6a_decision_gate
        ),
        "v6b_search_decision": (
            search_decision_gate
        ),
        "selected_candidate_rank_1": (
            selected_rank_gate
        ),
        "atom_identity_and_order": (
            identity_gate
        ),
        "atom_count": atom_count_gate,
        "composition": composition_gate,
        "removed_hydrogens_absent": (
            removed_atoms_absent_gate
        ),
        "only_authorized_atoms_moved": (
            authorized_motion_only_gate
        ),
        "authorized_shifts_within_limit": (
            authorized_shift_gate
        ),
        "cap_translated_rigidly": (
            rigid_cap_translation_gate
        ),
        "closure_edges_present": (
            closure_inventory_gate
        ),
        "nominal_degree_exact": degree_gate,
        "single_connected_component": (
            connected_gate
        ),
        "all_nominal_bonds_in_range": (
            bond_gate
        ),
        "modified_region_bond_margin": (
            modified_margin_gate
        ),
        "closure_BN_geometry": (
            closure_geometry_gate
        ),
        "no_geometric_reconnectivity": (
            reconnectivity_gate
        ),
        "modified_region_nonnominal_pair_clearance": (
            nonnominal_clearance_gate
        ),
        "no_overcoordinated_atoms": (
            overcoordination_gate
        ),
        "no_hard_contacts": (
            hard_contact_gate
        ),
    }

    passed = all(gates.values())

    if passed:
        decision = (
            "QM_F06_UPPER_V6B_PRE_QM_GATE_PASS_"
            "ORCA_INPUT_DESIGN_AUTHORIZED"
        )
    else:
        decision = (
            "QM_F06_UPPER_V6B_PRE_QM_GATE_FAIL_"
            "STRUCTURAL_REVIEW_REQUIRED"
        )

    with FORMAL_XYZ.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            f"{len(atom_ids)}\n"
        )

        handle.write(
            "QM_F06 UPPER V6-B formally audited "
            "pre-QM geometry; "
            f"gate_pass={passed}\n"
        )

        for atom_id in atom_ids:
            x_value, y_value, z_value = (
                candidate[atom_id]
            )

            handle.write(
                f"{element_by_id[atom_id]:2s} "
                f"{x_value: .12f} "
                f"{y_value: .12f} "
                f"{z_value: .12f}\n"
            )

    files_for_hash = {
        "v6a_xyz": V6A_XYZ,
        "v6a_map": V6A_MAP,
        "v6a_edges": V6A_EDGES,
        "v6a_report": V6A_REPORT,
        "candidate_xyz": CANDIDATE_XYZ,
        "search_ranking": SEARCH_RANKING,
        "search_report": SEARCH_REPORT,
        "formal_xyz": FORMAL_XYZ,
        "atom_audit": ATOM_AUDIT_CSV,
        "bond_audit": BOND_AUDIT_CSV,
        "contact_audit": CONTACT_AUDIT_CSV,
    }

    report = {
        "decision": decision,
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "selected_candidate": {
            "rank": int(
                selected_row["rank"]
            ),
            "sample_index": int(
                selected_row["sample_index"]
            ),
            "reported_score": float(
                selected_row["score"]
            ),
        },
        "geometry": {
            "atom_count": len(atom_ids),
            "composition": dict(
                sorted(composition.items())
            ),
            "authorized_moved_atoms": sorted(
                AUTHORIZED_MOVED_ATOMS
            ),
            "maximum_authorized_displacement_A": (
                max(
                    record["displacement_A"]
                    for record in atom_records
                    if record[
                        "authorized_to_move"
                    ]
                )
            ),
            "cap_vector_difference_A": (
                cap_vector_difference_A
            ),
        },
        "topology": {
            "nominal_edge_count": len(
                nominal_edges
            ),
            "geometric_edge_count": len(
                geometric_edges
            ),
            "degree_failure_count": len(
                degree_failures
            ),
            "gained_edge_count": len(
                gained_edges
            ),
            "gained_edges": [
                f"{first}--{second}"
                for first, second in gained_edges
            ],
            "lost_edge_count": len(
                lost_edges
            ),
            "lost_edges": [
                f"{first}--{second}"
                for first, second in lost_edges
            ],
            "overcoordinated_count": len(
                overcoordinated
            ),
            "overcoordinated_atoms": (
                overcoordinated
            ),
        },
        "geometry_margins": {
            "minimum_modified_bond_margin_A": (
                minimum_modified_margin
            ),
            "minimum_required_modified_margin_A": (
                MINIMUM_MODIFIED_BOND_MARGIN_A
            ),
            "minimum_modified_region_nonnominal_clearance_A": (
                minimum_nonnominal_clearance
            ),
            "minimum_required_modified_region_clearance_A": (
                MINIMUM_NONNOMINAL_CLEARANCE_A
            ),
            "clearance_scope": (
                "NONNOMINAL_PAIRS_INVOLVING_AT_LEAST_ONE_"
                "V6B_AUTHORIZED_MOVED_ATOM"
            ),
            "limiting_nonnominal_pair": (
                limiting_nonnominal_pair
            ),
            "closure_records": (
                closure_records
            ),
        },
        "failure_counts": {
            "unchanged_atom_failures": len(
                unchanged_failures
            ),
            "authorized_shift_failures": len(
                authorized_shift_failures
            ),
            "degree_failures": len(
                degree_failures
            ),
            "bond_failures": len(
                bond_failures
            ),
            "closure_failures": len(
                closure_failures
            ),
            "hard_contacts": len(
                hard_contacts
            ),
        },
        "gates": gates,
        "structural_pass": passed,
        "files": {
            key: str(path.relative_to(ROOT))
            for key, path in files_for_hash.items()
        },
        "sha256": {
            key: sha256(path)
            for key, path in files_for_hash.items()
        },
        "authorization": {
            "ORCA_input_design_authorized": (
                passed
            ),
            "ORCA_execution_authorized": False,
            "RESP_input_preparation_authorized": (
                False
            ),
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": (
                False
            ),
            "MD_authorized": False,
        },
    }

    REPORT_PATH.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 104)
    print("QM_F06 UPPER V6-B FORMAL PRE-QM AUDIT")
    print("=" * 104)

    for name, value in gates.items():
        print(
            f"{name:48s}: "
            f"{'PASS' if value else 'FAIL'}"
        )

    print()
    print("Atoms:", len(atom_ids))
    print(
        "Composition:",
        dict(sorted(composition.items())),
    )
    print("Nominal edges:", len(
        nominal_edges
    ))
    print("Geometric edges:", len(
        geometric_edges
    ))

    print()
    print(
        "Minimum modified-region bond margin A:",
        minimum_modified_margin,
    )
    print(
        "Minimum modified-region nonnominal clearance A:",
        minimum_nonnominal_clearance,
    )
    print(
        "Limiting nonnominal pair:",
        limiting_nonnominal_pair,
    )
    print(
        "Cap rigid-vector difference A:",
        cap_vector_difference_A,
    )

    print()
    print("Closure B-N records:")

    for record in closure_records:
        print(
            f"  {record['first_atom']:28s} -- "
            f"{record['second_atom']:28s} "
            f"{record['distance_A']:.6f} Å | "
            f"margin={record['margin_A']:.6f} Å | "
            f"{'PASS' if record['pass'] else 'FAIL'}"
        )

    print()
    print(
        "Unchanged-atom failures:",
        len(unchanged_failures),
    )
    print(
        "Authorized-shift failures:",
        len(authorized_shift_failures),
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

    if gained_edges:
        print()
        print("Gained geometric edges:")

        for first, second in gained_edges:
            print(
                f"  {first} -- {second} "
                f"{geometric_distances[(first, second)]:.6f} Å"
            )

    if lost_edges:
        print()
        print("Lost nominal edges:")

        for first, second in lost_edges:
            print(
                f"  {first} -- {second} "
                f"{distance(candidate[first], candidate[second]):.6f} Å"
            )

    print()
    print("Decision:", decision)
    print("Formal XYZ:", FORMAL_XYZ)
    print("Report:", REPORT_PATH)
    print("Atom audit:", ATOM_AUDIT_CSV)
    print("Bond audit:", BOND_AUDIT_CSV)
    print(
        "Contact audit:",
        CONTACT_AUDIT_CSV,
    )

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
