#!/usr/bin/env python3
"""
Full pre-QM structural audit for QM_F06 UPPER V4.

The audit validates:
- construction report authorization;
- atom count and composition;
- deterministic atom-map consistency;
- connectivity and nominal valence;
- restored canonical R2 bonds;
- artificial cap identities and bond lengths;
- B-H2 geometry at S:1738;
- hard-contact screening with topological exclusions;
- cap ownership and nearest-heavy-center identity;
- preliminary constraint-design classification.

No ORCA input is generated and no ORCA execution is authorized.
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
    "runs/phase1A/day030_qm_f06_upper_v4_construction"
)

XYZ_PATH = (
    CONSTRUCTION_DIR
    / "QM_F06_UPPER_V4_start.xyz"
)

MAP_PATH = (
    CONSTRUCTION_DIR
    / "QM_F06_UPPER_V4_atom_role_provenance_map.csv"
)

CAPS_PATH = (
    CONSTRUCTION_DIR
    / "QM_F06_UPPER_V4_artificial_caps.csv"
)

CONSTRUCTION_REPORT = (
    CONSTRUCTION_DIR
    / "QM_F06_UPPER_V4_CONSTRUCTION_REPORT.json"
)

GRAPH_ROOT = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph"
)

GRAPH_NODES = (
    GRAPH_ROOT
    / "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

GRAPH_EDGES = (
    GRAPH_ROOT
    / "r2_selected_longer_bn_bridge_graph_edges.csv"
)

V3_VALIDATED_BOND_AUDIT = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V3/pre_qm_audit/"
    "QM_F06_UPPER_BOUNDARY_V3_bond_audit.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_pre_qm_audit"
)

REPORT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_PRE_QM_AUDIT.json"
)

CONTACTS_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_contacts.csv"
)

VALENCE_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_valence.csv"
)

CAP_AUDIT_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_cap_audit.csv"
)

EXPECTED_COMPOSITION = Counter({
    "B": 15,
    "N": 11,
    "H": 20,
})

EXPECTED_ATOM_COUNT = 46

BN_MIN_A = 1.30
BN_MAX_A = 1.70

BH_MIN_A = 1.10
BH_MAX_A = 1.30

NH_MIN_A = 0.90
NH_MAX_A = 1.15

HH_HARD_CONTACT_A = 1.20
HX_HARD_CONTACT_A = 0.95
HEAVY_HEAVY_HARD_CONTACT_A = 1.20

BH2_MIN_ANGLE_DEG = 90.0
BH2_MAX_ANGLE_DEG = 150.0

CAP_NEAREST_CENTER_MARGIN_A = 0.15


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

    for index, line in enumerate(lines[2:2 + declared]):
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


def distance(first, second):
    return math.sqrt(sum(
        (a - b) ** 2
        for a, b in zip(first, second)
    ))


def angle(first, center, second):
    vector_1 = tuple(
        a - b
        for a, b in zip(first, center)
    )

    vector_2 = tuple(
        a - b
        for a, b in zip(second, center)
    )

    norm_1 = math.sqrt(sum(
        value * value
        for value in vector_1
    ))

    norm_2 = math.sqrt(sum(
        value * value
        for value in vector_2
    ))

    cosine = sum(
        a * b
        for a, b in zip(vector_1, vector_2)
    ) / (norm_1 * norm_2)

    cosine = max(-1.0, min(1.0, cosine))

    return math.degrees(math.acos(cosine))


def canonical_pair(first: str, second: str):
    return tuple(sorted((first, second)))


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for path in (
        XYZ_PATH,
        MAP_PATH,
        CAPS_PATH,
        CONSTRUCTION_REPORT,
        GRAPH_NODES,
        GRAPH_EDGES,
        V3_VALIDATED_BOND_AUDIT,
    ):
        require_file(path)

    atoms = read_xyz(XYZ_PATH)
    map_rows = read_csv(MAP_PATH)
    cap_rows = read_csv(CAPS_PATH)
    node_rows = read_csv(GRAPH_NODES)
    edge_rows = read_csv(GRAPH_EDGES)
    v3_bond_rows = read_csv(
        V3_VALIDATED_BOND_AUDIT
    )

    construction = json.loads(
        CONSTRUCTION_REPORT.read_text(
            encoding="utf-8"
        )
    )

    if not construction["basic_construction_pass"]:
        raise RuntimeError(
            "Construction report did not pass."
        )

    if not construction["authorization"][
        "pre_qm_structural_audit_authorized"
    ]:
        raise RuntimeError(
            "Pre-QM structural audit is not authorized."
        )

    if len(atoms) != len(map_rows):
        raise RuntimeError(
            "XYZ and atom-map row counts differ."
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

    node_data = {
        row["node_id"]: row
        for row in node_rows
    }

    canonical_adjacency = defaultdict(set)

    for row in edge_rows:
        first = row["source_node"]
        second = row["target_node"]

        canonical_adjacency[first].add(second)
        canonical_adjacency[second].add(first)

    artificial_caps = {
        row["cap_id"]: row
        for row in cap_rows
    }

    if len(artificial_caps) != 7:
        raise RuntimeError(
            "Expected seven artificial caps."
        )

    bonded_pairs = set()
    adjacency = defaultdict(set)

    # Import the complete connectivity already validated for V3.
    # Only retain edges whose two endpoints survive in V4.
    if not v3_bond_rows:
        raise RuntimeError(
            "The validated V3 bond audit contains no rows."
        )

    v3_bond_headers = set(v3_bond_rows[0])

    def resolve_column(
        headers,
        candidates,
        label,
    ):
        matches = [
            candidate
            for candidate in candidates
            if candidate in headers
        ]

        if len(matches) != 1:
            raise RuntimeError(
                f"Could not uniquely resolve {label}. "
                f"Candidates found: {matches}. "
                f"Available headers: {sorted(headers)}"
            )

        return matches[0]

    v3_first_key = resolve_column(
        v3_bond_headers,
        (
            "first_atom",
            "atom_id_1",
            "atom_1",
            "atom1",
            "atom_a",
            "node_a",
            "source_node",
            "inside_atom",
        ),
        "first V3 bond endpoint",
    )

    v3_second_key = resolve_column(
        v3_bond_headers,
        (
            "second_atom",
            "atom_id_2",
            "atom_2",
            "atom2",
            "atom_b",
            "node_b",
            "target_node",
            "outside_atom",
        ),
        "second V3 bond endpoint",
    )

    pass_candidates = [
        candidate
        for candidate in (
            "pass",
            "bond_pass",
            "distance_pass",
            "overall_pass",
            "valid",
        )
        if candidate in v3_bond_headers
    ]

    if len(pass_candidates) > 1:
        raise RuntimeError(
            "Multiple possible V3 bond pass columns found: "
            f"{pass_candidates}"
        )

    v3_pass_key = (
        pass_candidates[0]
        if pass_candidates
        else None
    )

    for row in v3_bond_rows:
        first = row[v3_first_key].strip()
        second = row[v3_second_key].strip()

        if not first or not second:
            raise RuntimeError(
                "Empty endpoint in validated V3 bond audit."
            )

        if first not in mapped or second not in mapped:
            continue

        if v3_pass_key is not None:
            pass_value = (
                row[v3_pass_key]
                .strip()
                .lower()
            )

            if pass_value not in {
                "true",
                "1",
                "yes",
                "pass",
            }:
                raise RuntimeError(
                    "A retained V3 edge did not pass its "
                    f"original audit: {first} -- {second}; "
                    f"{v3_pass_key}={row[v3_pass_key]!r}"
                )

        bonded_pairs.add(
            canonical_pair(first, second)
        )

    # Add canonical retained/restored real-atom bonds.
    for atom_id, record in mapped.items():
        if (
            record["artificial_cap"]
            .strip()
            .lower()
            == "true"
        ):
            continue

        if atom_id not in canonical_adjacency:
            continue

        for neighbor in canonical_adjacency[atom_id]:
            if neighbor not in mapped:
                continue

            neighbor_record = mapped[neighbor]

            if (
                neighbor_record["artificial_cap"]
                .strip()
                .lower()
                == "true"
            ):
                continue

            pair = canonical_pair(
                atom_id,
                neighbor,
            )

            bonded_pairs.add(pair)

    # Artificial cap bonds.
    for cap_id, cap_record in artificial_caps.items():
        center = cap_record["center_atom"]

        if cap_id not in mapped:
            raise RuntimeError(
                f"Cap absent from map: {cap_id}"
            )

        if center not in mapped:
            raise RuntimeError(
                f"Cap center absent from map: {center}"
            )

        bonded_pairs.add(
            canonical_pair(cap_id, center)
        )

    for first, second in bonded_pairs:
        adjacency[first].add(second)
        adjacency[second].add(first)

    # Graph connectivity.
    atom_ids = set(mapped)

    visited = set()
    queue = deque([next(iter(atom_ids))])

    while queue:
        atom_id = queue.popleft()

        if atom_id in visited:
            continue

        visited.add(atom_id)

        for neighbor in adjacency[atom_id]:
            if neighbor not in visited:
                queue.append(neighbor)

    connected_component_gate = (
        visited == atom_ids
    )

    # Valence audit.
    valence_records = []

    for atom_id, record in sorted(
        mapped.items(),
        key=lambda item: item[1]["index"],
    ):
        element = record["element"]
        degree = len(adjacency[atom_id])

        if element in {"B", "N"}:
            expected_degree = 3
        elif element == "H":
            expected_degree = 1
        else:
            expected_degree = None

        valence_pass = (
            expected_degree is not None
            and degree == expected_degree
        )

        valence_records.append({
            "index_0based": record["index"],
            "atom_id": atom_id,
            "element": element,
            "degree": degree,
            "expected_degree": expected_degree,
            "neighbors": "|".join(
                sorted(adjacency[atom_id])
            ),
            "pass": valence_pass,
        })

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

    valence_gate = all(
        row["pass"]
        for row in valence_records
    )

    # Cap audit.
    cap_audit_records = []

    for cap_id, cap_spec in artificial_caps.items():
        cap = mapped[cap_id]
        center_id = cap_spec["center_atom"]
        center = mapped[center_id]

        target_length = float(
            cap_spec["target_cap_bond_length_A"]
        )

        realized_length = distance(
            cap["xyz_A"],
            center["xyz_A"],
        )

        heavy_neighbors = sorted(
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

        nearest_distance, nearest_id, nearest_element = (
            heavy_neighbors[0]
        )

        second_distance = (
            heavy_neighbors[1][0]
            if len(heavy_neighbors) > 1
            else None
        )

        ownership_gate = (
            nearest_id == center_id
            and (
                second_distance is None
                or second_distance
                - nearest_distance
                >= CAP_NEAREST_CENTER_MARGIN_A
            )
        )

        if center["element"] == "B":
            bond_range_gate = (
                BH_MIN_A
                <= realized_length
                <= BH_MAX_A
            )
        elif center["element"] == "N":
            bond_range_gate = (
                NH_MIN_A
                <= realized_length
                <= NH_MAX_A
            )
        else:
            bond_range_gate = False

        cap_audit_records.append({
            "cap_id": cap_id,
            "center_atom": center_id,
            "center_element": center["element"],
            "target_length_A": target_length,
            "realized_length_A": realized_length,
            "nearest_heavy_atom": nearest_id,
            "nearest_heavy_element": nearest_element,
            "nearest_heavy_distance_A": nearest_distance,
            "second_heavy_distance_A": second_distance,
            "ownership_gate": ownership_gate,
            "bond_range_gate": bond_range_gate,
            "pass": (
                ownership_gate
                and bond_range_gate
                and abs(
                    realized_length
                    - target_length
                )
                <= 1.0e-8
            ),
        })

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

    cap_gate = all(
        row["pass"]
        for row in cap_audit_records
    )

    # S:1738 B-H2 geometry.
    s1738_caps = sorted(
        cap_id
        for cap_id, record in artificial_caps.items()
        if record["center_atom"] == "S:1738"
    )

    if len(s1738_caps) != 2:
        raise RuntimeError(
            "Expected exactly two caps on S:1738."
        )

    bh2_angle = angle(
        mapped[s1738_caps[0]]["xyz_A"],
        mapped["S:1738"]["xyz_A"],
        mapped[s1738_caps[1]]["xyz_A"],
    )

    bh2_geometry_gate = (
        BH2_MIN_ANGLE_DEG
        <= bh2_angle
        <= BH2_MAX_ANGLE_DEG
    )

    # Topological separations.
    graph_distance = {}

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

        for target, value in distances.items():
            graph_distance[
                canonical_pair(origin, target)
            ] = value

    contact_records = []
    hard_contacts = []

    ordered_atoms = sorted(
        mapped.items(),
        key=lambda item: item[1]["index"],
    )

    for first_position, (
        first_id,
        first,
    ) in enumerate(ordered_atoms):
        for second_id, second in ordered_atoms[
            first_position + 1:
        ]:
            value = distance(
                first["xyz_A"],
                second["xyz_A"],
            )

            separation = graph_distance.get(
                canonical_pair(first_id, second_id)
            )

            if separation in {1, 2}:
                classification = (
                    "TOPOLOGICAL_EXCLUSION"
                )
                is_hard = False
            else:
                first_element = first["element"]
                second_element = second["element"]

                if (
                    first_element == "H"
                    and second_element == "H"
                ):
                    threshold = HH_HARD_CONTACT_A
                elif (
                    first_element == "H"
                    or second_element == "H"
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

    # Constraint-design classification.
    proposed_fixed_ids = sorted(
        atom_id
        for atom_id, record in mapped.items()
        if (
            record["coordinate_source"]
            == "RETAINED_VALIDATED_V3A2_START"
            and record["v4_mobility_basis"]
            == "PENDING_PRE_QM_DESIGN"
            and (
                record["atom_role"]
                != "ARTIFICIAL_BOUNDARY_CAP"
            )
        )
    )

    restored_ids = sorted(
        atom_id
        for atom_id, record in mapped.items()
        if (
            record["atom_role"]
            == "RESTORED_CANONICAL_R2_ATOM"
        )
    )

    new_cap_ids = sorted(
        artificial_caps
    )

    constraint_design = {
        "fixed_core_candidates": proposed_fixed_ids,
        "restored_real_atom_candidates_for_mobility": restored_ids,
        "new_cap_candidates_for_mobility": new_cap_ids,
        "constraint_assignment_finalized": False,
    }

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
        "artificial_cap_geometry": cap_gate,
        "S1738_BH2_geometry": bh2_geometry_gate,
        "no_unresolved_hard_contacts": (
            hard_contact_gate
        ),
    }

    overall_pass = all(gates.values())

    decision = (
        "QM_F06_UPPER_V4_PRE_QM_STRUCTURAL_GATE_PASS_"
        "CONSTRAINT_DESIGN_AUTHORIZED"
        if overall_pass
        else
        "QM_F06_UPPER_V4_PRE_QM_STRUCTURAL_GATE_FAIL_"
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
        "S1738_BH2": {
            "caps": s1738_caps,
            "H_B_H_angle_deg": bh2_angle,
            "minimum_allowed_deg": (
                BH2_MIN_ANGLE_DEG
            ),
            "maximum_allowed_deg": (
                BH2_MAX_ANGLE_DEG
            ),
            "pass": bh2_geometry_gate,
        },
        "hard_contact_count": len(
            hard_contacts
        ),
        "hard_contacts": hard_contacts,
        "constraint_design": constraint_design,
        "files": {
            "valence_csv": str(
                VALENCE_CSV.relative_to(ROOT)
            ),
            "cap_audit_csv": str(
                CAP_AUDIT_CSV.relative_to(ROOT)
            ),
            "contacts_csv": str(
                CONTACTS_CSV.relative_to(ROOT)
            ),
        },
        "files_sha256": {
            "xyz": sha256(XYZ_PATH),
            "map": sha256(MAP_PATH),
            "caps": sha256(CAPS_PATH),
            "construction_report": sha256(
                CONSTRUCTION_REPORT
            ),
            "valence_csv": sha256(
                VALENCE_CSV
            ),
            "cap_audit_csv": sha256(
                CAP_AUDIT_CSV
            ),
            "contacts_csv": sha256(
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

    print("=" * 78)
    print("QM_F06 UPPER V4 PRE-QM STRUCTURAL AUDIT")
    print("=" * 78)

    for gate, passed in gates.items():
        print(
            f"{gate:40s}: "
            f"{'PASS' if passed else 'FAIL'}"
        )

    print()
    print("Composition:", dict(composition))
    print("Atom count:", len(atoms))
    print(
        "S:1738 H-B-H angle deg:",
        bh2_angle,
    )
    print(
        "Hard contacts:",
        len(hard_contacts),
    )
    print(
        "Valence failures:",
        sum(
            not row["pass"]
            for row in valence_records
        ),
    )
    print(
        "Cap failures:",
        sum(
            not row["pass"]
            for row in cap_audit_records
        ),
    )

    print()
    print("Decision:", decision)
    print("Report:", REPORT_PATH)
    print("Valence:", VALENCE_CSV)
    print("Caps:", CAP_AUDIT_CSV)
    print("Contacts:", CONTACTS_CSV)
    print(
        "Constraint design authorized:",
        overall_pass,
    )
    print("ORCA execution authorized: False")


if __name__ == "__main__":
    main()
