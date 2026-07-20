#!/usr/bin/env python3
"""
Pre-QM audit for QM_F06 UPPER Boundary V3-A2.

The audit verifies:
- expected composition B8 N7 H15;
- exact atom count;
- identity and mobility of P:1641, HCAP:UPPER:01 and HCAP:UPPER:04;
- repaired N-H distances, H-H distance and H-N-H angle;
- preservation of all prior V3-A constraints except indices 14 and 17;
- charge and multiplicity;
- ORCA resources;
- absence of GBW/MO reuse directives;
- absence of severe nonbonded contacts;
- topological exclusion of bonded and 1-3 pairs from generic hard-contact logic.

This script does not execute ORCA.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import Counter, deque
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

WORKFLOW = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V3/orca_v3a2_workflow"
)

XYZ_PATH = WORKFLOW / "v3a2_start.xyz"
INPUT_PATH = WORKFLOW / "v3a2.inp"
MAP_PATH = WORKFLOW / "v3a2_atom_role_constraint_map.csv"
SUMMARY_PATH = WORKFLOW / "v3a2_preparation_summary.json"

OUTPUT_DIR = WORKFLOW / "pre_qm_audit"
OUTPUT_JSON = OUTPUT_DIR / "QM_F06_UPPER_V3A2_PRE_QM_AUDIT.json"
OUTPUT_CONTACTS = OUTPUT_DIR / "QM_F06_UPPER_V3A2_contacts.csv"

EXPECTED_COMPOSITION = {
    "B": 8,
    "N": 7,
    "H": 15,
}

EXPECTED_ATOM_COUNT = 30
EXPECTED_FIXED_COUNT = 18
EXPECTED_MOBILE_COUNT = 12

CENTER_INDEX = 11
HEAVY_NEIGHBOR_INDEX = 13
CAP_INDICES = (14, 17)

EXPECTED_IDS = {
    11: "P:1641",
    13: "S:1710",
    14: "HCAP:UPPER:01",
    17: "HCAP:UPPER:04",
}

BOND_LIMITS = {
    frozenset(("B", "N")): (1.20, 1.85),
    frozenset(("B", "H")): (0.90, 1.45),
    frozenset(("N", "H")): (0.80, 1.30),
}

COVALENT_CONNECTIVITY_CUTOFFS = {
    frozenset(("B", "N")): 1.85,
    frozenset(("B", "H")): 1.45,
    frozenset(("N", "H")): 1.30,
}

HARD_CONTACT_CUTOFFS = {
    frozenset(("H", "H")): 1.20,
    frozenset(("B", "B")): 1.60,
    frozenset(("N", "N")): 1.60,
    frozenset(("B", "N")): 1.20,
    frozenset(("B", "H")): 0.90,
    frozenset(("N", "H")): 0.80,
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    return digest.hexdigest()


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def norm(vector):
    return math.sqrt(dot(vector, vector))


def distance(a, b):
    return norm(subtract(a, b))


def angle_deg(a, vertex, c):
    u = subtract(a, vertex)
    v = subtract(c, vertex)

    denominator = norm(u) * norm(v)

    if denominator < 1.0e-12:
        raise RuntimeError("Undefined angle from zero-length vector.")

    cosine = max(
        -1.0,
        min(1.0, dot(u, v) / denominator),
    )

    return math.degrees(math.acos(cosine))


def read_xyz(path: Path):
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    atom_count = int(lines[0].strip())
    atoms = []

    for index, line in enumerate(lines[2:2 + atom_count]):
        fields = line.split()

        if len(fields) < 4:
            raise RuntimeError(
                f"Malformed XYZ line for atom {index}: {line}"
            )

        atoms.append(
            {
                "index": index,
                "element": fields[0],
                "xyz": tuple(
                    float(value)
                    for value in fields[1:4]
                ),
            }
        )

    if len(atoms) != atom_count:
        raise RuntimeError(
            f"XYZ declares {atom_count} atoms but "
            f"contains {len(atoms)}."
        )

    return atoms


def read_map(path: Path):
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        headers = reader.fieldnames or []

    required = {
        "index_0based",
        "atom_id",
        "element",
        "v3a_fixed",
        "v3a2_fixed",
        "v3a2_mobile",
        "v3a2_mobility_basis",
    }

    missing = required - set(headers)

    if missing:
        raise RuntimeError(
            f"Constraint map missing columns: {sorted(missing)}"
        )

    return rows


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()

    if normalized == "true":
        return True

    if normalized == "false":
        return False

    raise RuntimeError(f"Invalid boolean value: {value!r}")


def extract_constraints(input_text: str) -> set[int]:
    return {
        int(value)
        for value in re.findall(
            r"\{\s*C\s+(\d+)\s+C\s*\}",
            input_text,
        )
    }


def graph_separation(
    graph: dict[int, set[int]],
    start: int,
    target: int,
    max_depth: int = 2,
):
    if start == target:
        return 0

    queue = deque([(start, 0)])
    visited = {start}

    while queue:
        node, depth = queue.popleft()

        if depth >= max_depth:
            continue

        for neighbor in graph[node]:
            if neighbor == target:
                return depth + 1

            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, depth + 1))

    return None


def build_geometry_graph(atoms):
    graph = {
        atom["index"]: set()
        for atom in atoms
    }

    edges = []

    for i, atom_i in enumerate(atoms):
        for j in range(i + 1, len(atoms)):
            atom_j = atoms[j]
            pair = frozenset(
                (atom_i["element"], atom_j["element"])
            )

            cutoff = COVALENT_CONNECTIVITY_CUTOFFS.get(pair)

            if cutoff is None:
                continue

            value = distance(
                atom_i["xyz"],
                atom_j["xyz"],
            )

            if value <= cutoff:
                graph[i].add(j)
                graph[j].add(i)
                edges.append(
                    {
                        "i": i,
                        "j": j,
                        "elements": sorted(pair),
                        "distance_A": value,
                    }
                )

    return graph, edges


def connected_components(graph):
    unseen = set(graph)
    components = []

    while unseen:
        seed = next(iter(unseen))
        stack = [seed]
        component = set()

        while stack:
            node = stack.pop()

            if node in component:
                continue

            component.add(node)
            unseen.discard(node)
            stack.extend(graph[node] - component)

        components.append(sorted(component))

    return components


def main() -> None:
    for path in (
        XYZ_PATH,
        INPUT_PATH,
        MAP_PATH,
        SUMMARY_PATH,
    ):
        require_file(path)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    atoms = read_xyz(XYZ_PATH)
    map_rows = read_map(MAP_PATH)
    input_text = INPUT_PATH.read_text(
        encoding="utf-8",
        errors="replace",
    )
    preparation_summary = json.loads(
        SUMMARY_PATH.read_text(encoding="utf-8")
    )

    map_by_index = {
        int(row["index_0based"]): row
        for row in map_rows
    }

    gates = {}

    composition = Counter(
        atom["element"]
        for atom in atoms
    )

    gates["atom_count"] = (
        len(atoms) == EXPECTED_ATOM_COUNT
    )
    gates["composition"] = (
        dict(composition) == EXPECTED_COMPOSITION
    )

    identity_checks = {}

    for index, expected_id in EXPECTED_IDS.items():
        actual_id = map_by_index[index]["atom_id"]
        identity_checks[str(index)] = {
            "expected": expected_id,
            "actual": actual_id,
            "pass": actual_id == expected_id,
        }

    gates["critical_atom_identity"] = all(
        item["pass"]
        for item in identity_checks.values()
    )

    v3a_fixed = {
        int(row["index_0based"])
        for row in map_rows
        if parse_bool(row["v3a_fixed"])
    }

    v3a2_fixed_from_map = {
        int(row["index_0based"])
        for row in map_rows
        if parse_bool(row["v3a2_fixed"])
    }

    v3a2_mobile_from_map = {
        int(row["index_0based"])
        for row in map_rows
        if parse_bool(row["v3a2_mobile"])
    }

    fixed_from_input = extract_constraints(input_text)

    expected_v3a2_fixed = (
        v3a_fixed - set(CAP_INDICES)
    )

    gates["only_caps_14_17_released"] = (
        v3a2_fixed_from_map == expected_v3a2_fixed
        and fixed_from_input == expected_v3a2_fixed
    )

    gates["center_remains_fixed"] = (
        CENTER_INDEX in v3a2_fixed_from_map
        and CENTER_INDEX in fixed_from_input
    )

    gates["caps_are_mobile"] = all(
        index in v3a2_mobile_from_map
        and index not in v3a2_fixed_from_map
        and index not in fixed_from_input
        for index in CAP_INDICES
    )

    gates["fixed_count"] = (
        len(v3a2_fixed_from_map)
        == EXPECTED_FIXED_COUNT
        and len(fixed_from_input)
        == EXPECTED_FIXED_COUNT
    )

    gates["mobile_count"] = (
        len(v3a2_mobile_from_map)
        == EXPECTED_MOBILE_COUNT
    )

    center = atoms[CENTER_INDEX]["xyz"]
    heavy = atoms[HEAVY_NEIGHBOR_INDEX]["xyz"]
    h1 = atoms[CAP_INDICES[0]]["xyz"]
    h2 = atoms[CAP_INDICES[1]]["xyz"]

    local_metrics = {
        "N_H14_A": distance(center, h1),
        "N_H17_A": distance(center, h2),
        "H14_H17_A": distance(h1, h2),
        "H14_N11_H17_deg": angle_deg(
            h1,
            center,
            h2,
        ),
        "N11_B13_A": distance(center, heavy),
    }

    gates["local_NH_distances"] = (
        0.98 <= local_metrics["N_H14_A"] <= 1.04
        and 0.98 <= local_metrics["N_H17_A"] <= 1.04
    )

    gates["geminal_HH_resolved"] = (
        local_metrics["H14_H17_A"] >= 1.50
    )

    gates["local_HNH_angle"] = (
        110.0
        <= local_metrics["H14_N11_H17_deg"]
        <= 130.0
    )

    gates["center_heavy_bond"] = (
        1.20 <= local_metrics["N11_B13_A"] <= 1.85
    )

    graph, inferred_edges = build_geometry_graph(atoms)
    components = connected_components(graph)

    gates["single_connected_component"] = (
        len(components) == 1
    )

    contacts = []
    hard_contacts = []

    for i, atom_i in enumerate(atoms):
        for j in range(i + 1, len(atoms)):
            atom_j = atoms[j]

            pair = frozenset(
                (atom_i["element"], atom_j["element"])
            )

            cutoff = HARD_CONTACT_CUTOFFS.get(pair)

            if cutoff is None:
                continue

            value = distance(
                atom_i["xyz"],
                atom_j["xyz"],
            )

            separation = graph_separation(
                graph,
                i,
                j,
                max_depth=2,
            )

            excluded_topologically = (
                separation in {1, 2}
            )

            is_hard = (
                value < cutoff
                and not excluded_topologically
            )

            if value < cutoff + 0.50:
                record = {
                    "index_1": i,
                    "atom_id_1": map_by_index[i]["atom_id"],
                    "element_1": atom_i["element"],
                    "index_2": j,
                    "atom_id_2": map_by_index[j]["atom_id"],
                    "element_2": atom_j["element"],
                    "distance_A": value,
                    "hard_cutoff_A": cutoff,
                    "graph_separation": separation,
                    "topologically_excluded": excluded_topologically,
                    "hard_contact": is_hard,
                }

                contacts.append(record)

                if is_hard:
                    hard_contacts.append(record)

    gates["no_unresolved_hard_contacts"] = (
        len(hard_contacts) == 0
    )

    charge_mult_match = re.search(
        r"(?im)^\s*\*\s+xyz\s+(-?\d+)\s+(\d+)\s*$",
        input_text,
    )

    charge = (
        int(charge_mult_match.group(1))
        if charge_mult_match
        else None
    )

    multiplicity = (
        int(charge_mult_match.group(2))
        if charge_mult_match
        else None
    )

    gates["charge_multiplicity"] = (
        charge == 0
        and multiplicity == 1
    )

    nprocs_match = re.search(
        r"(?im)^\s*nprocs\s+(\d+)\s*$",
        input_text,
    )

    maxcore_match = re.search(
        r"(?im)^\s*%maxcore\s+(\d+)\s*$",
        input_text,
    )

    nprocs = (
        int(nprocs_match.group(1))
        if nprocs_match
        else None
    )

    maxcore = (
        int(maxcore_match.group(1))
        if maxcore_match
        else None
    )

    gates["resources"] = (
        nprocs == 4
        and maxcore == 2500
    )

    forbidden_reuse_patterns = {
        "moread": r"(?i)\bmoread\b",
        "moinp": r"(?i)\bmoinp\b",
        "gbw_reference": r"(?i)\.gbw\b",
    }

    reuse_hits = {
        name: bool(re.search(pattern, input_text))
        for name, pattern
        in forbidden_reuse_patterns.items()
    }

    gates["fresh_scf"] = not any(
        reuse_hits.values()
    )

    gates["input_has_geometry"] = (
        len(
            re.findall(
                r"(?m)^\s*(?:B|N|H)\s+"
                r"[-+]?\d+\.\d+\s+"
                r"[-+]?\d+\.\d+\s+"
                r"[-+]?\d+\.\d+\s*$",
                input_text,
            )
        )
        == EXPECTED_ATOM_COUNT
    )

    summary_hh_distance = preparation_summary[
        "repaired_metrics"
    ]["H1_H2_A"]

    hh_consistency_tolerance_A = 1.0e-8

    gates["preparation_summary_consistent"] = (
        abs(
            summary_hh_distance
            - local_metrics["H14_H17_A"]
        )
        <= hh_consistency_tolerance_A
        and preparation_summary[
            "v3a2_fixed_count"
        ]
        == EXPECTED_FIXED_COUNT
        and preparation_summary[
            "v3a2_mobile_count"
        ]
        == EXPECTED_MOBILE_COUNT
    )

    overall_pass = all(gates.values())

    decision = (
        "QM_F06_UPPER_V3A2_PRE_QM_GATE_PASS_"
        "ORCA_EXECUTION_AUTHORIZED"
        if overall_pass
        else
        "QM_F06_UPPER_V3A2_PRE_QM_GATE_FAIL_"
        "ORCA_EXECUTION_BLOCKED"
    )

    contact_fields = [
        "index_1",
        "atom_id_1",
        "element_1",
        "index_2",
        "atom_id_2",
        "element_2",
        "distance_A",
        "hard_cutoff_A",
        "graph_separation",
        "topologically_excluded",
        "hard_contact",
    ]

    with OUTPUT_CONTACTS.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=contact_fields,
        )
        writer.writeheader()
        writer.writerows(contacts)

    report = {
        "decision": decision,
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "overall_pass": overall_pass,
        "gates": gates,
        "composition": dict(composition),
        "atom_count": len(atoms),
        "critical_atom_identity": identity_checks,
        "local_metrics": local_metrics,
        "constraints": {
            "v3a_fixed_count": len(v3a_fixed),
            "expected_v3a2_fixed": sorted(
                expected_v3a2_fixed
            ),
            "v3a2_fixed_from_map": sorted(
                v3a2_fixed_from_map
            ),
            "v3a2_fixed_from_input": sorted(
                fixed_from_input
            ),
            "v3a2_mobile_from_map": sorted(
                v3a2_mobile_from_map
            ),
        },
        "graph": {
            "inferred_edge_count": len(inferred_edges),
            "component_count": len(components),
            "components": components,
        },
        "contacts": {
            "reviewed_contact_count": len(contacts),
            "hard_contact_count": len(hard_contacts),
            "hard_contacts": hard_contacts,
        },
        "orca_input": {
            "charge": charge,
            "multiplicity": multiplicity,
            "nprocs": nprocs,
            "maxcore_mb_per_process": maxcore,
            "reuse_hits": reuse_hits,
        },
        "files_sha256": {
            "xyz": sha256(XYZ_PATH),
            "input": sha256(INPUT_PATH),
            "constraint_map": sha256(MAP_PATH),
            "preparation_summary": sha256(
                SUMMARY_PATH
            ),
            "contacts_csv": sha256(
                OUTPUT_CONTACTS
            ),
        },
        "authorization": {
            "orca_execution_authorized": overall_pass,
            "geometry_reference_accepted": False,
            "electronic_reference_accepted": False,
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    OUTPUT_JSON.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 78)
    print("QM_F06 UPPER V3-A2 PRE-QM AUDIT")
    print("=" * 78)

    for gate, passed in gates.items():
        print(
            f"{gate:40s}: "
            f"{'PASS' if passed else 'FAIL'}"
        )

    print()
    print("Composition:", dict(composition))
    print("Atom count:", len(atoms))
    print("Fixed atoms:", len(v3a2_fixed_from_map))
    print("Mobile atoms:", len(v3a2_mobile_from_map))
    print("N-H14:", local_metrics["N_H14_A"])
    print("N-H17:", local_metrics["N_H17_A"])
    print("H14-H17:", local_metrics["H14_H17_A"])
    print(
        "H14-N11-H17:",
        local_metrics["H14_N11_H17_deg"],
    )
    print("Hard contacts:", len(hard_contacts))
    print()
    print("Decision:", decision)
    print("Report:", OUTPUT_JSON)
    print("Contacts:", OUTPUT_CONTACTS)


if __name__ == "__main__":
    main()
