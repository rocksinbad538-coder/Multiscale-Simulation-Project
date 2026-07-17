#!/usr/bin/env python3
"""
Post-optimization validation of QM_F06 LOWER Boundary V2-A.

Compares the repaired V2-A input geometry with the converged ORCA V2-A
geometry and evaluates:

- ORCA workflow state and optimized XYZ;
- immobility of the 21 constrained atoms;
- displacement of the seven mobile atoms;
- preservation of all 27 intended bonds;
- graph degree and connectivity;
- topology-aware nonbonded contacts;
- artificial-cap and cap–bridge contacts;
- possible unintended covalent contacts.

No coordinates are modified and no QM calculation is executed.
"""

from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

V2_DIR = ROOT / (
    "runs/phase1A/day027_qm_f06_lower_boundary_redesign/"
    "QM_F06_LOWER_BOUNDARY_V2"
)

WORKFLOW_DIR = V2_DIR / "orca_v2_workflow"
STATE_PATH = WORKFLOW_DIR / "v2_workflow_state.json"

INITIAL_ATOMS_PATH = V2_DIR / (
    "cap02_repair/"
    "QM_F06_LOWER_BOUNDARY_V2_REPAIRED_atoms.csv"
)

FULL_EDGES_PATH = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_edges.csv"
)

OUTPUT_DIR = V2_DIR / "v2a_postprocessing"

BOND_RANGES = {
    tuple(sorted(("B", "N"))): (1.20, 1.85),
    tuple(sorted(("B", "H"))): (0.90, 1.45),
    tuple(sorted(("N", "H"))): (0.80, 1.30),
}

EXPECTED_DEGREES = {"B": 3, "N": 3, "H": 1}

VDW_RADII = {"H": 1.20, "B": 1.92, "N": 1.55}

COVALENT_THRESHOLDS = {
    tuple(sorted(("B", "B"))): 1.90,
    tuple(sorted(("B", "N"))): 1.85,
    tuple(sorted(("N", "N"))): 1.70,
    tuple(sorted(("B", "H"))): 1.45,
    tuple(sorted(("N", "H"))): 1.30,
    tuple(sorted(("H", "H"))): 0.90,
}

BRIDGE_ATOMS = {
    "BR4:LOWER:00:1",
    "BR4:LOWER:00:2",
    "BR4:LOWER:00:3",
    "BR4:LOWER:00:4",
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"No rows in {path}")
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return

    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(
            {field: row.get(field, "") for field in fields}
            for row in rows
        )


def xyz_from_row(row: dict[str, str]) -> tuple[float, float, float]:
    return (
        float(row["x_angstrom"]),
        float(row["y_angstrom"]),
        float(row["z_angstrom"]),
    )


def read_xyz(path: Path) -> list[tuple[str, float, float, float]]:
    require_file(path)
    lines = path.read_text(encoding="utf-8").splitlines()
    expected = int(lines[0].strip())

    rows = [
        (
            parts[0],
            float(parts[1]),
            float(parts[2]),
            float(parts[3]),
        )
        for line in lines[2:]
        if line.strip()
        for parts in [line.split()]
    ]

    if len(rows) != expected:
        raise RuntimeError(
            f"XYZ count mismatch in {path}: {expected} vs {len(rows)}"
        )

    return rows


def distance(a, b) -> float:
    return math.sqrt(
        sum((x - y) ** 2 for x, y in zip(a, b, strict=True))
    )


def canonical_edge(a: str, b: str) -> tuple[str, str]:
    return tuple(sorted((a, b)))


def shortest_path_length(adjacency, source, target):
    if source == target:
        return 0

    visited = {source}
    queue = deque([(source, 0)])

    while queue:
        node, depth = queue.popleft()

        for neighbor in adjacency[node]:
            if neighbor == target:
                return depth + 1
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, depth + 1))

    return None


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    require_file(STATE_PATH)
    state = json.loads(STATE_PATH.read_text(encoding="utf-8"))

    required_state = {
        "v2a_executed": state.get("v2a_executed") is True,
        "v2a_validation_pass": state.get("v2a_validation_pass") is True,
        "v2a_optimized_xyz_present": bool(state.get("v2a_optimized_xyz")),
    }

    if not all(required_state.values()):
        raise RuntimeError(f"Incomplete V2-A workflow state: {required_state}")

    optimized_xyz_path = ROOT / state["v2a_optimized_xyz"]
    require_file(optimized_xyz_path)

    atom_rows = read_csv(INITIAL_ATOMS_PATH)
    optimized_xyz = read_xyz(optimized_xyz_path)

    if len(atom_rows) != 28 or len(optimized_xyz) != 28:
        raise RuntimeError("Expected 28 atoms in manifest and optimized XYZ.")

    atom_ids = [row["atom_id"] for row in atom_rows]
    atom_set = set(atom_ids)

    elements = {row["atom_id"]: row["element"] for row in atom_rows}
    roles = {row["atom_id"]: row["atom_role"] for row in atom_rows}

    initial_coords = {
        row["atom_id"]: xyz_from_row(row)
        for row in atom_rows
    }

    optimized_coords = {}

    for atom_id, manifest_row, xyz_row in zip(
        atom_ids, atom_rows, optimized_xyz, strict=True
    ):
        if manifest_row["element"] != xyz_row[0]:
            raise RuntimeError(f"Element mismatch for {atom_id}")
        optimized_coords[atom_id] = xyz_row[1:]

    artificial_caps = {
        row["atom_id"]
        for row in atom_rows
        if row["artificial_cap"].lower() == "true"
    }

    edge_records = {}

    for row in read_csv(FULL_EDGES_PATH):
        a = row["source_node"]
        b = row["target_node"]

        if a in atom_set and b in atom_set:
            edge_records[canonical_edge(a, b)] = {
                "origin": "REAL_R2_GRAPH_EDGE",
                "source_edge_id": row["edge_id"],
                "edge_type": row["edge_type"],
            }

    for row in atom_rows:
        if row["artificial_cap"].lower() == "true":
            cap = row["atom_id"]
            parent = row["parent_inside_node"]

            if parent not in atom_set:
                raise RuntimeError(f"Missing cap parent: {cap} -> {parent}")

            edge_records[canonical_edge(cap, parent)] = {
                "origin": "ARTIFICIAL_CAP_EDGE",
                "source_edge_id": row["source_edge_id"],
                "edge_type": "QM_BOUNDARY_CAP",
            }

    adjacency = defaultdict(set)

    for a, b in edge_records:
        adjacency[a].add(b)
        adjacency[b].add(a)

    displacement_rows = []
    fixed_motion_failures = 0
    fixed_tolerance = 1.0e-6

    for atom_id in atom_ids:
        displacement = distance(
            initial_coords[atom_id],
            optimized_coords[atom_id],
        )

        fixed = roles[atom_id] not in {
            "REAL_R2_BOUNDARY_EXPANSION_ATOM",
            "ARTIFICIAL_BOUNDARY_CAP_V2",
        }

        fixed_pass = (not fixed) or displacement <= fixed_tolerance

        if not fixed_pass:
            fixed_motion_failures += 1

        displacement_rows.append(
            {
                "atom_id": atom_id,
                "element": elements[atom_id],
                "atom_role": roles[atom_id],
                "expected_fixed": fixed,
                "displacement_angstrom": f"{displacement:.12f}",
                "fixed_atom_motion_gate_pass": fixed_pass,
            }
        )

    bond_rows = []
    bond_failures = 0
    cap_bond_failures = 0

    for edge in sorted(edge_records):
        a, b = edge
        pair = tuple(sorted((elements[a], elements[b])))
        limits = BOND_RANGES[pair]

        measured = distance(
            optimized_coords[a],
            optimized_coords[b],
        )

        passed = limits[0] <= measured <= limits[1]
        is_cap = edge_records[edge]["origin"] == "ARTIFICIAL_CAP_EDGE"

        if not passed:
            bond_failures += 1
        if is_cap and not passed:
            cap_bond_failures += 1

        bond_rows.append(
            {
                "atom_1": a,
                "element_1": elements[a],
                "atom_2": b,
                "element_2": elements[b],
                "edge_origin": edge_records[edge]["origin"],
                "distance_angstrom": f"{measured:.10f}",
                "allowed_min_angstrom": limits[0],
                "allowed_max_angstrom": limits[1],
                "artificial_cap_edge": is_cap,
                "bond_gate_pass": passed,
            }
        )

    valence_rows = []
    valence_failures = 0

    for atom_id in atom_ids:
        observed = len(adjacency[atom_id])
        expected = EXPECTED_DEGREES[elements[atom_id]]
        passed = observed == expected

        if not passed:
            valence_failures += 1

        valence_rows.append(
            {
                "atom_id": atom_id,
                "element": elements[atom_id],
                "observed_degree": observed,
                "expected_degree": expected,
                "neighbors": "|".join(sorted(adjacency[atom_id])),
                "degree_gate_pass": passed,
            }
        )

    close_rows = []
    hard_contacts = 0
    hard_cap_contacts = 0
    bridge_cap_hard_contacts = 0
    unintended_covalent_contacts = 0

    for i, a in enumerate(atom_ids):
        for b in atom_ids[i + 1:]:
            separation = shortest_path_length(adjacency, a, b)

            if separation in {1, 2, 3}:
                continue

            measured = distance(
                optimized_coords[a],
                optimized_coords[b],
            )

            pair = tuple(sorted((elements[a], elements[b])))
            vdw_sum = VDW_RADII[elements[a]] + VDW_RADII[elements[b]]
            ratio = measured / vdw_sum

            involves_cap = a in artificial_caps or b in artificial_caps
            involves_bridge = a in BRIDGE_ATOMS or b in BRIDGE_ATOMS
            possible_new_bond = measured <= COVALENT_THRESHOLDS[pair]
            hard = ratio < 0.70

            if possible_new_bond:
                unintended_covalent_contacts += 1
            if hard:
                hard_contacts += 1
            if hard and involves_cap:
                hard_cap_contacts += 1
            if hard and involves_cap and involves_bridge:
                bridge_cap_hard_contacts += 1

            if ratio >= 0.90 and not possible_new_bond:
                continue

            close_rows.append(
                {
                    "atom_1": a,
                    "element_1": elements[a],
                    "atom_2": b,
                    "element_2": elements[b],
                    "graph_separation": (
                        separation if separation is not None else "DISCONNECTED"
                    ),
                    "distance_angstrom": f"{measured:.10f}",
                    "distance_over_vdw_sum": f"{ratio:.10f}",
                    "involves_artificial_cap": involves_cap,
                    "involves_bridge": involves_bridge,
                    "hard_contact_below_0p70": hard,
                    "possible_unintended_covalent_contact": possible_new_bond,
                }
            )

    visited = {atom_ids[0]}
    queue = deque([atom_ids[0]])

    while queue:
        node = queue.popleft()
        for neighbor in adjacency[node]:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)

    graph_connected = len(visited) == len(atom_ids)

    gate_pass = all(
        (
            graph_connected,
            fixed_motion_failures == 0,
            valence_failures == 0,
            bond_failures == 0,
            cap_bond_failures == 0,
            hard_cap_contacts == 0,
            bridge_cap_hard_contacts == 0,
            unintended_covalent_contacts == 0,
        )
    )

    decision = (
        "QM_F06_LOWER_BOUNDARY_V2A_OPTIMIZED_GEOMETRY_PASS"
        if gate_pass
        else "QM_F06_LOWER_BOUNDARY_V2A_OPTIMIZED_GEOMETRY_FAIL"
    )

    write_csv(
        OUTPUT_DIR / "QM_F06_LOWER_BOUNDARY_V2A_displacements.csv",
        displacement_rows,
    )
    write_csv(
        OUTPUT_DIR / "QM_F06_LOWER_BOUNDARY_V2A_bond_audit.csv",
        bond_rows,
    )
    write_csv(
        OUTPUT_DIR / "QM_F06_LOWER_BOUNDARY_V2A_valence_audit.csv",
        valence_rows,
    )
    write_csv(
        OUTPUT_DIR / "QM_F06_LOWER_BOUNDARY_V2A_nonbonded_audit.csv",
        close_rows,
    )

    mobile_rows = [
        row for row in displacement_rows
        if row["expected_fixed"] is False
    ]

    max_mobile = max(
        mobile_rows,
        key=lambda row: float(row["displacement_angstrom"]),
    )

    report_path = (
        OUTPUT_DIR / "QM_F06_LOWER_BOUNDARY_V2A_VALIDATION_REPORT.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 LOWER Boundary V2-A Validation — Day028",
                "",
                f"## Decision: **{decision}**",
                "",
                "## ORCA execution",
                "",
                "- Normal termination: **YES**",
                "- Geometry convergence: **YES**",
                f"- Optimized XYZ: `{optimized_xyz_path.relative_to(ROOT)}`",
                "",
                "## Constraint integrity",
                "",
                "- Expected fixed atoms: **21**",
                "- Expected mobile atoms: **7**",
                f"- Fixed-atom motion failures: **{fixed_motion_failures}**",
                (
                    "- Maximum mobile displacement: "
                    f"**{float(max_mobile['displacement_angstrom']):.6f} Å** "
                    f"for `{max_mobile['atom_id']}`"
                ),
                "",
                "## Connectivity and bonding",
                "",
                f"- Graph connected: **{graph_connected}**",
                f"- Intended edges: **{len(edge_records)}**",
                f"- Valence failures: **{valence_failures}**",
                f"- Bond-range failures: **{bond_failures}**",
                f"- Artificial-cap bond failures: **{cap_bond_failures}**",
                "",
                "## Topology-aware contacts",
                "",
                f"- Reported close contacts: **{len(close_rows)}**",
                f"- Hard contacts below 0.70: **{hard_contacts}**",
                f"- Hard cap contacts: **{hard_cap_contacts}**",
                f"- Hard bridge–cap contacts: **{bridge_cap_hard_contacts}**",
                (
                    "- Possible unintended covalent contacts: "
                    f"**{unintended_covalent_contacts}**"
                ),
                "",
                "## Authorization state",
                "",
                (
                    "- V2-A optimized structural gate: "
                    f"**{'PASS' if gate_pass else 'FAIL'}**"
                ),
                "- V2-B preparation: **PENDING THIS RESULT**",
                "- Final ESP/RESP calculation: **NOT AUTHORIZED**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = {
        "decision": decision,
        "optimized_xyz": str(optimized_xyz_path.relative_to(ROOT)),
        "fixed_atom_motion_failures": fixed_motion_failures,
        "graph_connected": graph_connected,
        "valence_failures": valence_failures,
        "bond_failures": bond_failures,
        "cap_bond_failures": cap_bond_failures,
        "hard_contacts": hard_contacts,
        "hard_cap_contacts": hard_cap_contacts,
        "bridge_cap_hard_contacts": bridge_cap_hard_contacts,
        "unintended_covalent_contacts": unintended_covalent_contacts,
        "v2a_structural_gate_pass": gate_pass,
        "v2b_preparation_authorized": gate_pass,
        "qm_executed_by_this_script": False,
    }

    (
        OUTPUT_DIR / "QM_F06_LOWER_BOUNDARY_V2A_validation_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("QM_F06 LOWER Boundary V2-A optimized-geometry audit completed.")
    print("Decision:", decision)
    print("Fixed-atom motion failures:", fixed_motion_failures)
    print("Bond failures:", bond_failures)
    print("Hard cap contacts:", hard_cap_contacts)
    print("Bridge-cap hard contacts:", bridge_cap_hard_contacts)
    print("Unintended covalent contacts:", unintended_covalent_contacts)
    print("V2-B preparation authorized:", gate_pass)
    print("Report:", report_path)


if __name__ == "__main__":
    main()
