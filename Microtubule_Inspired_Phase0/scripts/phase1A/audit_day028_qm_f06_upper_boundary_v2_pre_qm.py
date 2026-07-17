#!/usr/bin/env python3
"""
Pre-QM audit of QM_F06 UPPER Boundary V2.

No geometry is modified and no QM calculation is executed.
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
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V2"
)

ATOMS_PATH = V2_DIR / "QM_F06_UPPER_BOUNDARY_V2_atoms.csv"
CAPS_PATH = V2_DIR / "QM_F06_UPPER_BOUNDARY_V2_caps.csv"

FULL_EDGES_PATH = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_edges.csv"
)

OUTPUT_DIR = V2_DIR / "pre_qm_audit"

EXPECTED_DEGREES = {
    "B": 3,
    "N": 3,
    "H": 1,
}

BOND_RANGES = {
    tuple(sorted(("B", "N"))): (1.20, 1.85),
    tuple(sorted(("B", "H"))): (0.90, 1.45),
    tuple(sorted(("N", "H"))): (0.80, 1.30),
}

VDW_RADII = {
    "H": 1.20,
    "B": 1.92,
    "N": 1.55,
}

COVALENT_THRESHOLDS = {
    tuple(sorted(("B", "B"))): 1.90,
    tuple(sorted(("B", "N"))): 1.85,
    tuple(sorted(("N", "N"))): 1.70,
    tuple(sorted(("B", "H"))): 1.45,
    tuple(sorted(("N", "H"))): 1.30,
    tuple(sorted(("H", "H"))): 0.90,
}

BRIDGE_ATOMS = {
    "BR4:UPPER:00:1",
    "BR4:UPPER:00:2",
    "BR4:UPPER:00:3",
    "BR4:UPPER:00:4",
}

RESTORED_ATOMS = {
    "A:UPPER:11:3",
    "A:UPPER:13:3",
    "A:UPPER:14:2",
    "H4:UPPER:0046:0",
}

CRITICAL_SEAM_EDGE = tuple(
    sorted(("A:UPPER:14:4", "A:UPPER:13:3"))
)


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(f"No records in {path}")

    return rows


def write_csv(
    path: Path,
    rows: list[dict[str, Any]],
) -> None:
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

        for row in rows:
            writer.writerow(
                {
                    field: row.get(field, "")
                    for field in fields
                }
            )


def distance(first, second) -> float:
    return math.sqrt(
        sum(
            (a - b) ** 2
            for a, b in zip(first, second, strict=True)
        )
    )


def canonical_edge(
    first: str,
    second: str,
) -> tuple[str, str]:
    return tuple(sorted((first, second)))


def shortest_path_length(
    adjacency: dict[str, set[str]],
    source: str,
    target: str,
) -> int | None:
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

    atom_rows = read_csv(ATOMS_PATH)
    cap_rows = read_csv(CAPS_PATH)
    full_edge_rows = read_csv(FULL_EDGES_PATH)

    if len(atom_rows) != 28:
        raise RuntimeError(
            f"Expected 28 atoms; found {len(atom_rows)}"
        )

    atom_ids = [row["atom_id"] for row in atom_rows]

    if len(atom_ids) != len(set(atom_ids)):
        raise RuntimeError("Duplicate atom IDs detected.")

    atom_set = set(atom_ids)

    elements = {
        row["atom_id"]: row["element"]
        for row in atom_rows
    }

    roles = {
        row["atom_id"]: row["atom_role"]
        for row in atom_rows
    }

    coords = {
        row["atom_id"]: (
            float(row["x_angstrom"]),
            float(row["y_angstrom"]),
            float(row["z_angstrom"]),
        )
        for row in atom_rows
    }

    artificial_caps = {
        row["atom_id"]
        for row in atom_rows
        if row["artificial_cap"].lower() == "true"
    }

    edge_records: dict[
        tuple[str, str],
        dict[str, Any],
    ] = {}

    for row in full_edge_rows:
        first = row["source_node"]
        second = row["target_node"]

        if first in atom_set and second in atom_set:
            edge_records[
                canonical_edge(first, second)
            ] = {
                "edge_origin": "REAL_R2_GRAPH_EDGE",
                "source_edge_id": row["edge_id"],
                "edge_type": row["edge_type"],
            }

    for row in atom_rows:
        if row["artificial_cap"].lower() != "true":
            continue

        cap = row["atom_id"]
        parent = row["parent_inside_node"]

        if parent not in atom_set:
            raise RuntimeError(
                f"Artificial-cap parent missing: {cap} -> {parent}"
            )

        edge_records[
            canonical_edge(cap, parent)
        ] = {
            "edge_origin": "ARTIFICIAL_CAP_EDGE",
            "source_edge_id": row["source_edge_id"],
            "edge_type": "QM_BOUNDARY_CAP",
        }

    adjacency: dict[str, set[str]] = defaultdict(set)

    for first, second in edge_records:
        adjacency[first].add(second)
        adjacency[second].add(first)

    # Connectivity
    visited = {atom_ids[0]}
    queue = deque([atom_ids[0]])

    while queue:
        node = queue.popleft()

        for neighbor in adjacency[node]:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)

    graph_connected = len(visited) == len(atom_ids)

    # Valence
    valence_rows = []
    valence_failures = 0

    for atom_id in atom_ids:
        expected = EXPECTED_DEGREES[elements[atom_id]]
        observed = len(adjacency[atom_id])
        passed = observed == expected

        if not passed:
            valence_failures += 1

        valence_rows.append(
            {
                "atom_id": atom_id,
                "element": elements[atom_id],
                "atom_role": roles[atom_id],
                "expected_degree": expected,
                "observed_degree": observed,
                "neighbors": "|".join(
                    sorted(adjacency[atom_id])
                ),
                "degree_gate_pass": passed,
            }
        )

    # Bond geometry
    bond_rows = []
    bond_failures = 0
    cap_bond_failures = 0
    seam_bond_failures = 0

    for edge in sorted(edge_records):
        first, second = edge

        pair = tuple(
            sorted((elements[first], elements[second]))
        )

        if pair not in BOND_RANGES:
            raise RuntimeError(
                f"Unexpected bonded pair: {first} — {second}: {pair}"
            )

        lower, upper = BOND_RANGES[pair]
        measured = distance(coords[first], coords[second])
        passed = lower <= measured <= upper

        is_cap_edge = (
            edge_records[edge]["edge_origin"]
            == "ARTIFICIAL_CAP_EDGE"
        )

        touches_restored_region = (
            first in RESTORED_ATOMS
            or second in RESTORED_ATOMS
        )

        crosses_coordinate_seam = (
            touches_restored_region
            and not (
                first in RESTORED_ATOMS
                and second in RESTORED_ATOMS
            )
        )

        critical_restored_edge = edge == CRITICAL_SEAM_EDGE

        if not passed:
            bond_failures += 1

        if is_cap_edge and not passed:
            cap_bond_failures += 1

        if crosses_coordinate_seam and not passed:
            seam_bond_failures += 1

        bond_rows.append(
            {
                "atom_1": first,
                "element_1": elements[first],
                "atom_2": second,
                "element_2": elements[second],
                "edge_origin": edge_records[edge]["edge_origin"],
                "source_edge_id": edge_records[edge]["source_edge_id"],
                "distance_angstrom": f"{measured:.10f}",
                "allowed_min_angstrom": lower,
                "allowed_max_angstrom": upper,
                "artificial_cap_edge": is_cap_edge,
                "touches_restored_region": touches_restored_region,
                "crosses_coordinate_seam": crosses_coordinate_seam,
                "critical_restored_edge": critical_restored_edge,
                "bond_distance_gate_pass": passed,
            }
        )

    # Topology-aware nonbonded audit
    nonbonded_rows = []
    hard_contacts = 0
    hard_cap_contacts = 0
    bridge_cap_hard_contacts = 0
    unintended_covalent_contacts = 0

    for index, first in enumerate(atom_ids):
        for second in atom_ids[index + 1:]:
            separation = shortest_path_length(
                adjacency,
                first,
                second,
            )

            if separation in {1, 2, 3}:
                continue

            measured = distance(
                coords[first],
                coords[second],
            )

            vdw_sum = (
                VDW_RADII[elements[first]]
                + VDW_RADII[elements[second]]
            )

            ratio = measured / vdw_sum

            pair = tuple(
                sorted((elements[first], elements[second]))
            )

            possible_new_bond = (
                measured <= COVALENT_THRESHOLDS[pair]
            )

            involves_cap = (
                first in artificial_caps
                or second in artificial_caps
            )

            involves_bridge = (
                first in BRIDGE_ATOMS
                or second in BRIDGE_ATOMS
            )

            hard = ratio < 0.70

            if hard:
                hard_contacts += 1

            if hard and involves_cap:
                hard_cap_contacts += 1

            if hard and involves_cap and involves_bridge:
                bridge_cap_hard_contacts += 1

            if possible_new_bond:
                unintended_covalent_contacts += 1

            if ratio >= 0.90 and not possible_new_bond:
                continue

            nonbonded_rows.append(
                {
                    "atom_1": first,
                    "element_1": elements[first],
                    "atom_2": second,
                    "element_2": elements[second],
                    "graph_separation": (
                        separation
                        if separation is not None
                        else "DISCONNECTED"
                    ),
                    "distance_angstrom": f"{measured:.10f}",
                    "distance_over_vdw_sum": f"{ratio:.10f}",
                    "involves_artificial_cap": involves_cap,
                    "involves_bridge": involves_bridge,
                    "hard_contact_below_0p70": hard,
                    "possible_unintended_covalent_contact": (
                        possible_new_bond
                    ),
                }
            )

    element_counts = Counter(elements.values())

    formula = "".join(
        f"{element}{element_counts[element]}"
        for element in ("B", "N", "H")
        if element_counts[element]
    )

    neutral_valence_electrons = (
        3 * element_counts["B"]
        + 5 * element_counts["N"]
        + element_counts["H"]
    )

    provisional_charge = 0
    provisional_multiplicity = (
        1 if neutral_valence_electrons % 2 == 0 else 2
    )

    gate_pass = all(
        (
            graph_connected,
            len(edge_records) == 27,
            valence_failures == 0,
            bond_failures == 0,
            seam_bond_failures == 0,
            cap_bond_failures == 0,
            hard_cap_contacts == 0,
            bridge_cap_hard_contacts == 0,
            unintended_covalent_contacts == 0,
            neutral_valence_electrons % 2 == 0,
        )
    )

    decision = (
        "QM_F06_UPPER_BOUNDARY_V2_PASS_PRE_QM_GATE"
        if gate_pass
        else
        "QM_F06_UPPER_BOUNDARY_V2_FAIL_PRE_QM_GATE"
    )

    write_csv(
        OUTPUT_DIR
        / "QM_F06_UPPER_BOUNDARY_V2_valence_audit.csv",
        valence_rows,
    )

    write_csv(
        OUTPUT_DIR
        / "QM_F06_UPPER_BOUNDARY_V2_bond_audit.csv",
        bond_rows,
    )

    write_csv(
        OUTPUT_DIR
        / "QM_F06_UPPER_BOUNDARY_V2_nonbonded_audit.csv",
        nonbonded_rows,
    )

    summary = {
        "decision": decision,
        "formula": formula,
        "atom_count": len(atom_rows),
        "element_counts": dict(element_counts),
        "neutral_valence_electrons": neutral_valence_electrons,
        "provisional_charge": provisional_charge,
        "provisional_multiplicity": provisional_multiplicity,
        "connectivity_edges": len(edge_records),
        "graph_connected": graph_connected,
        "valence_failures": valence_failures,
        "bond_failures": bond_failures,
        "seam_bond_failures": seam_bond_failures,
        "cap_bond_failures": cap_bond_failures,
        "hard_contacts": hard_contacts,
        "hard_cap_contacts": hard_cap_contacts,
        "bridge_cap_hard_contacts": bridge_cap_hard_contacts,
        "unintended_covalent_contacts": (
            unintended_covalent_contacts
        ),
        "pre_qm_gate_pass": gate_pass,
        "orca_input_preparation_authorized": gate_pass,
        "orca_execution_authorized": False,
    }

    (
        OUTPUT_DIR
        / "QM_F06_UPPER_BOUNDARY_V2_pre_qm_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    report_path = (
        OUTPUT_DIR
        / "QM_F06_UPPER_BOUNDARY_V2_PRE_QM_AUDIT.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 UPPER Boundary V2 Pre-QM Audit — Day028",
                "",
                f"## Decision: **{decision}**",
                "",
                "## Composition",
                "",
                f"- Formula: **{formula}**",
                f"- Atoms: **{len(atom_rows)}**",
                (
                    "- Neutral valence-electron count: "
                    f"**{neutral_valence_electrons}**"
                ),
                (
                    "- Provisional charge/multiplicity: "
                    f"**{provisional_charge} / "
                    f"{provisional_multiplicity}**"
                ),
                "",
                "## Connectivity",
                "",
                f"- Reconstructed edges: **{len(edge_records)}**",
                (
                    "- Graph connected: "
                    f"**{'YES' if graph_connected else 'NO'}**"
                ),
                f"- Valence failures: **{valence_failures}**",
                "",
                "## Bond geometry",
                "",
                f"- Bond-range failures: **{bond_failures}**",
                (
                    "- Coordinate-seam bond failures: "
                    f"**{seam_bond_failures}**"
                ),
                (
                    "- Artificial-cap bond failures: "
                    f"**{cap_bond_failures}**"
                ),
                "",
                "## Topology-aware contacts",
                "",
                f"- Hard contacts: **{hard_contacts}**",
                (
                    "- Hard contacts involving caps: "
                    f"**{hard_cap_contacts}**"
                ),
                (
                    "- Hard bridge–cap contacts: "
                    f"**{bridge_cap_hard_contacts}**"
                ),
                (
                    "- Possible unintended covalent contacts: "
                    f"**{unintended_covalent_contacts}**"
                ),
                "",
                "## Authorization state",
                "",
                (
                    "- UPPER V2 pre-QM gate: "
                    f"**{'PASS' if gate_pass else 'FAIL'}**"
                ),
                (
                    "- ORCA input preparation: "
                    f"**{'AUTHORIZED' if gate_pass else 'NOT AUTHORIZED'}**"
                ),
                "- ORCA execution: **NOT AUTHORIZED**",
                "- ESP/RESP execution: **NOT AUTHORIZED**",
                "- Force-field parameter adoption: **NOT AUTHORIZED**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print("QM_F06 UPPER Boundary V2 pre-QM audit completed.")
    print("Decision:", decision)
    print("Formula:", formula)
    print("Valence electrons:", neutral_valence_electrons)
    print("Edges:", len(edge_records))
    print("Valence failures:", valence_failures)
    print("Bond failures:", bond_failures)
    print("Seam bond failures:", seam_bond_failures)
    print("Hard cap contacts:", hard_cap_contacts)
    print(
        "Bridge-cap hard contacts:",
        bridge_cap_hard_contacts,
    )
    print(
        "Unintended covalent contacts:",
        unintended_covalent_contacts,
    )
    print(
        "ORCA input preparation authorized:",
        gate_pass,
    )
    print("QM execution authorized: False")
    print("Report:", report_path)


if __name__ == "__main__":
    main()
