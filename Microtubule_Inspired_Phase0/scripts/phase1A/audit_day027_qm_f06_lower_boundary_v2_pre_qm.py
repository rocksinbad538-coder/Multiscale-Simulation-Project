#!/usr/bin/env python3
"""
Pre-QM audit of QM_F06 LOWER Boundary V2.

The V2 fragment combines:
- converged Stage-2 coordinates for retained V1 atoms;
- validated Day024 coordinates for restored real R2 atoms;
- three new artificial peripheral caps.

This audit reconstructs the intended topology and evaluates:

1. atom count, formula and electronic parity;
2. complete intended connectivity;
3. graph degree/valence;
4. all B-N, B-H and N-H bond distances;
5. the Stage2/Day024 coordinate seam;
6. topology-aware nonbonded contacts;
7. artificial-cap proximity to the bridge;
8. possible unintended covalent contacts.

No atom is moved and no QM calculation is executed.
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

ATOMS_PATH = V2_DIR / "QM_F06_LOWER_BOUNDARY_V2_atoms.csv"
CAPS_PATH = V2_DIR / "QM_F06_LOWER_BOUNDARY_V2_caps.csv"

FULL_EDGES_PATH = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_edges.csv"
)

OUTPUT_DIR = V2_DIR / "pre_qm_audit"

BOND_RANGES_ANGSTROM = {
    tuple(sorted(("B", "N"))): (1.20, 1.85),
    tuple(sorted(("B", "H"))): (0.90, 1.45),
    tuple(sorted(("N", "H"))): (0.80, 1.30),
}

EXPECTED_DEGREES = {
    "B": 3,
    "N": 3,
    "H": 1,
}

VDW_RADII = {
    "H": 1.20,
    "B": 1.92,
    "N": 1.55,
}

UNINTENDED_COVALENT_THRESHOLDS = {
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

RESTORED_REAL_ATOMS = {
    "A:LOWER:11:-3",
    "A:LOWER:13:-3",
    "A:LOWER:14:-2",
    "H4:LOWER:0016:0",
}

CRITICAL_RESTORED_EDGE = tuple(
    sorted(("A:LOWER:14:-4", "A:LOWER:13:-3"))
)


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty required file: {path}")


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(f"No rows found in: {path}")

    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows available for: {path}")

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


def canonical_edge(
    first: str,
    second: str,
) -> tuple[str, str]:
    return tuple(sorted((first, second)))


def coordinates(
    row: dict[str, str],
) -> tuple[float, float, float]:
    return (
        float(row["x_angstrom"]),
        float(row["y_angstrom"]),
        float(row["z_angstrom"]),
    )


def distance(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> float:
    return math.sqrt(
        sum(
            (a - b) ** 2
            for a, b in zip(first, second, strict=True)
        )
    )


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

    atoms = {
        row["atom_id"]: row
        for row in atom_rows
    }

    atom_ids = list(atoms)
    atom_set = set(atom_ids)

    if len(atom_rows) != 28:
        raise RuntimeError(
            f"Expected 28 atoms; found {len(atom_rows)}"
        )

    if len(atom_set) != len(atom_rows):
        raise RuntimeError("Duplicate atom identifiers detected.")

    elements = {
        atom_id: row["element"]
        for atom_id, row in atoms.items()
    }

    coords = {
        atom_id: coordinates(row)
        for atom_id, row in atoms.items()
    }

    artificial_caps = {
        atom_id
        for atom_id, row in atoms.items()
        if row["artificial_cap"].strip().lower() == "true"
    }

    # Reconstruct all real R2 edges whose endpoints are both present.
    edge_records: dict[tuple[str, str], dict[str, Any]] = {}

    for row in full_edge_rows:
        first = row["source_node"]
        second = row["target_node"]

        if first in atom_set and second in atom_set:
            key = canonical_edge(first, second)

            edge_records[key] = {
                "atom_1": key[0],
                "atom_2": key[1],
                "edge_origin": "REAL_R2_GRAPH_EDGE",
                "source_edge_id": row["edge_id"],
                "original_edge_type": row["edge_type"],
            }

    # Reconstruct every artificial cap edge from atom metadata.
    for atom_id in sorted(artificial_caps):
        parent = atoms[atom_id]["parent_inside_node"]

        if not parent:
            raise RuntimeError(
                f"Artificial cap has no parent: {atom_id}"
            )

        if parent not in atom_set:
            raise RuntimeError(
                f"Artificial cap parent absent: {atom_id} -> {parent}"
            )

        key = canonical_edge(atom_id, parent)

        edge_records[key] = {
            "atom_1": key[0],
            "atom_2": key[1],
            "edge_origin": "ARTIFICIAL_CAP_EDGE",
            "source_edge_id": atoms[atom_id]["source_edge_id"],
            "original_edge_type": "QM_BOUNDARY_CAP",
        }

    adjacency: dict[str, set[str]] = defaultdict(set)

    for first, second in edge_records:
        adjacency[first].add(second)
        adjacency[second].add(first)

    # Formula and provisional electronic state.
    element_counts = Counter(elements.values())

    valence_electrons = (
        3 * element_counts["B"]
        + 5 * element_counts["N"]
        + element_counts["H"]
    )

    provisional_charge = 0
    provisional_multiplicity = (
        1 if valence_electrons % 2 == 0 else 2
    )

    # Degree and valence audit.
    valence_rows: list[dict[str, Any]] = []
    valence_failures = 0

    for atom_id in atom_ids:
        degree = len(adjacency[atom_id])
        expected = EXPECTED_DEGREES[elements[atom_id]]
        passed = degree == expected

        if not passed:
            valence_failures += 1

        valence_rows.append(
            {
                "atom_id": atom_id,
                "element": elements[atom_id],
                "atom_role": atoms[atom_id]["atom_role"],
                "node_type": atoms[atom_id]["node_type"],
                "graph_degree": degree,
                "expected_degree": expected,
                "neighbors": "|".join(
                    sorted(adjacency[atom_id])
                ),
                "degree_gate_pass": passed,
            }
        )

    # Bond-distance and seam audit.
    bond_rows: list[dict[str, Any]] = []
    bond_failures = 0
    seam_bond_failures = 0
    cap_bond_failures = 0

    for edge in sorted(edge_records):
        first, second = edge
        record = edge_records[edge]

        pair = tuple(sorted((elements[first], elements[second])))
        limits = BOND_RANGES_ANGSTROM.get(pair)

        measured = distance(
            coords[first],
            coords[second],
        )

        passed = (
            limits is not None
            and limits[0] <= measured <= limits[1]
        )

        if not passed:
            bond_failures += 1

        is_cap_edge = (
            record["edge_origin"] == "ARTIFICIAL_CAP_EDGE"
        )

        if is_cap_edge and not passed:
            cap_bond_failures += 1

        coordinate_source_1 = atoms[first]["coordinate_source"]
        coordinate_source_2 = atoms[second]["coordinate_source"]

        crosses_coordinate_seam = (
            coordinate_source_1 != coordinate_source_2
            and not is_cap_edge
        )

        if crosses_coordinate_seam and not passed:
            seam_bond_failures += 1

        touches_restored_region = (
            first in RESTORED_REAL_ATOMS
            or second in RESTORED_REAL_ATOMS
        )

        bond_rows.append(
            {
                "atom_1": first,
                "element_1": elements[first],
                "atom_2": second,
                "element_2": elements[second],
                "edge_origin": record["edge_origin"],
                "source_edge_id": record["source_edge_id"],
                "original_edge_type": record["original_edge_type"],
                "coordinate_source_1": coordinate_source_1,
                "coordinate_source_2": coordinate_source_2,
                "crosses_coordinate_seam": crosses_coordinate_seam,
                "touches_restored_region": touches_restored_region,
                "critical_restored_edge": edge == CRITICAL_RESTORED_EDGE,
                "distance_angstrom": f"{measured:.10f}",
                "allowed_min_angstrom": (
                    limits[0] if limits else ""
                ),
                "allowed_max_angstrom": (
                    limits[1] if limits else ""
                ),
                "bond_distance_gate_pass": passed,
            }
        )

    # Topology-aware nonbonded audit.
    nonbonded_rows: list[dict[str, Any]] = []
    hard_contacts = 0
    hard_cap_contacts = 0
    unintended_covalent_contacts = 0
    bridge_cap_hard_contacts = 0

    for index, first in enumerate(atom_ids):
        for second in atom_ids[index + 1:]:
            graph_separation = shortest_path_length(
                adjacency,
                first,
                second,
            )

            # Exclude bonded, angular and torsional pairs.
            if graph_separation in {1, 2, 3}:
                continue

            measured = distance(
                coords[first],
                coords[second],
            )

            pair = tuple(sorted((elements[first], elements[second])))
            vdw_sum = (
                VDW_RADII[elements[first]]
                + VDW_RADII[elements[second]]
            )
            ratio = measured / vdw_sum

            covalent_threshold = (
                UNINTENDED_COVALENT_THRESHOLDS[pair]
            )
            possible_new_bond = measured <= covalent_threshold

            involves_cap = (
                first in artificial_caps
                or second in artificial_caps
            )
            involves_bridge = (
                first in BRIDGE_ATOMS
                or second in BRIDGE_ATOMS
            )
            hard = ratio < 0.70
            warning = 0.70 <= ratio < 0.80

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

            nonbonded_rows.append(
                {
                    "atom_1": first,
                    "element_1": elements[first],
                    "atom_role_1": atoms[first]["atom_role"],
                    "atom_2": second,
                    "element_2": elements[second],
                    "atom_role_2": atoms[second]["atom_role"],
                    "graph_separation": (
                        graph_separation
                        if graph_separation is not None
                        else "DISCONNECTED"
                    ),
                    "distance_angstrom": f"{measured:.10f}",
                    "vdw_sum_angstrom": f"{vdw_sum:.10f}",
                    "distance_over_vdw_sum": f"{ratio:.10f}",
                    "covalent_threshold_angstrom": (
                        f"{covalent_threshold:.10f}"
                    ),
                    "possible_unintended_covalent_contact": (
                        possible_new_bond
                    ),
                    "involves_artificial_cap": involves_cap,
                    "involves_bridge": involves_bridge,
                    "hard_contact_below_0p70": hard,
                    "strong_compression_0p70_to_0p80": warning,
                }
            )

    restored_bonds = [
        row
        for row in bond_rows
        if row["touches_restored_region"]
    ]

    seam_bonds = [
        row
        for row in bond_rows
        if row["crosses_coordinate_seam"]
    ]

    critical_edge_rows = [
        row
        for row in bond_rows
        if row["critical_restored_edge"]
    ]

    if len(critical_edge_rows) != 1:
        raise RuntimeError(
            "Critical restored edge was not reconstructed exactly once."
        )

    critical_edge = critical_edge_rows[0]

    formula = (
        f"B{element_counts['B']}"
        f"N{element_counts['N']}"
        f"H{element_counts['H']}"
    )

    graph_connected = False

    if atom_ids:
        visited = {atom_ids[0]}
        queue = deque([atom_ids[0]])

        while queue:
            node = queue.popleft()

            for neighbor in adjacency[node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)

        graph_connected = len(visited) == len(atom_ids)

    pre_qm_pass = all(
        (
            graph_connected,
            valence_failures == 0,
            bond_failures == 0,
            seam_bond_failures == 0,
            cap_bond_failures == 0,
            unintended_covalent_contacts == 0,
            hard_cap_contacts == 0,
            bridge_cap_hard_contacts == 0,
        )
    )

    decision = (
        "QM_F06_LOWER_BOUNDARY_V2_PASS_PRE_QM_GATE"
        if pre_qm_pass
        else "QM_F06_LOWER_BOUNDARY_V2_FAIL_PRE_QM_GATE"
    )

    write_csv(
        OUTPUT_DIR
        / "QM_F06_LOWER_BOUNDARY_V2_valence_audit.csv",
        valence_rows,
    )

    write_csv(
        OUTPUT_DIR
        / "QM_F06_LOWER_BOUNDARY_V2_bond_audit.csv",
        bond_rows,
    )

    if nonbonded_rows:
        write_csv(
            OUTPUT_DIR
            / "QM_F06_LOWER_BOUNDARY_V2_nonbonded_audit.csv",
            nonbonded_rows,
        )

    gate_rows = [
        {
            "gate": "ATOM_COUNT",
            "required": 28,
            "observed": len(atom_rows),
            "pass": len(atom_rows) == 28,
        },
        {
            "gate": "GRAPH_CONNECTED",
            "required": True,
            "observed": graph_connected,
            "pass": graph_connected,
        },
        {
            "gate": "VALENCE_FAILURES",
            "required": 0,
            "observed": valence_failures,
            "pass": valence_failures == 0,
        },
        {
            "gate": "BOND_DISTANCE_FAILURES",
            "required": 0,
            "observed": bond_failures,
            "pass": bond_failures == 0,
        },
        {
            "gate": "COORDINATE_SEAM_BOND_FAILURES",
            "required": 0,
            "observed": seam_bond_failures,
            "pass": seam_bond_failures == 0,
        },
        {
            "gate": "ARTIFICIAL_CAP_BOND_FAILURES",
            "required": 0,
            "observed": cap_bond_failures,
            "pass": cap_bond_failures == 0,
        },
        {
            "gate": "UNINTENDED_COVALENT_CONTACTS",
            "required": 0,
            "observed": unintended_covalent_contacts,
            "pass": unintended_covalent_contacts == 0,
        },
        {
            "gate": "HARD_CAP_CONTACTS",
            "required": 0,
            "observed": hard_cap_contacts,
            "pass": hard_cap_contacts == 0,
        },
        {
            "gate": "BRIDGE_CAP_HARD_CONTACTS",
            "required": 0,
            "observed": bridge_cap_hard_contacts,
            "pass": bridge_cap_hard_contacts == 0,
        },
    ]

    write_csv(
        OUTPUT_DIR
        / "QM_F06_LOWER_BOUNDARY_V2_pre_qm_gates.csv",
        gate_rows,
    )

    report_path = (
        OUTPUT_DIR
        / "QM_F06_LOWER_BOUNDARY_V2_PRE_QM_AUDIT.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 LOWER Boundary V2 Pre-QM Audit — Day027",
                "",
                f"## Decision: **{decision}**",
                "",
                "## Composition and electronic parity",
                "",
                f"- Formula: **{formula}**",
                f"- Atoms: **{len(atom_rows)}**",
                f"- B atoms: **{element_counts['B']}**",
                f"- N atoms: **{element_counts['N']}**",
                f"- H atoms: **{element_counts['H']}**",
                (
                    "- Neutral valence-electron count: "
                    f"**{valence_electrons}**"
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
                    "- Graph fully connected: "
                    f"**{'YES' if graph_connected else 'NO'}**"
                ),
                f"- Degree/valence failures: **{valence_failures}**",
                "",
                "## Bond geometry",
                "",
                f"- Bond-distance failures: **{bond_failures}**",
                (
                    "- Coordinate-seam bond failures: "
                    f"**{seam_bond_failures}**"
                ),
                (
                    "- Artificial-cap bond failures: "
                    f"**{cap_bond_failures}**"
                ),
                (
                    "- Bonds touching restored region: "
                    f"**{len(restored_bonds)}**"
                ),
                (
                    "- Stage2/Day024 seam bonds: "
                    f"**{len(seam_bonds)}**"
                ),
                (
                    "- Critical restored edge "
                    "`A:LOWER:14:-4 — A:LOWER:13:-3`: "
                    f"**{float(critical_edge['distance_angstrom']):.6f} Å**"
                ),
                "",
                "## Topology-aware nonbonded audit",
                "",
                (
                    "- Contacts reported below 0.90 vdW ratio "
                    "or covalent threshold: "
                    f"**{len(nonbonded_rows)}**"
                ),
                f"- Hard contacts below 0.70: **{hard_contacts}**",
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
                    "- V2 pre-QM geometry gate: "
                    f"**{'PASS' if pre_qm_pass else 'FAIL'}**"
                ),
                "- ORCA input preparation: **NOT YET AUTHORIZED**",
                "- QM execution: **NOT AUTHORIZED**",
                "",
                "## Interpretation boundary",
                "",
                (
                    "The provisional neutral singlet assignment follows "
                    "only from stoichiometric electron parity. SCF "
                    "stability and the absence of a lower-energy "
                    "open-shell state remain to be tested after the "
                    "geometric gate passes."
                ),
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = {
        "decision": decision,
        "formula": formula,
        "atom_count": len(atom_rows),
        "valence_electron_count": valence_electrons,
        "provisional_charge": provisional_charge,
        "provisional_multiplicity": provisional_multiplicity,
        "edge_count": len(edge_records),
        "graph_connected": graph_connected,
        "valence_failures": valence_failures,
        "bond_distance_failures": bond_failures,
        "coordinate_seam_bond_failures": seam_bond_failures,
        "artificial_cap_bond_failures": cap_bond_failures,
        "hard_contacts": hard_contacts,
        "hard_cap_contacts": hard_cap_contacts,
        "bridge_cap_hard_contacts": bridge_cap_hard_contacts,
        "unintended_covalent_contacts": (
            unintended_covalent_contacts
        ),
        "pre_qm_gate_pass": pre_qm_pass,
        "orca_input_preparation_authorized": False,
        "qm_execution_authorized": False,
        "required_next_step": (
            "PREPARE_BOUNDARY_V2_CONSTRAINED_OPTIMIZATION"
            if pre_qm_pass
            else "REVIEW_BOUNDARY_V2_FAILED_GATES"
        ),
    }

    (
        OUTPUT_DIR
        / "QM_F06_LOWER_BOUNDARY_V2_pre_qm_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("QM_F06 LOWER Boundary V2 pre-QM audit completed.")
    print("Decision:", decision)
    print("Formula:", formula)
    print("Valence electrons:", valence_electrons)
    print("Edges:", len(edge_records))
    print("Valence failures:", valence_failures)
    print("Bond failures:", bond_failures)
    print("Seam bond failures:", seam_bond_failures)
    print("Hard cap contacts:", hard_cap_contacts)
    print("Bridge-cap hard contacts:", bridge_cap_hard_contacts)
    print(
        "Unintended covalent contacts:",
        unintended_covalent_contacts,
    )
    print("QM execution authorized: False")
    print("Report:", report_path)


if __name__ == "__main__":
    main()
