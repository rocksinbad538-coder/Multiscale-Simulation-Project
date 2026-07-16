#!/usr/bin/env python3
"""
Connectivity-aware structural audit of the converged QM_F06 LOWER Stage 1.

Compares the repaired initial fragment against the converged ORCA geometry.

Checks:
- original bonded distances before/after optimization;
- B-N, B-H and N-H bond-range gates;
- bond stretching and compression;
- bridge-specific bonds;
- bridge and attachment-center angles;
- bridge torsions;
- unexpected short nonbonded contacts;
- preservation of fragment connectivity.

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

F06_DIR = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/QM_F06"
)

WORKFLOW_DIR = F06_DIR / (
    "orca_workflow/QM_F06_LOWER_CAPPED_REPAIRED"
)

INITIAL_ATOMS = F06_DIR / (
    "QM_F06_LOWER_CAPPED_REPAIRED_atoms.csv"
)

ORIGINAL_FRAGMENT_ATOMS = F06_DIR / (
    "QM_F06_LOWER_atoms.csv"
)

INTERNAL_EDGES = F06_DIR / (
    "QM_F06_LOWER_internal_edges.csv"
)

BOUNDARY_AUDIT = F06_DIR / (
    "QM_F06_LOWER_boundary_edge_audit.csv"
)

CAPS = F06_DIR / (
    "QM_F06_LOWER_CAPPED_caps.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day027_qm_f06_lower_stage1_validation/"
    "local_structure_audit"
)

BOND_RANGE = {
    frozenset(("B", "N")): (1.20, 1.85),
    frozenset(("B", "H")): (0.90, 1.45),
    frozenset(("N", "H")): (0.80, 1.30),
}

VDW_RADII = {
    "H": 1.20,
    "B": 1.92,
    "N": 1.55,
}

HARD_NONBONDED_RATIO = 0.70
WARNING_NONBONDED_RATIO = 0.80


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(f"No rows found in {path}")

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
            {
                field: row.get(field, "")
                for field in fields
            }
            for row in rows
        )


def xyz_from_row(
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


def vector(
    origin: tuple[float, float, float],
    target: tuple[float, float, float],
) -> tuple[float, float, float]:
    return tuple(
        b - a
        for a, b in zip(origin, target, strict=True)
    )


def norm(v: tuple[float, float, float]) -> float:
    return math.sqrt(sum(component * component for component in v))


def dot(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> float:
    return sum(
        a * b
        for a, b in zip(first, second, strict=True)
    )


def cross(
    first: tuple[float, float, float],
    second: tuple[float, float, float],
) -> tuple[float, float, float]:
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def angle_deg(
    atom_1: tuple[float, float, float],
    center: tuple[float, float, float],
    atom_3: tuple[float, float, float],
) -> float:
    first = vector(center, atom_1)
    second = vector(center, atom_3)

    denominator = norm(first) * norm(second)
    if denominator <= 1.0e-14:
        raise RuntimeError("Undefined angle due to zero-length vector.")

    cosine = max(
        -1.0,
        min(1.0, dot(first, second) / denominator),
    )

    return math.degrees(math.acos(cosine))


def torsion_deg(
    atom_1: tuple[float, float, float],
    atom_2: tuple[float, float, float],
    atom_3: tuple[float, float, float],
    atom_4: tuple[float, float, float],
) -> float:
    b1 = vector(atom_2, atom_1)
    b2 = vector(atom_2, atom_3)
    b3 = vector(atom_3, atom_4)

    n1 = cross(b1, b2)
    n2 = cross(b2, b3)

    if norm(n1) <= 1.0e-14 or norm(n2) <= 1.0e-14:
        raise RuntimeError("Undefined torsion due to collinearity.")

    n1u = tuple(value / norm(n1) for value in n1)
    n2u = tuple(value / norm(n2) for value in n2)
    b2u = tuple(value / norm(b2) for value in b2)

    m1 = cross(n1u, b2u)

    return math.degrees(
        math.atan2(
            dot(m1, n2u),
            dot(n1u, n2u),
        )
    )


def angular_difference(first: float, second: float) -> float:
    return abs((second - first + 180.0) % 360.0 - 180.0)


def shortest_path_length(
    adjacency: dict[str, set[str]],
    source: str,
    target: str,
) -> int | None:
    if source == target:
        return 0

    queue = deque([(source, 0)])
    visited = {source}

    while queue:
        node, depth = queue.popleft()

        for neighbor in adjacency[node]:
            if neighbor == target:
                return depth + 1

            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, depth + 1))

    return None


def locate_optimized_xyz() -> Path:
    state_path = WORKFLOW_DIR / "workflow_state.json"
    require_file(state_path)

    state = json.loads(
        state_path.read_text(encoding="utf-8")
    )

    path_text = state.get("stage1_optimized_xyz")

    if not path_text:
        raise RuntimeError(
            "stage1_optimized_xyz is absent from workflow state."
        )

    path = ROOT / path_text
    require_file(path)
    return path


def read_xyz(path: Path) -> list[tuple[str, float, float, float]]:
    lines = path.read_text(encoding="utf-8").splitlines()

    count = int(lines[0].strip())
    coordinates = [
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

    if len(coordinates) != count:
        raise RuntimeError(
            f"XYZ count mismatch: header={count}, rows={len(coordinates)}"
        )

    return coordinates


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    initial_rows = read_csv(INITIAL_ATOMS)
    original_rows = read_csv(ORIGINAL_FRAGMENT_ATOMS)
    internal_edges = read_csv(INTERNAL_EDGES)
    boundary_rows = read_csv(BOUNDARY_AUDIT)
    cap_rows = read_csv(CAPS)

    optimized_xyz = locate_optimized_xyz()
    optimized_rows = read_xyz(optimized_xyz)

    if len(initial_rows) != len(optimized_rows):
        raise RuntimeError(
            "Initial and optimized atom counts differ."
        )

    atom_ids = [row["atom_id"] for row in initial_rows]

    initial_coordinates = {
        row["atom_id"]: xyz_from_row(row)
        for row in initial_rows
    }

    optimized_coordinates: dict[str, tuple[float, float, float]] = {}
    elements: dict[str, str] = {}

    for atom_id, initial, optimized in zip(
        atom_ids,
        initial_rows,
        optimized_rows,
        strict=True,
    ):
        element, x, y, z = optimized

        if element != initial["element"]:
            raise RuntimeError(
                f"Element mismatch for {atom_id}: "
                f"{initial['element']} vs {element}"
            )

        optimized_coordinates[atom_id] = (x, y, z)
        elements[atom_id] = element

    original_flags = {
        row["node_id"]: {
            "is_bridge_atom": (
                row["is_bridge_atom"].lower() == "true"
            ),
            "is_attachment_center": (
                row["is_attachment_center"].lower() == "true"
            ),
        }
        for row in original_rows
    }

    adjacency: dict[str, set[str]] = defaultdict(set)
    edge_origin: dict[tuple[str, str], str] = {}

    def add_edge(
        first: str,
        second: str,
        origin: str,
    ) -> None:
        edge = tuple(sorted((first, second)))
        adjacency[first].add(second)
        adjacency[second].add(first)
        edge_origin[edge] = origin

    for row in internal_edges:
        add_edge(
            row["source_node"],
            row["target_node"],
            "ORIGINAL_INTERNAL_EDGE",
        )

    for row in boundary_rows:
        if row["preliminary_action"] == "INCLUDE_EXISTING_HYDROGEN":
            add_edge(
                row["inside_node"],
                row["outside_node"],
                "RESTORED_EXISTING_HYDROGEN_EDGE",
            )

    for row in cap_rows:
        add_edge(
            row["parent_inside_node"],
            row["cap_id"],
            "ARTIFICIAL_CAP_EDGE",
        )

    bond_rows: list[dict[str, Any]] = []
    bond_failures: list[dict[str, Any]] = []

    for first, second in sorted(edge_origin):
        pair = frozenset((elements[first], elements[second]))
        initial_distance = distance(
            initial_coordinates[first],
            initial_coordinates[second],
        )
        optimized_distance = distance(
            optimized_coordinates[first],
            optimized_coordinates[second],
        )

        limits = BOND_RANGE.get(pair)
        bond_gate = (
            limits is not None
            and limits[0] <= optimized_distance <= limits[1]
        )

        relative_change = (
            (optimized_distance - initial_distance)
            / initial_distance
        )

        touches_bridge = (
            original_flags.get(first, {}).get("is_bridge_atom", False)
            or original_flags.get(second, {}).get("is_bridge_atom", False)
        )

        touches_attachment = (
            original_flags.get(first, {}).get(
                "is_attachment_center",
                False,
            )
            or original_flags.get(second, {}).get(
                "is_attachment_center",
                False,
            )
        )

        row = {
            "atom_1": first,
            "element_1": elements[first],
            "atom_2": second,
            "element_2": elements[second],
            "edge_origin": edge_origin[(first, second)],
            "initial_distance_angstrom": f"{initial_distance:.10f}",
            "optimized_distance_angstrom": f"{optimized_distance:.10f}",
            "absolute_change_angstrom": (
                f"{optimized_distance - initial_distance:.10f}"
            ),
            "relative_change_percent": f"{100.0 * relative_change:.6f}",
            "allowed_min_angstrom": limits[0] if limits else "",
            "allowed_max_angstrom": limits[1] if limits else "",
            "touches_bridge": touches_bridge,
            "touches_attachment": touches_attachment,
            "bond_range_gate_pass": bond_gate,
        }

        bond_rows.append(row)

        if not bond_gate:
            bond_failures.append(row)

    angle_rows: list[dict[str, Any]] = []

    for center in sorted(adjacency):
        neighbors = sorted(adjacency[center])

        for first_index, first in enumerate(neighbors):
            for third in neighbors[first_index + 1:]:
                initial_angle = angle_deg(
                    initial_coordinates[first],
                    initial_coordinates[center],
                    initial_coordinates[third],
                )
                optimized_angle = angle_deg(
                    optimized_coordinates[first],
                    optimized_coordinates[center],
                    optimized_coordinates[third],
                )

                angle_rows.append(
                    {
                        "atom_1": first,
                        "center_atom": center,
                        "atom_3": third,
                        "element_pattern": (
                            f"{elements[first]}-"
                            f"{elements[center]}-"
                            f"{elements[third]}"
                        ),
                        "initial_angle_deg": f"{initial_angle:.8f}",
                        "optimized_angle_deg": f"{optimized_angle:.8f}",
                        "change_deg": (
                            f"{optimized_angle - initial_angle:.8f}"
                        ),
                        "center_is_bridge": (
                            original_flags.get(center, {}).get(
                                "is_bridge_atom",
                                False,
                            )
                        ),
                        "center_is_attachment": (
                            original_flags.get(center, {}).get(
                                "is_attachment_center",
                                False,
                            )
                        ),
                    }
                )

    torsion_rows: list[dict[str, Any]] = []
    seen_torsions: set[tuple[str, str, str, str]] = set()

    for second in adjacency:
        for third in adjacency[second]:
            if second >= third:
                continue

            for first in adjacency[second] - {third}:
                for fourth in adjacency[third] - {second}:
                    torsion = (first, second, third, fourth)
                    reverse = tuple(reversed(torsion))
                    key = min(torsion, reverse)

                    if key in seen_torsions:
                        continue

                    seen_torsions.add(key)

                    initial_torsion = torsion_deg(
                        initial_coordinates[first],
                        initial_coordinates[second],
                        initial_coordinates[third],
                        initial_coordinates[fourth],
                    )
                    optimized_torsion = torsion_deg(
                        optimized_coordinates[first],
                        optimized_coordinates[second],
                        optimized_coordinates[third],
                        optimized_coordinates[fourth],
                    )

                    involves_bridge = any(
                        original_flags.get(atom, {}).get(
                            "is_bridge_atom",
                            False,
                        )
                        for atom in torsion
                    )

                    torsion_rows.append(
                        {
                            "atom_1": first,
                            "atom_2": second,
                            "atom_3": third,
                            "atom_4": fourth,
                            "element_pattern": "-".join(
                                elements[atom]
                                for atom in torsion
                            ),
                            "initial_torsion_deg": (
                                f"{initial_torsion:.8f}"
                            ),
                            "optimized_torsion_deg": (
                                f"{optimized_torsion:.8f}"
                            ),
                            "absolute_change_deg": (
                                f"{angular_difference(initial_torsion, optimized_torsion):.8f}"
                            ),
                            "involves_bridge": involves_bridge,
                        }
                    )

    nonbonded_rows: list[dict[str, Any]] = []
    blocking_nonbonded = 0

    for first_index, first in enumerate(atom_ids):
        for second in atom_ids[first_index + 1:]:
            separation = shortest_path_length(
                adjacency,
                first,
                second,
            )

            if separation in {1, 2, 3}:
                continue

            measured = distance(
                optimized_coordinates[first],
                optimized_coordinates[second],
            )

            vdw_sum = (
                VDW_RADII[elements[first]]
                + VDW_RADII[elements[second]]
            )
            ratio = measured / vdw_sum

            if ratio >= 0.90:
                continue

            involves_cap = (
                initial_rows[atom_ids.index(first)]["artificial_cap"].lower()
                == "true"
                or initial_rows[atom_ids.index(second)]["artificial_cap"].lower()
                == "true"
            )

            hard = ratio < HARD_NONBONDED_RATIO
            warning = (
                HARD_NONBONDED_RATIO <= ratio
                < WARNING_NONBONDED_RATIO
            )

            if hard:
                blocking_nonbonded += 1

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
                    "optimized_distance_angstrom": f"{measured:.10f}",
                    "vdw_sum_angstrom": f"{vdw_sum:.10f}",
                    "distance_over_vdw_sum": f"{ratio:.10f}",
                    "involves_artificial_cap": involves_cap,
                    "hard_clash": hard,
                    "strong_compression": warning,
                }
            )

    bridge_bonds = [
        row for row in bond_rows
        if row["touches_bridge"]
    ]
    attachment_bonds = [
        row for row in bond_rows
        if row["touches_attachment"]
    ]
    bridge_angles = [
        row for row in angle_rows
        if row["center_is_bridge"]
    ]
    bridge_torsions = [
        row for row in torsion_rows
        if row["involves_bridge"]
    ]

    maximum_bond_change = max(
        bond_rows,
        key=lambda row: abs(
            float(row["absolute_change_angstrom"])
        ),
    )
    maximum_angle_change = max(
        angle_rows,
        key=lambda row: abs(float(row["change_deg"])),
    )
    maximum_bridge_torsion_change = max(
        bridge_torsions,
        key=lambda row: float(row["absolute_change_deg"]),
    )

    connectivity_pass = len(bond_failures) == 0
    nonbonded_pass = blocking_nonbonded == 0
    stage2_authorized = connectivity_pass and nonbonded_pass

    write_csv(
        OUTPUT_DIR / "QM_F06_LOWER_STAGE1_bond_audit.csv",
        bond_rows,
    )
    write_csv(
        OUTPUT_DIR / "QM_F06_LOWER_STAGE1_angle_audit.csv",
        angle_rows,
    )
    write_csv(
        OUTPUT_DIR / "QM_F06_LOWER_STAGE1_torsion_audit.csv",
        torsion_rows,
    )
    write_csv(
        OUTPUT_DIR / "QM_F06_LOWER_STAGE1_nonbonded_audit.csv",
        nonbonded_rows,
    )

    decision = (
        "QM_F06_LOWER_STAGE1_LOCAL_STRUCTURE_PASS_STAGE2_AUTHORIZED"
        if stage2_authorized
        else
        "QM_F06_LOWER_STAGE1_LOCAL_STRUCTURE_FAIL_STAGE2_BLOCKED"
    )

    report = OUTPUT_DIR / (
        "QM_F06_LOWER_STAGE1_LOCAL_STRUCTURE_AUDIT.md"
    )

    report.write_text(
        "\n".join(
            [
                "# QM_F06 LOWER Stage-1 Local Structure Audit — Day027",
                "",
                f"## Decision: **{decision}**",
                "",
                "## Connectivity",
                "",
                f"- Reconstructed bonds: **{len(bond_rows)}**",
                f"- Bond-range failures: **{len(bond_failures)}**",
                f"- Bridge-touching bonds: **{len(bridge_bonds)}**",
                (
                    "- Attachment-touching bonds: "
                    f"**{len(attachment_bonds)}**"
                ),
                (
                    "- Maximum bond-length change: "
                    f"**{abs(float(maximum_bond_change['absolute_change_angstrom'])):.6f} Å** "
                    f"for `{maximum_bond_change['atom_1']} — "
                    f"{maximum_bond_change['atom_2']}`"
                ),
                "",
                "## Angular response",
                "",
                f"- Total angles: **{len(angle_rows)}**",
                f"- Bridge-centered angles: **{len(bridge_angles)}**",
                (
                    "- Maximum angular change: "
                    f"**{abs(float(maximum_angle_change['change_deg'])):.4f}°** "
                    f"for `{maximum_angle_change['atom_1']} — "
                    f"{maximum_angle_change['center_atom']} — "
                    f"{maximum_angle_change['atom_3']}`"
                ),
                "",
                "## Torsional response",
                "",
                f"- Total torsions: **{len(torsion_rows)}**",
                (
                    "- Bridge-involving torsions: "
                    f"**{len(bridge_torsions)}**"
                ),
                (
                    "- Maximum bridge torsional change: "
                    f"**{float(maximum_bridge_torsion_change['absolute_change_deg']):.4f}°** "
                    f"for `{maximum_bridge_torsion_change['atom_1']} — "
                    f"{maximum_bridge_torsion_change['atom_2']} — "
                    f"{maximum_bridge_torsion_change['atom_3']} — "
                    f"{maximum_bridge_torsion_change['atom_4']}`"
                ),
                "",
                "## Long-range contacts",
                "",
                (
                    "- Contacts below 0.90 of vdW sum: "
                    f"**{len(nonbonded_rows)}**"
                ),
                (
                    "- Hard clashes below 0.70 of vdW sum: "
                    f"**{blocking_nonbonded}**"
                ),
                "",
                "## Gate state",
                "",
                (
                    "- Connectivity/bond-range gate: "
                    f"**{'PASS' if connectivity_pass else 'FAIL'}**"
                ),
                (
                    "- Long-range steric gate: "
                    f"**{'PASS' if nonbonded_pass else 'FAIL'}**"
                ),
                (
                    "- Stage-2 execution: "
                    f"**{'AUTHORIZED' if stage2_authorized else 'NOT AUTHORIZED'}**"
                ),
                "",
                (
                    "Large Cartesian displacements are acceptable only if "
                    "connectivity, bond ranges and nonbonded topology remain "
                    "chemically valid. This report applies that criterion."
                ),
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = {
        "decision": decision,
        "optimized_xyz": str(optimized_xyz.relative_to(ROOT)),
        "bond_count": len(bond_rows),
        "bond_failures": len(bond_failures),
        "hard_nonbonded_clashes": blocking_nonbonded,
        "stage2_execution_authorized": stage2_authorized,
        "qm_execution_performed_by_this_script": False,
        "required_next_step": (
            "RUN_STAGE2_LOWER_PREFLIGHT"
            if stage2_authorized
            else "REVIEW_FAILED_LOCAL_STRUCTURE_GATES"
        ),
    }

    (
        OUTPUT_DIR
        / "QM_F06_LOWER_STAGE1_local_structure_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("QM_F06 LOWER Stage-1 local-structure audit completed.")
    print(f"Decision: {decision}")
    print("Bond failures:", len(bond_failures))
    print("Hard nonbonded clashes:", blocking_nonbonded)
    print("Stage-2 execution authorized:", stage2_authorized)
    print(f"Report: {report}")


if __name__ == "__main__":
    main()
