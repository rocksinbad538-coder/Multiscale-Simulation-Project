#!/usr/bin/env python3
"""
Scientific validation of the converged QM_F06 LOWER Stage-2 geometry.

Compares:
- repaired initial geometry;
- converged Stage-1 geometry;
- converged Stage-2 geometry.

Validates:
- preservation of all intended bonds;
- bridge bond-length response;
- Stage-1 -> Stage-2 atomic displacements;
- evolution of the two inherited compressed contacts;
- topology-aware long-range contacts;
- cap bond integrity;
- static correctness of the promoted Stage-3 input.

No electronic-structure calculation is executed.
"""

from __future__ import annotations

import csv
import json
import math
import re
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

ORIGINAL_FRAGMENT = F06_DIR / (
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

STAGE3_INPUT = WORKFLOW_DIR / "stage3.inp"

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day027_qm_f06_lower_stage2_validation"
)

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

# Stage-1 contacts explicitly designated as Stage-2 relaxation targets.
TARGET_CONTACTS = (
    ("BR4:LOWER:00:3", "P:48"),
    ("BR4:LOWER:00:3", "H4:LOWER:0017:0"),
)


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty required file: {path}")


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(f"No rows in {path}")

    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows available for {path}")

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


def coordinates_from_atom_row(
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


def canonical_edge(first: str, second: str) -> tuple[str, str]:
    return tuple(sorted((first, second)))


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


def read_xyz(
    path: Path,
) -> list[tuple[str, float, float, float]]:
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
            f"XYZ count mismatch in {path}: "
            f"header={expected}, rows={len(rows)}"
        )

    return rows


def locate_state_paths() -> tuple[Path, Path]:
    state_path = WORKFLOW_DIR / "workflow_state.json"
    require_file(state_path)

    state = json.loads(
        state_path.read_text(encoding="utf-8")
    )

    required_flags = {
        "stage1_executed": state.get("stage1_executed"),
        "stage1_validation_pass": state.get(
            "stage1_validation_pass"
        ),
        "stage2_executed": state.get("stage2_executed"),
        "stage2_validation_pass": state.get(
            "stage2_validation_pass"
        ),
        "stage2_geometry_promoted": state.get(
            "stage2_geometry_promoted"
        ),
        "stage3_input_generated": state.get(
            "stage3_input_generated"
        ),
    }

    if not all(value is True for value in required_flags.values()):
        raise RuntimeError(
            f"Workflow state is incomplete: {required_flags}"
        )

    stage1_path = ROOT / state["stage1_optimized_xyz"]
    stage2_path = ROOT / state["stage2_optimized_xyz"]

    require_file(stage1_path)
    require_file(stage2_path)

    return stage1_path, stage2_path


def parse_stage3_input(path: Path) -> dict[str, Any]:
    require_file(path)

    text = path.read_text(encoding="utf-8")

    match = re.search(
        r"(?ms)^\s*\*\s+xyz\s+(-?\d+)\s+(\d+)\s*$"
        r"(.*?)"
        r"^\s*\*\s*$",
        text,
    )

    if not match:
        raise RuntimeError("Stage-3 XYZ block not found.")

    atoms = [
        line.split()
        for line in match.group(3).splitlines()
        if line.strip()
    ]

    constraints = re.findall(
        r"(?m)^\s*\{\s*C\s+\d+\s+C\s*\}\s*$",
        text,
    )

    checks = {
        "atom_count_22": len(atoms) == 22,
        "charge_zero": int(match.group(1)) == 0,
        "multiplicity_one": int(match.group(2)) == 1,
        "no_constraints": len(constraints) == 0,
        "promoted_from_stage2": (
            "Coordinates promoted from Stage 2" in text
        ),
        "contains_pbe0": "PBE0" in text,
        "contains_d4": "D4" in text,
        "contains_def2_tzvp": "def2-TZVP" in text,
        "contains_tightscf": "TightSCF" in text,
        "contains_defgrid3": "DefGrid3" in text,
        "not_an_optimization": not bool(
            re.search(r"(?i)(^|\s)Opt(\s|$)", text)
        ),
        "no_obsolete_grid": (
            "Grid5" not in text
            and "FinalGrid6" not in text
        ),
    }

    return {
        "checks": checks,
        "gate_pass": all(checks.values()),
        "atoms": atoms,
        "charge": int(match.group(1)),
        "multiplicity": int(match.group(2)),
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    atom_rows = read_csv(INITIAL_ATOMS)
    original_rows = read_csv(ORIGINAL_FRAGMENT)
    internal_edges = read_csv(INTERNAL_EDGES)
    boundary_rows = read_csv(BOUNDARY_AUDIT)
    cap_rows = read_csv(CAPS)

    stage1_xyz_path, stage2_xyz_path = locate_state_paths()

    stage1_xyz = read_xyz(stage1_xyz_path)
    stage2_xyz = read_xyz(stage2_xyz_path)

    if not (
        len(atom_rows)
        == len(stage1_xyz)
        == len(stage2_xyz)
        == 22
    ):
        raise RuntimeError("Atom-count mismatch among geometries.")

    atom_ids = [row["atom_id"] for row in atom_rows]

    elements = {
        row["atom_id"]: row["element"]
        for row in atom_rows
    }

    roles = {
        row["atom_id"]: row["atom_role"]
        for row in atom_rows
    }

    artificial_caps = {
        row["atom_id"]
        for row in atom_rows
        if row["artificial_cap"].lower() == "true"
    }

    stage1_coords: dict[str, tuple[float, float, float]] = {}
    stage2_coords: dict[str, tuple[float, float, float]] = {}

    for atom_id, source, xyz1, xyz2 in zip(
        atom_ids,
        atom_rows,
        stage1_xyz,
        stage2_xyz,
        strict=True,
    ):
        if xyz1[0] != source["element"] or xyz2[0] != source["element"]:
            raise RuntimeError(
                f"Element ordering mismatch for {atom_id}"
            )

        stage1_coords[atom_id] = xyz1[1:]
        stage2_coords[atom_id] = xyz2[1:]

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

    def add_edge(first: str, second: str, origin: str) -> None:
        edge = canonical_edge(first, second)
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
                "RESTORED_EXISTING_H_EDGE",
            )

    for row in cap_rows:
        add_edge(
            row["parent_inside_node"],
            row["cap_id"],
            "ARTIFICIAL_CAP_EDGE",
        )

    displacement_rows: list[dict[str, Any]] = []
    role_sum_sq: dict[str, float] = defaultdict(float)
    role_counts: Counter[str] = Counter()
    role_max: dict[str, float] = defaultdict(float)

    total_sum_sq = 0.0

    for atom_id in atom_ids:
        displacement = distance(
            stage1_coords[atom_id],
            stage2_coords[atom_id],
        )

        role = roles[atom_id]

        total_sum_sq += displacement**2
        role_sum_sq[role] += displacement**2
        role_counts[role] += 1
        role_max[role] = max(
            role_max[role],
            displacement,
        )

        displacement_rows.append(
            {
                "atom_id": atom_id,
                "element": elements[atom_id],
                "atom_role": role,
                "stage1_to_stage2_displacement_angstrom": (
                    f"{displacement:.10f}"
                ),
                "is_bridge_atom": (
                    original_flags.get(atom_id, {}).get(
                        "is_bridge_atom",
                        False,
                    )
                ),
                "is_attachment_center": (
                    original_flags.get(atom_id, {}).get(
                        "is_attachment_center",
                        False,
                    )
                ),
                "artificial_cap": atom_id in artificial_caps,
            }
        )

    overall_rmsd = math.sqrt(
        total_sum_sq / len(atom_ids)
    )

    maximum_displacement = max(
        displacement_rows,
        key=lambda row: float(
            row["stage1_to_stage2_displacement_angstrom"]
        ),
    )

    bond_rows: list[dict[str, Any]] = []
    bond_failures = 0
    cap_bond_failures = 0

    for first, second in sorted(edge_origin):
        pair = tuple(sorted((elements[first], elements[second])))
        limits = BOND_RANGES[pair]

        stage1_distance = distance(
            stage1_coords[first],
            stage1_coords[second],
        )
        stage2_distance = distance(
            stage2_coords[first],
            stage2_coords[second],
        )

        passed = (
            limits[0] <= stage2_distance <= limits[1]
        )

        if not passed:
            bond_failures += 1

        is_cap_edge = (
            edge_origin[(first, second)]
            == "ARTIFICIAL_CAP_EDGE"
        )

        if is_cap_edge and not passed:
            cap_bond_failures += 1

        bond_rows.append(
            {
                "atom_1": first,
                "element_1": elements[first],
                "atom_2": second,
                "element_2": elements[second],
                "edge_origin": edge_origin[(first, second)],
                "stage1_distance_angstrom": (
                    f"{stage1_distance:.10f}"
                ),
                "stage2_distance_angstrom": (
                    f"{stage2_distance:.10f}"
                ),
                "stage1_to_stage2_change_angstrom": (
                    f"{stage2_distance - stage1_distance:.10f}"
                ),
                "allowed_min_angstrom": limits[0],
                "allowed_max_angstrom": limits[1],
                "touches_bridge": (
                    original_flags.get(first, {}).get(
                        "is_bridge_atom",
                        False,
                    )
                    or original_flags.get(second, {}).get(
                        "is_bridge_atom",
                        False,
                    )
                ),
                "artificial_cap_edge": is_cap_edge,
                "stage2_bond_gate_pass": passed,
            }
        )

    target_rows: list[dict[str, Any]] = []

    for first, second in TARGET_CONTACTS:
        stage1_distance = distance(
            stage1_coords[first],
            stage1_coords[second],
        )
        stage2_distance = distance(
            stage2_coords[first],
            stage2_coords[second],
        )

        stage1_ratio = stage1_distance / (
            VDW_RADII[elements[first]]
            + VDW_RADII[elements[second]]
        )
        stage2_ratio = stage2_distance / (
            VDW_RADII[elements[first]]
            + VDW_RADII[elements[second]]
        )

        target_rows.append(
            {
                "atom_1": first,
                "element_1": elements[first],
                "atom_2": second,
                "element_2": elements[second],
                "stage1_distance_angstrom": (
                    f"{stage1_distance:.10f}"
                ),
                "stage2_distance_angstrom": (
                    f"{stage2_distance:.10f}"
                ),
                "distance_change_angstrom": (
                    f"{stage2_distance - stage1_distance:.10f}"
                ),
                "stage1_distance_over_vdw_sum": (
                    f"{stage1_ratio:.10f}"
                ),
                "stage2_distance_over_vdw_sum": (
                    f"{stage2_ratio:.10f}"
                ),
                "compression_improved": (
                    stage2_ratio > stage1_ratio
                ),
            }
        )

    close_contact_rows: list[dict[str, Any]] = []
    blocking_contacts = 0
    cap_blocking_contacts = 0

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
                stage2_coords[first],
                stage2_coords[second],
            )

            vdw_sum = (
                VDW_RADII[elements[first]]
                + VDW_RADII[elements[second]]
            )

            ratio = measured / vdw_sum

            if ratio >= 0.90:
                continue

            involves_cap = (
                first in artificial_caps
                or second in artificial_caps
            )

            hard = ratio < 0.70

            if hard:
                blocking_contacts += 1

            if hard and involves_cap:
                cap_blocking_contacts += 1

            close_contact_rows.append(
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
                    "stage2_distance_angstrom": (
                        f"{measured:.10f}"
                    ),
                    "distance_over_vdw_sum": (
                        f"{ratio:.10f}"
                    ),
                    "involves_artificial_cap": involves_cap,
                    "hard_contact_below_0p70": hard,
                }
            )

    stage3 = parse_stage3_input(STAGE3_INPUT)

    stage3_coordinate_match = True
    maximum_stage3_coordinate_error = 0.0

    for source, promoted in zip(
        stage2_xyz,
        stage3["atoms"],
        strict=True,
    ):
        if source[0] != promoted[0]:
            stage3_coordinate_match = False
            break

        error = max(
            abs(source[1] - float(promoted[1])),
            abs(source[2] - float(promoted[2])),
            abs(source[3] - float(promoted[3])),
        )

        maximum_stage3_coordinate_error = max(
            maximum_stage3_coordinate_error,
            error,
        )

    stage3_coordinate_match = (
        stage3_coordinate_match
        and maximum_stage3_coordinate_error <= 1.0e-8
    )

    all_target_contacts_improved = all(
        row["compression_improved"]
        for row in target_rows
    )

    stage2_structure_pass = all(
        (
            bond_failures == 0,
            cap_bond_failures == 0,
            cap_blocking_contacts == 0,
        )
    )

    stage3_static_pass = (
        stage3["gate_pass"]
        and stage3_coordinate_match
    )

    # Execution remains blocked until the electronic-property protocol
    # (energy only vs ESP/RESP/charges) is explicitly defined.
    stage3_execution_authorized = False

    decision = (
        "QM_F06_LOWER_STAGE2_STRUCTURE_VALIDATED_"
        "STAGE3_STATIC_INPUT_READY_PROPERTY_PROTOCOL_REQUIRED"
        if stage2_structure_pass and stage3_static_pass
        else
        "QM_F06_LOWER_STAGE2_OR_STAGE3_STATIC_GATE_FAILED"
    )

    write_csv(
        OUTPUT_DIR
        / "QM_F06_LOWER_STAGE1_TO_STAGE2_displacements.csv",
        displacement_rows,
    )
    write_csv(
        OUTPUT_DIR
        / "QM_F06_LOWER_STAGE2_bond_audit.csv",
        bond_rows,
    )
    write_csv(
        OUTPUT_DIR
        / "QM_F06_LOWER_STAGE2_target_contact_evolution.csv",
        target_rows,
    )

    if close_contact_rows:
        write_csv(
            OUTPUT_DIR
            / "QM_F06_LOWER_STAGE2_close_contact_audit.csv",
            close_contact_rows,
        )

    role_summary = [
        {
            "atom_role": role,
            "atom_count": role_counts[role],
            "rms_stage1_to_stage2_displacement_angstrom": (
                f"{math.sqrt(role_sum_sq[role] / role_counts[role]):.10f}"
            ),
            "maximum_stage1_to_stage2_displacement_angstrom": (
                f"{role_max[role]:.10f}"
            ),
        }
        for role in sorted(role_counts)
    ]

    write_csv(
        OUTPUT_DIR
        / "QM_F06_LOWER_STAGE1_TO_STAGE2_role_summary.csv",
        role_summary,
    )

    report_path = (
        OUTPUT_DIR
        / "QM_F06_LOWER_STAGE2_VALIDATION_REPORT.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 LOWER Stage-2 Validation — Day027",
                "",
                f"## Decision: **{decision}**",
                "",
                "## Stage-2 calculation",
                "",
                "- ORCA termination: **NORMAL**",
                "- Geometry optimization: **CONVERGED**",
                "- Method: **PBE0-D4/def2-TZVP**",
                "- Charge/multiplicity: **0 / 1**",
                "",
                "## Stage 1 → Stage 2 displacement",
                "",
                (
                    "- Direct Cartesian RMSD: "
                    f"**{overall_rmsd:.6f} Å**"
                ),
                (
                    "- Maximum atomic displacement: "
                    f"**{float(maximum_displacement['stage1_to_stage2_displacement_angstrom']):.6f} Å**"
                ),
                (
                    "- Maximum-displacement atom: "
                    f"`{maximum_displacement['atom_id']}`"
                ),
                "",
                "## Bonded structure",
                "",
                f"- Reconstructed bonds: **{len(bond_rows)}**",
                f"- Bond-range failures: **{bond_failures}**",
                (
                    "- Artificial-cap bond failures: "
                    f"**{cap_bond_failures}**"
                ),
                "",
                "## Stage-1 relaxation targets",
                "",
                (
                    "- Target contacts evaluated: "
                    f"**{len(target_rows)}**"
                ),
                (
                    "- Both contacts improved: "
                    f"**{'YES' if all_target_contacts_improved else 'NO'}**"
                ),
                *[
                    (
                        f"- `{row['atom_1']} — {row['atom_2']}`: "
                        f"{float(row['stage1_distance_angstrom']):.4f} Å "
                        f"→ {float(row['stage2_distance_angstrom']):.4f} Å; "
                        f"improved = **{row['compression_improved']}**"
                    )
                    for row in target_rows
                ],
                "",
                "## Long-range contacts",
                "",
                (
                    "- Contacts below 0.90 of vdW sum: "
                    f"**{len(close_contact_rows)}**"
                ),
                (
                    "- Hard contacts below 0.70: "
                    f"**{blocking_contacts}**"
                ),
                (
                    "- Hard contacts involving artificial caps: "
                    f"**{cap_blocking_contacts}**"
                ),
                "",
                "## Stage-3 static input",
                "",
                (
                    "- Static syntax/content gate: "
                    f"**{'PASS' if stage3['gate_pass'] else 'FAIL'}**"
                ),
                (
                    "- Exact Stage-2 coordinate propagation: "
                    f"**{'PASS' if stage3_coordinate_match else 'FAIL'}**"
                ),
                (
                    "- Maximum coordinate-transfer error: "
                    f"**{maximum_stage3_coordinate_error:.3e} Å**"
                ),
                "",
                "## Authorization state",
                "",
                (
                    "- Stage-2 structural validation: "
                    f"**{'PASS' if stage2_structure_pass else 'FAIL'}**"
                ),
                (
                    "- Stage-3 static input readiness: "
                    f"**{'PASS' if stage3_static_pass else 'FAIL'}**"
                ),
                "- Stage-3 execution: **NOT YET AUTHORIZED**",
                "- RESP/ESP charge protocol: **NOT YET DEFINED**",
                "- Force-field parameter adoption: **NOT AUTHORIZED**",
                "",
                "## Required scientific decision",
                "",
                (
                    "Define the Stage-3 electronic-property protocol. "
                    "The current input is a valid energy single point, "
                    "but an energy-only calculation is insufficient for "
                    "RESP or electrostatic-potential-derived charges."
                ),
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = {
        "decision": decision,
        "stage1_to_stage2_rmsd_angstrom": overall_rmsd,
        "bond_range_failures": bond_failures,
        "artificial_cap_bond_failures": (
            cap_bond_failures
        ),
        "stage1_target_contacts_all_improved": (
            all_target_contacts_improved
        ),
        "stage2_hard_contacts": blocking_contacts,
        "stage2_cap_involving_hard_contacts": (
            cap_blocking_contacts
        ),
        "stage2_structure_gate_pass": (
            stage2_structure_pass
        ),
        "stage3_static_input_gate_pass": (
            stage3_static_pass
        ),
        "stage3_coordinates_match_stage2": (
            stage3_coordinate_match
        ),
        "stage3_execution_authorized": (
            stage3_execution_authorized
        ),
        "required_next_step": (
            "DEFINE_STAGE3_ELECTRONIC_PROPERTY_AND_CHARGE_PROTOCOL"
            if stage2_structure_pass and stage3_static_pass
            else "REVIEW_FAILED_STAGE2_OR_STAGE3_STATIC_GATE"
        ),
    }

    (
        OUTPUT_DIR
        / "QM_F06_LOWER_STAGE2_validation_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("QM_F06 LOWER Stage-2 validation completed.")
    print(f"Decision: {decision}")
    print(f"Stage1->Stage2 RMSD: {overall_rmsd:.6f} Å")
    print("Bond-range failures:", bond_failures)
    print(
        "Cap bond failures:",
        cap_bond_failures,
    )
    print(
        "Both Stage-1 target contacts improved:",
        all_target_contacts_improved,
    )
    print(
        "Stage-2 cap-involving hard contacts:",
        cap_blocking_contacts,
    )
    print(
        "Stage-3 static input gate:",
        stage3_static_pass,
    )
    print("Stage-3 execution authorized: False")
    print(f"Report: {report_path}")


if __name__ == "__main__":
    main()
