#!/usr/bin/env python3
"""
Final topology-aware steric audit of repaired QM_F06 capped geometries.

This script:
- reads the REPAIRED atom manifests;
- reconstructs the original fragment topology;
- excludes 1-2, 1-3 and 1-4 pairs from the hard vdW gate;
- separates cap-induced failures from contacts inherited from R2;
- verifies X-H bond lengths after cap repositioning;
- does not modify coordinates or execute QM calculations.
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

ENDS = ("LOWER", "UPPER")

VDW_RADII = {
    "H": 1.20,
    "B": 1.92,
    "N": 1.55,
}

HARD_FAILURE_CLASSES = {
    "SEVERE_CLASH",
    "STRONG_COMPRESSION",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(f"No rows found in {path}")

    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "fragment",
        "atom_1",
        "element_1",
        "role_1",
        "atom_2",
        "element_2",
        "role_2",
        "graph_separation",
        "distance_angstrom",
        "vdw_sum_angstrom",
        "distance_over_vdw_sum",
        "classification",
        "involves_artificial_cap",
        "contact_origin",
        "blocks_qm_input_preparation",
    ]

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()

        for row in rows:
            writer.writerow(row)


def coordinates(
    atom: dict[str, str],
) -> tuple[float, float, float]:
    return (
        float(atom["x_angstrom"]),
        float(atom["y_angstrom"]),
        float(atom["z_angstrom"]),
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


def classify(ratio: float) -> str:
    if ratio < 0.70:
        return "SEVERE_CLASH"

    if ratio < 0.80:
        return "STRONG_COMPRESSION"

    if ratio < 0.90:
        return "CLOSE_CONTACT"

    return "ACCEPTABLE"


def main() -> None:
    combined_rows: list[dict[str, Any]] = []
    report_sections: list[str] = []

    all_cap_gates_pass = True
    all_bond_length_gates_pass = True

    for end in ENDS:
        original_label = f"QM_F06_{end}_CAPPED"
        repaired_label = f"{original_label}_REPAIRED"

        atoms = read_csv(
            F06_DIR / f"{repaired_label}_atoms.csv"
        )

        internal_edges = read_csv(
            F06_DIR / f"QM_F06_{end}_internal_edges.csv"
        )

        boundary_audit = read_csv(
            F06_DIR / f"QM_F06_{end}_boundary_edge_audit.csv"
        )

        caps = read_csv(
            F06_DIR / f"{original_label}_caps.csv"
        )

        atom_lookup = {
            row["atom_id"]: row
            for row in atoms
        }

        adjacency: dict[str, set[str]] = defaultdict(set)

        def add_edge(first: str, second: str) -> None:
            adjacency[first].add(second)
            adjacency[second].add(first)

        for edge in internal_edges:
            add_edge(
                edge["source_node"],
                edge["target_node"],
            )

        for edge in boundary_audit:
            if (
                edge["preliminary_action"]
                == "INCLUDE_EXISTING_HYDROGEN"
            ):
                add_edge(
                    edge["inside_node"],
                    edge["outside_node"],
                )

        for cap in caps:
            add_edge(
                cap["parent_inside_node"],
                cap["cap_id"],
            )

        # Verify repaired X-H distances.
        cap_distance_rows: list[dict[str, Any]] = []
        cap_distance_failures = 0

        for cap in caps:
            cap_id = cap["cap_id"]
            parent_id = cap["parent_inside_node"]

            measured = distance(
                coordinates(atom_lookup[parent_id]),
                coordinates(atom_lookup[cap_id]),
            )

            target = float(
                cap["target_XH_distance_angstrom"]
            )
            error = abs(measured - target)
            passed = error <= 1.0e-8

            if not passed:
                cap_distance_failures += 1

            cap_distance_rows.append(
                {
                    "fragment": repaired_label,
                    "cap_id": cap_id,
                    "parent_atom": parent_id,
                    "target_distance_angstrom": f"{target:.10f}",
                    "measured_distance_angstrom": f"{measured:.10f}",
                    "absolute_error_angstrom": f"{error:.12e}",
                    "distance_gate_pass": passed,
                }
            )

        all_bond_length_gates_pass = (
            all_bond_length_gates_pass
            and cap_distance_failures == 0
        )

        cap_distance_path = (
            F06_DIR
            / f"{repaired_label}_cap_distance_validation.csv"
        )

        with cap_distance_path.open(
            "w",
            newline="",
            encoding="utf-8",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(cap_distance_rows[0].keys()),
            )
            writer.writeheader()
            writer.writerows(cap_distance_rows)

        atom_ids = sorted(atom_lookup)
        audit_rows: list[dict[str, Any]] = []

        excluded_counts = Counter()

        for index, atom_1_id in enumerate(atom_ids):
            for atom_2_id in atom_ids[index + 1:]:
                separation = shortest_path_length(
                    adjacency,
                    atom_1_id,
                    atom_2_id,
                )

                if separation in {1, 2, 3}:
                    excluded_counts[separation] += 1
                    continue

                atom_1 = atom_lookup[atom_1_id]
                atom_2 = atom_lookup[atom_2_id]

                measured = distance(
                    coordinates(atom_1),
                    coordinates(atom_2),
                )

                vdw_sum = (
                    VDW_RADII[atom_1["element"]]
                    + VDW_RADII[atom_2["element"]]
                )

                ratio = measured / vdw_sum
                classification = classify(ratio)

                involves_cap = (
                    atom_1["artificial_cap"].lower() == "true"
                    or atom_2["artificial_cap"].lower() == "true"
                )

                contact_origin = (
                    "CAP_INTRODUCED_OR_CAP_INVOLVING"
                    if involves_cap
                    else "INHERITED_FROM_R2_GEOMETRY"
                )

                blocks_qm = (
                    involves_cap
                    and classification in HARD_FAILURE_CLASSES
                )

                if ratio < 1.0:
                    row = {
                        "fragment": repaired_label,
                        "atom_1": atom_1_id,
                        "element_1": atom_1["element"],
                        "role_1": atom_1["atom_role"],
                        "atom_2": atom_2_id,
                        "element_2": atom_2["element"],
                        "role_2": atom_2["atom_role"],
                        "graph_separation": (
                            separation
                            if separation is not None
                            else "DISCONNECTED"
                        ),
                        "distance_angstrom": f"{measured:.10f}",
                        "vdw_sum_angstrom": f"{vdw_sum:.10f}",
                        "distance_over_vdw_sum": f"{ratio:.10f}",
                        "classification": classification,
                        "involves_artificial_cap": involves_cap,
                        "contact_origin": contact_origin,
                        "blocks_qm_input_preparation": blocks_qm,
                    }

                    audit_rows.append(row)
                    combined_rows.append(row)

        class_counts = Counter(
            row["classification"]
            for row in audit_rows
        )

        inherited_hard = [
            row
            for row in audit_rows
            if (
                row["classification"] in HARD_FAILURE_CLASSES
                and not row["involves_artificial_cap"]
            )
        ]

        cap_hard = [
            row
            for row in audit_rows
            if row["blocks_qm_input_preparation"]
        ]

        cap_gate_pass = len(cap_hard) == 0

        all_cap_gates_pass = (
            all_cap_gates_pass and cap_gate_pass
        )

        output_path = (
            F06_DIR
            / f"{repaired_label}_final_steric_audit.csv"
        )
        write_csv(output_path, audit_rows)

        report_sections.extend(
            [
                f"## {repaired_label}",
                "",
                f"- Excluded 1–2 pairs: **{excluded_counts[1]}**",
                f"- Excluded 1–3 pairs: **{excluded_counts[2]}**",
                f"- Excluded 1–4 pairs: **{excluded_counts[3]}**",
                (
                    "- Long-range contacts below vdW sum: "
                    f"**{len(audit_rows)}**"
                ),
                (
                    "- Classification counts: "
                    f"`{dict(sorted(class_counts.items()))}`"
                ),
                (
                    "- Hard contacts inherited from R2: "
                    f"**{len(inherited_hard)}**"
                ),
                (
                    "- Hard contacts involving artificial caps: "
                    f"**{len(cap_hard)}**"
                ),
                (
                    "- Repaired cap-distance failures: "
                    f"**{cap_distance_failures}**"
                ),
                (
                    "- Artificial-cap steric gate: "
                    f"**{'PASS' if cap_gate_pass else 'FAIL'}**"
                ),
                "",
            ]
        )

    qm_input_authorized = (
        all_cap_gates_pass
        and all_bond_length_gates_pass
    )

    decision = (
        "QM_F06_REPAIRED_FRAGMENTS_READY_FOR_QM_INPUT_PREPARATION"
        if qm_input_authorized
        else
        "QM_F06_REPAIRED_FRAGMENTS_REQUIRE_FURTHER_BOUNDARY_REDESIGN"
    )

    combined_path = (
        F06_DIR
        / "QM_F06_repaired_fragment_final_steric_audit_combined.csv"
    )
    write_csv(combined_path, combined_rows)

    report_path = (
        F06_DIR
        / "QM_F06_REPAIRED_FRAGMENT_FINAL_AUDIT.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 Repaired-Fragment Final Audit — Day026",
                "",
                "## Decision",
                "",
                f"**{decision}**",
                "",
                "## Gate logic",
                "",
                (
                    "Hard long-range contacts involving artificial caps "
                    "are blocking. Hard contacts composed only of original "
                    "R2 atoms are retained as inherited conformational "
                    "strain and must be monitored during QM optimization."
                ),
                "",
                *report_sections,
                "## Authorization state",
                "",
                (
                    "- Artificial-cap geometry gate: "
                    f"**{'PASSED' if all_cap_gates_pass else 'FAILED'}**"
                ),
                (
                    "- Repaired X–H distance gate: "
                    f"**{'PASSED' if all_bond_length_gates_pass else 'FAILED'}**"
                ),
                (
                    "- QM input preparation: "
                    f"**{'AUTHORIZED' if qm_input_authorized else 'NOT AUTHORIZED'}**"
                ),
                "- QM calculation execution: **NOT AUTHORIZED**",
                "",
                "## Required next step",
                "",
                (
                    "Prepare reproducible electronic-structure input files "
                    "for the LOWER and UPPER repaired fragments if both "
                    "boundary gates pass. Execution remains a separate "
                    "authorization decision."
                ),
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = {
        "decision": decision,
        "artificial_cap_steric_gate_pass": all_cap_gates_pass,
        "repaired_xh_distance_gate_pass": (
            all_bond_length_gates_pass
        ),
        "qm_input_preparation_authorized": (
            qm_input_authorized
        ),
        "qm_calculation_execution_authorized": False,
        "required_next_step": (
            "PREPARE_QM_INPUT_FILES"
            if qm_input_authorized
            else "REDESIGN_FRAGMENT_BOUNDARY"
        ),
    }

    (
        F06_DIR
        / "QM_F06_repaired_fragment_final_audit_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("Final repaired-fragment audit completed.")
    print(f"Decision: {decision}")
    print(
        "QM input preparation authorized:",
        qm_input_authorized,
    )
    print("QM calculation execution authorized: False")


if __name__ == "__main__":
    main()
