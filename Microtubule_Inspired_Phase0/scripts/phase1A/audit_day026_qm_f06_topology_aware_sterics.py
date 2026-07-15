#!/usr/bin/env python3
"""
Topology-aware steric audit for capped QM_F06 fragments.

Pairs separated by:
- 1 bond (1-2),
- 2 bonds (1-3),
- 3 bonds (1-4)

are excluded from the hard van der Waals clash gate.

Only pairs with graph separation >= 4, or disconnected pairs, are
evaluated using vdW-radius normalization.

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

F06_DIR = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/QM_F06"
)

ENDS = ("LOWER", "UPPER")

VDW_RADII = {
    "H": 1.20,
    "B": 1.92,
    "N": 1.55,
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
        "requires_geometry_repair",
    ]

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()

        for row in rows:
            writer.writerow(row)


def xyz(atom: dict[str, str]) -> tuple[float, float, float]:
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
    fragment_passes: list[bool] = []

    for end in ENDS:
        label = f"QM_F06_{end}_CAPPED"

        atoms = read_csv(
            F06_DIR / f"{label}_atoms.csv"
        )
        internal_edges = read_csv(
            F06_DIR / f"QM_F06_{end}_internal_edges.csv"
        )
        boundary_audit = read_csv(
            F06_DIR / f"QM_F06_{end}_boundary_edge_audit.csv"
        )
        caps = read_csv(
            F06_DIR / f"{label}_caps.csv"
        )

        atom_lookup = {
            row["atom_id"]: row
            for row in atoms
        }

        adjacency: dict[str, set[str]] = defaultdict(set)

        def add_edge(first: str, second: str) -> None:
            adjacency[first].add(second)
            adjacency[second].add(first)

        for row in internal_edges:
            add_edge(row["source_node"], row["target_node"])

        for row in boundary_audit:
            if (
                row["preliminary_action"]
                == "INCLUDE_EXISTING_HYDROGEN"
            ):
                add_edge(
                    row["inside_node"],
                    row["outside_node"],
                )

        for row in caps:
            add_edge(
                row["parent_inside_node"],
                row["cap_id"],
            )

        atom_ids = sorted(atom_lookup)
        audit_rows: list[dict[str, Any]] = []

        excluded_12 = 0
        excluded_13 = 0
        excluded_14 = 0

        for index, atom_1_id in enumerate(atom_ids):
            for atom_2_id in atom_ids[index + 1:]:
                separation = shortest_path_length(
                    adjacency,
                    atom_1_id,
                    atom_2_id,
                )

                if separation == 1:
                    excluded_12 += 1
                    continue
                if separation == 2:
                    excluded_13 += 1
                    continue
                if separation == 3:
                    excluded_14 += 1
                    continue

                atom_1 = atom_lookup[atom_1_id]
                atom_2 = atom_lookup[atom_2_id]

                measured = distance(
                    xyz(atom_1),
                    xyz(atom_2),
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

                requires_repair = (
                    classification
                    in {
                        "SEVERE_CLASH",
                        "STRONG_COMPRESSION",
                    }
                )

                if ratio < 1.0:
                    row = {
                        "fragment": label,
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
                        "requires_geometry_repair": requires_repair,
                    }

                    audit_rows.append(row)
                    combined_rows.append(row)

        class_counts = Counter(
            row["classification"]
            for row in audit_rows
        )

        hard_failures = [
            row
            for row in audit_rows
            if row["requires_geometry_repair"]
        ]

        cap_failures = [
            row
            for row in hard_failures
            if row["involves_artificial_cap"]
        ]

        fragment_pass = len(hard_failures) == 0
        fragment_passes.append(fragment_pass)

        write_csv(
            F06_DIR
            / f"{label}_topology_aware_steric_audit.csv",
            audit_rows,
        )

        report_sections.extend(
            [
                f"## {label}",
                "",
                f"- Excluded 1–2 pairs: **{excluded_12}**",
                f"- Excluded 1–3 pairs: **{excluded_13}**",
                f"- Excluded 1–4 pairs: **{excluded_14}**",
                (
                    "- Long-range contacts below vdW sum: "
                    f"**{len(audit_rows)}**"
                ),
                (
                    "- Classification counts: "
                    f"`{dict(sorted(class_counts.items()))}`"
                ),
                f"- Hard steric failures: **{len(hard_failures)}**",
                (
                    "- Hard failures involving artificial caps: "
                    f"**{len(cap_failures)}**"
                ),
                (
                    "- Topology-aware steric gate: "
                    f"**{'PASS' if fragment_pass else 'FAIL'}**"
                ),
                "",
            ]
        )

    overall_pass = all(fragment_passes)

    decision = (
        "QM_F06_CAPPED_FRAGMENTS_PASS_TOPOLOGY_AWARE_STERIC_GATE"
        if overall_pass
        else
        "QM_F06_CAPPED_FRAGMENTS_REQUIRE_TARGETED_GEOMETRY_REPAIR"
    )

    combined_path = (
        F06_DIR
        / "QM_F06_topology_aware_steric_audit_combined.csv"
    )
    write_csv(combined_path, combined_rows)

    report_path = (
        F06_DIR
        / "QM_F06_TOPOLOGY_AWARE_STERIC_AUDIT.md"
    )
    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 Topology-Aware Steric Audit — Day026",
                "",
                "## Correction to the previous audit",
                "",
                (
                    "The previous pair-specific audit incorrectly applied "
                    "a van der Waals overlap criterion to bonded 1–2, "
                    "angular 1–3 and torsional 1–4 intramolecular pairs."
                ),
                "",
                (
                    "The present audit excludes all pairs with graph "
                    "separation ≤3 and applies the hard steric gate only "
                    "to longer-range or disconnected atom pairs."
                ),
                "",
                f"## Decision: **{decision}**",
                "",
                *report_sections,
                "## Authorization state",
                "",
                "- Graph and valence gate: **PASSED**",
                "- Topology-aware steric gate: "
                f"**{'PASSED' if overall_pass else 'FAILED'}**",
                "- QM input preparation: "
                f"**{'AUTHORIZED' if overall_pass else 'NOT AUTHORIZED'}**",
                "- QM calculation executed: **NO**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = {
        "decision": decision,
        "topology_aware_steric_gate_pass": overall_pass,
        "qm_input_preparation_authorized": overall_pass,
        "qm_calculation_executed": False,
        "required_next_step": (
            "PREPARE_QM_INPUTS"
            if overall_pass
            else "REPAIR_ONLY_REMAINING_LONG_RANGE_STERIC_FAILURES"
        ),
    }

    (
        F06_DIR
        / "QM_F06_topology_aware_steric_audit_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("Topology-aware QM_F06 steric audit completed.")
    print(f"Decision: {decision}")
    print("QM input preparation authorized:", overall_pass)


if __name__ == "__main__":
    main()
