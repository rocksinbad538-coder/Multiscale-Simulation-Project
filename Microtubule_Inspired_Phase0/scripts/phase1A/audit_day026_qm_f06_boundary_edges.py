#!/usr/bin/env python3
"""
Audit boundary edges of the extracted QM_F06 bridge fragments.

The script does not modify geometries and does not add capping atoms.
It identifies:

- the chemical identity of every cut edge;
- the local graph degrees of the inside and outside atoms;
- the original B–N bond distance;
- whether a cut touches the bridge core or an attachment center;
- whether fragment expansion must be considered before hydrogen capping.
"""

from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

COORDINATES = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "28_r2_inner_h_reflected_direction_refinement/"
    "r2_selected_four_atom_refined_full_coordinates.csv"
)

NODES = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

EDGES = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_edges.csv"
)

F06_DIR = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/QM_F06"
)

ENDS = ("LOWER", "UPPER")


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")

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
        writer.writerows(
            {
                field: row.get(field, "")
                for field in fields
            }
            for row in rows
        )


def distance_angstrom(
    first: dict[str, str],
    second: dict[str, str],
) -> float:
    dx = float(first["x_nm"]) - float(second["x_nm"])
    dy = float(first["y_nm"]) - float(second["y_nm"])
    dz = float(first["z_nm"]) - float(second["z_nm"])

    return 10.0 * math.sqrt(dx * dx + dy * dy + dz * dz)


def main() -> None:
    coordinate_rows = read_csv(COORDINATES)
    node_rows = read_csv(NODES)
    edge_rows = read_csv(EDGES)

    coordinates = {
        row["node_id"]: row
        for row in coordinate_rows
    }
    nodes = {
        row["node_id"]: row
        for row in node_rows
    }

    adjacency: dict[str, set[str]] = defaultdict(set)

    for edge in edge_rows:
        source = edge["source_node"]
        target = edge["target_node"]
        adjacency[source].add(target)
        adjacency[target].add(source)

    combined_rows: list[dict[str, Any]] = []
    report_sections: list[str] = []
    expansion_required = False

    for end in ENDS:
        atom_rows = read_csv(
            F06_DIR / f"QM_F06_{end}_atoms.csv"
        )
        boundary_rows = read_csv(
            F06_DIR / f"QM_F06_{end}_boundary_edges.csv"
        )

        fragment_atoms = {
            row["node_id"]: row
            for row in atom_rows
        }

        audited_rows: list[dict[str, Any]] = []

        for row in boundary_rows:
            inside = row["inside_node"]
            outside = row["outside_node"]

            inside_atom = fragment_atoms[inside]
            inside_meta = nodes[inside]
            outside_meta = nodes[outside]

            touches_bridge_core = (
                inside_atom["is_bridge_atom"].strip().lower() == "true"
            )
            touches_attachment_center = (
                inside_atom["is_attachment_center"].strip().lower()
                == "true"
            )

            if touches_bridge_core or touches_attachment_center:
                preliminary_action = (
                    "EXPAND_FRAGMENT_BEFORE_CAPPING"
                )
                expansion_required = True
            elif outside_meta["element"] == "H":
                preliminary_action = (
                    "INCLUDE_EXISTING_HYDROGEN"
                )
                expansion_required = True
            elif (
                inside_meta["element"] in {"B", "N"}
                and outside_meta["element"] in {"B", "N"}
            ):
                preliminary_action = (
                    "CANDIDATE_BN_CUT_FOR_HYDROGEN_CAPPING"
                )
            else:
                preliminary_action = (
                    "MANUAL_CHEMICAL_REVIEW_REQUIRED"
                )
                expansion_required = True

            audited = {
                "fragment": f"QM_F06_{end}",
                "edge_id": row["edge_id"],
                "inside_node": inside,
                "inside_element": inside_meta["element"],
                "inside_node_type": inside_meta["node_type"],
                "inside_graph_degree": len(adjacency[inside]),
                "outside_node": outside,
                "outside_element": outside_meta["element"],
                "outside_node_type": outside_meta["node_type"],
                "outside_graph_degree": len(adjacency[outside]),
                "original_edge_type": row["edge_type"],
                "original_distance_angstrom": (
                    f"{distance_angstrom(coordinates[inside], coordinates[outside]):.8f}"
                ),
                "touches_bridge_core": touches_bridge_core,
                "touches_attachment_center": touches_attachment_center,
                "preliminary_action": preliminary_action,
                "capping_authorized": False,
            }

            audited_rows.append(audited)
            combined_rows.append(audited)

        output = F06_DIR / (
            f"QM_F06_{end}_boundary_edge_audit.csv"
        )
        write_csv(output, audited_rows)

        action_counts = Counter(
            row["preliminary_action"]
            for row in audited_rows
        )
        pair_counts = Counter(
            (
                row["inside_element"],
                row["outside_element"],
            )
            for row in audited_rows
        )

        report_sections.extend(
            [
                f"## QM_F06_{end}",
                "",
                f"- Boundary edges audited: **{len(audited_rows)}**",
                f"- Element-pair counts: `{dict(sorted(pair_counts.items()))}`",
                f"- Preliminary actions: `{dict(sorted(action_counts.items()))}`",
                "",
            ]
        )

    combined_path = F06_DIR / (
        "QM_F06_boundary_edge_audit_combined.csv"
    )
    write_csv(combined_path, combined_rows)

    action_counts = Counter(
        row["preliminary_action"]
        for row in combined_rows
    )
    pair_counts = Counter(
        (
            row["inside_element"],
            row["outside_element"],
        )
        for row in combined_rows
    )

    decision = (
        "FRAGMENT_EXPANSION_REQUIRED_BEFORE_CAPPING"
        if expansion_required
        else "BN_BOUNDARY_CAPPING_DESIGN_CAN_PROCEED"
    )

    report = F06_DIR / "QM_F06_BOUNDARY_EDGE_AUDIT.md"
    report.write_text(
        "\n".join(
            [
                "# QM_F06 Boundary-Edge Audit — Day026",
                "",
                "## Scope",
                "",
                (
                    "All graph edges cut during extraction of the LOWER "
                    "and UPPER QM_F06 fragments were chemically classified."
                ),
                "",
                (
                    "No atoms were added, removed or moved. Hydrogen "
                    "capping remains unauthorized pending this audit."
                ),
                "",
                "## Combined result",
                "",
                f"- Total cut edges: **{len(combined_rows)}**",
                f"- Element-pair counts: `{dict(sorted(pair_counts.items()))}`",
                f"- Preliminary actions: `{dict(sorted(action_counts.items()))}`",
                f"- Decision: **{decision}**",
                "",
                *report_sections,
                "## Authorization state",
                "",
                "- Artificial capping authorized: **NO**",
                "- Geometry optimization authorized: **NO**",
                "- QM calculation executed: **NO**",
                "",
                "## Required next step",
                "",
                (
                    "Review every preliminary action and either enlarge "
                    "the graph fragment or define chemically valid B–H/N–H "
                    "cap placement along the original cut-bond vectors."
                ),
                "",
            ]
        ),
        encoding="utf-8",
    )

    summary = {
        "decision": decision,
        "boundary_edges_audited": len(combined_rows),
        "element_pair_counts": {
            f"{a}-{b}": count
            for (a, b), count in sorted(pair_counts.items())
        },
        "preliminary_action_counts": dict(
            sorted(action_counts.items())
        ),
        "artificial_capping_authorized": False,
        "geometry_optimization_authorized": False,
        "qm_calculation_executed": False,
        "required_next_step": (
            "REVIEW_BOUNDARY_AUDIT_AND_SELECT_EXPANSION_OR_CAP_DESIGN"
        ),
    }

    (
        F06_DIR / "QM_F06_boundary_edge_audit_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("QM_F06 boundary-edge audit completed.")
    print(f"Boundary edges audited: {len(combined_rows)}")
    print(f"Decision: {decision}")
    print(f"Output: {combined_path}")


if __name__ == "__main__":
    main()
