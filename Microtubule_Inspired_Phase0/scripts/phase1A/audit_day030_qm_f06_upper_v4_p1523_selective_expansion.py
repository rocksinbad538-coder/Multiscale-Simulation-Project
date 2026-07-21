#!/usr/bin/env python3
"""
Compare direct B-H closure against selective inclusion of P:1523.

P:1523 is shared by P:1580 and P:1582. The audit determines:
- cuts closed by adding P:1523;
- new cuts introduced around P:1523;
- artificial H count before and after inclusion;
- whether the selective expansion reduces geminal BH2 groups;
- whether additional N-H or heavy-shell closure would be required.

No geometry is constructed.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

GRAPH_ROOT = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph"
)

NODES_PATH = (
    GRAPH_ROOT
    / "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

EDGES_PATH = (
    GRAPH_ROOT
    / "r2_selected_longer_bn_bridge_graph_edges.csv"
)

V3_MAP = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V3/orca_v3a2_workflow/"
    "v3a2_atom_role_constraint_map.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day030_qm_f06_upper_v4_p1523_selective_expansion"
)

REPORT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_P1523_SELECTIVE_EXPANSION.json"
)

COMPARISON_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_P1523_boundary_comparison.csv"
)

DEFECTIVE_CAPS = {
    "HCAP:UPPER:01",
    "HCAP:UPPER:04",
}

FIRST_SHELL = {
    "P:1640",
    "P:1581",
    "P:1583",
    "S:1739",
    "P:1639",
    "H4:UPPER:0203:0",
}

SECOND_SHELL = {
    "P:1580",
    "P:1582",
    "P:1638",
    "P:1642",
    "S:1738",
}

SELECTIVE_ATOM = "P:1523"


def read_csv(path: Path):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    node_rows = read_csv(NODES_PATH)
    edge_rows = read_csv(EDGES_PATH)
    v3_rows = read_csv(V3_MAP)

    nodes = {
        row["node_id"]: {
            "element": row["element"],
            "node_type": row["node_type"],
        }
        for row in node_rows
    }

    adjacency = defaultdict(set)

    for row in edge_rows:
        a = row["source_node"]
        b = row["target_node"]
        adjacency[a].add(b)
        adjacency[b].add(a)

    retained_v3 = {
        row["atom_id"]
        for row in v3_rows
    } - DEFECTIVE_CAPS

    direct_model = (
        retained_v3
        | FIRST_SHELL
        | SECOND_SHELL
    )

    selective_model = (
        direct_model
        | {SELECTIVE_ATOM}
    )

    if SELECTIVE_ATOM not in nodes:
        raise RuntimeError(
            f"Missing selective atom: {SELECTIVE_ATOM}"
        )

    def boundary_summary(
        included,
        candidate_centers,
    ):
        records = []

        for atom_id in sorted(candidate_centers):
            if atom_id not in included:
                continue

            if atom_id not in nodes:
                continue

            heavy_neighbors = {
                neighbor
                for neighbor in adjacency[atom_id]
                if nodes[neighbor]["element"] != "H"
            }

            retained_heavy = (
                heavy_neighbors & included
            )

            outside_heavy = (
                heavy_neighbors - included
            )

            if not outside_heavy:
                continue

            element = nodes[atom_id]["element"]
            h_caps = len(outside_heavy)

            records.append({
                "center": atom_id,
                "element": element,
                "retained_heavy_degree": len(
                    retained_heavy
                ),
                "outside_heavy_degree": len(
                    outside_heavy
                ),
                "outside_heavy_neighbors": "|".join(
                    sorted(outside_heavy)
                ),
                "artificial_h_caps_required": h_caps,
                "nominal_coordination_after_caps": (
                    len(retained_heavy) + h_caps
                ),
                "geminal_double_H_cap": h_caps == 2,
            })

        return records

    direct_candidate_centers = set(
        SECOND_SHELL
    )

    selective_candidate_centers = (
        set(SECOND_SHELL)
        | {SELECTIVE_ATOM}
    )

    direct_records = boundary_summary(
        direct_model,
        direct_candidate_centers,
    )

    selective_records = boundary_summary(
        selective_model,
        selective_candidate_centers,
    )

    def summarize(records):
        return {
            "boundary_center_count": len(records),
            "artificial_h_caps": sum(
                row["artificial_h_caps_required"]
                for row in records
            ),
            "double_H_boundary_centers": sorted(
                row["center"]
                for row in records
                if row["geminal_double_H_cap"]
            ),
            "nitrogen_boundary_centers": sorted(
                row["center"]
                for row in records
                if row["element"] == "N"
            ),
        }

    direct_summary = summarize(
        direct_records
    )

    selective_summary = summarize(
        selective_records
    )

    p1523_record = next(
        (
            row for row in selective_records
            if row["center"] == SELECTIVE_ATOM
        ),
        None,
    )

    comparison_rows = []

    for strategy, records in (
        ("DIRECT_BH_CLOSURE", direct_records),
        ("P1523_SELECTIVE_EXPANSION", selective_records),
    ):
        for row in records:
            comparison_rows.append({
                "strategy": strategy,
                **row,
            })

    fields = [
        "strategy",
        "center",
        "element",
        "retained_heavy_degree",
        "outside_heavy_degree",
        "outside_heavy_neighbors",
        "artificial_h_caps_required",
        "nominal_coordination_after_caps",
        "geminal_double_H_cap",
    ]

    with COMPARISON_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
        )
        writer.writeheader()
        writer.writerows(comparison_rows)

    artificial_h_reduction = (
        direct_summary["artificial_h_caps"]
        - selective_summary["artificial_h_caps"]
    )

    report = {
        "decision": (
            "QM_F06_UPPER_V4_P1523_SELECTIVE_"
            "EXPANSION_AUDITED"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "selective_atom": SELECTIVE_ATOM,
        "selective_atom_element": nodes[
            SELECTIVE_ATOM
        ]["element"],
        "selective_atom_canonical_neighbors": sorted(
            adjacency[SELECTIVE_ATOM]
        ),
        "direct_closure": direct_summary,
        "selective_expansion": selective_summary,
        "p1523_boundary_record": p1523_record,
        "artificial_h_reduction": (
            artificial_h_reduction
        ),
        "authorization": {
            "preferred_strategy_selected": False,
            "v4_geometry_construction_authorized": False,
            "orca_execution_authorized": False,
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
    print("QM_F06 UPPER V4 P:1523 SELECTIVE-EXPANSION AUDIT")
    print("=" * 78)
    print(
        "P:1523 element:",
        nodes[SELECTIVE_ATOM]["element"],
    )
    print(
        "P:1523 canonical neighbors:",
        sorted(adjacency[SELECTIVE_ATOM]),
    )

    print()
    print("DIRECT B-H CLOSURE")
    print(
        "Artificial H:",
        direct_summary["artificial_h_caps"],
    )
    print(
        "Double-H centers:",
        direct_summary[
            "double_H_boundary_centers"
        ],
    )
    print(
        "N boundary centers:",
        direct_summary[
            "nitrogen_boundary_centers"
        ],
    )

    print()
    print("P:1523 SELECTIVE EXPANSION")
    print(
        "Artificial H:",
        selective_summary["artificial_h_caps"],
    )
    print(
        "Double-H centers:",
        selective_summary[
            "double_H_boundary_centers"
        ],
    )
    print(
        "N boundary centers:",
        selective_summary[
            "nitrogen_boundary_centers"
        ],
    )
    print(
        "P:1523 boundary record:",
        p1523_record,
    )
    print(
        "Artificial-H reduction:",
        artificial_h_reduction,
    )

    print()
    print("Decision:", report["decision"])
    print("Report:", REPORT_PATH)
    print("Comparison:", COMPARISON_CSV)
    print("Strategy selected: False")
    print("V4 construction authorized: False")


if __name__ == "__main__":
    main()
