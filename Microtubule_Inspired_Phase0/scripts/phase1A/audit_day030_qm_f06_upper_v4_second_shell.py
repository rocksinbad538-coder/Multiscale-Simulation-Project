#!/usr/bin/env python3
"""
Audit the second heavy-atom shell required for QM_F06 UPPER V4.

Adds the five heavy atoms exposed by the mandatory first-shell
restoration and determines the resulting external heavy-heavy cuts.

No geometry is constructed and no QM execution is authorized.
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
    "day030_qm_f06_upper_v4_second_shell"
)

REPORT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_SECOND_SHELL_AUDIT.json"
)

CUTS_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_second_shell_cuts.csv"
)

DEFECTIVE_CAPS = {
    "HCAP:UPPER:01",
    "HCAP:UPPER:04",
}

FIRST_SHELL_RESTORATION = {
    "P:1640",
    "P:1581",
    "P:1583",
    "S:1739",
    "P:1639",
    "H4:UPPER:0203:0",
}

SECOND_HEAVY_SHELL = {
    "P:1580",
    "P:1582",
    "P:1638",
    "P:1642",
    "S:1738",
}


def read_csv(path: Path):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle)
        return list(reader)


def pair(a: str, b: str) -> tuple[str, str]:
    return tuple(sorted((a, b)))


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

    adjacency: dict[str, set[str]] = defaultdict(set)
    edge_data = {}

    for row in edge_rows:
        a = row["source_node"]
        b = row["target_node"]

        adjacency[a].add(b)
        adjacency[b].add(a)
        edge_data[pair(a, b)] = row

    current_v3 = {
        row["atom_id"]
        for row in v3_rows
    }

    retained_v3 = current_v3 - DEFECTIVE_CAPS

    first_shell_model = (
        retained_v3
        | FIRST_SHELL_RESTORATION
    )

    second_shell_model = (
        first_shell_model
        | SECOND_HEAVY_SHELL
    )

    missing = (
        FIRST_SHELL_RESTORATION
        | SECOND_HEAVY_SHELL
    ) - set(nodes)

    if missing:
        raise RuntimeError(
            "Required atoms missing from canonical graph: "
            f"{sorted(missing)}"
        )

    def heavy_cuts(included: set[str]):
        records = {}

        for atom_id in sorted(included):
            if atom_id not in nodes:
                continue

            if nodes[atom_id]["element"] == "H":
                continue

            for neighbor in sorted(adjacency[atom_id]):
                if neighbor in included:
                    continue

                if nodes[neighbor]["element"] == "H":
                    continue

                edge_pair = pair(atom_id, neighbor)
                edge = edge_data[edge_pair]

                records[edge_pair] = {
                    "included_atom": atom_id,
                    "included_element": (
                        nodes[atom_id]["element"]
                    ),
                    "outside_atom": neighbor,
                    "outside_element": (
                        nodes[neighbor]["element"]
                    ),
                    "outside_node_type": (
                        nodes[neighbor]["node_type"]
                    ),
                    "edge_id": edge["edge_id"],
                    "edge_type": edge["edge_type"],
                }

        return records

    first_cuts = heavy_cuts(first_shell_model)
    second_cuts = heavy_cuts(second_shell_model)

    first_pairs = set(first_cuts)
    second_pairs = set(second_cuts)

    closed_by_second_shell = (
        first_pairs - second_pairs
    )

    newly_created_by_second_shell = (
        second_pairs - first_pairs
    )

    persistent = first_pairs & second_pairs

    rows = []

    for edge_pair in sorted(
        first_pairs | second_pairs
    ):
        if edge_pair in closed_by_second_shell:
            classification = "CLOSED_BY_SECOND_SHELL"
            record = first_cuts[edge_pair]
        elif edge_pair in newly_created_by_second_shell:
            classification = "NEW_AFTER_SECOND_SHELL"
            record = second_cuts[edge_pair]
        else:
            classification = "PERSISTENT_PREEXISTING_CUT"
            record = second_cuts[edge_pair]

        rows.append({
            **record,
            "classification": classification,
        })

    fieldnames = [
        "included_atom",
        "included_element",
        "outside_atom",
        "outside_element",
        "outside_node_type",
        "edge_id",
        "edge_type",
        "classification",
    ]

    with CUTS_PATH.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(rows)

    new_records = [
        second_cuts[value]
        for value in sorted(
            newly_created_by_second_shell
        )
    ]

    new_outside_atoms = sorted({
        row["outside_atom"]
        for row in new_records
    })

    included_boundary_atoms = sorted({
        row["included_atom"]
        for row in new_records
    })

    report = {
        "decision": (
            "QM_F06_UPPER_V4_SECOND_HEAVY_SHELL_"
            "AUDITED_FINAL_CLOSURE_PENDING"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "defective_caps_removed": sorted(
            DEFECTIVE_CAPS
        ),
        "first_shell_restoration": sorted(
            FIRST_SHELL_RESTORATION
        ),
        "second_heavy_shell": sorted(
            SECOND_HEAVY_SHELL
        ),
        "counts": {
            "retained_v3_atoms": len(
                retained_v3
            ),
            "first_shell_model_atoms": len(
                first_shell_model
            ),
            "second_shell_model_atoms": len(
                second_shell_model
            ),
            "first_shell_heavy_cuts": len(
                first_pairs
            ),
            "second_shell_heavy_cuts": len(
                second_pairs
            ),
            "cuts_closed_by_second_shell": len(
                closed_by_second_shell
            ),
            "new_cuts_after_second_shell": len(
                newly_created_by_second_shell
            ),
            "persistent_preexisting_cuts": len(
                persistent
            ),
        },
        "new_cuts_after_second_shell": (
            new_records
        ),
        "newly_exposed_heavy_atoms": (
            new_outside_atoms
        ),
        "new_boundary_atoms_inside_model": (
            included_boundary_atoms
        ),
        "authorization": {
            "v4_geometry_construction_authorized": False,
            "orca_execution_authorized": False,
            "geometric_reference_accepted": False,
            "electronic_reference_accepted": False,
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
    print("QM_F06 UPPER V4 SECOND-SHELL AUDIT")
    print("=" * 78)

    print(
        "Second-shell atoms:",
        sorted(SECOND_HEAVY_SHELL),
    )
    print(
        "First-shell heavy cuts:",
        len(first_pairs),
    )
    print(
        "Second-shell heavy cuts:",
        len(second_pairs),
    )
    print(
        "Cuts closed by second shell:",
        len(closed_by_second_shell),
    )
    print(
        "New cuts after second shell:",
        len(newly_created_by_second_shell),
    )
    print(
        "New boundary atoms inside model:",
        included_boundary_atoms,
    )
    print(
        "Newly exposed outside atoms:",
        new_outside_atoms,
    )

    print()
    print("NEW HEAVY-HEAVY CUTS AFTER SECOND SHELL")

    for row in new_records:
        print(
            f"{row['included_atom']:24s} -- "
            f"{row['outside_atom']:24s} "
            f"[{row['edge_id']}]"
        )

    print()
    print("Decision:", report["decision"])
    print("Report:", REPORT_PATH)
    print("Cuts:", CUTS_PATH)
    print("V4 construction authorized: False")


if __name__ == "__main__":
    main()
