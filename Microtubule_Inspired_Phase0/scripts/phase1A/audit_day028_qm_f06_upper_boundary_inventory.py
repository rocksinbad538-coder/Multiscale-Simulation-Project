#!/usr/bin/env python3
"""
Inventory and boundary audit for the existing QM_F06 UPPER fragment.

The script:

- locates the current UPPER repaired atom manifest;
- reconstructs its atom set and graph boundary;
- classifies artificial caps and their real parent atoms;
- identifies real R2 atoms immediately outside the fragment;
- reports first-shell boundary-expansion candidates;
- does not add, remove or move atoms;
- does not prepare or execute QM calculations.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

F06_DIR = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/QM_F06"
)

FULL_NODES = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

FULL_EDGES = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_edges.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "upper_boundary_inventory"
)


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
        writer.writerows(
            {
                field: row.get(field, "")
                for field in fields
            }
            for row in rows
        )


def find_upper_manifest() -> Path:
    candidates = sorted(
        path
        for path in F06_DIR.glob("*UPPER*CAPPED*REPAIRED*atoms.csv")
        if path.is_file()
    )

    if not candidates:
        candidates = sorted(
            path
            for path in F06_DIR.glob("*UPPER*atoms.csv")
            if path.is_file()
        )

    if not candidates:
        raise RuntimeError(
            "No UPPER atom manifest found in "
            f"{F06_DIR}"
        )

    return candidates[-1]


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    manifest_path = find_upper_manifest()

    fragment_rows = read_csv(manifest_path)
    node_rows = read_csv(FULL_NODES)
    edge_rows = read_csv(FULL_EDGES)

    node_map = {
        row["node_id"]: row
        for row in node_rows
    }

    fragment_ids = {
        row["atom_id"]
        for row in fragment_rows
    }

    real_fragment_ids = {
        row["atom_id"]
        for row in fragment_rows
        if row.get("artificial_cap", "").lower() != "true"
    }

    artificial_caps = [
        row
        for row in fragment_rows
        if row.get("artificial_cap", "").lower() == "true"
    ]

    adjacency: dict[str, set[str]] = defaultdict(set)

    for row in edge_rows:
        source = row["source_node"]
        target = row["target_node"]

        adjacency[source].add(target)
        adjacency[target].add(source)

    boundary_rows: list[dict[str, Any]] = []

    for inside in sorted(real_fragment_ids):
        for outside in sorted(adjacency.get(inside, set())):
            if outside in real_fragment_ids:
                continue

            if outside.startswith("HCAP"):
                continue

            outside_record = node_map.get(outside, {})

            boundary_rows.append(
                {
                    "inside_node": inside,
                    "inside_element": (
                        node_map.get(inside, {}).get(
                            "element",
                            next(
                                (
                                    row["element"]
                                    for row in fragment_rows
                                    if row["atom_id"] == inside
                                ),
                                "",
                            ),
                        )
                    ),
                    "inside_node_type": (
                        node_map.get(inside, {}).get(
                            "node_type",
                            next(
                                (
                                    row.get("node_type", "")
                                    for row in fragment_rows
                                    if row["atom_id"] == inside
                                ),
                                "",
                            ),
                        )
                    ),
                    "outside_node": outside,
                    "outside_element": outside_record.get(
                        "element",
                        "",
                    ),
                    "outside_node_type": outside_record.get(
                        "node_type",
                        "",
                    ),
                    "outside_end": outside_record.get(
                        "end",
                        "",
                    ),
                    "outside_already_in_manifest": (
                        outside in fragment_ids
                    ),
                    "candidate_real_boundary_expansion": (
                        outside in node_map
                        and outside not in fragment_ids
                    ),
                }
            )

    candidate_ids = sorted(
        {
            row["outside_node"]
            for row in boundary_rows
            if row["candidate_real_boundary_expansion"] is True
        }
    )

    candidate_rows: list[dict[str, Any]] = []

    for candidate in candidate_ids:
        record = node_map[candidate]

        internal_neighbors = sorted(
            neighbor
            for neighbor in adjacency[candidate]
            if neighbor in real_fragment_ids
        )

        external_neighbors = sorted(
            neighbor
            for neighbor in adjacency[candidate]
            if neighbor not in real_fragment_ids
        )

        candidate_rows.append(
            {
                "candidate_node": candidate,
                "element": record.get("element", ""),
                "node_type": record.get("node_type", ""),
                "end": record.get("end", ""),
                "graph_degree": len(adjacency[candidate]),
                "neighbors_inside_current_fragment": "|".join(
                    internal_neighbors
                ),
                "neighbors_outside_current_fragment": "|".join(
                    external_neighbors
                ),
                "inside_neighbor_count": len(internal_neighbors),
                "outside_neighbor_count": len(external_neighbors),
            }
        )

    cap_rows = []

    for cap in artificial_caps:
        cap_rows.append(
            {
                "cap_id": cap["atom_id"],
                "element": cap["element"],
                "parent_inside_node": cap.get(
                    "parent_inside_node",
                    "",
                ),
                "source_edge_id": cap.get(
                    "source_edge_id",
                    "",
                ),
                "atom_role": cap.get("atom_role", ""),
                "node_type": cap.get("node_type", ""),
            }
        )

    element_counts = Counter(
        row["element"]
        for row in fragment_rows
    )

    role_counts = Counter(
        row.get("atom_role", "")
        for row in fragment_rows
    )

    node_type_counts = Counter(
        row.get("node_type", "")
        for row in fragment_rows
    )

    write_csv(
        OUTPUT_DIR / "QM_F06_UPPER_boundary_edges.csv",
        boundary_rows,
    )

    write_csv(
        OUTPUT_DIR / "QM_F06_UPPER_expansion_candidates.csv",
        candidate_rows,
    )

    write_csv(
        OUTPUT_DIR / "QM_F06_UPPER_artificial_caps.csv",
        cap_rows,
    )

    summary = {
        "decision": (
            "QM_F06_UPPER_BOUNDARY_INVENTORY_COMPLETED_"
            "CONSTRUCTION_NOT_AUTHORIZED"
        ),
        "manifest": str(manifest_path.relative_to(ROOT)),
        "fragment_atom_count": len(fragment_rows),
        "element_counts": dict(element_counts),
        "role_counts": dict(role_counts),
        "node_type_counts": dict(node_type_counts),
        "artificial_cap_count": len(artificial_caps),
        "boundary_edge_count": len(boundary_rows),
        "unique_real_expansion_candidates": len(candidate_rows),
        "upper_boundary_construction_authorized": False,
        "upper_qm_execution_authorized": False,
    }

    (
        OUTPUT_DIR / "QM_F06_UPPER_boundary_inventory_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    report_path = (
        OUTPUT_DIR / "QM_F06_UPPER_BOUNDARY_INVENTORY_REPORT.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 UPPER Boundary Inventory — Day028",
                "",
                f"- Manifest: `{manifest_path.relative_to(ROOT)}`",
                f"- Current atoms: **{len(fragment_rows)}**",
                f"- Element counts: `{dict(element_counts)}`",
                f"- Artificial caps: **{len(artificial_caps)}**",
                f"- Boundary edges: **{len(boundary_rows)}**",
                (
                    "- Unique real expansion candidates: "
                    f"**{len(candidate_rows)}**"
                ),
                "",
                "## Decision",
                "",
                (
                    "**QM_F06_UPPER_BOUNDARY_INVENTORY_COMPLETED_"
                    "CONSTRUCTION_NOT_AUTHORIZED**"
                ),
                "",
                "The inventory identifies possible real R2 atoms that "
                "could replace or displace artificial boundary caps. "
                "No candidate is selected solely because it is the "
                "geometric counterpart of a LOWER atom."
                "",
                "## Authorization state",
                "",
                "- UPPER inventory: **COMPLETED**",
                "- UPPER boundary construction: **NOT AUTHORIZED**",
                "- UPPER QM input preparation: **NOT AUTHORIZED**",
                "- UPPER QM execution: **NOT AUTHORIZED**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print("QM_F06 UPPER boundary inventory completed.")
    print("Manifest:", manifest_path)
    print("Current atoms:", len(fragment_rows))
    print("Artificial caps:", len(artificial_caps))
    print("Boundary edges:", len(boundary_rows))
    print("Expansion candidates:", len(candidate_rows))
    print("UPPER construction authorized: False")
    print("Report:", report_path)


if __name__ == "__main__":
    main()
