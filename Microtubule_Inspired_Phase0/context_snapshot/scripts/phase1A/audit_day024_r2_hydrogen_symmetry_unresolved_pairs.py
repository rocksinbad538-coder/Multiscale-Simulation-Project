#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/phase1A/day024_chemical_end_rim_design"

GATE3M = BASE / "16_r2_selected_full_density_longer_bn_bridge_graph"
GATE3P1 = BASE / "24_r2_hydrogen_symmetry_refinement_preflight"
OUT = BASE / "25_r2_hydrogen_symmetry_unresolved_pair_audit"

GRAPH_NODES = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

GRAPH_EDGES = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_edges.csv"
)

GRAPH_PATHS = (
    GATE3M
    / "r2_selected_longer_bn_bridge_paths.csv"
)

HEAVY_PAIRS = (
    GATE3P1
    / "r2_heavy_lower_upper_pair_candidates.csv"
)

H_PAIRS = (
    GATE3P1
    / "r2_hydrogen_lower_upper_pair_candidates.csv"
)

UNRESOLVED_HEAVY = (
    OUT
    / "r2_unresolved_heavy_pair_metadata.csv"
)

UNRESOLVED_H = (
    OUT
    / "r2_unresolved_hydrogen_pair_metadata.csv"
)

PATH_METADATA = (
    OUT
    / "r2_bridge_path_pairing_metadata.csv"
)

SUMMARY = (
    OUT
    / "r2_unresolved_symmetry_pair_audit_summary.csv"
)

JSON_OUT = (
    OUT
    / "r2_unresolved_symmetry_pair_audit.json"
)

REPORT = (
    OUT
    / "R2_UNRESOLVED_SYMMETRY_PAIR_AUDIT_DAY024.md"
)


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty required file: {path}"
        )


def read_rows(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(
            f"No rows found in {path}"
        )

    return rows


def write_rows(
    path: Path,
    rows: list[dict[str, Any]],
) -> None:
    if not rows:
        raise RuntimeError(
            f"No rows available for {path}"
        )

    fields: list[str] = []

    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(
                {
                    field: row.get(field, "")
                    for field in fields
                }
            )


def nonempty_metadata(
    row: dict[str, str],
) -> dict[str, str]:
    excluded = {
        "x_nm",
        "y_nm",
        "z_nm",
    }

    return {
        key: value
        for key, value in row.items()
        if (
            key not in excluded
            and value not in {
                "",
                None,
            }
        )
    }


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        GRAPH_NODES,
        GRAPH_EDGES,
        GRAPH_PATHS,
        HEAVY_PAIRS,
        H_PAIRS,
    ):
        require_file(required)

    node_rows = read_rows(
        GRAPH_NODES
    )

    edge_rows = read_rows(
        GRAPH_EDGES
    )

    path_rows = read_rows(
        GRAPH_PATHS
    )

    heavy_pair_rows = read_rows(
        HEAVY_PAIRS
    )

    H_pair_rows = read_rows(
        H_PAIRS
    )

    nodes = {
        row["node_id"]: row
        for row in node_rows
    }

    adjacency = {
        node_id: set()
        for node_id in nodes
    }

    for row in edge_rows:
        first = row["source_node"]
        second = row["target_node"]

        adjacency[first].add(second)
        adjacency[second].add(first)

    unresolved_heavy_ids = [
        row["lower_node"]
        for row in heavy_pair_rows
        if row.get(
            "upper_node",
            "",
        ) == ""
    ]

    unresolved_H_ids = [
        row["lower_H"]
        for row in H_pair_rows
        if row.get(
            "pair_status"
        ) != "UNIQUE"
    ]

    unresolved_heavy_rows = []

    for node_id in unresolved_heavy_ids:
        row = nodes[node_id]

        heavy_neighbors = sorted(
            neighbor
            for neighbor in adjacency[node_id]
            if nodes[neighbor]["element"] != "H"
        )

        hydrogen_neighbors = sorted(
            neighbor
            for neighbor in adjacency[node_id]
            if nodes[neighbor]["element"] == "H"
        )

        unresolved_heavy_rows.append(
            {
                **nonempty_metadata(row),
                "heavy_degree": len(
                    heavy_neighbors
                ),
                "H_degree": len(
                    hydrogen_neighbors
                ),
                "heavy_neighbors": (
                    " | ".join(
                        heavy_neighbors
                    )
                ),
                "H_neighbors": (
                    " | ".join(
                        hydrogen_neighbors
                    )
                ),
            }
        )

    unresolved_H_rows = []

    for node_id in unresolved_H_ids:
        row = nodes[node_id]

        heavy_neighbors = [
            neighbor
            for neighbor in adjacency[node_id]
            if nodes[neighbor]["element"] != "H"
        ]

        attached_heavy = (
            heavy_neighbors[0]
            if len(heavy_neighbors) == 1
            else ""
        )

        attached_row = (
            nodes[attached_heavy]
            if attached_heavy
            else {}
        )

        unresolved_H_rows.append(
            {
                **nonempty_metadata(row),
                "attached_heavy_node": (
                    attached_heavy
                ),
                "attached_heavy_element": (
                    attached_row.get(
                        "element",
                        "",
                    )
                ),
                "attached_heavy_type": (
                    attached_row.get(
                        "node_type",
                        "",
                    )
                ),
                "attached_heavy_end": (
                    attached_row.get(
                        "end",
                        "",
                    )
                ),
                "attached_heavy_metadata": (
                    json.dumps(
                        nonempty_metadata(
                            attached_row
                        ),
                        sort_keys=True,
                    )
                    if attached_row
                    else ""
                ),
            }
        )

    write_rows(
        UNRESOLVED_HEAVY,
        unresolved_heavy_rows,
    )

    write_rows(
        UNRESOLVED_H,
        unresolved_H_rows,
    )

    path_metadata_rows = []

    for row in path_rows:
        path_metadata_rows.append(
            {
                **row,
                "seed_node_type": nodes.get(
                    row.get(
                        "seed_node",
                        "",
                    ),
                    {},
                ).get(
                    "node_type",
                    "",
                ),
                "annulus_node_type": nodes.get(
                    row.get(
                        "annulus_node",
                        "",
                    ),
                    {},
                ).get(
                    "node_type",
                    "",
                ),
            }
        )

    write_rows(
        PATH_METADATA,
        path_metadata_rows,
    )

    heavy_type_counts = Counter(
        row["node_type"]
        for row in unresolved_heavy_rows
    )

    H_role_counts = Counter(
        row["node_type"]
        for row in unresolved_H_rows
    )

    H_parent_type_counts = Counter(
        row["attached_heavy_type"]
        for row in unresolved_H_rows
    )

    lower_path_rows = [
        row
        for row in path_rows
        if row.get("end") == "LOWER"
    ]

    upper_path_rows = [
        row
        for row in path_rows
        if row.get("end") == "UPPER"
    ]

    summary = {
        "unresolved_heavy_count": len(
            unresolved_heavy_ids
        ),
        "unresolved_H_count": len(
            unresolved_H_ids
        ),
        "unresolved_heavy_types": (
            " | ".join(
                f"{key}:{value}"
                for key, value
                in sorted(
                    heavy_type_counts.items()
                )
            )
        ),
        "unresolved_H_roles": (
            " | ".join(
                f"{key}:{value}"
                for key, value
                in sorted(
                    H_role_counts.items()
                )
            )
        ),
        "unresolved_H_attached_heavy_types": (
            " | ".join(
                f"{key}:{value}"
                for key, value
                in sorted(
                    H_parent_type_counts.items()
                )
            )
        ),
        "lower_bridge_paths": len(
            lower_path_rows
        ),
        "upper_bridge_paths": len(
            upper_path_rows
        ),
        "path_based_correspondence_required": any(
            "BRIDGE"
            in key
            for key in heavy_type_counts
        ),
        "coordinates_modified": False,
        "molecular_topology_generated": False,
        "energy_minimized": False,
        "MD_performed": False,
        "required_next_step": (
            "BUILD_TOPOLOGY_AWARE_LOWER_UPPER_CORRESPONDENCE"
        ),
    }

    write_rows(
        SUMMARY,
        [summary],
    )

    JSON_OUT.write_text(
        json.dumps(
            {
                "summary": summary,
                "unresolved_heavy_type_counts": dict(
                    heavy_type_counts
                ),
                "unresolved_H_role_counts": dict(
                    H_role_counts
                ),
                "unresolved_H_parent_type_counts": dict(
                    H_parent_type_counts
                ),
                "limitations": [
                    (
                        "This audit does not modify coordinates."
                    ),
                    (
                        "No topology, parameterization, minimization, "
                        "MD or QM calculation is generated."
                    ),
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    heavy_lines = "\n".join(
        f"- {key}: **{value}**"
        for key, value
        in sorted(
            heavy_type_counts.items()
        )
    )

    H_lines = "\n".join(
        f"- {key}: **{value}**"
        for key, value
        in sorted(
            H_role_counts.items()
        )
    )

    parent_lines = "\n".join(
        f"- {key}: **{value}**"
        for key, value
        in sorted(
            H_parent_type_counts.items()
        )
    )

    REPORT.write_text(
        f"""# R2 Unresolved Symmetry-Pair Audit

## Unresolved heavy nodes

- Total: **{len(unresolved_heavy_ids)}**

{heavy_lines}

## Unresolved hydrogen nodes

- Total: **{len(unresolved_H_ids)}**

{H_lines}

## Attached-heavy classes for unresolved H

{parent_lines}

## Bridge paths

- Lower paths: **{len(lower_path_rows)}**
- Upper paths: **{len(upper_path_rows)}**
- Path-based correspondence required:
  **{summary['path_based_correspondence_required']}**

## Restrictions

- Coordinates modified: **NO**
- Molecular topology generated: **NO**
- Energy minimization performed: **NO**
- MD performed: **NO**

## Required next step

`{summary['required_next_step']}`
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 unresolved symmetry-pair audit completed."
    )

    print(
        "Unresolved heavy count: "
        f"{len(unresolved_heavy_ids)}"
    )

    print(
        "Unresolved heavy types: "
        + (
            "NONE"
            if not heavy_type_counts
            else " | ".join(
                f"{key}:{value}"
                for key, value
                in sorted(
                    heavy_type_counts.items()
                )
            )
        )
    )

    print(
        "Unresolved H count: "
        f"{len(unresolved_H_ids)}"
    )

    print(
        "Unresolved H roles: "
        + (
            "NONE"
            if not H_role_counts
            else " | ".join(
                f"{key}:{value}"
                for key, value
                in sorted(
                    H_role_counts.items()
                )
            )
        )
    )

    print(
        "Unresolved H attached-heavy types: "
        + (
            "NONE"
            if not H_parent_type_counts
            else " | ".join(
                f"{key}:{value}"
                for key, value
                in sorted(
                    H_parent_type_counts.items()
                )
            )
        )
    )

    print(
        "Bridge paths lower/upper: "
        f"{len(lower_path_rows)}/"
        f"{len(upper_path_rows)}"
    )

    print(
        "Path-based correspondence required: "
        f"{summary['path_based_correspondence_required']}"
    )

    print(
        "Coordinates modified: NO"
    )

    print(
        "Molecular topology generated: NO"
    )

    print(
        "Energy minimization performed: NO"
    )

    print(
        "MD performed: NO"
    )

    print(
        "Required next step: "
        f"{summary['required_next_step']}"
    )

    for path in (
        UNRESOLVED_HEAVY,
        UNRESOLVED_H,
        PATH_METADATA,
        SUMMARY,
        JSON_OUT,
        REPORT,
    ):
        print(
            f"Wrote: {path.relative_to(ROOT)}"
        )


if __name__ == "__main__":
    main()
