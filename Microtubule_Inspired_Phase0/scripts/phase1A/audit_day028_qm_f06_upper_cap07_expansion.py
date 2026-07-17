#!/usr/bin/env python3
"""
Targeted topology audit for the QM_F06 UPPER HCAP07 boundary.

The audit determines whether HCAP:UPPER:07 can be replaced by the real
R2 atom A:UPPER:13:3 and identifies the minimum chemically complete
first-shell expansion around that restored atom.

No atoms are added, removed, moved or optimized.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict, deque
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

F06_DIR = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/QM_F06"
)

UPPER_MANIFEST = F06_DIR / (
    "QM_F06_UPPER_CAPPED_REPAIRED_atoms.csv"
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

FULL_COORDINATES = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "28_r2_inner_h_reflected_direction_refinement/"
    "r2_selected_four_atom_refined_full_coordinates.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "upper_cap07_expansion_audit"
)

CAP_ID = "HCAP:UPPER:07"
INSIDE_PARENT = "A:UPPER:14:4"
RESTORED_CENTER = "A:UPPER:13:3"


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"No rows in {path}")
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
            {field: row.get(field, "") for field in fields}
            for row in rows
        )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fragment_rows = read_csv(UPPER_MANIFEST)
    node_rows = read_csv(FULL_NODES)
    edge_rows = read_csv(FULL_EDGES)
    coordinate_rows = read_csv(FULL_COORDINATES)

    fragment_map = {
        row["atom_id"]: row
        for row in fragment_rows
    }

    node_map = {
        row["node_id"]: row
        for row in node_rows
    }

    coordinate_map = {
        row["node_id"]: row
        for row in coordinate_rows
    }

    if CAP_ID not in fragment_map:
        raise RuntimeError(f"Missing cap: {CAP_ID}")

    if INSIDE_PARENT not in fragment_map:
        raise RuntimeError(
            f"Missing inside parent: {INSIDE_PARENT}"
        )

    if RESTORED_CENTER not in node_map:
        raise RuntimeError(
            f"Missing restored center in graph: {RESTORED_CENTER}"
        )

    adjacency: dict[str, set[str]] = defaultdict(set)
    edge_map: dict[tuple[str, str], dict[str, str]] = {}

    for row in edge_rows:
        source = row["source_node"]
        target = row["target_node"]

        adjacency[source].add(target)
        adjacency[target].add(source)

        edge_map[tuple(sorted((source, target)))] = row

    expected_cut = tuple(
        sorted((INSIDE_PARENT, RESTORED_CENTER))
    )

    if expected_cut not in edge_map:
        raise RuntimeError(
            "Expected real R2 cut edge is absent: "
            f"{INSIDE_PARENT} -- {RESTORED_CENTER}"
        )

    fragment_real_ids = {
        row["atom_id"]
        for row in fragment_rows
        if row.get("artificial_cap", "").lower() != "true"
    }

    # Shell 0: restored center.
    # Shell 1: all direct real graph neighbors.
    shell_0 = {RESTORED_CENTER}
    shell_1 = set(adjacency[RESTORED_CENTER])

    proposed_real_additions = (
        shell_0 | shell_1
    ) - fragment_real_ids

    # Include hydrogen passivants directly connected to the proposed
    # heavy atoms. These are already contained in shell_1 when present,
    # but this explicit classification is retained for traceability.
    hydrogen_additions = {
        atom_id
        for atom_id in proposed_real_additions
        if node_map.get(atom_id, {}).get("element") == "H"
    }

    heavy_additions = proposed_real_additions - hydrogen_additions

    # New cuts created after including shell 0 + shell 1.
    proposed_real_set = (
        fragment_real_ids | proposed_real_additions
    )

    new_boundary_rows: list[dict[str, Any]] = []

    for inside in sorted(proposed_real_set):
        for outside in sorted(adjacency.get(inside, set())):
            if outside in proposed_real_set:
                continue

            edge = edge_map[tuple(sorted((inside, outside)))]

            new_boundary_rows.append(
                {
                    "inside_node": inside,
                    "inside_element": (
                        node_map.get(inside, {}).get("element", "")
                    ),
                    "inside_node_type": (
                        node_map.get(inside, {}).get("node_type", "")
                    ),
                    "outside_node": outside,
                    "outside_element": (
                        node_map.get(outside, {}).get("element", "")
                    ),
                    "outside_node_type": (
                        node_map.get(outside, {}).get("node_type", "")
                    ),
                    "edge_id": edge["edge_id"],
                    "edge_type": edge["edge_type"],
                    "cut_touches_new_region": (
                        inside in proposed_real_additions
                    ),
                    "candidate_for_new_cap": (
                        node_map.get(inside, {}).get("element")
                        in {"B", "N"}
                        and node_map.get(outside, {}).get("element")
                        in {"B", "N"}
                    ),
                }
            )

    addition_rows: list[dict[str, Any]] = []

    for atom_id in sorted(proposed_real_additions):
        node = node_map.get(atom_id, {})
        coordinate = coordinate_map.get(atom_id)

        if coordinate is None:
            raise RuntimeError(
                f"No validated coordinate for proposed atom: {atom_id}"
            )

        addition_rows.append(
            {
                "atom_id": atom_id,
                "element": node.get("element", ""),
                "node_type": node.get("node_type", ""),
                "end": node.get("end", ""),
                "graph_degree": len(adjacency[atom_id]),
                "neighbors": "|".join(sorted(adjacency[atom_id])),
                "shell": (
                    0 if atom_id == RESTORED_CENTER else 1
                ),
                "heavy_atom": (
                    node.get("element", "") in {"B", "N"}
                ),
                "hydrogen": (
                    node.get("element", "") == "H"
                ),
                "coordinate_available": True,
                "x_nm": coordinate["x_nm"],
                "y_nm": coordinate["y_nm"],
                "z_nm": coordinate["z_nm"],
            }
        )

    target_rows = []

    for atom_id in (
        CAP_ID,
        INSIDE_PARENT,
        RESTORED_CENTER,
        *sorted(shell_1),
    ):
        source = fragment_map.get(atom_id, node_map.get(atom_id, {}))

        target_rows.append(
            {
                "atom_id": atom_id,
                "element": source.get("element", ""),
                "node_type": source.get("node_type", ""),
                "currently_in_fragment": atom_id in fragment_map,
                "artificial_cap": (
                    fragment_map.get(
                        atom_id,
                        {},
                    ).get("artificial_cap", "").lower()
                    == "true"
                ),
                "graph_neighbors": "|".join(
                    sorted(adjacency.get(atom_id, set()))
                ),
            }
        )

    write_csv(
        OUTPUT_DIR / "QM_F06_UPPER_CAP07_target_environment.csv",
        target_rows,
    )

    write_csv(
        OUTPUT_DIR / "QM_F06_UPPER_CAP07_proposed_real_additions.csv",
        addition_rows,
    )

    write_csv(
        OUTPUT_DIR / "QM_F06_UPPER_CAP07_new_boundary_edges.csv",
        new_boundary_rows,
    )

    candidate_new_cap_edges = [
        row
        for row in new_boundary_rows
        if (
            row["cut_touches_new_region"] is True
            and row["candidate_for_new_cap"] is True
        )
    ]

    summary = {
        "decision": (
            "QM_F06_UPPER_CAP07_REAL_COORDINATION_"
            "EXPANSION_DEFINED_CONSTRUCTION_PENDING"
        ),
        "cap_to_remove": CAP_ID,
        "inside_parent": INSIDE_PARENT,
        "restored_center": RESTORED_CENTER,
        "proposed_real_additions": sorted(
            proposed_real_additions
        ),
        "heavy_additions": sorted(heavy_additions),
        "hydrogen_additions": sorted(hydrogen_additions),
        "proposed_real_addition_count": len(
            proposed_real_additions
        ),
        "candidate_new_cap_edge_count": len(
            candidate_new_cap_edges
        ),
        "all_coordinates_available": True,
        "upper_v2_construction_authorized": True,
        "upper_qm_execution_authorized": False,
    }

    (
        OUTPUT_DIR / "QM_F06_UPPER_CAP07_expansion_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    report_path = (
        OUTPUT_DIR
        / "QM_F06_UPPER_CAP07_EXPANSION_AUDIT.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 UPPER HCAP07 Expansion Audit — Day028",
                "",
                f"- Artificial cap selected for removal: `{CAP_ID}`",
                f"- Current real parent: `{INSIDE_PARENT}`",
                f"- Restored real center: `{RESTORED_CENTER}`",
                "",
                "## Proposed real additions",
                "",
                *[
                    f"- `{atom_id}`"
                    for atom_id in sorted(proposed_real_additions)
                ],
                "",
                "## Counts",
                "",
                (
                    f"- Real atoms added: "
                    f"**{len(proposed_real_additions)}**"
                ),
                f"- Heavy atoms added: **{len(heavy_additions)}**",
                (
                    f"- Existing real hydrogens added: "
                    f"**{len(hydrogen_additions)}**"
                ),
                (
                    "- New peripheral B–N cuts requiring review: "
                    f"**{len(candidate_new_cap_edges)}**"
                ),
                "",
                "## Decision",
                "",
                (
                    "**QM_F06_UPPER_CAP07_REAL_COORDINATION_"
                    "EXPANSION_DEFINED_CONSTRUCTION_PENDING**"
                ),
                "",
                "## Authorization state",
                "",
                "- UPPER V2 construction: **AUTHORIZED**",
                "- Artificial-cap placement: **PENDING CUT AUDIT**",
                "- ORCA input preparation: **NOT AUTHORIZED**",
                "- ORCA execution: **NOT AUTHORIZED**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print("QM_F06 UPPER HCAP07 expansion audit completed.")
    print("Cap to remove:", CAP_ID)
    print("Restored center:", RESTORED_CENTER)
    print(
        "Proposed real additions:",
        len(proposed_real_additions),
    )
    print("Heavy additions:", len(heavy_additions))
    print(
        "Existing real H additions:",
        len(hydrogen_additions),
    )
    print(
        "Candidate new cap edges:",
        len(candidate_new_cap_edges),
    )
    print("UPPER V2 construction authorized: True")
    print("QM execution authorized: False")
    print("Report:", report_path)


if __name__ == "__main__":
    main()
