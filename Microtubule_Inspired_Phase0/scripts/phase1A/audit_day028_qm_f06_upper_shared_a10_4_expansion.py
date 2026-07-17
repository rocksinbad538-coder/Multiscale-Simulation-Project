#!/usr/bin/env python3
"""
Audit the shared omitted atom A:UPPER:10:4 responsible for the
HCAP:UPPER:05 / HCAPV2:UPPER:02 overlap.

No coordinates are modified and no QM calculation is executed.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

V2_DIR = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V2"
)

ATOMS_PATH = V2_DIR / "QM_F06_UPPER_BOUNDARY_V2_atoms.csv"

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

FULL_COORDS = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "28_r2_inner_h_reflected_direction_refinement/"
    "r2_selected_four_atom_refined_full_coordinates.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "upper_shared_a10_4_expansion_audit"
)

SHARED_ATOM = "A:UPPER:10:4"

CAPS_TO_REMOVE = {
    "HCAP:UPPER:05",
    "HCAPV2:UPPER:02",
}

KNOWN_INSIDE_NEIGHBORS = {
    "A:UPPER:11:3",
    "A:UPPER:11:5",
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(f"No records in {path}")

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

        for row in rows:
            writer.writerow(
                {
                    field: row.get(field, "")
                    for field in fields
                }
            )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fragment_rows = read_csv(ATOMS_PATH)
    node_rows = read_csv(FULL_NODES)
    edge_rows = read_csv(FULL_EDGES)
    coordinate_rows = read_csv(FULL_COORDS)

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

    adjacency: dict[str, set[str]] = defaultdict(set)
    edge_map: dict[tuple[str, str], dict[str, str]] = {}

    for row in edge_rows:
        first = row["source_node"]
        second = row["target_node"]

        adjacency[first].add(second)
        adjacency[second].add(first)

        edge_map[tuple(sorted((first, second)))] = row

    if SHARED_ATOM not in node_map:
        raise RuntimeError(
            f"Shared atom missing from graph: {SHARED_ATOM}"
        )

    if SHARED_ATOM not in coordinate_map:
        raise RuntimeError(
            f"Shared atom lacks validated coordinates: {SHARED_ATOM}"
        )

    observed_neighbors = set(adjacency[SHARED_ATOM])

    if not KNOWN_INSIDE_NEIGHBORS.issubset(observed_neighbors):
        raise RuntimeError(
            "Expected inside neighbors absent. "
            f"Observed: {sorted(observed_neighbors)}"
        )

    third_neighbors = sorted(
        observed_neighbors - KNOWN_INSIDE_NEIGHBORS
    )

    if len(third_neighbors) != 1:
        raise RuntimeError(
            "Expected exactly one third neighbor for "
            f"{SHARED_ATOM}; found {third_neighbors}"
        )

    third_neighbor = third_neighbors[0]

    additions = {
        SHARED_ATOM,
    }

    if third_neighbor not in fragment_map:
        additions.add(third_neighbor)

    # Include an existing real hydrogen if the third neighbor is H.
    third_neighbor_is_real_h = (
        node_map[third_neighbor]["element"] == "H"
    )

    proposed_real_set = {
        row["atom_id"]
        for row in fragment_rows
        if row["atom_id"] not in CAPS_TO_REMOVE
        and row["artificial_cap"].lower() != "true"
    } | additions

    new_cut_rows = []

    for inside in sorted(additions):
        for outside in sorted(adjacency[inside]):
            if outside in proposed_real_set:
                continue

            edge = edge_map[
                tuple(sorted((inside, outside)))
            ]

            inside_element = node_map[inside]["element"]
            outside_element = node_map[outside]["element"]

            new_cut_rows.append(
                {
                    "inside_node": inside,
                    "inside_element": inside_element,
                    "inside_node_type": node_map[inside]["node_type"],
                    "outside_node": outside,
                    "outside_element": outside_element,
                    "outside_node_type": node_map[outside]["node_type"],
                    "edge_id": edge["edge_id"],
                    "edge_type": edge["edge_type"],
                    "candidate_for_artificial_cap": (
                        {inside_element, outside_element}
                        == {"B", "N"}
                    ),
                    "outside_is_real_hydrogen": (
                        outside_element == "H"
                    ),
                }
            )

    environment_rows = []

    for atom_id in [
        SHARED_ATOM,
        *sorted(observed_neighbors),
    ]:
        node = node_map[atom_id]
        coordinate = coordinate_map.get(atom_id)

        environment_rows.append(
            {
                "atom_id": atom_id,
                "element": node["element"],
                "node_type": node["node_type"],
                "end": node["end"],
                "currently_in_fragment": atom_id in fragment_map,
                "selected_real_addition": atom_id in additions,
                "graph_degree": len(adjacency[atom_id]),
                "neighbors": "|".join(
                    sorted(adjacency[atom_id])
                ),
                "coordinate_available": coordinate is not None,
                "x_nm": coordinate["x_nm"] if coordinate else "",
                "y_nm": coordinate["y_nm"] if coordinate else "",
                "z_nm": coordinate["z_nm"] if coordinate else "",
            }
        )

    write_csv(
        OUTPUT_DIR
        / "QM_F06_UPPER_A10_4_environment.csv",
        environment_rows,
    )

    write_csv(
        OUTPUT_DIR
        / "QM_F06_UPPER_A10_4_new_cut_edges.csv",
        new_cut_rows,
    )

    candidate_cap_edges = [
        row
        for row in new_cut_rows
        if row["candidate_for_artificial_cap"] is True
    ]

    summary = {
        "decision": (
            "QM_F06_UPPER_SHARED_A10_4_EXPANSION_AUDITED_"
            "BOUNDARY_V3_CONSTRUCTION_PENDING"
        ),
        "shared_atom": SHARED_ATOM,
        "caps_to_remove": sorted(CAPS_TO_REMOVE),
        "shared_atom_neighbors": sorted(observed_neighbors),
        "third_neighbor": third_neighbor,
        "third_neighbor_element": (
            node_map[third_neighbor]["element"]
        ),
        "third_neighbor_node_type": (
            node_map[third_neighbor]["node_type"]
        ),
        "third_neighbor_already_in_fragment": (
            third_neighbor in fragment_map
        ),
        "third_neighbor_is_real_hydrogen": (
            third_neighbor_is_real_h
        ),
        "proposed_real_additions": sorted(additions),
        "new_candidate_cap_edges": candidate_cap_edges,
        "boundary_v3_construction_authorized": True,
        "orca_input_preparation_authorized": False,
        "orca_execution_authorized": False,
    }

    (
        OUTPUT_DIR
        / "QM_F06_UPPER_A10_4_expansion_summary.json"
    ).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    report_path = (
        OUTPUT_DIR
        / "QM_F06_UPPER_A10_4_EXPANSION_AUDIT.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 UPPER Shared A10:4 Expansion Audit — Day028",
                "",
                f"- Shared omitted atom: `{SHARED_ATOM}`",
                "",
                "## Existing real neighbors",
                "",
                *[
                    f"- `{atom_id}`"
                    for atom_id in sorted(KNOWN_INSIDE_NEIGHBORS)
                ],
                "",
                "## Third coordination partner",
                "",
                f"- Atom: `{third_neighbor}`",
                (
                    f"- Element: "
                    f"`{node_map[third_neighbor]['element']}`"
                ),
                (
                    f"- Node type: "
                    f"`{node_map[third_neighbor]['node_type']}`"
                ),
                (
                    "- Already present in fragment: "
                    f"**{'YES' if third_neighbor in fragment_map else 'NO'}**"
                ),
                "",
                "## Selected correction",
                "",
                "- Remove `HCAP:UPPER:05`",
                "- Remove `HCAPV2:UPPER:02`",
                f"- Restore `{SHARED_ATOM}`",
                *(
                    [f"- Restore `{third_neighbor}`"]
                    if third_neighbor not in fragment_map
                    else []
                ),
                "",
                "## New cut edges requiring treatment",
                "",
                *[
                    (
                        f"- `{row['edge_id']}`: "
                        f"`{row['inside_node']} — {row['outside_node']}`"
                    )
                    for row in candidate_cap_edges
                ],
                "",
                "## Decision",
                "",
                (
                    "**QM_F06_UPPER_SHARED_A10_4_EXPANSION_"
                    "AUDITED_BOUNDARY_V3_CONSTRUCTION_PENDING**"
                ),
                "",
                "## Authorization state",
                "",
                "- Boundary V3 construction: **AUTHORIZED**",
                "- Pre-QM audit: **REQUIRED AFTER CONSTRUCTION**",
                "- ORCA input preparation: **NOT AUTHORIZED**",
                "- ORCA execution: **NOT AUTHORIZED**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print("QM_F06 UPPER shared A10:4 audit completed.")
    print("Shared atom:", SHARED_ATOM)
    print("Neighbors:", sorted(observed_neighbors))
    print("Third neighbor:", third_neighbor)
    print(
        "Third neighbor already present:",
        third_neighbor in fragment_map,
    )
    print(
        "Third neighbor is real H:",
        third_neighbor_is_real_h,
    )
    print("Proposed real additions:", sorted(additions))
    print("Candidate new cap edges:", len(candidate_cap_edges))
    print("Boundary V3 construction authorized: True")
    print("ORCA execution authorized: False")
    print("Report:", report_path)


if __name__ == "__main__":
    main()
