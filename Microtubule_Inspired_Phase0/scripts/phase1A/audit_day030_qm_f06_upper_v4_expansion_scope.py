#!/usr/bin/env python3
"""
Determine the chemically required expansion scope for QM_F06 UPPER V4.

The audit:
- loads the canonical R2 selected graph;
- loads the current V3 atom inventory;
- removes the two defective artificial caps conceptually;
- restores the real first shell around P:1641;
- identifies all newly exposed heavy-atom cuts;
- classifies real R2 hydrogens separately from artificial QM caps;
- produces reproducible JSON and CSV reports;
- does not construct a geometry or authorize QM execution.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

GRAPH_ROOT = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph"
)

GRAPH_NODES = (
    GRAPH_ROOT
    / "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

GRAPH_EDGES = (
    GRAPH_ROOT
    / "r2_selected_longer_bn_bridge_graph_edges.csv"
)

V3_MAP = ROOT / (
    "runs/phase1A/day028_qm_f06_upper_transferability/"
    "QM_F06_UPPER_BOUNDARY_V3/orca_v3a2_workflow/"
    "v3a2_atom_role_constraint_map.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_expansion_scope"
)

OUTPUT_JSON = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_EXPANSION_SCOPE.json"
)

OUTPUT_NODES = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_candidate_nodes.csv"
)

OUTPUT_EDGES = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_candidate_edges.csv"
)

CENTRAL_ATOM = "P:1641"

DEFECTIVE_CAPS = {
    "HCAP:UPPER:01",
    "HCAP:UPPER:04",
}

MANDATORY_REAL_RESTORATION = {
    "P:1640",
    "P:1581",
    "P:1583",
    "S:1739",
    "P:1639",
    "H4:UPPER:0203:0",
}

EXPECTED_CENTRAL_NEIGHBORS = {
    "P:1640",
    "S:1710",
    "S:1739",
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def read_csv(path: Path):
    require_file(path)

    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        headers = reader.fieldnames or []

    return headers, rows


def find_column(
    headers: list[str],
    candidates: tuple[str, ...],
) -> str:
    lowered = {
        header.strip().lower(): header
        for header in headers
    }

    for candidate in candidates:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]

    raise RuntimeError(
        f"Could not find any of {candidates}; "
        f"headers={headers}"
    )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    node_headers, node_rows = read_csv(GRAPH_NODES)
    edge_headers, edge_rows = read_csv(GRAPH_EDGES)
    _, v3_rows = read_csv(V3_MAP)

    node_id_key = find_column(
        node_headers,
        ("node_id", "atom_id", "id"),
    )

    element_key = find_column(
        node_headers,
        ("element",),
    )

    role_key = find_column(
        node_headers,
        ("node_role", "atom_role", "role", "node_type"),
    )

    edge_a_key = find_column(
        edge_headers,
        (
            "node_a",
            "atom_id_1",
            "source",
            "source_node",
            "atom1",
        ),
    )

    edge_b_key = find_column(
        edge_headers,
        (
            "node_b",
            "atom_id_2",
            "target",
            "target_node",
            "atom2",
        ),
    )

    edge_id_key = find_column(
        edge_headers,
        ("edge_id", "id"),
    )

    nodes = {
        row[node_id_key]: {
            "atom_id": row[node_id_key],
            "element": row[element_key],
            "role": row[role_key],
            "raw": row,
        }
        for row in node_rows
    }

    adjacency: dict[str, set[str]] = defaultdict(set)
    edge_by_pair = {}

    for row in edge_rows:
        a = row[edge_a_key]
        b = row[edge_b_key]

        adjacency[a].add(b)
        adjacency[b].add(a)

        edge_by_pair[frozenset((a, b))] = row

    current_v3_ids = {
        row["atom_id"]
        for row in v3_rows
    }

    retained_v3_ids = (
        current_v3_ids - DEFECTIVE_CAPS
    )

    missing_mandatory = (
        MANDATORY_REAL_RESTORATION - set(nodes)
    )

    if missing_mandatory:
        raise RuntimeError(
            "Mandatory restoration atoms absent from canonical graph: "
            f"{sorted(missing_mandatory)}"
        )

    proposed_v4_ids = (
        retained_v3_ids
        | MANDATORY_REAL_RESTORATION
    )

    central_neighbors = adjacency[CENTRAL_ATOM]

    central_neighbor_gate = (
        EXPECTED_CENTRAL_NEIGHBORS
        .issubset(central_neighbors)
    )

    exposed_edges = []

    for atom_id in sorted(proposed_v4_ids):
        if atom_id not in nodes:
            # Current QM-only artificial atoms are not in the
            # canonical graph and are ignored here.
            continue

        for neighbor in sorted(adjacency[atom_id]):
            if neighbor in proposed_v4_ids:
                continue

            pair = frozenset((atom_id, neighbor))
            edge = edge_by_pair[pair]

            exposed_edges.append(
                {
                    "included_atom": atom_id,
                    "included_element": nodes[atom_id]["element"],
                    "outside_atom": neighbor,
                    "outside_element": nodes[neighbor]["element"],
                    "outside_role": nodes[neighbor]["role"],
                    "edge_id": edge[edge_id_key],
                    "heavy_heavy_cut": (
                        nodes[atom_id]["element"] != "H"
                        and nodes[neighbor]["element"] != "H"
                    ),
                }
            )

    exposed_heavy_atoms = sorted({
        record["outside_atom"]
        for record in exposed_edges
        if record["heavy_heavy_cut"]
    })

    first_shell = {
        CENTRAL_ATOM,
        *adjacency[CENTRAL_ATOM],
    }

    distance_from_center = {
        CENTRAL_ATOM: 0,
    }

    queue = deque([CENTRAL_ATOM])

    while queue:
        atom = queue.popleft()
        depth = distance_from_center[atom]

        if depth >= 3:
            continue

        for neighbor in adjacency[atom]:
            if neighbor not in distance_from_center:
                distance_from_center[neighbor] = depth + 1
                queue.append(neighbor)

    node_records = []

    all_relevant = sorted(
        proposed_v4_ids
        | set(exposed_heavy_atoms)
    )

    for atom_id in all_relevant:
        canonical = nodes.get(atom_id)

        if canonical is None:
            node_records.append(
                {
                    "atom_id": atom_id,
                    "element": "",
                    "canonical_role": "QM_ONLY_OR_ARTIFICIAL",
                    "in_current_v3": atom_id in current_v3_ids,
                    "retained_in_proposed_v4": atom_id in proposed_v4_ids,
                    "mandatory_restoration": (
                        atom_id in MANDATORY_REAL_RESTORATION
                    ),
                    "exposed_next_shell": (
                        atom_id in exposed_heavy_atoms
                    ),
                    "graph_distance_from_P1641": "",
                    "canonical_degree": "",
                }
            )
            continue

        node_records.append(
            {
                "atom_id": atom_id,
                "element": canonical["element"],
                "canonical_role": canonical["role"],
                "in_current_v3": atom_id in current_v3_ids,
                "retained_in_proposed_v4": atom_id in proposed_v4_ids,
                "mandatory_restoration": (
                    atom_id in MANDATORY_REAL_RESTORATION
                ),
                "exposed_next_shell": (
                    atom_id in exposed_heavy_atoms
                ),
                "graph_distance_from_P1641": (
                    distance_from_center.get(atom_id, "")
                ),
                "canonical_degree": len(adjacency[atom_id]),
            }
        )

    node_fields = [
        "atom_id",
        "element",
        "canonical_role",
        "in_current_v3",
        "retained_in_proposed_v4",
        "mandatory_restoration",
        "exposed_next_shell",
        "graph_distance_from_P1641",
        "canonical_degree",
    ]

    with OUTPUT_NODES.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=node_fields,
        )
        writer.writeheader()
        writer.writerows(node_records)

    edge_fields = [
        "included_atom",
        "included_element",
        "outside_atom",
        "outside_element",
        "outside_role",
        "edge_id",
        "heavy_heavy_cut",
    ]

    with OUTPUT_EDGES.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=edge_fields,
        )
        writer.writeheader()
        writer.writerows(exposed_edges)

    report = {
        "decision": (
            "QM_F06_UPPER_V4_FIRST_SHELL_DEFINED_"
            "OUTER_BOUNDARY_CLOSURE_PENDING"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "central_atom": CENTRAL_ATOM,
        "defective_caps_to_remove": sorted(
            DEFECTIVE_CAPS
        ),
        "mandatory_real_atoms_to_restore": sorted(
            MANDATORY_REAL_RESTORATION
        ),
        "expected_central_neighbors": sorted(
            EXPECTED_CENTRAL_NEIGHBORS
        ),
        "canonical_central_neighbors": sorted(
            central_neighbors
        ),
        "central_neighbor_gate_pass": (
            central_neighbor_gate
        ),
        "current_v3_atom_count": len(current_v3_ids),
        "proposed_v4_first_shell_atom_count": len(
            proposed_v4_ids
        ),
        "newly_exposed_heavy_atoms": exposed_heavy_atoms,
        "exposed_edge_count": len(exposed_edges),
        "exposed_heavy_edge_count": sum(
            record["heavy_heavy_cut"]
            for record in exposed_edges
        ),
        "first_shell": sorted(first_shell),
        "files": {
            "candidate_nodes_csv": str(
                OUTPUT_NODES.relative_to(ROOT)
            ),
            "candidate_edges_csv": str(
                OUTPUT_EDGES.relative_to(ROOT)
            ),
        },
        "files_sha256": {
            "canonical_nodes": sha256(GRAPH_NODES),
            "canonical_edges": sha256(GRAPH_EDGES),
            "v3_map": sha256(V3_MAP),
            "candidate_nodes_csv": sha256(OUTPUT_NODES),
            "candidate_edges_csv": sha256(OUTPUT_EDGES),
        },
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

    OUTPUT_JSON.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 78)
    print("QM_F06 UPPER V4 EXPANSION-SCOPE AUDIT")
    print("=" * 78)
    print(
        "Central-neighbor gate:",
        "PASS" if central_neighbor_gate else "FAIL",
    )
    print(
        "Defective caps removed:",
        sorted(DEFECTIVE_CAPS),
    )
    print(
        "Mandatory restoration:",
        sorted(MANDATORY_REAL_RESTORATION),
    )
    print(
        "Newly exposed heavy atoms:",
        exposed_heavy_atoms,
    )
    print(
        "Exposed heavy-heavy cuts:",
        report["exposed_heavy_edge_count"],
    )
    print()
    print("Decision:", report["decision"])
    print("Report:", OUTPUT_JSON)
    print("Nodes:", OUTPUT_NODES)
    print("Edges:", OUTPUT_EDGES)
    print("V4 construction authorized: False")


if __name__ == "__main__":
    main()
