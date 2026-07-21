#!/usr/bin/env python3
"""
Evaluate final outer-boundary closure alternatives for QM_F06 UPPER V4.

The audit determines:
- heavy coordination retained for each second-shell boundary atom;
- number of cut heavy bonds;
- number of artificial H caps required for direct closure;
- whether the resulting B coordination is nominally trivalent;
- outside atoms shared by multiple boundary centers;
- whether a selective additional expansion merits consideration.

No geometry is constructed and no QM execution is authorized.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
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
    "day030_qm_f06_upper_v4_final_closure"
)

REPORT_PATH = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_FINAL_CLOSURE_AUDIT.json"
)

BOUNDARY_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_boundary_centers.csv"
)

OUTSIDE_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_outside_atom_leverage.csv"
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

EXPECTED_BOUNDARY_CENTERS = SECOND_SHELL


def read_csv(path: Path) -> list[dict[str, str]]:
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
            "raw": row,
        }
        for row in node_rows
    }

    adjacency: dict[str, set[str]] = defaultdict(set)

    for row in edge_rows:
        a = row["source_node"]
        b = row["target_node"]

        adjacency[a].add(b)
        adjacency[b].add(a)

    current_v3 = {
        row["atom_id"]
        for row in v3_rows
    }

    retained_v3 = current_v3 - DEFECTIVE_CAPS

    proposed_v4 = (
        retained_v3
        | FIRST_SHELL
        | SECOND_SHELL
    )

    missing = (
        FIRST_SHELL | SECOND_SHELL
    ) - set(nodes)

    if missing:
        raise RuntimeError(
            "Required canonical atoms are missing: "
            f"{sorted(missing)}"
        )

    boundary_records = []
    outside_to_centers: dict[str, set[str]] = defaultdict(set)

    for center in sorted(EXPECTED_BOUNDARY_CENTERS):
        element = nodes[center]["element"]

        heavy_neighbors = sorted(
            neighbor
            for neighbor in adjacency[center]
            if nodes[neighbor]["element"] != "H"
        )

        real_h_neighbors = sorted(
            neighbor
            for neighbor in adjacency[center]
            if nodes[neighbor]["element"] == "H"
        )

        retained_heavy_neighbors = sorted(
            neighbor
            for neighbor in heavy_neighbors
            if neighbor in proposed_v4
        )

        outside_heavy_neighbors = sorted(
            neighbor
            for neighbor in heavy_neighbors
            if neighbor not in proposed_v4
        )

        retained_real_h_neighbors = sorted(
            neighbor
            for neighbor in real_h_neighbors
            if neighbor in proposed_v4
        )

        artificial_h_caps_required = len(
            outside_heavy_neighbors
        )

        resulting_coordination = (
            len(retained_heavy_neighbors)
            + len(retained_real_h_neighbors)
            + artificial_h_caps_required
        )

        for outside_atom in outside_heavy_neighbors:
            outside_to_centers[outside_atom].add(center)

        direct_cap_topology_pass = (
            element == "B"
            and 1 <= len(retained_heavy_neighbors) <= 2
            and artificial_h_caps_required in {1, 2}
            and resulting_coordination == 3
        )

        boundary_records.append({
            "center_atom": center,
            "element": element,
            "canonical_node_type": nodes[
                center
            ]["node_type"],
            "canonical_heavy_degree": len(
                heavy_neighbors
            ),
            "retained_heavy_degree": len(
                retained_heavy_neighbors
            ),
            "retained_heavy_neighbors": "|".join(
                retained_heavy_neighbors
            ),
            "outside_heavy_degree": len(
                outside_heavy_neighbors
            ),
            "outside_heavy_neighbors": "|".join(
                outside_heavy_neighbors
            ),
            "canonical_real_h_neighbors": "|".join(
                real_h_neighbors
            ),
            "retained_real_h_count": len(
                retained_real_h_neighbors
            ),
            "artificial_h_caps_required": (
                artificial_h_caps_required
            ),
            "resulting_nominal_coordination": (
                resulting_coordination
            ),
            "direct_BH_cap_topology_pass": (
                direct_cap_topology_pass
            ),
        })

    boundary_fields = [
        "center_atom",
        "element",
        "canonical_node_type",
        "canonical_heavy_degree",
        "retained_heavy_degree",
        "retained_heavy_neighbors",
        "outside_heavy_degree",
        "outside_heavy_neighbors",
        "canonical_real_h_neighbors",
        "retained_real_h_count",
        "artificial_h_caps_required",
        "resulting_nominal_coordination",
        "direct_BH_cap_topology_pass",
    ]

    with BOUNDARY_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=boundary_fields,
        )
        writer.writeheader()
        writer.writerows(boundary_records)

    outside_records = []

    for outside_atom, centers in sorted(
        outside_to_centers.items()
    ):
        outside_records.append({
            "outside_atom": outside_atom,
            "element": nodes[outside_atom]["element"],
            "node_type": nodes[outside_atom]["node_type"],
            "connected_boundary_center_count": len(centers),
            "connected_boundary_centers": "|".join(
                sorted(centers)
            ),
            "high_leverage_shared_atom": len(centers) >= 2,
            "canonical_heavy_degree": sum(
                nodes[neighbor]["element"] != "H"
                for neighbor in adjacency[outside_atom]
            ),
        })

    outside_fields = [
        "outside_atom",
        "element",
        "node_type",
        "connected_boundary_center_count",
        "connected_boundary_centers",
        "high_leverage_shared_atom",
        "canonical_heavy_degree",
    ]

    with OUTSIDE_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=outside_fields,
        )
        writer.writeheader()
        writer.writerows(outside_records)

    all_boundary_boron = all(
        row["element"] == "B"
        for row in boundary_records
    )

    all_direct_caps_topologically_valid = all(
        row["direct_BH_cap_topology_pass"]
        for row in boundary_records
    )

    total_artificial_h_required = sum(
        row["artificial_h_caps_required"]
        for row in boundary_records
    )

    shared_outside_atoms = [
        row["outside_atom"]
        for row in outside_records
        if row["high_leverage_shared_atom"]
    ]

    cap_count_distribution = dict(Counter(
        row["artificial_h_caps_required"]
        for row in boundary_records
    ))

    if (
        all_boundary_boron
        and all_direct_caps_topologically_valid
    ):
        decision = (
            "QM_F06_UPPER_V4_DIRECT_BH_CLOSURE_"
            "TOPOLOGICALLY_POSSIBLE_"
            "CHEMICAL_PRECEDENT_COMPARISON_PENDING"
        )
    else:
        decision = (
            "QM_F06_UPPER_V4_DIRECT_BH_CLOSURE_FAIL_"
            "ADDITIONAL_HEAVY_EXPANSION_REQUIRED"
        )

    report = {
        "decision": decision,
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "proposed_v4_atom_count_before_new_caps": len(
            proposed_v4
        ),
        "boundary_centers": sorted(
            EXPECTED_BOUNDARY_CENTERS
        ),
        "all_boundary_centers_are_boron": (
            all_boundary_boron
        ),
        "all_direct_caps_topologically_valid": (
            all_direct_caps_topologically_valid
        ),
        "total_artificial_h_caps_required": (
            total_artificial_h_required
        ),
        "artificial_h_count_distribution": (
            cap_count_distribution
        ),
        "shared_high_leverage_outside_atoms": (
            shared_outside_atoms
        ),
        "boundary_records": boundary_records,
        "outside_records": outside_records,
        "authorization": {
            "direct_BH_closure_chemically_accepted": False,
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
    print("QM_F06 UPPER V4 FINAL-CLOSURE AUDIT")
    print("=" * 78)
    print(
        "All boundary centers are B:",
        all_boundary_boron,
    )
    print(
        "All direct B-H closures topologically valid:",
        all_direct_caps_topologically_valid,
    )
    print(
        "Total artificial H caps required:",
        total_artificial_h_required,
    )
    print(
        "Cap-count distribution:",
        cap_count_distribution,
    )
    print(
        "Shared high-leverage outside atoms:",
        shared_outside_atoms,
    )

    print()
    print("BOUNDARY CENTERS")

    for row in boundary_records:
        print(
            f"{row['center_atom']:12s} "
            f"element={row['element']} "
            f"retained_heavy={row['retained_heavy_degree']} "
            f"cut_heavy={row['outside_heavy_degree']} "
            f"H_caps={row['artificial_h_caps_required']} "
            f"coord={row['resulting_nominal_coordination']} "
            f"pass={row['direct_BH_cap_topology_pass']}"
        )

    print()
    print("OUTSIDE ATOM LEVERAGE")

    for row in outside_records:
        print(
            f"{row['outside_atom']:24s} "
            f"centers={row['connected_boundary_center_count']} "
            f"{row['connected_boundary_centers']}"
        )

    print()
    print("Decision:", decision)
    print("Report:", REPORT_PATH)
    print("Boundary CSV:", BOUNDARY_CSV)
    print("Outside CSV:", OUTSIDE_CSV)
    print("V4 construction authorized: False")


if __name__ == "__main__":
    main()
