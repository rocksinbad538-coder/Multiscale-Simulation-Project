#!/usr/bin/env python3
"""
Compute the incremental heavy-atom boundary cuts introduced by the
mandatory QM_F06 UPPER V4 restoration.

The analysis separates:
1. heavy-heavy cuts already present in retained V3;
2. heavy-heavy cuts in the proposed V4 first-shell model;
3. genuinely new cuts caused by restoring the real P:1641 environment.

No geometry is constructed and no QM execution is authorized.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
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
    "runs/phase1A/day030_qm_f06_upper_v4_boundary_delta"
)

OUTPUT_JSON = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_BOUNDARY_DELTA.json"
)

OUTPUT_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V4_boundary_cut_delta.csv"
)

DEFECTIVE_CAPS = {
    "HCAP:UPPER:01",
    "HCAP:UPPER:04",
}

MANDATORY_RESTORATION = {
    "P:1640",
    "P:1581",
    "P:1583",
    "S:1739",
    "P:1639",
    "H4:UPPER:0203:0",
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


def canonical_pair(a: str, b: str) -> tuple[str, str]:
    return tuple(sorted((a, b)))


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    _, node_rows = read_csv(GRAPH_NODES)
    _, edge_rows = read_csv(GRAPH_EDGES)
    _, v3_rows = read_csv(V3_MAP)

    nodes = {
        row["node_id"]: {
            "element": row["element"],
            "role": row["node_type"],
        }
        for row in node_rows
    }

    adjacency: dict[str, set[str]] = defaultdict(set)
    edge_metadata = {}

    for row in edge_rows:
        a = row["source_node"]
        b = row["target_node"]
        pair = canonical_pair(a, b)

        adjacency[a].add(b)
        adjacency[b].add(a)
        edge_metadata[pair] = row

    current_v3_ids = {
        row["atom_id"]
        for row in v3_rows
    }

    retained_v3_ids = (
        current_v3_ids - DEFECTIVE_CAPS
    )

    proposed_v4_ids = (
        retained_v3_ids | MANDATORY_RESTORATION
    )

    def heavy_boundary_cuts(
        included_ids: set[str],
    ) -> dict[tuple[str, str], dict]:
        records = {}

        for included_atom in sorted(included_ids):
            if included_atom not in nodes:
                continue

            if nodes[included_atom]["element"] == "H":
                continue

            for outside_atom in sorted(
                adjacency[included_atom]
            ):
                if outside_atom in included_ids:
                    continue

                if nodes[outside_atom]["element"] == "H":
                    continue

                pair = canonical_pair(
                    included_atom,
                    outside_atom,
                )

                edge = edge_metadata[pair]

                records[pair] = {
                    "included_atom": included_atom,
                    "included_element": nodes[
                        included_atom
                    ]["element"],
                    "outside_atom": outside_atom,
                    "outside_element": nodes[
                        outside_atom
                    ]["element"],
                    "outside_role": nodes[
                        outside_atom
                    ]["role"],
                    "edge_id": edge["edge_id"],
                    "edge_type": edge["edge_type"],
                }

        return records

    baseline_cuts = heavy_boundary_cuts(
        retained_v3_ids
    )

    proposed_cuts = heavy_boundary_cuts(
        proposed_v4_ids
    )

    baseline_pairs = set(baseline_cuts)
    proposed_pairs = set(proposed_cuts)

    new_pairs = proposed_pairs - baseline_pairs
    removed_pairs = baseline_pairs - proposed_pairs
    persistent_pairs = baseline_pairs & proposed_pairs

    records = []

    for pair in sorted(
        baseline_pairs | proposed_pairs
    ):
        if pair in new_pairs:
            classification = "NEW_V4_CUT"
            record = proposed_cuts[pair]
        elif pair in removed_pairs:
            classification = "CLOSED_BY_V4"
            record = baseline_cuts[pair]
        else:
            classification = "PREEXISTING_V3_CUT"
            record = proposed_cuts[pair]

        records.append({
            **record,
            "classification": classification,
            "in_baseline_v3": pair in baseline_pairs,
            "in_proposed_v4": pair in proposed_pairs,
            "touches_mandatory_restoration": bool(
                set(pair) & MANDATORY_RESTORATION
            ),
        })

    fieldnames = [
        "included_atom",
        "included_element",
        "outside_atom",
        "outside_element",
        "outside_role",
        "edge_id",
        "edge_type",
        "classification",
        "in_baseline_v3",
        "in_proposed_v4",
        "touches_mandatory_restoration",
    ]

    with OUTPUT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(records)

    new_cut_records = [
        proposed_cuts[pair]
        for pair in sorted(new_pairs)
    ]

    new_outside_atoms = sorted({
        record["outside_atom"]
        for record in new_cut_records
    })

    closed_cut_records = [
        baseline_cuts[pair]
        for pair in sorted(removed_pairs)
    ]

    report = {
        "decision": (
            "QM_F06_UPPER_V4_BOUNDARY_DELTA_DEFINED_"
            "CLOSURE_CLASSIFICATION_PENDING"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "defective_caps_removed": sorted(
            DEFECTIVE_CAPS
        ),
        "mandatory_restoration": sorted(
            MANDATORY_RESTORATION
        ),
        "counts": {
            "current_v3_atoms": len(
                current_v3_ids
            ),
            "retained_v3_atoms": len(
                retained_v3_ids
            ),
            "proposed_v4_atoms": len(
                proposed_v4_ids
            ),
            "baseline_v3_heavy_cuts": len(
                baseline_pairs
            ),
            "proposed_v4_heavy_cuts": len(
                proposed_pairs
            ),
            "preexisting_v3_cuts": len(
                persistent_pairs
            ),
            "cuts_closed_by_v4": len(
                removed_pairs
            ),
            "genuinely_new_v4_cuts": len(
                new_pairs
            ),
        },
        "new_v4_cuts": new_cut_records,
        "newly_exposed_heavy_atoms": (
            new_outside_atoms
        ),
        "cuts_closed_by_v4": (
            closed_cut_records
        ),
        "files": {
            "delta_csv": str(
                OUTPUT_CSV.relative_to(ROOT)
            ),
        },
        "files_sha256": {
            "canonical_nodes": sha256(
                GRAPH_NODES
            ),
            "canonical_edges": sha256(
                GRAPH_EDGES
            ),
            "v3_map": sha256(V3_MAP),
            "delta_csv": sha256(OUTPUT_CSV),
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
    print("QM_F06 UPPER V4 BOUNDARY-DELTA AUDIT")
    print("=" * 78)
    print(
        "Baseline V3 heavy cuts:",
        len(baseline_pairs),
    )
    print(
        "Proposed V4 heavy cuts:",
        len(proposed_pairs),
    )
    print(
        "Pre-existing V3 cuts:",
        len(persistent_pairs),
    )
    print(
        "Cuts closed by V4:",
        len(removed_pairs),
    )
    print(
        "Genuinely new V4 cuts:",
        len(new_pairs),
    )
    print(
        "Newly exposed heavy atoms:",
        new_outside_atoms,
    )

    print()
    print("NEW V4 HEAVY-HEAVY CUTS")

    for record in new_cut_records:
        print(
            f"{record['included_atom']:24s} -- "
            f"{record['outside_atom']:24s} "
            f"[{record['edge_id']}]"
        )

    print()
    print("CUTS CLOSED BY V4")

    for record in closed_cut_records:
        print(
            f"{record['included_atom']:24s} -- "
            f"{record['outside_atom']:24s} "
            f"[{record['edge_id']}]"
        )

    print()
    print("Decision:", report["decision"])
    print("Report:", OUTPUT_JSON)
    print("CSV:", OUTPUT_CSV)
    print("V4 construction authorized: False")


if __name__ == "__main__":
    main()
