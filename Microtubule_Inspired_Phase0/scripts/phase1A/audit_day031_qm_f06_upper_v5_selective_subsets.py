#!/usr/bin/env python3
"""
Exhaustively evaluate selective V5 expansions.

All subsets of the seven first/second-shell candidates are tested
against the canonical heavy-atom graph.

Ranking priorities:
1. close the two V4 cuts implicated in the ambiguous S:1738 caps;
2. minimize total heavy-heavy boundary cuts;
3. minimize newly created cuts relative to V4;
4. minimize added real atoms;
5. minimize boundary centers with more than one external heavy neighbor.

This script selects no model and authorizes no construction.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

CANONICAL_EDGES = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_edges.csv"
)

CANONICAL_NODES = ROOT / (
    "runs/phase1A/day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

V4_MAP = ROOT / (
    "runs/phase1A/day030_qm_f06_upper_v4_construction/"
    "QM_F06_UPPER_V4_atom_role_provenance_map.csv"
)

SCOPE_REPORT = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5_expansion_scope/"
    "QM_F06_UPPER_V5_EXPANSION_SCOPE.json"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day031_qm_f06_upper_v5_selective_subsets"
)

REPORT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5_SELECTIVE_SUBSET_AUDIT.json"
)

SUBSET_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5_selective_subset_ranking.csv"
)

CUT_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V5_top_subset_boundary_cuts.csv"
)

TARGET_CUTS = {
    tuple(sorted((
        "S:1738",
        "BR4:UPPER:14:1",
    ))),
    tuple(sorted((
        "S:1738",
        "P:1637",
    ))),
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
        return reader.fieldnames or [], list(reader)


def find_column(headers, candidates):
    lookup = {
        header.lower(): header
        for header in headers
    }

    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]

    raise RuntimeError(
        f"Could not resolve column from {candidates}; "
        f"headers={headers}"
    )


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {
        "true",
        "1",
        "yes",
        "pass",
    }


def pair(first: str, second: str):
    return tuple(sorted((first, second)))


def boundary_cuts(model, adjacency):
    cuts = set()

    for atom_id in model:
        for neighbor in adjacency.get(atom_id, set()):
            if neighbor not in model:
                cuts.add(pair(atom_id, neighbor))

    return cuts


def main() -> None:
    for path in (
        CANONICAL_EDGES,
        CANONICAL_NODES,
        V4_MAP,
        SCOPE_REPORT,
    ):
        require_file(path)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    edge_headers, edge_rows = read_csv(
        CANONICAL_EDGES
    )

    node_headers, node_rows = read_csv(
        CANONICAL_NODES
    )

    _, v4_rows = read_csv(V4_MAP)

    scope = json.loads(
        SCOPE_REPORT.read_text(
            encoding="utf-8"
        )
    )

    first_key = find_column(
        edge_headers,
        (
            "node_a",
            "atom_id_1",
            "source",
            "source_node",
            "atom1",
            "first_atom",
        ),
    )

    second_key = find_column(
        edge_headers,
        (
            "node_b",
            "atom_id_2",
            "target",
            "target_node",
            "atom2",
            "second_atom",
        ),
    )

    node_id_key = find_column(
        node_headers,
        ("node_id", "atom_id", "id"),
    )

    node_element_key = find_column(
        node_headers,
        ("element", "atom_element"),
    )

    node_elements = {
        row[node_id_key]: row[node_element_key]
        for row in node_rows
    }

    adjacency = defaultdict(set)
    edge_ids = {}

    edge_id_key = next(
        (
            key
            for key in ("edge_id", "id", "bond_id")
            if key in edge_headers
        ),
        None,
    )

    for row in edge_rows:
        first = row[first_key]
        second = row[second_key]

        if not first or not second:
            continue

        if (
            node_elements.get(first) == "H"
            or node_elements.get(second) == "H"
        ):
            continue

        adjacency[first].add(second)
        adjacency[second].add(first)

        edge_ids[pair(first, second)] = (
            row[edge_id_key]
            if edge_id_key
            else ""
        )

    v4_real = {
        row["atom_id"]
        for row in v4_rows
        if not parse_bool(row["artificial_cap"])
    }

    candidates = sorted(
        set(scope["first_shell_new_atoms"])
        | set(scope["second_shell_new_atoms"])
    )

    if len(candidates) != 7:
        raise RuntimeError(
            "Expected exactly seven V5 candidates; "
            f"found {candidates}"
        )

    baseline_cuts = boundary_cuts(
        v4_real,
        adjacency,
    )

    records = []
    subset_cut_map = {}

    for subset_size in range(len(candidates) + 1):
        for subset_tuple in itertools.combinations(
            candidates,
            subset_size,
        ):
            subset = set(subset_tuple)
            model = set(v4_real) | subset

            cuts = boundary_cuts(
                model,
                adjacency,
            )

            closed = baseline_cuts - cuts
            new = cuts - baseline_cuts

            target_closed = (
                TARGET_CUTS
                & closed
            )

            target_remaining = (
                TARGET_CUTS
                & cuts
            )

            external_degree = Counter()

            for first, second in cuts:
                if first in model:
                    external_degree[first] += 1
                if second in model:
                    external_degree[second] += 1

            multi_cut_centers = sorted(
                atom_id
                for atom_id, count in external_degree.items()
                if count > 1
            )

            maximum_external_degree = (
                max(external_degree.values())
                if external_degree
                else 0
            )

            # Lower score is better. The target-cut penalty
            # dominates all secondary criteria.
            score = (
                1000 * len(target_remaining)
                + 100 * len(cuts)
                + 20 * len(new)
                + 5 * len(multi_cut_centers)
                + len(subset)
            )

            subset_id = (
                "BASELINE"
                if not subset
                else "+".join(sorted(subset))
            )

            record = {
                "rank": None,
                "subset_id": subset_id,
                "added_atom_count": len(subset),
                "added_atoms": "|".join(
                    sorted(subset)
                ),
                "total_real_atoms": len(model),
                "boundary_cut_count": len(cuts),
                "cuts_closed_vs_v4": len(closed),
                "new_cuts_vs_v4": len(new),
                "target_cuts_closed": len(
                    target_closed
                ),
                "target_cuts_remaining": len(
                    target_remaining
                ),
                "multi_cut_boundary_center_count": len(
                    multi_cut_centers
                ),
                "maximum_external_heavy_degree": (
                    maximum_external_degree
                ),
                "multi_cut_boundary_centers": "|".join(
                    multi_cut_centers
                ),
                "score": score,
            }

            records.append(record)
            subset_cut_map[subset_id] = cuts

    records.sort(
        key=lambda row: (
            row["target_cuts_remaining"],
            row["boundary_cut_count"],
            row["new_cuts_vs_v4"],
            row["multi_cut_boundary_center_count"],
            row["added_atom_count"],
            row["added_atoms"],
        )
    )

    for rank, row in enumerate(records, start=1):
        row["rank"] = rank

    with SUBSET_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(records[0]),
        )
        writer.writeheader()
        writer.writerows(records)

    top_records = records[:10]
    top_ids = {
        row["subset_id"]
        for row in top_records
    }

    cut_records = []

    for subset_id in top_ids:
        model = set(v4_real)

        if subset_id != "BASELINE":
            model.update(subset_id.split("+"))

        cuts = subset_cut_map[subset_id]

        for first, second in sorted(cuts):
            inside = (
                first if first in model else second
            )

            outside = (
                second if first in model else first
            )

            cut_records.append({
                "subset_id": subset_id,
                "edge_id": edge_ids.get(
                    pair(first, second),
                    "",
                ),
                "inside_atom": inside,
                "inside_element": node_elements.get(
                    inside,
                    "",
                ),
                "outside_atom": outside,
                "outside_element": node_elements.get(
                    outside,
                    "",
                ),
                "preexisting_v4_cut": (
                    pair(first, second)
                    in baseline_cuts
                ),
                "target_cut": (
                    pair(first, second)
                    in TARGET_CUTS
                ),
            })

    with CUT_CSV.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(cut_records[0]),
        )
        writer.writeheader()
        writer.writerows(cut_records)

    best = records[0]

    report = {
        "decision": (
            "QM_F06_UPPER_V5_SELECTIVE_SUBSETS_AUDITED_"
            "CHEMICAL_MODEL_SELECTION_PENDING"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "candidate_atoms": candidates,
        "subset_count": len(records),
        "v4_baseline": {
            "real_atom_count": len(v4_real),
            "boundary_cut_count": len(
                baseline_cuts
            ),
        },
        "target_cuts": [
            list(value)
            for value in sorted(TARGET_CUTS)
        ],
        "best_ranked_subset": best,
        "top_ten": top_records,
        "files": {
            "ranking": str(
                SUBSET_CSV.relative_to(ROOT)
            ),
            "top_subset_cuts": str(
                CUT_CSV.relative_to(ROOT)
            ),
        },
        "files_sha256": {
            "canonical_edges": sha256(
                CANONICAL_EDGES
            ),
            "canonical_nodes": sha256(
                CANONICAL_NODES
            ),
            "v4_map": sha256(V4_MAP),
            "scope_report": sha256(
                SCOPE_REPORT
            ),
            "ranking": sha256(SUBSET_CSV),
            "top_subset_cuts": sha256(
                CUT_CSV
            ),
        },
        "authorization": {
            "v5_model_selected": False,
            "v5_geometry_construction_authorized": False,
            "orca_input_preparation_authorized": False,
            "orca_execution_authorized": False,
            "RESP_execution_authorized": False,
            "force_field_adoption_authorized": False,
            "MD_authorized": False,
        },
    }

    REPORT.write_text(
        json.dumps(
            report,
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    print("=" * 110)
    print("QM_F06 UPPER V5 SELECTIVE-SUBSET AUDIT")
    print("=" * 110)

    print("Candidate atoms:", candidates)
    print("Subsets evaluated:", len(records))
    print(
        "V4 baseline cuts:",
        len(baseline_cuts),
    )
    print()

    print("TOP 15 SUBSETS")

    for row in records[:15]:
        print(
            f"rank={row['rank']:3d} "
            f"added={row['added_atom_count']:1d} "
            f"cuts={row['boundary_cut_count']:2d} "
            f"closed={row['cuts_closed_vs_v4']:2d} "
            f"new={row['new_cuts_vs_v4']:2d} "
            f"target_remaining="
            f"{row['target_cuts_remaining']:1d} "
            f"multi_centers="
            f"{row['multi_cut_boundary_center_count']:2d} "
            f"{row['subset_id']}"
        )

    print()
    print(
        "Best ranked subset:",
        best["subset_id"],
    )
    print(
        "Best boundary cuts:",
        best["boundary_cut_count"],
    )
    print(
        "Best target cuts remaining:",
        best["target_cuts_remaining"],
    )
    print(
        "Best new cuts versus V4:",
        best["new_cuts_vs_v4"],
    )
    print()

    print(
        "Decision:",
        report["decision"],
    )
    print("Report:", REPORT)
    print("Ranking:", SUBSET_CSV)
    print("Top cuts:", CUT_CSV)
    print()
    print("V5 model selected: False")
    print("V5 construction authorized: False")
    print("ORCA execution authorized: False")


if __name__ == "__main__":
    main()
