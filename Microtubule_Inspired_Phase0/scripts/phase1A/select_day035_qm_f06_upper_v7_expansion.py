#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict, deque
from itertools import combinations
import csv
import json

ROOT = Path(__file__).resolve().parents[2]

CANONICAL_EDGES = ROOT / (
    "runs/phase1A/"
    "day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_edges.csv"
)

CANONICAL_NODES = ROOT / (
    "runs/phase1A/"
    "day024_chemical_end_rim_design/"
    "16_r2_selected_full_density_longer_bn_bridge_graph/"
    "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

V6B_MAP = ROOT / (
    "runs/phase1A/"
    "day033_qm_f06_upper_v6b_orca_input/"
    "QM_F06_UPPER_V6B_constraint_map.csv"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/"
    "day035_qm_f06_upper_v7_expansion_selection"
)

OUTPUT_RANKING = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7_expansion_ranking.csv"
)

OUTPUT_CUTS = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7_selected_boundary_cuts.csv"
)

OUTPUT_REPORT = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7_EXPANSION_SELECTION.json"
)

REQUIRED_ATOM = "A:UPPER:13:1"
REMOVED_CAP = "HCAPV2:UPPER:03"

# Search to graph radius 2 around the missing canonical atom.
SEARCH_RADIUS = 2

# Keep the exhaustive search manageable and chemically local.
MAX_OPTIONAL_ATOMS = 12


def load_csv(path):
    with path.open(
        encoding="utf-8",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle))


def subsets(values):
    values = sorted(values)

    for size in range(len(values) + 1):
        for subset in combinations(values, size):
            yield set(subset)


def main():
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    edge_rows = load_csv(CANONICAL_EDGES)
    node_rows = load_csv(CANONICAL_NODES)
    map_rows = load_csv(V6B_MAP)

    element = {
        row["node_id"]: row["element"]
        for row in node_rows
    }

    node_type = {
        row["node_id"]: row["node_type"]
        for row in node_rows
    }

    heavy_graph = defaultdict(set)
    heavy_edges = set()

    for row in edge_rows:
        first = row["source_node"]
        second = row["target_node"]

        if (
            row["source_element"] == "H"
            or row["target_element"] == "H"
        ):
            continue

        pair = tuple(sorted((first, second)))
        heavy_edges.add(pair)

        heavy_graph[first].add(second)
        heavy_graph[second].add(first)

    retained_real = {
        row["atom_id"]
        for row in map_rows
        if (
            row["element"] != "H"
            and row["atom_id"] != REMOVED_CAP
        )
    }

    if REQUIRED_ATOM in retained_real:
        raise RuntimeError(
            f"{REQUIRED_ATOM} is already retained."
        )

    # Breadth-first local neighborhood.
    distance = {REQUIRED_ATOM: 0}
    queue = deque([REQUIRED_ATOM])

    while queue:
        center = queue.popleft()

        if distance[center] >= SEARCH_RADIUS:
            continue

        for neighbor in heavy_graph[center]:
            if neighbor not in distance:
                distance[neighbor] = (
                    distance[center] + 1
                )
                queue.append(neighbor)

    local_candidates = {
        atom_id
        for atom_id, radius in distance.items()
        if (
            radius <= SEARCH_RADIUS
            and atom_id not in retained_real
            and element.get(atom_id) in {"B", "N"}
        )
    }

    local_candidates.add(REQUIRED_ATOM)

    optional = (
        local_candidates
        - {REQUIRED_ATOM}
    )

    if len(optional) > MAX_OPTIONAL_ATOMS:
        raise RuntimeError(
            "Local candidate space is unexpectedly large: "
            f"{len(optional)} optional atoms."
        )

    records = []
    cut_records_by_subset = {}

    for optional_subset in subsets(optional):
        added = {
            REQUIRED_ATOM,
            *optional_subset,
        }

        inside = retained_real | added

        boundary_cuts = []
        external_count = defaultdict(int)

        for first in inside:
            for second in heavy_graph[first]:
                if second in inside:
                    continue

                boundary_cuts.append(
                    tuple(sorted((first, second)))
                )
                external_count[first] += 1

        boundary_cuts = sorted(
            set(boundary_cuts)
        )

        multi_cut_centers = sorted(
            atom_id
            for atom_id, count
            in external_count.items()
            if count > 1
        )

        maximum_external_degree = max(
            external_count.values(),
            default=0,
        )

        required_internal = (
            tuple(sorted((
                REQUIRED_ATOM,
                "A:UPPER:14:2",
            )))
            in heavy_edges
            and REQUIRED_ATOM in inside
            and "A:UPPER:14:2" in inside
        )

        required_atom_external_degree = (
            external_count.get(
                REQUIRED_ATOM,
                0,
            )
        )

        # Primary scientific gate:
        # the newly restored atom may have at most one cut.
        required_center_gate = (
            required_atom_external_degree <= 1
        )

        # Prefer no new center with more than one cut.
        no_multicut_gate = (
            len(multi_cut_centers) == 0
        )

        # Score lexicographically encoded as an integer.
        score = (
            1_000_000
            * int(not required_internal)
            + 100_000
            * int(not required_center_gate)
            + 10_000
            * len(multi_cut_centers)
            + 1_000
            * maximum_external_degree
            + 50
            * len(boundary_cuts)
            + len(added)
        )

        subset_id = "+".join(
            sorted(added)
        )

        record = {
            "subset_id": subset_id,
            "added_atom_count": len(added),
            "added_atoms": "|".join(
                sorted(added)
            ),
            "total_real_atom_count": (
                len(inside)
            ),
            "boundary_cut_count": (
                len(boundary_cuts)
            ),
            "multi_cut_boundary_center_count": (
                len(multi_cut_centers)
            ),
            "maximum_external_heavy_degree": (
                maximum_external_degree
            ),
            "multi_cut_boundary_centers": "|".join(
                multi_cut_centers
            ),
            "required_A13_1_A14_2_internal": (
                required_internal
            ),
            "A13_1_external_heavy_degree": (
                required_atom_external_degree
            ),
            "A13_1_single_cut_gate": (
                required_center_gate
            ),
            "no_multi_cut_center_gate": (
                no_multicut_gate
            ),
            "score": score,
        }

        records.append(record)
        cut_records_by_subset[
            subset_id
        ] = boundary_cuts

    records.sort(
        key=lambda row: (
            int(row["score"]),
            int(row["added_atom_count"]),
            row["subset_id"],
        )
    )

    for rank, row in enumerate(
        records,
        start=1,
    ):
        row["rank"] = rank

    selected = records[0]

    fieldnames = [
        "rank",
        "subset_id",
        "added_atom_count",
        "added_atoms",
        "total_real_atom_count",
        "boundary_cut_count",
        "multi_cut_boundary_center_count",
        "maximum_external_heavy_degree",
        "multi_cut_boundary_centers",
        "required_A13_1_A14_2_internal",
        "A13_1_external_heavy_degree",
        "A13_1_single_cut_gate",
        "no_multi_cut_center_gate",
        "score",
    ]

    with OUTPUT_RANKING.open(
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

    selected_cuts = (
        cut_records_by_subset[
            selected["subset_id"]
        ]
    )

    with OUTPUT_CUTS.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "inside_atom",
                "inside_element",
                "outside_atom",
                "outside_element",
            ],
        )
        writer.writeheader()

        for first, second in selected_cuts:
            if first in retained_real or (
                first in
                set(
                    selected[
                        "added_atoms"
                    ].split("|")
                )
            ):
                inside_atom = first
                outside_atom = second
            else:
                inside_atom = second
                outside_atom = first

            writer.writerow({
                "inside_atom": inside_atom,
                "inside_element": element[
                    inside_atom
                ],
                "outside_atom": outside_atom,
                "outside_element": element[
                    outside_atom
                ],
            })

    selected_pass = (
        selected[
            "required_A13_1_A14_2_internal"
        ]
        and selected[
            "A13_1_single_cut_gate"
        ]
    )

    report = {
        "model": "QM_F06_UPPER_V7",
        "decision": (
            "QM_F06_UPPER_V7_EXPANSION_"
            "SUBSET_SELECTED_CONSTRUCTION_REQUIRED"
            if selected_pass
            else
            "QM_F06_UPPER_V7_NO_ACCEPTABLE_"
            "LOCAL_EXPANSION_SUBSET"
        ),
        "required_atom": REQUIRED_ATOM,
        "removed_cap": REMOVED_CAP,
        "search_radius": SEARCH_RADIUS,
        "retained_real_atom_count": (
            len(retained_real)
        ),
        "local_candidate_atoms": sorted(
            local_candidates
        ),
        "optional_atom_count": len(optional),
        "subset_count": len(records),
        "selected_candidate": selected,
        "construction_authorized": (
            selected_pass
        ),
        "orca_authorized": False,
        "RESP_authorized": False,
        "force_field_adoption_authorized": False,
        "MD_authorized": False,
    }

    OUTPUT_REPORT.write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=" * 108)
    print("QM_F06 UPPER V7 EXPANSION SUBSET SELECTION")
    print("=" * 108)
    print(
        "Retained V6-B real atoms:",
        len(retained_real),
    )
    print(
        "Local candidate atoms:",
        len(local_candidates),
    )
    print(
        "Optional atoms:",
        len(optional),
    )
    print(
        "Subsets evaluated:",
        len(records),
    )

    print()
    print("Local candidate inventory:")

    for atom_id in sorted(
        local_candidates
    ):
        print(
            f"  {atom_id:28s} "
            f"{element.get(atom_id, '?'):2s} "
            f"radius={distance.get(atom_id)} "
            f"type={node_type.get(atom_id, '')}"
        )

    print()
    print("Top 15 subsets:")

    for row in records[:15]:
        print(
            f"rank={row['rank']:3d} | "
            f"added={row['added_atom_count']:2d} | "
            f"cuts={row['boundary_cut_count']:2d} | "
            f"multi={row['multi_cut_boundary_center_count']:2d} | "
            f"A13:1 external="
            f"{row['A13_1_external_heavy_degree']} | "
            f"score={row['score']:7d} | "
            f"{row['added_atoms']}"
        )

    print()
    print("Selected candidate:")

    for key, value in selected.items():
        print(f"  {key}: {value}")

    print()
    print("Decision:", report["decision"])
    print("Ranking:", OUTPUT_RANKING)
    print("Selected cuts:", OUTPUT_CUTS)
    print("Report:", OUTPUT_REPORT)
    print()
    print(
        "V7 construction authorized:",
        selected_pass,
    )
    print("ORCA authorized: False")
    print("RESP authorized: False")


if __name__ == "__main__":
    main()
