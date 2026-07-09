#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter, deque
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

GATE3A_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "01_r2_parent_rim_chemical_audit"
)

GATE3C_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "03_r2_graded_heteropolar_bn_collar_connectivity_blueprint"
)

CYCLE_AUDIT_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "04_r2_graded_collar_cycle_topology_audit"
)

OUTPUT_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "05_r2_hexagonal_edge_completion_seed"
)

PARENT_ATOMS_CSV = (
    GATE3A_ROOT
    / "r2_parent_hbn_atoms.csv"
)

PARENT_BONDS_CSV = (
    GATE3A_ROOT
    / "r2_parent_hbn_geometry_derived_bonds.csv"
)

TERMINAL_ATOMS_CSV = (
    GATE3A_ROOT
    / "r2_parent_terminal_rim_atoms.csv"
)

PARENT_SUMMARY_CSV = (
    GATE3A_ROOT
    / "r2_parent_rim_chemical_audit_summary.csv"
)

BLUEPRINT_SUMMARY_CSV = (
    GATE3C_ROOT
    / "r2_graded_collar_connectivity_blueprint_summary.csv"
)

CYCLE_AUDIT_SUMMARY_CSV = (
    CYCLE_AUDIT_ROOT
    / "r2_collar_cycle_topology_audit_summary.csv"
)

STEP_SEARCH_CSV = (
    OUTPUT_ROOT
    / "r2_hexagonal_edge_completion_step_search.csv"
)

ADDED_NODES_CSV = (
    OUTPUT_ROOT
    / "r2_hexagonal_edge_completion_nodes.csv"
)

ADDED_EDGES_CSV = (
    OUTPUT_ROOT
    / "r2_hexagonal_edge_completion_edges.csv"
)

END_SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_hexagonal_edge_completion_end_summary.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_hexagonal_edge_completion_seed_summary.csv"
)

GATES_CSV = (
    OUTPUT_ROOT
    / "r2_hexagonal_edge_completion_seed_gates.csv"
)

SEED_JSON = (
    OUTPUT_ROOT
    / "r2_hexagonal_edge_completion_seed.json"
)

SOURCE_MANIFEST_CSV = (
    OUTPUT_ROOT
    / "r2_hexagonal_edge_completion_source_manifest.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_HEXAGONAL_EDGE_COMPLETION_SEED_DAY024.md"
)

EXPECTED_PARENT_DECISION = (
    "R2_PARENT_RIM_AND_CHEMICAL_CONSTRAINTS_AUDITED"
)

EXPECTED_BLUEPRINT_DECISION = (
    "R2_GRADED_HETEROPOLAR_BN_COLLAR_CONNECTIVITY_BLUEPRINT_VALIDATED"
)

EXPECTED_CYCLE_AUDIT_DECISION = (
    "R2_GRADED_HETEROPOLAR_BN_COLLAR_REQUIRES_HEXAGONAL_GRAPH_REDESIGN"
)

PASS_DECISION = (
    "R2_HEXAGONAL_EDGE_COMPLETION_SEED_VALIDATED"
)

EXPECTED_PARENT_ATOMS = 1680
EXPECTED_PARENT_BONDS = 2460

EXPECTED_TERMINALS_TOTAL = 60
EXPECTED_TERMINALS_PER_END = 30

EXPECTED_INTERIOR_DEGREE3 = 1620

EXPECTED_ADDED_ATOMS_PER_END = 30
EXPECTED_ADDED_ATOMS_TOTAL = 60

EXPECTED_NEW_EDGES_PER_END = 60
EXPECTED_NEW_EDGES_TOTAL = 120

EXPECTED_SHORTEST_PARENT_PATH = 4
EXPECTED_CLOSED_CYCLE_LENGTH = 6

MAX_STEP = 15


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
    )


def require_file(path: Path) -> None:
    if (
        not path.exists()
        or not path.is_file()
        or path.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing or empty required file: {path}"
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def read_csv_rows(
    path: Path,
) -> list[dict[str, str]]:
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


def read_single_csv_row(
    path: Path,
) -> dict[str, str]:
    rows = read_csv_rows(path)

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected one row in {path}; found {len(rows)}"
        )

    return rows[0]


def write_csv(
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


def parse_int(
    row: dict[str, str],
    key: str,
) -> int:
    try:
        return int(float(row[key]))
    except (
        KeyError,
        TypeError,
        ValueError,
    ) as exc:
        raise RuntimeError(
            f"Could not parse integer field {key!r}"
        ) from exc


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


def shortest_path_length(
    adjacency: dict[int, set[int]],
    source: int,
    target: int,
) -> int | None:
    if source == target:
        return 0

    visited = {
        source
    }

    queue: deque[
        tuple[int, int]
    ] = deque(
        [
            (
                source,
                0,
            )
        ]
    )

    while queue:
        node, distance = queue.popleft()

        for neighbor in adjacency[node]:
            if neighbor == target:
                return distance + 1

            if neighbor in visited:
                continue

            visited.add(neighbor)

            queue.append(
                (
                    neighbor,
                    distance + 1,
                )
            )

    return None


def count_simple_paths_exact_length(
    adjacency: dict[int, set[int]],
    source: int,
    target: int,
    length: int,
) -> int:
    count = 0
    visited = {
        source
    }

    def walk(
        node: int,
        depth: int,
    ) -> None:
        nonlocal count

        if depth == length:
            if node == target:
                count += 1
            return

        remaining = (
            length - depth
        )

        if remaining <= 0:
            return

        for neighbor in adjacency[node]:
            if neighbor in visited:
                continue

            visited.add(neighbor)

            walk(
                neighbor,
                depth + 1,
            )

            visited.remove(neighbor)

    walk(
        source,
        0,
    )

    return count


def connected_components(
    adjacency: dict[int, set[int]],
) -> list[set[int]]:
    remaining = set(adjacency)
    components: list[
        set[int]
    ] = []

    while remaining:
        start = min(remaining)

        component: set[int] = set()

        queue: deque[int] = deque(
            [start]
        )

        while queue:
            node = queue.popleft()

            if node in component:
                continue

            component.add(node)

            queue.extend(
                adjacency[node]
                - component
            )

        components.append(component)

        remaining -= component

    return components


def count_four_cycles(
    adjacency: dict[int, set[int]],
) -> int:
    nodes = sorted(adjacency)

    raw_count = 0

    for first_index, first in enumerate(
        nodes
    ):
        first_neighbors = adjacency[first]

        for second in nodes[
            first_index + 1:
        ]:
            common_count = len(
                first_neighbors
                & adjacency[second]
            )

            if common_count >= 2:
                raw_count += (
                    common_count
                    * (
                        common_count - 1
                    )
                    // 2
                )

    if raw_count % 2 != 0:
        raise RuntimeError(
            "Four-cycle common-neighbor count "
            "was not divisible by two."
        )

    return raw_count // 2


def is_bipartite(
    adjacency: dict[int, set[int]],
) -> bool:
    colors: dict[int, int] = {}

    for start in sorted(adjacency):
        if start in colors:
            continue

        colors[start] = 0

        queue: deque[int] = deque(
            [start]
        )

        while queue:
            node = queue.popleft()

            for neighbor in adjacency[node]:
                if neighbor not in colors:
                    colors[neighbor] = (
                        1 - colors[node]
                    )

                    queue.append(neighbor)

                elif (
                    colors[neighbor]
                    == colors[node]
                ):
                    return False

    return True


def pairing_for_step(
    terminal_indices: list[int],
    step: int,
) -> list[tuple[int, int]]:
    count = len(
        terminal_indices
    )

    pairs = []

    for index in range(count):
        first = terminal_indices[
            index
        ]

        second = terminal_indices[
            (
                index + step
            )
            % count
        ]

        pairs.append(
            tuple(
                sorted(
                    (
                        first,
                        second,
                    )
                )
            )
        )

    return pairs


def add_candidate_repair_nodes(
    parent_adjacency: dict[int, set[int]],
    pairs: list[tuple[int, int]],
    first_new_node: int,
) -> tuple[
    dict[int, set[int]],
    list[int],
]:
    adjacency = {
        node: set(neighbors)
        for node, neighbors
        in parent_adjacency.items()
    }

    new_nodes = []

    for pair_index, (
        first,
        second,
    ) in enumerate(pairs):
        new_node = (
            first_new_node
            + pair_index
        )

        new_nodes.append(
            new_node
        )

        adjacency[new_node] = {
            first,
            second,
        }

        adjacency[first].add(
            new_node
        )

        adjacency[second].add(
            new_node
        )

    return adjacency, new_nodes


def evaluate_step(
    end: str,
    terminal_indices: list[int],
    step: int,
    parent_adjacency: dict[int, set[int]],
    parent_four_cycles: int,
) -> dict[str, Any]:
    pairs = pairing_for_step(
        terminal_indices,
        step,
    )

    unique_pairs = set(
        pairs
    )

    incidence = Counter()

    for first, second in pairs:
        incidence[first] += 1
        incidence[second] += 1

    distances = []

    path4_counts = []

    for first, second in pairs:
        distance = shortest_path_length(
            parent_adjacency,
            first,
            second,
        )

        distances.append(
            -1
            if distance is None
            else distance
        )

        if distance == EXPECTED_SHORTEST_PARENT_PATH:
            path4_counts.append(
                count_simple_paths_exact_length(
                    parent_adjacency,
                    first,
                    second,
                    EXPECTED_SHORTEST_PARENT_PATH,
                )
            )
        else:
            path4_counts.append(0)

    candidate_adjacency, new_nodes = (
        add_candidate_repair_nodes(
            parent_adjacency,
            pairs,
            EXPECTED_PARENT_ATOMS,
        )
    )

    candidate_four_cycles = (
        count_four_cycles(
            candidate_adjacency
        )
    )

    all_distances_four = (
        distances
        and all(
            value
            == EXPECTED_SHORTEST_PARENT_PATH
            for value in distances
        )
    )

    all_have_six_cycle = (
        path4_counts
        and all(
            value >= 1
            for value in path4_counts
        )
    )

    pair_uniqueness = (
        len(unique_pairs)
        == EXPECTED_TERMINALS_PER_END
    )

    incidence_valid = (
        len(incidence)
        == EXPECTED_TERMINALS_PER_END
        and min(
            incidence.values()
        )
        == 2
        and max(
            incidence.values()
        )
        == 2
    )

    no_new_four_cycles = (
        candidate_four_cycles
        == parent_four_cycles
    )

    score = (
        1000
        * int(
            all_distances_four
        )
        + 500
        * int(
            all_have_six_cycle
        )
        + 250
        * int(
            pair_uniqueness
        )
        + 250
        * int(
            incidence_valid
        )
        + 500
        * int(
            no_new_four_cycles
        )
        - 100
        * abs(
            candidate_four_cycles
            - parent_four_cycles
        )
        - sum(
            abs(
                value
                - EXPECTED_SHORTEST_PARENT_PATH
            )
            for value in distances
            if value >= 0
        )
    )

    return {
        "end": end,
        "step": step,
        "pair_count": len(
            pairs
        ),
        "unique_pair_count": len(
            unique_pairs
        ),
        "terminal_incidence_minimum": (
            min(
                incidence.values()
            )
            if incidence
            else 0
        ),
        "terminal_incidence_maximum": (
            max(
                incidence.values()
            )
            if incidence
            else 0
        ),
        "shortest_parent_path_minimum": (
            min(distances)
        ),
        "shortest_parent_path_maximum": (
            max(distances)
        ),
        "pairings_with_parent_path_length4": sum(
            value
            == EXPECTED_SHORTEST_PARENT_PATH
            for value in distances
        ),
        "minimum_length4_paths_per_pair": (
            min(path4_counts)
        ),
        "maximum_length4_paths_per_pair": (
            max(path4_counts)
        ),
        "parent_four_cycles": (
            parent_four_cycles
        ),
        "candidate_four_cycles": (
            candidate_four_cycles
        ),
        "new_four_cycles": (
            candidate_four_cycles
            - parent_four_cycles
        ),
        "all_pairings_close_six_member_cycles": (
            all_distances_four
            and all_have_six_cycle
        ),
        "pair_uniqueness_pass": (
            pair_uniqueness
        ),
        "terminal_incidence_pass": (
            incidence_valid
        ),
        "no_new_four_cycles_pass": (
            no_new_four_cycles
        ),
        "candidate_score": (
            score
        ),
        "new_node_count": len(
            new_nodes
        ),
    }


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        PARENT_ATOMS_CSV,
        PARENT_BONDS_CSV,
        TERMINAL_ATOMS_CSV,
        PARENT_SUMMARY_CSV,
        BLUEPRINT_SUMMARY_CSV,
        CYCLE_AUDIT_SUMMARY_CSV,
    ):
        require_file(required)

    parent_atoms = read_csv_rows(
        PARENT_ATOMS_CSV
    )

    parent_bonds = read_csv_rows(
        PARENT_BONDS_CSV
    )

    terminal_rows = read_csv_rows(
        TERMINAL_ATOMS_CSV
    )

    parent_summary = read_single_csv_row(
        PARENT_SUMMARY_CSV
    )

    blueprint_summary = read_single_csv_row(
        BLUEPRINT_SUMMARY_CSV
    )

    cycle_summary = read_single_csv_row(
        CYCLE_AUDIT_SUMMARY_CSV
    )

    if parent_summary.get(
        "decision"
    ) != EXPECTED_PARENT_DECISION:
        raise RuntimeError(
            "Gate 3A parent audit is not accepted."
        )

    if blueprint_summary.get(
        "decision"
    ) != EXPECTED_BLUEPRINT_DECISION:
        raise RuntimeError(
            "Gate 3C blueprint record is not in "
            "the expected accepted combinatorial state."
        )

    if cycle_summary.get(
        "decision"
    ) != EXPECTED_CYCLE_AUDIT_DECISION:
        raise RuntimeError(
            "Gate 3C.1 did not authorize the "
            "hexagonal-graph redesign path."
        )

    if len(parent_atoms) != EXPECTED_PARENT_ATOMS:
        raise RuntimeError(
            "Unexpected parent atom count: "
            f"{len(parent_atoms)}/"
            f"{EXPECTED_PARENT_ATOMS}"
        )

    if len(parent_bonds) != EXPECTED_PARENT_BONDS:
        raise RuntimeError(
            "Unexpected parent bond count: "
            f"{len(parent_bonds)}/"
            f"{EXPECTED_PARENT_BONDS}"
        )

    elements: dict[int, str] = {}

    for row in parent_atoms:
        index = int(
            float(
                row[
                    "hbn_local_index_0based"
                ]
            )
        )

        elements[index] = row[
            "element"
        ]

    if len(elements) != EXPECTED_PARENT_ATOMS:
        raise RuntimeError(
            "Parent atom indices are incomplete or duplicated."
        )

    parent_adjacency: dict[
        int,
        set[int],
    ] = {
        index: set()
        for index in range(
            EXPECTED_PARENT_ATOMS
        )
    }

    parent_edge_set: set[
        tuple[int, int]
    ] = set()

    for row in parent_bonds:
        first = int(
            float(
                row[
                    "atom_i_local_0based"
                ]
            )
        )

        second = int(
            float(
                row[
                    "atom_j_local_0based"
                ]
            )
        )

        pair = tuple(
            sorted(
                (
                    first,
                    second,
                )
            )
        )

        if pair in parent_edge_set:
            raise RuntimeError(
                "Duplicate parent edge detected."
            )

        parent_edge_set.add(
            pair
        )

        parent_adjacency[first].add(
            second
        )

        parent_adjacency[second].add(
            first
        )

    parent_degrees = {
        node: len(neighbors)
        for node, neighbors
        in parent_adjacency.items()
    }

    degree_counts = Counter(
        parent_degrees.values()
    )

    parent_four_cycles = (
        count_four_cycles(
            parent_adjacency
        )
    )

    parent_bipartite = is_bipartite(
        parent_adjacency
    )

    parent_components = connected_components(
        parent_adjacency
    )

    terminals_by_end: dict[
        str,
        list[dict[str, str]],
    ] = {}

    for end in (
        "LOWER",
        "UPPER",
    ):
        rows = [
            row
            for row in terminal_rows
            if row.get("end") == end
        ]

        rows.sort(
            key=lambda row: int(
                float(
                    row[
                        "circumferential_order"
                    ]
                )
            )
        )

        terminals_by_end[end] = rows

    step_search_rows = []

    selected_by_end: dict[
        str,
        dict[str, Any],
    ] = {}

    for end in (
        "LOWER",
        "UPPER",
    ):
        end_rows = terminals_by_end[
            end
        ]

        if len(end_rows) != EXPECTED_TERMINALS_PER_END:
            raise RuntimeError(
                f"{end}: unexpected terminal count "
                f"{len(end_rows)}/"
                f"{EXPECTED_TERMINALS_PER_END}"
            )

        terminal_indices = [
            int(
                float(
                    row[
                        "hbn_local_index_0based"
                    ]
                )
            )
            for row in end_rows
        ]

        terminal_elements = {
            elements[index]
            for index in terminal_indices
        }

        if len(terminal_elements) != 1:
            raise RuntimeError(
                f"{end}: terminal population is not "
                "elementally pure."
            )

        candidate_rows = []

        for step in range(
            1,
            MAX_STEP + 1,
        ):
            candidate = evaluate_step(
                end,
                terminal_indices,
                step,
                parent_adjacency,
                parent_four_cycles,
            )

            candidate_rows.append(
                candidate
            )

            step_search_rows.append(
                candidate
            )

        candidate_rows.sort(
            key=lambda row: (
                -int(
                    row[
                        "candidate_score"
                    ]
                ),
                int(
                    row[
                        "step"
                    ]
                ),
            )
        )

        selected = candidate_rows[0]

        selected[
            "selected"
        ] = True

        selected_by_end[end] = (
            selected
        )

    for row in step_search_rows:
        selected = selected_by_end[
            row["end"]
        ]

        row["selected"] = (
            int(
                row["step"]
            )
            == int(
                selected["step"]
            )
        )

    write_csv(
        STEP_SEARCH_CSV,
        step_search_rows,
    )

    combined_adjacency = {
        node: set(neighbors)
        for node, neighbors
        in parent_adjacency.items()
    }

    added_node_rows = []
    added_edge_rows = []
    end_summary_rows = []

    next_new_node = EXPECTED_PARENT_ATOMS

    for end in (
        "LOWER",
        "UPPER",
    ):
        selected = selected_by_end[
            end
        ]

        step = int(
            selected[
                "step"
            ]
        )

        end_rows = terminals_by_end[
            end
        ]

        terminal_indices = [
            int(
                float(
                    row[
                        "hbn_local_index_0based"
                    ]
                )
            )
            for row in end_rows
        ]

        parent_element = elements[
            terminal_indices[0]
        ]

        added_element = (
            "N"
            if parent_element == "B"
            else "B"
        )

        pairs = pairing_for_step(
            terminal_indices,
            step,
        )

        end_new_nodes = []

        path_lengths = []
        path4_multiplicities = []

        for pair_index, (
            first,
            second,
        ) in enumerate(pairs):
            new_node = (
                next_new_node
                + pair_index
            )

            end_new_nodes.append(
                new_node
            )

            combined_adjacency[
                new_node
            ] = {
                first,
                second,
            }

            combined_adjacency[
                first
            ].add(
                new_node
            )

            combined_adjacency[
                second
            ].add(
                new_node
            )

            path_length = shortest_path_length(
                parent_adjacency,
                first,
                second,
            )

            if path_length is None:
                raise RuntimeError(
                    f"{end}: disconnected parent pairing."
                )

            path_multiplicity = (
                count_simple_paths_exact_length(
                    parent_adjacency,
                    first,
                    second,
                    EXPECTED_SHORTEST_PARENT_PATH,
                )
                if path_length
                == EXPECTED_SHORTEST_PARENT_PATH
                else 0
            )

            path_lengths.append(
                path_length
            )

            path4_multiplicities.append(
                path_multiplicity
            )

            added_node_rows.append(
                {
                    "end": end,
                    "added_node_global_index_0based": (
                        new_node
                    ),
                    "added_node_id": (
                        f"{end}:HEX_EDGE:"
                        f"{pair_index}"
                    ),
                    "circumferential_index": (
                        pair_index
                    ),
                    "element": (
                        added_element
                    ),
                    "parent_element": (
                        parent_element
                    ),
                    "heavy_degree_in_seed": 2,
                    "reserved_future_valence": 1,
                    "target_final_coordination": 3,
                    "coordinates_assigned": False,
                    "formal_charge_assigned": False,
                    "force_field_type_assigned": False,
                }
            )

            for edge_local_index, parent_node in enumerate(
                (
                    first,
                    second,
                ),
                start=1,
            ):
                added_edge_rows.append(
                    {
                        "end": end,
                        "edge_id": (
                            f"{end}:HEX_EDGE:"
                            f"{pair_index}:"
                            f"{edge_local_index}"
                        ),
                        "parent_node_0based": (
                            parent_node
                        ),
                        "added_node_0based": (
                            new_node
                        ),
                        "parent_element": (
                            elements[
                                parent_node
                            ]
                        ),
                        "added_element": (
                            added_element
                        ),
                        "heteropolar_BN_edge": (
                            elements[
                                parent_node
                            ]
                            != added_element
                        ),
                        "parent_pair_shortest_path": (
                            path_length
                        ),
                        "closed_cycle_length": (
                            path_length + 2
                        ),
                        "length4_parent_path_multiplicity": (
                            path_multiplicity
                        ),
                        "formal_bond_order_assigned": False,
                        "coordinates_assigned": False,
                    }
                )

        incidence = Counter()

        for first, second in pairs:
            incidence[first] += 1
            incidence[second] += 1

        end_summary_rows.append(
            {
                "end": end,
                "selected_circumferential_step": (
                    step
                ),
                "parent_terminal_element": (
                    parent_element
                ),
                "added_complementary_element": (
                    added_element
                ),
                "parent_terminal_atoms": len(
                    terminal_indices
                ),
                "added_edge_completion_atoms": len(
                    end_new_nodes
                ),
                "new_parent_to_added_edges": (
                    2
                    * len(
                        end_new_nodes
                    )
                ),
                "terminal_new_edge_incidence_minimum": (
                    min(
                        incidence.values()
                    )
                ),
                "terminal_new_edge_incidence_maximum": (
                    max(
                        incidence.values()
                    )
                ),
                "parent_pair_shortest_path_minimum": (
                    min(
                        path_lengths
                    )
                ),
                "parent_pair_shortest_path_maximum": (
                    max(
                        path_lengths
                    )
                ),
                "closed_cycle_length_minimum": (
                    min(
                        value + 2
                        for value in path_lengths
                    )
                ),
                "closed_cycle_length_maximum": (
                    max(
                        value + 2
                        for value in path_lengths
                    )
                ),
                "length4_parent_paths_per_pair_minimum": (
                    min(
                        path4_multiplicities
                    )
                ),
                "length4_parent_paths_per_pair_maximum": (
                    max(
                        path4_multiplicities
                    )
                ),
            }
        )

        next_new_node += len(
            end_new_nodes
        )

    write_csv(
        ADDED_NODES_CSV,
        added_node_rows,
    )

    write_csv(
        ADDED_EDGES_CSV,
        added_edge_rows,
    )

    write_csv(
        END_SUMMARY_CSV,
        end_summary_rows,
    )

    augmented_four_cycles = (
        count_four_cycles(
            combined_adjacency
        )
    )

    augmented_components = (
        connected_components(
            combined_adjacency
        )
    )

    augmented_bipartite = (
        is_bipartite(
            combined_adjacency
        )
    )

    augmented_degrees = {
        node: len(neighbors)
        for node, neighbors
        in combined_adjacency.items()
    }

    terminal_indices_all = {
        int(
            float(
                row[
                    "hbn_local_index_0based"
                ]
            )
        )
        for row in terminal_rows
    }

    added_indices = {
        int(
            row[
                "added_node_global_index_0based"
            ]
        )
        for row in added_node_rows
    }

    parent_terminal_degree_failures = [
        node
        for node in terminal_indices_all
        if augmented_degrees[
            node
        ]
        != 3
    ]

    parent_nonterminal_degree_failures = [
        node
        for node in range(
            EXPECTED_PARENT_ATOMS
        )
        if node not in terminal_indices_all
        and augmented_degrees[
            node
        ]
        != 3
    ]

    added_degree_failures = [
        node
        for node in added_indices
        if augmented_degrees[
            node
        ]
        != 2
    ]

    nonheteropolar_new_edges = [
        row
        for row in added_edge_rows
        if not bool(
            row[
                "heteropolar_BN_edge"
            ]
        )
    ]

    all_closed_cycle_lengths = [
        int(
            row[
                "closed_cycle_length"
            ]
        )
        for row in added_edge_rows[
            ::2
        ]
    ]

    selected_steps = {
        row["end"]: int(
            row[
                "selected_circumferential_step"
            ]
        )
        for row in end_summary_rows
    }

    lower_summary = next(
        row
        for row in end_summary_rows
        if row[
            "end"
        ]
        == "LOWER"
    )

    upper_summary = next(
        row
        for row in end_summary_rows
        if row[
            "end"
        ]
        == "UPPER"
    )

    gates = {
        "Gate3A_parent_audit_is_accepted": (
            parent_summary.get(
                "decision"
            )
            == EXPECTED_PARENT_DECISION
        ),
        "Gate3C_combinatorial_blueprint_record_is_preserved": (
            blueprint_summary.get(
                "decision"
            )
            == EXPECTED_BLUEPRINT_DECISION
        ),
        "Gate3C_cycle_audit_requires_hexagonal_redesign": (
            cycle_summary.get(
                "decision"
            )
            == EXPECTED_CYCLE_AUDIT_DECISION
        ),
        "parent_graph_has_1680_atoms": (
            len(
                parent_adjacency
            )
            == EXPECTED_PARENT_ATOMS
        ),
        "parent_graph_has_2460_bonds": (
            len(
                parent_edge_set
            )
            == EXPECTED_PARENT_BONDS
        ),
        "parent_graph_has_60_degree1_terminals": (
            degree_counts.get(
                1,
                0,
            )
            == EXPECTED_TERMINALS_TOTAL
        ),
        "parent_graph_has_1620_degree3_atoms": (
            degree_counts.get(
                3,
                0,
            )
            == EXPECTED_INTERIOR_DEGREE3
        ),
        "parent_graph_has_no_four_member_cycles": (
            parent_four_cycles == 0
        ),
        "parent_graph_is_bipartite": (
            parent_bipartite
        ),
        "parent_graph_is_connected": (
            len(
                parent_components
            )
            == 1
        ),
        "lower_selected_pairing_closes_only_six_member_cycles": (
            int(
                lower_summary[
                    "closed_cycle_length_minimum"
                ]
            )
            == EXPECTED_CLOSED_CYCLE_LENGTH
            and int(
                lower_summary[
                    "closed_cycle_length_maximum"
                ]
            )
            == EXPECTED_CLOSED_CYCLE_LENGTH
        ),
        "upper_selected_pairing_closes_only_six_member_cycles": (
            int(
                upper_summary[
                    "closed_cycle_length_minimum"
                ]
            )
            == EXPECTED_CLOSED_CYCLE_LENGTH
            and int(
                upper_summary[
                    "closed_cycle_length_maximum"
                ]
            )
            == EXPECTED_CLOSED_CYCLE_LENGTH
        ),
        "each_terminal_receives_exactly_two_new_edges": (
            all(
                int(
                    row[
                        "terminal_new_edge_incidence_minimum"
                    ]
                )
                == 2
                and int(
                    row[
                        "terminal_new_edge_incidence_maximum"
                    ]
                )
                == 2
                for row in end_summary_rows
            )
        ),
        "30_complementary_atoms_are_added_per_end": (
            all(
                int(
                    row[
                        "added_edge_completion_atoms"
                    ]
                )
                == EXPECTED_ADDED_ATOMS_PER_END
                for row in end_summary_rows
            )
        ),
        "60_complementary_atoms_are_added_total": (
            len(
                added_node_rows
            )
            == EXPECTED_ADDED_ATOMS_TOTAL
        ),
        "60_new_edges_are_added_per_end": (
            all(
                int(
                    row[
                        "new_parent_to_added_edges"
                    ]
                )
                == EXPECTED_NEW_EDGES_PER_END
                for row in end_summary_rows
            )
        ),
        "120_new_edges_are_added_total": (
            len(
                added_edge_rows
            )
            == EXPECTED_NEW_EDGES_TOTAL
        ),
        "lower_B_end_receives_N_completion_row": (
            lower_summary[
                "parent_terminal_element"
            ]
            == "B"
            and lower_summary[
                "added_complementary_element"
            ]
            == "N"
        ),
        "upper_N_end_receives_B_completion_row": (
            upper_summary[
                "parent_terminal_element"
            ]
            == "N"
            and upper_summary[
                "added_complementary_element"
            ]
            == "B"
        ),
        "all_new_edges_are_heteropolar_BN": (
            len(
                nonheteropolar_new_edges
            )
            == 0
        ),
        "all_parent_terminal_atoms_reach_degree3": (
            len(
                parent_terminal_degree_failures
            )
            == 0
        ),
        "all_parent_nonterminal_atoms_remain_degree3": (
            len(
                parent_nonterminal_degree_failures
            )
            == 0
        ),
        "all_added_completion_atoms_have_degree2": (
            len(
                added_degree_failures
            )
            == 0
        ),
        "augmented_graph_has_no_four_member_cycles": (
            augmented_four_cycles == 0
        ),
        "augmented_graph_is_bipartite": (
            augmented_bipartite
        ),
        "augmented_graph_is_connected": (
            len(
                augmented_components
            )
            == 1
        ),
        "every_added_atom_closes_at_least_one_six_member_cycle": (
            all(
                int(
                    row[
                        "length4_parent_path_multiplicity"
                    ]
                )
                >= 1
                for row in added_edge_rows[
                    ::2
                ]
            )
        ),
        "no_coordinates_were_assigned": all(
            not bool(
                row[
                    "coordinates_assigned"
                ]
            )
            for row in added_node_rows
        )
        and all(
            not bool(
                row[
                    "coordinates_assigned"
                ]
            )
            for row in added_edge_rows
        ),
        "no_formal_charges_were_assigned": all(
            not bool(
                row[
                    "formal_charge_assigned"
                ]
            )
            for row in added_node_rows
        ),
        "no_force_field_types_were_assigned": all(
            not bool(
                row[
                    "force_field_type_assigned"
                ]
            )
            for row in added_node_rows
        ),
    }

    failed_gates = [
        name
        for name, passed
        in gates.items()
        if not passed
    ]

    accepted = (
        len(
            failed_gates
        )
        == 0
    )

    decision = (
        PASS_DECISION
        if accepted
        else
        "R2_HEXAGONAL_EDGE_COMPLETION_SEED_REQUIRES_REVIEW"
    )

    required_next_step = (
        "DESIGN_AND_VALIDATE_R2_ANNULAR_CAP_ATTACHMENT_TO_STANDARDIZED_DEGREE2_RIM"
        if accepted
        else
        "REVIEW_R2_HEXAGONAL_EDGE_COMPLETION_SEED_FAILURES"
    )

    summary = {
        "decision": decision,
        "parent_atoms": (
            EXPECTED_PARENT_ATOMS
        ),
        "parent_bonds": (
            EXPECTED_PARENT_BONDS
        ),
        "parent_degree1_atoms_before_completion": (
            degree_counts.get(
                1,
                0,
            )
        ),
        "parent_degree3_atoms_before_completion": (
            degree_counts.get(
                3,
                0,
            )
        ),
        "parent_four_member_cycles": (
            parent_four_cycles
        ),
        "lower_selected_step": (
            selected_steps[
                "LOWER"
            ]
        ),
        "upper_selected_step": (
            selected_steps[
                "UPPER"
            ]
        ),
        "added_atoms_lower": (
            parse_int(
                lower_summary,
                "added_edge_completion_atoms",
            )
        ),
        "added_atoms_upper": (
            parse_int(
                upper_summary,
                "added_edge_completion_atoms",
            )
        ),
        "added_atoms_total": (
            len(
                added_node_rows
            )
        ),
        "new_edges_lower": (
            parse_int(
                lower_summary,
                "new_parent_to_added_edges",
            )
        ),
        "new_edges_upper": (
            parse_int(
                upper_summary,
                "new_parent_to_added_edges",
            )
        ),
        "new_edges_total": (
            len(
                added_edge_rows
            )
        ),
        "lower_parent_element": (
            lower_summary[
                "parent_terminal_element"
            ]
        ),
        "lower_added_element": (
            lower_summary[
                "added_complementary_element"
            ]
        ),
        "upper_parent_element": (
            upper_summary[
                "parent_terminal_element"
            ]
        ),
        "upper_added_element": (
            upper_summary[
                "added_complementary_element"
            ]
        ),
        "closed_cycle_length_minimum": (
            min(
                all_closed_cycle_lengths
            )
        ),
        "closed_cycle_length_maximum": (
            max(
                all_closed_cycle_lengths
            )
        ),
        "augmented_four_member_cycles": (
            augmented_four_cycles
        ),
        "augmented_graph_bipartite": (
            augmented_bipartite
        ),
        "augmented_graph_components": (
            len(
                augmented_components
            )
        ),
        "parent_terminal_degree_failures": (
            len(
                parent_terminal_degree_failures
            )
        ),
        "parent_nonterminal_degree_failures": (
            len(
                parent_nonterminal_degree_failures
            )
        ),
        "added_degree2_failures": (
            len(
                added_degree_failures
            )
        ),
        "standardized_degree2_rim_created": (
            accepted
        ),
        "candidate_is_final_cap_chemistry": False,
        "annular_cap_graph_authorized": (
            accepted
        ),
        "coordinate_generation_authorized": False,
        "molecular_topology_generation_authorized": False,
        "formal_charge_assignment_authorized": False,
        "force_field_parameterization_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "QM_authorized": False,
        "failed_gates": (
            " | ".join(
                failed_gates
            )
        ),
        "required_next_step": (
            required_next_step
        ),
    }

    write_csv(
        SUMMARY_CSV,
        [
            summary
        ],
    )

    write_csv(
        GATES_CSV,
        [
            {
                "gate": name,
                "pass": passed,
            }
            for name, passed
            in gates.items()
        ],
    )

    SEED_JSON.write_text(
        json.dumps(
            {
                "summary": summary,
                "end_summaries": (
                    end_summary_rows
                ),
                "gates": gates,
                "interpretation": {
                    "edge_motif": (
                        "The accepted parent exposes degree-1 "
                        "polar termini. The added complementary "
                        "row converts those termini to degree 3 "
                        "and creates a new degree-2 attachment rim."
                    ),
                    "cycle_topology": (
                        "Each added completion atom connects a "
                        "pair of parent termini separated by a "
                        "four-edge parent path, thereby closing "
                        "a six-member heavy-atom cycle."
                    ),
                    "scope": (
                        "This is a coordinate-free graph seed, "
                        "not a final cap or molecular model."
                    ),
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    source_rows = [
        {
            "role": (
                "Gate3A_parent_atoms"
            ),
            "file": relative(
                PARENT_ATOMS_CSV
            ),
            "sha256": sha256(
                PARENT_ATOMS_CSV
            ),
        },
        {
            "role": (
                "Gate3A_parent_bonds"
            ),
            "file": relative(
                PARENT_BONDS_CSV
            ),
            "sha256": sha256(
                PARENT_BONDS_CSV
            ),
        },
        {
            "role": (
                "Gate3A_terminal_atoms"
            ),
            "file": relative(
                TERMINAL_ATOMS_CSV
            ),
            "sha256": sha256(
                TERMINAL_ATOMS_CSV
            ),
        },
        {
            "role": (
                "Gate3C_blueprint_summary"
            ),
            "file": relative(
                BLUEPRINT_SUMMARY_CSV
            ),
            "sha256": sha256(
                BLUEPRINT_SUMMARY_CSV
            ),
        },
        {
            "role": (
                "Gate3C_cycle_audit_summary"
            ),
            "file": relative(
                CYCLE_AUDIT_SUMMARY_CSV
            ),
            "sha256": sha256(
                CYCLE_AUDIT_SUMMARY_CSV
            ),
        },
    ]

    write_csv(
        SOURCE_MANIFEST_CSV,
        source_rows,
    )

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed
        in gates.items()
    )

    REPORT_MD.write_text(
        f"""# R2 Hexagonal Edge-Completion Seed

## Scope

This stage replaces the rejected unconstrained bipartite collar
interface with an edge-completion motif derived from the explicit
parent BN graph.

No coordinates, formal bond orders, partial charges, molecular
topology, force-field parameters, minimization, MD, or QM calculation
were generated.

## Parent graph

- Atoms:
  **{EXPECTED_PARENT_ATOMS}**
- Bonds:
  **{EXPECTED_PARENT_BONDS}**
- Degree-1 terminal atoms:
  **{degree_counts.get(1, 0)}**
- Degree-3 interior atoms:
  **{degree_counts.get(3, 0)}**
- Four-member cycles:
  **{parent_four_cycles}**
- Bipartite:
  **{parent_bipartite}**

## Lower end

- Parent termination:
  **{lower_summary['parent_terminal_element']}**
- Added complementary row:
  **{lower_summary['added_complementary_element']}**
- Selected circumferential step:
  **{lower_summary['selected_circumferential_step']}**
- Parent terminals:
  **{lower_summary['parent_terminal_atoms']}**
- Added atoms:
  **{lower_summary['added_edge_completion_atoms']}**
- New edges:
  **{lower_summary['new_parent_to_added_edges']}**
- Closed cycle length:
  **{lower_summary['closed_cycle_length_minimum']}–
  {lower_summary['closed_cycle_length_maximum']}**

## Upper end

- Parent termination:
  **{upper_summary['parent_terminal_element']}**
- Added complementary row:
  **{upper_summary['added_complementary_element']}**
- Selected circumferential step:
  **{upper_summary['selected_circumferential_step']}**
- Parent terminals:
  **{upper_summary['parent_terminal_atoms']}**
- Added atoms:
  **{upper_summary['added_edge_completion_atoms']}**
- New edges:
  **{upper_summary['new_parent_to_added_edges']}**
- Closed cycle length:
  **{upper_summary['closed_cycle_length_minimum']}–
  {upper_summary['closed_cycle_length_maximum']}**

## Resulting graph

- Added B/N atoms total:
  **{len(added_node_rows)}**
- Added parent-to-completion edges:
  **{len(added_edge_rows)}**
- Parent-terminal degree failures:
  **{len(parent_terminal_degree_failures)}**
- Parent-nonterminal degree failures:
  **{len(parent_nonterminal_degree_failures)}**
- Added degree-2-rim failures:
  **{len(added_degree_failures)}**
- Four-member cycles:
  **{augmented_four_cycles}**
- Bipartite:
  **{augmented_bipartite}**
- Connected components:
  **{len(augmented_components)}**

The new complementary row is not the final cap. It is a validated
hexagonal edge-completion seed that converts the degree-1 polar
parent termini into a conventional degree-2 attachment rim.

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Standardized degree-2 rim created:
  **{'YES' if accepted else 'NO'}**
- Annular-cap graph design authorized:
  **{'YES' if accepted else 'NO'}**
- Coordinate generation authorized:
  **NO**
- Molecular topology generation authorized:
  **NO**
- Energy minimization authorized:
  **NO**
- MD authorized:
  **NO**
- QM authorized:
  **NO**
- Required next step:
  `{required_next_step}`
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 hexagonal edge-completion "
        "seed completed."
    )

    print(
        "Parent atoms / bonds / d1 / d3 / "
        "four-cycles: "
        f"{EXPECTED_PARENT_ATOMS}/"
        f"{EXPECTED_PARENT_BONDS}/"
        f"{degree_counts.get(1, 0)}/"
        f"{degree_counts.get(3, 0)}/"
        f"{parent_four_cycles}"
    )

    for row in end_summary_rows:
        print(
            f"{row['end']} selected step / parent / "
            "added / atoms / edges: "
            f"{row['selected_circumferential_step']}/"
            f"{row['parent_terminal_element']}/"
            f"{row['added_complementary_element']}/"
            f"{row['added_edge_completion_atoms']}/"
            f"{row['new_parent_to_added_edges']}"
        )

        print(
            f"{row['end']} parent path / closed cycle "
            "min-max: "
            f"{row['parent_pair_shortest_path_minimum']}-"
            f"{row['parent_pair_shortest_path_maximum']} / "
            f"{row['closed_cycle_length_minimum']}-"
            f"{row['closed_cycle_length_maximum']}"
        )

    print(
        "Added atoms / new edges total: "
        f"{len(added_node_rows)}/"
        f"{len(added_edge_rows)}"
    )

    print(
        "Parent terminal / nonterminal / added-rim "
        "degree failures: "
        f"{len(parent_terminal_degree_failures)}/"
        f"{len(parent_nonterminal_degree_failures)}/"
        f"{len(added_degree_failures)}"
    )

    print(
        "Augmented four-member cycles: "
        f"{augmented_four_cycles}"
    )

    print(
        "Augmented bipartite / connected components: "
        f"{augmented_bipartite}/"
        f"{len(augmented_components)}"
    )

    print(
        f"Decision: {decision}"
    )

    print(
        "Failed gates: "
        + (
            "NONE"
            if not failed_gates
            else " | ".join(
                failed_gates
            )
        )
    )

    print(
        "Standardized degree-2 rim created: "
        f"{'YES' if accepted else 'NO'}"
    )

    print(
        "Annular-cap graph design authorized: "
        f"{'YES' if accepted else 'NO'}"
    )

    print(
        "Coordinate generation authorized: NO"
    )

    print(
        "Molecular topology generation authorized: NO"
    )

    print(
        "Formal charge assignment authorized: NO"
    )

    print(
        "Force-field parameterization authorized: NO"
    )

    print(
        "Energy minimization authorized: NO"
    )

    print(
        "MD authorized: NO"
    )

    print(
        "QM authorized: NO"
    )

    print(
        f"Required next step: {required_next_step}"
    )

    for path in (
        STEP_SEARCH_CSV,
        ADDED_NODES_CSV,
        ADDED_EDGES_CSV,
        END_SUMMARY_CSV,
        SUMMARY_CSV,
        GATES_CSV,
        SEED_JSON,
        SOURCE_MANIFEST_CSV,
        REPORT_MD,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not accepted:
        raise RuntimeError(
            "R2 hexagonal edge-completion seed "
            "requires review."
        )


if __name__ == "__main__":
    main()
