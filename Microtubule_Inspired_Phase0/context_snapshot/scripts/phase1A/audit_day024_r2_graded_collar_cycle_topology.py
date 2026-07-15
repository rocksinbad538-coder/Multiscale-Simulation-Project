#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, deque
from itertools import combinations
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

INPUT_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "03_r2_graded_heteropolar_bn_collar_connectivity_blueprint"
)

OUTPUT_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "04_r2_graded_collar_cycle_topology_audit"
)

NODES_CSV = (
    INPUT_ROOT
    / "r2_selected_collar_connectivity_nodes.csv"
)

EDGES_CSV = (
    INPUT_ROOT
    / "r2_selected_collar_connectivity_edges.csv"
)

BLUEPRINT_SUMMARY_CSV = (
    INPUT_ROOT
    / "r2_graded_collar_connectivity_blueprint_summary.csv"
)

BLUEPRINT_GATES_CSV = (
    INPUT_ROOT
    / "r2_graded_collar_connectivity_blueprint_gates.csv"
)

CYCLES_CSV = (
    OUTPUT_ROOT
    / "r2_collar_simple_cycles_up_to_length10.csv"
)

NODE_PARTICIPATION_CSV = (
    OUTPUT_ROOT
    / "r2_collar_node_cycle_participation.csv"
)

EDGE_PARTICIPATION_CSV = (
    OUTPUT_ROOT
    / "r2_collar_edge_cycle_participation.csv"
)

END_SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_collar_cycle_topology_end_summary.csv"
)

AUDIT_SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_collar_cycle_topology_audit_summary.csv"
)

GATES_CSV = (
    OUTPUT_ROOT
    / "r2_collar_cycle_topology_audit_gates.csv"
)

AUDIT_JSON = (
    OUTPUT_ROOT
    / "r2_collar_cycle_topology_audit.json"
)

SOURCE_MANIFEST_CSV = (
    OUTPUT_ROOT
    / "r2_collar_cycle_topology_source_manifest.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_GRADED_COLLAR_CYCLE_TOPOLOGY_AUDIT_DAY024.md"
)

EXPECTED_BLUEPRINT_DECISION = (
    "R2_GRADED_HETEROPOLAR_BN_COLLAR_CONNECTIVITY_BLUEPRINT_VALIDATED"
)

PASS_DECISION = (
    "R2_GRADED_HETEROPOLAR_BN_COLLAR_CYCLE_TOPOLOGY_VALIDATED"
)

FAIL_DECISION = (
    "R2_GRADED_HETEROPOLAR_BN_COLLAR_REQUIRES_HEXAGONAL_GRAPH_REDESIGN"
)

MAX_ENUMERATED_CYCLE_LENGTH = 10


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


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


def canonical_cycle(
    cycle: list[int],
) -> tuple[int, ...]:
    count = len(cycle)

    candidates: list[
        tuple[int, ...]
    ] = []

    forward = list(cycle)
    reverse = list(reversed(cycle))

    for sequence in (
        forward,
        reverse,
    ):
        for offset in range(count):
            candidates.append(
                tuple(
                    sequence[offset:]
                    + sequence[:offset]
                )
            )

    return min(candidates)


def enumerate_simple_cycles(
    adjacency: dict[int, set[int]],
    maximum_length: int,
) -> set[tuple[int, ...]]:
    """
    Enumerate undirected simple cycles up to maximum_length.

    The start node is required to be the minimum-numbered node in each
    explored cycle. Canonicalization removes reversal duplicates.
    """

    cycles: set[
        tuple[int, ...]
    ] = set()

    node_ids = sorted(adjacency)

    for start in node_ids:
        path = [start]
        visited = {start}

        def walk(current: int) -> None:
            for neighbor in sorted(
                adjacency[current]
            ):
                if neighbor == start:
                    if len(path) >= 3:
                        cycles.add(
                            canonical_cycle(
                                path
                            )
                        )

                    continue

                if len(path) >= maximum_length:
                    continue

                if neighbor in visited:
                    continue

                if neighbor < start:
                    continue

                visited.add(neighbor)
                path.append(neighbor)

                walk(neighbor)

                path.pop()
                visited.remove(neighbor)

        walk(start)

    return cycles


def count_four_cycles_by_common_neighbors(
    adjacency: dict[int, set[int]],
) -> int:
    raw_count = 0

    for first, second in combinations(
        sorted(adjacency),
        2,
    ):
        common = (
            adjacency[first]
            & adjacency[second]
        )

        if len(common) >= 2:
            raw_count += (
                len(common)
                * (len(common) - 1)
                // 2
            )

    if raw_count % 2 != 0:
        raise RuntimeError(
            "Common-neighbor four-cycle count "
            "was not divisible by two."
        )

    return raw_count // 2


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


def bipartite_coloring(
    adjacency: dict[int, set[int]],
) -> tuple[
    bool,
    dict[int, int],
]:
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

                    queue.append(
                        neighbor
                    )

                elif (
                    colors[neighbor]
                    == colors[node]
                ):
                    return False, colors

    return True, colors


def cycle_edges(
    cycle: tuple[int, ...],
) -> list[tuple[int, int]]:
    return [
        tuple(
            sorted(
                (
                    cycle[index],
                    cycle[
                        (index + 1)
                        % len(cycle)
                    ],
                )
            )
        )
        for index in range(
            len(cycle)
        )
    ]


def audit_end(
    end: str,
    node_rows: list[dict[str, str]],
    edge_rows: list[dict[str, str]],
) -> tuple[
    dict[str, Any],
    dict[str, bool],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
]:
    end_nodes = [
        row
        for row in node_rows
        if row.get("end") == end
        and row.get("node_type")
        != "INNER_PASSIVANT_H"
    ]

    end_edges = [
        row
        for row in edge_rows
        if row.get("end") == end
        and parse_bool(
            row.get(
                "heavy_atom_edge",
                "false",
            )
        )
    ]

    if not end_nodes:
        raise RuntimeError(
            f"No heavy nodes found for {end}."
        )

    if not end_edges:
        raise RuntimeError(
            f"No heavy edges found for {end}."
        )

    nodes_by_id = {
        row["node_id"]: row
        for row in end_nodes
    }

    if len(nodes_by_id) != len(
        end_nodes
    ):
        raise RuntimeError(
            f"{end}: duplicate node identifiers."
        )

    sorted_node_ids = sorted(
        nodes_by_id
    )

    integer_by_node = {
        node_id: index
        for index, node_id in enumerate(
            sorted_node_ids
        )
    }

    node_by_integer = {
        index: node_id
        for node_id, index
        in integer_by_node.items()
    }

    adjacency: dict[
        int,
        set[int],
    ] = {
        index: set()
        for index in node_by_integer
    }

    edge_records: dict[
        tuple[int, int],
        dict[str, str],
    ] = {}

    duplicate_edge_count = 0
    self_edge_count = 0

    for row in end_edges:
        first_id = row[
            "source_node"
        ]

        second_id = row[
            "target_node"
        ]

        if (
            first_id not in integer_by_node
            or second_id not in integer_by_node
        ):
            raise RuntimeError(
                f"{end}: heavy edge references "
                "a missing heavy node."
            )

        first = integer_by_node[
            first_id
        ]

        second = integer_by_node[
            second_id
        ]

        if first == second:
            self_edge_count += 1
            continue

        pair = tuple(
            sorted(
                (
                    first,
                    second,
                )
            )
        )

        if pair in edge_records:
            duplicate_edge_count += 1
            continue

        edge_records[pair] = row

        adjacency[first].add(
            second
        )

        adjacency[second].add(
            first
        )

    components = connected_components(
        adjacency
    )

    is_bipartite, colors = (
        bipartite_coloring(
            adjacency
        )
    )

    cycles = enumerate_simple_cycles(
        adjacency,
        MAX_ENUMERATED_CYCLE_LENGTH,
    )

    cycle_counts = Counter(
        len(cycle)
        for cycle in cycles
    )

    enumerated_four_cycles = (
        cycle_counts.get(
            4,
            0,
        )
    )

    common_neighbor_four_cycles = (
        count_four_cycles_by_common_neighbors(
            adjacency
        )
    )

    girth = (
        min(
            len(cycle)
            for cycle in cycles
        )
        if cycles
        else 0
    )

    node_cycle_counts: dict[
        int,
        Counter[int],
    ] = {
        node: Counter()
        for node in adjacency
    }

    edge_cycle_counts: dict[
        tuple[int, int],
        Counter[int],
    ] = {
        pair: Counter()
        for pair in edge_records
    }

    cycle_rows: list[
        dict[str, Any]
    ] = []

    for cycle_index, cycle in enumerate(
        sorted(
            cycles,
            key=lambda item: (
                len(item),
                item,
            ),
        ),
        start=1,
    ):
        length = len(cycle)

        for node in cycle:
            node_cycle_counts[
                node
            ][
                length
            ] += 1

        for pair in cycle_edges(cycle):
            if pair not in edge_cycle_counts:
                raise RuntimeError(
                    f"{end}: cycle references "
                    "an unavailable edge."
                )

            edge_cycle_counts[
                pair
            ][
                length
            ] += 1

        node_ids = [
            node_by_integer[
                node
            ]
            for node in cycle
        ]

        layers = [
            int(
                float(
                    nodes_by_id[
                        node_id
                    ][
                        "layer"
                    ]
                )
            )
            for node_id in node_ids
        ]

        elements = [
            nodes_by_id[
                node_id
            ][
                "element"
            ]
            for node_id in node_ids
        ]

        cycle_rows.append(
            {
                "end": end,
                "cycle_id": (
                    f"{end}_C{cycle_index:05d}"
                ),
                "cycle_length": length,
                "node_ids": " | ".join(
                    node_ids
                ),
                "layer_sequence": " | ".join(
                    str(value)
                    for value in layers
                ),
                "element_sequence": "".join(
                    elements
                ),
            }
        )

    node_participation_rows: list[
        dict[str, Any]
    ] = []

    for integer_node in sorted(
        adjacency
    ):
        node_id = node_by_integer[
            integer_node
        ]

        source = nodes_by_id[
            node_id
        ]

        counts = node_cycle_counts[
            integer_node
        ]

        node_participation_rows.append(
            {
                "end": end,
                "node_id": node_id,
                "node_type": source[
                    "node_type"
                ],
                "layer": source[
                    "layer"
                ],
                "element": source[
                    "element"
                ],
                "heavy_degree": len(
                    adjacency[
                        integer_node
                    ]
                ),
                "cycles_length3": counts.get(
                    3,
                    0,
                ),
                "cycles_length4": counts.get(
                    4,
                    0,
                ),
                "cycles_length5": counts.get(
                    5,
                    0,
                ),
                "cycles_length6": counts.get(
                    6,
                    0,
                ),
                "cycles_length7": counts.get(
                    7,
                    0,
                ),
                "cycles_length8": counts.get(
                    8,
                    0,
                ),
                "cycles_length9": counts.get(
                    9,
                    0,
                ),
                "cycles_length10": counts.get(
                    10,
                    0,
                ),
            }
        )

    edge_participation_rows: list[
        dict[str, Any]
    ] = []

    for pair in sorted(
        edge_records
    ):
        source_row = edge_records[
            pair
        ]

        counts = edge_cycle_counts[
            pair
        ]

        edge_participation_rows.append(
            {
                "end": end,
                "source_node": (
                    node_by_integer[
                        pair[0]
                    ]
                ),
                "target_node": (
                    node_by_integer[
                        pair[1]
                    ]
                ),
                "interface": source_row.get(
                    "interface",
                    "",
                ),
                "edge_type": source_row.get(
                    "edge_type",
                    "",
                ),
                "cycles_length3": counts.get(
                    3,
                    0,
                ),
                "cycles_length4": counts.get(
                    4,
                    0,
                ),
                "cycles_length5": counts.get(
                    5,
                    0,
                ),
                "cycles_length6": counts.get(
                    6,
                    0,
                ),
                "cycles_length7": counts.get(
                    7,
                    0,
                ),
                "cycles_length8": counts.get(
                    8,
                    0,
                ),
                "cycles_length9": counts.get(
                    9,
                    0,
                ),
                "cycles_length10": counts.get(
                    10,
                    0,
                ),
            }
        )

    parent_nodes = [
        integer_by_node[
            row["node_id"]
        ]
        for row in end_nodes
        if row[
            "node_type"
        ]
        == "PARENT_TERMINAL"
    ]

    added_bn_nodes = [
        integer_by_node[
            row["node_id"]
        ]
        for row in end_nodes
        if row[
            "node_type"
        ]
        in {
            "COLLAR_BN",
            "INNER_BOUNDARY_BN",
        }
    ]

    parent_degree_failures = sum(
        len(adjacency[node]) != 2
        for node in parent_nodes
    )

    added_degree_failures = sum(
        len(adjacency[node]) != 3
        for node in added_bn_nodes
    )

    nodes_in_four_cycles = {
        node
        for cycle in cycles
        if len(cycle) == 4
        for node in cycle
    }

    added_nodes_in_four_cycles = sum(
        node in nodes_in_four_cycles
        for node in added_bn_nodes
    )

    parent_collar_edges_in_four_cycles = sum(
        row.get(
            "edge_type",
            "",
        )
        == "PARENT_TO_COLLAR"
        and int(
            row[
                "cycles_length4"
            ]
        )
        > 0
        for row in edge_participation_rows
    )

    odd_cycle_count = sum(
        count
        for length, count
        in cycle_counts.items()
        if length % 2 == 1
    )

    cycle_rank = (
        len(edge_records)
        - len(adjacency)
        + len(components)
    )

    gates = {
        f"{end}_heavy_graph_is_connected": (
            len(components) == 1
        ),
        f"{end}_heavy_graph_is_bipartite": (
            is_bipartite
        ),
        f"{end}_has_no_self_edges": (
            self_edge_count == 0
        ),
        f"{end}_has_no_duplicate_edges": (
            duplicate_edge_count == 0
        ),
        f"{end}_parent_blueprint_degree_is_two": (
            parent_degree_failures == 0
        ),
        f"{end}_added_BN_degree_is_three": (
            added_degree_failures == 0
        ),
        f"{end}_four_cycle_enumerators_agree": (
            enumerated_four_cycles
            == common_neighbor_four_cycles
        ),
        f"{end}_contains_no_triangles": (
            cycle_counts.get(
                3,
                0,
            )
            == 0
        ),
        f"{end}_contains_no_four_member_cycles": (
            enumerated_four_cycles == 0
        ),
        f"{end}_contains_no_five_member_cycles": (
            cycle_counts.get(
                5,
                0,
            )
            == 0
        ),
        f"{end}_contains_no_odd_cycles_up_to_length10": (
            odd_cycle_count == 0
        ),
        f"{end}_girth_is_at_least_six": (
            girth >= 6
        ),
        f"{end}_contains_hexagonal_cycles": (
            cycle_counts.get(
                6,
                0,
            )
            > 0
        ),
        f"{end}_added_BN_nodes_avoid_four_cycles": (
            added_nodes_in_four_cycles == 0
        ),
        f"{end}_parent_collar_edges_avoid_four_cycles": (
            parent_collar_edges_in_four_cycles
            == 0
        ),
    }

    summary = {
        "end": end,
        "heavy_nodes": len(
            adjacency
        ),
        "heavy_edges": len(
            edge_records
        ),
        "connected_components": len(
            components
        ),
        "cycle_rank": cycle_rank,
        "bipartite": is_bipartite,
        "girth": girth,
        "cycles_length3": (
            cycle_counts.get(
                3,
                0,
            )
        ),
        "cycles_length4": (
            cycle_counts.get(
                4,
                0,
            )
        ),
        "cycles_length5": (
            cycle_counts.get(
                5,
                0,
            )
        ),
        "cycles_length6": (
            cycle_counts.get(
                6,
                0,
            )
        ),
        "cycles_length7": (
            cycle_counts.get(
                7,
                0,
            )
        ),
        "cycles_length8": (
            cycle_counts.get(
                8,
                0,
            )
        ),
        "cycles_length9": (
            cycle_counts.get(
                9,
                0,
            )
        ),
        "cycles_length10": (
            cycle_counts.get(
                10,
                0,
            )
        ),
        "enumerated_four_cycles": (
            enumerated_four_cycles
        ),
        "common_neighbor_four_cycles": (
            common_neighbor_four_cycles
        ),
        "added_BN_nodes_in_four_cycles": (
            added_nodes_in_four_cycles
        ),
        "parent_collar_edges_in_four_cycles": (
            parent_collar_edges_in_four_cycles
        ),
        "parent_degree_failures": (
            parent_degree_failures
        ),
        "added_BN_degree_failures": (
            added_degree_failures
        ),
        "self_edges": (
            self_edge_count
        ),
        "duplicate_edges": (
            duplicate_edge_count
        ),
        "maximum_enumerated_cycle_length": (
            MAX_ENUMERATED_CYCLE_LENGTH
        ),
    }

    return (
        summary,
        gates,
        cycle_rows,
        node_participation_rows,
        edge_participation_rows,
    )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        NODES_CSV,
        EDGES_CSV,
        BLUEPRINT_SUMMARY_CSV,
        BLUEPRINT_GATES_CSV,
    ):
        require_file(required)

    nodes = read_csv_rows(
        NODES_CSV
    )

    edges = read_csv_rows(
        EDGES_CSV
    )

    blueprint = read_single_csv_row(
        BLUEPRINT_SUMMARY_CSV
    )

    blueprint_gates = read_csv_rows(
        BLUEPRINT_GATES_CSV
    )

    if blueprint.get(
        "decision"
    ) != EXPECTED_BLUEPRINT_DECISION:
        raise RuntimeError(
            "Gate 3C is not in the accepted state."
        )

    failed_upstream_gates = [
        row.get("gate", "")
        for row in blueprint_gates
        if not parse_bool(
            row.get(
                "pass",
                "false",
            )
        )
    ]

    if failed_upstream_gates:
        raise RuntimeError(
            "Gate 3C contains failed upstream gates: "
            + " | ".join(
                failed_upstream_gates
            )
        )

    all_end_summaries: list[
        dict[str, Any]
    ] = []

    all_gates: dict[
        str,
        bool
    ] = {}

    all_cycle_rows: list[
        dict[str, Any]
    ] = []

    all_node_rows: list[
        dict[str, Any]
    ] = []

    all_edge_rows: list[
        dict[str, Any]
    ] = []

    for end in (
        "LOWER",
        "UPPER",
    ):
        (
            end_summary,
            end_gates,
            cycle_rows,
            node_rows,
            edge_rows,
        ) = audit_end(
            end,
            nodes,
            edges,
        )

        all_end_summaries.append(
            end_summary
        )

        all_gates.update(
            end_gates
        )

        all_cycle_rows.extend(
            cycle_rows
        )

        all_node_rows.extend(
            node_rows
        )

        all_edge_rows.extend(
            edge_rows
        )

    all_gates[
        "Gate3C_blueprint_was_previously_accepted"
    ] = (
        blueprint.get(
            "decision"
        )
        == EXPECTED_BLUEPRINT_DECISION
    )

    all_gates[
        "Gate3C_had_no_failed_upstream_gates"
    ] = (
        len(
            failed_upstream_gates
        )
        == 0
    )

    lower = next(
        row
        for row in all_end_summaries
        if row["end"] == "LOWER"
    )

    upper = next(
        row
        for row in all_end_summaries
        if row["end"] == "UPPER"
    )

    all_gates[
        "lower_and_upper_cycle_statistics_are_symmetric"
    ] = all(
        lower[key] == upper[key]
        for key in (
            "heavy_nodes",
            "heavy_edges",
            "cycle_rank",
            "girth",
            "cycles_length3",
            "cycles_length4",
            "cycles_length5",
            "cycles_length6",
            "cycles_length7",
            "cycles_length8",
            "cycles_length9",
            "cycles_length10",
        )
    )

    failed_gates = [
        name
        for name, passed
        in all_gates.items()
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
        else FAIL_DECISION
    )

    required_next_step = (
        "BUILD_AND_VALIDATE_R2_GRADED_HETEROPOLAR_COLLAR_STATIC_COORDINATE_EMBEDDING"
        if accepted
        else
        "REDESIGN_R2_COLLAR_GRAPH_USING_HEXAGONAL_LATTICE_TEMPLATE_AND_REAUDIT"
    )

    write_csv(
        CYCLES_CSV,
        all_cycle_rows,
    )

    write_csv(
        NODE_PARTICIPATION_CSV,
        all_node_rows,
    )

    write_csv(
        EDGE_PARTICIPATION_CSV,
        all_edge_rows,
    )

    write_csv(
        END_SUMMARY_CSV,
        all_end_summaries,
    )

    write_csv(
        GATES_CSV,
        [
            {
                "gate": name,
                "pass": passed,
            }
            for name, passed
            in all_gates.items()
        ],
    )

    summary = {
        "decision": decision,
        "Gate3C_blueprint_decision": (
            blueprint[
                "decision"
            ]
        ),
        "lower_girth": (
            lower[
                "girth"
            ]
        ),
        "upper_girth": (
            upper[
                "girth"
            ]
        ),
        "lower_four_member_cycles": (
            lower[
                "cycles_length4"
            ]
        ),
        "upper_four_member_cycles": (
            upper[
                "cycles_length4"
            ]
        ),
        "lower_six_member_cycles": (
            lower[
                "cycles_length6"
            ]
        ),
        "upper_six_member_cycles": (
            upper[
                "cycles_length6"
            ]
        ),
        "lower_eight_member_cycles": (
            lower[
                "cycles_length8"
            ]
        ),
        "upper_eight_member_cycles": (
            upper[
                "cycles_length8"
            ]
        ),
        "lower_added_BN_nodes_in_four_cycles": (
            lower[
                "added_BN_nodes_in_four_cycles"
            ]
        ),
        "upper_added_BN_nodes_in_four_cycles": (
            upper[
                "added_BN_nodes_in_four_cycles"
            ]
        ),
        "failed_gates": (
            " | ".join(
                failed_gates
            )
        ),
        "cycle_topology_accepted": (
            accepted
        ),
        "previous_coordinate_embedding_authorization_superseded": (
            not accepted
        ),
        "static_coordinate_embedding_authorized": (
            accepted
        ),
        "molecular_topology_generation_authorized": False,
        "formal_charge_assignment_authorized": False,
        "force_field_parameterization_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "QM_authorized": False,
        "required_next_step": (
            required_next_step
        ),
    }

    write_csv(
        AUDIT_SUMMARY_CSV,
        [
            summary
        ],
    )

    payload = {
        "summary": summary,
        "end_summaries": (
            all_end_summaries
        ),
        "gates": all_gates,
        "interpretation": {
            "four_member_cycles": (
                "Four-member B-N-B-N cycles are not accepted "
                "for the current low-strain h-BN-like collar "
                "hypothesis."
            ),
            "six_member_cycles": (
                "Six-member cycles are compatible with the "
                "intended hexagonal-network motif."
            ),
            "scope": (
                "This is a graph-topology audit. It does not "
                "assign coordinates, energies, charges, or "
                "force-field parameters."
            ),
        },
    }

    AUDIT_JSON.write_text(
        json.dumps(
            payload,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    source_rows = [
        {
            "role": (
                "Gate3C_nodes"
            ),
            "file": relative(
                NODES_CSV
            ),
            "sha256": sha256(
                NODES_CSV
            ),
        },
        {
            "role": (
                "Gate3C_edges"
            ),
            "file": relative(
                EDGES_CSV
            ),
            "sha256": sha256(
                EDGES_CSV
            ),
        },
        {
            "role": (
                "Gate3C_summary"
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
                "Gate3C_gates"
            ),
            "file": relative(
                BLUEPRINT_GATES_CSV
            ),
            "sha256": sha256(
                BLUEPRINT_GATES_CSV
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
        in all_gates.items()
    )

    REPORT_MD.write_text(
        f"""# R2 Graded-Collar Cycle-Topology Audit

## Scope

This stage audits the heavy-atom graph produced by Gate 3C.

No coordinates, formal bond orders, charges, force-field parameters,
minimization, MD, or QM calculation were generated.

## Lower end

- Heavy nodes/edges:
  **{lower['heavy_nodes']}/{lower['heavy_edges']}**
- Connected components:
  **{lower['connected_components']}**
- Bipartite:
  **{lower['bipartite']}**
- Cycle rank:
  **{lower['cycle_rank']}**
- Girth:
  **{lower['girth']}**
- Cycles of length 3/4/5/6/7/8:
  **{lower['cycles_length3']}/
  {lower['cycles_length4']}/
  {lower['cycles_length5']}/
  {lower['cycles_length6']}/
  {lower['cycles_length7']}/
  {lower['cycles_length8']}**
- Added BN nodes participating in four-member cycles:
  **{lower['added_BN_nodes_in_four_cycles']}**

## Upper end

- Heavy nodes/edges:
  **{upper['heavy_nodes']}/{upper['heavy_edges']}**
- Connected components:
  **{upper['connected_components']}**
- Bipartite:
  **{upper['bipartite']}**
- Cycle rank:
  **{upper['cycle_rank']}**
- Girth:
  **{upper['girth']}**
- Cycles of length 3/4/5/6/7/8:
  **{upper['cycles_length3']}/
  {upper['cycles_length4']}/
  {upper['cycles_length5']}/
  {upper['cycles_length6']}/
  {upper['cycles_length7']}/
  {upper['cycles_length8']}**
- Added BN nodes participating in four-member cycles:
  **{upper['added_BN_nodes_in_four_cycles']}**

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Previous coordinate-embedding authorization superseded:
  **{'YES' if not accepted else 'NO'}**
- Static coordinate embedding authorized:
  **{'YES' if accepted else 'NO'}**
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

## Interpretation

Gate 3C established an abstract trivalent, connected and heteropolar
graph. This audit determines whether that graph also has a
low-strain hexagonal-network cycle topology.

The presence of four-member B-N-B-N cycles blocks direct promotion to
the coordinate-embedding stage. In that event, the graph must be
rebuilt from an explicit hexagonal-lattice or nanotube-junction
template rather than from unconstrained bipartite interface flows.
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 graded-collar cycle-topology "
        "audit completed."
    )

    for end_summary in all_end_summaries:
        print(
            f"{end_summary['end']} heavy nodes/edges/"
            "components/bipartite: "
            f"{end_summary['heavy_nodes']}/"
            f"{end_summary['heavy_edges']}/"
            f"{end_summary['connected_components']}/"
            f"{end_summary['bipartite']}"
        )

        print(
            f"{end_summary['end']} girth and cycles "
            "L3/L4/L5/L6/L7/L8: "
            f"{end_summary['girth']} | "
            f"{end_summary['cycles_length3']}/"
            f"{end_summary['cycles_length4']}/"
            f"{end_summary['cycles_length5']}/"
            f"{end_summary['cycles_length6']}/"
            f"{end_summary['cycles_length7']}/"
            f"{end_summary['cycles_length8']}"
        )

        print(
            f"{end_summary['end']} added BN nodes in "
            "four-member cycles: "
            f"{end_summary['added_BN_nodes_in_four_cycles']}"
        )

        print(
            f"{end_summary['end']} parent-collar edges in "
            "four-member cycles: "
            f"{end_summary['parent_collar_edges_in_four_cycles']}"
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
        "Previous coordinate-embedding authorization "
        "superseded: "
        f"{'YES' if not accepted else 'NO'}"
    )

    print(
        "Static coordinate embedding authorized: "
        f"{'YES' if accepted else 'NO'}"
    )

    print(
        "Molecular topology generation authorized: NO"
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
        CYCLES_CSV,
        NODE_PARTICIPATION_CSV,
        EDGE_PARTICIPATION_CSV,
        END_SUMMARY_CSV,
        AUDIT_SUMMARY_CSV,
        GATES_CSV,
        AUDIT_JSON,
        SOURCE_MANIFEST_CSV,
        REPORT_MD,
    ):
        print(
            f"Wrote: {relative(path)}"
        )


if __name__ == "__main__":
    main()
