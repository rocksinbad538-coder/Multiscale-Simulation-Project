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

GATE3F_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "07_r2_reconstruction_vs_partial_attachment_contingency"
)

GATE3G1_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "09_r2_direct_junction_geometric_lower_bound"
)

GATE3H_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "10_r2_alternating_bn_oligomer_bridge_feasibility"
)

OUTPUT_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "11_r2_alternating_bn_trimer_bridge_graph"
)

DESIGN_NODES_CSV = (
    GATE3F_ROOT
    / "r2_partial_attachment_passivated_annulus_nodes.csv"
)

DESIGN_EDGES_CSV = (
    GATE3F_ROOT
    / "r2_partial_attachment_passivated_annulus_edges.csv"
)

DESIGN_SUMMARY_CSV = (
    GATE3F_ROOT
    / "r2_reconstruction_vs_partial_attachment_summary.csv"
)

DIRECT_LOWER_BOUND_SUMMARY_CSV = (
    GATE3G1_ROOT
    / "r2_direct_junction_geometric_lower_bound_summary.csv"
)

BRIDGE_FEASIBILITY_SUMMARY_CSV = (
    GATE3H_ROOT
    / "r2_bn_oligomer_bridge_feasibility_summary.csv"
)

SELECTED_BRIDGE_CANDIDATE_CSV = (
    GATE3H_ROOT
    / "r2_bn_oligomer_bridge_selected_candidate.csv"
)

GRAPH_NODES_CSV = (
    OUTPUT_ROOT
    / "r2_alternating_bn_trimer_bridge_graph_nodes.csv"
)

GRAPH_EDGES_CSV = (
    OUTPUT_ROOT
    / "r2_alternating_bn_trimer_bridge_graph_edges.csv"
)

BRIDGE_PATHS_CSV = (
    OUTPUT_ROOT
    / "r2_alternating_bn_trimer_bridge_paths.csv"
)

END_SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_alternating_bn_trimer_bridge_end_summary.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_alternating_bn_trimer_bridge_graph_summary.csv"
)

GATES_CSV = (
    OUTPUT_ROOT
    / "r2_alternating_bn_trimer_bridge_graph_gates.csv"
)

GRAPH_JSON = (
    OUTPUT_ROOT
    / "r2_alternating_bn_trimer_bridge_graph.json"
)

SOURCE_MANIFEST_CSV = (
    OUTPUT_ROOT
    / "r2_alternating_bn_trimer_bridge_source_manifest.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_ALTERNATING_BN_TRIMER_BRIDGE_GRAPH_DAY024.md"
)

EXPECTED_DESIGN_DECISION = (
    "R2_PARTIAL_HETEROPOLAR_ANNULUS_ATTACHMENT_AND_"
    "COMPLEMENTARY_PASSIVATION_GRAPH_VALIDATED"
)

EXPECTED_DIRECT_REJECTION = (
    "R2_PARTIAL_ATTACHMENT_DIRECT_BN_JUNCTION_GEOMETRIC_"
    "LOWER_BOUND_FAILED"
)

EXPECTED_BRIDGE_DECISION = (
    "R2_SHORTEST_ALTERNATING_BN_OLIGOMER_BRIDGE_CLASS_IDENTIFIED"
)

PASS_DECISION = (
    "R2_ALTERNATING_BN_TRIMER_BRIDGE_GRAPH_VALIDATED"
)

EXPECTED_BRIDGE_ATOMS_PER_PATH = 3
EXPECTED_BRIDGE_PATHS_PER_END = 15
EXPECTED_BRIDGE_PATHS_TOTAL = 30

EXPECTED_BRIDGE_ATOMS_PER_END = (
    EXPECTED_BRIDGE_ATOMS_PER_PATH
    * EXPECTED_BRIDGE_PATHS_PER_END
)

EXPECTED_BRIDGE_ATOMS_TOTAL = (
    2
    * EXPECTED_BRIDGE_ATOMS_PER_END
)

EXPECTED_BRIDGE_EDGES_PER_PATH = 4
EXPECTED_BRIDGE_EDGES_PER_END = (
    EXPECTED_BRIDGE_PATHS_PER_END
    * EXPECTED_BRIDGE_EDGES_PER_PATH
)

EXPECTED_BRIDGE_EDGES_TOTAL = (
    2
    * EXPECTED_BRIDGE_EDGES_PER_END
)

EXPECTED_PARENT_ATOMS = 1680
EXPECTED_SEED_ATOMS_TOTAL = 60
EXPECTED_ANNULUS_ATOMS_TOTAL = 252

EXPECTED_HEAVY_ATOMS_PER_END = 201
EXPECTED_ADDED_HEAVY_ATOMS_TOTAL = 402
EXPECTED_TOTAL_HEAVY_ATOMS = (
    EXPECTED_PARENT_ATOMS
    + EXPECTED_ADDED_HEAVY_ATOMS_TOTAL
)

EXPECTED_H_PER_END = 87
EXPECTED_H_TOTAL = 174

EXPECTED_TOTAL_NODES = (
    EXPECTED_TOTAL_HEAVY_ATOMS
    + EXPECTED_H_TOTAL
)

EXPECTED_SEED_H_PER_END = 15
EXPECTED_OUTER_H_PER_END = 15
EXPECTED_INNER_H_PER_END = 12
EXPECTED_BRIDGE_H_PER_END = 45

MINIMUM_ACCEPTED_HEAVY_GIRTH = 6
MINIMUM_BRIDGE_CONTAINING_CYCLE_LENGTH = 6


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


def parse_float(
    row: dict[str, str],
    key: str,
) -> float:
    try:
        value = float(row[key])
    except (
        KeyError,
        TypeError,
        ValueError,
    ) as exc:
        raise RuntimeError(
            f"Could not parse numeric field {key!r}"
        ) from exc

    if not math.isfinite(value):
        raise RuntimeError(
            f"Non-finite field {key!r}"
        )

    return value


def opposite_element(element: str) -> str:
    if element == "B":
        return "N"

    if element == "N":
        return "B"

    raise RuntimeError(
        f"Unexpected BN element: {element}"
    )


def connected_components(
    adjacency: dict[str, set[str]],
) -> list[set[str]]:
    remaining = set(adjacency)
    components: list[set[str]] = []

    while remaining:
        start = min(remaining)
        component: set[str] = set()
        queue: deque[str] = deque([start])

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
    adjacency: dict[str, set[str]],
) -> tuple[
    bool,
    dict[str, int],
]:
    colors: dict[str, int] = {}

    for start in sorted(adjacency):
        if start in colors:
            continue

        colors[start] = 0
        queue: deque[str] = deque([start])

        while queue:
            node = queue.popleft()

            for neighbor in adjacency[node]:
                if neighbor not in colors:
                    colors[neighbor] = (
                        1 - colors[node]
                    )

                    queue.append(neighbor)

                elif colors[neighbor] == colors[node]:
                    return False, colors

    return True, colors


def count_four_cycles(
    adjacency: dict[str, set[str]],
) -> int:
    nodes = sorted(adjacency)
    raw_count = 0

    for first_index, first in enumerate(nodes):
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


def shortest_path_length(
    adjacency: dict[str, set[str]],
    source: str,
    target: str,
    excluded_edges: set[
        tuple[str, str]
    ] | None = None,
) -> int | None:
    if source == target:
        return 0

    excluded = (
        excluded_edges
        if excluded_edges is not None
        else set()
    )

    visited = {source}

    queue: deque[
        tuple[str, int]
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
            pair = tuple(
                sorted(
                    (
                        node,
                        neighbor,
                    )
                )
            )

            if pair in excluded:
                continue

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


def contains_six_cycle(
    adjacency: dict[str, set[str]],
) -> bool:
    for first in adjacency:
        for second in adjacency[first]:
            if first >= second:
                continue

            excluded = {
                tuple(
                    sorted(
                        (
                            first,
                            second,
                        )
                    )
                )
            }

            alternative = shortest_path_length(
                adjacency,
                first,
                second,
                excluded,
            )

            if alternative == 5:
                return True

    return False


def add_node(
    node_rows: list[dict[str, Any]],
    node_elements: dict[str, str],
    node_types: dict[str, str],
    adjacency: dict[str, set[str]],
    *,
    node_id: str,
    element: str,
    node_type: str,
    end: str,
    metadata: dict[str, Any],
) -> None:
    if node_id in node_elements:
        raise RuntimeError(
            f"Duplicate node identifier: {node_id}"
        )

    node_elements[node_id] = element
    node_types[node_id] = node_type
    adjacency[node_id] = set()

    node_rows.append(
        {
            "node_id": node_id,
            "element": element,
            "node_type": node_type,
            "end": end,
            **metadata,
            "coordinates_assigned": False,
            "formal_charge_assigned": False,
            "force_field_type_assigned": False,
        }
    )


def add_edge(
    edge_rows: list[dict[str, Any]],
    edge_pairs: set[
        tuple[str, str]
    ],
    adjacency: dict[str, set[str]],
    node_elements: dict[str, str],
    *,
    source: str,
    target: str,
    edge_type: str,
    end: str,
    heavy_atom_edge: bool,
    metadata: dict[str, Any] | None = None,
) -> None:
    if source == target:
        raise RuntimeError(
            f"Self-edge requested for {source}"
        )

    if (
        source not in adjacency
        or target not in adjacency
    ):
        raise RuntimeError(
            f"Edge references missing nodes: "
            f"{source} | {target}"
        )

    pair = tuple(
        sorted(
            (
                source,
                target,
            )
        )
    )

    if pair in edge_pairs:
        raise RuntimeError(
            f"Duplicate edge requested: {pair}"
        )

    edge_pairs.add(pair)

    adjacency[source].add(target)
    adjacency[target].add(source)

    source_element = node_elements[source]
    target_element = node_elements[target]

    edge_rows.append(
        {
            "edge_id": (
                f"E:{len(edge_rows) + 1}"
            ),
            "source_node": source,
            "target_node": target,
            "source_element": source_element,
            "target_element": target_element,
            "edge_type": edge_type,
            "end": end,
            "heavy_atom_edge": heavy_atom_edge,
            "heteropolar_BN_edge": (
                heavy_atom_edge
                and {
                    source_element,
                    target_element,
                }
                == {
                    "B",
                    "N",
                }
            ),
            **(
                metadata
                if metadata is not None
                else {}
            ),
            "coordinates_assigned": False,
            "formal_bond_order_assigned": False,
        }
    )


def add_hydrogen(
    *,
    node_rows: list[dict[str, Any]],
    edge_rows: list[dict[str, Any]],
    edge_pairs: set[
        tuple[str, str]
    ],
    adjacency: dict[str, set[str]],
    node_elements: dict[str, str],
    node_types: dict[str, str],
    heavy_node: str,
    end: str,
    hydrogen_role: str,
    hydrogen_index: int,
) -> str:
    hydrogen_id = (
        f"H:{end}:{hydrogen_role}:"
        f"{hydrogen_index:04d}"
    )

    add_node(
        node_rows,
        node_elements,
        node_types,
        adjacency,
        node_id=hydrogen_id,
        element="H",
        node_type=hydrogen_role,
        end=end,
        metadata={
            "attached_to": heavy_node,
        },
    )

    add_edge(
        edge_rows,
        edge_pairs,
        adjacency,
        node_elements,
        source=heavy_node,
        target=hydrogen_id,
        edge_type=(
            f"{hydrogen_role}_BOND"
        ),
        end=end,
        heavy_atom_edge=False,
    )

    return hydrogen_id


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        DESIGN_NODES_CSV,
        DESIGN_EDGES_CSV,
        DESIGN_SUMMARY_CSV,
        DIRECT_LOWER_BOUND_SUMMARY_CSV,
        BRIDGE_FEASIBILITY_SUMMARY_CSV,
        SELECTED_BRIDGE_CANDIDATE_CSV,
    ):
        require_file(required)

    design_nodes = read_csv_rows(
        DESIGN_NODES_CSV
    )

    design_edges = read_csv_rows(
        DESIGN_EDGES_CSV
    )

    design_summary = read_single_csv_row(
        DESIGN_SUMMARY_CSV
    )

    direct_summary = read_single_csv_row(
        DIRECT_LOWER_BOUND_SUMMARY_CSV
    )

    bridge_summary = read_single_csv_row(
        BRIDGE_FEASIBILITY_SUMMARY_CSV
    )

    selected_candidate_rows = read_csv_rows(
        SELECTED_BRIDGE_CANDIDATE_CSV
    )

    if design_summary.get(
        "decision"
    ) != EXPECTED_DESIGN_DECISION:
        raise RuntimeError(
            "Gate 3F graph design is not accepted."
        )

    if direct_summary.get(
        "decision"
    ) != EXPECTED_DIRECT_REJECTION:
        raise RuntimeError(
            "Gate 3G.1 does not contain the "
            "expected direct-junction rejection."
        )

    if bridge_summary.get(
        "decision"
    ) != EXPECTED_BRIDGE_DECISION:
        raise RuntimeError(
            "Gate 3H bridge screen is not accepted."
        )

    selected_bridge_atoms = parse_int(
        bridge_summary,
        "selected_shortest_bridge_atoms_per_attachment",
    )

    if (
        selected_bridge_atoms
        != EXPECTED_BRIDGE_ATOMS_PER_PATH
    ):
        raise RuntimeError(
            "Unexpected selected bridge class: "
            f"{selected_bridge_atoms}/"
            f"{EXPECTED_BRIDGE_ATOMS_PER_PATH}"
        )

    selected_mappings = {
        row["end"]: row
        for row in selected_candidate_rows
        if row.get(
            "classification"
        )
        == (
            "SELECTED_SHORTEST_FEASIBLE_"
            "BRIDGE_CLASS_MAPPING"
        )
    }

    if set(selected_mappings) != {
        "LOWER",
        "UPPER",
    }:
        raise RuntimeError(
            "Could not resolve both selected end mappings."
        )

    original_nodes_by_id = {
        row[
            "node_id"
        ]: row
        for row in design_nodes
    }

    if len(original_nodes_by_id) != len(
        design_nodes
    ):
        raise RuntimeError(
            "Duplicate node identifiers in Gate 3F."
        )

    node_rows: list[
        dict[str, Any]
    ] = []

    edge_rows: list[
        dict[str, Any]
    ] = []

    node_elements: dict[
        str,
        str
    ] = {}

    node_types: dict[
        str,
        str
    ] = {}

    adjacency: dict[
        str,
        set[str]
    ] = {}

    edge_pairs: set[
        tuple[str, str]
    ] = set()

    retained_heavy_node_ids = []

    for row in design_nodes:
        if row[
            "element"
        ] == "H":
            continue

        node_id = row[
            "node_id"
        ]

        retained_heavy_node_ids.append(
            node_id
        )

        metadata = {
            key: value
            for key, value
            in row.items()
            if key
            not in {
                "node_id",
                "element",
                "node_type",
                "end",
                "coordinates_assigned",
                "formal_charge_assigned",
                "force_field_type_assigned",
            }
            and value != ""
        }

        add_node(
            node_rows,
            node_elements,
            node_types,
            adjacency,
            node_id=node_id,
            element=row[
                "element"
            ],
            node_type=row[
                "node_type"
            ],
            end=row[
                "end"
            ],
            metadata={
                **metadata,
                "graph_source": (
                    "GATE3F_RETAINED_HEAVY_NODE"
                ),
            },
        )

    direct_edges_removed = 0
    passivation_edges_removed = 0
    retained_heavy_edges = 0

    for row in design_edges:
        edge_type = row[
            "edge_type"
        ]

        heavy_atom_edge = parse_bool(
            row[
                "heavy_atom_edge"
            ]
        )

        if (
            edge_type
            == "PARTIAL_HETEROPOLAR_SEED_TO_ANNULUS"
        ):
            direct_edges_removed += 1
            continue

        if not heavy_atom_edge:
            passivation_edges_removed += 1
            continue

        source = row[
            "source_node"
        ]

        target = row[
            "target_node"
        ]

        if (
            source not in adjacency
            or target not in adjacency
        ):
            raise RuntimeError(
                "Retained heavy edge references "
                "a removed node."
            )

        retained_heavy_edges += 1

        add_edge(
            edge_rows,
            edge_pairs,
            adjacency,
            node_elements,
            source=source,
            target=target,
            edge_type=edge_type,
            end=row[
                "end"
            ],
            heavy_atom_edge=True,
            metadata={
                "graph_source": (
                    "GATE3F_RETAINED_HEAVY_EDGE"
                ),
            },
        )

    if direct_edges_removed != 30:
        raise RuntimeError(
            "Unexpected number of rejected direct edges: "
            f"{direct_edges_removed}/30"
        )

    bridge_path_rows: list[
        dict[str, Any]
    ] = []

    selected_seed_ids_by_end: dict[
        str,
        set[str]
    ] = {}

    selected_outer_ids_by_end: dict[
        str,
        set[str]
    ] = {}

    bridge_node_ids_by_end: dict[
        str,
        list[str]
    ] = {
        "LOWER": [],
        "UPPER": [],
    }

    bridge_path_edge_sets: dict[
        str,
        set[tuple[str, str]]
    ] = {}

    end_preliminary_rows = []

    for end in (
        "LOWER",
        "UPPER",
    ):
        mapping = selected_mappings[
            end
        ]

        bridge_atoms = parse_int(
            mapping,
            "bridge_atoms_per_attachment",
        )

        if (
            bridge_atoms
            != EXPECTED_BRIDGE_ATOMS_PER_PATH
        ):
            raise RuntimeError(
                f"{end}: selected mapping uses "
                f"{bridge_atoms} bridge atoms."
            )

        seed_parity = parse_int(
            mapping,
            "seed_parity",
        )

        orientation = parse_int(
            mapping,
            "orientation",
        )

        rotation_index = parse_int(
            mapping,
            "rotation_index",
        )

        chirality = parse_int(
            mapping,
            "chirality",
        )

        selected_gap_nm = parse_float(
            mapping,
            "selected_gap_nm",
        )

        feasible_gap_minimum_nm = parse_float(
            mapping,
            "feasible_gap_minimum_nm",
        )

        feasible_gap_maximum_nm = parse_float(
            mapping,
            "feasible_gap_maximum_nm",
        )

        seed_rows = [
            row
            for row in design_nodes
            if row[
                "end"
            ]
            == end
            and row[
                "node_type"
            ]
            == "HEXAGONAL_EDGE_COMPLETION_SEED"
        ]

        seed_rows.sort(
            key=lambda row: parse_int(
                row,
                "circumferential_index",
            )
        )

        if len(seed_rows) != 30:
            raise RuntimeError(
                f"{end}: seed population "
                f"{len(seed_rows)}/30"
            )

        seed_elements = {
            row[
                "element"
            ]
            for row in seed_rows
        }

        if len(seed_elements) != 1:
            raise RuntimeError(
                f"{end}: seed is not elementally homogeneous."
            )

        seed_element = next(
            iter(
                seed_elements
            )
        )

        bridge_sequence = mapping[
            "bridge_element_sequence"
        ].split("-")

        if len(bridge_sequence) != 3:
            raise RuntimeError(
                f"{end}: unexpected bridge sequence "
                f"{bridge_sequence}"
            )

        expected_sequence = []

        current_element = seed_element

        for _ in range(3):
            current_element = opposite_element(
                current_element
            )

            expected_sequence.append(
                current_element
            )

        if bridge_sequence != expected_sequence:
            raise RuntimeError(
                f"{end}: bridge sequence "
                f"{bridge_sequence} does not match "
                f"{expected_sequence}"
            )

        required_annulus_element = mapping[
            "required_annulus_element"
        ]

        if (
            required_annulus_element
            != seed_element
        ):
            raise RuntimeError(
                f"{end}: a three-atom bridge should "
                "terminate at the same endpoint element "
                "as the seed."
            )

        selected_seed_rows = [
            row
            for row in seed_rows
            if parse_int(
                row,
                "circumferential_index",
            )
            % 2
            == seed_parity
        ]

        if (
            len(
                selected_seed_rows
            )
            != EXPECTED_BRIDGE_PATHS_PER_END
        ):
            raise RuntimeError(
                f"{end}: parity selection produced "
                f"{len(selected_seed_rows)}/"
                f"{EXPECTED_BRIDGE_PATHS_PER_END}"
            )

        outer_rows = [
            row
            for row in design_nodes
            if row[
                "end"
            ]
            == end
            and row[
                "node_type"
            ]
            == "ANNULUS_OUTER_BOUNDARY"
            and row[
                "element"
            ]
            == required_annulus_element
        ]

        outer_rows.sort(
            key=lambda row: parse_float(
                row,
                "angle_turns",
            )
        )

        if (
            len(
                outer_rows
            )
            != EXPECTED_BRIDGE_PATHS_PER_END
        ):
            raise RuntimeError(
                f"{end}: required annulus endpoint "
                f"population {len(outer_rows)}/"
                f"{EXPECTED_BRIDGE_PATHS_PER_END}"
            )

        mapped_outer_rows = [
            outer_rows[
                (
                    orientation
                    * index
                    + rotation_index
                )
                % EXPECTED_BRIDGE_PATHS_PER_END
            ]
            for index in range(
                EXPECTED_BRIDGE_PATHS_PER_END
            )
        ]

        selected_seed_ids = {
            row[
                "node_id"
            ]
            for row in selected_seed_rows
        }

        selected_outer_ids = {
            row[
                "node_id"
            ]
            for row in mapped_outer_rows
        }

        if len(selected_seed_ids) != 15:
            raise RuntimeError(
                f"{end}: duplicate selected seed sites."
            )

        if len(selected_outer_ids) != 15:
            raise RuntimeError(
                f"{end}: duplicate selected annulus sites."
            )

        selected_seed_ids_by_end[
            end
        ] = selected_seed_ids

        selected_outer_ids_by_end[
            end
        ] = selected_outer_ids

        for bridge_index, (
            seed_row,
            outer_row,
        ) in enumerate(
            zip(
                selected_seed_rows,
                mapped_outer_rows,
            )
        ):
            seed_id = seed_row[
                "node_id"
            ]

            annulus_id = outer_row[
                "node_id"
            ]

            path_id = (
                f"{end}:BRIDGE:"
                f"{bridge_index:02d}"
            )

            bridge_ids = []

            for position_index, element in enumerate(
                bridge_sequence,
                start=1,
            ):
                bridge_id = (
                    f"BR:{end}:"
                    f"{bridge_index:02d}:"
                    f"{position_index}"
                )

                bridge_ids.append(
                    bridge_id
                )

                bridge_node_ids_by_end[
                    end
                ].append(
                    bridge_id
                )

                add_node(
                    node_rows,
                    node_elements,
                    node_types,
                    adjacency,
                    node_id=bridge_id,
                    element=element,
                    node_type=(
                        "ALTERNATING_BN_TRIMER_BRIDGE"
                    ),
                    end=end,
                    metadata={
                        "bridge_path_id": (
                            path_id
                        ),
                        "bridge_index": (
                            bridge_index
                        ),
                        "bridge_position": (
                            position_index
                        ),
                        "attached_seed_node": (
                            seed_id
                        ),
                        "attached_annulus_node": (
                            annulus_id
                        ),
                        "selected_gap_nm": (
                            selected_gap_nm
                        ),
                        "feasible_gap_minimum_nm": (
                            feasible_gap_minimum_nm
                        ),
                        "feasible_gap_maximum_nm": (
                            feasible_gap_maximum_nm
                        ),
                        "mapping_seed_parity": (
                            seed_parity
                        ),
                        "mapping_orientation": (
                            orientation
                        ),
                        "mapping_rotation": (
                            rotation_index
                        ),
                        "mapping_chirality": (
                            chirality
                        ),
                        "graph_source": (
                            "GATE3H_SELECTED_TRIMER_BRIDGE"
                        ),
                    },
                )

            path_nodes = [
                seed_id,
                *bridge_ids,
                annulus_id,
            ]

            path_edge_pairs = set()

            for path_edge_index in range(
                len(path_nodes) - 1
            ):
                first = path_nodes[
                    path_edge_index
                ]

                second = path_nodes[
                    path_edge_index + 1
                ]

                pair = tuple(
                    sorted(
                        (
                            first,
                            second,
                        )
                    )
                )

                path_edge_pairs.add(pair)

                add_edge(
                    edge_rows,
                    edge_pairs,
                    adjacency,
                    node_elements,
                    source=first,
                    target=second,
                    edge_type=(
                        "ALTERNATING_BN_TRIMER_BRIDGE"
                    ),
                    end=end,
                    heavy_atom_edge=True,
                    metadata={
                        "bridge_path_id": (
                            path_id
                        ),
                        "bridge_edge_position": (
                            path_edge_index + 1
                        ),
                        "graph_source": (
                            "GATE3H_SELECTED_TRIMER_BRIDGE"
                        ),
                    },
                )

            bridge_path_edge_sets[
                path_id
            ] = path_edge_pairs

            bridge_path_rows.append(
                {
                    "bridge_path_id": path_id,
                    "end": end,
                    "bridge_index": bridge_index,
                    "seed_node": seed_id,
                    "seed_element": (
                        node_elements[
                            seed_id
                        ]
                    ),
                    "bridge_node_1": (
                        bridge_ids[0]
                    ),
                    "bridge_element_1": (
                        node_elements[
                            bridge_ids[0]
                        ]
                    ),
                    "bridge_node_2": (
                        bridge_ids[1]
                    ),
                    "bridge_element_2": (
                        node_elements[
                            bridge_ids[1]
                        ]
                    ),
                    "bridge_node_3": (
                        bridge_ids[2]
                    ),
                    "bridge_element_3": (
                        node_elements[
                            bridge_ids[2]
                        ]
                    ),
                    "annulus_node": (
                        annulus_id
                    ),
                    "annulus_element": (
                        node_elements[
                            annulus_id
                        ]
                    ),
                    "heavy_edges_in_path": (
                        EXPECTED_BRIDGE_EDGES_PER_PATH
                    ),
                    "selected_gap_nm": (
                        selected_gap_nm
                    ),
                    "feasible_gap_minimum_nm": (
                        feasible_gap_minimum_nm
                    ),
                    "feasible_gap_maximum_nm": (
                        feasible_gap_maximum_nm
                    ),
                    "mapping_seed_parity": (
                        seed_parity
                    ),
                    "mapping_orientation": (
                        orientation
                    ),
                    "mapping_rotation": (
                        rotation_index
                    ),
                    "mapping_chirality": (
                        chirality
                    ),
                }
            )

        end_preliminary_rows.append(
            {
                "end": end,
                "seed_element": (
                    seed_element
                ),
                "bridge_sequence": (
                    "-".join(
                        bridge_sequence
                    )
                ),
                "annulus_endpoint_element": (
                    required_annulus_element
                ),
                "selected_seed_sites": (
                    len(
                        selected_seed_ids
                    )
                ),
                "selected_annulus_sites": (
                    len(
                        selected_outer_ids
                    )
                ),
                "bridge_paths": (
                    EXPECTED_BRIDGE_PATHS_PER_END
                ),
                "bridge_atoms": (
                    len(
                        bridge_node_ids_by_end[
                            end
                        ]
                    )
                ),
                "bridge_heavy_edges": (
                    EXPECTED_BRIDGE_EDGES_PER_END
                ),
                "selected_gap_nm": (
                    selected_gap_nm
                ),
                "feasible_gap_minimum_nm": (
                    feasible_gap_minimum_nm
                ),
                "feasible_gap_maximum_nm": (
                    feasible_gap_maximum_nm
                ),
                "mapping_seed_parity": (
                    seed_parity
                ),
                "mapping_orientation": (
                    orientation
                ),
                "mapping_rotation": (
                    rotation_index
                ),
                "mapping_chirality": (
                    chirality
                ),
            }
        )

    hydrogen_counts_by_end_role: dict[
        str,
        Counter[str]
    ] = {
        "LOWER": Counter(),
        "UPPER": Counter(),
    }

    hydrogen_global_index = 0

    for end in (
        "LOWER",
        "UPPER",
    ):
        seed_ids = sorted(
            node_id
            for node_id in adjacency
            if (
                node_types[
                    node_id
                ]
                == "HEXAGONAL_EDGE_COMPLETION_SEED"
                and original_nodes_by_id[
                    node_id
                ][
                    "end"
                ]
                == end
            )
        )

        outer_ids = sorted(
            node_id
            for node_id in adjacency
            if (
                node_types[
                    node_id
                ]
                == "ANNULUS_OUTER_BOUNDARY"
                and original_nodes_by_id[
                    node_id
                ][
                    "end"
                ]
                == end
            )
        )

        inner_ids = sorted(
            node_id
            for node_id in adjacency
            if (
                node_types[
                    node_id
                ]
                == "ANNULUS_INNER_BOUNDARY"
                and original_nodes_by_id[
                    node_id
                ][
                    "end"
                ]
                == end
            )
        )

        unselected_seed_ids = [
            node_id
            for node_id in seed_ids
            if node_id
            not in selected_seed_ids_by_end[
                end
            ]
        ]

        unselected_outer_ids = [
            node_id
            for node_id in outer_ids
            if node_id
            not in selected_outer_ids_by_end[
                end
            ]
        ]

        if len(unselected_seed_ids) != 15:
            raise RuntimeError(
                f"{end}: unselected seed count "
                f"{len(unselected_seed_ids)}/15"
            )

        if len(unselected_outer_ids) != 15:
            raise RuntimeError(
                f"{end}: unselected outer count "
                f"{len(unselected_outer_ids)}/15"
            )

        if len(inner_ids) != 12:
            raise RuntimeError(
                f"{end}: inner boundary count "
                f"{len(inner_ids)}/12"
            )

        for node_id in unselected_seed_ids:
            add_hydrogen(
                node_rows=node_rows,
                edge_rows=edge_rows,
                edge_pairs=edge_pairs,
                adjacency=adjacency,
                node_elements=node_elements,
                node_types=node_types,
                heavy_node=node_id,
                end=end,
                hydrogen_role=(
                    "SEED_PASSIVANT_H"
                ),
                hydrogen_index=(
                    hydrogen_global_index
                ),
            )

            hydrogen_global_index += 1

            hydrogen_counts_by_end_role[
                end
            ][
                "SEED_PASSIVANT_H"
            ] += 1

        for node_id in unselected_outer_ids:
            add_hydrogen(
                node_rows=node_rows,
                edge_rows=edge_rows,
                edge_pairs=edge_pairs,
                adjacency=adjacency,
                node_elements=node_elements,
                node_types=node_types,
                heavy_node=node_id,
                end=end,
                hydrogen_role=(
                    "ANNULUS_OUTER_PASSIVANT_H"
                ),
                hydrogen_index=(
                    hydrogen_global_index
                ),
            )

            hydrogen_global_index += 1

            hydrogen_counts_by_end_role[
                end
            ][
                "ANNULUS_OUTER_PASSIVANT_H"
            ] += 1

        for node_id in inner_ids:
            add_hydrogen(
                node_rows=node_rows,
                edge_rows=edge_rows,
                edge_pairs=edge_pairs,
                adjacency=adjacency,
                node_elements=node_elements,
                node_types=node_types,
                heavy_node=node_id,
                end=end,
                hydrogen_role=(
                    "ANNULUS_INNER_PASSIVANT_H"
                ),
                hydrogen_index=(
                    hydrogen_global_index
                ),
            )

            hydrogen_global_index += 1

            hydrogen_counts_by_end_role[
                end
            ][
                "ANNULUS_INNER_PASSIVANT_H"
            ] += 1

        for node_id in bridge_node_ids_by_end[
            end
        ]:
            add_hydrogen(
                node_rows=node_rows,
                edge_rows=edge_rows,
                edge_pairs=edge_pairs,
                adjacency=adjacency,
                node_elements=node_elements,
                node_types=node_types,
                heavy_node=node_id,
                end=end,
                hydrogen_role=(
                    "BRIDGE_PASSIVANT_H"
                ),
                hydrogen_index=(
                    hydrogen_global_index
                ),
            )

            hydrogen_global_index += 1

            hydrogen_counts_by_end_role[
                end
            ][
                "BRIDGE_PASSIVANT_H"
            ] += 1

    heavy_nodes = {
        node_id
        for node_id, element
        in node_elements.items()
        if element != "H"
    }

    hydrogen_nodes = {
        node_id
        for node_id, element
        in node_elements.items()
        if element == "H"
    }

    heavy_adjacency = {
        node_id: {
            neighbor
            for neighbor in adjacency[
                node_id
            ]
            if neighbor in heavy_nodes
        }
        for node_id in heavy_nodes
    }

    full_components = connected_components(
        adjacency
    )

    heavy_components = connected_components(
        heavy_adjacency
    )

    heavy_bipartite, heavy_colors = (
        bipartite_coloring(
            heavy_adjacency
        )
    )

    heavy_four_cycles = count_four_cycles(
        heavy_adjacency
    )

    heavy_contains_six_cycle = (
        contains_six_cycle(
            heavy_adjacency
        )
    )

    if heavy_four_cycles > 0:
        heavy_girth = 4
    elif (
        heavy_bipartite
        and heavy_contains_six_cycle
    ):
        heavy_girth = 6
    else:
        heavy_girth = 0

    heavy_degree_failures = [
        node_id
        for node_id in heavy_nodes
        if len(
            adjacency[
                node_id
            ]
        )
        != 3
    ]

    hydrogen_degree_failures = [
        node_id
        for node_id in hydrogen_nodes
        if len(
            adjacency[
                node_id
            ]
        )
        != 1
    ]

    bridge_heavy_degree_failures = [
        node_id
        for node_id
        in (
            bridge_node_ids_by_end[
                "LOWER"
            ]
            + bridge_node_ids_by_end[
                "UPPER"
            ]
        )
        if len(
            heavy_adjacency[
                node_id
            ]
        )
        != 2
    ]

    nonheteropolar_heavy_edges = [
        row
        for row in edge_rows
        if parse_bool(
            row[
                "heavy_atom_edge"
            ]
        )
        and not parse_bool(
            row[
                "heteropolar_BN_edge"
            ]
        )
    ]

    direct_seed_annulus_edges_remaining = [
        row
        for row in edge_rows
        if row[
            "edge_type"
        ]
        == "PARTIAL_HETEROPOLAR_SEED_TO_ANNULUS"
    ]

    bridge_graph_edges = [
        row
        for row in edge_rows
        if row[
            "edge_type"
        ]
        == "ALTERNATING_BN_TRIMER_BRIDGE"
    ]

    bridge_path_cycles = []

    for row in bridge_path_rows:
        path_id = row[
            "bridge_path_id"
        ]

        excluded_edges = (
            bridge_path_edge_sets[
                path_id
            ]
        )

        alternative_length = (
            shortest_path_length(
                heavy_adjacency,
                row[
                    "seed_node"
                ],
                row[
                    "annulus_node"
                ],
                excluded_edges,
            )
        )

        cycle_length = (
            None
            if alternative_length is None
            else (
                alternative_length
                + EXPECTED_BRIDGE_EDGES_PER_PATH
            )
        )

        row[
            "alternative_heavy_path_length"
        ] = (
            ""
            if alternative_length is None
            else alternative_length
        )

        row[
            "shortest_cycle_containing_bridge_path"
        ] = (
            ""
            if cycle_length is None
            else cycle_length
        )

        bridge_path_cycles.append(
            cycle_length
        )

    disconnected_bridge_paths = sum(
        value is None
        for value in bridge_path_cycles
    )

    finite_bridge_cycles = [
        int(value)
        for value in bridge_path_cycles
        if value is not None
    ]

    if not finite_bridge_cycles:
        raise RuntimeError(
            "No bridge-containing cycle lengths "
            "could be resolved."
        )

    bridge_cycle_minimum = min(
        finite_bridge_cycles
    )

    bridge_cycle_maximum = max(
        finite_bridge_cycles
    )

    composition_by_end = {
        "LOWER": Counter(),
        "UPPER": Counter(),
    }

    bridge_composition_by_end = {
        "LOWER": Counter(),
        "UPPER": Counter(),
    }

    node_row_by_id = {
        row[
            "node_id"
        ]: row
        for row in node_rows
    }

    for node_id, element in node_elements.items():
        row = node_row_by_id[
            node_id
        ]

        end = row[
            "end"
        ]

        if end in composition_by_end:
            composition_by_end[
                end
            ][
                element
            ] += 1

        if (
            row[
                "node_type"
            ]
            == "ALTERNATING_BN_TRIMER_BRIDGE"
        ):
            bridge_composition_by_end[
                end
            ][
                element
            ] += 1

    end_summary_rows = []

    for preliminary in end_preliminary_rows:
        end = preliminary[
            "end"
        ]

        role_counts = (
            hydrogen_counts_by_end_role[
                end
            ]
        )

        end_bridge_rows = [
            row
            for row in bridge_path_rows
            if row[
                "end"
            ]
            == end
        ]

        end_cycle_lengths = [
            int(
                row[
                    "shortest_cycle_containing_bridge_path"
                ]
            )
            for row in end_bridge_rows
            if row[
                "shortest_cycle_containing_bridge_path"
            ]
            != ""
        ]

        end_summary_rows.append(
            {
                **preliminary,
                "bridge_B_atoms": (
                    bridge_composition_by_end[
                        end
                    ][
                        "B"
                    ]
                ),
                "bridge_N_atoms": (
                    bridge_composition_by_end[
                        end
                    ][
                        "N"
                    ]
                ),
                "seed_H_passivants": (
                    role_counts[
                        "SEED_PASSIVANT_H"
                    ]
                ),
                "outer_H_passivants": (
                    role_counts[
                        "ANNULUS_OUTER_PASSIVANT_H"
                    ]
                ),
                "inner_H_passivants": (
                    role_counts[
                        "ANNULUS_INNER_PASSIVANT_H"
                    ]
                ),
                "bridge_H_passivants": (
                    role_counts[
                        "BRIDGE_PASSIVANT_H"
                    ]
                ),
                "total_H_passivants": sum(
                    role_counts.values()
                ),
                "total_added_heavy_atoms": (
                    EXPECTED_HEAVY_ATOMS_PER_END
                ),
                "shortest_bridge_containing_cycle": (
                    min(
                        end_cycle_lengths
                    )
                ),
                "longest_bridge_containing_cycle": (
                    max(
                        end_cycle_lengths
                    )
                ),
            }
        )

    combined_bridge_composition = (
        bridge_composition_by_end[
            "LOWER"
        ]
        + bridge_composition_by_end[
            "UPPER"
        ]
    )

    node_type_counts = Counter(
        node_types.values()
    )

    heavy_edge_count = sum(
        parse_bool(
            row[
                "heavy_atom_edge"
            ]
        )
        for row in edge_rows
    )

    hydrogen_edge_count = (
        len(edge_rows)
        - heavy_edge_count
    )

    gates = {
        "Gate3F_graph_design_is_accepted": (
            design_summary.get(
                "decision"
            )
            == EXPECTED_DESIGN_DECISION
        ),
        "Gate3G1_direct_junction_is_rejected": (
            direct_summary.get(
                "decision"
            )
            == EXPECTED_DIRECT_REJECTION
        ),
        "Gate3H_trimer_bridge_class_is_selected": (
            bridge_summary.get(
                "decision"
            )
            == EXPECTED_BRIDGE_DECISION
            and selected_bridge_atoms
            == EXPECTED_BRIDGE_ATOMS_PER_PATH
        ),
        "30_rejected_direct_seed_annulus_edges_were_removed": (
            direct_edges_removed
            == EXPECTED_BRIDGE_PATHS_TOTAL
        ),
        "no_rejected_direct_seed_annulus_edges_remain": (
            len(
                direct_seed_annulus_edges_remaining
            )
            == 0
        ),
        "15_trimer_bridge_paths_were_built_per_end": (
            all(
                int(
                    row[
                        "bridge_paths"
                    ]
                )
                == EXPECTED_BRIDGE_PATHS_PER_END
                for row in end_summary_rows
            )
        ),
        "30_trimer_bridge_paths_were_built_total": (
            len(
                bridge_path_rows
            )
            == EXPECTED_BRIDGE_PATHS_TOTAL
        ),
        "45_bridge_atoms_were_added_per_end": (
            all(
                int(
                    row[
                        "bridge_atoms"
                    ]
                )
                == EXPECTED_BRIDGE_ATOMS_PER_END
                for row in end_summary_rows
            )
        ),
        "90_bridge_atoms_were_added_total": (
            node_type_counts[
                "ALTERNATING_BN_TRIMER_BRIDGE"
            ]
            == EXPECTED_BRIDGE_ATOMS_TOTAL
        ),
        "60_bridge_heavy_edges_were_added_per_end": (
            all(
                int(
                    row[
                        "bridge_heavy_edges"
                    ]
                )
                == EXPECTED_BRIDGE_EDGES_PER_END
                for row in end_summary_rows
            )
        ),
        "120_bridge_heavy_edges_were_added_total": (
            len(
                bridge_graph_edges
            )
            == EXPECTED_BRIDGE_EDGES_TOTAL
        ),
        "lower_bridge_sequence_is_B_N_B": (
            next(
                row
                for row in end_summary_rows
                if row[
                    "end"
                ]
                == "LOWER"
            )[
                "bridge_sequence"
            ]
            == "B-N-B"
        ),
        "upper_bridge_sequence_is_N_B_N": (
            next(
                row
                for row in end_summary_rows
                if row[
                    "end"
                ]
                == "UPPER"
            )[
                "bridge_sequence"
            ]
            == "N-B-N"
        ),
        "combined_bridge_composition_is_45B_45N": (
            combined_bridge_composition[
                "B"
            ]
            == 45
            and combined_bridge_composition[
                "N"
            ]
            == 45
        ),
        "all_heavy_atoms_have_total_coordination3": (
            len(
                heavy_degree_failures
            )
            == 0
        ),
        "all_bridge_atoms_have_two_heavy_neighbors": (
            len(
                bridge_heavy_degree_failures
            )
            == 0
        ),
        "all_H_atoms_have_coordination1": (
            len(
                hydrogen_degree_failures
            )
            == 0
        ),
        "all_heavy_edges_are_heteropolar_BN": (
            len(
                nonheteropolar_heavy_edges
            )
            == 0
        ),
        "full_graph_is_connected": (
            len(
                full_components
            )
            == 1
        ),
        "heavy_graph_is_connected": (
            len(
                heavy_components
            )
            == 1
        ),
        "heavy_graph_is_bipartite": (
            heavy_bipartite
        ),
        "heavy_graph_contains_no_four_member_cycles": (
            heavy_four_cycles == 0
        ),
        "heavy_graph_girth_is_at_least6": (
            heavy_girth
            >= MINIMUM_ACCEPTED_HEAVY_GIRTH
        ),
        "every_bridge_path_participates_in_a_cycle": (
            disconnected_bridge_paths == 0
        ),
        "every_bridge_containing_cycle_has_length_at_least6": (
            bridge_cycle_minimum
            >= MINIMUM_BRIDGE_CONTAINING_CYCLE_LENGTH
        ),
        "87_H_passivants_were_added_per_end": (
            all(
                int(
                    row[
                        "total_H_passivants"
                    ]
                )
                == EXPECTED_H_PER_END
                for row in end_summary_rows
            )
        ),
        "174_H_passivants_were_added_total": (
            len(
                hydrogen_nodes
            )
            == EXPECTED_H_TOTAL
        ),
        "passivation_partition_is_15_15_12_45_per_end": (
            all(
                int(
                    row[
                        "seed_H_passivants"
                    ]
                )
                == EXPECTED_SEED_H_PER_END
                and int(
                    row[
                        "outer_H_passivants"
                    ]
                )
                == EXPECTED_OUTER_H_PER_END
                and int(
                    row[
                        "inner_H_passivants"
                    ]
                )
                == EXPECTED_INNER_H_PER_END
                and int(
                    row[
                        "bridge_H_passivants"
                    ]
                )
                == EXPECTED_BRIDGE_H_PER_END
                for row in end_summary_rows
            )
        ),
        "201_added_heavy_atoms_are_present_per_end": (
            all(
                int(
                    row[
                        "total_added_heavy_atoms"
                    ]
                )
                == EXPECTED_HEAVY_ATOMS_PER_END
                for row in end_summary_rows
            )
        ),
        "2082_total_heavy_atoms_are_present": (
            len(
                heavy_nodes
            )
            == EXPECTED_TOTAL_HEAVY_ATOMS
        ),
        "2256_total_nodes_are_present": (
            len(
                adjacency
            )
            == EXPECTED_TOTAL_NODES
        ),
        "no_coordinates_were_assigned": all(
            not parse_bool(
                row.get(
                    "coordinates_assigned",
                    "false",
                )
            )
            for row in node_rows
        )
        and all(
            not parse_bool(
                row.get(
                    "coordinates_assigned",
                    "false",
                )
            )
            for row in edge_rows
        ),
        "no_formal_charges_were_assigned": all(
            not parse_bool(
                row.get(
                    "formal_charge_assigned",
                    "false",
                )
            )
            for row in node_rows
        ),
        "no_force_field_types_were_assigned": all(
            not parse_bool(
                row.get(
                    "force_field_type_assigned",
                    "false",
                )
            )
            for row in node_rows
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
        "R2_ALTERNATING_BN_TRIMER_BRIDGE_GRAPH_REQUIRES_REVIEW"
    )

    required_next_step = (
        "BUILD_AND_VALIDATE_R2_ALTERNATING_BN_TRIMER_BRIDGE_"
        "STATIC_COORDINATE_EMBEDDING"
        if accepted
        else
        "REVIEW_R2_ALTERNATING_BN_TRIMER_BRIDGE_GRAPH_FAILURES"
    )

    summary = {
        "decision": decision,
        "direct_seed_annulus_edges_removed": (
            direct_edges_removed
        ),
        "direct_seed_annulus_edges_remaining": (
            len(
                direct_seed_annulus_edges_remaining
            )
        ),
        "bridge_paths_lower": (
            sum(
                row[
                    "end"
                ]
                == "LOWER"
                for row in bridge_path_rows
            )
        ),
        "bridge_paths_upper": (
            sum(
                row[
                    "end"
                ]
                == "UPPER"
                for row in bridge_path_rows
            )
        ),
        "bridge_paths_total": (
            len(
                bridge_path_rows
            )
        ),
        "bridge_atoms_lower": (
            len(
                bridge_node_ids_by_end[
                    "LOWER"
                ]
            )
        ),
        "bridge_atoms_upper": (
            len(
                bridge_node_ids_by_end[
                    "UPPER"
                ]
            )
        ),
        "bridge_atoms_total": (
            node_type_counts[
                "ALTERNATING_BN_TRIMER_BRIDGE"
            ]
        ),
        "bridge_heavy_edges_total": (
            len(
                bridge_graph_edges
            )
        ),
        "bridge_B_atoms_total": (
            combined_bridge_composition[
                "B"
            ]
        ),
        "bridge_N_atoms_total": (
            combined_bridge_composition[
                "N"
            ]
        ),
        "total_heavy_atoms": (
            len(
                heavy_nodes
            )
        ),
        "total_H_atoms": (
            len(
                hydrogen_nodes
            )
        ),
        "total_nodes": (
            len(
                adjacency
            )
        ),
        "total_heavy_edges": (
            heavy_edge_count
        ),
        "total_H_edges": (
            hydrogen_edge_count
        ),
        "heavy_degree_failures": (
            len(
                heavy_degree_failures
            )
        ),
        "bridge_heavy_degree_failures": (
            len(
                bridge_heavy_degree_failures
            )
        ),
        "H_degree_failures": (
            len(
                hydrogen_degree_failures
            )
        ),
        "nonheteropolar_heavy_edges": (
            len(
                nonheteropolar_heavy_edges
            )
        ),
        "full_graph_components": (
            len(
                full_components
            )
        ),
        "heavy_graph_components": (
            len(
                heavy_components
            )
        ),
        "heavy_graph_bipartite": (
            heavy_bipartite
        ),
        "heavy_graph_four_member_cycles": (
            heavy_four_cycles
        ),
        "heavy_graph_girth": (
            heavy_girth
        ),
        "shortest_bridge_containing_cycle": (
            bridge_cycle_minimum
        ),
        "longest_bridge_containing_cycle": (
            bridge_cycle_maximum
        ),
        "bridge_paths_without_alternative_cycle": (
            disconnected_bridge_paths
        ),
        "candidate_is_final_chemistry": False,
        "static_coordinate_embedding_authorized": (
            accepted
        ),
        "molecular_topology_generation_authorized": False,
        "formal_charge_assignment_authorized": False,
        "force_field_parameterization_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "mobile_MD_authorized": False,
        "multitemperature_MD_authorized": False,
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
        GRAPH_NODES_CSV,
        node_rows,
    )

    write_csv(
        GRAPH_EDGES_CSV,
        edge_rows,
    )

    write_csv(
        BRIDGE_PATHS_CSV,
        bridge_path_rows,
    )

    write_csv(
        END_SUMMARY_CSV,
        end_summary_rows,
    )

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

    GRAPH_JSON.write_text(
        json.dumps(
            {
                "summary": summary,
                "end_summaries": (
                    end_summary_rows
                ),
                "gates": gates,
                "selected_mappings": (
                    selected_mappings
                ),
                "limitations": [
                    (
                        "The graph assigns connectivity and "
                        "coordination only. It contains no "
                        "three-dimensional bridge conformers."
                    ),
                    (
                        "Each bridge atom is provisionally "
                        "monohydrogenated to complete coordination."
                    ),
                    (
                        "The selected 0.1771 nm axial gap comes "
                        "from a conformational-envelope screen, "
                        "not an energy calculation."
                    ),
                    (
                        "Graph validity does not establish "
                        "energetic stability or synthetic feasibility."
                    ),
                    (
                        "No molecular topology, formal charges, "
                        "force-field parameters, minimization, "
                        "MD, or QM calculation was generated."
                    ),
                ],
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
                "Gate3F_design_nodes"
            ),
            "file": relative(
                DESIGN_NODES_CSV
            ),
            "sha256": sha256(
                DESIGN_NODES_CSV
            ),
        },
        {
            "role": (
                "Gate3F_design_edges"
            ),
            "file": relative(
                DESIGN_EDGES_CSV
            ),
            "sha256": sha256(
                DESIGN_EDGES_CSV
            ),
        },
        {
            "role": (
                "Gate3F_design_summary"
            ),
            "file": relative(
                DESIGN_SUMMARY_CSV
            ),
            "sha256": sha256(
                DESIGN_SUMMARY_CSV
            ),
        },
        {
            "role": (
                "Gate3G1_direct_lower_bound_summary"
            ),
            "file": relative(
                DIRECT_LOWER_BOUND_SUMMARY_CSV
            ),
            "sha256": sha256(
                DIRECT_LOWER_BOUND_SUMMARY_CSV
            ),
        },
        {
            "role": (
                "Gate3H_bridge_feasibility_summary"
            ),
            "file": relative(
                BRIDGE_FEASIBILITY_SUMMARY_CSV
            ),
            "sha256": sha256(
                BRIDGE_FEASIBILITY_SUMMARY_CSV
            ),
        },
        {
            "role": (
                "Gate3H_selected_bridge_candidate"
            ),
            "file": relative(
                SELECTED_BRIDGE_CANDIDATE_CSV
            ),
            "sha256": sha256(
                SELECTED_BRIDGE_CANDIDATE_CSV
            ),
        },
    ]

    write_csv(
        SOURCE_MANIFEST_CSV,
        source_rows,
    )

    end_lines = "\n".join(
        (
            f"### {row['end']}\n\n"
            f"- Bridge sequence: "
            f"**{row['bridge_sequence']}**\n"
            f"- Bridge paths/atoms/heavy edges: "
            f"**{row['bridge_paths']}/"
            f"{row['bridge_atoms']}/"
            f"{row['bridge_heavy_edges']}**\n"
            f"- Bridge B/N atoms: "
            f"**{row['bridge_B_atoms']}/"
            f"{row['bridge_N_atoms']}**\n"
            f"- H partition seed/outer/inner/bridge: "
            f"**{row['seed_H_passivants']}/"
            f"{row['outer_H_passivants']}/"
            f"{row['inner_H_passivants']}/"
            f"{row['bridge_H_passivants']}**\n"
            f"- Total H: "
            f"**{row['total_H_passivants']}**\n"
            f"- Selected gap: "
            f"**{float(row['selected_gap_nm']):.6f} nm**\n"
            f"- Shortest/longest bridge-containing cycle: "
            f"**{row['shortest_bridge_containing_cycle']}/"
            f"{row['longest_bridge_containing_cycle']}**"
        )
        for row in end_summary_rows
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
        f"""# R2 Alternating BN Trimer-Bridge Graph

## Scope

This gate replaces the rejected direct seed–annulus bonds with
alternating three-atom BN bridge paths.

No coordinates, molecular topology, formal charges, force-field
parameters, minimization, MD, or QM calculation were generated.

## Graph transformation

- Rejected direct seed–annulus edges removed:
  **{direct_edges_removed}**
- Direct seed–annulus edges remaining:
  **{len(direct_seed_annulus_edges_remaining)}**
- Bridge paths:
  **{len(bridge_path_rows)}**
- Bridge atoms:
  **{node_type_counts['ALTERNATING_BN_TRIMER_BRIDGE']}**
- Bridge heavy edges:
  **{len(bridge_graph_edges)}**

{end_lines}

## Combined graph

- Total heavy atoms:
  **{len(heavy_nodes)}**
- Total H atoms:
  **{len(hydrogen_nodes)}**
- Total nodes:
  **{len(adjacency)}**
- Heavy/H edges:
  **{heavy_edge_count}/{hydrogen_edge_count}**
- Heavy-degree failures:
  **{len(heavy_degree_failures)}**
- Bridge heavy-degree failures:
  **{len(bridge_heavy_degree_failures)}**
- H-degree failures:
  **{len(hydrogen_degree_failures)}**
- Nonheteropolar heavy edges:
  **{len(nonheteropolar_heavy_edges)}**
- Heavy connected components:
  **{len(heavy_components)}**
- Bipartite:
  **{heavy_bipartite}**
- Four-member heavy cycles:
  **{heavy_four_cycles}**
- Heavy-graph girth:
  **{heavy_girth}**

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Candidate is final chemistry:
  **NO**
- Static coordinate embedding authorized:
  **{'YES' if accepted else 'NO'}**
- Molecular topology generation authorized:
  **NO**
- Formal charge assignment authorized:
  **NO**
- Force-field parameterization authorized:
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

This gate establishes graph-level valence completion for the selected
three-atom bridge class. It does not prove that 30 simultaneous bridge
conformers can be embedded without steric clashes, unacceptable bond
angles or excessive strain.

The next gate must construct explicit bridge conformers for both ends
while preserving the parent scaffold, annulus aperture, H passivation,
and the selected axial separation.
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 alternating BN trimer-bridge "
        "graph completed."
    )

    print(
        "Rejected direct seed-annulus edges "
        "removed / remaining: "
        f"{direct_edges_removed}/"
        f"{len(direct_seed_annulus_edges_remaining)}"
    )

    print(
        "Bridge paths lower/upper/total: "
        f"{sum(row['end'] == 'LOWER' for row in bridge_path_rows)}/"
        f"{sum(row['end'] == 'UPPER' for row in bridge_path_rows)}/"
        f"{len(bridge_path_rows)}"
    )

    print(
        "Bridge atoms lower/upper/total: "
        f"{len(bridge_node_ids_by_end['LOWER'])}/"
        f"{len(bridge_node_ids_by_end['UPPER'])}/"
        f"{node_type_counts['ALTERNATING_BN_TRIMER_BRIDGE']}"
    )

    print(
        "Bridge heavy edges total: "
        f"{len(bridge_graph_edges)}"
    )

    for row in end_summary_rows:
        print(
            f"{row['end']} sequence / bridge B/N / "
            "H seed/outer/inner/bridge/total: "
            f"{row['bridge_sequence']}/"
            f"{row['bridge_B_atoms']}/"
            f"{row['bridge_N_atoms']}/"
            f"{row['seed_H_passivants']}/"
            f"{row['outer_H_passivants']}/"
            f"{row['inner_H_passivants']}/"
            f"{row['bridge_H_passivants']}/"
            f"{row['total_H_passivants']}"
        )

        print(
            f"{row['end']} gap range / selected / "
            "bridge-cycle min-max: "
            f"{float(row['feasible_gap_minimum_nm']):.6f}-"
            f"{float(row['feasible_gap_maximum_nm']):.6f}/"
            f"{float(row['selected_gap_nm']):.6f}/"
            f"{row['shortest_bridge_containing_cycle']}-"
            f"{row['longest_bridge_containing_cycle']}"
        )

    print(
        "Total heavy / H / all nodes: "
        f"{len(heavy_nodes)}/"
        f"{len(hydrogen_nodes)}/"
        f"{len(adjacency)}"
    )

    print(
        "Heavy / H edges: "
        f"{heavy_edge_count}/"
        f"{hydrogen_edge_count}"
    )

    print(
        "Heavy / bridge-heavy / H degree failures: "
        f"{len(heavy_degree_failures)}/"
        f"{len(bridge_heavy_degree_failures)}/"
        f"{len(hydrogen_degree_failures)}"
    )

    print(
        "Nonheteropolar heavy edges: "
        f"{len(nonheteropolar_heavy_edges)}"
    )

    print(
        "Full/heavy components / bipartite / "
        "girth / four-cycles: "
        f"{len(full_components)}/"
        f"{len(heavy_components)}/"
        f"{heavy_bipartite}/"
        f"{heavy_girth}/"
        f"{heavy_four_cycles}"
    )

    print(
        "Bridge paths without cycle / "
        "bridge-cycle min-max: "
        f"{disconnected_bridge_paths}/"
        f"{bridge_cycle_minimum}/"
        f"{bridge_cycle_maximum}"
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
        "Candidate is final chemistry: NO"
    )

    print(
        "Static coordinate embedding authorized: "
        f"{'YES' if accepted else 'NO'}"
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
        GRAPH_NODES_CSV,
        GRAPH_EDGES_CSV,
        BRIDGE_PATHS_CSV,
        END_SUMMARY_CSV,
        SUMMARY_CSV,
        GATES_CSV,
        GRAPH_JSON,
        SOURCE_MANIFEST_CSV,
        REPORT_MD,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not accepted:
        raise RuntimeError(
            "R2 alternating BN trimer-bridge "
            "graph requires review."
        )


if __name__ == "__main__":
    main()
