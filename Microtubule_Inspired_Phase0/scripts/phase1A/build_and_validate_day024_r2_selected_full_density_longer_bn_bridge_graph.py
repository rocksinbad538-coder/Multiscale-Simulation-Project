#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, deque
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

BASE = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design"
)

GATE3I = (
    BASE
    / "11_r2_alternating_bn_trimer_bridge_graph"
)

GATE3L = (
    BASE
    / "15_r2_full_density_longer_bn_bridge_screen"
)

OUTPUT = (
    BASE
    / "16_r2_selected_full_density_longer_bn_bridge_graph"
)

SOURCE_NODES = (
    GATE3I
    / "r2_alternating_bn_trimer_bridge_graph_nodes.csv"
)

SOURCE_EDGES = (
    GATE3I
    / "r2_alternating_bn_trimer_bridge_graph_edges.csv"
)

SOURCE_GRAPH_SUMMARY = (
    GATE3I
    / "r2_alternating_bn_trimer_bridge_graph_summary.csv"
)

SELECTED_CANDIDATE = (
    GATE3L
    / "r2_full_density_longer_bridge_selected_candidate.csv"
)

SCREEN_SUMMARY = (
    GATE3L
    / "r2_full_density_longer_bridge_screen_summary.csv"
)

GRAPH_NODES = (
    OUTPUT
    / "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

GRAPH_EDGES = (
    OUTPUT
    / "r2_selected_longer_bn_bridge_graph_edges.csv"
)

BRIDGE_PATHS = (
    OUTPUT
    / "r2_selected_longer_bn_bridge_paths.csv"
)

END_SUMMARY = (
    OUTPUT
    / "r2_selected_longer_bn_bridge_end_summary.csv"
)

SUMMARY = (
    OUTPUT
    / "r2_selected_longer_bn_bridge_graph_summary.csv"
)

GATES = (
    OUTPUT
    / "r2_selected_longer_bn_bridge_graph_gates.csv"
)

GRAPH_JSON = (
    OUTPUT
    / "r2_selected_longer_bn_bridge_graph.json"
)

MANIFEST = (
    OUTPUT
    / "r2_selected_longer_bn_bridge_source_manifest.csv"
)

REPORT = (
    OUTPUT
    / "R2_SELECTED_FULL_DENSITY_LONGER_BN_BRIDGE_GRAPH_DAY024.md"
)

EXPECTED_SOURCE_GRAPH_DECISION = (
    "R2_ALTERNATING_BN_TRIMER_BRIDGE_GRAPH_VALIDATED"
)

EXPECTED_SCREEN_DECISION = (
    "R2_FULL_DENSITY_LONGER_BN_BRIDGE_CLASS_IDENTIFIED"
)

PASS_DECISION = (
    "R2_SELECTED_FULL_DENSITY_FOUR_ATOM_BN_BRIDGE_GRAPH_VALIDATED"
)

REVIEW_DECISION = (
    "R2_SELECTED_FULL_DENSITY_FOUR_ATOM_BN_BRIDGE_GRAPH_REQUIRES_REVIEW"
)

EXPECTED_BRIDGE_ATOMS_PER_PATH = 4
EXPECTED_BONDS_PER_PATH = 5
EXPECTED_ATTACHMENTS_PER_END = 15
EXPECTED_PATHS_TOTAL = 30

EXPECTED_BRIDGE_ATOMS_PER_END = 60
EXPECTED_BRIDGE_ATOMS_TOTAL = 120

EXPECTED_PARENT_ATOMS = 1680
EXPECTED_SEED_ATOMS = 60
EXPECTED_ANNULUS_ATOMS = 252

EXPECTED_HEAVY_ATOMS = 2112
EXPECTED_H_ATOMS = 204
EXPECTED_TOTAL_NODES = 2316

EXPECTED_HEAVY_EDGES = 3066
EXPECTED_H_EDGES = 204

EXPECTED_CYCLE_LENGTH = 16


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
    )


def require_file(path: Path) -> None:
    if (
        not path.is_file()
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


def read_rows(
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


def read_one(
    path: Path,
) -> dict[str, str]:
    rows = read_rows(path)

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected one row in {path}; found {len(rows)}"
        )

    return rows[0]


def write_rows(
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


def split_nodes(value: str) -> list[str]:
    nodes = [
        item.strip()
        for item in value.split("|")
        if item.strip()
    ]

    if not nodes:
        raise RuntimeError(
            "Could not parse node-list field."
        )

    return nodes


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

            for neighbor in adjacency[node]:
                if neighbor not in component:
                    queue.append(neighbor)

        components.append(component)
        remaining -= component

    return components


def bipartite(
    adjacency: dict[str, set[str]],
) -> bool:
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
                    colors[neighbor] = 1 - colors[node]
                    queue.append(neighbor)

                elif colors[neighbor] == colors[node]:
                    return False

    return True


def count_four_cycles(
    adjacency: dict[str, set[str]],
) -> int:
    nodes = sorted(adjacency)
    raw = 0

    for first_index, first in enumerate(nodes):
        first_neighbors = adjacency[first]

        for second in nodes[first_index + 1:]:
            common = len(
                first_neighbors
                & adjacency[second]
            )

            if common >= 2:
                raw += (
                    common
                    * (
                        common - 1
                    )
                    // 2
                )

    if raw % 2 != 0:
        raise RuntimeError(
            "Invalid four-cycle common-neighbor count."
        )

    return raw // 2


def shortest_path_length(
    adjacency: dict[str, set[str]],
    source: str,
    target: str,
    excluded_edges: set[tuple[str, str]],
) -> int | None:
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

    visited = {source}

    while queue:
        node, distance = queue.popleft()

        for neighbor in adjacency[node]:
            edge = tuple(
                sorted(
                    (
                        node,
                        neighbor,
                    )
                )
            )

            if edge in excluded_edges:
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


def main() -> None:
    OUTPUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        SOURCE_NODES,
        SOURCE_EDGES,
        SOURCE_GRAPH_SUMMARY,
        SELECTED_CANDIDATE,
        SCREEN_SUMMARY,
    ):
        require_file(required)

    source_nodes = read_rows(
        SOURCE_NODES
    )

    source_edges = read_rows(
        SOURCE_EDGES
    )

    source_graph_summary = read_one(
        SOURCE_GRAPH_SUMMARY
    )

    candidate_rows = read_rows(
        SELECTED_CANDIDATE
    )

    screen_summary = read_one(
        SCREEN_SUMMARY
    )

    if source_graph_summary.get(
        "decision"
    ) != EXPECTED_SOURCE_GRAPH_DECISION:
        raise RuntimeError(
            "Gate 3I source graph is not accepted."
        )

    if screen_summary.get(
        "decision"
    ) != EXPECTED_SCREEN_DECISION:
        raise RuntimeError(
            "Gate 3L longer-bridge screen is not accepted."
        )

    selected_bridge_atoms = parse_int(
        screen_summary,
        "selected_bridge_atoms_per_attachment",
    )

    selected_attachments = parse_int(
        screen_summary,
        "selected_attachments_per_end",
    )

    if (
        selected_bridge_atoms
        != EXPECTED_BRIDGE_ATOMS_PER_PATH
        or selected_attachments
        != EXPECTED_ATTACHMENTS_PER_END
    ):
        raise RuntimeError(
            "Unexpected selected bridge class or attachment count."
        )

    selected_end_rows = {
        row[
            "end"
        ]: row
        for row in candidate_rows
        if row.get(
            "classification"
        )
        == "SELECTED_END_MAPPING"
    }

    if set(selected_end_rows) != {
        "LOWER",
        "UPPER",
    }:
        raise RuntimeError(
            "Could not resolve both selected end mappings."
        )

    selected_class_rows = [
        row
        for row in candidate_rows
        if row.get(
            "classification"
        )
        == "SELECTED_CLASS"
    ]

    if len(selected_class_rows) != 1:
        raise RuntimeError(
            "Expected one selected-class row."
        )

    selected_class = selected_class_rows[0]

    if (
        parse_int(
            selected_class,
            "bridge_atoms_per_attachment",
        )
        != EXPECTED_BRIDGE_ATOMS_PER_PATH
        or parse_int(
            selected_class,
            "bonds_per_path",
        )
        != EXPECTED_BONDS_PER_PATH
    ):
        raise RuntimeError(
            "Selected-class metadata are inconsistent."
        )

    source_node_by_id = {
        row[
            "node_id"
        ]: row
        for row in source_nodes
    }

    node_rows: list[dict[str, Any]] = []
    edge_rows: list[dict[str, Any]] = []
    path_rows: list[dict[str, Any]] = []

    elements: dict[str, str] = {}
    node_types: dict[str, str] = {}
    node_ends: dict[str, str] = {}

    adjacency: dict[str, set[str]] = {}
    edge_pairs: set[tuple[str, str]] = set()

    def add_node(
        node_id: str,
        element: str,
        node_type: str,
        end: str,
        metadata: dict[str, Any] | None = None,
    ) -> None:
        if node_id in elements:
            raise RuntimeError(
                f"Duplicate node: {node_id}"
            )

        elements[node_id] = element
        node_types[node_id] = node_type
        node_ends[node_id] = end
        adjacency[node_id] = set()

        node_rows.append(
            {
                "node_id": node_id,
                "element": element,
                "node_type": node_type,
                "end": end,
                **(
                    metadata
                    if metadata is not None
                    else {}
                ),
                "coordinates_assigned": False,
                "formal_charge_assigned": False,
                "force_field_type_assigned": False,
            }
        )

    def add_edge(
        source: str,
        target: str,
        edge_type: str,
        end: str,
        heavy_atom_edge: bool,
        metadata: dict[str, Any] | None = None,
    ) -> None:
        if source == target:
            raise RuntimeError(
                f"Self edge requested for {source}"
            )

        if (
            source not in adjacency
            or target not in adjacency
        ):
            raise RuntimeError(
                f"Missing edge endpoint: {source} | {target}"
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
                f"Duplicate edge: {pair}"
            )

        edge_pairs.add(pair)
        adjacency[source].add(target)
        adjacency[target].add(source)

        source_element = elements[source]
        target_element = elements[target]

        heteropolar = (
            heavy_atom_edge
            and {
                source_element,
                target_element,
            }
            == {
                "B",
                "N",
            }
        )

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
                "heteropolar_BN_edge": heteropolar,
                **(
                    metadata
                    if metadata is not None
                    else {}
                ),
                "coordinates_assigned": False,
                "formal_bond_order_assigned": False,
            }
        )

    retained_heavy_ids = []

    for row in source_nodes:
        if row[
            "element"
        ] == "H":
            continue

        if row[
            "node_type"
        ] == "ALTERNATING_BN_TRIMER_BRIDGE":
            continue

        node_id = row[
            "node_id"
        ]

        retained_heavy_ids.append(
            node_id
        )

        metadata = {
            key: value
            for key, value in row.items()
            if (
                key
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
            )
        }

        metadata[
            "graph_source"
        ] = (
            "GATE3I_RETAINED_NONBRIDGE_HEAVY_NODE"
        )

        add_node(
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
            metadata=metadata,
        )

    removed_old_bridge_edges = 0
    removed_old_H_edges = 0
    retained_heavy_edges = 0

    for row in source_edges:
        if not parse_bool(
            row[
                "heavy_atom_edge"
            ]
        ):
            removed_old_H_edges += 1
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
            removed_old_bridge_edges += 1
            continue

        add_edge(
            source=source,
            target=target,
            edge_type=row[
                "edge_type"
            ],
            end=row[
                "end"
            ],
            heavy_atom_edge=True,
            metadata={
                "graph_source": (
                    "GATE3I_RETAINED_NONBRIDGE_HEAVY_EDGE"
                ),
            },
        )

        retained_heavy_edges += 1

    if removed_old_bridge_edges != 120:
        raise RuntimeError(
            "Unexpected removed trimer-edge count: "
            f"{removed_old_bridge_edges}/120"
        )

    selected_seed_ids_by_end: dict[
        str,
        set[str]
    ] = {}

    selected_annulus_ids_by_end: dict[
        str,
        set[str]
    ] = {}

    bridge_ids_by_end: dict[
        str,
        list[str]
    ] = {
        "LOWER": [],
        "UPPER": [],
    }

    path_edge_sets: dict[
        str,
        set[tuple[str, str]]
    ] = {}

    for end in (
        "LOWER",
        "UPPER",
    ):
        mapping = selected_end_rows[
            end
        ]

        seed_ids = split_nodes(
            mapping[
                "seed_nodes"
            ]
        )

        annulus_ids = split_nodes(
            mapping[
                "annulus_nodes"
            ]
        )

        if (
            len(seed_ids)
            != EXPECTED_ATTACHMENTS_PER_END
            or len(annulus_ids)
            != EXPECTED_ATTACHMENTS_PER_END
        ):
            raise RuntimeError(
                f"{end}: expected 15 seed and annulus endpoints."
            )

        if (
            len(set(seed_ids))
            != EXPECTED_ATTACHMENTS_PER_END
            or len(set(annulus_ids))
            != EXPECTED_ATTACHMENTS_PER_END
        ):
            raise RuntimeError(
                f"{end}: duplicate selected endpoint."
            )

        for node_id in (
            seed_ids
            + annulus_ids
        ):
            if node_id not in adjacency:
                raise RuntimeError(
                    f"{end}: selected endpoint missing from base graph: "
                    f"{node_id}"
                )

        seed_element = mapping[
            "seed_element"
        ]

        annulus_element = mapping[
            "annulus_endpoint_element"
        ]

        sequence = mapping[
            "bridge_element_sequence"
        ].split("-")

        if len(sequence) != EXPECTED_BRIDGE_ATOMS_PER_PATH:
            raise RuntimeError(
                f"{end}: invalid bridge sequence {sequence}"
            )

        expected_sequence = []
        current = seed_element

        for _ in range(
            EXPECTED_BRIDGE_ATOMS_PER_PATH
        ):
            current = opposite_element(
                current
            )

            expected_sequence.append(
                current
            )

        if sequence != expected_sequence:
            raise RuntimeError(
                f"{end}: sequence mismatch "
                f"{sequence} != {expected_sequence}"
            )

        if (
            opposite_element(
                sequence[-1]
            )
            != annulus_element
        ):
            raise RuntimeError(
                f"{end}: terminal bridge element does not "
                "connect heteropolarly to annulus endpoint."
            )

        if any(
            elements[node_id]
            != seed_element
            for node_id in seed_ids
        ):
            raise RuntimeError(
                f"{end}: selected seed element mismatch."
            )

        if any(
            elements[node_id]
            != annulus_element
            for node_id in annulus_ids
        ):
            raise RuntimeError(
                f"{end}: selected annulus element mismatch."
            )

        selected_seed_ids_by_end[
            end
        ] = set(seed_ids)

        selected_annulus_ids_by_end[
            end
        ] = set(annulus_ids)

        for bridge_index, (
            seed_id,
            annulus_id,
        ) in enumerate(
            zip(
                seed_ids,
                annulus_ids,
            )
        ):
            path_id = (
                f"{end}:BRIDGE4:"
                f"{bridge_index:02d}"
            )

            bridge_node_ids = []

            for bridge_position, element in enumerate(
                sequence,
                start=1,
            ):
                node_id = (
                    f"BR4:{end}:"
                    f"{bridge_index:02d}:"
                    f"{bridge_position}"
                )

                bridge_node_ids.append(
                    node_id
                )

                bridge_ids_by_end[
                    end
                ].append(
                    node_id
                )

                add_node(
                    node_id=node_id,
                    element=element,
                    node_type=(
                        "ALTERNATING_BN_FOUR_ATOM_BRIDGE"
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
                            bridge_position
                        ),
                        "attached_seed_node": (
                            seed_id
                        ),
                        "attached_annulus_node": (
                            annulus_id
                        ),
                        "mapping_id": (
                            mapping[
                                "mapping_id"
                            ]
                        ),
                        "mapping_seed_parity": (
                            mapping[
                                "seed_parity"
                            ]
                        ),
                        "mapping_orientation": (
                            mapping[
                                "orientation"
                            ]
                        ),
                        "mapping_rotation": (
                            mapping[
                                "rotation"
                            ]
                        ),
                        "graph_source": (
                            "GATE3L_SELECTED_FOUR_ATOM_BRIDGE"
                        ),
                    },
                )

            path_nodes = [
                seed_id,
                *bridge_node_ids,
                annulus_id,
            ]

            path_edges = set()

            for edge_position in range(
                len(path_nodes) - 1
            ):
                first = path_nodes[
                    edge_position
                ]

                second = path_nodes[
                    edge_position + 1
                ]

                pair = tuple(
                    sorted(
                        (
                            first,
                            second,
                        )
                    )
                )

                path_edges.add(pair)

                add_edge(
                    source=first,
                    target=second,
                    edge_type=(
                        "ALTERNATING_BN_FOUR_ATOM_BRIDGE"
                    ),
                    end=end,
                    heavy_atom_edge=True,
                    metadata={
                        "bridge_path_id": (
                            path_id
                        ),
                        "bridge_edge_position": (
                            edge_position + 1
                        ),
                        "graph_source": (
                            "GATE3L_SELECTED_FOUR_ATOM_BRIDGE"
                        ),
                    },
                )

            path_edge_sets[
                path_id
            ] = path_edges

            path_rows.append(
                {
                    "bridge_path_id": (
                        path_id
                    ),
                    "end": end,
                    "bridge_index": (
                        bridge_index
                    ),
                    "seed_node": (
                        seed_id
                    ),
                    "seed_element": (
                        elements[
                            seed_id
                        ]
                    ),
                    "bridge_node_1": (
                        bridge_node_ids[0]
                    ),
                    "bridge_element_1": (
                        elements[
                            bridge_node_ids[0]
                        ]
                    ),
                    "bridge_node_2": (
                        bridge_node_ids[1]
                    ),
                    "bridge_element_2": (
                        elements[
                            bridge_node_ids[1]
                        ]
                    ),
                    "bridge_node_3": (
                        bridge_node_ids[2]
                    ),
                    "bridge_element_3": (
                        elements[
                            bridge_node_ids[2]
                        ]
                    ),
                    "bridge_node_4": (
                        bridge_node_ids[3]
                    ),
                    "bridge_element_4": (
                        elements[
                            bridge_node_ids[3]
                        ]
                    ),
                    "annulus_node": (
                        annulus_id
                    ),
                    "annulus_element": (
                        elements[
                            annulus_id
                        ]
                    ),
                    "heavy_edges_in_path": (
                        EXPECTED_BONDS_PER_PATH
                    ),
                    "mapping_id": (
                        mapping[
                            "mapping_id"
                        ]
                    ),
                    "mapping_seed_parity": (
                        mapping[
                            "seed_parity"
                        ]
                    ),
                    "mapping_orientation": (
                        mapping[
                            "orientation"
                        ]
                    ),
                    "mapping_rotation": (
                        mapping[
                            "rotation"
                        ]
                    ),
                }
            )

    heavy_nodes = {
        node_id
        for node_id, element
        in elements.items()
        if element != "H"
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

    H_count_by_end_role: dict[
        str,
        Counter[str]
    ] = {
        "LOWER": Counter(),
        "UPPER": Counter(),
        "PARENT": Counter(),
    }

    hydrogen_index = 0

    for heavy_id in sorted(
        heavy_nodes
    ):
        heavy_degree = len(
            heavy_adjacency[
                heavy_id
            ]
        )

        if heavy_degree > 3:
            raise RuntimeError(
                f"Heavy degree exceeds 3: "
                f"{heavy_id} -> {heavy_degree}"
            )

        required_H = (
            3
            - heavy_degree
        )

        if required_H < 0:
            raise RuntimeError(
                f"Negative H requirement for {heavy_id}"
            )

        for local_H_index in range(
            required_H
        ):
            end = node_ends[
                heavy_id
            ]

            node_type = node_types[
                heavy_id
            ]

            if node_type == "HEXAGONAL_EDGE_COMPLETION_SEED":
                role = "SEED_PASSIVANT_H"

            elif node_type == "ANNULUS_OUTER_BOUNDARY":
                role = "ANNULUS_OUTER_PASSIVANT_H"

            elif node_type == "ANNULUS_INNER_BOUNDARY":
                role = "ANNULUS_INNER_PASSIVANT_H"

            elif node_type == "ALTERNATING_BN_FOUR_ATOM_BRIDGE":
                role = "BRIDGE_PASSIVANT_H"

            else:
                role = "PARENT_OR_OTHER_PASSIVANT_H"

            hydrogen_id = (
                f"H4:{end}:"
                f"{hydrogen_index:04d}:"
                f"{local_H_index}"
            )

            hydrogen_index += 1

            add_node(
                node_id=hydrogen_id,
                element="H",
                node_type=role,
                end=end,
                metadata={
                    "attached_to": (
                        heavy_id
                    ),
                    "graph_source": (
                        "COORDINATION_DERIVED_PASSIVATION"
                    ),
                },
            )

            add_edge(
                source=heavy_id,
                target=hydrogen_id,
                edge_type=(
                    f"{role}_BOND"
                ),
                end=end,
                heavy_atom_edge=False,
                metadata={
                    "graph_source": (
                        "COORDINATION_DERIVED_PASSIVATION"
                    ),
                },
            )

            if end not in H_count_by_end_role:
                H_count_by_end_role[
                    end
                ] = Counter()

            H_count_by_end_role[
                end
            ][role] += 1

    all_heavy_nodes = {
        node_id
        for node_id, element
        in elements.items()
        if element != "H"
    }

    hydrogen_nodes = {
        node_id
        for node_id, element
        in elements.items()
        if element == "H"
    }

    heavy_adjacency = {
        node_id: {
            neighbor
            for neighbor in adjacency[
                node_id
            ]
            if neighbor in all_heavy_nodes
        }
        for node_id in all_heavy_nodes
    }

    full_components = connected_components(
        adjacency
    )

    heavy_components = connected_components(
        heavy_adjacency
    )

    heavy_bipartite = bipartite(
        heavy_adjacency
    )

    four_cycles = count_four_cycles(
        heavy_adjacency
    )

    heavy_degree_over3 = [
        node_id
        for node_id in all_heavy_nodes
        if len(
            heavy_adjacency[
                node_id
            ]
        )
        > 3
    ]

    heavy_total_degree_failures = [
        node_id
        for node_id in all_heavy_nodes
        if len(
            adjacency[
                node_id
            ]
        )
        != 3
    ]

    H_degree_failures = [
        node_id
        for node_id in hydrogen_nodes
        if len(
            adjacency[
                node_id
            ]
        )
        != 1
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

    cycle_lengths = []
    paths_without_cycle = 0

    for row in path_rows:
        path_id = row[
            "bridge_path_id"
        ]

        alternative = shortest_path_length(
            heavy_adjacency,
            row[
                "seed_node"
            ],
            row[
                "annulus_node"
            ],
            path_edge_sets[
                path_id
            ],
        )

        if alternative is None:
            paths_without_cycle += 1
            cycle_length: int | None = None

        else:
            cycle_length = (
                alternative
                + EXPECTED_BONDS_PER_PATH
            )

            cycle_lengths.append(
                cycle_length
            )

        row[
            "alternative_heavy_path_length"
        ] = (
            ""
            if alternative is None
            else alternative
        )

        row[
            "shortest_cycle_containing_bridge_path"
        ] = (
            ""
            if cycle_length is None
            else cycle_length
        )

    if not cycle_lengths:
        raise RuntimeError(
            "No bridge-containing cycles were resolved."
        )

    cycle_minimum = min(
        cycle_lengths
    )

    cycle_maximum = max(
        cycle_lengths
    )

    heavy_edge_count = sum(
        parse_bool(
            row[
                "heavy_atom_edge"
            ]
        )
        for row in edge_rows
    )

    H_edge_count = (
        len(edge_rows)
        - heavy_edge_count
    )

    type_counts = Counter(
        node_types.values()
    )

    bridge_composition_by_end = {
        "LOWER": Counter(),
        "UPPER": Counter(),
    }

    for node_id in all_heavy_nodes:
        if node_types[
            node_id
        ] != "ALTERNATING_BN_FOUR_ATOM_BRIDGE":
            continue

        bridge_composition_by_end[
            node_ends[
                node_id
            ]
        ][
            elements[
                node_id
            ]
        ] += 1

    end_summary_rows = []

    for end in (
        "LOWER",
        "UPPER",
    ):
        mapping = selected_end_rows[
            end
        ]

        end_paths = [
            row
            for row in path_rows
            if row[
                "end"
            ]
            == end
        ]

        end_cycles = [
            int(
                row[
                    "shortest_cycle_containing_bridge_path"
                ]
            )
            for row in end_paths
            if row[
                "shortest_cycle_containing_bridge_path"
            ]
            != ""
        ]

        H_roles = H_count_by_end_role[
            end
        ]

        end_summary_rows.append(
            {
                "end": end,
                "mapping_id": (
                    mapping[
                        "mapping_id"
                    ]
                ),
                "seed_element": (
                    mapping[
                        "seed_element"
                    ]
                ),
                "bridge_sequence": (
                    mapping[
                        "bridge_element_sequence"
                    ]
                ),
                "annulus_endpoint_element": (
                    mapping[
                        "annulus_endpoint_element"
                    ]
                ),
                "bridge_paths": (
                    len(
                        end_paths
                    )
                ),
                "bridge_atoms": (
                    len(
                        bridge_ids_by_end[
                            end
                        ]
                    )
                ),
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
                    H_roles[
                        "SEED_PASSIVANT_H"
                    ]
                ),
                "outer_H_passivants": (
                    H_roles[
                        "ANNULUS_OUTER_PASSIVANT_H"
                    ]
                ),
                "inner_H_passivants": (
                    H_roles[
                        "ANNULUS_INNER_PASSIVANT_H"
                    ]
                ),
                "bridge_H_passivants": (
                    H_roles[
                        "BRIDGE_PASSIVANT_H"
                    ]
                ),
                "other_H_passivants": (
                    H_roles[
                        "PARENT_OR_OTHER_PASSIVANT_H"
                    ]
                ),
                "total_H_passivants": sum(
                    H_roles.values()
                ),
                "cycle_minimum": min(
                    end_cycles
                ),
                "cycle_maximum": max(
                    end_cycles
                ),
                "screen_minimum_angle_deg": (
                    mapping[
                        "minimum_angle_deg"
                    ]
                ),
                "screen_minimum_local_clearance_nm": (
                    mapping[
                        "minimum_local_clearance_nm"
                    ]
                ),
                "screen_minimum_interbridge_clearance_nm": (
                    mapping[
                        "minimum_interbridge_clearance_nm"
                    ]
                ),
                "screen_maximum_bond_deviation_nm": (
                    mapping[
                        "maximum_bond_deviation_nm"
                    ]
                ),
            }
        )

    gates = {
        "Gate3I_source_graph_is_accepted": (
            source_graph_summary.get(
                "decision"
            )
            == EXPECTED_SOURCE_GRAPH_DECISION
        ),
        "Gate3L_longer_bridge_screen_is_accepted": (
            screen_summary.get(
                "decision"
            )
            == EXPECTED_SCREEN_DECISION
        ),
        "selected_bridge_class_has_four_atoms": (
            selected_bridge_atoms
            == EXPECTED_BRIDGE_ATOMS_PER_PATH
        ),
        "selected_architecture_has_15_attachments_per_end": (
            selected_attachments
            == EXPECTED_ATTACHMENTS_PER_END
        ),
        "30_old_trimer_paths_were_replaced": (
            removed_old_bridge_edges
            == 120
        ),
        "30_new_four_atom_bridge_paths_were_built": (
            len(
                path_rows
            )
            == EXPECTED_PATHS_TOTAL
        ),
        "15_bridge_paths_were_built_per_end": all(
            int(
                row[
                    "bridge_paths"
                ]
            )
            == EXPECTED_ATTACHMENTS_PER_END
            for row in end_summary_rows
        ),
        "60_bridge_atoms_were_added_per_end": all(
            int(
                row[
                    "bridge_atoms"
                ]
            )
            == EXPECTED_BRIDGE_ATOMS_PER_END
            for row in end_summary_rows
        ),
        "120_bridge_atoms_were_added_total": (
            type_counts[
                "ALTERNATING_BN_FOUR_ATOM_BRIDGE"
            ]
            == EXPECTED_BRIDGE_ATOMS_TOTAL
        ),
        "2112_heavy_atoms_are_present": (
            len(
                all_heavy_nodes
            )
            == EXPECTED_HEAVY_ATOMS
        ),
        "204_H_atoms_are_present": (
            len(
                hydrogen_nodes
            )
            == EXPECTED_H_ATOMS
        ),
        "2316_total_nodes_are_present": (
            len(
                adjacency
            )
            == EXPECTED_TOTAL_NODES
        ),
        "3066_heavy_edges_are_present": (
            heavy_edge_count
            == EXPECTED_HEAVY_EDGES
        ),
        "204_H_edges_are_present": (
            H_edge_count
            == EXPECTED_H_EDGES
        ),
        "all_heavy_atoms_have_total_coordination3": (
            len(
                heavy_total_degree_failures
            )
            == 0
        ),
        "all_H_atoms_have_coordination1": (
            len(
                H_degree_failures
            )
            == 0
        ),
        "no_heavy_atom_has_heavy_degree_above3": (
            len(
                heavy_degree_over3
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
            four_cycles == 0
        ),
        "every_bridge_path_participates_in_a_cycle": (
            paths_without_cycle == 0
        ),
        "all_bridge_containing_cycles_have_length16": (
            cycle_minimum
            == EXPECTED_CYCLE_LENGTH
            and cycle_maximum
            == EXPECTED_CYCLE_LENGTH
        ),
        "no_coordinates_were_assigned": all(
            not parse_bool(
                row.get(
                    "coordinates_assigned",
                    "false",
                )
            )
            for row in node_rows
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
        else REVIEW_DECISION
    )

    required_next_step = (
        "BUILD_AND_VALIDATE_R2_SELECTED_FOUR_ATOM_BN_"
        "BRIDGE_STATIC_COORDINATE_EMBEDDING"
        if accepted
        else
        "REVIEW_R2_SELECTED_FOUR_ATOM_BN_BRIDGE_GRAPH_FAILURES"
    )

    summary = {
        "decision": decision,
        "bridge_atoms_per_path": (
            EXPECTED_BRIDGE_ATOMS_PER_PATH
        ),
        "bonds_per_path": (
            EXPECTED_BONDS_PER_PATH
        ),
        "attachments_per_end": (
            EXPECTED_ATTACHMENTS_PER_END
        ),
        "bridge_paths_total": (
            len(
                path_rows
            )
        ),
        "bridge_atoms_lower": (
            len(
                bridge_ids_by_end[
                    "LOWER"
                ]
            )
        ),
        "bridge_atoms_upper": (
            len(
                bridge_ids_by_end[
                    "UPPER"
                ]
            )
        ),
        "bridge_atoms_total": (
            type_counts[
                "ALTERNATING_BN_FOUR_ATOM_BRIDGE"
            ]
        ),
        "total_heavy_atoms": (
            len(
                all_heavy_nodes
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
        "heavy_edges": (
            heavy_edge_count
        ),
        "H_edges": (
            H_edge_count
        ),
        "heavy_total_degree_failures": (
            len(
                heavy_total_degree_failures
            )
        ),
        "H_degree_failures": (
            len(
                H_degree_failures
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
        "heavy_graph_four_cycles": (
            four_cycles
        ),
        "bridge_paths_without_cycle": (
            paths_without_cycle
        ),
        "bridge_cycle_minimum": (
            cycle_minimum
        ),
        "bridge_cycle_maximum": (
            cycle_maximum
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

    write_rows(
        GRAPH_NODES,
        node_rows,
    )

    write_rows(
        GRAPH_EDGES,
        edge_rows,
    )

    write_rows(
        BRIDGE_PATHS,
        path_rows,
    )

    write_rows(
        END_SUMMARY,
        end_summary_rows,
    )

    write_rows(
        SUMMARY,
        [
            summary
        ],
    )

    write_rows(
        GATES,
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
                "limitations": [
                    (
                        "This gate validates graph connectivity and "
                        "coordination only."
                    ),
                    (
                        "The selected conformers from Gate 3L are not "
                        "yet applied as coordinates."
                    ),
                    (
                        "Hydrogen atoms are added from graph coordination "
                        "requirements, not from protonation energetics."
                    ),
                    (
                        "No molecular topology, charges, force-field "
                        "parameters, minimization, MD or QM calculation "
                        "was generated."
                    ),
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    manifest_rows = [
        {
            "role": (
                "Gate3I_source_nodes"
            ),
            "file": relative(
                SOURCE_NODES
            ),
            "sha256": sha256(
                SOURCE_NODES
            ),
        },
        {
            "role": (
                "Gate3I_source_edges"
            ),
            "file": relative(
                SOURCE_EDGES
            ),
            "sha256": sha256(
                SOURCE_EDGES
            ),
        },
        {
            "role": (
                "Gate3I_source_summary"
            ),
            "file": relative(
                SOURCE_GRAPH_SUMMARY
            ),
            "sha256": sha256(
                SOURCE_GRAPH_SUMMARY
            ),
        },
        {
            "role": (
                "Gate3L_selected_candidate"
            ),
            "file": relative(
                SELECTED_CANDIDATE
            ),
            "sha256": sha256(
                SELECTED_CANDIDATE
            ),
        },
        {
            "role": (
                "Gate3L_screen_summary"
            ),
            "file": relative(
                SCREEN_SUMMARY
            ),
            "sha256": sha256(
                SCREEN_SUMMARY
            ),
        },
    ]

    write_rows(
        MANIFEST,
        manifest_rows,
    )

    end_lines = "\n".join(
        (
            f"### {row['end']}\n\n"
            f"- Mapping: **{row['mapping_id']}**\n"
            f"- Bridge sequence: **{row['bridge_sequence']}**\n"
            f"- Paths/bridge atoms: "
            f"**{row['bridge_paths']}/{row['bridge_atoms']}**\n"
            f"- Bridge B/N: "
            f"**{row['bridge_B_atoms']}/{row['bridge_N_atoms']}**\n"
            f"- H seed/outer/inner/bridge/other/total: "
            f"**{row['seed_H_passivants']}/"
            f"{row['outer_H_passivants']}/"
            f"{row['inner_H_passivants']}/"
            f"{row['bridge_H_passivants']}/"
            f"{row['other_H_passivants']}/"
            f"{row['total_H_passivants']}**\n"
            f"- Bridge-containing cycles: "
            f"**{row['cycle_minimum']}–{row['cycle_maximum']}**"
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

    REPORT.write_text(
        f"""# R2 Selected Full-Density Four-Atom BN Bridge Graph

## Scope

This gate replaces the previous three-atom BN bridges with the shortest
full-density longer bridge class selected by Gate 3L.

No coordinates, molecular topology, formal charges, force-field
parameters, minimization, MD or QM calculation were generated.

## Architecture

- Bridge atoms per path:
  **{EXPECTED_BRIDGE_ATOMS_PER_PATH}**
- Bonds per path:
  **{EXPECTED_BONDS_PER_PATH}**
- Attachments per end:
  **{EXPECTED_ATTACHMENTS_PER_END}**
- Total paths:
  **{len(path_rows)}**
- Total bridge atoms:
  **{type_counts['ALTERNATING_BN_FOUR_ATOM_BRIDGE']}**

{end_lines}

## Complete graph

- Heavy atoms:
  **{len(all_heavy_nodes)}**
- H atoms:
  **{len(hydrogen_nodes)}**
- Total nodes:
  **{len(adjacency)}**
- Heavy/H edges:
  **{heavy_edge_count}/{H_edge_count}**
- Heavy/H degree failures:
  **{len(heavy_total_degree_failures)}/
  {len(H_degree_failures)}**
- Nonheteropolar heavy edges:
  **{len(nonheteropolar_heavy_edges)}**
- Full/heavy connected components:
  **{len(full_components)}/{len(heavy_components)}**
- Bipartite:
  **{heavy_bipartite}**
- Four-member cycles:
  **{four_cycles}**
- Bridge cycle range:
  **{cycle_minimum}–{cycle_maximum}**

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
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
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 selected full-density four-atom "
        "BN bridge graph completed."
    )

    print(
        "Removed old trimer heavy edges / old H edges: "
        f"{removed_old_bridge_edges}/"
        f"{removed_old_H_edges}"
    )

    print(
        "Bridge paths lower/upper/total: "
        f"{sum(row['end'] == 'LOWER' for row in path_rows)}/"
        f"{sum(row['end'] == 'UPPER' for row in path_rows)}/"
        f"{len(path_rows)}"
    )

    print(
        "Bridge atoms lower/upper/total: "
        f"{len(bridge_ids_by_end['LOWER'])}/"
        f"{len(bridge_ids_by_end['UPPER'])}/"
        f"{type_counts['ALTERNATING_BN_FOUR_ATOM_BRIDGE']}"
    )

    for row in end_summary_rows:
        print(
            f"{row['end']} mapping / sequence / "
            "bridge B-N / H seed-outer-inner-bridge-other-total: "
            f"{row['mapping_id']}/"
            f"{row['bridge_sequence']}/"
            f"{row['bridge_B_atoms']}-"
            f"{row['bridge_N_atoms']}/"
            f"{row['seed_H_passivants']}-"
            f"{row['outer_H_passivants']}-"
            f"{row['inner_H_passivants']}-"
            f"{row['bridge_H_passivants']}-"
            f"{row['other_H_passivants']}-"
            f"{row['total_H_passivants']}"
        )

    print(
        "Total heavy / H / all nodes: "
        f"{len(all_heavy_nodes)}/"
        f"{len(hydrogen_nodes)}/"
        f"{len(adjacency)}"
    )

    print(
        "Heavy / H edges: "
        f"{heavy_edge_count}/"
        f"{H_edge_count}"
    )

    print(
        "Heavy total-degree / H-degree failures: "
        f"{len(heavy_total_degree_failures)}/"
        f"{len(H_degree_failures)}"
    )

    print(
        "Nonheteropolar heavy edges: "
        f"{len(nonheteropolar_heavy_edges)}"
    )

    print(
        "Full/heavy components / bipartite / four-cycles: "
        f"{len(full_components)}/"
        f"{len(heavy_components)}/"
        f"{heavy_bipartite}/"
        f"{four_cycles}"
    )

    print(
        "Paths without cycle / bridge-cycle min-max: "
        f"{paths_without_cycle}/"
        f"{cycle_minimum}/"
        f"{cycle_maximum}"
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
        GRAPH_NODES,
        GRAPH_EDGES,
        BRIDGE_PATHS,
        END_SUMMARY,
        SUMMARY,
        GATES,
        GRAPH_JSON,
        MANIFEST,
        REPORT,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not accepted:
        raise RuntimeError(
            "Selected four-atom BN bridge graph requires review."
        )


if __name__ == "__main__":
    main()
