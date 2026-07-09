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

GATE3B_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "02_r2_polar_end_specific_candidate_ranking"
)

OUTPUT_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "03_r2_graded_heteropolar_bn_collar_connectivity_blueprint"
)

PARENT_SUMMARY = (
    GATE3A_ROOT
    / "r2_parent_rim_chemical_audit_summary.csv"
)

TERMINAL_ATOMS = (
    GATE3A_ROOT
    / "r2_parent_terminal_rim_atoms.csv"
)

CANDIDATE_SELECTION = (
    GATE3B_ROOT
    / "r2_polar_end_specific_candidate_selection_summary.csv"
)

DESIGN_CONTRACT = (
    GATE3B_ROOT
    / "r2_primary_graded_bn_collar_design_contract.json"
)

SEQUENCES_CSV = (
    OUTPUT_ROOT
    / "r2_feasible_collar_ring_population_sequences.csv"
)

NODES_CSV = (
    OUTPUT_ROOT
    / "r2_selected_collar_connectivity_nodes.csv"
)

EDGES_CSV = (
    OUTPUT_ROOT
    / "r2_selected_collar_connectivity_edges.csv"
)

INTERFACES_CSV = (
    OUTPUT_ROOT
    / "r2_selected_collar_interface_summary.csv"
)

BLUEPRINT_SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_graded_collar_connectivity_blueprint_summary.csv"
)

GATES_CSV = (
    OUTPUT_ROOT
    / "r2_graded_collar_connectivity_blueprint_gates.csv"
)

BLUEPRINT_JSON = (
    OUTPUT_ROOT
    / "r2_graded_collar_connectivity_blueprint.json"
)

SOURCE_MANIFEST_CSV = (
    OUTPUT_ROOT
    / "r2_graded_collar_connectivity_source_manifest.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_GRADED_HETEROPOLAR_BN_COLLAR_CONNECTIVITY_BLUEPRINT_DAY024.md"
)

EXPECTED_GATE3A_DECISION = (
    "R2_PARENT_RIM_AND_CHEMICAL_CONSTRAINTS_AUDITED"
)

EXPECTED_GATE3B_DECISION = (
    "R2_POLAR_END_SPECIFIC_VALENCE_COMPLETION_CANDIDATES_RANKED"
)

EXPECTED_PRIMARY_CANDIDATE = (
    "C2_GRADED_HETEROPOLAR_BN_COLLAR_ANNULUS"
)

PASS_DECISION = (
    "R2_GRADED_HETEROPOLAR_BN_COLLAR_CONNECTIVITY_BLUEPRINT_VALIDATED"
)

EXPECTED_RING_POPULATIONS = (
    30,
    30,
    30,
    22,
    20,
    12,
)

EXPECTED_PARENT_TERMINALS_PER_END = 30
EXPECTED_PARENT_NEW_BONDS_PER_END = 60
EXPECTED_PARENT_NEW_BONDS_TOTAL = 120

EXPECTED_ADDED_HEAVY_ATOMS_PER_END = 144
EXPECTED_ADDED_HEAVY_ATOMS_TOTAL = 288

EXPECTED_INNER_PASSIVANTS_PER_END = 12
EXPECTED_INNER_PASSIVANTS_TOTAL = 24

EXPECTED_STERIC_BEADS_PER_END = 144
EXPECTED_STERIC_BEADS_TOTAL = 288

MAX_INTERFACE_ANGULAR_SPAN_TURNS = 0.10
MAX_HEAVY_ATOM_ESTIMATE_RELATIVE_ERROR = 0.02
MAX_INNER_SPACING_RELATIVE_ERROR = 0.15

MIN_SEQUENCE_LAYERS = 4
MAX_SEQUENCE_LAYERS = 8
MIN_RING_POPULATION = 6

EXPECTED_SELECTED_LAYER_COUNT = 6

SEQUENCE_SCORE_WEIGHTS = {
    "heavy_atom_error": 35.0,
    "inner_spacing_error": 45.0,
    "layer_count_error": 4.0,
    "abrupt_contraction": 10.0,
}


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
            f"Expected one row in {path}; "
            f"found {len(rows)}"
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


class FlowEdge:
    def __init__(
        self,
        to: int,
        reverse: int,
        capacity: int,
        original_capacity: int,
        tag: tuple[int, int] | None = None,
    ) -> None:
        self.to = to
        self.reverse = reverse
        self.capacity = capacity
        self.original_capacity = original_capacity
        self.tag = tag


class Dinic:
    def __init__(
        self,
        node_count: int,
    ) -> None:
        self.graph: list[
            list[FlowEdge]
        ] = [
            []
            for _ in range(node_count)
        ]

    def add_edge(
        self,
        source: int,
        target: int,
        capacity: int,
        tag: tuple[int, int] | None = None,
    ) -> None:
        forward = FlowEdge(
            to=target,
            reverse=len(
                self.graph[target]
            ),
            capacity=capacity,
            original_capacity=capacity,
            tag=tag,
        )

        reverse = FlowEdge(
            to=source,
            reverse=len(
                self.graph[source]
            ),
            capacity=0,
            original_capacity=0,
            tag=None,
        )

        self.graph[source].append(
            forward
        )

        self.graph[target].append(
            reverse
        )

    def max_flow(
        self,
        source: int,
        sink: int,
    ) -> int:
        total_flow = 0
        node_count = len(self.graph)

        while True:
            levels = [-1] * node_count
            levels[source] = 0

            queue: deque[int] = deque(
                [source]
            )

            while queue:
                node = queue.popleft()

                for edge in self.graph[node]:
                    if (
                        edge.capacity > 0
                        and levels[edge.to] < 0
                    ):
                        levels[edge.to] = (
                            levels[node] + 1
                        )

                        queue.append(
                            edge.to
                        )

            if levels[sink] < 0:
                return total_flow

            iterators = [0] * node_count

            def send_flow(
                node: int,
                available: int,
            ) -> int:
                if node == sink:
                    return available

                while (
                    iterators[node]
                    < len(
                        self.graph[node]
                    )
                ):
                    edge_index = (
                        iterators[node]
                    )

                    edge = self.graph[
                        node
                    ][
                        edge_index
                    ]

                    if (
                        edge.capacity > 0
                        and levels[edge.to]
                        == levels[node] + 1
                    ):
                        sent = send_flow(
                            edge.to,
                            min(
                                available,
                                edge.capacity,
                            ),
                        )

                        if sent > 0:
                            edge.capacity -= sent

                            reverse_edge = (
                                self.graph[
                                    edge.to
                                ][
                                    edge.reverse
                                ]
                            )

                            reverse_edge.capacity += sent

                            return sent

                    iterators[node] += 1

                return 0

            while True:
                sent = send_flow(
                    source,
                    10**9,
                )

                if sent == 0:
                    break

                total_flow += sent


def circular_distance_turns(
    left_index: int,
    left_count: int,
    right_index: int,
    right_count: int,
) -> float:
    left_angle = (
        left_index
        / left_count
    )

    right_angle = (
        right_index
        / right_count
    )

    difference = abs(
        left_angle
        - right_angle
    )

    return min(
        difference,
        1.0 - difference,
    )


def balanced_degree_sequence(
    node_count: int,
    total_degree: int,
) -> list[int]:
    if not (
        node_count
        <= total_degree
        <= 2 * node_count
    ):
        raise RuntimeError(
            "Balanced degree sequence requires each "
            "node to have degree 1 or 2: "
            f"nodes={node_count}, total={total_degree}"
        )

    degrees = [1] * node_count
    extras = (
        total_degree
        - node_count
    )

    if extras == 0:
        return degrees

    if extras == node_count:
        return [2] * node_count

    selected: set[int] = set()

    for index in range(extras):
        candidate = int(
            math.floor(
                (
                    index + 0.5
                )
                * node_count
                / extras
            )
        ) % node_count

        while candidate in selected:
            candidate = (
                candidate + 1
            ) % node_count

        selected.add(candidate)

    for index in selected:
        degrees[index] = 2

    if sum(degrees) != total_degree:
        raise RuntimeError(
            "Balanced degree sequence total mismatch."
        )

    return degrees


def realize_local_interface(
    left_degrees: list[int],
    right_degrees: list[int],
) -> tuple[
    list[tuple[int, int]],
    float,
]:
    if sum(left_degrees) != sum(
        right_degrees
    ):
        raise RuntimeError(
            "Interface degree sums do not match."
        )

    left_count = len(
        left_degrees
    )

    right_count = len(
        right_degrees
    )

    required_flow = sum(
        left_degrees
    )

    thresholds = sorted(
        {
            circular_distance_turns(
                left_index,
                left_count,
                right_index,
                right_count,
            )
            for left_index in range(
                left_count
            )
            for right_index in range(
                right_count
            )
        }
    )

    for threshold in thresholds:
        source = 0
        left_offset = 1
        right_offset = (
            left_offset
            + left_count
        )

        sink = (
            right_offset
            + right_count
        )

        flow = Dinic(
            sink + 1
        )

        for left_index, degree in enumerate(
            left_degrees
        ):
            flow.add_edge(
                source,
                left_offset
                + left_index,
                degree,
            )

        for left_index in range(
            left_count
        ):
            for right_index in range(
                right_count
            ):
                distance = (
                    circular_distance_turns(
                        left_index,
                        left_count,
                        right_index,
                        right_count,
                    )
                )

                if (
                    distance
                    <= threshold
                    + 1.0e-12
                ):
                    flow.add_edge(
                        left_offset
                        + left_index,
                        right_offset
                        + right_index,
                        1,
                        tag=(
                            left_index,
                            right_index,
                        ),
                    )

        for right_index, degree in enumerate(
            right_degrees
        ):
            flow.add_edge(
                right_offset
                + right_index,
                sink,
                degree,
            )

        observed_flow = flow.max_flow(
            source,
            sink,
        )

        if observed_flow != required_flow:
            continue

        pairs: list[
            tuple[int, int]
        ] = []

        for left_index in range(
            left_count
        ):
            node = (
                left_offset
                + left_index
            )

            for edge in flow.graph[node]:
                if (
                    edge.tag is not None
                    and edge.original_capacity == 1
                    and edge.capacity == 0
                ):
                    pairs.append(
                        edge.tag
                    )

        if len(pairs) != required_flow:
            raise RuntimeError(
                "Extracted interface-edge count mismatch."
            )

        if len(set(pairs)) != len(
            pairs
        ):
            raise RuntimeError(
                "Duplicate interface edges were produced."
            )

        observed_left = [0] * left_count
        observed_right = [0] * right_count

        for left_index, right_index in pairs:
            observed_left[
                left_index
            ] += 1

            observed_right[
                right_index
            ] += 1

        if observed_left != left_degrees:
            raise RuntimeError(
                "Left interface degrees were not reproduced."
            )

        if observed_right != right_degrees:
            raise RuntimeError(
                "Right interface degrees were not reproduced."
            )

        maximum_span = max(
            circular_distance_turns(
                left_index,
                left_count,
                right_index,
                right_count,
            )
            for left_index, right_index
            in pairs
        )

        return pairs, maximum_span

    raise RuntimeError(
        "No local simple bipartite interface realization "
        "was found."
    )


def enumerate_population_sequences(
    target_heavy_atoms: float,
    parent_radius_nm: float,
    aperture_radius_nm: float,
) -> list[dict[str, Any]]:
    parent_spacing_nm = (
        2.0
        * math.pi
        * parent_radius_nm
        / EXPECTED_PARENT_TERMINALS_PER_END
    )

    rows: list[
        dict[str, Any]
    ] = []

    def record_sequence(
        sequence: list[int],
    ) -> None:
        total_heavy_atoms = sum(
            sequence
        )

        inner_population = (
            sequence[-1]
        )

        inner_spacing_nm = (
            2.0
            * math.pi
            * aperture_radius_nm
            / inner_population
        )

        heavy_atom_relative_error = (
            abs(
                total_heavy_atoms
                - target_heavy_atoms
            )
            / target_heavy_atoms
        )

        inner_spacing_relative_error = (
            abs(
                inner_spacing_nm
                - parent_spacing_nm
            )
            / parent_spacing_nm
        )

        layer_count_relative_error = (
            abs(
                len(sequence)
                - EXPECTED_SELECTED_LAYER_COUNT
            )
            / EXPECTED_SELECTED_LAYER_COUNT
        )

        abrupt_contraction_penalty = sum(
            max(
                0.0,
                (
                    sequence[index]
                    - sequence[index + 1]
                )
                / sequence[index]
                - 0.35,
            )
            for index in range(
                len(sequence) - 1
            )
        )

        score = (
            100.0
            - SEQUENCE_SCORE_WEIGHTS[
                "heavy_atom_error"
            ]
            * heavy_atom_relative_error
            - SEQUENCE_SCORE_WEIGHTS[
                "inner_spacing_error"
            ]
            * inner_spacing_relative_error
            - SEQUENCE_SCORE_WEIGHTS[
                "layer_count_error"
            ]
            * layer_count_relative_error
            - SEQUENCE_SCORE_WEIGHTS[
                "abrupt_contraction"
            ]
            * abrupt_contraction_penalty
        )

        rows.append(
            {
                "ring_populations": (
                    "-".join(
                        str(value)
                        for value in sequence
                    )
                ),
                "layer_count": len(
                    sequence
                ),
                "total_added_heavy_atoms_per_end": (
                    total_heavy_atoms
                ),
                "inner_boundary_population": (
                    inner_population
                ),
                "parent_terminal_spacing_nm": (
                    parent_spacing_nm
                ),
                "inner_boundary_spacing_nm": (
                    inner_spacing_nm
                ),
                "heavy_atom_relative_error": (
                    heavy_atom_relative_error
                ),
                "inner_spacing_relative_error": (
                    inner_spacing_relative_error
                ),
                "layer_count_relative_error": (
                    layer_count_relative_error
                ),
                "abrupt_contraction_penalty": (
                    abrupt_contraction_penalty
                ),
                "sequence_screening_score": score,
            }
        )

    def recurse(
        sequence: list[int],
        incoming_edges: int,
    ) -> None:
        current_population = (
            sequence[-1]
        )

        if (
            len(sequence)
            >= MIN_SEQUENCE_LAYERS
            and incoming_edges
            == 2 * current_population
        ):
            record_sequence(
                sequence
            )

        if len(sequence) >= MAX_SEQUENCE_LAYERS:
            return

        outgoing_edges = (
            3
            * current_population
            - incoming_edges
        )

        if outgoing_edges <= 0:
            return

        for next_population in range(
            current_population,
            MIN_RING_POPULATION - 1,
            -1,
        ):
            if not (
                max(
                    current_population,
                    next_population,
                )
                <= outgoing_edges
                <= 2
                * min(
                    current_population,
                    next_population,
                )
            ):
                continue

            recurse(
                sequence
                + [
                    next_population
                ],
                outgoing_edges,
            )

    recurse(
        [
            EXPECTED_PARENT_TERMINALS_PER_END
        ],
        EXPECTED_PARENT_NEW_BONDS_PER_END,
    )

    filtered = [
        row
        for row in rows
        if (
            0.75
            * target_heavy_atoms
            <= int(
                row[
                    "total_added_heavy_atoms_per_end"
                ]
            )
            <= 1.25
            * target_heavy_atoms
        )
        and 6
        <= int(
            row[
                "inner_boundary_population"
            ]
        )
        <= 18
    ]

    filtered.sort(
        key=lambda row: (
            -float(
                row[
                    "sequence_screening_score"
                ]
            ),
            float(
                row[
                    "heavy_atom_relative_error"
                ]
            ),
            float(
                row[
                    "inner_spacing_relative_error"
                ]
            ),
            abs(
                int(
                    row[
                        "layer_count"
                    ]
                )
                - EXPECTED_SELECTED_LAYER_COUNT
            ),
            str(
                row[
                    "ring_populations"
                ]
            ),
        )
    )

    for rank, row in enumerate(
        filtered,
        start=1,
    ):
        row["rank"] = rank

    return filtered


def added_element(
    parent_element: str,
    layer_index_1based: int,
) -> str:
    if parent_element not in {
        "B",
        "N",
    }:
        raise RuntimeError(
            f"Unexpected parent element: {parent_element}"
        )

    if layer_index_1based % 2 == 1:
        return (
            "N"
            if parent_element == "B"
            else "B"
        )

    return parent_element


def connected_component_size(
    adjacency: dict[
        str,
        set[str],
    ],
    first_node: str,
) -> int:
    visited: set[str] = set()
    stack = [first_node]

    while stack:
        node = stack.pop()

        if node in visited:
            continue

        visited.add(node)

        stack.extend(
            adjacency.get(
                node,
                set(),
            )
            - visited
        )

    return len(visited)


def build_end_blueprint(
    end: str,
    parent_rows: list[
        dict[str, str]
    ],
    sequence: tuple[int, ...],
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, Any],
]:
    ordered_parent = sorted(
        parent_rows,
        key=lambda row: int(
            float(
                row[
                    "circumferential_order"
                ]
            )
        ),
    )

    if len(
        ordered_parent
    ) != EXPECTED_PARENT_TERMINALS_PER_END:
        raise RuntimeError(
            f"{end}: unexpected parent terminal count."
        )

    parent_elements = {
        row["element"]
        for row in ordered_parent
    }

    if len(parent_elements) != 1:
        raise RuntimeError(
            f"{end}: parent end is not elementally pure."
        )

    parent_element = next(
        iter(
            parent_elements
        )
    )

    expected_parent_element = (
        "B"
        if end == "LOWER"
        else "N"
    )

    if parent_element != expected_parent_element:
        raise RuntimeError(
            f"{end}: parent element "
            f"{parent_element}/{expected_parent_element}"
        )

    node_rows: list[
        dict[str, Any]
    ] = []

    edge_rows: list[
        dict[str, Any]
    ] = []

    interface_rows: list[
        dict[str, Any]
    ] = []

    node_element: dict[
        str,
        str
    ] = {}

    node_type: dict[
        str,
        str
    ] = {}

    parent_ids: list[str] = []

    for order, row in enumerate(
        ordered_parent
    ):
        parent_local_index = int(
            float(
                row[
                    "hbn_local_index_0based"
                ]
            )
        )

        node_id = (
            f"{end}:P:"
            f"{parent_local_index}"
        )

        parent_ids.append(
            node_id
        )

        node_element[node_id] = (
            parent_element
        )

        node_type[node_id] = (
            "PARENT_TERMINAL"
        )

        node_rows.append(
            {
                "end": end,
                "node_id": node_id,
                "node_type": (
                    "PARENT_TERMINAL"
                ),
                "layer": 0,
                "circumferential_index": (
                    order
                ),
                "element": parent_element,
                "parent_hbn_local_index_0based": (
                    parent_local_index
                ),
                "parent_existing_coordination": 1,
                "blueprint_new_coordination": 2,
                "target_total_coordination": 3,
                "formal_charge_assigned": False,
                "coordinates_assigned": False,
            }
        )

    previous_ids = parent_ids
    previous_outgoing_degrees = [
        2
        for _ in previous_ids
    ]

    for layer_index, population in enumerate(
        sequence,
        start=1,
    ):
        element = added_element(
            parent_element,
            layer_index,
        )

        current_ids = [
            (
                f"{end}:L{layer_index}:"
                f"{index}"
            )
            for index in range(
                population
            )
        ]

        total_interface_edges = sum(
            previous_outgoing_degrees
        )

        current_incoming_degrees = (
            balanced_degree_sequence(
                population,
                total_interface_edges,
            )
        )

        (
            pairs,
            maximum_span,
        ) = realize_local_interface(
            previous_outgoing_degrees,
            current_incoming_degrees,
        )

        observed_left = [0] * len(
            previous_ids
        )

        observed_right = [0] * len(
            current_ids
        )

        interface_name = (
            "PARENT_TO_L1"
            if layer_index == 1
            else (
                f"L{layer_index - 1}"
                f"_TO_L{layer_index}"
            )
        )

        for left_index, right_index in pairs:
            source_id = previous_ids[
                left_index
            ]

            target_id = current_ids[
                right_index
            ]

            source_element = (
                node_element[
                    source_id
                ]
            )

            target_element = element

            observed_left[
                left_index
            ] += 1

            observed_right[
                right_index
            ] += 1

            edge_rows.append(
                {
                    "end": end,
                    "edge_id": (
                        f"{end}:E:"
                        f"{len(edge_rows) + 1}"
                    ),
                    "interface": (
                        interface_name
                    ),
                    "source_node": (
                        source_id
                    ),
                    "target_node": (
                        target_id
                    ),
                    "source_element": (
                        source_element
                    ),
                    "target_element": (
                        target_element
                    ),
                    "edge_type": (
                        "PARENT_TO_COLLAR"
                        if layer_index == 1
                        else "COLLAR_BN"
                    ),
                    "heavy_atom_edge": True,
                    "heteropolar_BN_edge": (
                        {
                            source_element,
                            target_element,
                        }
                        == {
                            "B",
                            "N",
                        }
                    ),
                    "angular_span_turns": (
                        circular_distance_turns(
                            left_index,
                            len(
                                previous_ids
                            ),
                            right_index,
                            len(
                                current_ids
                            ),
                        )
                    ),
                    "formal_bond_order_assigned": False,
                    "coordinates_assigned": False,
                }
            )

        if (
            observed_left
            != previous_outgoing_degrees
        ):
            raise RuntimeError(
                f"{end} {interface_name}: "
                "left-degree mismatch."
            )

        if (
            observed_right
            != current_incoming_degrees
        ):
            raise RuntimeError(
                f"{end} {interface_name}: "
                "right-degree mismatch."
            )

        interface_rows.append(
            {
                "end": end,
                "interface": interface_name,
                "left_nodes": len(
                    previous_ids
                ),
                "right_nodes": len(
                    current_ids
                ),
                "edge_count": len(
                    pairs
                ),
                "left_degree_minimum": min(
                    previous_outgoing_degrees
                ),
                "left_degree_maximum": max(
                    previous_outgoing_degrees
                ),
                "right_degree_minimum": min(
                    current_incoming_degrees
                ),
                "right_degree_maximum": max(
                    current_incoming_degrees
                ),
                "maximum_angular_span_turns": (
                    maximum_span
                ),
                "locality_threshold_pass": (
                    maximum_span
                    <= MAX_INTERFACE_ANGULAR_SPAN_TURNS
                ),
            }
        )

        is_inner_layer = (
            layer_index
            == len(sequence)
        )

        if is_inner_layer:
            outgoing_degrees = [
                0
                for _ in current_ids
            ]
        else:
            outgoing_degrees = [
                3 - degree
                for degree
                in current_incoming_degrees
            ]

            if not all(
                degree in {
                    1,
                    2,
                }
                for degree in outgoing_degrees
            ):
                raise RuntimeError(
                    f"{end} layer {layer_index}: "
                    "invalid outgoing degree."
                )

        for index, node_id in enumerate(
            current_ids
        ):
            node_element[node_id] = (
                element
            )

            node_type[node_id] = (
                "INNER_BOUNDARY_BN"
                if is_inner_layer
                else "COLLAR_BN"
            )

            node_rows.append(
                {
                    "end": end,
                    "node_id": node_id,
                    "node_type": (
                        node_type[node_id]
                    ),
                    "layer": layer_index,
                    "circumferential_index": (
                        index
                    ),
                    "element": element,
                    "parent_hbn_local_index_0based": (
                        ""
                    ),
                    "incoming_blueprint_degree": (
                        current_incoming_degrees[
                            index
                        ]
                    ),
                    "outgoing_blueprint_degree": (
                        outgoing_degrees[
                            index
                        ]
                    ),
                    "target_total_coordination": 3,
                    "formal_charge_assigned": False,
                    "coordinates_assigned": False,
                }
            )

        previous_ids = current_ids
        previous_outgoing_degrees = (
            outgoing_degrees
        )

    inner_ids = previous_ids

    inner_bn_degrees = Counter()

    for row in edge_rows:
        if row[
            "heavy_atom_edge"
        ]:
            inner_bn_degrees[
                row[
                    "source_node"
                ]
            ] += 1

            inner_bn_degrees[
                row[
                    "target_node"
                ]
            ] += 1

    for inner_index, inner_id in enumerate(
        inner_ids
    ):
        if inner_bn_degrees[
            inner_id
        ] != 2:
            raise RuntimeError(
                f"{end}: inner BN node {inner_id} "
                "does not have heavy-atom degree 2."
            )

        hydrogen_id = (
            f"{end}:H:{inner_index}"
        )

        node_element[
            hydrogen_id
        ] = "H"

        node_type[
            hydrogen_id
        ] = "INNER_PASSIVANT_H"

        node_rows.append(
            {
                "end": end,
                "node_id": hydrogen_id,
                "node_type": (
                    "INNER_PASSIVANT_H"
                ),
                "layer": (
                    len(sequence) + 1
                ),
                "circumferential_index": (
                    inner_index
                ),
                "element": "H",
                "parent_hbn_local_index_0based": (
                    ""
                ),
                "incoming_blueprint_degree": 1,
                "outgoing_blueprint_degree": 0,
                "target_total_coordination": 1,
                "formal_charge_assigned": False,
                "coordinates_assigned": False,
            }
        )

        edge_rows.append(
            {
                "end": end,
                "edge_id": (
                    f"{end}:E:"
                    f"{len(edge_rows) + 1}"
                ),
                "interface": (
                    "INNER_BOUNDARY_TO_H"
                ),
                "source_node": inner_id,
                "target_node": (
                    hydrogen_id
                ),
                "source_element": (
                    node_element[
                        inner_id
                    ]
                ),
                "target_element": "H",
                "edge_type": (
                    "INNER_EDGE_PASSIVATION"
                ),
                "heavy_atom_edge": False,
                "heteropolar_BN_edge": False,
                "angular_span_turns": 0.0,
                "formal_bond_order_assigned": False,
                "coordinates_assigned": False,
            }
        )

    interface_rows.append(
        {
            "end": end,
            "interface": (
                "INNER_BOUNDARY_TO_H"
            ),
            "left_nodes": len(
                inner_ids
            ),
            "right_nodes": len(
                inner_ids
            ),
            "edge_count": len(
                inner_ids
            ),
            "left_degree_minimum": 1,
            "left_degree_maximum": 1,
            "right_degree_minimum": 1,
            "right_degree_maximum": 1,
            "maximum_angular_span_turns": 0.0,
            "locality_threshold_pass": True,
        }
    )

    adjacency: dict[
        str,
        set[str]
    ] = {
        row["node_id"]: set()
        for row in node_rows
    }

    graph_degrees = Counter()

    edge_pairs: list[
        tuple[str, str]
    ] = []

    for edge in edge_rows:
        first = str(
            edge[
                "source_node"
            ]
        )

        second = str(
            edge[
                "target_node"
            ]
        )

        pair = tuple(
            sorted(
                (
                    first,
                    second,
                )
            )
        )

        edge_pairs.append(
            pair
        )

        adjacency[first].add(
            second
        )

        adjacency[second].add(
            first
        )

        graph_degrees[first] += 1
        graph_degrees[second] += 1

    duplicate_edges = (
        len(edge_pairs)
        - len(
            set(
                edge_pairs
            )
        )
    )

    self_edges = sum(
        first == second
        for first, second
        in edge_pairs
    )

    added_bn_ids = [
        row["node_id"]
        for row in node_rows
        if row[
            "node_type"
        ]
        in {
            "COLLAR_BN",
            "INNER_BOUNDARY_BN",
        }
    ]

    hydrogen_ids = [
        row["node_id"]
        for row in node_rows
        if row[
            "node_type"
        ]
        == "INNER_PASSIVANT_H"
    ]

    parent_final_degree_failures = [
        node_id
        for node_id in parent_ids
        if (
            1
            + graph_degrees[
                node_id
            ]
        )
        != 3
    ]

    added_bn_degree_failures = [
        node_id
        for node_id in added_bn_ids
        if graph_degrees[
            node_id
        ]
        != 3
    ]

    hydrogen_degree_failures = [
        node_id
        for node_id in hydrogen_ids
        if graph_degrees[
            node_id
        ]
        != 1
    ]

    heavy_edges = [
        row
        for row in edge_rows
        if bool(
            row[
                "heavy_atom_edge"
            ]
        )
    ]

    same_element_heavy_edges = [
        row
        for row in heavy_edges
        if row[
            "source_element"
        ]
        == row[
            "target_element"
        ]
    ]

    nonheteropolar_heavy_edges = [
        row
        for row in heavy_edges
        if not bool(
            row[
                "heteropolar_BN_edge"
            ]
        )
    ]

    component_size = (
        connected_component_size(
            adjacency,
            node_rows[0][
                "node_id"
            ],
        )
    )

    composition = Counter(
        row["element"]
        for row in node_rows
        if row[
            "node_type"
        ]
        != "PARENT_TERMINAL"
    )

    maximum_interface_span = max(
        float(
            row[
                "maximum_angular_span_turns"
            ]
        )
        for row in interface_rows
    )

    parent_to_collar_edges = sum(
        row[
            "edge_type"
        ]
        == "PARENT_TO_COLLAR"
        for row in edge_rows
    )

    heavy_edge_count = len(
        heavy_edges
    )

    heavy_vertex_count = (
        len(parent_ids)
        + len(added_bn_ids)
    )

    heavy_cycle_rank = (
        heavy_edge_count
        - heavy_vertex_count
        + 1
    )

    metrics = {
        "end": end,
        "parent_element": (
            parent_element
        ),
        "first_layer_element": (
            added_element(
                parent_element,
                1,
            )
        ),
        "ring_populations": (
            "-".join(
                str(value)
                for value in sequence
            )
        ),
        "layer_count": len(
            sequence
        ),
        "parent_terminal_nodes": len(
            parent_ids
        ),
        "added_B_atoms": (
            composition.get(
                "B",
                0,
            )
        ),
        "added_N_atoms": (
            composition.get(
                "N",
                0,
            )
        ),
        "added_heavy_atoms": (
            composition.get(
                "B",
                0,
            )
            + composition.get(
                "N",
                0,
            )
        ),
        "added_H_atoms": (
            composition.get(
                "H",
                0,
            )
        ),
        "parent_to_collar_edges": (
            parent_to_collar_edges
        ),
        "heavy_atom_edges": (
            heavy_edge_count
        ),
        "passivation_edges": (
            len(edge_rows)
            - heavy_edge_count
        ),
        "total_blueprint_edges": len(
            edge_rows
        ),
        "same_element_heavy_edges": (
            len(
                same_element_heavy_edges
            )
        ),
        "nonheteropolar_heavy_edges": (
            len(
                nonheteropolar_heavy_edges
            )
        ),
        "duplicate_edges": (
            duplicate_edges
        ),
        "self_edges": self_edges,
        "connected_nodes": (
            component_size
        ),
        "total_nodes_in_end_graph": (
            len(node_rows)
        ),
        "is_connected": (
            component_size
            == len(node_rows)
        ),
        "parent_final_degree_failures": (
            len(
                parent_final_degree_failures
            )
        ),
        "added_BN_degree_failures": (
            len(
                added_bn_degree_failures
            )
        ),
        "hydrogen_degree_failures": (
            len(
                hydrogen_degree_failures
            )
        ),
        "maximum_interface_angular_span_turns": (
            maximum_interface_span
        ),
        "inner_boundary_population": (
            len(inner_ids)
        ),
        "heavy_graph_cycle_rank": (
            heavy_cycle_rank
        ),
    }

    return (
        node_rows,
        edge_rows,
        interface_rows,
        metrics,
    )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        PARENT_SUMMARY,
        TERMINAL_ATOMS,
        CANDIDATE_SELECTION,
        DESIGN_CONTRACT,
    ):
        require_file(required)

    parent = read_single_csv_row(
        PARENT_SUMMARY
    )

    selection = read_single_csv_row(
        CANDIDATE_SELECTION
    )

    terminal_rows = read_csv_rows(
        TERMINAL_ATOMS
    )

    contract = json.loads(
        DESIGN_CONTRACT.read_text(
            encoding="utf-8",
        )
    )

    if parent.get(
        "decision"
    ) != EXPECTED_GATE3A_DECISION:
        raise RuntimeError(
            "Gate 3A is not in the accepted state."
        )

    if selection.get(
        "decision"
    ) != EXPECTED_GATE3B_DECISION:
        raise RuntimeError(
            "Gate 3B is not in the accepted state."
        )

    if selection.get(
        "primary_candidate"
    ) != EXPECTED_PRIMARY_CANDIDATE:
        raise RuntimeError(
            "Unexpected Gate 3B primary candidate."
        )

    if contract.get(
        "primary_candidate_id"
    ) != EXPECTED_PRIMARY_CANDIDATE:
        raise RuntimeError(
            "Design contract and selection disagree."
        )

    if bool(
        contract.get(
            "explicit_coordinate_generation_authorized",
            True,
        )
    ):
        raise RuntimeError(
            "Upstream contract unexpectedly authorizes "
            "coordinate generation."
        )

    aperture_diameter_nm = parse_float(
        parent,
        "target_aperture_diameter_nm",
    )

    aperture_radius_nm = parse_float(
        parent,
        "target_aperture_radius_nm",
    )

    open_area_fraction = parse_float(
        parent,
        "target_open_area_fraction",
    )

    parent_rim_radius_nm = parse_float(
        parent,
        "parent_rim_mean_radius_nm",
    )

    target_heavy_atoms_per_end = parse_float(
        parent,
        "estimated_monolayer_hBN_atoms_per_end",
    )

    sequence_rows = (
        enumerate_population_sequences(
            target_heavy_atoms_per_end,
            parent_rim_radius_nm,
            aperture_radius_nm,
        )
    )

    if not sequence_rows:
        raise RuntimeError(
            "No admissible ring-population sequences "
            "were found."
        )

    selected_sequence_row = (
        sequence_rows[0]
    )

    selected_sequence = tuple(
        int(value)
        for value in str(
            selected_sequence_row[
                "ring_populations"
            ]
        ).split("-")
    )

    write_csv(
        SEQUENCES_CSV,
        sequence_rows,
    )

    lower_parent_rows = [
        row
        for row in terminal_rows
        if row.get(
            "end"
        )
        == "LOWER"
    ]

    upper_parent_rows = [
        row
        for row in terminal_rows
        if row.get(
            "end"
        )
        == "UPPER"
    ]

    (
        lower_nodes,
        lower_edges,
        lower_interfaces,
        lower_metrics,
    ) = build_end_blueprint(
        "LOWER",
        lower_parent_rows,
        selected_sequence,
    )

    (
        upper_nodes,
        upper_edges,
        upper_interfaces,
        upper_metrics,
    ) = build_end_blueprint(
        "UPPER",
        upper_parent_rows,
        selected_sequence,
    )

    all_nodes = (
        lower_nodes
        + upper_nodes
    )

    all_edges = (
        lower_edges
        + upper_edges
    )

    all_interfaces = (
        lower_interfaces
        + upper_interfaces
    )

    write_csv(
        NODES_CSV,
        all_nodes,
    )

    write_csv(
        EDGES_CSV,
        all_edges,
    )

    write_csv(
        INTERFACES_CSV,
        all_interfaces,
    )

    added_composition = Counter(
        row["element"]
        for row in all_nodes
        if row[
            "node_type"
        ]
        != "PARENT_TERMINAL"
    )

    total_added_heavy_atoms = (
        added_composition.get(
            "B",
            0,
        )
        + added_composition.get(
            "N",
            0,
        )
    )

    total_added_hydrogens = (
        added_composition.get(
            "H",
            0,
        )
    )

    heavy_atom_relative_error = float(
        selected_sequence_row[
            "heavy_atom_relative_error"
        ]
    )

    inner_spacing_relative_error = float(
        selected_sequence_row[
            "inner_spacing_relative_error"
        ]
    )

    maximum_interface_span = max(
        float(
            row[
                "maximum_angular_span_turns"
            ]
        )
        for row in all_interfaces
    )

    total_parent_to_collar_edges = (
        int(
            lower_metrics[
                "parent_to_collar_edges"
            ]
        )
        + int(
            upper_metrics[
                "parent_to_collar_edges"
            ]
        )
    )

    total_same_element_heavy_edges = (
        int(
            lower_metrics[
                "same_element_heavy_edges"
            ]
        )
        + int(
            upper_metrics[
                "same_element_heavy_edges"
            ]
        )
    )

    total_degree_failures = sum(
        int(
            metrics[field]
        )
        for metrics in (
            lower_metrics,
            upper_metrics,
        )
        for field in (
            "parent_final_degree_failures",
            "added_BN_degree_failures",
            "hydrogen_degree_failures",
        )
    )

    gates = {
        "Gate3A_parent_audit_is_accepted": (
            parent.get(
                "decision"
            )
            == EXPECTED_GATE3A_DECISION
        ),
        "Gate3B_candidate_ranking_is_accepted": (
            selection.get(
                "decision"
            )
            == EXPECTED_GATE3B_DECISION
        ),
        "primary_candidate_is_C2_graded_heteropolar_BN_collar": (
            selection.get(
                "primary_candidate"
            )
            == EXPECTED_PRIMARY_CANDIDATE
        ),
        "selected_ring_population_sequence_is_expected": (
            selected_sequence
            == EXPECTED_RING_POPULATIONS
        ),
        "selected_blueprint_has_six_layers": (
            len(
                selected_sequence
            )
            == EXPECTED_SELECTED_LAYER_COUNT
        ),
        "selected_blueprint_has_144_added_heavy_atoms_per_end": (
            int(
                lower_metrics[
                    "added_heavy_atoms"
                ]
            )
            == EXPECTED_ADDED_HEAVY_ATOMS_PER_END
            and int(
                upper_metrics[
                    "added_heavy_atoms"
                ]
            )
            == EXPECTED_ADDED_HEAVY_ATOMS_PER_END
        ),
        "total_added_heavy_atoms_match_288_R2_steric_beads": (
            total_added_heavy_atoms
            == EXPECTED_ADDED_HEAVY_ATOMS_TOTAL
            == EXPECTED_STERIC_BEADS_TOTAL
        ),
        "heavy_atom_estimate_relative_error_is_within_2_percent": (
            heavy_atom_relative_error
            <= MAX_HEAVY_ATOM_ESTIMATE_RELATIVE_ERROR
        ),
        "inner_boundary_spacing_error_is_within_15_percent": (
            inner_spacing_relative_error
            <= MAX_INNER_SPACING_RELATIVE_ERROR
        ),
        "lower_parent_first_layer_is_B_to_N": (
            lower_metrics[
                "parent_element"
            ]
            == "B"
            and lower_metrics[
                "first_layer_element"
            ]
            == "N"
        ),
        "upper_parent_first_layer_is_N_to_B": (
            upper_metrics[
                "parent_element"
            ]
            == "N"
            and upper_metrics[
                "first_layer_element"
            ]
            == "B"
        ),
        "parent_to_collar_bond_count_is_60_per_end": (
            int(
                lower_metrics[
                    "parent_to_collar_edges"
                ]
            )
            == EXPECTED_PARENT_NEW_BONDS_PER_END
            and int(
                upper_metrics[
                    "parent_to_collar_edges"
                ]
            )
            == EXPECTED_PARENT_NEW_BONDS_PER_END
        ),
        "parent_to_collar_bond_count_is_120_total": (
            total_parent_to_collar_edges
            == EXPECTED_PARENT_NEW_BONDS_TOTAL
        ),
        "all_parent_and_added_atom_coordination_targets_are_met": (
            total_degree_failures == 0
        ),
        "all_heavy_atom_edges_are_heteropolar_BN": (
            total_same_element_heavy_edges
            == 0
            and int(
                lower_metrics[
                    "nonheteropolar_heavy_edges"
                ]
            )
            == 0
            and int(
                upper_metrics[
                    "nonheteropolar_heavy_edges"
                ]
            )
            == 0
        ),
        "lower_end_graph_is_connected": (
            bool(
                lower_metrics[
                    "is_connected"
                ]
            )
        ),
        "upper_end_graph_is_connected": (
            bool(
                upper_metrics[
                    "is_connected"
                ]
            )
        ),
        "blueprint_has_no_duplicate_edges": (
            int(
                lower_metrics[
                    "duplicate_edges"
                ]
            )
            == 0
            and int(
                upper_metrics[
                    "duplicate_edges"
                ]
            )
            == 0
        ),
        "blueprint_has_no_self_edges": (
            int(
                lower_metrics[
                    "self_edges"
                ]
            )
            == 0
            and int(
                upper_metrics[
                    "self_edges"
                ]
            )
            == 0
        ),
        "interface_angular_locality_is_within_threshold": (
            maximum_interface_span
            <= MAX_INTERFACE_ANGULAR_SPAN_TURNS
        ),
        "inner_boundary_has_12_H_passivants_per_end": (
            int(
                lower_metrics[
                    "added_H_atoms"
                ]
            )
            == EXPECTED_INNER_PASSIVANTS_PER_END
            and int(
                upper_metrics[
                    "added_H_atoms"
                ]
            )
            == EXPECTED_INNER_PASSIVANTS_PER_END
        ),
        "total_inner_passivants_are_24": (
            total_added_hydrogens
            == EXPECTED_INNER_PASSIVANTS_TOTAL
        ),
        "combined_added_BN_composition_is_balanced": (
            added_composition.get(
                "B",
                0,
            )
            == 144
            and added_composition.get(
                "N",
                0,
            )
            == 144
        ),
        "no_formal_charges_were_assigned": all(
            not bool(
                row[
                    "formal_charge_assigned"
                ]
            )
            for row in all_nodes
        ),
        "no_coordinates_were_assigned": all(
            not bool(
                row[
                    "coordinates_assigned"
                ]
            )
            for row in all_nodes
        )
        and all(
            not bool(
                row[
                    "coordinates_assigned"
                ]
            )
            for row in all_edges
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
        "R2_GRADED_HETEROPOLAR_BN_COLLAR_CONNECTIVITY_BLUEPRINT_REQUIRES_REVIEW"
    )

    required_next_step = (
        "BUILD_AND_VALIDATE_R2_GRADED_HETEROPOLAR_COLLAR_STATIC_COORDINATE_EMBEDDING"
        if accepted
        else
        "REVIEW_R2_GRADED_COLLAR_CONNECTIVITY_BLUEPRINT_FAILURES"
    )

    summary = {
        "decision": decision,
        "primary_candidate": (
            EXPECTED_PRIMARY_CANDIDATE
        ),
        "primary_candidate_is_final_chemistry": False,
        "selected_sequence_rank": (
            selected_sequence_row[
                "rank"
            ]
        ),
        "selected_ring_populations": (
            selected_sequence_row[
                "ring_populations"
            ]
        ),
        "selected_layer_count": len(
            selected_sequence
        ),
        "selected_sequence_score": (
            selected_sequence_row[
                "sequence_screening_score"
            ]
        ),
        "target_heavy_atoms_per_end": (
            target_heavy_atoms_per_end
        ),
        "selected_heavy_atoms_per_end": (
            sum(
                selected_sequence
            )
        ),
        "heavy_atom_relative_error": (
            heavy_atom_relative_error
        ),
        "inner_boundary_population": (
            selected_sequence[-1]
        ),
        "inner_boundary_spacing_nm": (
            selected_sequence_row[
                "inner_boundary_spacing_nm"
            ]
        ),
        "parent_terminal_spacing_nm": (
            selected_sequence_row[
                "parent_terminal_spacing_nm"
            ]
        ),
        "inner_spacing_relative_error": (
            inner_spacing_relative_error
        ),
        "target_aperture_diameter_nm": (
            aperture_diameter_nm
        ),
        "target_aperture_radius_nm": (
            aperture_radius_nm
        ),
        "target_open_area_fraction": (
            open_area_fraction
        ),
        "lower_added_B_atoms": (
            lower_metrics[
                "added_B_atoms"
            ]
        ),
        "lower_added_N_atoms": (
            lower_metrics[
                "added_N_atoms"
            ]
        ),
        "lower_added_H_atoms": (
            lower_metrics[
                "added_H_atoms"
            ]
        ),
        "upper_added_B_atoms": (
            upper_metrics[
                "added_B_atoms"
            ]
        ),
        "upper_added_N_atoms": (
            upper_metrics[
                "added_N_atoms"
            ]
        ),
        "upper_added_H_atoms": (
            upper_metrics[
                "added_H_atoms"
            ]
        ),
        "combined_added_B_atoms": (
            added_composition.get(
                "B",
                0,
            )
        ),
        "combined_added_N_atoms": (
            added_composition.get(
                "N",
                0,
            )
        ),
        "combined_added_H_atoms": (
            added_composition.get(
                "H",
                0,
            )
        ),
        "total_parent_to_collar_edges": (
            total_parent_to_collar_edges
        ),
        "total_same_element_heavy_edges": (
            total_same_element_heavy_edges
        ),
        "total_coordination_failures": (
            total_degree_failures
        ),
        "maximum_interface_angular_span_turns": (
            maximum_interface_span
        ),
        "lower_end_connected": (
            lower_metrics[
                "is_connected"
            ]
        ),
        "upper_end_connected": (
            upper_metrics[
                "is_connected"
            ]
        ),
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
        "QM_recalculation_authorized": False,
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
        BLUEPRINT_SUMMARY_CSV,
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

    blueprint_payload = {
        "summary": summary,
        "lower_end_metrics": (
            lower_metrics
        ),
        "upper_end_metrics": (
            upper_metrics
        ),
        "gates": gates,
        "selected_sequence": list(
            selected_sequence
        ),
        "limitations": [
            (
                "The graph is not a three-dimensional "
                "coordinate embedding."
            ),
            (
                "No bond lengths, bond angles, dihedrals, "
                "formal bond orders, or partial charges "
                "have been assigned."
            ),
            (
                "Graph feasibility does not establish "
                "energetic stability or synthetic realizability."
            ),
            (
                "The target aperture is represented by "
                "the selected inner-boundary population "
                "and inherited geometric constraints only."
            ),
        ],
    }

    BLUEPRINT_JSON.write_text(
        json.dumps(
            blueprint_payload,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    source_rows = [
        {
            "role": (
                "Gate3A_parent_summary"
            ),
            "file": relative(
                PARENT_SUMMARY
            ),
            "sha256": sha256(
                PARENT_SUMMARY
            ),
        },
        {
            "role": (
                "Gate3A_terminal_atoms"
            ),
            "file": relative(
                TERMINAL_ATOMS
            ),
            "sha256": sha256(
                TERMINAL_ATOMS
            ),
        },
        {
            "role": (
                "Gate3B_candidate_selection"
            ),
            "file": relative(
                CANDIDATE_SELECTION
            ),
            "sha256": sha256(
                CANDIDATE_SELECTION
            ),
        },
        {
            "role": (
                "Gate3B_design_contract"
            ),
            "file": relative(
                DESIGN_CONTRACT
            ),
            "sha256": sha256(
                DESIGN_CONTRACT
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
        for name, passed in gates.items()
    )

    interface_lines = "\n".join(
        (
            f"- {row['end']} `{row['interface']}`: "
            f"{row['left_nodes']} → {row['right_nodes']} nodes; "
            f"{row['edge_count']} edges; "
            f"maximum angular span "
            f"{float(row['maximum_angular_span_turns']):.6f} turns"
        )
        for row in all_interfaces
    )

    REPORT_MD.write_text(
        f"""# R2 Graded Heteropolar BN Collar Connectivity Blueprint

## Scope

This gate constructs and validates a coordinate-free graph blueprint
for the selected C2 graded heteropolar BN collar–annulus candidate.

No three-dimensional coordinates, molecular topology, bond-order
model, partial charges, force-field parameters, minimization, MD, or
QM calculation were generated.

## Selected ring-population sequence

- Sequence:
  **{selected_sequence_row['ring_populations']}**
- Layers:
  **{len(selected_sequence)}**
- Added B/N atoms per end:
  **{sum(selected_sequence)}**
- Screening estimate:
  **{target_heavy_atoms_per_end:.3f} atoms/end**
- Relative heavy-atom error:
  **{heavy_atom_relative_error:.6f}**
- Inner-boundary population:
  **{selected_sequence[-1]}**
- Inner-boundary spacing proxy:
  **{float(selected_sequence_row['inner_boundary_spacing_nm']):.6f} nm**
- Parent terminal-site spacing:
  **{float(selected_sequence_row['parent_terminal_spacing_nm']):.6f} nm**
- Relative spacing error:
  **{inner_spacing_relative_error:.6f}**

The selected graph contains exactly 144 added B/N atoms per end,
matching the 144 steric R2 beads per end.

## End-specific composition

### Lower B-terminated end

- First added layer:
  **N**
- Added B/N/H:
  **{lower_metrics['added_B_atoms']}/
  {lower_metrics['added_N_atoms']}/
  {lower_metrics['added_H_atoms']}**
- Parent-to-collar edges:
  **{lower_metrics['parent_to_collar_edges']}**
- Same-element heavy edges:
  **{lower_metrics['same_element_heavy_edges']}**
- Connected:
  **{lower_metrics['is_connected']}**

### Upper N-terminated end

- First added layer:
  **B**
- Added B/N/H:
  **{upper_metrics['added_B_atoms']}/
  {upper_metrics['added_N_atoms']}/
  {upper_metrics['added_H_atoms']}**
- Parent-to-collar edges:
  **{upper_metrics['parent_to_collar_edges']}**
- Same-element heavy edges:
  **{upper_metrics['same_element_heavy_edges']}**
- Connected:
  **{upper_metrics['is_connected']}**

### Combined added structure

- B/N/H:
  **{added_composition.get('B', 0)}/
  {added_composition.get('N', 0)}/
  {added_composition.get('H', 0)}**
- Added heavy atoms:
  **{total_added_heavy_atoms}**
- Parent-to-collar bonds:
  **{total_parent_to_collar_edges}**
- Coordination failures:
  **{total_degree_failures}**
- Same-element heavy edges:
  **{total_same_element_heavy_edges}**

## Interfaces

{interface_lines}

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

## Interpretation limitation

The graph demonstrates that the parent valence deficits, end
asymmetry, heavy-atom population, connectivity, and heteropolar
bonding constraints can be satisfied simultaneously at the abstract
topological level. It does not establish that the graph can be
embedded in three dimensions with chemically acceptable B–N and X–H
bond lengths, bond angles, strain, planarity, aperture size, or
energetic stability.
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 graded heteropolar BN collar "
        "connectivity blueprint completed."
    )

    print(
        "Feasible population sequences retained: "
        f"{len(sequence_rows)}"
    )

    print(
        "Selected sequence rank / populations / score: "
        f"{selected_sequence_row['rank']}/"
        f"{selected_sequence_row['ring_populations']}/"
        f"{float(selected_sequence_row['sequence_screening_score']):.6f}"
    )

    print(
        "Selected layers / added heavy atoms per end: "
        f"{len(selected_sequence)}/"
        f"{sum(selected_sequence)}"
    )

    print(
        "Target / selected heavy atoms per end / "
        "relative error: "
        f"{target_heavy_atoms_per_end:.3f}/"
        f"{sum(selected_sequence)}/"
        f"{heavy_atom_relative_error:.6f}"
    )

    print(
        "Inner population / parent spacing / "
        "inner spacing / relative error: "
        f"{selected_sequence[-1]}/"
        f"{float(selected_sequence_row['parent_terminal_spacing_nm']):.6f}/"
        f"{float(selected_sequence_row['inner_boundary_spacing_nm']):.6f}/"
        f"{inner_spacing_relative_error:.6f}"
    )

    print(
        "Lower parent / first layer / added B/N/H: "
        f"{lower_metrics['parent_element']}/"
        f"{lower_metrics['first_layer_element']}/"
        f"{lower_metrics['added_B_atoms']}/"
        f"{lower_metrics['added_N_atoms']}/"
        f"{lower_metrics['added_H_atoms']}"
    )

    print(
        "Upper parent / first layer / added B/N/H: "
        f"{upper_metrics['parent_element']}/"
        f"{upper_metrics['first_layer_element']}/"
        f"{upper_metrics['added_B_atoms']}/"
        f"{upper_metrics['added_N_atoms']}/"
        f"{upper_metrics['added_H_atoms']}"
    )

    print(
        "Combined added B/N/H: "
        f"{added_composition.get('B', 0)}/"
        f"{added_composition.get('N', 0)}/"
        f"{added_composition.get('H', 0)}"
    )

    print(
        "Parent-to-collar bonds lower/upper/total: "
        f"{lower_metrics['parent_to_collar_edges']}/"
        f"{upper_metrics['parent_to_collar_edges']}/"
        f"{total_parent_to_collar_edges}"
    )

    print(
        "Coordination failures: "
        f"{total_degree_failures}"
    )

    print(
        "Same-element heavy edges: "
        f"{total_same_element_heavy_edges}"
    )

    print(
        "Lower / upper connected: "
        f"{lower_metrics['is_connected']}/"
        f"{upper_metrics['is_connected']}"
    )

    print(
        "Maximum interface angular span: "
        f"{maximum_interface_span:.6f} turns"
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
        SEQUENCES_CSV,
        NODES_CSV,
        EDGES_CSV,
        INTERFACES_CSV,
        BLUEPRINT_SUMMARY_CSV,
        GATES_CSV,
        BLUEPRINT_JSON,
        SOURCE_MANIFEST_CSV,
        REPORT_MD,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not accepted:
        raise RuntimeError(
            "R2 graded collar connectivity blueprint "
            "requires review."
        )


if __name__ == "__main__":
    main()
