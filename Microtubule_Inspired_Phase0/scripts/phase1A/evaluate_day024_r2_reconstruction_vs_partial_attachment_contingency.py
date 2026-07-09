#!/usr/bin/env python3

from __future__ import annotations

import cmath
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

GATE3D_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "05_r2_hexagonal_edge_completion_seed"
)

GATE3E_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "06_r2_hexagonal_annulus_attachment_feasibility"
)

OUTPUT_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "07_r2_reconstruction_vs_partial_attachment_contingency"
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

SEED_NODES_CSV = (
    GATE3D_ROOT
    / "r2_hexagonal_edge_completion_nodes.csv"
)

SEED_EDGES_CSV = (
    GATE3D_ROOT
    / "r2_hexagonal_edge_completion_edges.csv"
)

SEED_SUMMARY_CSV = (
    GATE3D_ROOT
    / "r2_hexagonal_edge_completion_seed_summary.csv"
)

ANNULUS_FEASIBILITY_SUMMARY_CSV = (
    GATE3E_ROOT
    / "r2_hexagonal_annulus_attachment_feasibility_summary.csv"
)

RECONSTRUCTION_SCREEN_CSV = (
    OUTPUT_ROOT
    / "r2_parent_full_shell_reconstruction_screen.csv"
)

MAPPING_SEARCH_CSV = (
    OUTPUT_ROOT
    / "r2_partial_attachment_mapping_search.csv"
)

SELECTED_NODES_CSV = (
    OUTPUT_ROOT
    / "r2_partial_attachment_passivated_annulus_nodes.csv"
)

SELECTED_EDGES_CSV = (
    OUTPUT_ROOT
    / "r2_partial_attachment_passivated_annulus_edges.csv"
)

END_SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_partial_attachment_passivated_annulus_end_summary.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_reconstruction_vs_partial_attachment_summary.csv"
)

GATES_CSV = (
    OUTPUT_ROOT
    / "r2_reconstruction_vs_partial_attachment_gates.csv"
)

DESIGN_JSON = (
    OUTPUT_ROOT
    / "r2_partial_attachment_passivated_annulus_design.json"
)

SOURCE_MANIFEST_CSV = (
    OUTPUT_ROOT
    / "r2_reconstruction_vs_partial_attachment_source_manifest.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_RECONSTRUCTION_VS_PARTIAL_ATTACHMENT_CONTINGENCY_DAY024.md"
)

EXPECTED_PARENT_DECISION = (
    "R2_PARENT_RIM_AND_CHEMICAL_CONSTRAINTS_AUDITED"
)

EXPECTED_SEED_DECISION = (
    "R2_HEXAGONAL_EDGE_COMPLETION_SEED_VALIDATED"
)

EXPECTED_ANNULUS_DECISION = (
    "R2_PURE_HEXAGONAL_BN_ANNULUS_DIRECT_ATTACHMENT_"
    "NOT_FEASIBLE_WITH_CURRENT_HOMOPOLAR_DEGREE2_RIM"
)

PASS_DECISION = (
    "R2_PARTIAL_HETEROPOLAR_ANNULUS_ATTACHMENT_AND_"
    "COMPLEMENTARY_PASSIVATION_GRAPH_VALIDATED"
)

EXPECTED_PARENT_ATOMS = 1680
EXPECTED_PARENT_BONDS = 2460

EXPECTED_TERMINALS_PER_END = 30
EXPECTED_SEED_ATOMS_PER_END = 30

EXPECTED_DIRECT_ATTACHMENTS_PER_END = 15
EXPECTED_UNATTACHED_SEED_SITES_PER_END = 15
EXPECTED_UNATTACHED_OUTER_SITES_PER_END = 15

EXPECTED_ANNULUS_HEAVY_ATOMS = 126
EXPECTED_TOTAL_HEAVY_ATOMS_PER_END = 156

MAX_HEAVY_ATOM_RELATIVE_ERROR = 0.15
MAX_OUTER_RADIUS_RELATIVE_ERROR = 0.15
MAX_INNER_RADIUS_RELATIVE_ERROR = 0.30

MAX_MAPPING_ANGULAR_RESIDUAL_TURNS = 0.10
MIN_ATTACHMENT_CREATED_CYCLE_LENGTH = 6

MAX_RECONSTRUCTION_DEPTH = 4

VERTEX_OFFSETS = (
    (2, 0),
    (1, 1),
    (-1, 1),
    (-2, 0),
    (-1, -1),
    (1, -1),
)


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


def cell_centers(
    shell_index: int,
) -> list[tuple[int, int]]:
    radius = shell_index - 1
    cells = []

    for q in range(
        -radius,
        radius + 1,
    ):
        for r in range(
            -radius,
            radius + 1,
        ):
            s = -q - r

            if max(
                abs(q),
                abs(r),
                abs(s),
            ) <= radius:
                cells.append(
                    (
                        q,
                        r,
                    )
                )

    return cells


def build_hexagonal_flake(
    shell_index: int,
) -> tuple[
    set[tuple[int, int]],
    set[
        tuple[
            tuple[int, int],
            tuple[int, int],
        ]
    ],
]:
    vertices: set[
        tuple[int, int]
    ] = set()

    edges: set[
        tuple[
            tuple[int, int],
            tuple[int, int],
        ]
    ] = set()

    for q, r in cell_centers(
        shell_index
    ):
        center_x = 3 * q
        center_y = q + 2 * r

        ring = [
            (
                center_x + dx,
                center_y + dy,
            )
            for dx, dy in VERTEX_OFFSETS
        ]

        vertices.update(
            ring
        )

        for index in range(6):
            first = ring[index]
            second = ring[
                (index + 1) % 6
            ]

            edges.add(
                tuple(
                    sorted(
                        (
                            first,
                            second,
                        )
                    )
                )
            )

    return vertices, edges


def build_adjacency(
    nodes: set[Any],
    edges: set[tuple[Any, Any]],
) -> dict[Any, set[Any]]:
    adjacency = {
        node: set()
        for node in nodes
    }

    for first, second in edges:
        adjacency[first].add(
            second
        )

        adjacency[second].add(
            first
        )

    return adjacency


def connected_components(
    adjacency: dict[Any, set[Any]],
) -> list[set[Any]]:
    remaining = set(adjacency)
    components = []

    while remaining:
        start = next(iter(remaining))
        component = set()
        queue: deque[Any] = deque(
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

        components.append(
            component
        )

        remaining -= component

    return components


def bipartite_coloring(
    adjacency: dict[Any, set[Any]],
) -> tuple[
    bool,
    dict[Any, int],
]:
    colors = {}

    for start in adjacency:
        if start in colors:
            continue

        colors[start] = 0

        queue: deque[Any] = deque(
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


def count_four_cycles(
    adjacency: dict[Any, set[Any]],
) -> int:
    nodes = list(adjacency)
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
            "Four-cycle count was not divisible by two."
        )

    return raw_count // 2


def graph_girth(
    adjacency: dict[Any, set[Any]],
) -> int:
    best = math.inf

    for start in adjacency:
        distances = {
            start: 0
        }

        parents = {
            start: None
        }

        queue: deque[Any] = deque(
            [start]
        )

        while queue:
            node = queue.popleft()

            if (
                distances[node] * 2 + 1
                >= best
            ):
                continue

            for neighbor in adjacency[node]:
                if neighbor not in distances:
                    distances[neighbor] = (
                        distances[node] + 1
                    )

                    parents[neighbor] = node

                    queue.append(
                        neighbor
                    )

                elif (
                    parents[node] != neighbor
                    and parents.get(
                        neighbor
                    )
                    != node
                ):
                    cycle_length = (
                        distances[node]
                        + distances[neighbor]
                        + 1
                    )

                    best = min(
                        best,
                        cycle_length,
                    )

    return (
        0
        if math.isinf(best)
        else int(best)
    )


def shortest_path_lengths(
    adjacency: dict[Any, set[Any]],
    source: Any,
) -> dict[Any, int]:
    distances = {
        source: 0
    }

    queue: deque[Any] = deque(
        [source]
    )

    while queue:
        node = queue.popleft()

        for neighbor in adjacency[node]:
            if neighbor in distances:
                continue

            distances[neighbor] = (
                distances[node] + 1
            )

            queue.append(
                neighbor
            )

    return distances


def multi_source_distances(
    adjacency: dict[int, set[int]],
    sources: list[int],
) -> dict[int, int]:
    distances = {
        source: 0
        for source in sources
    }

    queue: deque[int] = deque(
        sources
    )

    while queue:
        node = queue.popleft()

        for neighbor in adjacency[node]:
            if neighbor in distances:
                continue

            distances[neighbor] = (
                distances[node] + 1
            )

            queue.append(
                neighbor
            )

    return distances


def circular_distance_turns(
    first: float,
    second: float,
) -> float:
    difference = abs(
        first - second
    ) % 1.0

    return min(
        difference,
        1.0 - difference,
    )


def physical_xy(
    vertex: tuple[int, int],
    bond_length_nm: float,
) -> tuple[float, float]:
    integer_x, integer_y = vertex

    return (
        integer_x
        * bond_length_nm
        / 2.0,
        integer_y
        * math.sqrt(3.0)
        * bond_length_nm
        / 2.0,
    )


def angle_turns(
    vertex: tuple[int, int],
    bond_length_nm: float,
) -> float:
    x_nm, y_nm = physical_xy(
        vertex,
        bond_length_nm,
    )

    return (
        math.atan2(
            y_nm,
            x_nm,
        )
        / (
            2.0
            * math.pi
        )
    ) % 1.0


def radius_nm(
    vertex: tuple[int, int],
    bond_length_nm: float,
) -> float:
    x_nm, y_nm = physical_xy(
        vertex,
        bond_length_nm,
    )

    return math.hypot(
        x_nm,
        y_nm,
    )


def optimal_phase_residual(
    source_angles: list[float],
    target_angles: list[float],
) -> tuple[
    float,
    float,
]:
    if len(source_angles) != len(
        target_angles
    ):
        raise RuntimeError(
            "Angular lists have different sizes."
        )

    phase_vector = sum(
        cmath.exp(
            2j
            * math.pi
            * (
                target
                - source
            )
        )
        for source, target
        in zip(
            source_angles,
            target_angles,
        )
    )

    phase = (
        cmath.phase(
            phase_vector
        )
        / (
            2.0
            * math.pi
        )
    ) % 1.0

    residuals = [
        circular_distance_turns(
            (
                source + phase
            ) % 1.0,
            target,
        )
        for source, target
        in zip(
            source_angles,
            target_angles,
        )
    ]

    return (
        phase,
        max(residuals),
    )


def add_graph_edge(
    adjacency: dict[str, set[str]],
    edge_set: set[tuple[str, str]],
    first: str,
    second: str,
) -> None:
    if first == second:
        raise RuntimeError(
            f"Self-edge requested for {first}"
        )

    pair = tuple(
        sorted(
            (
                first,
                second,
            )
        )
    )

    if pair in edge_set:
        raise RuntimeError(
            f"Duplicate edge requested: {pair}"
        )

    edge_set.add(pair)

    adjacency.setdefault(
        first,
        set(),
    ).add(second)

    adjacency.setdefault(
        second,
        set(),
    ).add(first)


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    required_files = (
        PARENT_ATOMS_CSV,
        PARENT_BONDS_CSV,
        TERMINAL_ATOMS_CSV,
        PARENT_SUMMARY_CSV,
        SEED_NODES_CSV,
        SEED_EDGES_CSV,
        SEED_SUMMARY_CSV,
        ANNULUS_FEASIBILITY_SUMMARY_CSV,
    )

    for required in required_files:
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

    seed_nodes = read_csv_rows(
        SEED_NODES_CSV
    )

    seed_edges = read_csv_rows(
        SEED_EDGES_CSV
    )

    parent_summary = read_single_csv_row(
        PARENT_SUMMARY_CSV
    )

    seed_summary = read_single_csv_row(
        SEED_SUMMARY_CSV
    )

    annulus_summary = read_single_csv_row(
        ANNULUS_FEASIBILITY_SUMMARY_CSV
    )

    if parent_summary.get(
        "decision"
    ) != EXPECTED_PARENT_DECISION:
        raise RuntimeError(
            "Gate 3A parent audit is not accepted."
        )

    if seed_summary.get(
        "decision"
    ) != EXPECTED_SEED_DECISION:
        raise RuntimeError(
            "Gate 3D edge-completion seed is not accepted."
        )

    if annulus_summary.get(
        "decision"
    ) != EXPECTED_ANNULUS_DECISION:
        raise RuntimeError(
            "Gate 3E does not contain the expected "
            "negative direct-attachment decision."
        )

    if len(parent_atoms) != EXPECTED_PARENT_ATOMS:
        raise RuntimeError(
            "Unexpected parent atom count."
        )

    if len(parent_bonds) != EXPECTED_PARENT_BONDS:
        raise RuntimeError(
            "Unexpected parent bond count."
        )

    bond_length_nm = parse_float(
        parent_summary,
        "BN_bond_mean_nm",
    )

    parent_rim_radius_nm = parse_float(
        parent_summary,
        "parent_rim_mean_radius_nm",
    )

    target_aperture_radius_nm = parse_float(
        parent_summary,
        "target_aperture_radius_nm",
    )

    target_aperture_diameter_nm = parse_float(
        parent_summary,
        "target_aperture_diameter_nm",
    )

    target_heavy_atoms_per_end = parse_float(
        parent_summary,
        "estimated_monolayer_hBN_atoms_per_end",
    )

    outer_shell = parse_int(
        annulus_summary,
        "best_geometry_outer_shell",
    )

    inner_shell = parse_int(
        annulus_summary,
        "best_geometry_inner_shell",
    )

    if (
        outer_shell != 5
        or inner_shell != 2
    ):
        raise RuntimeError(
            "Unexpected best geometry template: "
            f"{outer_shell}/{inner_shell}"
        )

    parent_elements: dict[
        int,
        str
    ] = {}

    parent_adjacency: dict[
        int,
        set[int]
    ] = {
        index: set()
        for index in range(
            EXPECTED_PARENT_ATOMS
        )
    }

    for row in parent_atoms:
        index = int(
            float(
                row[
                    "hbn_local_index_0based"
                ]
            )
        )

        parent_elements[index] = (
            row[
                "element"
            ]
        )

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
                "Duplicate parent edge."
            )

        parent_edge_set.add(pair)

        parent_adjacency[first].add(
            second
        )

        parent_adjacency[second].add(
            first
        )

    terminal_by_end: dict[
        str,
        list[int]
    ] = {}

    for end in (
        "LOWER",
        "UPPER",
    ):
        rows = [
            row
            for row in terminal_rows
            if row.get(
                "end"
            )
            == end
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

        terminal_by_end[end] = [
            int(
                float(
                    row[
                        "hbn_local_index_0based"
                    ]
                )
            )
            for row in rows
        ]

    reconstruction_rows: list[
        dict[str, Any]
    ] = []

    reconstruction_candidates = []

    for end in (
        "LOWER",
        "UPPER",
    ):
        distances = multi_source_distances(
            parent_adjacency,
            terminal_by_end[end],
        )

        for depth in range(
            0,
            MAX_RECONSTRUCTION_DEPTH + 1,
        ):
            removed = {
                node
                for node, distance
                in distances.items()
                if distance < depth
            }

            frontier = sorted(
                node
                for node, distance
                in distances.items()
                if distance == depth
            )

            remaining_degrees = {
                node: sum(
                    neighbor not in removed
                    for neighbor in parent_adjacency[
                        node
                    ]
                )
                for node in frontier
            }

            b_count = sum(
                parent_elements[node] == "B"
                for node in frontier
            )

            n_count = sum(
                parent_elements[node] == "N"
                for node in frontier
            )

            degree_counts = Counter(
                remaining_degrees.values()
            )

            simple_mixed_degree2_rim = (
                len(frontier) == 30
                and b_count == 15
                and n_count == 15
                and degree_counts.get(
                    2,
                    0,
                )
                == 30
            )

            row = {
                "end": end,
                "removed_graph_shells": (
                    depth
                ),
                "removed_atoms": len(
                    removed
                ),
                "frontier_atoms": len(
                    frontier
                ),
                "frontier_B_atoms": (
                    b_count
                ),
                "frontier_N_atoms": (
                    n_count
                ),
                "frontier_degree0": (
                    degree_counts.get(
                        0,
                        0,
                    )
                ),
                "frontier_degree1": (
                    degree_counts.get(
                        1,
                        0,
                    )
                ),
                "frontier_degree2": (
                    degree_counts.get(
                        2,
                        0,
                    )
                ),
                "frontier_degree3": (
                    degree_counts.get(
                        3,
                        0,
                    )
                ),
                "simple_30_site_mixed_degree2_rim": (
                    simple_mixed_degree2_rim
                ),
            }

            reconstruction_rows.append(
                row
            )

            if simple_mixed_degree2_rim:
                reconstruction_candidates.append(
                    row
                )

    write_csv(
        RECONSTRUCTION_SCREEN_CSV,
        reconstruction_rows,
    )

    (
        outer_vertices,
        outer_edges,
    ) = build_hexagonal_flake(
        outer_shell
    )

    (
        removed_vertices,
        _,
    ) = build_hexagonal_flake(
        inner_shell
    )

    annulus_vertices = (
        outer_vertices
        - removed_vertices
    )

    annulus_edges = {
        edge
        for edge in outer_edges
        if (
            edge[0] in annulus_vertices
            and edge[1] in annulus_vertices
        )
    }

    outer_full_adjacency = build_adjacency(
        outer_vertices,
        outer_edges,
    )

    annulus_adjacency = build_adjacency(
        annulus_vertices,
        annulus_edges,
    )

    annulus_bipartite, annulus_colors = (
        bipartite_coloring(
            annulus_adjacency
        )
    )

    if not annulus_bipartite:
        raise RuntimeError(
            "Selected annulus is not bipartite."
        )

    outer_boundary = sorted(
        (
            vertex
            for vertex in annulus_vertices
            if (
                len(
                    outer_full_adjacency[
                        vertex
                    ]
                )
                == 2
                and len(
                    annulus_adjacency[
                        vertex
                    ]
                )
                == 2
            )
        ),
        key=lambda vertex: angle_turns(
            vertex,
            bond_length_nm,
        ),
    )

    inner_boundary = sorted(
        (
            vertex
            for vertex in annulus_vertices
            if (
                len(
                    outer_full_adjacency[
                        vertex
                    ]
                )
                == 3
                and len(
                    annulus_adjacency[
                        vertex
                    ]
                )
                == 2
            )
        ),
        key=lambda vertex: angle_turns(
            vertex,
            bond_length_nm,
        ),
    )

    annulus_degree_counts = Counter(
        len(neighbors)
        for neighbors
        in annulus_adjacency.values()
    )

    if len(annulus_vertices) != EXPECTED_ANNULUS_HEAVY_ATOMS:
        raise RuntimeError(
            "Unexpected selected-annulus atom count: "
            f"{len(annulus_vertices)}/"
            f"{EXPECTED_ANNULUS_HEAVY_ATOMS}"
        )

    if len(outer_boundary) != 30:
        raise RuntimeError(
            "Unexpected outer-boundary population."
        )

    if set(outer_boundary) & set(
        inner_boundary
    ):
        raise RuntimeError(
            "Outer and inner boundaries overlap."
        )

    seed_nodes_by_end: dict[
        str,
        list[dict[str, str]]
    ] = {}

    seed_node_by_global: dict[
        int,
        dict[str, str]
    ] = {}

    for row in seed_nodes:
        global_index = int(
            float(
                row[
                    "added_node_global_index_0based"
                ]
            )
        )

        seed_node_by_global[
            global_index
        ] = row

    for end in (
        "LOWER",
        "UPPER",
    ):
        rows = [
            row
            for row in seed_nodes
            if row.get(
                "end"
            )
            == end
        ]

        rows.sort(
            key=lambda row: int(
                float(
                    row[
                        "circumferential_index"
                    ]
                )
            )
        )

        if len(rows) != EXPECTED_SEED_ATOMS_PER_END:
            raise RuntimeError(
                f"{end}: unexpected seed-node count."
            )

        seed_nodes_by_end[
            end
        ] = rows

    parent_seed_adjacency: dict[
        str,
        set[str]
    ] = {}

    parent_seed_edge_set: set[
        tuple[str, str]
    ] = set()

    for index in range(
        EXPECTED_PARENT_ATOMS
    ):
        node_id = (
            f"P:{index}"
        )

        parent_seed_adjacency[
            node_id
        ] = set()

    for first, second in parent_edge_set:
        add_graph_edge(
            parent_seed_adjacency,
            parent_seed_edge_set,
            f"P:{first}",
            f"P:{second}",
        )

    for row in seed_nodes:
        global_index = int(
            float(
                row[
                    "added_node_global_index_0based"
                ]
            )
        )

        parent_seed_adjacency[
            f"S:{global_index}"
        ] = set()

    for row in seed_edges:
        parent_index = int(
            float(
                row[
                    "parent_node_0based"
                ]
            )
        )

        seed_index = int(
            float(
                row[
                    "added_node_0based"
                ]
            )
        )

        add_graph_edge(
            parent_seed_adjacency,
            parent_seed_edge_set,
            f"P:{parent_index}",
            f"S:{seed_index}",
        )

    mapping_rows: list[
        dict[str, Any]
    ] = []

    selected_mapping_by_end: dict[
        str,
        dict[str, Any]
    ] = {}

    for end in (
        "LOWER",
        "UPPER",
    ):
        end_seed_rows = (
            seed_nodes_by_end[
                end
            ]
        )

        seed_element = (
            end_seed_rows[0][
                "element"
            ]
        )

        complementary_element = (
            "B"
            if seed_element == "N"
            else "N"
        )

        color_element = {
            0: complementary_element,
            1: seed_element,
        }

        attachable_outer = sorted(
            (
                vertex
                for vertex in outer_boundary
                if color_element[
                    annulus_colors[
                        vertex
                    ]
                ]
                == complementary_element
            ),
            key=lambda vertex: angle_turns(
                vertex,
                bond_length_nm,
            ),
        )

        if len(
            attachable_outer
        ) != EXPECTED_DIRECT_ATTACHMENTS_PER_END:
            raise RuntimeError(
                f"{end}: expected 15 complementary "
                "outer-boundary sites."
            )

        annulus_distance_maps = {
            vertex: shortest_path_lengths(
                annulus_adjacency,
                vertex,
            )
            for vertex in attachable_outer
        }

        candidate_rows = []

        for seed_parity in (
            0,
            1,
        ):
            selected_seed_rows = [
                row
                for row in end_seed_rows
                if int(
                    float(
                        row[
                            "circumferential_index"
                        ]
                    )
                )
                % 2
                == seed_parity
            ]

            if len(
                selected_seed_rows
            ) != EXPECTED_DIRECT_ATTACHMENTS_PER_END:
                raise RuntimeError(
                    f"{end}: alternating seed selection "
                    "did not produce 15 sites."
                )

            selected_seed_ids = [
                "S:"
                + str(
                    int(
                        float(
                            row[
                                "added_node_global_index_0based"
                            ]
                        )
                    )
                )
                for row in selected_seed_rows
            ]

            seed_distance_maps = {
                node_id: shortest_path_lengths(
                    parent_seed_adjacency,
                    node_id,
                )
                for node_id in selected_seed_ids
            }

            seed_angles = [
                (
                    int(
                        float(
                            row[
                                "circumferential_index"
                            ]
                        )
                    )
                    / EXPECTED_SEED_ATOMS_PER_END
                )
                % 1.0
                for row in selected_seed_rows
            ]

            for orientation in (
                1,
                -1,
            ):
                for rotation in range(
                    EXPECTED_DIRECT_ATTACHMENTS_PER_END
                ):
                    mapped_outer = [
                        attachable_outer[
                            (
                                orientation
                                * index
                                + rotation
                            )
                            % EXPECTED_DIRECT_ATTACHMENTS_PER_END
                        ]
                        for index in range(
                            EXPECTED_DIRECT_ATTACHMENTS_PER_END
                        )
                    ]

                    outer_angles = [
                        angle_turns(
                            vertex,
                            bond_length_nm,
                        )
                        for vertex in mapped_outer
                    ]

                    (
                        phase_turns,
                        max_angular_residual,
                    ) = optimal_phase_residual(
                        seed_angles,
                        outer_angles,
                    )

                    created_cycle_lengths = []

                    for first_index in range(
                        EXPECTED_DIRECT_ATTACHMENTS_PER_END
                    ):
                        for second_index in range(
                            first_index + 1,
                            EXPECTED_DIRECT_ATTACHMENTS_PER_END,
                        ):
                            first_seed = (
                                selected_seed_ids[
                                    first_index
                                ]
                            )

                            second_seed = (
                                selected_seed_ids[
                                    second_index
                                ]
                            )

                            seed_distance = (
                                seed_distance_maps[
                                    first_seed
                                ][
                                    second_seed
                                ]
                            )

                            first_outer = (
                                mapped_outer[
                                    first_index
                                ]
                            )

                            second_outer = (
                                mapped_outer[
                                    second_index
                                ]
                            )

                            annulus_distance = (
                                annulus_distance_maps[
                                    first_outer
                                ][
                                    second_outer
                                ]
                            )

                            created_cycle_lengths.append(
                                seed_distance
                                + annulus_distance
                                + 2
                            )

                    minimum_created_cycle = min(
                        created_cycle_lengths
                    )

                    mean_created_cycle = sum(
                        created_cycle_lengths
                    ) / len(
                        created_cycle_lengths
                    )

                    candidate = {
                        "end": end,
                        "seed_element": (
                            seed_element
                        ),
                        "annulus_attachment_element": (
                            complementary_element
                        ),
                        "seed_parity": (
                            seed_parity
                        ),
                        "orientation": (
                            orientation
                        ),
                        "rotation": (
                            rotation
                        ),
                        "global_phase_turns": (
                            phase_turns
                        ),
                        "maximum_angular_residual_turns": (
                            max_angular_residual
                        ),
                        "minimum_attachment_created_cycle_length": (
                            minimum_created_cycle
                        ),
                        "mean_attachment_created_cycle_length": (
                            mean_created_cycle
                        ),
                        "angular_locality_pass": (
                            max_angular_residual
                            <= MAX_MAPPING_ANGULAR_RESIDUAL_TURNS
                        ),
                        "cycle_length_pass": (
                            minimum_created_cycle
                            >= MIN_ATTACHMENT_CREATED_CYCLE_LENGTH
                        ),
                    }

                    candidate_rows.append(
                        candidate
                    )

                    mapping_rows.append(
                        candidate
                    )

        feasible_candidates = [
            row
            for row in candidate_rows
            if bool(
                row[
                    "angular_locality_pass"
                ]
            )
            and bool(
                row[
                    "cycle_length_pass"
                ]
            )
        ]

        if not feasible_candidates:
            raise RuntimeError(
                f"{end}: no feasible partial-attachment "
                "mapping was found."
            )

        feasible_candidates.sort(
            key=lambda row: (
                float(
                    row[
                        "maximum_angular_residual_turns"
                    ]
                ),
                -int(
                    row[
                        "minimum_attachment_created_cycle_length"
                    ]
                ),
                float(
                    row[
                        "mean_attachment_created_cycle_length"
                    ]
                ),
                int(
                    row[
                        "seed_parity"
                    ]
                ),
                -int(
                    row[
                        "orientation"
                    ]
                ),
                int(
                    row[
                        "rotation"
                    ]
                ),
            )
        )

        selected_mapping_by_end[
            end
        ] = dict(
            feasible_candidates[0]
        )

    for row in mapping_rows:
        selected = (
            selected_mapping_by_end[
                row["end"]
            ]
        )

        row["selected"] = all(
            row[key]
            == selected[key]
            for key in (
                "seed_parity",
                "orientation",
                "rotation",
            )
        )

    write_csv(
        MAPPING_SEARCH_CSV,
        mapping_rows,
    )

    graph_adjacency: dict[
        str,
        set[str]
    ] = {}

    graph_edges: set[
        tuple[str, str]
    ] = set()

    node_elements: dict[
        str,
        str
    ] = {}

    node_types: dict[
        str,
        str
    ] = {}

    node_rows: list[
        dict[str, Any]
    ] = []

    edge_rows: list[
        dict[str, Any]
    ] = []

    def register_node(
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

        node_elements[node_id] = (
            element
        )

        node_types[node_id] = (
            node_type
        )

        graph_adjacency[node_id] = set()

        node_rows.append(
            {
                "node_id": node_id,
                "element": element,
                "node_type": (
                    node_type
                ),
                "end": end,
                **metadata,
                "coordinates_assigned": False,
                "formal_charge_assigned": False,
                "force_field_type_assigned": False,
            }
        )

    def register_edge(
        first: str,
        second: str,
        edge_type: str,
        end: str,
        heavy_atom_edge: bool,
    ) -> None:
        add_graph_edge(
            graph_adjacency,
            graph_edges,
            first,
            second,
        )

        first_element = node_elements[
            first
        ]

        second_element = node_elements[
            second
        ]

        edge_rows.append(
            {
                "edge_id": (
                    f"E:{len(edge_rows) + 1}"
                ),
                "source_node": first,
                "target_node": second,
                "source_element": (
                    first_element
                ),
                "target_element": (
                    second_element
                ),
                "edge_type": (
                    edge_type
                ),
                "end": end,
                "heavy_atom_edge": (
                    heavy_atom_edge
                ),
                "heteropolar_BN_edge": (
                    heavy_atom_edge
                    and {
                        first_element,
                        second_element,
                    }
                    == {
                        "B",
                        "N",
                    }
                ),
                "coordinates_assigned": False,
                "formal_bond_order_assigned": False,
            }
        )

    for index in range(
        EXPECTED_PARENT_ATOMS
    ):
        register_node(
            f"P:{index}",
            parent_elements[
                index
            ],
            "PARENT_HBN",
            "PARENT",
            {
                "source_index_0based": (
                    index
                ),
            },
        )

    for first, second in sorted(
        parent_edge_set
    ):
        register_edge(
            f"P:{first}",
            f"P:{second}",
            "PARENT_BN",
            "PARENT",
            True,
        )

    for row in seed_nodes:
        global_index = int(
            float(
                row[
                    "added_node_global_index_0based"
                ]
            )
        )

        register_node(
            f"S:{global_index}",
            row[
                "element"
            ],
            "HEXAGONAL_EDGE_COMPLETION_SEED",
            row[
                "end"
            ],
            {
                "source_index_0based": (
                    global_index
                ),
                "circumferential_index": (
                    int(
                        float(
                            row[
                                "circumferential_index"
                            ]
                        )
                    )
                ),
            },
        )

    for row in seed_edges:
        parent_index = int(
            float(
                row[
                    "parent_node_0based"
                ]
            )
        )

        seed_index = int(
            float(
                row[
                    "added_node_0based"
                ]
            )
        )

        register_edge(
            f"P:{parent_index}",
            f"S:{seed_index}",
            "PARENT_TO_COMPLETION_SEED",
            row[
                "end"
            ],
            True,
        )

    end_summary_rows = []

    for end in (
        "LOWER",
        "UPPER",
    ):
        end_seed_rows = (
            seed_nodes_by_end[
                end
            ]
        )

        seed_element = (
            end_seed_rows[0][
                "element"
            ]
        )

        complementary_element = (
            "B"
            if seed_element == "N"
            else "N"
        )

        color_element = {
            0: complementary_element,
            1: seed_element,
        }

        annulus_node_id: dict[
            tuple[int, int],
            str
        ] = {}

        for vertex in sorted(
            annulus_vertices
        ):
            node_id = (
                f"A:{end}:"
                f"{vertex[0]}:{vertex[1]}"
            )

            annulus_node_id[
                vertex
            ] = node_id

            node_type = (
                "ANNULUS_INNER_BOUNDARY"
                if vertex in inner_boundary
                else (
                    "ANNULUS_OUTER_BOUNDARY"
                    if vertex in outer_boundary
                    else "ANNULUS_INTERIOR"
                )
            )

            register_node(
                node_id,
                color_element[
                    annulus_colors[
                        vertex
                    ]
                ],
                node_type,
                end,
                {
                    "lattice_x": (
                        vertex[0]
                    ),
                    "lattice_y": (
                        vertex[1]
                    ),
                    "radius_proxy_nm": (
                        radius_nm(
                            vertex,
                            bond_length_nm,
                        )
                    ),
                    "angle_turns": (
                        angle_turns(
                            vertex,
                            bond_length_nm,
                        )
                    ),
                },
            )

        for first, second in sorted(
            annulus_edges
        ):
            register_edge(
                annulus_node_id[
                    first
                ],
                annulus_node_id[
                    second
                ],
                "ANNULUS_BN",
                end,
                True,
            )

        selected_mapping = (
            selected_mapping_by_end[
                end
            ]
        )

        seed_parity = int(
            selected_mapping[
                "seed_parity"
            ]
        )

        orientation = int(
            selected_mapping[
                "orientation"
            ]
        )

        rotation = int(
            selected_mapping[
                "rotation"
            ]
        )

        selected_seed_rows = [
            row
            for row in end_seed_rows
            if int(
                float(
                    row[
                        "circumferential_index"
                    ]
                )
            )
            % 2
            == seed_parity
        ]

        unselected_seed_rows = [
            row
            for row in end_seed_rows
            if row not in selected_seed_rows
        ]

        attachable_outer = sorted(
            (
                vertex
                for vertex in outer_boundary
                if color_element[
                    annulus_colors[
                        vertex
                    ]
                ]
                == complementary_element
            ),
            key=lambda vertex: angle_turns(
                vertex,
                bond_length_nm,
            ),
        )

        mapped_outer = [
            attachable_outer[
                (
                    orientation
                    * index
                    + rotation
                )
                % EXPECTED_DIRECT_ATTACHMENTS_PER_END
            ]
            for index in range(
                EXPECTED_DIRECT_ATTACHMENTS_PER_END
            )
        ]

        attached_outer = set(
            mapped_outer
        )

        for seed_row, outer_vertex in zip(
            selected_seed_rows,
            mapped_outer,
        ):
            seed_index = int(
                float(
                    seed_row[
                        "added_node_global_index_0based"
                    ]
                )
            )

            register_edge(
                f"S:{seed_index}",
                annulus_node_id[
                    outer_vertex
                ],
                "PARTIAL_HETEROPOLAR_SEED_TO_ANNULUS",
                end,
                True,
            )

        passivant_count = 0

        for seed_row in unselected_seed_rows:
            seed_index = int(
                float(
                    seed_row[
                        "added_node_global_index_0based"
                    ]
                )
            )

            hydrogen_id = (
                f"H:{end}:SEED:"
                f"{passivant_count}"
            )

            passivant_count += 1

            register_node(
                hydrogen_id,
                "H",
                "SEED_PASSIVANT_H",
                end,
                {
                    "attached_to": (
                        f"S:{seed_index}"
                    ),
                },
            )

            register_edge(
                f"S:{seed_index}",
                hydrogen_id,
                "SEED_H_PASSIVATION",
                end,
                False,
            )

        unconnected_outer = [
            vertex
            for vertex in outer_boundary
            if vertex not in attached_outer
        ]

        for outer_vertex in unconnected_outer:
            hydrogen_id = (
                f"H:{end}:OUTER:"
                f"{passivant_count}"
            )

            passivant_count += 1

            register_node(
                hydrogen_id,
                "H",
                "ANNULUS_OUTER_PASSIVANT_H",
                end,
                {
                    "attached_to": (
                        annulus_node_id[
                            outer_vertex
                        ]
                    ),
                },
            )

            register_edge(
                annulus_node_id[
                    outer_vertex
                ],
                hydrogen_id,
                "ANNULUS_OUTER_H_PASSIVATION",
                end,
                False,
            )

        for inner_vertex in inner_boundary:
            hydrogen_id = (
                f"H:{end}:INNER:"
                f"{passivant_count}"
            )

            passivant_count += 1

            register_node(
                hydrogen_id,
                "H",
                "ANNULUS_INNER_PASSIVANT_H",
                end,
                {
                    "attached_to": (
                        annulus_node_id[
                            inner_vertex
                        ]
                    ),
                },
            )

            register_edge(
                annulus_node_id[
                    inner_vertex
                ],
                hydrogen_id,
                "ANNULUS_INNER_H_PASSIVATION",
                end,
                False,
            )

        outer_radii = [
            radius_nm(
                vertex,
                bond_length_nm,
            )
            for vertex in outer_boundary
        ]

        inner_radii = [
            radius_nm(
                vertex,
                bond_length_nm,
            )
            for vertex in inner_boundary
        ]

        outer_radius_mean_nm = (
            sum(
                outer_radii
            )
            / len(
                outer_radii
            )
        )

        inner_radius_mean_nm = (
            sum(
                inner_radii
            )
            / len(
                inner_radii
            )
        )

        end_summary_rows.append(
            {
                "end": end,
                "seed_element": (
                    seed_element
                ),
                "annulus_attachment_element": (
                    complementary_element
                ),
                "seed_sites": len(
                    end_seed_rows
                ),
                "covalently_attached_seed_sites": (
                    len(
                        selected_seed_rows
                    )
                ),
                "H_passivated_seed_sites": (
                    len(
                        unselected_seed_rows
                    )
                ),
                "annulus_heavy_atoms": (
                    len(
                        annulus_vertices
                    )
                ),
                "annulus_outer_boundary_atoms": (
                    len(
                        outer_boundary
                    )
                ),
                "annulus_inner_boundary_atoms": (
                    len(
                        inner_boundary
                    )
                ),
                "covalently_attached_outer_sites": (
                    len(
                        attached_outer
                    )
                ),
                "H_passivated_outer_sites": (
                    len(
                        unconnected_outer
                    )
                ),
                "H_passivated_inner_sites": (
                    len(
                        inner_boundary
                    )
                ),
                "total_H_passivants": (
                    passivant_count
                ),
                "total_added_heavy_atoms": (
                    EXPECTED_SEED_ATOMS_PER_END
                    + len(
                        annulus_vertices
                    )
                ),
                "mapping_seed_parity": (
                    seed_parity
                ),
                "mapping_orientation": (
                    orientation
                ),
                "mapping_rotation": (
                    rotation
                ),
                "mapping_maximum_angular_residual_turns": (
                    selected_mapping[
                        "maximum_angular_residual_turns"
                    ]
                ),
                "minimum_attachment_created_cycle_length": (
                    selected_mapping[
                        "minimum_attachment_created_cycle_length"
                    ]
                ),
                "outer_radius_mean_nm": (
                    outer_radius_mean_nm
                ),
                "inner_radius_mean_nm": (
                    inner_radius_mean_nm
                ),
            }
        )

    write_csv(
        SELECTED_NODES_CSV,
        node_rows,
    )

    write_csv(
        SELECTED_EDGES_CSV,
        edge_rows,
    )

    write_csv(
        END_SUMMARY_CSV,
        end_summary_rows,
    )

    total_degrees = {
        node: len(neighbors)
        for node, neighbors
        in graph_adjacency.items()
    }

    heavy_nodes = {
        node
        for node, element
        in node_elements.items()
        if element != "H"
    }

    hydrogen_nodes = {
        node
        for node, element
        in node_elements.items()
        if element == "H"
    }

    bn_degree_failures = [
        node
        for node in heavy_nodes
        if total_degrees[
            node
        ]
        != 3
    ]

    hydrogen_degree_failures = [
        node
        for node in hydrogen_nodes
        if total_degrees[
            node
        ]
        != 1
    ]

    heavy_edge_rows = [
        row
        for row in edge_rows
        if bool(
            row[
                "heavy_atom_edge"
            ]
        )
    ]

    nonheteropolar_heavy_edges = [
        row
        for row in heavy_edge_rows
        if not bool(
            row[
                "heteropolar_BN_edge"
            ]
        )
    ]

    heavy_adjacency = {
        node: {
            neighbor
            for neighbor in graph_adjacency[
                node
            ]
            if neighbor in heavy_nodes
        }
        for node in heavy_nodes
    }

    heavy_components = connected_components(
        heavy_adjacency
    )

    heavy_bipartite, _ = (
        bipartite_coloring(
            heavy_adjacency
        )
    )

    heavy_four_cycles = count_four_cycles(
        heavy_adjacency
    )

    heavy_girth = graph_girth(
        heavy_adjacency
    )

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

    combined_heavy_added = (
        int(
            lower_summary[
                "total_added_heavy_atoms"
            ]
        )
        + int(
            upper_summary[
                "total_added_heavy_atoms"
            ]
        )
    )

    combined_hydrogens = (
        int(
            lower_summary[
                "total_H_passivants"
            ]
        )
        + int(
            upper_summary[
                "total_H_passivants"
            ]
        )
    )

    heavy_atom_relative_error = (
        abs(
            EXPECTED_TOTAL_HEAVY_ATOMS_PER_END
            - target_heavy_atoms_per_end
        )
        / target_heavy_atoms_per_end
    )

    outer_radius_relative_error = (
        abs(
            float(
                lower_summary[
                    "outer_radius_mean_nm"
                ]
            )
            - parent_rim_radius_nm
        )
        / parent_rim_radius_nm
    )

    inner_radius_relative_error = (
        abs(
            float(
                lower_summary[
                    "inner_radius_mean_nm"
                ]
            )
            - target_aperture_radius_nm
        )
        / target_aperture_radius_nm
    )

    reconstruction_feasible = (
        len(
            reconstruction_candidates
        )
        > 0
    )

    gates = {
        "Gate3A_parent_audit_is_accepted": (
            parent_summary.get(
                "decision"
            )
            == EXPECTED_PARENT_DECISION
        ),
        "Gate3D_hexagonal_edge_completion_seed_is_accepted": (
            seed_summary.get(
                "decision"
            )
            == EXPECTED_SEED_DECISION
        ),
        "Gate3E_direct_30_site_attachment_was_rejected": (
            annulus_summary.get(
                "decision"
            )
            == EXPECTED_ANNULUS_DECISION
        ),
        "simple_full_shell_reconstruction_does_not_create_30_site_15B_15N_rim": (
            not reconstruction_feasible
        ),
        "selected_annulus_is_n5_m2": (
            outer_shell == 5
            and inner_shell == 2
        ),
        "selected_annulus_has_126_heavy_atoms": (
            len(
                annulus_vertices
            )
            == EXPECTED_ANNULUS_HEAVY_ATOMS
        ),
        "selected_annulus_has_30_outer_boundary_sites": (
            len(
                outer_boundary
            )
            == 30
        ),
        "selected_annulus_is_connected": (
            len(
                connected_components(
                    annulus_adjacency
                )
            )
            == 1
        ),
        "selected_annulus_is_bipartite": (
            annulus_bipartite
        ),
        "selected_annulus_has_no_four_member_cycles": (
            count_four_cycles(
                annulus_adjacency
            )
            == 0
        ),
        "selected_annulus_has_only_degree2_or_degree3_atoms": (
            all(
                degree in {
                    2,
                    3,
                }
                for degree
                in annulus_degree_counts
            )
        ),
        "15_alternating_seed_sites_are_attached_per_end": (
            all(
                int(
                    row[
                        "covalently_attached_seed_sites"
                    ]
                )
                == EXPECTED_DIRECT_ATTACHMENTS_PER_END
                for row in end_summary_rows
            )
        ),
        "15_remaining_seed_sites_are_H_passivated_per_end": (
            all(
                int(
                    row[
                        "H_passivated_seed_sites"
                    ]
                )
                == EXPECTED_UNATTACHED_SEED_SITES_PER_END
                for row in end_summary_rows
            )
        ),
        "15_noncomplementary_outer_sites_are_H_passivated_per_end": (
            all(
                int(
                    row[
                        "H_passivated_outer_sites"
                    ]
                )
                == EXPECTED_UNATTACHED_OUTER_SITES_PER_END
                for row in end_summary_rows
            )
        ),
        "all_inner_boundary_sites_are_H_passivated": (
            all(
                int(
                    row[
                        "H_passivated_inner_sites"
                    ]
                )
                == len(
                    inner_boundary
                )
                for row in end_summary_rows
            )
        ),
        "all_BN_atoms_reach_total_coordination3": (
            len(
                bn_degree_failures
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
        "combined_heavy_graph_is_connected": (
            len(
                heavy_components
            )
            == 1
        ),
        "combined_heavy_graph_is_bipartite": (
            heavy_bipartite
        ),
        "combined_heavy_graph_has_no_four_member_cycles": (
            heavy_four_cycles == 0
        ),
        "combined_heavy_graph_girth_is_at_least6": (
            heavy_girth >= 6
        ),
        "attachment_mapping_is_angularly_local": (
            all(
                float(
                    row[
                        "mapping_maximum_angular_residual_turns"
                    ]
                )
                <= MAX_MAPPING_ANGULAR_RESIDUAL_TURNS
                for row in end_summary_rows
            )
        ),
        "attachment_created_cycles_have_length_at_least6": (
            all(
                int(
                    row[
                        "minimum_attachment_created_cycle_length"
                    ]
                )
                >= MIN_ATTACHMENT_CREATED_CYCLE_LENGTH
                for row in end_summary_rows
            )
        ),
        "heavy_atom_population_error_is_within15_percent": (
            heavy_atom_relative_error
            <= MAX_HEAVY_ATOM_RELATIVE_ERROR
        ),
        "outer_radius_error_is_within15_percent": (
            outer_radius_relative_error
            <= MAX_OUTER_RADIUS_RELATIVE_ERROR
        ),
        "inner_radius_proxy_error_is_within30_percent": (
            inner_radius_relative_error
            <= MAX_INNER_RADIUS_RELATIVE_ERROR
        ),
        "no_coordinates_were_assigned": all(
            not bool(
                row[
                    "coordinates_assigned"
                ]
            )
            for row in node_rows
        )
        and all(
            not bool(
                row[
                    "coordinates_assigned"
                ]
            )
            for row in edge_rows
        ),
        "no_formal_charges_were_assigned": all(
            not bool(
                row[
                    "formal_charge_assigned"
                ]
            )
            for row in node_rows
        ),
        "no_force_field_types_were_assigned": all(
            not bool(
                row[
                    "force_field_type_assigned"
                ]
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
        "R2_RECONSTRUCTION_VS_PARTIAL_ATTACHMENT_"
        "CONTINGENCY_REQUIRES_REVIEW"
    )

    required_next_step = (
        "BUILD_AND_VALIDATE_R2_PARTIAL_ATTACHMENT_ANNULUS_"
        "STATIC_COORDINATE_EMBEDDING"
        if accepted
        else
        "REVIEW_R2_RECONSTRUCTION_VS_PARTIAL_ATTACHMENT_FAILURES"
    )

    summary = {
        "decision": decision,
        "simple_full_shell_reconstruction_candidate_count": (
            len(
                reconstruction_candidates
            )
        ),
        "selected_architecture": (
            "N5_M2_HEXAGONAL_BN_ANNULUS_WITH_15_ALTERNATING_"
            "HETEROPOLAR_ATTACHMENTS_AND_COMPLEMENTARY_H_PASSIVATION"
        ),
        "annulus_outer_shell": (
            outer_shell
        ),
        "annulus_inner_shell": (
            inner_shell
        ),
        "annulus_heavy_atoms_per_end": (
            len(
                annulus_vertices
            )
        ),
        "edge_completion_seed_atoms_per_end": (
            EXPECTED_SEED_ATOMS_PER_END
        ),
        "total_added_heavy_atoms_per_end": (
            EXPECTED_TOTAL_HEAVY_ATOMS_PER_END
        ),
        "target_heavy_atoms_per_end": (
            target_heavy_atoms_per_end
        ),
        "heavy_atom_relative_error": (
            heavy_atom_relative_error
        ),
        "covalent_seed_annulus_attachments_per_end": (
            EXPECTED_DIRECT_ATTACHMENTS_PER_END
        ),
        "seed_H_passivants_per_end": (
            EXPECTED_UNATTACHED_SEED_SITES_PER_END
        ),
        "outer_annulus_H_passivants_per_end": (
            EXPECTED_UNATTACHED_OUTER_SITES_PER_END
        ),
        "inner_annulus_H_passivants_per_end": (
            len(
                inner_boundary
            )
        ),
        "total_H_passivants_per_end": (
            int(
                lower_summary[
                    "total_H_passivants"
                ]
            )
        ),
        "combined_added_heavy_atoms": (
            combined_heavy_added
        ),
        "combined_added_H_atoms": (
            combined_hydrogens
        ),
        "heavy_graph_components": (
            len(
                heavy_components
            )
        ),
        "heavy_graph_bipartite": (
            heavy_bipartite
        ),
        "heavy_graph_girth": (
            heavy_girth
        ),
        "heavy_graph_four_member_cycles": (
            heavy_four_cycles
        ),
        "BN_coordination_failures": (
            len(
                bn_degree_failures
            )
        ),
        "H_coordination_failures": (
            len(
                hydrogen_degree_failures
            )
        ),
        "nonheteropolar_heavy_edges": (
            len(
                nonheteropolar_heavy_edges
            )
        ),
        "outer_radius_mean_nm": (
            lower_summary[
                "outer_radius_mean_nm"
            ]
        ),
        "inner_radius_mean_nm": (
            lower_summary[
                "inner_radius_mean_nm"
            ]
        ),
        "target_parent_rim_radius_nm": (
            parent_rim_radius_nm
        ),
        "target_aperture_radius_nm": (
            target_aperture_radius_nm
        ),
        "target_aperture_diameter_nm": (
            target_aperture_diameter_nm
        ),
        "outer_radius_relative_error": (
            outer_radius_relative_error
        ),
        "inner_radius_relative_error": (
            inner_radius_relative_error
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

    DESIGN_JSON.write_text(
        json.dumps(
            {
                "summary": summary,
                "end_summaries": (
                    end_summary_rows
                ),
                "gates": gates,
                "selected_mappings": (
                    selected_mapping_by_end
                ),
                "limitations": [
                    (
                        "The annulus uses an abstract two-dimensional "
                        "hexagonal-lattice template. It has not been "
                        "embedded relative to the BNNT in three dimensions."
                    ),
                    (
                        "Hydrogen orientations, B-H and N-H bond lengths, "
                        "bond angles, partial charges, and force-field "
                        "parameters have not been assigned."
                    ),
                    (
                        "Graph-level valence completion does not establish "
                        "energetic stability or synthetic realizability."
                    ),
                    (
                        "The effective hydrated aperture must be measured "
                        "after static three-dimensional embedding."
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
                "Gate3D_seed_nodes"
            ),
            "file": relative(
                SEED_NODES_CSV
            ),
            "sha256": sha256(
                SEED_NODES_CSV
            ),
        },
        {
            "role": (
                "Gate3D_seed_edges"
            ),
            "file": relative(
                SEED_EDGES_CSV
            ),
            "sha256": sha256(
                SEED_EDGES_CSV
            ),
        },
        {
            "role": (
                "Gate3E_annulus_feasibility_summary"
            ),
            "file": relative(
                ANNULUS_FEASIBILITY_SUMMARY_CSV
            ),
            "sha256": sha256(
                ANNULUS_FEASIBILITY_SUMMARY_CSV
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

    reconstruction_lines = "\n".join(
        (
            f"- {row['end']} depth {row['removed_graph_shells']}: "
            f"frontier={row['frontier_atoms']}; "
            f"B/N={row['frontier_B_atoms']}/"
            f"{row['frontier_N_atoms']}; "
            f"d1/d2/d3={row['frontier_degree1']}/"
            f"{row['frontier_degree2']}/"
            f"{row['frontier_degree3']}; "
            f"mixed-30-site-rim="
            f"{row['simple_30_site_mixed_degree2_rim']}"
        )
        for row in reconstruction_rows
    )

    REPORT_MD.write_text(
        f"""# R2 Reconstruction versus Partial-Attachment Contingency

## Scope

This gate compares two responses to the negative Gate 3E result:

1. removal of complete graph shells from the BNNT end in an attempt to
   expose a mixed 15 B / 15 N coordination-two rim;
2. partial heteropolar attachment of the validated hexagonal edge seed
   to the n=5, m=2 annulus, with explicit hydrogen passivation of every
   unmatched coordination-two site.

No three-dimensional molecular coordinates, partial charges,
force-field parameters, minimization, MD, or QM calculation were
generated.

## Full-shell parent reconstruction screen

{reconstruction_lines}

Simple full-shell reconstruction candidates:

- **{len(reconstruction_candidates)}**

## Selected partial-attachment architecture

- Annulus:
  **n={outer_shell}, m={inner_shell}**
- Annulus B/N atoms:
  **{len(annulus_vertices)} per end**
- Edge-completion seed:
  **{EXPECTED_SEED_ATOMS_PER_END} atoms per end**
- Total added heavy atoms:
  **{EXPECTED_TOTAL_HEAVY_ATOMS_PER_END} per end**
- Target heavy atoms:
  **{target_heavy_atoms_per_end:.3f} per end**
- Relative heavy-atom error:
  **{heavy_atom_relative_error:.6f}**
- Direct heteropolar attachments:
  **{EXPECTED_DIRECT_ATTACHMENTS_PER_END} per end**
- H-passivated seed sites:
  **{EXPECTED_UNATTACHED_SEED_SITES_PER_END} per end**
- H-passivated outer-annulus sites:
  **{EXPECTED_UNATTACHED_OUTER_SITES_PER_END} per end**
- H-passivated inner-annulus sites:
  **{len(inner_boundary)} per end**
- Total H passivants:
  **{lower_summary['total_H_passivants']} per end**

## Final graph audit

- Heavy-graph components:
  **{len(heavy_components)}**
- Heavy-graph bipartite:
  **{heavy_bipartite}**
- Heavy-graph girth:
  **{heavy_girth}**
- Four-member heavy cycles:
  **{heavy_four_cycles}**
- B/N coordination failures:
  **{len(bn_degree_failures)}**
- H coordination failures:
  **{len(hydrogen_degree_failures)}**
- Nonheteropolar heavy edges:
  **{len(nonheteropolar_heavy_edges)}**

## Geometric proxies

- Outer radius:
  **{float(lower_summary['outer_radius_mean_nm']):.6f} nm**
- Parent-rim target:
  **{parent_rim_radius_nm:.6f} nm**
- Outer-radius relative error:
  **{outer_radius_relative_error:.6f}**
- Inner radius:
  **{float(lower_summary['inner_radius_mean_nm']):.6f} nm**
- Aperture-radius target:
  **{target_aperture_radius_nm:.6f} nm**
- Inner-radius relative error:
  **{inner_radius_relative_error:.6f}**

The inner-radius value is a nucleus-position lattice proxy. The
effective steric aperture after inward H termination must be measured
during the static coordinate-embedding gate.

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
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 reconstruction versus partial-attachment "
        "contingency audit completed."
    )

    print(
        "Simple full-shell reconstruction candidates: "
        f"{len(reconstruction_candidates)}"
    )

    for end in (
        "LOWER",
        "UPPER",
    ):
        end_rows = [
            row
            for row in reconstruction_rows
            if row[
                "end"
            ]
            == end
        ]

        print(
            f"{end} reconstruction frontier "
            "depth:B/N/d2:"
        )

        for row in end_rows:
            print(
                "  "
                f"{row['removed_graph_shells']}:"
                f"{row['frontier_B_atoms']}/"
                f"{row['frontier_N_atoms']}/"
                f"{row['frontier_degree2']}"
            )

    print(
        "Selected annulus n/m / annulus atoms / "
        "total heavy atoms per end: "
        f"{outer_shell}/"
        f"{inner_shell}/"
        f"{len(annulus_vertices)}/"
        f"{EXPECTED_TOTAL_HEAVY_ATOMS_PER_END}"
    )

    print(
        "Target heavy atoms per end / relative error: "
        f"{target_heavy_atoms_per_end:.3f}/"
        f"{heavy_atom_relative_error:.6f}"
    )

    for row in end_summary_rows:
        print(
            f"{row['end']} attached seed / H-seed / "
            "H-outer / H-inner / total H: "
            f"{row['covalently_attached_seed_sites']}/"
            f"{row['H_passivated_seed_sites']}/"
            f"{row['H_passivated_outer_sites']}/"
            f"{row['H_passivated_inner_sites']}/"
            f"{row['total_H_passivants']}"
        )

        print(
            f"{row['end']} mapping parity/orientation/"
            "rotation/angular residual/min cycle: "
            f"{row['mapping_seed_parity']}/"
            f"{row['mapping_orientation']}/"
            f"{row['mapping_rotation']}/"
            f"{float(row['mapping_maximum_angular_residual_turns']):.6f}/"
            f"{row['minimum_attachment_created_cycle_length']}"
        )

    print(
        "Combined added heavy / H atoms: "
        f"{combined_heavy_added}/"
        f"{combined_hydrogens}"
    )

    print(
        "B/N coordination failures / H failures / "
        "nonheteropolar heavy edges: "
        f"{len(bn_degree_failures)}/"
        f"{len(hydrogen_degree_failures)}/"
        f"{len(nonheteropolar_heavy_edges)}"
    )

    print(
        "Heavy graph components / bipartite / "
        "girth / four-cycles: "
        f"{len(heavy_components)}/"
        f"{heavy_bipartite}/"
        f"{heavy_girth}/"
        f"{heavy_four_cycles}"
    )

    print(
        "Outer radius / target / relative error: "
        f"{float(lower_summary['outer_radius_mean_nm']):.6f}/"
        f"{parent_rim_radius_nm:.6f}/"
        f"{outer_radius_relative_error:.6f}"
    )

    print(
        "Inner radius / target / relative error: "
        f"{float(lower_summary['inner_radius_mean_nm']):.6f}/"
        f"{target_aperture_radius_nm:.6f}/"
        f"{inner_radius_relative_error:.6f}"
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
        RECONSTRUCTION_SCREEN_CSV,
        MAPPING_SEARCH_CSV,
        SELECTED_NODES_CSV,
        SELECTED_EDGES_CSV,
        END_SUMMARY_CSV,
        SUMMARY_CSV,
        GATES_CSV,
        DESIGN_JSON,
        SOURCE_MANIFEST_CSV,
        REPORT_MD,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not accepted:
        raise RuntimeError(
            "R2 reconstruction versus partial-attachment "
            "contingency requires review."
        )


if __name__ == "__main__":
    main()
