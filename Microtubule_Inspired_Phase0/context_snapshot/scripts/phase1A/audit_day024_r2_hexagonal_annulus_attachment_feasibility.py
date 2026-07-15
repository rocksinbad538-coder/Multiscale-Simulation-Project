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

GATE3C1_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "04_r2_graded_collar_cycle_topology_audit"
)

GATE3D_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "05_r2_hexagonal_edge_completion_seed"
)

OUTPUT_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "06_r2_hexagonal_annulus_attachment_feasibility"
)

PARENT_SUMMARY_CSV = (
    GATE3A_ROOT
    / "r2_parent_rim_chemical_audit_summary.csv"
)

CYCLE_AUDIT_SUMMARY_CSV = (
    GATE3C1_ROOT
    / "r2_collar_cycle_topology_audit_summary.csv"
)

SEED_SUMMARY_CSV = (
    GATE3D_ROOT
    / "r2_hexagonal_edge_completion_seed_summary.csv"
)

SEED_GATES_CSV = (
    GATE3D_ROOT
    / "r2_hexagonal_edge_completion_seed_gates.csv"
)

TEMPLATES_CSV = (
    OUTPUT_ROOT
    / "r2_hexagonal_annulus_template_screen.csv"
)

BEST_CANDIDATES_CSV = (
    OUTPUT_ROOT
    / "r2_hexagonal_annulus_best_candidates.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_hexagonal_annulus_attachment_feasibility_summary.csv"
)

GATES_CSV = (
    OUTPUT_ROOT
    / "r2_hexagonal_annulus_attachment_audit_gates.csv"
)

AUDIT_JSON = (
    OUTPUT_ROOT
    / "r2_hexagonal_annulus_attachment_feasibility.json"
)

SOURCE_MANIFEST_CSV = (
    OUTPUT_ROOT
    / "r2_hexagonal_annulus_attachment_source_manifest.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_HEXAGONAL_ANNULUS_ATTACHMENT_FEASIBILITY_DAY024.md"
)

EXPECTED_PARENT_DECISION = (
    "R2_PARENT_RIM_AND_CHEMICAL_CONSTRAINTS_AUDITED"
)

EXPECTED_CYCLE_DECISION = (
    "R2_GRADED_HETEROPOLAR_BN_COLLAR_REQUIRES_HEXAGONAL_GRAPH_REDESIGN"
)

EXPECTED_SEED_DECISION = (
    "R2_HEXAGONAL_EDGE_COMPLETION_SEED_VALIDATED"
)

PASS_DECISION = (
    "R2_PURE_HEXAGONAL_BN_ANNULUS_DIRECT_ATTACHMENT_FEASIBLE"
)

NEGATIVE_DECISION = (
    "R2_PURE_HEXAGONAL_BN_ANNULUS_DIRECT_ATTACHMENT_NOT_FEASIBLE_"
    "WITH_CURRENT_HOMOPOLAR_DEGREE2_RIM"
)

EXPECTED_SEED_ATOMS_PER_END = 30
EXPECTED_ATTACHMENT_SITES_PER_END = 30

MIN_OUTER_SHELL = 2
MAX_OUTER_SHELL = 10

MAX_TOTAL_HEAVY_RELATIVE_ERROR = 0.15
MAX_OUTER_RADIUS_RELATIVE_ERROR = 0.15
MAX_INNER_RADIUS_RELATIVE_ERROR = 0.30

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
        rows = list(
            csv.DictReader(handle)
        )

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


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


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
    vertices: set[tuple[int, int]],
    edges: set[
        tuple[
            tuple[int, int],
            tuple[int, int],
        ]
    ],
) -> dict[
    tuple[int, int],
    set[tuple[int, int]],
]:
    adjacency = {
        vertex: set()
        for vertex in vertices
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
    adjacency: dict[
        tuple[int, int],
        set[tuple[int, int]],
    ],
) -> list[
    set[tuple[int, int]]
]:
    remaining = set(
        adjacency
    )

    components = []

    while remaining:
        start = min(
            remaining
        )

        component = set()

        queue: deque[
            tuple[int, int]
        ] = deque(
            [start]
        )

        while queue:
            node = queue.popleft()

            if node in component:
                continue

            component.add(
                node
            )

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
    adjacency: dict[
        tuple[int, int],
        set[tuple[int, int]],
    ],
) -> tuple[
    bool,
    dict[
        tuple[int, int],
        int,
    ],
]:
    colors: dict[
        tuple[int, int],
        int,
    ] = {}

    for start in sorted(
        adjacency
    ):
        if start in colors:
            continue

        colors[start] = 0

        queue: deque[
            tuple[int, int]
        ] = deque(
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
    adjacency: dict[
        tuple[int, int],
        set[tuple[int, int]],
    ],
) -> int:
    nodes = sorted(
        adjacency
    )

    raw_count = 0

    for first_index, first in enumerate(
        nodes
    ):
        for second in nodes[
            first_index + 1:
        ]:
            common = (
                adjacency[first]
                & adjacency[second]
            )

            if len(common) >= 2:
                raw_count += (
                    len(common)
                    * (
                        len(common) - 1
                    )
                    // 2
                )

    if raw_count % 2 != 0:
        raise RuntimeError(
            "Four-cycle count was not divisible by two."
        )

    return raw_count // 2


def physical_radius_nm(
    vertex: tuple[int, int],
    bond_length_nm: float,
) -> float:
    integer_x, integer_y = vertex

    x_nm = (
        integer_x
        * bond_length_nm
        / 2.0
    )

    y_nm = (
        integer_y
        * math.sqrt(3.0)
        * bond_length_nm
        / 2.0
    )

    return math.hypot(
        x_nm,
        y_nm,
    )


def mean_value(
    values: list[float],
) -> float:
    if not values:
        raise RuntimeError(
            "Cannot calculate the mean of an empty list."
        )

    return sum(values) / len(values)


def build_annulus_template(
    outer_shell: int,
    inner_shell: int,
    bond_length_nm: float,
    parent_rim_radius_nm: float,
    target_aperture_radius_nm: float,
    target_total_added_heavy_atoms_per_end: float,
) -> dict[str, Any]:
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
            edge[0]
            in annulus_vertices
            and edge[1]
            in annulus_vertices
        )
    }

    outer_full_adjacency = (
        build_adjacency(
            outer_vertices,
            outer_edges,
        )
    )

    annulus_adjacency = (
        build_adjacency(
            annulus_vertices,
            annulus_edges,
        )
    )

    degrees = Counter(
        len(neighbors)
        for neighbors
        in annulus_adjacency.values()
    )

    outer_boundary = [
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
    ]

    inner_boundary = [
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
    ]

    bipartite, colors = (
        bipartite_coloring(
            annulus_adjacency
        )
    )

    outer_color_counts = Counter(
        colors[
            vertex
        ]
        for vertex in outer_boundary
    )

    inner_color_counts = Counter(
        colors[
            vertex
        ]
        for vertex in inner_boundary
    )

    max_heteropolar_attachments = max(
        outer_color_counts.get(
            0,
            0,
        ),
        outer_color_counts.get(
            1,
            0,
        ),
    )

    homopolar_bonds_required_for_30_site_attachment = max(
        0,
        EXPECTED_ATTACHMENT_SITES_PER_END
        - max_heteropolar_attachments,
    )

    unused_outer_boundary_sites = max(
        0,
        len(outer_boundary)
        - EXPECTED_ATTACHMENT_SITES_PER_END,
    )

    outer_radii = [
        physical_radius_nm(
            vertex,
            bond_length_nm,
        )
        for vertex in outer_boundary
    ]

    inner_radii = [
        physical_radius_nm(
            vertex,
            bond_length_nm,
        )
        for vertex in inner_boundary
    ]

    outer_radius_mean_nm = (
        mean_value(
            outer_radii
        )
    )

    outer_radius_maximum_nm = max(
        outer_radii
    )

    inner_radius_mean_nm = (
        mean_value(
            inner_radii
        )
    )

    outer_radius_relative_error = abs(
        outer_radius_mean_nm
        - parent_rim_radius_nm
    ) / parent_rim_radius_nm

    inner_radius_relative_error = abs(
        inner_radius_mean_nm
        - target_aperture_radius_nm
    ) / target_aperture_radius_nm

    total_added_heavy_atoms_per_end = (
        EXPECTED_SEED_ATOMS_PER_END
        + len(
            annulus_vertices
        )
    )

    total_heavy_relative_error = abs(
        total_added_heavy_atoms_per_end
        - target_total_added_heavy_atoms_per_end
    ) / target_total_added_heavy_atoms_per_end

    direct_30_site_heteropolar_attachment_possible = (
        max_heteropolar_attachments
        >= EXPECTED_ATTACHMENT_SITES_PER_END
    )

    geometric_and_population_constraints_pass = (
        total_heavy_relative_error
        <= MAX_TOTAL_HEAVY_RELATIVE_ERROR
        and outer_radius_relative_error
        <= MAX_OUTER_RADIUS_RELATIVE_ERROR
        and inner_radius_relative_error
        <= MAX_INNER_RADIUS_RELATIVE_ERROR
    )

    graph_integrity_pass = (
        len(
            connected_components(
                annulus_adjacency
            )
        )
        == 1
        and bipartite
        and count_four_cycles(
            annulus_adjacency
        )
        == 0
        and degrees.get(
            1,
            0,
        )
        == 0
        and all(
            degree in {
                2,
                3,
            }
            for degree in degrees
        )
    )

    direct_candidate_pass = (
        direct_30_site_heteropolar_attachment_possible
        and geometric_and_population_constraints_pass
        and graph_integrity_pass
    )

    screening_penalty = (
        100.0
        * total_heavy_relative_error
        + 100.0
        * outer_radius_relative_error
        + 100.0
        * inner_radius_relative_error
        + 10.0
        * homopolar_bonds_required_for_30_site_attachment
    )

    return {
        "outer_shell_n": outer_shell,
        "inner_shell_m": inner_shell,
        "annulus_heavy_atoms": len(
            annulus_vertices
        ),
        "annulus_edges": len(
            annulus_edges
        ),
        "degree2_atoms": (
            degrees.get(
                2,
                0,
            )
        ),
        "degree3_atoms": (
            degrees.get(
                3,
                0,
            )
        ),
        "outer_boundary_atoms": len(
            outer_boundary
        ),
        "outer_boundary_color0": (
            outer_color_counts.get(
                0,
                0,
            )
        ),
        "outer_boundary_color1": (
            outer_color_counts.get(
                1,
                0,
            )
        ),
        "inner_boundary_atoms": len(
            inner_boundary
        ),
        "inner_boundary_color0": (
            inner_color_counts.get(
                0,
                0,
            )
        ),
        "inner_boundary_color1": (
            inner_color_counts.get(
                1,
                0,
            )
        ),
        "maximum_heteropolar_attachment_sites": (
            max_heteropolar_attachments
        ),
        "homopolar_bonds_required_for_30_site_attachment": (
            homopolar_bonds_required_for_30_site_attachment
        ),
        "unused_outer_boundary_sites_after_30_attachments": (
            unused_outer_boundary_sites
        ),
        "outer_radius_mean_nm": (
            outer_radius_mean_nm
        ),
        "outer_radius_maximum_nm": (
            outer_radius_maximum_nm
        ),
        "inner_radius_mean_nm": (
            inner_radius_mean_nm
        ),
        "outer_radius_relative_error": (
            outer_radius_relative_error
        ),
        "inner_radius_relative_error": (
            inner_radius_relative_error
        ),
        "seed_heavy_atoms_per_end": (
            EXPECTED_SEED_ATOMS_PER_END
        ),
        "total_added_heavy_atoms_per_end": (
            total_added_heavy_atoms_per_end
        ),
        "total_heavy_relative_error": (
            total_heavy_relative_error
        ),
        "connected_components": len(
            connected_components(
                annulus_adjacency
            )
        ),
        "bipartite": bipartite,
        "four_member_cycles": (
            count_four_cycles(
                annulus_adjacency
            )
        ),
        "graph_integrity_pass": (
            graph_integrity_pass
        ),
        "direct_30_site_heteropolar_attachment_possible": (
            direct_30_site_heteropolar_attachment_possible
        ),
        "geometric_and_population_constraints_pass": (
            geometric_and_population_constraints_pass
        ),
        "direct_candidate_pass": (
            direct_candidate_pass
        ),
        "screening_penalty": (
            screening_penalty
        ),
    }


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        PARENT_SUMMARY_CSV,
        CYCLE_AUDIT_SUMMARY_CSV,
        SEED_SUMMARY_CSV,
        SEED_GATES_CSV,
    ):
        require_file(required)

    parent = read_single_csv_row(
        PARENT_SUMMARY_CSV
    )

    cycle_audit = read_single_csv_row(
        CYCLE_AUDIT_SUMMARY_CSV
    )

    seed = read_single_csv_row(
        SEED_SUMMARY_CSV
    )

    seed_gates = read_csv_rows(
        SEED_GATES_CSV
    )

    if parent.get(
        "decision"
    ) != EXPECTED_PARENT_DECISION:
        raise RuntimeError(
            "Gate 3A is not in the accepted state."
        )

    if cycle_audit.get(
        "decision"
    ) != EXPECTED_CYCLE_DECISION:
        raise RuntimeError(
            "Gate 3C.1 does not contain the expected "
            "hexagonal-redesign decision."
        )

    if seed.get(
        "decision"
    ) != EXPECTED_SEED_DECISION:
        raise RuntimeError(
            "Gate 3D is not in the accepted state."
        )

    failed_seed_gates = [
        row.get(
            "gate",
            "",
        )
        for row in seed_gates
        if not parse_bool(
            row.get(
                "pass",
                "false",
            )
        )
    ]

    if failed_seed_gates:
        raise RuntimeError(
            "Gate 3D contains failed gates: "
            + " | ".join(
                failed_seed_gates
            )
        )

    bond_length_nm = parse_float(
        parent,
        "BN_bond_mean_nm",
    )

    parent_rim_radius_nm = parse_float(
        parent,
        "parent_rim_mean_radius_nm",
    )

    target_aperture_radius_nm = parse_float(
        parent,
        "target_aperture_radius_nm",
    )

    target_aperture_diameter_nm = parse_float(
        parent,
        "target_aperture_diameter_nm",
    )

    target_total_added_heavy_atoms_per_end = parse_float(
        parent,
        "estimated_monolayer_hBN_atoms_per_end",
    )

    if parse_int(
        seed,
        "added_atoms_lower",
    ) != EXPECTED_SEED_ATOMS_PER_END:
        raise RuntimeError(
            "Unexpected lower seed atom count."
        )

    if parse_int(
        seed,
        "added_atoms_upper",
    ) != EXPECTED_SEED_ATOMS_PER_END:
        raise RuntimeError(
            "Unexpected upper seed atom count."
        )

    template_rows = []

    for outer_shell in range(
        MIN_OUTER_SHELL,
        MAX_OUTER_SHELL + 1,
    ):
        for inner_shell in range(
            1,
            outer_shell,
        ):
            template_rows.append(
                build_annulus_template(
                    outer_shell,
                    inner_shell,
                    bond_length_nm,
                    parent_rim_radius_nm,
                    target_aperture_radius_nm,
                    target_total_added_heavy_atoms_per_end,
                )
            )

    template_rows.sort(
        key=lambda row: (
            not bool(
                row[
                    "direct_candidate_pass"
                ]
            ),
            float(
                row[
                    "screening_penalty"
                ]
            ),
            int(
                row[
                    "outer_shell_n"
                ]
            ),
            int(
                row[
                    "inner_shell_m"
                ]
            ),
        )
    )

    for rank, row in enumerate(
        template_rows,
        start=1,
    ):
        row[
            "overall_screening_rank"
        ] = rank

    write_csv(
        TEMPLATES_CSV,
        template_rows,
    )

    direct_candidates = [
        row
        for row in template_rows
        if bool(
            row[
                "direct_candidate_pass"
            ]
        )
    ]

    attachment_capable = [
        row
        for row in template_rows
        if bool(
            row[
                "direct_30_site_heteropolar_attachment_possible"
            ]
        )
    ]

    best_overall = min(
        template_rows,
        key=lambda row: float(
            row[
                "screening_penalty"
            ]
        ),
    )

    best_geometry_population = min(
        template_rows,
        key=lambda row: (
            float(
                row[
                    "total_heavy_relative_error"
                ]
            )
            + float(
                row[
                    "outer_radius_relative_error"
                ]
            )
            + float(
                row[
                    "inner_radius_relative_error"
                ]
            )
        ),
    )

    best_attachment_capable = (
        min(
            attachment_capable,
            key=lambda row: float(
                row[
                    "screening_penalty"
                ]
            ),
        )
        if attachment_capable
        else None
    )

    selected = (
        min(
            direct_candidates,
            key=lambda row: float(
                row[
                    "screening_penalty"
                ]
            ),
        )
        if direct_candidates
        else None
    )

    best_candidate_rows = [
        {
            "classification": (
                "BEST_OVERALL_PENALTY"
            ),
            **best_overall,
        },
        {
            "classification": (
                "BEST_GEOMETRY_AND_POPULATION"
            ),
            **best_geometry_population,
        },
    ]

    if best_attachment_capable is not None:
        best_candidate_rows.append(
            {
                "classification": (
                    "BEST_30_SITE_ATTACHMENT_CAPABLE"
                ),
                **best_attachment_capable,
            }
        )

    if selected is not None:
        best_candidate_rows.append(
            {
                "classification": (
                    "SELECTED_DIRECT_CANDIDATE"
                ),
                **selected,
            }
        )

    write_csv(
        BEST_CANDIDATES_CSV,
        best_candidate_rows,
    )

    audit_gates = {
        "Gate3A_parent_audit_is_accepted": (
            parent.get(
                "decision"
            )
            == EXPECTED_PARENT_DECISION
        ),
        "Gate3C1_hexagonal_redesign_path_is_active": (
            cycle_audit.get(
                "decision"
            )
            == EXPECTED_CYCLE_DECISION
        ),
        "Gate3D_edge_completion_seed_is_accepted": (
            seed.get(
                "decision"
            )
            == EXPECTED_SEED_DECISION
        ),
        "Gate3D_has_no_failed_gates": (
            len(
                failed_seed_gates
            )
            == 0
        ),
        "template_family_contains_45_annuli": (
            len(
                template_rows
            )
            == 45
        ),
        "all_templates_are_connected": all(
            int(
                row[
                    "connected_components"
                ]
            )
            == 1
            for row in template_rows
        ),
        "all_templates_are_bipartite": all(
            bool(
                row[
                    "bipartite"
                ]
            )
            for row in template_rows
        ),
        "all_templates_have_no_four_member_cycles": all(
            int(
                row[
                    "four_member_cycles"
                ]
            )
            == 0
            for row in template_rows
        ),
        "all_templates_have_only_degree2_or_degree3_atoms": all(
            bool(
                row[
                    "graph_integrity_pass"
                ]
            )
            for row in template_rows
        ),
        "screening_metrics_are_finite": all(
            all(
                math.isfinite(
                    float(
                        row[field]
                    )
                )
                for field in (
                    "outer_radius_relative_error",
                    "inner_radius_relative_error",
                    "total_heavy_relative_error",
                    "screening_penalty",
                )
            )
            for row in template_rows
        ),
    }

    failed_audit_gates = [
        name
        for name, passed
        in audit_gates.items()
        if not passed
    ]

    audit_integrity_pass = (
        len(
            failed_audit_gates
        )
        == 0
    )

    direct_attachment_feasible = (
        audit_integrity_pass
        and selected is not None
    )

    decision = (
        PASS_DECISION
        if direct_attachment_feasible
        else NEGATIVE_DECISION
    )

    required_next_step = (
        "BUILD_AND_VALIDATE_SELECTED_HEXAGONAL_BN_ANNULUS_ATTACHMENT_GRAPH"
        if direct_attachment_feasible
        else
        "EVALUATE_C3_PARENT_RIM_RECONSTRUCTION_AND_HYBRID_LINKER_CONTINGENCIES"
    )

    summary = {
        "decision": decision,
        "templates_screened": len(
            template_rows
        ),
        "direct_candidates_passing_all_constraints": (
            len(
                direct_candidates
            )
        ),
        "attachment_capable_templates": (
            len(
                attachment_capable
            )
        ),
        "target_BN_bond_length_nm": (
            bond_length_nm
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
        "target_total_added_heavy_atoms_per_end": (
            target_total_added_heavy_atoms_per_end
        ),
        "best_geometry_outer_shell": (
            best_geometry_population[
                "outer_shell_n"
            ]
        ),
        "best_geometry_inner_shell": (
            best_geometry_population[
                "inner_shell_m"
            ]
        ),
        "best_geometry_annulus_heavy_atoms": (
            best_geometry_population[
                "annulus_heavy_atoms"
            ]
        ),
        "best_geometry_total_added_heavy_atoms_per_end": (
            best_geometry_population[
                "total_added_heavy_atoms_per_end"
            ]
        ),
        "best_geometry_outer_boundary_atoms": (
            best_geometry_population[
                "outer_boundary_atoms"
            ]
        ),
        "best_geometry_outer_boundary_color0": (
            best_geometry_population[
                "outer_boundary_color0"
            ]
        ),
        "best_geometry_outer_boundary_color1": (
            best_geometry_population[
                "outer_boundary_color1"
            ]
        ),
        "best_geometry_maximum_heteropolar_attachments": (
            best_geometry_population[
                "maximum_heteropolar_attachment_sites"
            ]
        ),
        "best_geometry_homopolar_bonds_required": (
            best_geometry_population[
                "homopolar_bonds_required_for_30_site_attachment"
            ]
        ),
        "best_geometry_outer_radius_mean_nm": (
            best_geometry_population[
                "outer_radius_mean_nm"
            ]
        ),
        "best_geometry_inner_radius_mean_nm": (
            best_geometry_population[
                "inner_radius_mean_nm"
            ]
        ),
        "best_geometry_total_heavy_relative_error": (
            best_geometry_population[
                "total_heavy_relative_error"
            ]
        ),
        "best_geometry_outer_radius_relative_error": (
            best_geometry_population[
                "outer_radius_relative_error"
            ]
        ),
        "best_geometry_inner_radius_relative_error": (
            best_geometry_population[
                "inner_radius_relative_error"
            ]
        ),
        "best_attachment_capable_outer_shell": (
            ""
            if best_attachment_capable is None
            else best_attachment_capable[
                "outer_shell_n"
            ]
        ),
        "best_attachment_capable_inner_shell": (
            ""
            if best_attachment_capable is None
            else best_attachment_capable[
                "inner_shell_m"
            ]
        ),
        "best_attachment_capable_outer_radius_mean_nm": (
            ""
            if best_attachment_capable is None
            else best_attachment_capable[
                "outer_radius_mean_nm"
            ]
        ),
        "best_attachment_capable_inner_radius_mean_nm": (
            ""
            if best_attachment_capable is None
            else best_attachment_capable[
                "inner_radius_mean_nm"
            ]
        ),
        "best_attachment_capable_total_added_heavy_atoms": (
            ""
            if best_attachment_capable is None
            else best_attachment_capable[
                "total_added_heavy_atoms_per_end"
            ]
        ),
        "selected_outer_shell": (
            ""
            if selected is None
            else selected[
                "outer_shell_n"
            ]
        ),
        "selected_inner_shell": (
            ""
            if selected is None
            else selected[
                "inner_shell_m"
            ]
        ),
        "audit_integrity_pass": (
            audit_integrity_pass
        ),
        "direct_pure_BN_annulus_attachment_feasible": (
            direct_attachment_feasible
        ),
        "current_homopolar_degree2_seed_retained": True,
        "coordinate_generation_authorized": False,
        "molecular_topology_generation_authorized": False,
        "formal_charge_assignment_authorized": False,
        "force_field_parameterization_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "QM_authorized": False,
        "failed_audit_gates": (
            " | ".join(
                failed_audit_gates
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
            in audit_gates.items()
        ],
    )

    AUDIT_JSON.write_text(
        json.dumps(
            {
                "summary": summary,
                "audit_gates": (
                    audit_gates
                ),
                "best_geometry_and_population": (
                    best_geometry_population
                ),
                "best_attachment_capable": (
                    best_attachment_capable
                ),
                "selected_direct_candidate": (
                    selected
                ),
                "interpretation": {
                    "negative_result": (
                        "A negative feasibility decision means "
                        "that no pure hexagonal BN annulus in the "
                        "screened family simultaneously satisfies "
                        "30 heteropolar attachments, radius, pore, "
                        "and heavy-atom constraints."
                    ),
                    "not_a_stability_calculation": (
                        "This is a graph and lattice-size audit, "
                        "not an energy or synthetic-feasibility "
                        "calculation."
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
                "Gate3A_parent_summary"
            ),
            "file": relative(
                PARENT_SUMMARY_CSV
            ),
            "sha256": sha256(
                PARENT_SUMMARY_CSV
            ),
        },
        {
            "role": (
                "Gate3C1_cycle_audit_summary"
            ),
            "file": relative(
                CYCLE_AUDIT_SUMMARY_CSV
            ),
            "sha256": sha256(
                CYCLE_AUDIT_SUMMARY_CSV
            ),
        },
        {
            "role": (
                "Gate3D_seed_summary"
            ),
            "file": relative(
                SEED_SUMMARY_CSV
            ),
            "sha256": sha256(
                SEED_SUMMARY_CSV
            ),
        },
        {
            "role": (
                "Gate3D_seed_gates"
            ),
            "file": relative(
                SEED_GATES_CSV
            ),
            "sha256": sha256(
                SEED_GATES_CSV
            ),
        },
    ]

    write_csv(
        SOURCE_MANIFEST_CSV,
        source_rows,
    )

    audit_gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed
        in audit_gates.items()
    )

    attachment_text = (
        "No attachment-capable template was found."
        if best_attachment_capable is None
        else (
            f"Best attachment-capable template: "
            f"n={best_attachment_capable['outer_shell_n']}, "
            f"m={best_attachment_capable['inner_shell_m']}; "
            f"outer radius "
            f"{best_attachment_capable['outer_radius_mean_nm']:.6f} nm; "
            f"inner radius "
            f"{best_attachment_capable['inner_radius_mean_nm']:.6f} nm; "
            f"total added heavy atoms "
            f"{best_attachment_capable['total_added_heavy_atoms_per_end']}."
        )
    )

    REPORT_MD.write_text(
        f"""# R2 Hexagonal-Annulus Attachment Feasibility Audit

## Scope

This gate screens exact graph templates constructed as a hexagonal
benzenoid flake with a concentric inner flake removed.

No molecular coordinates, bond orders, partial charges, force-field
parameters, topology, minimization, MD, or QM calculation were
generated.

## Targets

- Parent-rim radius:
  **{parent_rim_radius_nm:.6f} nm**
- Aperture radius/diameter:
  **{target_aperture_radius_nm:.6f}/
  {target_aperture_diameter_nm:.6f} nm**
- BN bond-length reference:
  **{bond_length_nm:.6f} nm**
- Edge-completion atoms already added:
  **{EXPECTED_SEED_ATOMS_PER_END} per end**
- Total heavy-atom screening target:
  **{target_total_added_heavy_atoms_per_end:.3f} per end**
- Required heteropolar attachment sites:
  **{EXPECTED_ATTACHMENT_SITES_PER_END} per end**

## Template family

- Outer shell indices:
  **{MIN_OUTER_SHELL}–{MAX_OUTER_SHELL}**
- Templates screened:
  **{len(template_rows)}**
- Attachment-capable templates:
  **{len(attachment_capable)}**
- Templates satisfying every constraint:
  **{len(direct_candidates)}**

## Best geometry/population template

- Outer/inner shell:
  **{best_geometry_population['outer_shell_n']}/
  {best_geometry_population['inner_shell_m']}**
- Annulus heavy atoms:
  **{best_geometry_population['annulus_heavy_atoms']}**
- Total heavy atoms including the completion seed:
  **{best_geometry_population['total_added_heavy_atoms_per_end']}**
- Outer boundary:
  **{best_geometry_population['outer_boundary_atoms']} atoms**
- Outer sublattice populations:
  **{best_geometry_population['outer_boundary_color0']}/
  {best_geometry_population['outer_boundary_color1']}**
- Maximum heteropolar attachments:
  **{best_geometry_population['maximum_heteropolar_attachment_sites']}**
- Homopolar attachments required to reach 30:
  **{best_geometry_population['homopolar_bonds_required_for_30_site_attachment']}**
- Mean outer radius:
  **{best_geometry_population['outer_radius_mean_nm']:.6f} nm**
- Mean inner radius:
  **{best_geometry_population['inner_radius_mean_nm']:.6f} nm**
- Heavy/outer-radius/inner-radius relative errors:
  **{best_geometry_population['total_heavy_relative_error']:.6f}/
  {best_geometry_population['outer_radius_relative_error']:.6f}/
  {best_geometry_population['inner_radius_relative_error']:.6f}**

## Best 30-site attachment-capable template

{attachment_text}

## Audit gates

{audit_gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed audit-integrity gates:
  **{'NONE' if not failed_audit_gates else ' | '.join(failed_audit_gates)}**
- Direct pure-BN annulus attachment feasible:
  **{'YES' if direct_attachment_feasible else 'NO'}**
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

## Interpretation

A negative result does not invalidate the hexagonal edge-completion
seed. It shows that a closed, pure-hexagonal BN annulus cannot be
attached directly to the resulting elementally homogeneous degree-2
rim while simultaneously satisfying the 30-site heteropolar-junction,
radius, aperture, and atom-population constraints in the screened
template family.

The next comparison must therefore examine the previously retained
parent-rim reconstruction contingency and an explicitly end-specific
hybrid-linker architecture.
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 hexagonal-annulus attachment "
        "feasibility audit completed."
    )

    print(
        "Templates screened / attachment-capable / "
        "fully feasible: "
        f"{len(template_rows)}/"
        f"{len(attachment_capable)}/"
        f"{len(direct_candidates)}"
    )

    print(
        "Targets bond / outer radius / aperture radius / "
        "heavy atoms per end: "
        f"{bond_length_nm:.6f}/"
        f"{parent_rim_radius_nm:.6f}/"
        f"{target_aperture_radius_nm:.6f}/"
        f"{target_total_added_heavy_atoms_per_end:.3f}"
    )

    print(
        "Best geometry template n/m / annulus atoms / "
        "total atoms with seed: "
        f"{best_geometry_population['outer_shell_n']}/"
        f"{best_geometry_population['inner_shell_m']}/"
        f"{best_geometry_population['annulus_heavy_atoms']}/"
        f"{best_geometry_population['total_added_heavy_atoms_per_end']}"
    )

    print(
        "Best geometry outer boundary color0/color1 / "
        "maximum heteropolar attachments / "
        "required homopolar attachments: "
        f"{best_geometry_population['outer_boundary_color0']}/"
        f"{best_geometry_population['outer_boundary_color1']}/"
        f"{best_geometry_population['maximum_heteropolar_attachment_sites']}/"
        f"{best_geometry_population['homopolar_bonds_required_for_30_site_attachment']}"
    )

    print(
        "Best geometry outer/inner radius: "
        f"{best_geometry_population['outer_radius_mean_nm']:.6f}/"
        f"{best_geometry_population['inner_radius_mean_nm']:.6f} nm"
    )

    print(
        "Best geometry heavy/outer/inner relative errors: "
        f"{best_geometry_population['total_heavy_relative_error']:.6f}/"
        f"{best_geometry_population['outer_radius_relative_error']:.6f}/"
        f"{best_geometry_population['inner_radius_relative_error']:.6f}"
    )

    if best_attachment_capable is None:
        print(
            "Best 30-site attachment-capable template: NONE"
        )
    else:
        print(
            "Best 30-site attachment-capable template "
            "n/m / outer radius / inner radius / total atoms: "
            f"{best_attachment_capable['outer_shell_n']}/"
            f"{best_attachment_capable['inner_shell_m']}/"
            f"{best_attachment_capable['outer_radius_mean_nm']:.6f}/"
            f"{best_attachment_capable['inner_radius_mean_nm']:.6f}/"
            f"{best_attachment_capable['total_added_heavy_atoms_per_end']}"
        )

        print(
            "Best attachment-capable heavy/outer/inner "
            "relative errors: "
            f"{best_attachment_capable['total_heavy_relative_error']:.6f}/"
            f"{best_attachment_capable['outer_radius_relative_error']:.6f}/"
            f"{best_attachment_capable['inner_radius_relative_error']:.6f}"
        )

    print(
        f"Decision: {decision}"
    )

    print(
        "Failed audit-integrity gates: "
        + (
            "NONE"
            if not failed_audit_gates
            else " | ".join(
                failed_audit_gates
            )
        )
    )

    print(
        "Direct pure-BN annulus attachment feasible: "
        f"{'YES' if direct_attachment_feasible else 'NO'}"
    )

    print(
        "Current hexagonal edge-completion seed retained: YES"
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
        TEMPLATES_CSV,
        BEST_CANDIDATES_CSV,
        SUMMARY_CSV,
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
