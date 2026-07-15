#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/phase1A/day024_chemical_end_rim_design"

GATE3M = (
    BASE
    / "16_r2_selected_full_density_longer_bn_bridge_graph"
)

GATE3P2 = (
    BASE
    / "28_r2_inner_h_reflected_direction_refinement"
)

GRAPH_NODES = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

GRAPH_EDGES = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_edges.csv"
)

GRAPH_SUMMARY = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_summary.csv"
)

REFINED_COORDINATES = (
    GATE3P2
    / "r2_selected_four_atom_refined_full_coordinates.csv"
)

REFINEMENT_SUMMARY = (
    GATE3P2
    / "r2_inner_h_refinement_summary.csv"
)

OUT = (
    BASE
    / "29_r2_chemical_realizability_and_parameterization_scope"
)

ATOM_INVENTORY = (
    OUT
    / "r2_chemical_atom_inventory.csv"
)

BOND_INVENTORY = (
    OUT
    / "r2_chemical_bond_inventory.csv"
)

BOND_STATISTICS = (
    OUT
    / "r2_chemical_bond_statistics.csv"
)

LOCAL_ENVIRONMENTS = (
    OUT
    / "r2_local_chemical_environment_inventory.csv"
)

NODE_ENVIRONMENTS = (
    OUT
    / "r2_node_chemical_environment_assignments.csv"
)

PARAMETERIZATION_SCOPE = (
    OUT
    / "r2_parameterization_scope.csv"
)

CRITICAL_CENTERS = (
    OUT
    / "r2_parameterization_critical_centers.csv"
)

SUMMARY = (
    OUT
    / "r2_chemical_realizability_and_parameterization_scope_summary.csv"
)

GATES = (
    OUT
    / "r2_chemical_realizability_and_parameterization_scope_gates.csv"
)

JSON_OUT = (
    OUT
    / "r2_chemical_realizability_and_parameterization_scope.json"
)

MANIFEST = (
    OUT
    / "r2_chemical_realizability_and_parameterization_scope_manifest.csv"
)

REPORT = (
    OUT
    / "R2_CHEMICAL_REALIZABILITY_AND_PARAMETERIZATION_SCOPE_DAY024.md"
)

EXPECTED_GRAPH_DECISION = (
    "R2_SELECTED_FULL_DENSITY_FOUR_ATOM_BN_BRIDGE_GRAPH_VALIDATED"
)

EXPECTED_REFINEMENT_DECISION = (
    "R2_SELECTED_FOUR_ATOM_BN_BRIDGE_"
    "HYDROGEN_COORDINATES_VALIDATED_AFTER_SYMMETRY_REFINEMENT"
)

PASS_DECISION = (
    "R2_STATIC_CHEMICAL_REALIZABILITY_VALIDATED_"
    "PARAMETERIZATION_SCOPE_DEFINED"
)

REVIEW_DECISION = (
    "R2_CHEMICAL_REALIZABILITY_OR_PARAMETERIZATION_SCOPE_REQUIRES_REVIEW"
)

EXPECTED_HEAVY = 2112
EXPECTED_H = 204
EXPECTED_TOTAL = 2316
EXPECTED_BH = 102
EXPECTED_NH = 102

BN_TARGET_NM = 0.144973
BH_TARGET_NM = 0.119
NH_TARGET_NM = 0.101

MAX_BN_DEVIATION_NM = 0.003
MAX_XH_DEVIATION_NM = 0.002

BRIDGE_TYPE = "ALTERNATING_BN_FOUR_ATOM_BRIDGE"
SEED_TYPE = "HEXAGONAL_EDGE_COMPLETION_SEED"

EDGE_TYPES = {
    "HEXAGONAL_EDGE_COMPLETION_SEED",
    "ANNULUS_OUTER_BOUNDARY",
    "ANNULUS_INNER_BOUNDARY",
}

ANNULUS_TYPES = {
    "ANNULUS_INTERIOR",
    "ANNULUS_OUTER_BOUNDARY",
    "ANNULUS_INNER_BOUNDARY",
}

PASSIVANT_TYPES = {
    "BRIDGE_PASSIVANT_H",
    "SEED_PASSIVANT_H",
    "ANNULUS_OUTER_PASSIVANT_H",
    "ANNULUS_INNER_PASSIVANT_H",
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty required file: {path}"
        )


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
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


def read_rows(path: Path) -> list[dict[str, str]]:
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


def read_one(path: Path) -> dict[str, str]:
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
            f"Non-finite value in {key!r}"
        )

    return value


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


def bond_class(
    first_element: str,
    second_element: str,
) -> str:
    pair = tuple(
        sorted(
            (
                first_element,
                second_element,
            )
        )
    )

    if pair == (
        "B",
        "N",
    ):
        return "B-N"

    if pair == (
        "B",
        "H",
    ):
        return "B-H"

    if pair == (
        "H",
        "N",
    ):
        return "N-H"

    return (
        f"{first_element}-{second_element}"
    )


def target_bond_length(
    class_name: str,
) -> float | None:
    return {
        "B-N": BN_TARGET_NM,
        "B-H": BH_TARGET_NM,
        "N-H": NH_TARGET_NM,
    }.get(
        class_name
    )


def classify_region(
    row: dict[str, str],
    adjacency: dict[str, set[str]],
    nodes: dict[str, dict[str, str]],
) -> str:
    node_id = row["node_id"]
    node_type = row["node_type"]

    if row["element"] == "H":
        if node_type == "BRIDGE_PASSIVANT_H":
            return "BRIDGE_PASSIVANT"

        if node_type in PASSIVANT_TYPES:
            return "EDGE_PASSIVANT"

        return "OTHER_HYDROGEN"

    if node_type == BRIDGE_TYPE:
        return "FOUR_ATOM_BRIDGE"

    if node_type == SEED_TYPE:
        return "SEED_EDGE"

    if node_type in {
        "ANNULUS_OUTER_BOUNDARY",
        "ANNULUS_INNER_BOUNDARY",
    }:
        return "ANNULUS_EDGE"

    if node_type == "ANNULUS_INTERIOR":
        bridge_neighbor = any(
            nodes[neighbor]["node_type"]
            == BRIDGE_TYPE
            for neighbor in adjacency[node_id]
        )

        return (
            "ANNULUS_BRIDGE_ATTACHMENT"
            if bridge_neighbor
            else "ANNULUS_INTERIOR"
        )

    if node_type == "PARENT_HBN":
        seed_neighbor = any(
            nodes[neighbor]["node_type"]
            == SEED_TYPE
            for neighbor in adjacency[node_id]
        )

        return (
            "PARENT_SEED_ATTACHMENT"
            if seed_neighbor
            else "PARENT_HBN_BULK_LIKE"
        )

    return "OTHER_HEAVY"


def environment_signature(
    node_id: str,
    nodes: dict[str, dict[str, str]],
    adjacency: dict[str, set[str]],
    region: str,
) -> str:
    row = nodes[node_id]

    neighbor_elements = Counter(
        nodes[neighbor]["element"]
        for neighbor in adjacency[node_id]
    )

    neighbor_types = Counter(
        nodes[neighbor]["node_type"]
        for neighbor in adjacency[node_id]
    )

    element_part = ",".join(
        f"{key}:{value}"
        for key, value
        in sorted(
            neighbor_elements.items()
        )
    )

    type_part = ",".join(
        f"{key}:{value}"
        for key, value
        in sorted(
            neighbor_types.items()
        )
    )

    return (
        f"center={row['element']}|"
        f"type={row['node_type']}|"
        f"region={region}|"
        f"degree={len(adjacency[node_id])}|"
        f"neighbor_elements={element_part}|"
        f"neighbor_types={type_part}"
    )


def graph_is_connected(
    node_ids: set[str],
    adjacency: dict[str, set[str]],
) -> bool:
    if not node_ids:
        return False

    start = next(
        iter(
            node_ids
        )
    )

    visited = {
        start
    }

    queue = deque(
        [
            start
        ]
    )

    while queue:
        current = queue.popleft()

        for neighbor in adjacency[current]:
            if neighbor not in visited:
                visited.add(
                    neighbor
                )

                queue.append(
                    neighbor
                )

    return visited == node_ids


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        GRAPH_NODES,
        GRAPH_EDGES,
        GRAPH_SUMMARY,
        REFINED_COORDINATES,
        REFINEMENT_SUMMARY,
    ):
        require_file(required)

    graph_summary = read_one(
        GRAPH_SUMMARY
    )

    refinement_summary = read_one(
        REFINEMENT_SUMMARY
    )

    if graph_summary.get(
        "decision"
    ) != EXPECTED_GRAPH_DECISION:
        raise RuntimeError(
            "Gate 3M graph is not accepted."
        )

    if refinement_summary.get(
        "decision"
    ) != EXPECTED_REFINEMENT_DECISION:
        raise RuntimeError(
            "Gate 3P.2 refined coordinates are not accepted."
        )

    node_rows = read_rows(
        GRAPH_NODES
    )

    edge_rows = read_rows(
        GRAPH_EDGES
    )

    coordinate_rows = read_rows(
        REFINED_COORDINATES
    )

    nodes = {
        row["node_id"]: row
        for row in node_rows
    }

    positions = {
        row["node_id"]: np.asarray(
            [
                parse_float(
                    row,
                    "x_nm",
                ),
                parse_float(
                    row,
                    "y_nm",
                ),
                parse_float(
                    row,
                    "z_nm",
                ),
            ],
            dtype=float,
        )
        for row in coordinate_rows
    }

    node_ids = set(
        nodes
    )

    coordinate_ids = set(
        positions
    )

    adjacency = {
        node_id: set()
        for node_id in node_ids
    }

    seen_edges = set()
    duplicate_edges = 0
    self_edges = 0
    missing_edge_nodes = 0

    bond_rows = []

    for edge_index, row in enumerate(
        edge_rows,
        start=1,
    ):
        first = row["source_node"]
        second = row["target_node"]

        if (
            first not in nodes
            or second not in nodes
        ):
            missing_edge_nodes += 1
            continue

        if first == second:
            self_edges += 1
            continue

        edge_key = tuple(
            sorted(
                (
                    first,
                    second,
                )
            )
        )

        if edge_key in seen_edges:
            duplicate_edges += 1
            continue

        seen_edges.add(
            edge_key
        )

        adjacency[first].add(
            second
        )

        adjacency[second].add(
            first
        )

        first_element = nodes[
            first
        ][
            "element"
        ]

        second_element = nodes[
            second
        ][
            "element"
        ]

        class_name = bond_class(
            first_element,
            second_element,
        )

        distance = float(
            np.linalg.norm(
                positions[first]
                - positions[second]
            )
        )

        target = target_bond_length(
            class_name
        )

        deviation = (
            abs(
                distance
                - target
            )
            if target is not None
            else ""
        )

        bond_rows.append(
            {
                "edge_index": edge_index,
                "source_node": first,
                "target_node": second,
                "source_element": first_element,
                "target_element": second_element,
                "source_type": nodes[
                    first
                ][
                    "node_type"
                ],
                "target_type": nodes[
                    second
                ][
                    "node_type"
                ],
                "bond_class": class_name,
                "distance_nm": distance,
                "target_nm": (
                    ""
                    if target is None
                    else target
                ),
                "absolute_deviation_nm": (
                    deviation
                ),
            }
        )

    write_rows(
        BOND_INVENTORY,
        bond_rows,
    )

    element_counts = Counter(
        row["element"]
        for row in node_rows
    )

    node_type_counts = Counter(
        row["node_type"]
        for row in node_rows
    )

    inventory_rows = []

    for element in sorted(
        element_counts
    ):
        inventory_rows.append(
            {
                "inventory_class": "ELEMENT",
                "name": element,
                "count": element_counts[
                    element
                ],
            }
        )

    for node_type in sorted(
        node_type_counts
    ):
        inventory_rows.append(
            {
                "inventory_class": "NODE_TYPE",
                "name": node_type,
                "count": node_type_counts[
                    node_type
                ],
            }
        )

    write_rows(
        ATOM_INVENTORY,
        inventory_rows,
    )

    bond_class_values = defaultdict(
        list
    )

    for row in bond_rows:
        bond_class_values[
            row[
                "bond_class"
            ]
        ].append(
            float(
                row[
                    "distance_nm"
                ]
            )
        )

    bond_statistics_rows = []

    for class_name in sorted(
        bond_class_values
    ):
        values = np.asarray(
            bond_class_values[
                class_name
            ],
            dtype=float,
        )

        target = target_bond_length(
            class_name
        )

        bond_statistics_rows.append(
            {
                "bond_class": class_name,
                "count": int(
                    values.size
                ),
                "target_nm": (
                    ""
                    if target is None
                    else target
                ),
                "minimum_nm": float(
                    np.min(
                        values
                    )
                ),
                "mean_nm": float(
                    np.mean(
                        values
                    )
                ),
                "maximum_nm": float(
                    np.max(
                        values
                    )
                ),
                "maximum_absolute_deviation_nm": (
                    ""
                    if target is None
                    else float(
                        np.max(
                            np.abs(
                                values
                                - target
                            )
                        )
                    )
                ),
            }
        )

    write_rows(
        BOND_STATISTICS,
        bond_statistics_rows,
    )

    node_assignment_rows = []
    signature_counts = Counter()
    signature_example = {}

    region_counts = Counter()

    for node_id in sorted(
        nodes
    ):
        row = nodes[
            node_id
        ]

        region = classify_region(
            row,
            adjacency,
            nodes,
        )

        signature = environment_signature(
            node_id,
            nodes,
            adjacency,
            region,
        )

        signature_counts[
            signature
        ] += 1

        signature_example.setdefault(
            signature,
            node_id,
        )

        region_counts[
            region
        ] += 1

        neighbor_elements = Counter(
            nodes[neighbor][
                "element"
            ]
            for neighbor in adjacency[
                node_id
            ]
        )

        node_assignment_rows.append(
            {
                "node_id": node_id,
                "element": row[
                    "element"
                ],
                "node_type": row[
                    "node_type"
                ],
                "end": row[
                    "end"
                ],
                "region": region,
                "total_degree": len(
                    adjacency[
                        node_id
                    ]
                ),
                "B_neighbors": neighbor_elements[
                    "B"
                ],
                "N_neighbors": neighbor_elements[
                    "N"
                ],
                "H_neighbors": neighbor_elements[
                    "H"
                ],
                "environment_signature": signature,
            }
        )

    write_rows(
        NODE_ENVIRONMENTS,
        node_assignment_rows,
    )

    environment_rows = []

    for environment_index, signature in enumerate(
        sorted(
            signature_counts
        ),
        start=1,
    ):
        example_id = signature_example[
            signature
        ]

        example = next(
            row
            for row in node_assignment_rows
            if row["node_id"]
            == example_id
        )

        environment_rows.append(
            {
                "environment_id": (
                    f"ENV_{environment_index:04d}"
                ),
                "count": signature_counts[
                    signature
                ],
                "example_node": example_id,
                "element": example[
                    "element"
                ],
                "node_type": example[
                    "node_type"
                ],
                "region": example[
                    "region"
                ],
                "total_degree": example[
                    "total_degree"
                ],
                "environment_signature": signature,
            }
        )

    write_rows(
        LOCAL_ENVIRONMENTS,
        environment_rows,
    )

    bridge_ids = {
        node_id
        for node_id, row
        in nodes.items()
        if row["node_type"]
        == BRIDGE_TYPE
    }

    bridge_first_shell = set()

    for bridge_id in bridge_ids:
        for neighbor in adjacency[
            bridge_id
        ]:
            if nodes[
                neighbor
            ][
                "element"
            ] != "H":
                bridge_first_shell.add(
                    neighbor
                )

    critical_rows = []

    for node_id in sorted(
        nodes
    ):
        row = nodes[
            node_id
        ]

        reasons = []

        if row[
            "node_type"
        ] == BRIDGE_TYPE:
            reasons.append(
                "FOUR_ATOM_BRIDGE_CENTER"
            )

        if node_id in bridge_first_shell:
            reasons.append(
                "FIRST_SHELL_OF_BRIDGE"
            )

        if row[
            "node_type"
        ] in EDGE_TYPES:
            reasons.append(
                "EDGE_OR_RIM_CENTER"
            )

        if row[
            "node_type"
        ] in PASSIVANT_TYPES:
            reasons.append(
                "HYDROGEN_PASSIVATED_EDGE_OR_BRIDGE"
            )

        if reasons:
            assignment = next(
                item
                for item in node_assignment_rows
                if item["node_id"]
                == node_id
            )

            critical_rows.append(
                {
                    "node_id": node_id,
                    "element": row[
                        "element"
                    ],
                    "node_type": row[
                        "node_type"
                    ],
                    "end": row[
                        "end"
                    ],
                    "region": assignment[
                        "region"
                    ],
                    "reasons": (
                        " | ".join(
                            reasons
                        )
                    ),
                    "parameterization_status": (
                        "UNASSESSED"
                    ),
                    "reference_data_required": True,
                }
            )

    write_rows(
        CRITICAL_CENTERS,
        critical_rows,
    )

    scope_rows = [
        {
            "scope_region": "PARENT_HBN_BULK_LIKE",
            "node_count": region_counts[
                "PARENT_HBN_BULK_LIKE"
            ],
            "chemical_description": (
                "Three-coordinate curved h-BN parent scaffold "
                "away from reconstructed ends."
            ),
            "expected_parameterization_route": (
                "Audit transferable bulk/curved h-BN bonded and "
                "nonbonded parameters against published references."
            ),
            "novel_parameterization_expected": (
                "POSSIBLY_NO_FOR_BULK_INTERIOR"
            ),
            "QM_reference_priority": "MEDIUM",
            "force_field_assignment_authorized": False,
        },
        {
            "scope_region": "PARENT_AND_ANNULUS_ATTACHMENTS",
            "node_count": (
                region_counts[
                    "PARENT_SEED_ATTACHMENT"
                ]
                + region_counts[
                    "ANNULUS_BRIDGE_ATTACHMENT"
                ]
            ),
            "chemical_description": (
                "Parent/seed and annulus/bridge junction atoms."
            ),
            "expected_parameterization_route": (
                "Explicit local bonded-environment comparison and "
                "small-cluster QM validation."
            ),
            "novel_parameterization_expected": "LIKELY",
            "QM_reference_priority": "HIGH",
            "force_field_assignment_authorized": False,
        },
        {
            "scope_region": "ANNULUS_AND_SEED_EDGE",
            "node_count": (
                region_counts[
                    "ANNULUS_INTERIOR"
                ]
                + region_counts[
                    "ANNULUS_EDGE"
                ]
                + region_counts[
                    "SEED_EDGE"
                ]
            ),
            "chemical_description": (
                "Reconstructed end annulus, inner/outer edge and "
                "hexagonal completion seed."
            ),
            "expected_parameterization_route": (
                "Audit edge-specific B/N/H types, equilibrium angles, "
                "torsions and nonbonded terms."
            ),
            "novel_parameterization_expected": "LIKELY",
            "QM_reference_priority": "HIGH",
            "force_field_assignment_authorized": False,
        },
        {
            "scope_region": "FOUR_ATOM_BRIDGE",
            "node_count": region_counts[
                "FOUR_ATOM_BRIDGE"
            ],
            "chemical_description": (
                "Alternating B-N-B-N four-atom bridge paths."
            ),
            "expected_parameterization_route": (
                "Dedicated constrained-cluster QM scans or optimized "
                "reference fragments for bonds, angles and torsions."
            ),
            "novel_parameterization_expected": "YES",
            "QM_reference_priority": "HIGHEST",
            "force_field_assignment_authorized": False,
        },
        {
            "scope_region": "HYDROGEN_PASSIVANTS",
            "node_count": (
                region_counts[
                    "BRIDGE_PASSIVANT"
                ]
                + region_counts[
                    "EDGE_PASSIVANT"
                ]
            ),
            "chemical_description": (
                "B-H and N-H passivation at bridge and reconstructed "
                "end environments."
            ),
            "expected_parameterization_route": (
                "Audit B-H/N-H bonded terms, partial charges and "
                "water-interaction behavior for each local environment."
            ),
            "novel_parameterization_expected": "LIKELY",
            "QM_reference_priority": "HIGH",
            "force_field_assignment_authorized": False,
        },
    ]

    write_rows(
        PARAMETERIZATION_SCOPE,
        scope_rows,
    )

    heavy_ids = {
        node_id
        for node_id, row
        in nodes.items()
        if row[
            "element"
        ] != "H"
    }

    H_ids = {
        node_id
        for node_id, row
        in nodes.items()
        if row[
            "element"
        ] == "H"
    }

    heavy_valence_failures = [
        node_id
        for node_id in heavy_ids
        if len(
            adjacency[
                node_id
            ]
        )
        != 3
    ]

    H_valence_failures = [
        node_id
        for node_id in H_ids
        if len(
            adjacency[
                node_id
            ]
        )
        != 1
    ]

    nonheteropolar_heavy_edges = [
        row
        for row in bond_rows
        if (
            row[
                "source_element"
            ]
            != "H"
            and row[
                "target_element"
            ]
            != "H"
            and row[
                "bond_class"
            ]
            != "B-N"
        )
    ]

    invalid_H_edges = [
        row
        for row in bond_rows
        if (
            "H"
            in {
                row[
                    "source_element"
                ],
                row[
                    "target_element"
                ],
            }
            and row[
                "bond_class"
            ]
            not in {
                "B-H",
                "N-H",
            }
        )
    ]

    bond_stats_by_class = {
        row[
            "bond_class"
        ]: row
        for row in bond_statistics_rows
    }

    maximum_BN_deviation = float(
        bond_stats_by_class[
            "B-N"
        ][
            "maximum_absolute_deviation_nm"
        ]
    )

    maximum_BH_deviation = float(
        bond_stats_by_class[
            "B-H"
        ][
            "maximum_absolute_deviation_nm"
        ]
    )

    maximum_NH_deviation = float(
        bond_stats_by_class[
            "N-H"
        ][
            "maximum_absolute_deviation_nm"
        ]
    )

    gates = {
        "Gate3M_graph_is_accepted": (
            graph_summary.get(
                "decision"
            )
            == EXPECTED_GRAPH_DECISION
        ),
        "Gate3P2_refined_coordinates_are_accepted": (
            refinement_summary.get(
                "decision"
            )
            == EXPECTED_REFINEMENT_DECISION
        ),
        "graph_contains_2316_nodes": (
            len(
                nodes
            )
            == EXPECTED_TOTAL
        ),
        "coordinates_exist_for_all_2316_nodes": (
            coordinate_ids
            == node_ids
        ),
        "graph_contains_2112_heavy_and_204_H_nodes": (
            len(
                heavy_ids
            )
            == EXPECTED_HEAVY
            and len(
                H_ids
            )
            == EXPECTED_H
        ),
        "graph_has_no_missing_edge_nodes": (
            missing_edge_nodes
            == 0
        ),
        "graph_has_no_self_edges": (
            self_edges
            == 0
        ),
        "graph_has_no_duplicate_edges": (
            duplicate_edges
            == 0
        ),
        "graph_is_connected": (
            graph_is_connected(
                node_ids,
                adjacency,
            )
        ),
        "all_heavy_atoms_are_three_coordinate": (
            len(
                heavy_valence_failures
            )
            == 0
        ),
        "all_H_atoms_are_one_coordinate": (
            len(
                H_valence_failures
            )
            == 0
        ),
        "all_heavy_heavy_edges_are_BN": (
            len(
                nonheteropolar_heavy_edges
            )
            == 0
        ),
        "all_H_edges_are_BH_or_NH": (
            len(
                invalid_H_edges
            )
            == 0
        ),
        "102_BH_and_102_NH_bonds_are_present": (
            int(
                bond_stats_by_class[
                    "B-H"
                ][
                    "count"
                ]
            )
            == EXPECTED_BH
            and int(
                bond_stats_by_class[
                    "N-H"
                ][
                    "count"
                ]
            )
            == EXPECTED_NH
        ),
        "BN_bond_deviation_is_at_most0p003nm": (
            maximum_BN_deviation
            <= MAX_BN_DEVIATION_NM
        ),
        "BH_and_NH_deviations_are_at_most0p002nm": (
            maximum_BH_deviation
            <= MAX_XH_DEVIATION_NM
            and maximum_NH_deviation
            <= MAX_XH_DEVIATION_NM
        ),
        "local_chemical_environments_are_fully_enumerated": (
            sum(
                row[
                    "count"
                ]
                for row in environment_rows
            )
            == EXPECTED_TOTAL
        ),
        "parameterization_critical_centers_are_explicitly_identified": (
            len(
                critical_rows
            )
            > 0
        ),
        "force_field_coverage_is_not_assumed": all(
            row[
                "force_field_assignment_authorized"
            ]
            is False
            for row in scope_rows
        ),
        "no_topology_charges_parameters_minimization_MD_or_QM_generated": True,
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
        "AUDIT_BN_H_FORCE_FIELD_COVERAGE_AND_DEFINE_"
        "QM_REFERENCE_FRAGMENT_SET"
        if accepted
        else
        "REVIEW_R2_STATIC_CHEMICAL_REALIZABILITY_FAILURES"
    )

    summary = {
        "decision": decision,
        "total_nodes": len(
            nodes
        ),
        "heavy_nodes": len(
            heavy_ids
        ),
        "H_nodes": len(
            H_ids
        ),
        "B_nodes": element_counts[
            "B"
        ],
        "N_nodes": element_counts[
            "N"
        ],
        "graph_edges": len(
            seen_edges
        ),
        "B_N_bonds": int(
            bond_stats_by_class[
                "B-N"
            ][
                "count"
            ]
        ),
        "B_H_bonds": int(
            bond_stats_by_class[
                "B-H"
            ][
                "count"
            ]
        ),
        "N_H_bonds": int(
            bond_stats_by_class[
                "N-H"
            ][
                "count"
            ]
        ),
        "maximum_BN_deviation_nm": (
            maximum_BN_deviation
        ),
        "maximum_BH_deviation_nm": (
            maximum_BH_deviation
        ),
        "maximum_NH_deviation_nm": (
            maximum_NH_deviation
        ),
        "heavy_valence_failure_count": len(
            heavy_valence_failures
        ),
        "H_valence_failure_count": len(
            H_valence_failures
        ),
        "nonheteropolar_heavy_edge_count": len(
            nonheteropolar_heavy_edges
        ),
        "invalid_H_edge_count": len(
            invalid_H_edges
        ),
        "unique_local_environment_count": len(
            environment_rows
        ),
        "parameterization_critical_center_count": len(
            critical_rows
        ),
        "static_chemical_graph_realizable": (
            accepted
        ),
        "existing_force_field_coverage_established": False,
        "molecular_topology_generation_authorized": False,
        "formal_charge_assignment_authorized": False,
        "force_field_parameterization_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "QM_authorized": False,
        "force_field_and_QM_reference_audit_authorized": (
            accepted
        ),
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

    JSON_OUT.write_text(
        json.dumps(
            {
                "summary": summary,
                "gates": gates,
                "element_counts": dict(
                    element_counts
                ),
                "node_type_counts": dict(
                    node_type_counts
                ),
                "region_counts": dict(
                    region_counts
                ),
                "limitations": [
                    (
                        "This gate establishes static graph and "
                        "coordinate realizability only."
                    ),
                    (
                        "No force-field coverage is inferred from "
                        "elemental composition alone."
                    ),
                    (
                        "Bridge, attachment, reconstructed-edge and "
                        "passivated environments require explicit "
                        "reference-data assessment."
                    ),
                    (
                        "No molecular topology, formal charges, "
                        "force-field parameters, minimization, MD or QM "
                        "calculation is generated."
                    ),
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    write_rows(
        MANIFEST,
        [
            {
                "role": "Gate3M_graph_nodes",
                "file": relative(
                    GRAPH_NODES
                ),
                "sha256": sha256(
                    GRAPH_NODES
                ),
            },
            {
                "role": "Gate3M_graph_edges",
                "file": relative(
                    GRAPH_EDGES
                ),
                "sha256": sha256(
                    GRAPH_EDGES
                ),
            },
            {
                "role": "Gate3M_graph_summary",
                "file": relative(
                    GRAPH_SUMMARY
                ),
                "sha256": sha256(
                    GRAPH_SUMMARY
                ),
            },
            {
                "role": "Gate3P2_refined_coordinates",
                "file": relative(
                    REFINED_COORDINATES
                ),
                "sha256": sha256(
                    REFINED_COORDINATES
                ),
            },
            {
                "role": "Gate3P2_refinement_summary",
                "file": relative(
                    REFINEMENT_SUMMARY
                ),
                "sha256": sha256(
                    REFINEMENT_SUMMARY
                ),
            },
        ],
    )

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed
        in gates.items()
    )

    scope_lines = "\n".join(
        (
            f"- **{row['scope_region']}**: "
            f"{row['node_count']} nodes; "
            f"QM priority={row['QM_reference_priority']}; "
            f"novel parameterization="
            f"{row['novel_parameterization_expected']}."
        )
        for row in scope_rows
    )

    REPORT.write_text(
        f"""# R2 Chemical Realizability and Parameterization Scope

## Static chemical inventory

- Total nodes: **{len(nodes)}**
- Heavy/H: **{len(heavy_ids)}/{len(H_ids)}**
- B/N: **{element_counts['B']}/{element_counts['N']}**
- Graph edges: **{len(seen_edges)}**
- B-N/B-H/N-H bonds:
  **{summary['B_N_bonds']}/{summary['B_H_bonds']}/{summary['N_H_bonds']}**

## Bond geometry

- Maximum B-N deviation:
  **{maximum_BN_deviation:.9f} nm**
- Maximum B-H deviation:
  **{maximum_BH_deviation:.12e} nm**
- Maximum N-H deviation:
  **{maximum_NH_deviation:.12e} nm**

## Valence and graph chemistry

- Heavy valence failures:
  **{len(heavy_valence_failures)}**
- H valence failures:
  **{len(H_valence_failures)}**
- Nonheteropolar heavy-heavy edges:
  **{len(nonheteropolar_heavy_edges)}**
- Invalid H edges:
  **{len(invalid_H_edges)}**

## Chemical environments

- Unique local environments:
  **{len(environment_rows)}**
- Parameterization-critical centers:
  **{len(critical_rows)}**

## Parameterization scope

{scope_lines}

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Static chemical graph realizable:
  **{'YES' if accepted else 'NO'}**
- Existing force-field coverage established:
  **NO**
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
- Force-field and QM-reference audit authorized:
  **{'YES' if accepted else 'NO'}**
- Required next step:
  `{required_next_step}`
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 chemical realizability and "
        "parameterization-scope audit completed."
    )

    print(
        "Nodes heavy/H/total and B/N: "
        f"{len(heavy_ids)}/"
        f"{len(H_ids)}/"
        f"{len(nodes)} and "
        f"{element_counts['B']}/"
        f"{element_counts['N']}"
    )

    print(
        "B-N / B-H / N-H bond counts: "
        f"{summary['B_N_bonds']}/"
        f"{summary['B_H_bonds']}/"
        f"{summary['N_H_bonds']}"
    )

    print(
        "Maximum B-N / B-H / N-H deviations: "
        f"{maximum_BN_deviation:.9f}/"
        f"{maximum_BH_deviation:.12e}/"
        f"{maximum_NH_deviation:.12e} nm"
    )

    print(
        "Heavy/H valence failures: "
        f"{len(heavy_valence_failures)}/"
        f"{len(H_valence_failures)}"
    )

    print(
        "Nonheteropolar heavy edges / invalid H edges: "
        f"{len(nonheteropolar_heavy_edges)}/"
        f"{len(invalid_H_edges)}"
    )

    print(
        "Unique local environments / critical centers: "
        f"{len(environment_rows)}/"
        f"{len(critical_rows)}"
    )

    print(
        "Static chemical graph realizable: "
        f"{accepted}"
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
        "Existing force-field coverage established: NO"
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
        "Force-field and QM-reference audit authorized: "
        f"{'YES' if accepted else 'NO'}"
    )

    print(
        f"Required next step: {required_next_step}"
    )

    for path in (
        ATOM_INVENTORY,
        BOND_INVENTORY,
        BOND_STATISTICS,
        LOCAL_ENVIRONMENTS,
        NODE_ENVIRONMENTS,
        PARAMETERIZATION_SCOPE,
        CRITICAL_CENTERS,
        SUMMARY,
        GATES,
        JSON_OUT,
        MANIFEST,
        REPORT,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not accepted:
        raise RuntimeError(
            "Chemical realizability or parameterization scope "
            "requires review."
        )


if __name__ == "__main__":
    main()
