#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/phase1A/day024_chemical_end_rim_design"

GATE3M = BASE / "16_r2_selected_full_density_longer_bn_bridge_graph"
GATE3P = BASE / "23_r2_four_atom_hydrogen_coordinate_embedding"
GATE3P1 = BASE / "24_r2_hydrogen_symmetry_refinement_preflight"

OUT = BASE / "26_r2_complete_end_symmetry_correspondence"

GRAPH_NODES = GATE3M / "r2_selected_longer_bn_bridge_graph_nodes.csv"
GRAPH_EDGES = GATE3M / "r2_selected_longer_bn_bridge_graph_edges.csv"

FULL_COORDINATES = (
    GATE3P
    / "r2_selected_four_atom_full_coordinates.csv"
)

PREVIOUS_HEAVY_PAIRS = (
    GATE3P1
    / "r2_heavy_lower_upper_pair_candidates.csv"
)

PREVIOUS_H_PAIRS = (
    GATE3P1
    / "r2_hydrogen_lower_upper_pair_candidates.csv"
)

COMPLETE_HEAVY_PAIRS = (
    OUT
    / "r2_complete_lower_upper_heavy_pairs.csv"
)

COMPLETE_H_PAIRS = (
    OUT
    / "r2_complete_lower_upper_hydrogen_pairs.csv"
)

TRANSFORM_SUMMARY = (
    OUT
    / "r2_end_symmetry_transform_comparison.csv"
)

SELECTED_TRANSFORM = (
    OUT
    / "r2_selected_annulus_symmetry_transform.csv"
)

DEVIATION_SUMMARY = (
    OUT
    / "r2_symmetry_deviation_by_node_type.csv"
)

SUMMARY = (
    OUT
    / "r2_complete_end_symmetry_correspondence_summary.csv"
)

GATES = (
    OUT
    / "r2_complete_end_symmetry_correspondence_gates.csv"
)

JSON_OUT = (
    OUT
    / "r2_complete_end_symmetry_correspondence.json"
)

REPORT = (
    OUT
    / "R2_COMPLETE_END_SYMMETRY_CORRESPONDENCE_DAY024.md"
)

PASS_DECISION = (
    "R2_COMPLETE_END_SYMMETRY_CORRESPONDENCE_VALIDATED"
)

REVIEW_DECISION = (
    "R2_COMPLETE_END_SYMMETRY_CORRESPONDENCE_REQUIRES_REVIEW"
)

EXPECTED_HEAVY_PER_END = 216
EXPECTED_H_PER_END = 102

ANNULUS_TYPES = {
    "ANNULUS_INTERIOR",
    "ANNULUS_OUTER_BOUNDARY",
    "ANNULUS_INNER_BOUNDARY",
}

RIGID_RIM_TYPES = {
    *ANNULUS_TYPES,
    "HEXAGONAL_EDGE_COMPLETION_SEED",
}

INNER_BOUNDARY_TYPES = {
    "ANNULUS_INNER_BOUNDARY",
}

MAX_RIGID_RIM_RMSD_NM = 1.0e-8
MAX_RIGID_RIM_DEVIATION_NM = 1.0e-7


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty required file: {path}"
        )


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
    value = float(row[key])

    if not math.isfinite(value):
        raise RuntimeError(
            f"Non-finite field {key!r}"
        )

    return value


def coordinates(
    rows: list[dict[str, str]],
) -> dict[str, np.ndarray]:
    return {
        row["node_id"]: np.asarray(
            [
                parse_float(row, "x_nm"),
                parse_float(row, "y_nm"),
                parse_float(row, "z_nm"),
            ],
            dtype=float,
        )
        for row in rows
    }


def orthogonal_fit(
    source: np.ndarray,
    target: np.ndarray,
    allow_reflection: bool,
) -> tuple[
    np.ndarray,
    np.ndarray,
    float,
    float,
]:
    source_center = np.mean(
        source,
        axis=0,
    )

    target_center = np.mean(
        target,
        axis=0,
    )

    source_centered = (
        source
        - source_center
    )

    target_centered = (
        target
        - target_center
    )

    covariance = (
        source_centered.T
        @ target_centered
    )

    u_matrix, _, vt_matrix = np.linalg.svd(
        covariance
    )

    rotation = (
        vt_matrix.T
        @ u_matrix.T
    )

    if (
        not allow_reflection
        and np.linalg.det(rotation) < 0.0
    ):
        vt_matrix[-1, :] *= -1.0

        rotation = (
            vt_matrix.T
            @ u_matrix.T
        )

    translation = (
        target_center
        - rotation
        @ source_center
    )

    transformed = (
        source
        @ rotation.T
        + translation
    )

    deviations = np.linalg.norm(
        transformed
        - target,
        axis=1,
    )

    rmsd = float(
        np.sqrt(
            np.mean(
                deviations ** 2
            )
        )
    )

    maximum = float(
        np.max(
            deviations
        )
    )

    return (
        rotation,
        translation,
        rmsd,
        maximum,
    )


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        GRAPH_NODES,
        GRAPH_EDGES,
        FULL_COORDINATES,
        PREVIOUS_HEAVY_PAIRS,
        PREVIOUS_H_PAIRS,
    ):
        require_file(required)

    node_rows = read_rows(
        GRAPH_NODES
    )

    edge_rows = read_rows(
        GRAPH_EDGES
    )

    coordinate_rows = read_rows(
        FULL_COORDINATES
    )

    previous_heavy_pair_rows = read_rows(
        PREVIOUS_HEAVY_PAIRS
    )

    nodes = {
        row["node_id"]: row
        for row in node_rows
    }

    positions = coordinates(
        coordinate_rows
    )

    adjacency = {
        node_id: set()
        for node_id in nodes
    }

    for row in edge_rows:
        first = row["source_node"]
        second = row["target_node"]

        adjacency[first].add(second)
        adjacency[second].add(first)

    lower_heavy = [
        row
        for row in node_rows
        if (
            row["element"] != "H"
            and row["end"] == "LOWER"
        )
    ]

    upper_heavy = [
        row
        for row in node_rows
        if (
            row["element"] != "H"
            and row["end"] == "UPPER"
        )
    ]

    lower_H = [
        row
        for row in node_rows
        if (
            row["element"] == "H"
            and row["end"] == "LOWER"
        )
    ]

    upper_H = [
        row
        for row in node_rows
        if (
            row["element"] == "H"
            and row["end"] == "UPPER"
        )
    ]

    heavy_pair_map = {}

    for row in previous_heavy_pair_rows:
        lower_id = row.get(
            "lower_node",
            "",
        )

        upper_id = row.get(
            "upper_node",
            "",
        )

        if lower_id and upper_id:
            heavy_pair_map[
                lower_id
            ] = upper_id

    upper_seed_by_index = {
        int(
            row[
                "circumferential_index"
            ]
        ): row["node_id"]
        for row in upper_heavy
        if row["node_type"]
        == "HEXAGONAL_EDGE_COMPLETION_SEED"
    }

    lower_seed_rows = [
        row
        for row in lower_heavy
        if row["node_type"]
        == "HEXAGONAL_EDGE_COMPLETION_SEED"
    ]

    if (
        len(lower_seed_rows) != 30
        or len(upper_seed_by_index) != 30
    ):
        raise RuntimeError(
            "Expected 30 seed nodes per end."
        )

    for lower_row in lower_seed_rows:
        index = int(
            lower_row[
                "circumferential_index"
            ]
        )

        if index not in upper_seed_by_index:
            raise RuntimeError(
                f"No upper seed for circumferential index {index}."
            )

        heavy_pair_map[
            lower_row["node_id"]
        ] = upper_seed_by_index[
            index
        ]

    if len(heavy_pair_map) != EXPECTED_HEAVY_PER_END:
        raise RuntimeError(
            "Incomplete heavy correspondence: "
            f"{len(heavy_pair_map)}"
        )

    if len(set(heavy_pair_map.values())) != EXPECTED_HEAVY_PER_END:
        raise RuntimeError(
            "Upper heavy correspondence is not one-to-one."
        )

    heavy_pair_rows = []

    for lower_id in sorted(
        heavy_pair_map
    ):
        upper_id = heavy_pair_map[
            lower_id
        ]

        lower_row = nodes[
            lower_id
        ]

        upper_row = nodes[
            upper_id
        ]

        if (
            lower_row["node_type"]
            != upper_row["node_type"]
        ):
            raise RuntimeError(
                "Incompatible heavy node types: "
                f"{lower_id} ({lower_row['node_type']}) / "
                f"{upper_id} ({upper_row['node_type']})"
            )

        lower_element = lower_row["element"]
        upper_element = upper_row["element"]

        if lower_element == upper_element:
            element_relation = "SAME_ELEMENT"

        elif {
            lower_element,
            upper_element,
        } == {
            "B",
            "N",
        }:
            element_relation = "B_N_SUBLATTICE_SWAP"

        else:
            raise RuntimeError(
                "Incompatible heavy-element relation: "
                f"{lower_id} ({lower_element}) / "
                f"{upper_id} ({upper_element})"
            )

        method = (
            "CIRCUMFERENTIAL_INDEX"
            if lower_row["node_type"]
            == "HEXAGONAL_EDGE_COMPLETION_SEED"
            else "PREVIOUS_VALID_PAIR"
        )

        heavy_pair_rows.append(
            {
                "lower_node": lower_id,
                "upper_node": upper_id,
                "lower_element": lower_element,
                "upper_element": upper_element,
                "element_relation": element_relation,
                "node_type": lower_row["node_type"],
                "pair_method": method,
                "circumferential_index": lower_row.get(
                    "circumferential_index",
                    "",
                ),
            }
        )

    write_rows(
        COMPLETE_HEAVY_PAIRS,
        heavy_pair_rows,
    )

    attached_heavy = {}

    for row in lower_H + upper_H:
        hydrogen_id = row["node_id"]

        heavy_neighbors = [
            neighbor
            for neighbor in adjacency[
                hydrogen_id
            ]
            if nodes[
                neighbor
            ][
                "element"
            ]
            != "H"
        ]

        if len(heavy_neighbors) != 1:
            raise RuntimeError(
                f"{hydrogen_id}: expected one attached heavy."
            )

        attached_heavy[
            hydrogen_id
        ] = heavy_neighbors[0]

    upper_H_by_parent_and_role = defaultdict(list)

    for row in upper_H:
        hydrogen_id = row["node_id"]

        upper_H_by_parent_and_role[
            (
                attached_heavy[
                    hydrogen_id
                ],
                row["node_type"],
            )
        ].append(
            hydrogen_id
        )

    H_pair_rows = []
    H_pair_map = {}

    for lower_row in lower_H:
        lower_id = lower_row[
            "node_id"
        ]

        lower_parent = attached_heavy[
            lower_id
        ]

        upper_parent = heavy_pair_map[
            lower_parent
        ]

        candidates = upper_H_by_parent_and_role[
            (
                upper_parent,
                lower_row["node_type"],
            )
        ]

        if len(candidates) != 1:
            raise RuntimeError(
                f"{lower_id}: expected one upper H candidate; "
                f"found {len(candidates)}"
            )

        upper_id = candidates[0]

        H_pair_map[
            lower_id
        ] = upper_id

        H_pair_rows.append(
            {
                "lower_H": lower_id,
                "upper_H": upper_id,
                "hydrogen_role": lower_row["node_type"],
                "lower_attached_heavy": lower_parent,
                "upper_attached_heavy": upper_parent,
                "pair_method": (
                    "ATTACHED_HEAVY_PAIR_AND_ROLE"
                ),
            }
        )

    if len(H_pair_map) != EXPECTED_H_PER_END:
        raise RuntimeError(
            "Incomplete H correspondence: "
            f"{len(H_pair_map)}"
        )

    if len(set(H_pair_map.values())) != EXPECTED_H_PER_END:
        raise RuntimeError(
            "Upper H correspondence is not one-to-one."
        )

    write_rows(
        COMPLETE_H_PAIRS,
        H_pair_rows,
    )

    subset_definitions = {
        "ALL_END_HEAVY": (
            lambda row: True
        ),
        "RIGID_RIM_HEAVY": (
            lambda row: row[
                "node_type"
            ] in RIGID_RIM_TYPES
        ),
        "ANNULUS_HEAVY": (
            lambda row: row[
                "node_type"
            ] in ANNULUS_TYPES
        ),
        "INNER_BOUNDARY_HEAVY": (
            lambda row: row[
                "node_type"
            ] in INNER_BOUNDARY_TYPES
        ),
        "SEED_HEAVY": (
            lambda row: row[
                "node_type"
            ]
            == "HEXAGONAL_EDGE_COMPLETION_SEED"
        ),
        "FOUR_ATOM_BRIDGE_HEAVY": (
            lambda row: row[
                "node_type"
            ]
            == "ALTERNATING_BN_FOUR_ATOM_BRIDGE"
        ),
    }

    transform_rows = []
    transform_objects = {}

    for subset_name, selector in subset_definitions.items():
        selected_rows = [
            row
            for row in heavy_pair_rows
            if selector(
                nodes[
                    row["lower_node"]
                ]
            )
        ]

        if len(selected_rows) < 3:
            continue

        source = np.asarray(
            [
                positions[
                    row["lower_node"]
                ]
                for row in selected_rows
            ],
            dtype=float,
        )

        target = np.asarray(
            [
                positions[
                    row["upper_node"]
                ]
                for row in selected_rows
            ],
            dtype=float,
        )

        for transform_type, allow_reflection in (
            (
                "PROPER_ROTATION_ONLY",
                False,
            ),
            (
                "ORTHOGONAL_REFLECTION_ALLOWED",
                True,
            ),
        ):
            (
                rotation,
                translation,
                rmsd,
                maximum,
            ) = orthogonal_fit(
                source,
                target,
                allow_reflection,
            )

            determinant = float(
                np.linalg.det(
                    rotation
                )
            )

            transform_rows.append(
                {
                    "subset": subset_name,
                    "transform_type": transform_type,
                    "pair_count": len(
                        selected_rows
                    ),
                    "determinant": determinant,
                    "RMSD_nm": rmsd,
                    "maximum_deviation_nm": maximum,
                }
            )

            transform_objects[
                (
                    subset_name,
                    transform_type,
                )
            ] = (
                rotation,
                translation,
                rmsd,
                maximum,
            )

    write_rows(
        TRANSFORM_SUMMARY,
        transform_rows,
    )

    candidate_keys = [
        (
            "RIGID_RIM_HEAVY",
            "ORTHOGONAL_REFLECTION_ALLOWED",
        ),
        (
            "ANNULUS_HEAVY",
            "ORTHOGONAL_REFLECTION_ALLOWED",
        ),
        (
            "INNER_BOUNDARY_HEAVY",
            "ORTHOGONAL_REFLECTION_ALLOWED",
        ),
        (
            "RIGID_RIM_HEAVY",
            "PROPER_ROTATION_ONLY",
        ),
        (
            "ANNULUS_HEAVY",
            "PROPER_ROTATION_ONLY",
        ),
    ]

    selected_key = min(
        candidate_keys,
        key=lambda key: (
            transform_objects[
                key
            ][2],
            transform_objects[
                key
            ][3],
        ),
    )

    (
        selected_rotation,
        selected_translation,
        selected_rmsd,
        selected_maximum,
    ) = transform_objects[
        selected_key
    ]

    selected_subset, selected_type = selected_key

    transform_matrix_rows = []

    for index in range(3):
        transform_matrix_rows.append(
            {
                "selected_subset": selected_subset,
                "selected_transform_type": selected_type,
                "row": index,
                "R0": float(
                    selected_rotation[
                        index,
                        0,
                    ]
                ),
                "R1": float(
                    selected_rotation[
                        index,
                        1,
                    ]
                ),
                "R2": float(
                    selected_rotation[
                        index,
                        2,
                    ]
                ),
                "translation_nm": float(
                    selected_translation[
                        index
                    ]
                ),
                "determinant": float(
                    np.linalg.det(
                        selected_rotation
                    )
                ),
                "subset_RMSD_nm": selected_rmsd,
                "subset_maximum_deviation_nm": (
                    selected_maximum
                ),
            }
        )

    write_rows(
        SELECTED_TRANSFORM,
        transform_matrix_rows,
    )

    deviations_by_type = defaultdict(list)

    for row in heavy_pair_rows:
        lower_id = row[
            "lower_node"
        ]

        upper_id = row[
            "upper_node"
        ]

        transformed = (
            selected_rotation
            @ positions[
                lower_id
            ]
            + selected_translation
        )

        deviation = float(
            np.linalg.norm(
                transformed
                - positions[
                    upper_id
                ]
            )
        )

        deviations_by_type[
            nodes[
                lower_id
            ][
                "node_type"
            ]
        ].append(
            deviation
        )

    deviation_rows = []

    for node_type in sorted(
        deviations_by_type
    ):
        values = np.asarray(
            deviations_by_type[
                node_type
            ],
            dtype=float,
        )

        deviation_rows.append(
            {
                "node_type": node_type,
                "pair_count": int(
                    values.size
                ),
                "RMSD_nm": float(
                    np.sqrt(
                        np.mean(
                            values ** 2
                        )
                    )
                ),
                "mean_deviation_nm": float(
                    np.mean(
                        values
                    )
                ),
                "maximum_deviation_nm": float(
                    np.max(
                        values
                    )
                ),
            }
        )

    write_rows(
        DEVIATION_SUMMARY,
        deviation_rows,
    )

    H_deviations_by_role = defaultdict(list)

    for row in H_pair_rows:
        lower_id = row[
            "lower_H"
        ]

        upper_id = row[
            "upper_H"
        ]

        transformed = (
            selected_rotation
            @ positions[
                lower_id
            ]
            + selected_translation
        )

        deviation = float(
            np.linalg.norm(
                transformed
                - positions[
                    upper_id
                ]
            )
        )

        H_deviations_by_role[
            row[
                "hydrogen_role"
            ]
        ].append(
            deviation
        )

    inner_H_values = np.asarray(
        H_deviations_by_role[
            "ANNULUS_INNER_PASSIVANT_H"
        ],
        dtype=float,
    )

    rigid_transform_is_exact = (
        selected_rmsd
        <= MAX_RIGID_RIM_RMSD_NM
        and selected_maximum
        <= MAX_RIGID_RIM_DEVIATION_NM
    )

    gates = {
        "216_lower_heavy_nodes_are_paired": (
            len(
                heavy_pair_map
            )
            == EXPECTED_HEAVY_PER_END
        ),
        "216_upper_heavy_nodes_are_used_once": (
            len(
                set(
                    heavy_pair_map.values()
                )
            )
            == EXPECTED_HEAVY_PER_END
        ),
        "30_seed_pairs_are_resolved_by_circumferential_index": (
            sum(
                row[
                    "node_type"
                ]
                == "HEXAGONAL_EDGE_COMPLETION_SEED"
                for row in heavy_pair_rows
            )
            == 30
        ),
        "102_lower_H_nodes_are_paired": (
            len(
                H_pair_map
            )
            == EXPECTED_H_PER_END
        ),
        "102_upper_H_nodes_are_used_once": (
            len(
                set(
                    H_pair_map.values()
                )
            )
            == EXPECTED_H_PER_END
        ),
        "15_seed_H_pairs_are_resolved_through_seed_correspondence": (
            sum(
                row[
                    "hydrogen_role"
                ]
                == "SEED_PASSIVANT_H"
                for row in H_pair_rows
            )
            == 15
        ),
        "12_inner_H_pairs_are_resolved": (
            len(
                inner_H_values
            )
            == 12
        ),
        "selected_rigid_rim_transform_is_exact": (
            rigid_transform_is_exact
        ),
        "coordinates_are_not_modified": True,
        "no_topology_minimization_MD_or_QM_is_generated": True,
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
        "REFINE_R2_HYDROGEN_ORIENTATIONS_WITH_"
        "VALIDATED_END_SYMMETRY_TRANSFORM"
        if accepted
        else
        "REVIEW_R2_END_SYMMETRY_TRANSFORM_AND_CORRESPONDENCE"
    )

    summary = {
        "decision": decision,
        "heavy_pairs": len(
            heavy_pair_map
        ),
        "seed_heavy_pairs": sum(
            row[
                "node_type"
            ]
            == "HEXAGONAL_EDGE_COMPLETION_SEED"
            for row in heavy_pair_rows
        ),
        "same_element_heavy_pairs": sum(
            row[
                "element_relation"
            ]
            == "SAME_ELEMENT"
            for row in heavy_pair_rows
        ),
        "B_N_sublattice_swap_heavy_pairs": sum(
            row[
                "element_relation"
            ]
            == "B_N_SUBLATTICE_SWAP"
            for row in heavy_pair_rows
        ),
        "H_pairs": len(
            H_pair_map
        ),
        "seed_H_pairs": sum(
            row[
                "hydrogen_role"
            ]
            == "SEED_PASSIVANT_H"
            for row in H_pair_rows
        ),
        "inner_H_pairs": int(
            inner_H_values.size
        ),
        "selected_subset": selected_subset,
        "selected_transform_type": selected_type,
        "selected_transform_determinant": float(
            np.linalg.det(
                selected_rotation
            )
        ),
        "selected_transform_RMSD_nm": (
            selected_rmsd
        ),
        "selected_transform_maximum_deviation_nm": (
            selected_maximum
        ),
        "inner_H_existing_symmetry_RMSD_nm": float(
            np.sqrt(
                np.mean(
                    inner_H_values ** 2
                )
            )
        ),
        "inner_H_existing_symmetry_maximum_nm": float(
            np.max(
                inner_H_values
            )
        ),
        "coordinates_modified": False,
        "molecular_topology_generated": False,
        "energy_minimized": False,
        "MD_performed": False,
        "QM_performed": False,
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
                "selected_rotation": (
                    selected_rotation.tolist()
                ),
                "selected_translation_nm": (
                    selected_translation.tolist()
                ),
                "gates": gates,
                "limitations": [
                    (
                        "This gate resolves correspondence and "
                        "audits symmetry only."
                    ),
                    (
                        "No heavy or hydrogen coordinates are modified."
                    ),
                    (
                        "Bridge conformers are evaluated separately "
                        "from the rigid annulus/seed reference."
                    ),
                    (
                        "No topology, charges, parameterization, "
                        "minimization, MD or QM calculation is generated."
                    ),
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    transform_lines = "\n".join(
        (
            f"- {row['subset']} / "
            f"{row['transform_type']}: "
            f"det={float(row['determinant']):.6f}; "
            f"RMSD={float(row['RMSD_nm']):.12e} nm; "
            f"max={float(row['maximum_deviation_nm']):.12e} nm"
        )
        for row in transform_rows
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
        f"""# R2 Complete End-Symmetry Correspondence

## Complete correspondence

- Heavy pairs: **{len(heavy_pair_map)}**
- Seed-heavy pairs: **{summary['seed_heavy_pairs']}**
- H pairs: **{len(H_pair_map)}**
- Seed-H pairs: **{summary['seed_H_pairs']}**
- Inner-H pairs: **{summary['inner_H_pairs']}**

## Transform comparison

{transform_lines}

## Selected rigid-rim transformation

- Subset: **{selected_subset}**
- Type: **{selected_type}**
- Determinant:
  **{summary['selected_transform_determinant']:.12f}**
- RMSD:
  **{selected_rmsd:.12e} nm**
- Maximum deviation:
  **{selected_maximum:.12e} nm**

## Existing inner-H asymmetry

- RMSD:
  **{summary['inner_H_existing_symmetry_RMSD_nm']:.12e} nm**
- Maximum:
  **{summary['inner_H_existing_symmetry_maximum_nm']:.12e} nm**

## Gates

{gate_lines}

## Decision

- Decision: **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Coordinates modified: **NO**
- Molecular topology generated: **NO**
- Energy minimization performed: **NO**
- MD performed: **NO**
- QM performed: **NO**
- Required next step:
  `{required_next_step}`
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 complete end-symmetry "
        "correspondence audit completed."
    )

    print(
        "Heavy pairs / seed-heavy pairs: "
        f"{len(heavy_pair_map)}/"
        f"{summary['seed_heavy_pairs']}"
    )

    print(
        "Heavy element relations same/sublattice-swap: "
        f"{summary['same_element_heavy_pairs']}/"
        f"{summary['B_N_sublattice_swap_heavy_pairs']}"
    )

    print(
        "H pairs / seed-H pairs / inner-H pairs: "
        f"{len(H_pair_map)}/"
        f"{summary['seed_H_pairs']}/"
        f"{summary['inner_H_pairs']}"
    )

    for row in transform_rows:
        print(
            f"{row['subset']} {row['transform_type']} "
            "det/RMSD/max: "
            f"{float(row['determinant']):.6f}/"
            f"{float(row['RMSD_nm']):.12e}/"
            f"{float(row['maximum_deviation_nm']):.12e}"
        )

    print(
        "Selected subset / transform type: "
        f"{selected_subset}/"
        f"{selected_type}"
    )

    print(
        "Selected transform determinant / RMSD / maximum: "
        f"{float(np.linalg.det(selected_rotation)):.12f}/"
        f"{selected_rmsd:.12e}/"
        f"{selected_maximum:.12e} nm"
    )

    print(
        "Existing inner-H symmetry RMSD / maximum: "
        f"{summary['inner_H_existing_symmetry_RMSD_nm']:.12e}/"
        f"{summary['inner_H_existing_symmetry_maximum_nm']:.12e} nm"
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
        "Coordinates modified: NO"
    )

    print(
        "Molecular topology generated: NO"
    )

    print(
        "Energy minimization performed: NO"
    )

    print(
        "MD performed: NO"
    )

    print(
        "QM performed: NO"
    )

    print(
        f"Required next step: {required_next_step}"
    )

    for path in (
        COMPLETE_HEAVY_PAIRS,
        COMPLETE_H_PAIRS,
        TRANSFORM_SUMMARY,
        SELECTED_TRANSFORM,
        DEVIATION_SUMMARY,
        SUMMARY,
        GATES,
        JSON_OUT,
        REPORT,
    ):
        print(
            f"Wrote: {path.relative_to(ROOT)}"
        )

    if not accepted:
        raise RuntimeError(
            "Complete end-symmetry correspondence requires review."
        )


if __name__ == "__main__":
    main()
