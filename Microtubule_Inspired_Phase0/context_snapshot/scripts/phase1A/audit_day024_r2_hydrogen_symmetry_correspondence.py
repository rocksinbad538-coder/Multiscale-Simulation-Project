#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/phase1A/day024_chemical_end_rim_design"

GATE3M = BASE / "16_r2_selected_full_density_longer_bn_bridge_graph"
GATE3O = BASE / "19_r2_selected_four_atom_heavy_coordinate_embedding"
GATE3P = BASE / "23_r2_four_atom_hydrogen_coordinate_embedding"

GRAPH_NODES = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

GRAPH_EDGES = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_edges.csv"
)

HEAVY_COORDINATES = (
    GATE3O
    / "r2_selected_four_atom_heavy_coordinates.csv"
)

FULL_COORDINATES = (
    GATE3P
    / "r2_selected_four_atom_full_coordinates.csv"
)

H_ORIENTATIONS = (
    GATE3P
    / "r2_selected_four_atom_hydrogen_orientations.csv"
)

H_SUMMARY = (
    GATE3P
    / "r2_selected_four_atom_hydrogen_embedding_summary.csv"
)

OUT = (
    BASE
    / "24_r2_hydrogen_symmetry_refinement_preflight"
)

PAIR_CSV = (
    OUT
    / "r2_hydrogen_lower_upper_pair_candidates.csv"
)

HEAVY_PAIR_CSV = (
    OUT
    / "r2_heavy_lower_upper_pair_candidates.csv"
)

TRANSFORM_CSV = (
    OUT
    / "r2_lower_to_upper_symmetry_transform.csv"
)

SUMMARY_CSV = (
    OUT
    / "r2_hydrogen_symmetry_preflight_summary.csv"
)

JSON_OUT = (
    OUT
    / "r2_hydrogen_symmetry_preflight.json"
)

REPORT = (
    OUT
    / "R2_HYDROGEN_SYMMETRY_REFINEMENT_PREFLIGHT_DAY024.md"
)


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


def canonical_identifier(node_id: str) -> str:
    value = node_id

    replacements = (
        (":LOWER:", ":END:"),
        (":UPPER:", ":END:"),
        ("LOWER:", "END:"),
        ("UPPER:", "END:"),
        (":LOWER", ":END"),
        (":UPPER", ":END"),
        ("LOWER", "END"),
        ("UPPER", "END"),
    )

    for old, new in replacements:
        value = value.replace(old, new)

    return value


def metadata_signature(
    row: dict[str, str],
) -> tuple[str, ...]:
    ignored = {
        "node_id",
        "end",
        "x_nm",
        "y_nm",
        "z_nm",
    }

    preferred = (
        "element",
        "node_type",
        "circumferential_index",
        "angle_turns",
        "ring_index",
        "radial_index",
        "bridge_index",
        "bridge_position",
        "attachment_index",
        "parent_node",
        "attached_to",
    )

    values = []

    for key in preferred:
        if key in row and row[key] != "":
            values.append(
                f"{key}={row[key]}"
            )

    for key in sorted(row):
        if (
            key not in ignored
            and key not in preferred
            and row[key] != ""
        ):
            values.append(
                f"{key}={row[key]}"
            )

    return tuple(values)


def kabsch(
    source: np.ndarray,
    target: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
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

    if np.linalg.det(rotation) < 0.0:
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

    rmsd = float(
        np.sqrt(
            np.mean(
                np.sum(
                    (
                        transformed
                        - target
                    )
                    ** 2,
                    axis=1,
                )
            )
        )
    )

    return (
        rotation,
        translation,
        rmsd,
    )


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        GRAPH_NODES,
        GRAPH_EDGES,
        HEAVY_COORDINATES,
        FULL_COORDINATES,
        H_ORIENTATIONS,
        H_SUMMARY,
    ):
        require_file(required)

    node_rows = read_rows(
        GRAPH_NODES
    )

    edge_rows = read_rows(
        GRAPH_EDGES
    )

    heavy_rows = read_rows(
        HEAVY_COORDINATES
    )

    full_rows = read_rows(
        FULL_COORDINATES
    )

    orientation_rows = read_rows(
        H_ORIENTATIONS
    )

    h_summary = read_one(
        H_SUMMARY
    )

    nodes = {
        row["node_id"]: row
        for row in node_rows
    }

    heavy_coordinates = coordinates(
        heavy_rows
    )

    full_coordinates = coordinates(
        full_rows
    )

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

    adjacency = {
        node_id: set()
        for node_id in nodes
    }

    for row in edge_rows:
        first = row["source_node"]
        second = row["target_node"]

        adjacency[first].add(second)
        adjacency[second].add(first)

    attached_heavy = {}

    for row in lower_H + upper_H:
        node_id = row["node_id"]

        heavy_neighbors = [
            neighbor
            for neighbor in adjacency[node_id]
            if nodes[neighbor]["element"] != "H"
        ]

        if len(heavy_neighbors) != 1:
            raise RuntimeError(
                f"{node_id}: expected one attached heavy atom."
            )

        attached_heavy[node_id] = heavy_neighbors[0]

    upper_heavy_by_canonical = defaultdict(list)

    for row in upper_heavy:
        upper_heavy_by_canonical[
            canonical_identifier(
                row["node_id"]
            )
        ].append(row)

    heavy_pair_rows = []
    heavy_pairs = []

    for lower_row in lower_heavy:
        lower_id = lower_row["node_id"]

        candidates = upper_heavy_by_canonical[
            canonical_identifier(
                lower_id
            )
        ]

        if len(candidates) == 1:
            upper_row = candidates[0]
            method = "CANONICAL_NODE_ID"

        else:
            signature = metadata_signature(
                lower_row
            )

            candidates = [
                row
                for row in upper_heavy
                if (
                    row["element"]
                    == lower_row["element"]
                    and row["node_type"]
                    == lower_row["node_type"]
                    and metadata_signature(row)
                    == signature
                )
            ]

            if len(candidates) == 1:
                upper_row = candidates[0]
                method = "METADATA_SIGNATURE"

            else:
                heavy_pair_rows.append(
                    {
                        "lower_node": lower_id,
                        "upper_node": "",
                        "pair_method": "UNRESOLVED",
                        "candidate_count": len(candidates),
                        "node_type": lower_row["node_type"],
                        "element": lower_row["element"],
                    }
                )

                continue

        upper_id = upper_row["node_id"]

        heavy_pairs.append(
            (
                lower_id,
                upper_id,
            )
        )

        heavy_pair_rows.append(
            {
                "lower_node": lower_id,
                "upper_node": upper_id,
                "pair_method": method,
                "candidate_count": 1,
                "node_type": lower_row["node_type"],
                "element": lower_row["element"],
            }
        )

    write_rows(
        HEAVY_PAIR_CSV,
        heavy_pair_rows,
    )

    if len(heavy_pairs) < 3:
        raise RuntimeError(
            "Insufficient heavy pairs for rigid-transform audit."
        )

    source = np.asarray(
        [
            heavy_coordinates[lower_id]
            for lower_id, _
            in heavy_pairs
        ],
        dtype=float,
    )

    target = np.asarray(
        [
            heavy_coordinates[upper_id]
            for _, upper_id
            in heavy_pairs
        ],
        dtype=float,
    )

    rotation, translation, heavy_rmsd = kabsch(
        source,
        target,
    )

    transformed_source = (
        source
        @ rotation.T
        + translation
    )

    heavy_pair_deviations = np.linalg.norm(
        transformed_source
        - target,
        axis=1,
    )

    for row, deviation in zip(
        [
            row
            for row in heavy_pair_rows
            if row["upper_node"] != ""
        ],
        heavy_pair_deviations,
    ):
        row["symmetry_deviation_nm"] = float(
            deviation
        )

    write_rows(
        HEAVY_PAIR_CSV,
        heavy_pair_rows,
    )

    upper_H_by_role_and_parent = defaultdict(list)

    inverse_heavy_pair = {
        lower_id: upper_id
        for lower_id, upper_id
        in heavy_pairs
    }

    for row in upper_H:
        upper_H_by_role_and_parent[
            (
                row["node_type"],
                attached_heavy[
                    row["node_id"]
                ],
            )
        ].append(row)

    H_pair_rows = []

    for lower_row in lower_H:
        lower_id = lower_row["node_id"]
        lower_parent = attached_heavy[
            lower_id
        ]

        expected_upper_parent = inverse_heavy_pair.get(
            lower_parent
        )

        candidates = upper_H_by_role_and_parent[
            (
                lower_row["node_type"],
                expected_upper_parent,
            )
        ]

        if len(candidates) == 1:
            upper_row = candidates[0]
            upper_id = upper_row["node_id"]

            transformed_H = (
                rotation
                @ full_coordinates[
                    lower_id
                ]
                + translation
            )

            deviation = float(
                np.linalg.norm(
                    transformed_H
                    - full_coordinates[
                        upper_id
                    ]
                )
            )

            H_pair_rows.append(
                {
                    "lower_H": lower_id,
                    "upper_H": upper_id,
                    "hydrogen_role": lower_row["node_type"],
                    "lower_attached_heavy": lower_parent,
                    "upper_attached_heavy": expected_upper_parent,
                    "pair_status": "UNIQUE",
                    "symmetry_deviation_nm": deviation,
                    "lower_x_nm": full_coordinates[lower_id][0],
                    "lower_y_nm": full_coordinates[lower_id][1],
                    "lower_z_nm": full_coordinates[lower_id][2],
                    "upper_x_nm": full_coordinates[upper_id][0],
                    "upper_y_nm": full_coordinates[upper_id][1],
                    "upper_z_nm": full_coordinates[upper_id][2],
                    "transformed_lower_x_nm": transformed_H[0],
                    "transformed_lower_y_nm": transformed_H[1],
                    "transformed_lower_z_nm": transformed_H[2],
                }
            )

        else:
            H_pair_rows.append(
                {
                    "lower_H": lower_id,
                    "upper_H": "",
                    "hydrogen_role": lower_row["node_type"],
                    "lower_attached_heavy": lower_parent,
                    "upper_attached_heavy": (
                        expected_upper_parent
                        or ""
                    ),
                    "pair_status": "UNRESOLVED",
                    "candidate_count": len(candidates),
                }
            )

    write_rows(
        PAIR_CSV,
        H_pair_rows,
    )

    resolved_H = [
        row
        for row in H_pair_rows
        if row["pair_status"] == "UNIQUE"
    ]

    inner_H = [
        row
        for row in resolved_H
        if row["hydrogen_role"]
        == "ANNULUS_INNER_PASSIVANT_H"
    ]

    H_deviations = np.asarray(
        [
            float(
                row[
                    "symmetry_deviation_nm"
                ]
            )
            for row in resolved_H
        ],
        dtype=float,
    )

    inner_H_deviations = np.asarray(
        [
            float(
                row[
                    "symmetry_deviation_nm"
                ]
            )
            for row in inner_H
        ],
        dtype=float,
    )

    transform_rows = []

    for row_index in range(3):
        transform_rows.append(
            {
                "row": row_index,
                "R0": float(
                    rotation[row_index, 0]
                ),
                "R1": float(
                    rotation[row_index, 1]
                ),
                "R2": float(
                    rotation[row_index, 2]
                ),
                "translation_nm": float(
                    translation[row_index]
                ),
            }
        )

    write_rows(
        TRANSFORM_CSV,
        transform_rows,
    )

    unresolved_heavy = sum(
        row["upper_node"] == ""
        for row in heavy_pair_rows
    )

    unresolved_H = sum(
        row["pair_status"] != "UNIQUE"
        for row in H_pair_rows
    )

    summary = {
        "Gate3P_decision": h_summary.get(
            "decision",
            "",
        ),
        "lower_heavy_nodes": len(
            lower_heavy
        ),
        "upper_heavy_nodes": len(
            upper_heavy
        ),
        "resolved_heavy_pairs": len(
            heavy_pairs
        ),
        "unresolved_heavy_pairs": (
            unresolved_heavy
        ),
        "heavy_symmetry_RMSD_nm": (
            heavy_rmsd
        ),
        "heavy_symmetry_maximum_deviation_nm": float(
            np.max(
                heavy_pair_deviations
            )
        ),
        "lower_H_nodes": len(
            lower_H
        ),
        "upper_H_nodes": len(
            upper_H
        ),
        "resolved_H_pairs": len(
            resolved_H
        ),
        "unresolved_H_pairs": (
            unresolved_H
        ),
        "resolved_inner_H_pairs": len(
            inner_H
        ),
        "all_H_symmetry_RMSD_nm": (
            float(
                np.sqrt(
                    np.mean(
                        H_deviations ** 2
                    )
                )
            )
            if H_deviations.size
            else ""
        ),
        "all_H_symmetry_maximum_deviation_nm": (
            float(
                np.max(
                    H_deviations
                )
            )
            if H_deviations.size
            else ""
        ),
        "inner_H_symmetry_RMSD_nm": (
            float(
                np.sqrt(
                    np.mean(
                        inner_H_deviations ** 2
                    )
                )
            )
            if inner_H_deviations.size
            else ""
        ),
        "inner_H_symmetry_maximum_deviation_nm": (
            float(
                np.max(
                    inner_H_deviations
                )
            )
            if inner_H_deviations.size
            else ""
        ),
        "symmetry_constrained_refinement_ready": (
            unresolved_heavy == 0
            and unresolved_H == 0
            and len(resolved_H) == 102
            and len(inner_H) == 12
        ),
        "coordinates_modified": False,
        "topology_generated": False,
        "energy_minimized": False,
        "MD_performed": False,
    }

    write_rows(
        SUMMARY_CSV,
        [summary],
    )

    JSON_OUT.write_text(
        json.dumps(
            {
                "summary": summary,
                "rotation": rotation.tolist(),
                "translation_nm": translation.tolist(),
                "limitations": [
                    (
                        "This is a correspondence and symmetry audit only."
                    ),
                    (
                        "No hydrogen or heavy coordinates are modified."
                    ),
                    (
                        "No topology, charges, force-field parameters, "
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

    REPORT.write_text(
        f"""# R2 Hydrogen Symmetry Refinement Preflight

## Heavy correspondence

- Lower/upper heavy nodes:
  **{len(lower_heavy)}/{len(upper_heavy)}**
- Resolved/unresolved heavy pairs:
  **{len(heavy_pairs)}/{unresolved_heavy}**
- Heavy symmetry RMSD:
  **{heavy_rmsd:.12e} nm**
- Heavy maximum symmetry deviation:
  **{float(np.max(heavy_pair_deviations)):.12e} nm**

## Hydrogen correspondence

- Lower/upper H:
  **{len(lower_H)}/{len(upper_H)}**
- Resolved/unresolved H pairs:
  **{len(resolved_H)}/{unresolved_H}**
- Resolved inner-H pairs:
  **{len(inner_H)}**

## Existing H asymmetry

- All-H symmetry RMSD:
  **{summary['all_H_symmetry_RMSD_nm']} nm**
- All-H maximum deviation:
  **{summary['all_H_symmetry_maximum_deviation_nm']} nm**
- Inner-H symmetry RMSD:
  **{summary['inner_H_symmetry_RMSD_nm']} nm**
- Inner-H maximum deviation:
  **{summary['inner_H_symmetry_maximum_deviation_nm']} nm**

## Decision

- Symmetry-constrained refinement ready:
  **{summary['symmetry_constrained_refinement_ready']}**
- Coordinates modified:
  **NO**
- Molecular topology generated:
  **NO**
- Energy minimization performed:
  **NO**
- MD performed:
  **NO**
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 hydrogen symmetry correspondence "
        "preflight completed."
    )

    print(
        "Heavy lower/upper/resolved/unresolved: "
        f"{len(lower_heavy)}/"
        f"{len(upper_heavy)}/"
        f"{len(heavy_pairs)}/"
        f"{unresolved_heavy}"
    )

    print(
        "Heavy symmetry RMSD / maximum deviation: "
        f"{heavy_rmsd:.12e}/"
        f"{float(np.max(heavy_pair_deviations)):.12e} nm"
    )

    print(
        "H lower/upper/resolved/unresolved: "
        f"{len(lower_H)}/"
        f"{len(upper_H)}/"
        f"{len(resolved_H)}/"
        f"{unresolved_H}"
    )

    print(
        "Resolved inner-H pairs: "
        f"{len(inner_H)}"
    )

    print(
        "All-H symmetry RMSD / maximum: "
        f"{summary['all_H_symmetry_RMSD_nm']}/"
        f"{summary['all_H_symmetry_maximum_deviation_nm']} nm"
    )

    print(
        "Inner-H symmetry RMSD / maximum: "
        f"{summary['inner_H_symmetry_RMSD_nm']}/"
        f"{summary['inner_H_symmetry_maximum_deviation_nm']} nm"
    )

    print(
        "Symmetry-constrained refinement ready: "
        f"{summary['symmetry_constrained_refinement_ready']}"
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

    for path in (
        HEAVY_PAIR_CSV,
        PAIR_CSV,
        TRANSFORM_CSV,
        SUMMARY_CSV,
        JSON_OUT,
        REPORT,
    ):
        print(
            f"Wrote: {path.relative_to(ROOT)}"
        )


if __name__ == "__main__":
    main()
