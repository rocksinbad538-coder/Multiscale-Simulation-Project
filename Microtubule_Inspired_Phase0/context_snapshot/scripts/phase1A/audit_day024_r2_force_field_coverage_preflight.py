#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
from collections import Counter, defaultdict
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

GATE3Q = (
    BASE
    / "29_r2_chemical_realizability_and_parameterization_scope"
)

GRAPH_NODES = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_nodes.csv"
)

GRAPH_EDGES = (
    GATE3M
    / "r2_selected_longer_bn_bridge_graph_edges.csv"
)

COORDINATES = (
    GATE3P2
    / "r2_selected_four_atom_refined_full_coordinates.csv"
)

CHEMICAL_SUMMARY = (
    GATE3Q
    / "r2_chemical_realizability_and_parameterization_scope_summary.csv"
)

NODE_ENVIRONMENTS = (
    GATE3Q
    / "r2_node_chemical_environment_assignments.csv"
)

LOCAL_ENVIRONMENTS = (
    GATE3Q
    / "r2_local_chemical_environment_inventory.csv"
)

PARAMETERIZATION_SCOPE = (
    GATE3Q
    / "r2_parameterization_scope.csv"
)

CRITICAL_CENTERS = (
    GATE3Q
    / "r2_parameterization_critical_centers.csv"
)

OUT = (
    BASE
    / "30_r2_force_field_coverage_preflight"
)

BOND_REQUIREMENTS = (
    OUT
    / "r2_required_bond_terms.csv"
)

ANGLE_REQUIREMENTS = (
    OUT
    / "r2_required_angle_terms.csv"
)

TORSION_REQUIREMENTS = (
    OUT
    / "r2_required_torsion_terms.csv"
)

IMPROPER_REQUIREMENTS = (
    OUT
    / "r2_required_improper_centers.csv"
)

LOCAL_FF_FILES = (
    OUT
    / "r2_local_force_field_file_inventory.csv"
)

LOCAL_FF_HITS = (
    OUT
    / "r2_local_force_field_bn_h_text_hits.csv"
)

QM_FRAGMENT_CLASSES = (
    OUT
    / "r2_preliminary_qm_fragment_classes.csv"
)

SUMMARY = (
    OUT
    / "r2_force_field_coverage_preflight_summary.csv"
)

JSON_OUT = (
    OUT
    / "r2_force_field_coverage_preflight.json"
)

MANIFEST = (
    OUT
    / "r2_force_field_coverage_preflight_manifest.csv"
)

REPORT = (
    OUT
    / "R2_FORCE_FIELD_COVERAGE_PREFLIGHT_DAY024.md"
)

EXPECTED_CHEMICAL_DECISION = (
    "R2_STATIC_CHEMICAL_REALIZABILITY_VALIDATED_"
    "PARAMETERIZATION_SCOPE_DEFINED"
)

BOND_MAX_FOR_TOPOLOGY_NM = 0.25

FF_SUFFIXES = {
    ".itp",
    ".top",
    ".rtp",
    ".atp",
    ".ff",
    ".prm",
    ".par",
    ".str",
    ".lib",
    ".frc",
    ".xml",
    ".off",
    ".dat",
    ".txt",
}

EXCLUDED_DIR_NAMES = {
    ".git",
    "__pycache__",
    "runs",
}

SEARCH_PATTERNS = (
    r"\bBN\b",
    r"\bB-N\b",
    r"\bB_H\b",
    r"\bB-H\b",
    r"\bN_H\b",
    r"\bN-H\b",
    r"boron",
    r"nitride",
    r"h-?BN",
    r"hexagonal boron",
)


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty required file: {path}"
        )


def relative(path: Path) -> str:
    try:
        return str(
            path.resolve().relative_to(ROOT)
        )
    except ValueError:
        return str(
            path.resolve()
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
        rows = [
            {
                "status": "NO_ROWS",
            }
        ]

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
    value = float(
        row[key]
    )

    if not math.isfinite(value):
        raise RuntimeError(
            f"Non-finite numeric value in {key!r}"
        )

    return value


def canonical_bond(
    first: str,
    second: str,
) -> tuple[str, str]:
    return tuple(
        sorted(
            (
                first,
                second,
            )
        )
    )


def canonical_angle(
    first: str,
    center: str,
    third: str,
) -> tuple[str, str, str]:
    forward = (
        first,
        center,
        third,
    )

    reverse = (
        third,
        center,
        first,
    )

    return min(
        forward,
        reverse,
    )


def canonical_torsion(
    first: str,
    second: str,
    third: str,
    fourth: str,
) -> tuple[str, str, str, str]:
    forward = (
        first,
        second,
        third,
        fourth,
    )

    reverse = (
        fourth,
        third,
        second,
        first,
    )

    return min(
        forward,
        reverse,
    )


def angle_deg(
    first: np.ndarray,
    center: np.ndarray,
    third: np.ndarray,
) -> float:
    v1 = first - center
    v2 = third - center

    n1 = float(
        np.linalg.norm(v1)
    )

    n2 = float(
        np.linalg.norm(v2)
    )

    if n1 <= 1.0e-12 or n2 <= 1.0e-12:
        raise RuntimeError(
            "Zero-length vector in angle calculation."
        )

    cosine = float(
        np.clip(
            np.dot(v1, v2) / (n1 * n2),
            -1.0,
            1.0,
        )
    )

    return math.degrees(
        math.acos(cosine)
    )


def dihedral_deg(
    p0: np.ndarray,
    p1: np.ndarray,
    p2: np.ndarray,
    p3: np.ndarray,
) -> float:
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2

    b1_norm = float(
        np.linalg.norm(b1)
    )

    if b1_norm <= 1.0e-12:
        raise RuntimeError(
            "Zero central bond in torsion calculation."
        )

    b1_unit = b1 / b1_norm

    v = (
        b0
        - np.dot(b0, b1_unit)
        * b1_unit
    )

    w = (
        b2
        - np.dot(b2, b1_unit)
        * b1_unit
    )

    x_value = float(
        np.dot(v, w)
    )

    y_value = float(
        np.dot(
            np.cross(
                b1_unit,
                v,
            ),
            w,
        )
    )

    return math.degrees(
        math.atan2(
            y_value,
            x_value,
        )
    )


def node_label(
    row: dict[str, str],
) -> str:
    return (
        f"{row['element']}|"
        f"{row['node_type']}"
    )


def iter_candidate_ff_files() -> list[Path]:
    roots = [
        ROOT,
        ROOT / "forcefields",
        ROOT / "forcefield",
        ROOT / "ff",
        Path.home() / ".local/share/gromacs/top",
        Path("/usr/local/share/gromacs/top"),
        Path("/opt/homebrew/share/gromacs/top"),
    ]

    found = set()

    for search_root in roots:
        if not search_root.exists():
            continue

        if search_root == ROOT:
            for current_root, dir_names, file_names in os.walk(
                search_root
            ):
                dir_names[:] = [
                    name
                    for name in dir_names
                    if name not in EXCLUDED_DIR_NAMES
                ]

                current = Path(
                    current_root
                )

                for name in file_names:
                    path = current / name

                    if (
                        path.suffix.lower()
                        in FF_SUFFIXES
                    ):
                        found.add(
                            path.resolve()
                        )

        else:
            for path in search_root.rglob("*"):
                if (
                    path.is_file()
                    and path.suffix.lower()
                    in FF_SUFFIXES
                ):
                    found.add(
                        path.resolve()
                    )

    return sorted(
        found
    )


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        GRAPH_NODES,
        GRAPH_EDGES,
        COORDINATES,
        CHEMICAL_SUMMARY,
        NODE_ENVIRONMENTS,
        LOCAL_ENVIRONMENTS,
        PARAMETERIZATION_SCOPE,
        CRITICAL_CENTERS,
    ):
        require_file(required)

    chemical_summary = read_one(
        CHEMICAL_SUMMARY
    )

    if chemical_summary.get(
        "decision"
    ) != EXPECTED_CHEMICAL_DECISION:
        raise RuntimeError(
            "Gate 3Q chemical audit is not accepted."
        )

    node_rows = read_rows(
        GRAPH_NODES
    )

    edge_rows = read_rows(
        GRAPH_EDGES
    )

    coordinate_rows = read_rows(
        COORDINATES
    )

    environment_rows = read_rows(
        NODE_ENVIRONMENTS
    )

    local_environment_rows = read_rows(
        LOCAL_ENVIRONMENTS
    )

    scope_rows = read_rows(
        PARAMETERIZATION_SCOPE
    )

    critical_rows = read_rows(
        CRITICAL_CENTERS
    )

    nodes = {
        row["node_id"]: row
        for row in node_rows
    }

    positions = {
        row["node_id"]: np.asarray(
            [
                parse_float(row, "x_nm"),
                parse_float(row, "y_nm"),
                parse_float(row, "z_nm"),
            ],
            dtype=float,
        )
        for row in coordinate_rows
    }

    environment_by_node = {
        row["node_id"]: row
        for row in environment_rows
    }

    adjacency = {
        node_id: set()
        for node_id in nodes
    }

    edge_set = set()

    for row in edge_rows:
        first = row["source_node"]
        second = row["target_node"]

        key = tuple(
            sorted(
                (
                    first,
                    second,
                )
            )
        )

        edge_set.add(
            key
        )

        adjacency[first].add(
            second
        )

        adjacency[second].add(
            first
        )

    bond_counter = Counter()
    bond_examples = {}
    bond_lengths = defaultdict(list)

    for first, second in sorted(
        edge_set
    ):
        first_label = node_label(
            nodes[first]
        )

        second_label = node_label(
            nodes[second]
        )

        key = canonical_bond(
            first_label,
            second_label,
        )

        bond_counter[key] += 1

        bond_examples.setdefault(
            key,
            (
                first,
                second,
            ),
        )

        bond_lengths[key].append(
            float(
                np.linalg.norm(
                    positions[first]
                    - positions[second]
                )
            )
        )

    bond_requirement_rows = []

    for index, key in enumerate(
        sorted(
            bond_counter
        ),
        start=1,
    ):
        values = np.asarray(
            bond_lengths[key],
            dtype=float,
        )

        example = bond_examples[
            key
        ]

        bond_requirement_rows.append(
            {
                "bond_term_id": (
                    f"BOND_{index:04d}"
                ),
                "atom_type_1": key[0],
                "atom_type_2": key[1],
                "count": bond_counter[key],
                "minimum_nm": float(
                    np.min(values)
                ),
                "mean_nm": float(
                    np.mean(values)
                ),
                "maximum_nm": float(
                    np.max(values)
                ),
                "example_node_1": example[0],
                "example_node_2": example[1],
                "coverage_status": "UNASSESSED",
            }
        )

    write_rows(
        BOND_REQUIREMENTS,
        bond_requirement_rows,
    )

    angle_counter = Counter()
    angle_examples = {}
    angle_values = defaultdict(list)

    for center in sorted(
        nodes
    ):
        neighbors = sorted(
            adjacency[center]
        )

        for first_index in range(
            len(neighbors)
        ):
            for third_index in range(
                first_index + 1,
                len(neighbors),
            ):
                first = neighbors[
                    first_index
                ]

                third = neighbors[
                    third_index
                ]

                key = canonical_angle(
                    node_label(
                        nodes[first]
                    ),
                    node_label(
                        nodes[center]
                    ),
                    node_label(
                        nodes[third]
                    ),
                )

                angle_counter[key] += 1

                angle_examples.setdefault(
                    key,
                    (
                        first,
                        center,
                        third,
                    ),
                )

                angle_values[key].append(
                    angle_deg(
                        positions[first],
                        positions[center],
                        positions[third],
                    )
                )

    angle_requirement_rows = []

    for index, key in enumerate(
        sorted(
            angle_counter
        ),
        start=1,
    ):
        values = np.asarray(
            angle_values[key],
            dtype=float,
        )

        example = angle_examples[
            key
        ]

        angle_requirement_rows.append(
            {
                "angle_term_id": (
                    f"ANGLE_{index:04d}"
                ),
                "atom_type_1": key[0],
                "center_type": key[1],
                "atom_type_3": key[2],
                "count": angle_counter[key],
                "minimum_deg": float(
                    np.min(values)
                ),
                "mean_deg": float(
                    np.mean(values)
                ),
                "maximum_deg": float(
                    np.max(values)
                ),
                "example_node_1": example[0],
                "example_center": example[1],
                "example_node_3": example[2],
                "coverage_status": "UNASSESSED",
            }
        )

    write_rows(
        ANGLE_REQUIREMENTS,
        angle_requirement_rows,
    )

    torsion_counter = Counter()
    torsion_examples = {}
    torsion_values = defaultdict(list)

    for second, third in sorted(
        edge_set
    ):
        second_neighbors = sorted(
            adjacency[second]
            - {
                third
            }
        )

        third_neighbors = sorted(
            adjacency[third]
            - {
                second
            }
        )

        for first in second_neighbors:
            for fourth in third_neighbors:
                if first == fourth:
                    continue

                key = canonical_torsion(
                    node_label(
                        nodes[first]
                    ),
                    node_label(
                        nodes[second]
                    ),
                    node_label(
                        nodes[third]
                    ),
                    node_label(
                        nodes[fourth]
                    ),
                )

                torsion_counter[key] += 1

                torsion_examples.setdefault(
                    key,
                    (
                        first,
                        second,
                        third,
                        fourth,
                    ),
                )

                torsion_values[key].append(
                    dihedral_deg(
                        positions[first],
                        positions[second],
                        positions[third],
                        positions[fourth],
                    )
                )

    torsion_requirement_rows = []

    for index, key in enumerate(
        sorted(
            torsion_counter
        ),
        start=1,
    ):
        values = np.asarray(
            torsion_values[key],
            dtype=float,
        )

        example = torsion_examples[
            key
        ]

        torsion_requirement_rows.append(
            {
                "torsion_term_id": (
                    f"TORSION_{index:04d}"
                ),
                "atom_type_1": key[0],
                "atom_type_2": key[1],
                "atom_type_3": key[2],
                "atom_type_4": key[3],
                "count": torsion_counter[key],
                "minimum_deg": float(
                    np.min(values)
                ),
                "mean_deg": float(
                    np.mean(values)
                ),
                "maximum_deg": float(
                    np.max(values)
                ),
                "example_node_1": example[0],
                "example_node_2": example[1],
                "example_node_3": example[2],
                "example_node_4": example[3],
                "coverage_status": "UNASSESSED",
            }
        )

    write_rows(
        TORSION_REQUIREMENTS,
        torsion_requirement_rows,
    )

    improper_rows = []

    for node_id in sorted(
        nodes
    ):
        if nodes[
            node_id
        ][
            "element"
        ] == "H":
            continue

        neighbors = sorted(
            adjacency[
                node_id
            ]
        )

        if len(neighbors) != 3:
            continue

        neighbor_elements = Counter(
            nodes[neighbor][
                "element"
            ]
            for neighbor in neighbors
        )

        improper_rows.append(
            {
                "center_node": node_id,
                "center_element": nodes[
                    node_id
                ][
                    "element"
                ],
                "center_type": nodes[
                    node_id
                ][
                    "node_type"
                ],
                "region": environment_by_node[
                    node_id
                ][
                    "region"
                ],
                "neighbor_1": neighbors[0],
                "neighbor_2": neighbors[1],
                "neighbor_3": neighbors[2],
                "neighbor_elements": (
                    " | ".join(
                        f"{key}:{value}"
                        for key, value
                        in sorted(
                            neighbor_elements.items()
                        )
                    )
                ),
                "improper_or_planarity_term_status": (
                    "UNASSESSED"
                ),
            }
        )

    write_rows(
        IMPROPER_REQUIREMENTS,
        improper_rows,
    )

    candidate_ff_files = iter_candidate_ff_files()

    ff_file_rows = []
    hit_rows = []

    combined_pattern = re.compile(
        "|".join(
            SEARCH_PATTERNS
        ),
        re.IGNORECASE,
    )

    for path in candidate_ff_files:
        size_bytes = path.stat().st_size

        ff_file_rows.append(
            {
                "file": relative(
                    path
                ),
                "suffix": path.suffix.lower(),
                "size_bytes": size_bytes,
                "sha256": (
                    sha256(path)
                    if size_bytes
                    <= 50 * 1024 * 1024
                    else "SKIPPED_GT_50MB"
                ),
            }
        )

        if size_bytes > 20 * 1024 * 1024:
            continue

        try:
            lines = path.read_text(
                encoding="utf-8",
                errors="ignore",
            ).splitlines()
        except OSError:
            continue

        for line_number, line in enumerate(
            lines,
            start=1,
        ):
            if combined_pattern.search(
                line
            ):
                hit_rows.append(
                    {
                        "file": relative(
                            path
                        ),
                        "line_number": line_number,
                        "text": line.strip()[
                            :500
                        ],
                    }
                )

    write_rows(
        LOCAL_FF_FILES,
        ff_file_rows,
    )

    write_rows(
        LOCAL_FF_HITS,
        hit_rows,
    )

    critical_region_counts = Counter(
        row["region"]
        for row in critical_rows
    )

    fragment_rows = [
        {
            "fragment_class": "QM_F01_PARENT_BULK_HBN",
            "priority": "MEDIUM",
            "target_environment": (
                "PARENT_HBN_BULK_LIKE"
            ),
            "purpose": (
                "Reference B-N equilibrium geometry, curvature response, "
                "charges and nonbonded behavior for the parent scaffold."
            ),
            "minimum_content": (
                "Hydrogen-terminated h-BN patch or short nanotube segment."
            ),
            "geometry_source": (
                "Extract representative parent-scaffold environment."
            ),
            "calculation_authorized": False,
        },
        {
            "fragment_class": "QM_F02_PARENT_SEED_JUNCTION",
            "priority": "HIGH",
            "target_environment": (
                "PARENT_SEED_ATTACHMENT"
            ),
            "purpose": (
                "Validate reconstructed junction bonds, angles, "
                "charge redistribution and planarity."
            ),
            "minimum_content": (
                "Seed atom plus parent first and second coordination shells, "
                "hydrogen terminated."
            ),
            "geometry_source": (
                "Extract representative junction from each B/N sublattice."
            ),
            "calculation_authorized": False,
        },
        {
            "fragment_class": "QM_F03_ANNULUS_INNER_EDGE_BH",
            "priority": "HIGH",
            "target_environment": (
                "ANNULUS_INNER_BOUNDARY_WITH_B_H"
            ),
            "purpose": (
                "Validate inner-rim B-H bond, angles, partial charge "
                "and water-facing electrostatics."
            ),
            "minimum_content": (
                "Inner annulus B-H center with two coordination shells."
            ),
            "geometry_source": (
                "Extract representative lower/upper inner-rim environment."
            ),
            "calculation_authorized": False,
        },
        {
            "fragment_class": "QM_F04_ANNULUS_INNER_EDGE_NH",
            "priority": "HIGH",
            "target_environment": (
                "ANNULUS_INNER_BOUNDARY_WITH_N_H"
            ),
            "purpose": (
                "Validate inner-rim N-H bond, angles, partial charge "
                "and water-facing electrostatics."
            ),
            "minimum_content": (
                "Inner annulus N-H center with two coordination shells."
            ),
            "geometry_source": (
                "Extract representative lower/upper inner-rim environment."
            ),
            "calculation_authorized": False,
        },
        {
            "fragment_class": "QM_F05_ANNULUS_OUTER_EDGE",
            "priority": "HIGH",
            "target_environment": (
                "ANNULUS_OUTER_BOUNDARY"
            ),
            "purpose": (
                "Validate outer-rim B-H/N-H environments and local "
                "edge stiffness."
            ),
            "minimum_content": (
                "Representative outer-edge B-H and N-H fragments."
            ),
            "geometry_source": (
                "Extract one representative of each elemental environment."
            ),
            "calculation_authorized": False,
        },
        {
            "fragment_class": "QM_F06_FOUR_ATOM_BRIDGE",
            "priority": "HIGHEST",
            "target_environment": (
                "ALTERNATING_BN_FOUR_ATOM_BRIDGE"
            ),
            "purpose": (
                "Validate bridge bond, angle, torsional and charge terms."
            ),
            "minimum_content": (
                "Full B-N-B-N bridge plus both attachment centers and "
                "their first coordination shells."
            ),
            "geometry_source": (
                "Extract representative lower and upper bridge conformers."
            ),
            "calculation_authorized": False,
        },
        {
            "fragment_class": "QM_F07_BRIDGE_ATTACHMENT",
            "priority": "HIGHEST",
            "target_environment": (
                "ANNULUS_BRIDGE_ATTACHMENT_AND_SEED_BRIDGE_ATTACHMENT"
            ),
            "purpose": (
                "Resolve coupling between bridge parameters and attachment "
                "center parameters."
            ),
            "minimum_content": (
                "Bridge endpoint, attached rim/seed atom and two shells."
            ),
            "geometry_source": (
                "Extract both chemically distinct bridge endpoints."
            ),
            "calculation_authorized": False,
        },
    ]

    write_rows(
        QM_FRAGMENT_CLASSES,
        fragment_rows,
    )

    summary = {
        "decision": (
            "R2_FORCE_FIELD_COVERAGE_PREFLIGHT_COMPLETED"
        ),
        "required_bond_term_classes": len(
            bond_requirement_rows
        ),
        "required_angle_term_classes": len(
            angle_requirement_rows
        ),
        "required_torsion_term_classes": len(
            torsion_requirement_rows
        ),
        "three_coordinate_heavy_improper_centers": len(
            improper_rows
        ),
        "local_force_field_files_found": len(
            candidate_ff_files
        ),
        "local_BN_H_text_hits": len(
            hit_rows
        ),
        "unique_local_chemical_environments": len(
            local_environment_rows
        ),
        "parameterization_critical_centers": len(
            critical_rows
        ),
        "preliminary_QM_fragment_classes": len(
            fragment_rows
        ),
        "force_field_coverage_established": False,
        "topology_generation_authorized": False,
        "formal_charge_assignment_authorized": False,
        "force_field_parameterization_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "QM_calculation_authorized": False,
        "literature_and_force_field_source_review_required": True,
        "required_next_step": (
            "REVIEW_LOCAL_FORCE_FIELD_HITS_AND_PERFORM_"
            "PRIMARY_SOURCE_FORCE_FIELD_LITERATURE_AUDIT"
        ),
    }

    write_rows(
        SUMMARY,
        [
            summary
        ],
    )

    JSON_OUT.write_text(
        json.dumps(
            {
                "summary": summary,
                "critical_region_counts": dict(
                    critical_region_counts
                ),
                "parameterization_scope": scope_rows,
                "limitations": [
                    (
                        "Text hits do not establish scientifically valid "
                        "parameter coverage."
                    ),
                    (
                        "Atom-type naming compatibility is not assumed."
                    ),
                    (
                        "No topology, charges or parameter files are generated."
                    ),
                    (
                        "No minimization, MD or QM calculation is performed."
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
                "role": "Gate3P2_coordinates",
                "file": relative(
                    COORDINATES
                ),
                "sha256": sha256(
                    COORDINATES
                ),
            },
            {
                "role": "Gate3Q_summary",
                "file": relative(
                    CHEMICAL_SUMMARY
                ),
                "sha256": sha256(
                    CHEMICAL_SUMMARY
                ),
            },
            {
                "role": "Gate3Q_local_environments",
                "file": relative(
                    LOCAL_ENVIRONMENTS
                ),
                "sha256": sha256(
                    LOCAL_ENVIRONMENTS
                ),
            },
            {
                "role": "Gate3Q_critical_centers",
                "file": relative(
                    CRITICAL_CENTERS
                ),
                "sha256": sha256(
                    CRITICAL_CENTERS
                ),
            },
        ],
    )

    REPORT.write_text(
        f"""# R2 Force-Field Coverage Preflight

## Required bonded-term classes

- Bond classes: **{len(bond_requirement_rows)}**
- Angle classes: **{len(angle_requirement_rows)}**
- Torsion classes: **{len(torsion_requirement_rows)}**
- Three-coordinate heavy centers requiring planarity/improper review:
  **{len(improper_rows)}**

## Existing local assets

- Candidate force-field files found:
  **{len(candidate_ff_files)}**
- Text hits containing BN/B-H/N-H terminology:
  **{len(hit_rows)}**

Text matches are not treated as validated parameter coverage.

## Chemical scope

- Unique local environments:
  **{len(local_environment_rows)}**
- Parameterization-critical centers:
  **{len(critical_rows)}**

## Preliminary QM fragment set

- Fragment classes:
  **{len(fragment_rows)}**
- Highest-priority classes:
  four-atom bridge and bridge-attachment environments.

## Restrictions

- Force-field coverage established: **NO**
- Topology generation authorized: **NO**
- Charge assignment authorized: **NO**
- Parameterization authorized: **NO**
- Energy minimization authorized: **NO**
- MD authorized: **NO**
- QM calculation authorized: **NO**

## Required next step

`{summary['required_next_step']}`
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 force-field coverage preflight completed."
    )

    print(
        "Required bond / angle / torsion classes: "
        f"{len(bond_requirement_rows)}/"
        f"{len(angle_requirement_rows)}/"
        f"{len(torsion_requirement_rows)}"
    )

    print(
        "Three-coordinate heavy improper/planarity centers: "
        f"{len(improper_rows)}"
    )

    print(
        "Local force-field files / BN-H text hits: "
        f"{len(candidate_ff_files)}/"
        f"{len(hit_rows)}"
    )

    print(
        "Unique local environments / critical centers: "
        f"{len(local_environment_rows)}/"
        f"{len(critical_rows)}"
    )

    print(
        "Preliminary QM fragment classes: "
        f"{len(fragment_rows)}"
    )

    print(
        "Decision: "
        f"{summary['decision']}"
    )

    print(
        "Force-field coverage established: NO"
    )

    print(
        "Topology generation authorized: NO"
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
        "QM calculation authorized: NO"
    )

    print(
        "Required next step: "
        f"{summary['required_next_step']}"
    )

    for path in (
        BOND_REQUIREMENTS,
        ANGLE_REQUIREMENTS,
        TORSION_REQUIREMENTS,
        IMPROPER_REQUIREMENTS,
        LOCAL_FF_FILES,
        LOCAL_FF_HITS,
        QM_FRAGMENT_CLASSES,
        SUMMARY,
        JSON_OUT,
        MANIFEST,
        REPORT,
    ):
        print(
            f"Wrote: {relative(path)}"
        )


if __name__ == "__main__":
    main()
