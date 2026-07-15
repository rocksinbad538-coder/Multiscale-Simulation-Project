#!/usr/bin/env python3

from __future__ import annotations

import csv
import importlib.util
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

AUDIT_SCRIPT = (
    ROOT
    / "scripts/phase1A/"
    "audit_day024_r2_parent_rim_and_chemical_constraints.py"
)

SYSTEM_GRO = (
    ROOT
    / "runs/phase1A/day023_confinement_design/"
    "15_r2_frozen_solute_nvt_20ps_preparation/"
    "r2_frozen_solute_nvt_20ps_input.gro"
)

TPR_DUMP = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "01_r2_parent_rim_chemical_audit/"
    "r2_parent_system_tpr_dump.txt"
)

PREFERRED_TOP = (
    ROOT
    / "parameters/phase1A/accepted/"
    "hybrid_hbnBonded_kang2000_improperGeo100_validated/"
    "hbn_pyrene_4_hydratable_gap45_pyr5shift_clean032_"
    "hbnBonded_kang2000_improperGeo100.top"
)

OUTPUT_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "01_r2_parent_rim_chemical_audit/"
    "topology_terminal_classification"
)

TERMINAL_ATOMS_CSV = (
    OUTPUT_ROOT
    / "hbn_explicit_topology_terminal_atoms.csv"
)

END_SUMMARY_CSV = (
    OUTPUT_ROOT
    / "hbn_explicit_topology_terminal_end_summary.csv"
)

BOND_COMPARISON_CSV = (
    OUTPUT_ROOT
    / "hbn_explicit_vs_geometry_bond_comparison.csv"
)

TOPOLOGY_SOURCES_CSV = (
    OUTPUT_ROOT
    / "hbn_topology_include_sources.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "hbn_terminal_coordination_classification_summary.csv"
)

GATES_CSV = (
    OUTPUT_ROOT
    / "hbn_terminal_coordination_classification_gates.csv"
)

SUMMARY_JSON = (
    OUTPUT_ROOT
    / "hbn_terminal_coordination_classification.json"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "HBN_EXPLICIT_TOPOLOGY_AND_TERMINAL_CLASSIFICATION_DAY024.md"
)

EXPECTED_HBN_ATOMS = 1680
EXPECTED_B_ATOMS = 840
EXPECTED_N_ATOMS = 840
EXPECTED_BONDS = 2460

EXPECTED_TERMINAL_ATOMS = 60
EXPECTED_TERMINAL_ATOMS_PER_END = 30
EXPECTED_INTERIOR_ATOMS = 1620

BOND_MIN_NM = 0.115
BOND_MAX_NM = 0.175

VALID_BOND_MIN_NM = 0.125
VALID_BOND_MAX_NM = 0.165

MAX_END_RADIUS_DIFFERENCE_NM = 0.050

DECISION_PASS = (
    "HBN_EXPLICIT_TOPOLOGY_AND_TERMINAL_COORDINATION_CLASSIFIED"
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
                    field: row.get(
                        field,
                        "",
                    )
                    for field in fields
                }
            )


def load_audit_module():
    require_file(AUDIT_SCRIPT)

    spec = importlib.util.spec_from_file_location(
        "day024_parent_rim_audit",
        AUDIT_SCRIPT,
    )

    if (
        spec is None
        or spec.loader is None
    ):
        raise RuntimeError(
            "Could not load the parent-rim audit module."
        )

    module = importlib.util.module_from_spec(
        spec
    )

    spec.loader.exec_module(
        module
    )

    return module


def strip_comment(line: str) -> str:
    return line.split(
        ";",
        1,
    )[0].strip()


def resolve_include(
    parent: Path,
    include_name: str,
) -> Path | None:
    candidates = [
        parent.parent / include_name,
        ROOT / include_name,
        ROOT / "parameters" / include_name,
    ]

    for candidate in candidates:
        if (
            candidate.exists()
            and candidate.is_file()
        ):
            return candidate.resolve()

    return None


def collect_local_include_tree(
    top_path: Path,
) -> list[Path]:
    visited: set[Path] = set()
    ordered: list[Path] = []

    include_pattern = re.compile(
        r'^\s*#include\s+"([^"]+)"'
    )

    def visit(path: Path) -> None:
        resolved = path.resolve()

        if resolved in visited:
            return

        require_file(resolved)

        visited.add(resolved)
        ordered.append(resolved)

        text = resolved.read_text(
            encoding="utf-8",
            errors="replace",
        )

        for raw_line in text.splitlines():
            match = include_pattern.match(
                raw_line
            )

            if match is None:
                continue

            child = resolve_include(
                resolved,
                match.group(1),
            )

            if child is not None:
                visit(child)

    visit(top_path)

    return ordered


def parse_hbn_definition(
    files: list[Path],
) -> tuple[
    dict[int, dict[str, Any]],
    set[tuple[int, int]],
    list[dict[str, Any]],
]:
    atoms: dict[
        int,
        dict[str, Any]
    ] = {}

    bonds: set[
        tuple[int, int]
    ] = set()

    source_rows: list[
        dict[str, Any]
    ] = []

    connectivity_sections = {
        "bonds",
        "constraints",
    }

    for path in files:
        lines = path.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines()

        section = ""
        current_molecule = ""
        awaiting_molecule_name = False

        file_atom_records = 0
        file_connectivity_records = 0

        for line_number, raw_line in enumerate(
            lines,
            start=1,
        ):
            stripped = strip_comment(
                raw_line
            )

            if not stripped:
                continue

            if stripped.startswith("#"):
                continue

            section_match = re.match(
                r"^\[\s*([^\]]+?)\s*\]$",
                stripped,
            )

            if section_match is not None:
                section = (
                    section_match.group(1)
                    .strip()
                    .lower()
                )

                if section == "moleculetype":
                    current_molecule = ""
                    awaiting_molecule_name = True

                continue

            fields = stripped.split()

            if (
                section == "moleculetype"
                and awaiting_molecule_name
            ):
                current_molecule = fields[0]
                awaiting_molecule_name = False
                continue

            if current_molecule.upper() != "HBN":
                continue

            if section == "atoms":
                try:
                    atom_index_1based = int(
                        fields[0]
                    )
                except (
                    IndexError,
                    ValueError,
                ):
                    continue

                atoms[
                    atom_index_1based
                ] = {
                    "source_file": relative(
                        path
                    ),
                    "source_line": line_number,
                    "raw_record": stripped,
                }

                file_atom_records += 1

            elif section in connectivity_sections:
                if len(fields) < 2:
                    continue

                try:
                    first = int(
                        fields[0]
                    )
                    second = int(
                        fields[1]
                    )
                except ValueError:
                    continue

                if first == second:
                    raise RuntimeError(
                        "Self-connectivity detected in HBN topology: "
                        f"{relative(path)}:{line_number}"
                    )

                pair = tuple(
                    sorted(
                        (
                            first,
                            second,
                        )
                    )
                )

                bonds.add(pair)
                file_connectivity_records += 1

        if (
            file_atom_records > 0
            or file_connectivity_records > 0
        ):
            source_rows.append(
                {
                    "file": relative(
                        path
                    ),
                    "HBN_atom_records": (
                        file_atom_records
                    ),
                    "HBN_connectivity_records": (
                        file_connectivity_records
                    ),
                }
            )

    return atoms, bonds, source_rows


def find_accepted_topology() -> tuple[
    Path,
    list[Path],
    dict[int, dict[str, Any]],
    set[tuple[int, int]],
    list[dict[str, Any]],
]:
    candidates: list[Path] = []

    if PREFERRED_TOP.exists():
        candidates.append(
            PREFERRED_TOP
        )

    accepted_root = (
        ROOT
        / "parameters/phase1A/accepted"
    )

    if accepted_root.exists():
        for path in accepted_root.rglob(
            "*.top"
        ):
            if path not in candidates:
                candidates.append(
                    path
                )

    if not candidates:
        raise RuntimeError(
            "No accepted topology candidates were found."
        )

    evaluated = []

    for top_path in candidates:
        try:
            include_files = (
                collect_local_include_tree(
                    top_path
                )
            )

            (
                atoms,
                bonds,
                source_rows,
            ) = parse_hbn_definition(
                include_files
            )
        except Exception:
            continue

        score = 0

        if (
            top_path.resolve()
            == PREFERRED_TOP.resolve()
        ):
            score += 1000

        score -= abs(
            len(atoms)
            - EXPECTED_HBN_ATOMS
        )

        score -= abs(
            len(bonds)
            - EXPECTED_BONDS
        )

        if len(atoms) == EXPECTED_HBN_ATOMS:
            score += 500

        if len(bonds) == EXPECTED_BONDS:
            score += 500

        evaluated.append(
            (
                score,
                top_path,
                include_files,
                atoms,
                bonds,
                source_rows,
            )
        )

    if not evaluated:
        raise RuntimeError(
            "No topology candidate exposed a parsable HBN definition."
        )

    evaluated.sort(
        key=lambda item: item[0],
        reverse=True,
    )

    (
        _,
        top_path,
        include_files,
        atoms,
        bonds,
        source_rows,
    ) = evaluated[0]

    return (
        top_path,
        include_files,
        atoms,
        bonds,
        source_rows,
    )


def minimum_image(
    displacement: np.ndarray,
    box: np.ndarray,
) -> np.ndarray:
    return (
        displacement
        - box
        * np.round(
            displacement
            / box
        )
    )


def build_geometry_bond_set(
    positions: np.ndarray,
    elements: np.ndarray,
) -> set[tuple[int, int]]:
    b_indices = np.flatnonzero(
        elements == "B"
    )

    n_indices = np.flatnonzero(
        elements == "N"
    )

    displacement = (
        positions[
            b_indices,
            None,
            :,
        ]
        - positions[
            None,
            n_indices,
            :,
        ]
    )

    distances = np.linalg.norm(
        displacement,
        axis=2,
    )

    matches = np.argwhere(
        (
            distances >= BOND_MIN_NM
        )
        & (
            distances <= BOND_MAX_NM
        )
    )

    bonds: set[
        tuple[int, int]
    ] = set()

    for b_row, n_row in matches:
        first = int(
            b_indices[
                b_row
            ]
        ) + 1

        second = int(
            n_indices[
                n_row
            ]
        ) + 1

        bonds.add(
            tuple(
                sorted(
                    (
                        first,
                        second,
                    )
                )
            )
        )

    return bonds


def degree_array(
    bonds: set[tuple[int, int]],
) -> np.ndarray:
    degrees = np.zeros(
        EXPECTED_HBN_ATOMS,
        dtype=int,
    )

    for first, second in bonds:
        if not (
            1 <= first <= EXPECTED_HBN_ATOMS
            and 1 <= second <= EXPECTED_HBN_ATOMS
        ):
            raise RuntimeError(
                "HBN connectivity contains an out-of-range atom: "
                f"{first}-{second}"
            )

        degrees[
            first - 1
        ] += 1

        degrees[
            second - 1
        ] += 1

    return degrees


def bond_distances(
    bonds: set[tuple[int, int]],
    positions: np.ndarray,
) -> np.ndarray:
    values = []

    for first, second in sorted(
        bonds
    ):
        displacement = (
            positions[
                first - 1
            ]
            - positions[
                second - 1
            ]
        )

        values.append(
            float(
                np.linalg.norm(
                    displacement
                )
            )
        )

    return np.asarray(
        values,
        dtype=float,
    )


def end_metrics(
    label: str,
    indices: np.ndarray,
    positions: np.ndarray,
    elements: np.ndarray,
    center: np.ndarray,
    axis: np.ndarray,
    axial_coordinates: np.ndarray,
) -> tuple[
    dict[str, Any],
    list[dict[str, Any]],
]:
    if len(indices) == 0:
        raise RuntimeError(
            f"No terminal atoms found for {label}."
        )

    end_positions = positions[
        indices
    ]

    end_center = np.mean(
        end_positions,
        axis=0,
    )

    relative = (
        end_positions
        - end_center
    )

    axial_relative = (
        relative
        @ axis
    )

    perpendicular = (
        relative
        - np.outer(
            axial_relative,
            axis,
        )
    )

    radii = np.linalg.norm(
        perpendicular,
        axis=1,
    )

    rows = []

    for index, radius in zip(
        indices,
        radii,
    ):
        rows.append(
            {
                "end": label,
                "hbn_local_index_0based": int(
                    index
                ),
                "hbn_local_index_1based": int(
                    index + 1
                ),
                "element": str(
                    elements[
                        index
                    ]
                ),
                "coordination_number": 1,
                "x_nm": float(
                    positions[
                        index,
                        0
                    ]
                ),
                "y_nm": float(
                    positions[
                        index,
                        1
                    ]
                ),
                "z_nm": float(
                    positions[
                        index,
                        2
                    ]
                ),
                "axial_coordinate_nm": float(
                    axial_coordinates[
                        index
                    ]
                ),
                "radius_from_end_center_nm": float(
                    radius
                ),
            }
        )

    metrics = {
        "end": label,
        "terminal_atom_count": len(
            indices
        ),
        "B_count": int(
            np.count_nonzero(
                elements[
                    indices
                ]
                == "B"
            )
        ),
        "N_count": int(
            np.count_nonzero(
                elements[
                    indices
                ]
                == "N"
            )
        ),
        "axial_mean_nm": float(
            np.mean(
                axial_coordinates[
                    indices
                ]
            )
        ),
        "axial_minimum_nm": float(
            np.min(
                axial_coordinates[
                    indices
                ]
            )
        ),
        "axial_maximum_nm": float(
            np.max(
                axial_coordinates[
                    indices
                ]
            )
        ),
        "axial_standard_deviation_nm": float(
            np.std(
                axial_coordinates[
                    indices
                ]
            )
        ),
        "radius_mean_nm": float(
            np.mean(
                radii
            )
        ),
        "radius_standard_deviation_nm": float(
            np.std(
                radii
            )
        ),
        "radius_minimum_nm": float(
            np.min(
                radii
            )
        ),
        "radius_maximum_nm": float(
            np.max(
                radii
            )
        ),
    }

    return metrics, rows


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        AUDIT_SCRIPT,
        SYSTEM_GRO,
        TPR_DUMP,
    ):
        require_file(required)

    module = load_audit_module()

    atoms, box = module.read_gro(
        SYSTEM_GRO
    )

    dump_text = TPR_DUMP.read_text(
        encoding="utf-8",
        errors="replace",
    )

    hbn_tpr_atoms = (
        module.parse_hbn_tpr_atoms(
            dump_text
        )
    )

    if len(hbn_tpr_atoms) != EXPECTED_HBN_ATOMS:
        raise RuntimeError(
            "Unexpected HBN TPR atom count: "
            f"{len(hbn_tpr_atoms)}/"
            f"{EXPECTED_HBN_ATOMS}"
        )

    positions = np.asarray(
        [
            atoms[index][
                "position_nm"
            ]
            for index in range(
                EXPECTED_HBN_ATOMS
            )
        ],
        dtype=float,
    )

    elements = np.asarray(
        [
            row["element"]
            for row in hbn_tpr_atoms
        ],
        dtype=object,
    )

    (
        top_path,
        include_files,
        topology_atoms,
        explicit_bonds,
        source_rows,
    ) = find_accepted_topology()

    write_csv(
        TOPOLOGY_SOURCES_CSV,
        source_rows,
    )

    geometry_bonds = (
        build_geometry_bond_set(
            positions,
            elements,
        )
    )

    explicit_only = (
        explicit_bonds
        - geometry_bonds
    )

    geometry_only = (
        geometry_bonds
        - explicit_bonds
    )

    all_bonds = sorted(
        explicit_bonds
        | geometry_bonds
    )

    bond_comparison_rows = []

    for first, second in all_bonds:
        distance_nm = float(
            np.linalg.norm(
                positions[
                    first - 1
                ]
                - positions[
                    second - 1
                ]
            )
        )

        bond_comparison_rows.append(
            {
                "atom_i_1based": first,
                "atom_j_1based": second,
                "element_i": str(
                    elements[
                        first - 1
                    ]
                ),
                "element_j": str(
                    elements[
                        second - 1
                    ]
                ),
                "distance_nm": distance_nm,
                "in_explicit_topology": (
                    (
                        first,
                        second,
                    )
                    in explicit_bonds
                ),
                "in_geometry_graph": (
                    (
                        first,
                        second,
                    )
                    in geometry_bonds
                ),
                "classification": (
                    "MATCHED"
                    if (
                        first,
                        second,
                    )
                    in explicit_bonds
                    and (
                        first,
                        second,
                    )
                    in geometry_bonds
                    else (
                        "EXPLICIT_ONLY"
                        if (
                            first,
                            second,
                        )
                        in explicit_bonds
                        else "GEOMETRY_ONLY"
                    )
                ),
            }
        )

    write_csv(
        BOND_COMPARISON_CSV,
        bond_comparison_rows,
    )

    explicit_degrees = degree_array(
        explicit_bonds
    )

    geometry_degrees = degree_array(
        geometry_bonds
    )

    explicit_distances = bond_distances(
        explicit_bonds,
        positions,
    )

    (
        tube_center,
        tube_axis,
        pca_eigenvalues,
    ) = module.determine_tube_axis(
        positions
    )

    axial_coordinates = (
        positions
        - tube_center
    ) @ tube_axis

    terminal_indices = np.flatnonzero(
        explicit_degrees == 1
    )

    terminal_order = np.argsort(
        axial_coordinates[
            terminal_indices
        ]
    )

    ordered_terminal_indices = (
        terminal_indices[
            terminal_order
        ]
    )

    half = (
        len(
            ordered_terminal_indices
        )
        // 2
    )

    lower_indices = (
        ordered_terminal_indices[
            :half
        ]
    )

    upper_indices = (
        ordered_terminal_indices[
            half:
        ]
    )

    (
        lower_metrics,
        lower_rows,
    ) = end_metrics(
        "LOWER",
        lower_indices,
        positions,
        elements,
        tube_center,
        tube_axis,
        axial_coordinates,
    )

    (
        upper_metrics,
        upper_rows,
    ) = end_metrics(
        "UPPER",
        upper_indices,
        positions,
        elements,
        tube_center,
        tube_axis,
        axial_coordinates,
    )

    terminal_rows = (
        lower_rows
        + upper_rows
    )

    write_csv(
        TERMINAL_ATOMS_CSV,
        terminal_rows,
    )

    write_csv(
        END_SUMMARY_CSV,
        [
            lower_metrics,
            upper_metrics,
        ],
    )

    degree_counts = {
        degree: int(
            np.count_nonzero(
                explicit_degrees
                == degree
            )
        )
        for degree in range(
            int(
                np.max(
                    explicit_degrees
                )
            )
            + 1
        )
    }

    geometry_degree_counts = {
        degree: int(
            np.count_nonzero(
                geometry_degrees
                == degree
            )
        )
        for degree in range(
            int(
                np.max(
                    geometry_degrees
                )
            )
            + 1
        )
    }

    global_axial_min = float(
        np.min(
            axial_coordinates
        )
    )

    global_axial_max = float(
        np.max(
            axial_coordinates
        )
    )

    axial_span = (
        global_axial_max
        - global_axial_min
    )

    lower_localized_at_end = (
        float(
            lower_metrics[
                "axial_maximum_nm"
            ]
        )
        <= (
            global_axial_min
            + 0.10
            * axial_span
        )
    )

    upper_localized_at_end = (
        float(
            upper_metrics[
                "axial_minimum_nm"
            ]
        )
        >= (
            global_axial_max
            - 0.10
            * axial_span
        )
    )

    terminal_total_B = int(
        np.count_nonzero(
            elements[
                terminal_indices
            ]
            == "B"
        )
    )

    terminal_total_N = int(
        np.count_nonzero(
            elements[
                terminal_indices
            ]
            == "N"
        )
    )

    all_explicit_bonds_are_BN = all(
        elements[
            first - 1
        ]
        != elements[
            second - 1
        ]
        for first, second in explicit_bonds
    )

    gates = {
        "accepted_topology_contains_1680_HBN_atoms": (
            len(
                topology_atoms
            )
            == EXPECTED_HBN_ATOMS
        ),
        "accepted_topology_contains_2460_HBN_bonds": (
            len(
                explicit_bonds
            )
            == EXPECTED_BONDS
        ),
        "geometry_contains_2460_HBN_bonds": (
            len(
                geometry_bonds
            )
            == EXPECTED_BONDS
        ),
        "explicit_and_geometry_bond_sets_are_identical": (
            len(
                explicit_only
            )
            == 0
            and len(
                geometry_only
            )
            == 0
        ),
        "explicit_and_geometry_degree_arrays_are_identical": (
            np.array_equal(
                explicit_degrees,
                geometry_degrees,
            )
        ),
        "explicit_topology_has_60_degree1_atoms": (
            degree_counts.get(
                1,
                0,
            )
            == EXPECTED_TERMINAL_ATOMS
        ),
        "explicit_topology_has_zero_degree2_atoms": (
            degree_counts.get(
                2,
                0,
            )
            == 0
        ),
        "explicit_topology_has_1620_degree3_atoms": (
            degree_counts.get(
                3,
                0,
            )
            == EXPECTED_INTERIOR_ATOMS
        ),
        "explicit_topology_has_no_degree0_atoms": (
            degree_counts.get(
                0,
                0,
            )
            == 0
        ),
        "explicit_topology_has_no_degree4plus_atoms": (
            int(
                np.count_nonzero(
                    explicit_degrees
                    >= 4
                )
            )
            == 0
        ),
        "lower_end_has_30_terminal_atoms": (
            len(
                lower_indices
            )
            == EXPECTED_TERMINAL_ATOMS_PER_END
        ),
        "upper_end_has_30_terminal_atoms": (
            len(
                upper_indices
            )
            == EXPECTED_TERMINAL_ATOMS_PER_END
        ),
        "terminal_population_contains_30B_and_30N": (
            terminal_total_B == 30
            and terminal_total_N == 30
        ),
        "lower_terminal_cluster_is_at_axial_minimum": (
            lower_localized_at_end
        ),
        "upper_terminal_cluster_is_at_axial_maximum": (
            upper_localized_at_end
        ),
        "terminal_end_radii_are_symmetric": (
            abs(
                float(
                    lower_metrics[
                        "radius_mean_nm"
                    ]
                )
                - float(
                    upper_metrics[
                        "radius_mean_nm"
                    ]
                )
            )
            <= MAX_END_RADIUS_DIFFERENCE_NM
        ),
        "all_explicit_connectivity_is_BN": (
            all_explicit_bonds_are_BN
        ),
        "explicit_bond_lengths_are_chemically_plausible": (
            len(
                explicit_distances
            )
            == EXPECTED_BONDS
            and float(
                np.min(
                    explicit_distances
                )
            )
            >= VALID_BOND_MIN_NM
            and float(
                np.max(
                    explicit_distances
                )
            )
            <= VALID_BOND_MAX_NM
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
        DECISION_PASS
        if accepted
        else
        "HBN_TERMINAL_COORDINATION_CLASSIFICATION_REQUIRES_REVIEW"
    )

    required_next_step = (
        "REPAIR_R2_PARENT_RIM_AUDITOR_FOR_EXPLICIT_DEGREE1_TERMINI"
        if accepted
        else
        "REVIEW_HBN_EXPLICIT_TOPOLOGY_CLASSIFICATION_FAILURES"
    )

    summary = {
        "decision": decision,
        "selected_topology": relative(
            top_path
        ),
        "local_include_files": len(
            include_files
        ),
        "HBN_topology_atoms": len(
            topology_atoms
        ),
        "explicit_bonds": len(
            explicit_bonds
        ),
        "geometry_bonds": len(
            geometry_bonds
        ),
        "explicit_only_bonds": len(
            explicit_only
        ),
        "geometry_only_bonds": len(
            geometry_only
        ),
        "explicit_degree0_atoms": (
            degree_counts.get(
                0,
                0,
            )
        ),
        "explicit_degree1_atoms": (
            degree_counts.get(
                1,
                0,
            )
        ),
        "explicit_degree2_atoms": (
            degree_counts.get(
                2,
                0,
            )
        ),
        "explicit_degree3_atoms": (
            degree_counts.get(
                3,
                0,
            )
        ),
        "geometry_degree_counts": json.dumps(
            geometry_degree_counts,
            sort_keys=True,
        ),
        "terminal_B_atoms_total": (
            terminal_total_B
        ),
        "terminal_N_atoms_total": (
            terminal_total_N
        ),
        "lower_terminal_atoms": len(
            lower_indices
        ),
        "lower_terminal_B_atoms": (
            lower_metrics[
                "B_count"
            ]
        ),
        "lower_terminal_N_atoms": (
            lower_metrics[
                "N_count"
            ]
        ),
        "lower_terminal_axial_mean_nm": (
            lower_metrics[
                "axial_mean_nm"
            ]
        ),
        "lower_terminal_axial_std_nm": (
            lower_metrics[
                "axial_standard_deviation_nm"
            ]
        ),
        "lower_terminal_radius_mean_nm": (
            lower_metrics[
                "radius_mean_nm"
            ]
        ),
        "upper_terminal_atoms": len(
            upper_indices
        ),
        "upper_terminal_B_atoms": (
            upper_metrics[
                "B_count"
            ]
        ),
        "upper_terminal_N_atoms": (
            upper_metrics[
                "N_count"
            ]
        ),
        "upper_terminal_axial_mean_nm": (
            upper_metrics[
                "axial_mean_nm"
            ]
        ),
        "upper_terminal_axial_std_nm": (
            upper_metrics[
                "axial_standard_deviation_nm"
            ]
        ),
        "upper_terminal_radius_mean_nm": (
            upper_metrics[
                "radius_mean_nm"
            ]
        ),
        "explicit_bond_mean_nm": float(
            np.mean(
                explicit_distances
            )
        ),
        "explicit_bond_median_nm": float(
            np.median(
                explicit_distances
            )
        ),
        "explicit_bond_minimum_nm": float(
            np.min(
                explicit_distances
            )
        ),
        "explicit_bond_maximum_nm": float(
            np.max(
                explicit_distances
            )
        ),
        "parent_auditor_repair_authorized": (
            accepted
        ),
        "geometry_generation_authorized": False,
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

    SUMMARY_JSON.write_text(
        json.dumps(
            {
                "summary": summary,
                "gates": gates,
                "lower_end": (
                    lower_metrics
                ),
                "upper_end": (
                    upper_metrics
                ),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed in gates.items()
    )

    REPORT_MD.write_text(
        f"""# HBN Explicit Topology and Terminal Coordination Classification

## Purpose

This stage determines whether the 60 degree-1 atoms identified from
geometry are also present in the accepted explicit HBN topology.

No coordinates or topology were modified. No minimization, MD, or QM
calculation was executed.

## Accepted topology

- Selected topology:
  `{relative(top_path)}`
- Local files in include tree:
  **{len(include_files)}**
- HBN atom records:
  **{len(topology_atoms)}**
- Explicit HBN connections:
  **{len(explicit_bonds)}**

## Explicit versus geometric connectivity

- Geometric connections:
  **{len(geometry_bonds)}**
- Explicit-only connections:
  **{len(explicit_only)}**
- Geometry-only connections:
  **{len(geometry_only)}**
- Bond sets identical:
  **{'YES' if not explicit_only and not geometry_only else 'NO'}**

## Explicit coordination distribution

- Degree 0:
  **{degree_counts.get(0, 0)}**
- Degree 1:
  **{degree_counts.get(1, 0)}**
- Degree 2:
  **{degree_counts.get(2, 0)}**
- Degree 3:
  **{degree_counts.get(3, 0)}**
- Degree 4 or greater:
  **{np.count_nonzero(explicit_degrees >= 4)}**

The accepted scaffold therefore contains 60 singly coordinated
terminal sites and 1620 three-coordinate interior sites. The prior
assumption of 120 two-coordinate terminal sites is not applicable to
this structure.

## Terminal-site distribution

### Lower end

- Terminal atoms:
  **{len(lower_indices)}**
- B/N:
  **{lower_metrics['B_count']}/{lower_metrics['N_count']}**
- Axial mean/std:
  **{lower_metrics['axial_mean_nm']:.6f}/
  {lower_metrics['axial_standard_deviation_nm']:.6f} nm**
- Mean radius:
  **{lower_metrics['radius_mean_nm']:.6f} nm**

### Upper end

- Terminal atoms:
  **{len(upper_indices)}**
- B/N:
  **{upper_metrics['B_count']}/{upper_metrics['N_count']}**
- Axial mean/std:
  **{upper_metrics['axial_mean_nm']:.6f}/
  {upper_metrics['axial_standard_deviation_nm']:.6f} nm**
- Mean radius:
  **{upper_metrics['radius_mean_nm']:.6f} nm**

### Total terminal composition

- B/N:
  **{terminal_total_B}/{terminal_total_N}**

## Bond distances

- Mean/median:
  **{np.mean(explicit_distances):.6f}/
  {np.median(explicit_distances):.6f} nm**
- Minimum/maximum:
  **{np.min(explicit_distances):.6f}/
  {np.max(explicit_distances):.6f} nm**

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Parent-rim auditor repair authorized:
  **{'YES' if accepted else 'NO'}**
- Explicit geometry generation authorized:
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
        "Day024 HBN explicit-topology and terminal-site "
        "classification completed."
    )

    print(
        "Selected topology: "
        f"{relative(top_path)}"
    )

    print(
        "Topology include files / HBN atom records: "
        f"{len(include_files)}/"
        f"{len(topology_atoms)}"
    )

    print(
        "Explicit / geometry bonds: "
        f"{len(explicit_bonds)}/"
        f"{len(geometry_bonds)}"
    )

    print(
        "Explicit-only / geometry-only bonds: "
        f"{len(explicit_only)}/"
        f"{len(geometry_only)}"
    )

    print(
        "Explicit degree counts d0/d1/d2/d3/d4+: "
        f"{degree_counts.get(0, 0)}/"
        f"{degree_counts.get(1, 0)}/"
        f"{degree_counts.get(2, 0)}/"
        f"{degree_counts.get(3, 0)}/"
        f"{int(np.count_nonzero(explicit_degrees >= 4))}"
    )

    print(
        "Geometry degree counts d0/d1/d2/d3/d4+: "
        f"{geometry_degree_counts.get(0, 0)}/"
        f"{geometry_degree_counts.get(1, 0)}/"
        f"{geometry_degree_counts.get(2, 0)}/"
        f"{geometry_degree_counts.get(3, 0)}/"
        f"{int(np.count_nonzero(geometry_degrees >= 4))}"
    )

    print(
        "Terminal B / N atoms total: "
        f"{terminal_total_B}/"
        f"{terminal_total_N}"
    )

    print(
        "Lower terminal atoms / B / N / axial mean / "
        "axial std / radius: "
        f"{len(lower_indices)}/"
        f"{lower_metrics['B_count']}/"
        f"{lower_metrics['N_count']}/"
        f"{lower_metrics['axial_mean_nm']:.6f}/"
        f"{lower_metrics['axial_standard_deviation_nm']:.6f}/"
        f"{lower_metrics['radius_mean_nm']:.6f}"
    )

    print(
        "Upper terminal atoms / B / N / axial mean / "
        "axial std / radius: "
        f"{len(upper_indices)}/"
        f"{upper_metrics['B_count']}/"
        f"{upper_metrics['N_count']}/"
        f"{upper_metrics['axial_mean_nm']:.6f}/"
        f"{upper_metrics['axial_standard_deviation_nm']:.6f}/"
        f"{upper_metrics['radius_mean_nm']:.6f}"
    )

    print(
        "Explicit bond mean/median/min/max: "
        f"{np.mean(explicit_distances):.6f}/"
        f"{np.median(explicit_distances):.6f}/"
        f"{np.min(explicit_distances):.6f}/"
        f"{np.max(explicit_distances):.6f} nm"
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
        "Parent-rim auditor repair authorized: "
        f"{'YES' if accepted else 'NO'}"
    )

    print(
        "Explicit geometry generation authorized: NO"
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
        TERMINAL_ATOMS_CSV,
        END_SUMMARY_CSV,
        BOND_COMPARISON_CSV,
        TOPOLOGY_SOURCES_CSV,
        SUMMARY_CSV,
        GATES_CSV,
        SUMMARY_JSON,
        REPORT_MD,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not accepted:
        raise RuntimeError(
            "HBN explicit-topology classification "
            "requires review."
        )


if __name__ == "__main__":
    main()
