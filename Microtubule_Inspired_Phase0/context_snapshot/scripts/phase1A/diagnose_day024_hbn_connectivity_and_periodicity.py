#!/usr/bin/env python3

from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

AUDIT_SCRIPT = (
    ROOT
    / "scripts/phase1A/"
    "audit_day024_r2_parent_rim_and_chemical_constraints.py"
)

STAGE_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "01_r2_parent_rim_chemical_audit"
)

OUTPUT_ROOT = (
    STAGE_ROOT
    / "connectivity_diagnostic"
)

SYSTEM_GRO = (
    ROOT
    / "runs/phase1A/day023_confinement_design/"
    "15_r2_frozen_solute_nvt_20ps_preparation/"
    "r2_frozen_solute_nvt_20ps_input.gro"
)

TPR_DUMP = (
    STAGE_ROOT
    / "r2_parent_system_tpr_dump.txt"
)

SWEEP_CSV = (
    OUTPUT_ROOT
    / "hbn_connectivity_cutoff_sweep.csv"
)

NEIGHBOR_CSV = (
    OUTPUT_ROOT
    / "hbn_opposite_element_neighbor_ranks.csv"
)

BEST_CSV = (
    OUTPUT_ROOT
    / "hbn_connectivity_best_candidates.csv"
)

SUMMARY_JSON = (
    OUTPUT_ROOT
    / "hbn_connectivity_diagnostic_summary.json"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "HBN_CONNECTIVITY_AND_PERIODICITY_DIAGNOSTIC_DAY024.md"
)

EXPECTED_HBN_ATOMS = 1680
EXPECTED_B_ATOMS = 840
EXPECTED_N_ATOMS = 840
EXPECTED_EDGE_ATOMS = 120
EXPECTED_INTERIOR_ATOMS = 1560
EXPECTED_BONDS = 2460

ORIGINAL_LOWER_CUTOFF_NM = 0.115
ORIGINAL_UPPER_CUTOFF_NM = 0.175


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
                    field: row.get(field, "")
                    for field in fields
                }
            )


def load_audit_module():
    spec = importlib.util.spec_from_file_location(
        "day024_rim_audit",
        AUDIT_SCRIPT,
    )

    if (
        spec is None
        or spec.loader is None
    ):
        raise RuntimeError(
            "Could not load the Day024 audit module."
        )

    module = importlib.util.module_from_spec(
        spec
    )

    spec.loader.exec_module(
        module
    )

    return module


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


def summarize_values(
    values: np.ndarray,
) -> dict[str, float]:
    return {
        "minimum": float(
            np.min(values)
        ),
        "p01": float(
            np.percentile(
                values,
                1,
            )
        ),
        "p05": float(
            np.percentile(
                values,
                5,
            )
        ),
        "median": float(
            np.median(values)
        ),
        "p95": float(
            np.percentile(
                values,
                95,
            )
        ),
        "p99": float(
            np.percentile(
                values,
                99,
            )
        ),
        "maximum": float(
            np.max(values)
        ),
        "mean": float(
            np.mean(values)
        ),
    }


def degree_metrics(
    distances: np.ndarray,
    cutoff_nm: float,
) -> dict[str, Any]:
    bonded = (
        distances
        <= cutoff_nm
    )

    b_degrees = np.count_nonzero(
        bonded,
        axis=1,
    )

    n_degrees = np.count_nonzero(
        bonded,
        axis=0,
    )

    degrees = np.concatenate(
        (
            b_degrees,
            n_degrees,
        )
    )

    degree_0 = int(
        np.count_nonzero(
            degrees == 0
        )
    )

    degree_1 = int(
        np.count_nonzero(
            degrees == 1
        )
    )

    degree_2 = int(
        np.count_nonzero(
            degrees == 2
        )
    )

    degree_3 = int(
        np.count_nonzero(
            degrees == 3
        )
    )

    degree_4plus = int(
        np.count_nonzero(
            degrees >= 4
        )
    )

    bond_count = int(
        np.count_nonzero(
            bonded
        )
    )

    score = (
        abs(
            degree_2
            - EXPECTED_EDGE_ATOMS
        )
        + abs(
            degree_3
            - EXPECTED_INTERIOR_ATOMS
        )
        + 5
        * (
            degree_0
            + degree_1
            + degree_4plus
        )
        + abs(
            bond_count
            - EXPECTED_BONDS
        )
    )

    return {
        "cutoff_nm": cutoff_nm,
        "bond_count": bond_count,
        "degree_0": degree_0,
        "degree_1": degree_1,
        "degree_2": degree_2,
        "degree_3": degree_3,
        "degree_4plus": degree_4plus,
        "maximum_degree": int(
            np.max(degrees)
        ),
        "mean_degree": float(
            np.mean(degrees)
        ),
        "expected_graph_match": (
            bond_count == EXPECTED_BONDS
            and degree_2 == EXPECTED_EDGE_ATOMS
            and degree_3 == EXPECTED_INTERIOR_ATOMS
            and degree_0 == 0
            and degree_1 == 0
            and degree_4plus == 0
        ),
        "graph_score": int(
            score
        ),
    }


def ranked_neighbor_rows(
    distances: np.ndarray,
    b_indices: np.ndarray,
    n_indices: np.ndarray,
    mode: str,
) -> tuple[
    list[dict[str, Any]],
    dict[str, dict[str, float]],
]:
    b_sorted = np.sort(
        distances,
        axis=1,
    )

    n_sorted = np.sort(
        distances.T,
        axis=1,
    )

    rows: list[dict[str, Any]] = []

    rank_summary: dict[
        str,
        dict[str, float]
    ] = {}

    for element, indices, sorted_distances in (
        (
            "B",
            b_indices,
            b_sorted,
        ),
        (
            "N",
            n_indices,
            n_sorted,
        ),
    ):
        for rank in range(4):
            values = sorted_distances[
                :,
                rank
            ]

            key = (
                f"{mode}_{element}_rank_{rank + 1}"
            )

            rank_summary[key] = (
                summarize_values(
                    values
                )
            )

            for local_row, atom_index in enumerate(
                indices
            ):
                rows.append(
                    {
                        "distance_mode": mode,
                        "element": element,
                        "hbn_local_index_0based": int(
                            atom_index
                        ),
                        "neighbor_rank": (
                            rank + 1
                        ),
                        "distance_nm": float(
                            values[
                                local_row
                            ]
                        ),
                    }
                )

    return rows, rank_summary


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
            "Unexpected parsed HBN atom count: "
            f"{len(hbn_tpr_atoms)}/"
            f"{EXPECTED_HBN_ATOMS}"
        )

    hbn_positions = np.asarray(
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

    b_indices = np.flatnonzero(
        elements == "B"
    )

    n_indices = np.flatnonzero(
        elements == "N"
    )

    if (
        len(b_indices) != EXPECTED_B_ATOMS
        or len(n_indices) != EXPECTED_N_ATOMS
    ):
        raise RuntimeError(
            "Unexpected B/N composition."
        )

    b_positions = hbn_positions[
        b_indices
    ]

    n_positions = hbn_positions[
        n_indices
    ]

    raw_displacement = (
        b_positions[
            :,
            None,
            :,
        ]
        - n_positions[
            None,
            :,
            :,
        ]
    )

    pbc_displacement = minimum_image(
        raw_displacement,
        box,
    )

    raw_distances = np.linalg.norm(
        raw_displacement,
        axis=2,
    )

    pbc_distances = np.linalg.norm(
        pbc_displacement,
        axis=2,
    )

    neighbor_rows = []
    rank_summaries = {}

    for mode, distance_matrix in (
        (
            "RAW_NO_PBC",
            raw_distances,
        ),
        (
            "MINIMUM_IMAGE_XYZ",
            pbc_distances,
        ),
    ):
        rows, summaries = (
            ranked_neighbor_rows(
                distance_matrix,
                b_indices,
                n_indices,
                mode,
            )
        )

        neighbor_rows.extend(
            rows
        )

        rank_summaries.update(
            summaries
        )

    write_csv(
        NEIGHBOR_CSV,
        neighbor_rows,
    )

    cutoffs = np.round(
        np.arange(
            0.050,
            0.251,
            0.0025,
        ),
        4,
    )

    sweep_rows: list[
        dict[str, Any]
    ] = []

    for mode, distance_matrix in (
        (
            "RAW_NO_PBC",
            raw_distances,
        ),
        (
            "MINIMUM_IMAGE_XYZ",
            pbc_distances,
        ),
    ):
        for cutoff in cutoffs:
            metrics = degree_metrics(
                distance_matrix,
                float(cutoff),
            )

            metrics[
                "distance_mode"
            ] = mode

            sweep_rows.append(
                metrics
            )

    write_csv(
        SWEEP_CSV,
        sweep_rows,
    )

    sorted_candidates = sorted(
        sweep_rows,
        key=lambda row: (
            int(
                row[
                    "graph_score"
                ]
            ),
            abs(
                float(
                    row[
                        "cutoff_nm"
                    ]
                )
                - 0.145
            ),
        ),
    )

    best_rows = (
        sorted_candidates[:20]
    )

    write_csv(
        BEST_CSV,
        best_rows,
    )

    original_raw = {
        **degree_metrics(
            raw_distances,
            ORIGINAL_UPPER_CUTOFF_NM,
        ),
        "distance_mode": (
            "RAW_NO_PBC"
        ),
    }

    original_pbc = {
        **degree_metrics(
            pbc_distances,
            ORIGINAL_UPPER_CUTOFF_NM,
        ),
        "distance_mode": (
            "MINIMUM_IMAGE_XYZ"
        ),
    }

    original_band_raw = (
        (
            raw_distances
            >= ORIGINAL_LOWER_CUTOFF_NM
        )
        & (
            raw_distances
            <= ORIGINAL_UPPER_CUTOFF_NM
        )
    )

    original_band_pbc = (
        (
            pbc_distances
            >= ORIGINAL_LOWER_CUTOFF_NM
        )
        & (
            pbc_distances
            <= ORIGINAL_UPPER_CUTOFF_NM
        )
    )

    artificial_pbc_pairs = (
        (
            pbc_distances
            <= ORIGINAL_UPPER_CUTOFF_NM
        )
        & (
            raw_distances
            > ORIGINAL_UPPER_CUTOFF_NM
        )
    )

    artificial_pbc_pair_count = int(
        np.count_nonzero(
            artificial_pbc_pairs
        )
    )

    very_short_raw_pairs = int(
        np.count_nonzero(
            raw_distances
            < 0.100
        )
    )

    very_short_pbc_pairs = int(
        np.count_nonzero(
            pbc_distances
            < 0.100
        )
    )

    (
        tube_center,
        tube_axis,
        pca_eigenvalues,
    ) = module.determine_tube_axis(
        hbn_positions
    )

    axial_coordinates = (
        hbn_positions
        - tube_center
    ) @ tube_axis

    coordinate_minimum = np.min(
        hbn_positions,
        axis=0,
    )

    coordinate_maximum = np.max(
        hbn_positions,
        axis=0,
    )

    coordinate_span = (
        coordinate_maximum
        - coordinate_minimum
    )

    projected_box_extent = float(
        np.sum(
            np.abs(
                tube_axis
            )
            * box
        )
    )

    axial_span = float(
        np.max(
            axial_coordinates
        )
        - np.min(
            axial_coordinates
        )
    )

    exact_matches = [
        row
        for row in sweep_rows
        if bool(
            row[
                "expected_graph_match"
            ]
        )
    ]

    best = best_rows[0]

    summary = {
        "system_box_nm": (
            box.tolist()
        ),
        "HBN_coordinate_minimum_nm": (
            coordinate_minimum.tolist()
        ),
        "HBN_coordinate_maximum_nm": (
            coordinate_maximum.tolist()
        ),
        "HBN_coordinate_span_nm": (
            coordinate_span.tolist()
        ),
        "tube_center_nm": (
            tube_center.tolist()
        ),
        "tube_axis": (
            tube_axis.tolist()
        ),
        "PCA_eigenvalues": (
            pca_eigenvalues.tolist()
        ),
        "tube_axial_span_nm": (
            axial_span
        ),
        "box_extent_projected_on_axis_nm": (
            projected_box_extent
        ),
        "original_band_raw_bond_count": int(
            np.count_nonzero(
                original_band_raw
            )
        ),
        "original_band_pbc_bond_count": int(
            np.count_nonzero(
                original_band_pbc
            )
        ),
        "original_upper_cutoff_raw": (
            original_raw
        ),
        "original_upper_cutoff_pbc": (
            original_pbc
        ),
        "artificial_PBC_pairs_below_0p175nm": (
            artificial_pbc_pair_count
        ),
        "raw_pairs_below_0p100nm": (
            very_short_raw_pairs
        ),
        "PBC_pairs_below_0p100nm": (
            very_short_pbc_pairs
        ),
        "exact_expected_graph_matches": (
            exact_matches
        ),
        "best_candidate": (
            best
        ),
        "neighbor_rank_summaries": (
            rank_summaries
        ),
        "MD_executed": False,
        "QM_executed": False,
        "audit_repair_authorized": False,
        "required_next_step": (
            "CLASSIFY_CONNECTIVITY_FAILURE_AND_REPAIR_AUDITOR"
        ),
    }

    SUMMARY_JSON.write_text(
        json.dumps(
            summary,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    exact_match_lines = (
        "\n".join(
            (
                f"- {row['distance_mode']}: "
                f"cutoff={row['cutoff_nm']:.4f} nm; "
                f"bonds={row['bond_count']}; "
                f"degree2={row['degree_2']}; "
                f"degree3={row['degree_3']}"
            )
            for row in exact_matches
        )
        if exact_matches
        else "- NONE"
    )

    best_lines = "\n".join(
        (
            f"- {row['distance_mode']}: "
            f"cutoff={row['cutoff_nm']:.4f} nm; "
            f"score={row['graph_score']}; "
            f"bonds={row['bond_count']}; "
            f"d0/d1/d2/d3/d4+="
            f"{row['degree_0']}/"
            f"{row['degree_1']}/"
            f"{row['degree_2']}/"
            f"{row['degree_3']}/"
            f"{row['degree_4plus']}"
        )
        for row in best_rows[:10]
    )

    REPORT_MD.write_text(
        f"""# HBN Connectivity and Periodicity Diagnostic

## Scope

This diagnostic investigates why the original Day024 auditor found no
degree-2 terminal atoms.

No minimization, molecular dynamics, topology generation, or quantum
calculation was executed.

## Geometry

- Box:
  **{box[0]:.6f}, {box[1]:.6f}, {box[2]:.6f} nm**
- HBN coordinate span:
  **{coordinate_span[0]:.6f}, {coordinate_span[1]:.6f},
  {coordinate_span[2]:.6f} nm**
- PCA tube axis:
  **{tube_axis[0]:.8f}, {tube_axis[1]:.8f},
  {tube_axis[2]:.8f}**
- Axial tube span:
  **{axial_span:.6f} nm**
- Box extent projected on tube axis:
  **{projected_box_extent:.6f} nm**

## Original geometric bond rule

Original search band:

- minimum: **{ORIGINAL_LOWER_CUTOFF_NM:.6f} nm**
- maximum: **{ORIGINAL_UPPER_CUTOFF_NM:.6f} nm**

Raw, nonperiodic bond count in this band:

- **{np.count_nonzero(original_band_raw)}**

Minimum-image XYZ bond count in this band:

- **{np.count_nonzero(original_band_pbc)}**

Degree graph using all B–N distances up to 0.175 nm:

- Raw:
  bonds={original_raw['bond_count']},
  d0/d1/d2/d3/d4+ =
  {original_raw['degree_0']}/
  {original_raw['degree_1']}/
  {original_raw['degree_2']}/
  {original_raw['degree_3']}/
  {original_raw['degree_4plus']}
- Minimum-image XYZ:
  bonds={original_pbc['bond_count']},
  d0/d1/d2/d3/d4+ =
  {original_pbc['degree_0']}/
  {original_pbc['degree_1']}/
  {original_pbc['degree_2']}/
  {original_pbc['degree_3']}/
  {original_pbc['degree_4plus']}

Potential PBC-created pairs below 0.175 nm:

- **{artificial_pbc_pair_count}**

B–N pairs below 0.100 nm:

- Raw: **{very_short_raw_pairs}**
- Minimum-image XYZ: **{very_short_pbc_pairs}**

## Exact expected graph matches

{exact_match_lines}

## Ten best cutoff/mode candidates

{best_lines}

## Status

- Auditor repair authorized: **NO**
- MD executed: **NO**
- QM executed: **NO**
- Required next step:
  `CLASSIFY_CONNECTIVITY_FAILURE_AND_REPAIR_AUDITOR`
""",
        encoding="utf-8",
    )

    print(
        "Day024 HBN connectivity and periodicity "
        "diagnostic completed."
    )

    print(
        "System box: "
        f"{box[0]:.6f}/"
        f"{box[1]:.6f}/"
        f"{box[2]:.6f} nm"
    )

    print(
        "HBN coordinate span: "
        f"{coordinate_span[0]:.6f}/"
        f"{coordinate_span[1]:.6f}/"
        f"{coordinate_span[2]:.6f} nm"
    )

    print(
        "PCA tube axis: "
        f"{tube_axis[0]:.8f}/"
        f"{tube_axis[1]:.8f}/"
        f"{tube_axis[2]:.8f}"
    )

    print(
        "Tube axial span / projected box extent: "
        f"{axial_span:.6f}/"
        f"{projected_box_extent:.6f} nm"
    )

    print(
        "Original 0.115-0.175 nm band "
        "raw/PBC bond counts: "
        f"{np.count_nonzero(original_band_raw)}/"
        f"{np.count_nonzero(original_band_pbc)}"
    )

    print(
        "At cutoff <=0.175 nm, RAW "
        "bonds and d0/d1/d2/d3/d4+: "
        f"{original_raw['bond_count']} | "
        f"{original_raw['degree_0']}/"
        f"{original_raw['degree_1']}/"
        f"{original_raw['degree_2']}/"
        f"{original_raw['degree_3']}/"
        f"{original_raw['degree_4plus']}"
    )

    print(
        "At cutoff <=0.175 nm, PBC "
        "bonds and d0/d1/d2/d3/d4+: "
        f"{original_pbc['bond_count']} | "
        f"{original_pbc['degree_0']}/"
        f"{original_pbc['degree_1']}/"
        f"{original_pbc['degree_2']}/"
        f"{original_pbc['degree_3']}/"
        f"{original_pbc['degree_4plus']}"
    )

    print(
        "Potential PBC-created pairs below 0.175 nm: "
        f"{artificial_pbc_pair_count}"
    )

    print(
        "B-N pairs below 0.100 nm RAW/PBC: "
        f"{very_short_raw_pairs}/"
        f"{very_short_pbc_pairs}"
    )

    for mode in (
        "RAW_NO_PBC",
        "MINIMUM_IMAGE_XYZ",
    ):
        for element in (
            "B",
            "N",
        ):
            values = rank_summaries[
                f"{mode}_{element}_rank_1"
            ]

            print(
                f"{mode} {element} first-opposite-neighbor "
                "min/median/max: "
                f"{values['minimum']:.6f}/"
                f"{values['median']:.6f}/"
                f"{values['maximum']:.6f} nm"
            )

    print(
        "Exact expected graph matches: "
        f"{len(exact_matches)}"
    )

    for row in exact_matches[:10]:
        print(
            "EXACT "
            f"{row['distance_mode']} | "
            f"cutoff={row['cutoff_nm']:.4f} nm | "
            f"bonds={row['bond_count']} | "
            f"d2/d3="
            f"{row['degree_2']}/"
            f"{row['degree_3']}"
        )

    print(
        "Best graph candidates:"
    )

    for row in best_rows[:10]:
        print(
            f"  {row['distance_mode']} | "
            f"cutoff={row['cutoff_nm']:.4f} nm | "
            f"score={row['graph_score']} | "
            f"bonds={row['bond_count']} | "
            f"d0/d1/d2/d3/d4+="
            f"{row['degree_0']}/"
            f"{row['degree_1']}/"
            f"{row['degree_2']}/"
            f"{row['degree_3']}/"
            f"{row['degree_4plus']}"
        )

    print(
        "Auditor repair authorized: NO"
    )

    print(
        "MD executed: NO"
    )

    print(
        "QM executed: NO"
    )

    print(
        "Required next step: "
        "CLASSIFY_CONNECTIVITY_FAILURE_AND_REPAIR_AUDITOR"
    )

    for path in (
        SWEEP_CSV,
        NEIGHBOR_CSV,
        BEST_CSV,
        SUMMARY_JSON,
        REPORT_MD,
    ):
        print(
            f"Wrote: {relative(path)}"
        )


if __name__ == "__main__":
    main()
