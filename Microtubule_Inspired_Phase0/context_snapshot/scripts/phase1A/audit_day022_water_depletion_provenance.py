#!/usr/bin/env python3

from __future__ import annotations

import csv
import importlib.util
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

PROTOCOL = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

EXECUTION = PROTOCOL / "execution"

FROZEN_RUN = (
    ROOT
    / "runs/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032_"
    "nvt_100ps_frozenSolute"
)

FROZEN_ENDPOINT = (
    FROZEN_RUN
    / "nvt_100ps_frozenSolute.gro"
)

HBN_ITP = (
    PROTOCOL
    / "protocol_inputs/topology/"
    "hbn_bonded_mobile_release.itp"
)

VALIDATOR = (
    ROOT
    / "scripts/phase1A/"
    "validate_day022_stage08_mobile_100ps.py"
)

WATER_ANALYSIS = (
    EXECUTION
    / "08_nvt_mobile_100ps/"
    "mobile_vs_frozen_water"
)

TIMESERIES_CSV = (
    WATER_ANALYSIS
    / "mobile_frozen_water_timeseries.csv"
)

OUTPUT_ENDPOINT_CSV = (
    WATER_ANALYSIS
    / "water_depletion_endpoint_audit.csv"
)

OUTPUT_BLOCK_CSV = (
    WATER_ANALYSIS
    / "water_depletion_10ps_block_audit.csv"
)

OUTPUT_SUMMARY_CSV = (
    WATER_ANALYSIS
    / "water_depletion_provenance_summary.csv"
)

OUTPUT_REPORT = (
    WATER_ANALYSIS
    / "WATER_DEPLETION_PROVENANCE_AUDIT_DAY022.md"
)

HBN_COUNT = 1680
PYR_COUNT = 4
PYR_ATOMS = 26
SOLUTE_COUNT = HBN_COUNT + PYR_COUNT * PYR_ATOMS

TOTAL_ATOMS = 68320
WATER_ATOMS = TOTAL_ATOMS - SOLUTE_COUNT
WATER_ATOMS_PER_MOLECULE = 4
WATER_COUNT = WATER_ATOMS // WATER_ATOMS_PER_MOLECULE


STAGES = [
    {
        "label": "accepted_frozen_endpoint",
        "elapsed_mobile_ps": 0.0,
        "path": FROZEN_ENDPOINT,
        "required": True,
    },
    {
        "label": "02_nvt_k10000_1ps",
        "elapsed_mobile_ps": 1.0,
        "path": (
            EXECUTION
            / "02_nvt_k10000_1ps/"
            "02_nvt_k10000_1ps.gro"
        ),
        "required": False,
    },
    {
        "label": "03_nvt_k1000_2ps",
        "elapsed_mobile_ps": 3.0,
        "path": (
            EXECUTION
            / "03_nvt_k1000_2ps/"
            "03_nvt_k1000_2ps.gro"
        ),
        "required": False,
    },
    {
        "label": "03b_nvt_k1000_hold_2ps",
        "elapsed_mobile_ps": 5.0,
        "path": (
            EXECUTION
            / "03b_nvt_k1000_hold_2ps/"
            "03b_nvt_k1000_hold_2ps.gro"
        ),
        "required": False,
    },
    {
        "label": "04_nvt_k100_2ps",
        "elapsed_mobile_ps": 7.0,
        "path": (
            EXECUTION
            / "04_nvt_k100_2ps/"
            "04_nvt_k100_2ps.gro"
        ),
        "required": False,
    },
    {
        "label": "05_nvt_unrestrained_2ps",
        "elapsed_mobile_ps": 9.0,
        "path": (
            EXECUTION
            / "05_nvt_unrestrained_2ps/"
            "05_nvt_unrestrained_2ps.gro"
        ),
        "required": False,
    },
    {
        "label": "06_nvt_unrestrained_10ps",
        "elapsed_mobile_ps": 19.0,
        "path": (
            EXECUTION
            / "06_nvt_unrestrained_10ps/"
            "06_nvt_unrestrained_10ps.gro"
        ),
        "required": False,
    },
    {
        "label": "07_nvt_unrestrained_25ps",
        "elapsed_mobile_ps": 44.0,
        "path": (
            EXECUTION
            / "07_nvt_unrestrained_25ps/"
            "07_nvt_unrestrained_25ps.gro"
        ),
        "required": True,
    },
    {
        "label": "08_nvt_mobile_100ps",
        "elapsed_mobile_ps": 144.0,
        "path": (
            EXECUTION
            / "08_nvt_mobile_100ps/"
            "08_nvt_mobile_100ps.gro"
        ),
        "required": True,
    },
]


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
            f"Missing or empty file: {path}"
        )


def load_module(path: Path, name: str):
    specification = (
        importlib.util.spec_from_file_location(
            name,
            path,
        )
    )

    if (
        specification is None
        or specification.loader is None
    ):
        raise RuntimeError(
            f"Could not load module: {path}"
        )

    module = importlib.util.module_from_spec(
        specification
    )

    specification.loader.exec_module(module)

    return module


def read_gro(
    path: Path,
) -> tuple[np.ndarray, np.ndarray]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    if len(lines) < 3:
        raise RuntimeError(
            f"Malformed GRO file: {path}"
        )

    natoms = int(
        lines[1].strip()
    )

    if natoms != TOTAL_ATOMS:
        raise RuntimeError(
            f"Expected {TOTAL_ATOMS} atoms in {path}; "
            f"found {natoms}"
        )

    positions = np.empty(
        (natoms, 3),
        dtype=float,
    )

    for atom_index, line in enumerate(
        lines[2 : 2 + natoms]
    ):
        positions[
            atom_index
        ] = (
            float(line[20:28]),
            float(line[28:36]),
            float(line[36:44]),
        )

    box_values = [
        float(value)
        for value in lines[
            2 + natoms
        ].split()
    ]

    if len(box_values) not in {
        3,
        9,
    }:
        raise RuntimeError(
            f"Unexpected box format in {path}"
        )

    if len(box_values) == 9:
        off_diagonal = np.array(
            box_values[3:],
            dtype=float,
        )

        if not np.allclose(
            off_diagonal,
            0.0,
            atol=1.0e-10,
        ):
            raise RuntimeError(
                "Only an orthorhombic box is supported"
            )

    box = np.array(
        box_values[:3],
        dtype=float,
    )

    return positions, box


def minimum_image(
    vectors: np.ndarray,
    box: np.ndarray,
) -> np.ndarray:
    return (
        vectors
        - box
        * np.round(
            vectors / box
        )
    )


def tube_axis(
    coordinates: np.ndarray,
) -> np.ndarray:
    centered = (
        coordinates
        - coordinates.mean(axis=0)
    )

    covariance = np.cov(
        centered,
        rowvar=False,
    )

    eigenvalues, eigenvectors = (
        np.linalg.eigh(
            covariance
        )
    )

    axis = eigenvectors[
        :,
        int(
            np.argmax(
                eigenvalues
            )
        ),
    ]

    return (
        axis
        / np.linalg.norm(axis)
    )


def endpoint_metrics(
    path: Path,
    integrated,
    hbn_topology,
) -> dict[str, float]:
    positions, box = read_gro(
        path
    )

    hbn = integrated.unwrap_by_bonds(
        positions[:HBN_COUNT],
        box,
        hbn_topology["bonds"],
    )

    axis = tube_axis(
        hbn
    )

    center = hbn.mean(
        axis=0
    )

    hbn_relative = (
        hbn - center
    )

    axial = (
        hbn_relative @ axis
    )

    transverse = (
        hbn_relative
        - np.outer(
            axial,
            axis,
        )
    )

    radial = np.linalg.norm(
        transverse,
        axis=1,
    )

    axial_lower = float(
        np.quantile(
            axial,
            0.01,
        )
    )

    axial_upper = float(
        np.quantile(
            axial,
            0.99,
        )
    )

    central_lower = float(
        np.quantile(
            axial,
            0.10,
        )
    )

    central_upper = float(
        np.quantile(
            axial,
            0.90,
        )
    )

    central_mask = (
        (axial >= central_lower)
        & (axial <= central_upper)
    )

    wall_radius = float(
        np.median(
            radial[
                central_mask
            ]
        )
    )

    lumen_length = (
        axial_upper
        - axial_lower
    )

    water_oxygen = positions[
        SOLUTE_COUNT:
        TOTAL_ATOMS:
        WATER_ATOMS_PER_MOLECULE
    ]

    if len(water_oxygen) != WATER_COUNT:
        raise RuntimeError(
            "Unexpected water-oxygen count"
        )

    center_wrapped = np.mod(
        center,
        box,
    )

    water_near_tube = (
        center
        + minimum_image(
            water_oxygen
            - center_wrapped,
            box,
        )
    )

    water_relative = (
        water_near_tube
        - center
    )

    water_axial = (
        water_relative @ axis
    )

    water_transverse = (
        water_relative
        - np.outer(
            water_axial,
            axis,
        )
    )

    water_radial = np.linalg.norm(
        water_transverse,
        axis=1,
    )

    inside = (
        (water_axial >= axial_lower)
        & (water_axial <= axial_upper)
        & (water_radial <= wall_radius)
    )

    occupancy = int(
        np.count_nonzero(
            inside
        )
    )

    volume = (
        math.pi
        * wall_radius ** 2
        * lumen_length
    )

    return {
        "lumen_water_count": occupancy,
        "lumen_number_density_nm-3": (
            occupancy / volume
        ),
        "wall_radius_nm": wall_radius,
        "lumen_length_nm": lumen_length,
        "lumen_volume_nm3": volume,
    }


def read_timeseries() -> list[dict[str, str]]:
    require_file(
        TIMESERIES_CSV
    )

    with TIMESERIES_CSV.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(
            csv.DictReader(handle)
        )

    if len(rows) != 402:
        raise RuntimeError(
            "Expected 402 water-timeseries rows; "
            f"found {len(rows)}"
        )

    return rows


def dataset_rows(
    rows: list[dict[str, str]],
    dataset: str,
) -> list[dict[str, str]]:
    selected = [
        row
        for row in rows
        if row[
            "dataset"
        ].strip() == dataset
    ]

    selected.sort(
        key=lambda row: float(
            row["relative_time_ps"]
        )
    )

    if len(selected) != 201:
        raise RuntimeError(
            f"Expected 201 {dataset} rows; "
            f"found {len(selected)}"
        )

    return selected


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
) -> None:
    if not rows:
        raise RuntimeError(
            f"No rows available for {path}"
        )

    fields: list[str] = []
    seen: set[str] = set()

    for row in rows:
        for field in row:
            if field not in seen:
                seen.add(field)
                fields.append(field)

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


def block_rows(
    dataset: str,
    rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    output = []

    times = np.array(
        [
            float(
                row[
                    "relative_time_ps"
                ]
            )
            for row in rows
        ],
        dtype=float,
    )

    occupancy = np.array(
        [
            float(
                row[
                    "lumen_water_count"
                ]
            )
            for row in rows
        ],
        dtype=float,
    )

    density = np.array(
        [
            float(
                row[
                    "lumen_number_density_nm-3"
                ]
            )
            for row in rows
        ],
        dtype=float,
    )

    for block_index in range(10):
        lower = 10.0 * block_index
        upper = lower + 10.0

        if block_index == 9:
            mask = (
                (times >= lower)
                & (times <= upper)
            )
        else:
            mask = (
                (times >= lower)
                & (times < upper)
            )

        output.append(
            {
                "dataset": dataset,
                "block": block_index + 1,
                "start_ps": lower,
                "end_ps": upper,
                "frame_count": int(
                    np.count_nonzero(
                        mask
                    )
                ),
                "occupancy_mean": float(
                    occupancy[mask].mean()
                ),
                "occupancy_std": float(
                    occupancy[mask].std()
                ),
                "occupancy_min": float(
                    occupancy[mask].min()
                ),
                "occupancy_max": float(
                    occupancy[mask].max()
                ),
                "density_mean_nm-3": float(
                    density[mask].mean()
                ),
            }
        )

    return output


def main() -> None:
    require_file(
        VALIDATOR
    )

    require_file(
        HBN_ITP
    )

    integrated = load_module(
        VALIDATOR,
        "stage08_water_audit_validator",
    )

    hbn_topology = integrated.parse_itp(
        HBN_ITP
    )

    endpoint_rows = []

    for stage in STAGES:
        path = stage["path"]

        if not path.exists():
            if stage["required"]:
                raise RuntimeError(
                    f"Missing required endpoint: {path}"
                )

            endpoint_rows.append(
                {
                    "stage": stage["label"],
                    "elapsed_mobile_ps": (
                        stage[
                            "elapsed_mobile_ps"
                        ]
                    ),
                    "status": "MISSING_OPTIONAL",
                    "gro_path": relative(path),
                }
            )

            continue

        metrics = endpoint_metrics(
            path,
            integrated,
            hbn_topology,
        )

        endpoint_rows.append(
            {
                "stage": stage["label"],
                "elapsed_mobile_ps": (
                    stage[
                        "elapsed_mobile_ps"
                    ]
                ),
                "status": "ANALYZED",
                **metrics,
                "gro_path": relative(path),
            }
        )

    write_csv(
        OUTPUT_ENDPOINT_CSV,
        endpoint_rows,
    )

    timeseries = read_timeseries()

    frozen_rows = dataset_rows(
        timeseries,
        "frozen",
    )

    mobile_rows = dataset_rows(
        timeseries,
        "mobile",
    )

    blocks = (
        block_rows(
            "frozen",
            frozen_rows,
        )
        + block_rows(
            "mobile",
            mobile_rows,
        )
    )

    write_csv(
        OUTPUT_BLOCK_CSV,
        blocks,
    )

    endpoint_by_stage = {
        row["stage"]: row
        for row in endpoint_rows
        if row.get(
            "status"
        ) == "ANALYZED"
    }

    frozen_endpoint_count = int(
        round(
            float(
                endpoint_by_stage[
                    "accepted_frozen_endpoint"
                ][
                    "lumen_water_count"
                ]
            )
        )
    )

    stage07_endpoint_count = int(
        round(
            float(
                endpoint_by_stage[
                    "07_nvt_unrestrained_25ps"
                ][
                    "lumen_water_count"
                ]
            )
        )
    )

    stage08_endpoint_count = int(
        round(
            float(
                endpoint_by_stage[
                    "08_nvt_mobile_100ps"
                ][
                    "lumen_water_count"
                ]
            )
        )
    )

    frozen_last_count = int(
        round(
            float(
                frozen_rows[-1][
                    "lumen_water_count"
                ]
            )
        )
    )

    mobile_first_count = int(
        round(
            float(
                mobile_rows[0][
                    "lumen_water_count"
                ]
            )
        )
    )

    mobile_last_count = int(
        round(
            float(
                mobile_rows[-1][
                    "lumen_water_count"
                ]
            )
        )
    )

    frozen_endpoint_difference = abs(
        frozen_endpoint_count
        - frozen_last_count
    )

    stage07_stage08_difference = abs(
        stage07_endpoint_count
        - mobile_first_count
    )

    stage08_endpoint_difference = abs(
        stage08_endpoint_count
        - mobile_last_count
    )

    geometry_count_consistency = (
        frozen_endpoint_difference <= 2
        and stage07_stage08_difference <= 2
        and stage08_endpoint_difference <= 2
    )

    summary = {
        "existing_frozen_window_relative_to_branch_start_ps": (
            "-100_to_0"
        ),
        "stage08_mobile_window_relative_to_branch_start_ps": (
            "44_to_144"
        ),
        "windows_temporally_matched": False,
        "accepted_frozen_endpoint_occupancy": (
            frozen_endpoint_count
        ),
        "frozen_timeseries_last_occupancy": (
            frozen_last_count
        ),
        "stage07_endpoint_occupancy": (
            stage07_endpoint_count
        ),
        "stage08_timeseries_first_occupancy": (
            mobile_first_count
        ),
        "stage08_endpoint_occupancy": (
            stage08_endpoint_count
        ),
        "stage08_timeseries_last_occupancy": (
            mobile_last_count
        ),
        "frozen_endpoint_count_difference": (
            frozen_endpoint_difference
        ),
        "stage07_to_stage08_frame0_count_difference": (
            stage07_stage08_difference
        ),
        "stage08_endpoint_count_difference": (
            stage08_endpoint_difference
        ),
        "geometry_count_consistency": (
            "PASS"
            if geometry_count_consistency
            else "REVIEW"
        ),
        "direct_mobility_attribution": (
            "NOT_SUPPORTED"
        ),
        "matched_frozen_continuation_required": (
            True
        ),
        "recommended_matched_frozen_duration_ps": (
            144.0
        ),
        "recommended_comparison_window_ps": (
            "44_to_144"
        ),
        "electronic_snapshot_selection_authorized": (
            False
        ),
        "longer_mobile_production_authorized": (
            False
        ),
    }

    write_csv(
        OUTPUT_SUMMARY_CSV,
        [summary],
    )

    OUTPUT_REPORT.write_text(
        f"""# Water-Depletion Provenance Audit

## Temporal provenance

- Existing frozen trajectory relative to the mobile branch point:
  **-100 to 0 ps**
- Stage08 mobile trajectory relative to the branch point:
  **44 to 144 ps**
- Temporally matched windows: **NO**
- Direct attribution of the water difference to solute mobility:
  **NOT SUPPORTED**

## Count continuity

- Accepted frozen endpoint / frozen trajectory last frame:
  {frozen_endpoint_count} / {frozen_last_count}
- Stage07 endpoint / Stage08 first frame:
  {stage07_endpoint_count} / {mobile_first_count}
- Stage08 endpoint / Stage08 last frame:
  {stage08_endpoint_count} / {mobile_last_count}
- Geometry/count consistency: **{summary['geometry_count_consistency']}**

## Required control

A matched frozen continuation must begin from the same branch
coordinates and use the same water-velocity initialization as the
mobile protocol. It should span 144 ps. The scientifically matched
comparison window is frozen 44-144 ps versus mobile Stage08 0-100 ps.

No electronic snapshot selection or recalculation is authorized from
the current unmatched comparison.
""",
        encoding="utf-8",
    )

    print(
        "Day022 water-depletion provenance audit completed."
    )

    print(
        "Existing frozen window relative to branch: -100 to 0 ps"
    )

    print(
        "Stage08 mobile window relative to branch: 44 to 144 ps"
    )

    print(
        "Accepted frozen endpoint / frozen last-frame occupancy: "
        f"{frozen_endpoint_count}/{frozen_last_count}"
    )

    print(
        "Stage07 endpoint / Stage08 first-frame occupancy: "
        f"{stage07_endpoint_count}/{mobile_first_count}"
    )

    print(
        "Stage08 endpoint / Stage08 last-frame occupancy: "
        f"{stage08_endpoint_count}/{mobile_last_count}"
    )

    print(
        "Geometry/count consistency: "
        f"{summary['geometry_count_consistency']}"
    )

    print(
        "Direct mobility attribution: NOT_SUPPORTED"
    )

    print(
        "Matched frozen continuation required: YES"
    )

    print(
        "Recommended matched frozen duration: 144 ps"
    )

    print(
        "Matched comparison window: frozen 44-144 ps "
        "versus mobile Stage08 0-100 ps"
    )

    print(
        "Electronic snapshot selection authorized: NO"
    )

    print(
        f"Wrote: {relative(OUTPUT_REPORT)}"
    )


if __name__ == "__main__":
    main()
