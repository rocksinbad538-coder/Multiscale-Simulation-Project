#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import math
import os
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

ACCEPTED_R0 = (
    ROOT
    / "runs/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032_"
    "nvt_100ps_frozenSolute"
)

PROTOCOL_ROOT = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

STAGE08_ANALYSIS_ROOT = (
    PROTOCOL_ROOT
    / "execution/08_nvt_mobile_100ps"
)

DAY023_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design"
)

OUTPUT_ROOT = (
    DAY023_ROOT
    / "01_r0_t0_reference"
)

T0_GRO = (
    OUTPUT_ROOT
    / "r0_accepted_t0_hydrated_system.gro"
)

TRJCONV_LOG = (
    OUTPUT_ROOT
    / "r0_accepted_t0_hydrated_trjconv.log"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r0_t0_geometry_and_hydration_summary.csv"
)

END_ZONE_CSV = (
    OUTPUT_ROOT
    / "r0_t0_end_zone_water_counts.csv"
)

COMPOSITION_CSV = (
    OUTPUT_ROOT
    / "r0_t0_composition_audit.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R0_T0_HYDRATED_REFERENCE_AUDIT_DAY023.md"
)

CONTRACT_UPDATE_MD = (
    OUTPUT_ROOT
    / "R1_AUTHORITATIVE_START_STATE_DAY023.md"
)

EXPECTED_ATOMS = 68320
HBN_ATOMS = 1680
PYR_ATOMS = 104
SOLUTE_ATOMS = HBN_ATOMS + PYR_ATOMS
WATER_ATOMS = EXPECTED_ATOMS - SOLUTE_ATOMS
WATER_SITES = 4
EXPECTED_WATERS = WATER_ATOMS // WATER_SITES

FRAME_INTERVAL_PS = 0.5
TIME_TOLERANCE_PS = 1.0e-6
SOLUTE_MATCH_TOLERANCE_NM = 0.0020

END_RING_WIDTH_NM = 0.15
END_ZONE_WIDTH_NM = 0.30
PROVISIONAL_WALL_EXCLUSION_NM = 0.25


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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            block = handle.read(
                1024 * 1024
            )

            if not block:
                break

            digest.update(block)

    return digest.hexdigest()


def locate_gmx() -> str:
    executable = shutil.which("gmx")

    if executable:
        return executable

    fallback = Path(
        "/usr/local/gromacs/bin/gmx"
    )

    if fallback.exists():
        return str(fallback)

    raise RuntimeError(
        "Could not locate GROMACS."
    )


def discover_unique(
    directory: Path,
    suffix: str,
) -> Path:
    candidates = sorted(
        path
        for path in directory.glob(
            f"*{suffix}"
        )
        if path.is_file()
    )

    if len(candidates) != 1:
        raise RuntimeError(
            f"Expected exactly one {suffix} file in "
            f"{directory}; found {len(candidates)}:\n"
            + "\n".join(
                str(path)
                for path in candidates
            )
        )

    return candidates[0]


def parse_time_from_title(
    title: str,
) -> float | None:
    patterns = (
        r"\bt\s*=\s*([-+0-9.eE]+)",
        r"\btime\s*=\s*([-+0-9.eE]+)",
        r"\btime\s+([-+0-9.eE]+)",
    )

    for pattern in patterns:
        match = re.search(
            pattern,
            title,
            flags=re.IGNORECASE,
        )

        if match is not None:
            return float(
                match.group(1)
            )

    return None


def parse_box_matrix(
    values: list[float],
) -> np.ndarray:
    if len(values) == 3:
        return np.diag(
            np.asarray(
                values,
                dtype=float,
            )
        )

    if len(values) == 9:
        # GROMACS GRO triclinic order:
        # v1x v2y v3z v1y v1z v2x v2z v3x v3y
        return np.asarray(
            [
                [
                    values[0],
                    values[3],
                    values[4],
                ],
                [
                    values[5],
                    values[1],
                    values[6],
                ],
                [
                    values[7],
                    values[8],
                    values[2],
                ],
            ],
            dtype=float,
        )

    raise RuntimeError(
        "Unsupported GRO box field count: "
        f"{len(values)}"
    )


def read_gro(
    path: Path,
) -> tuple[
    str,
    list[dict[str, Any]],
    np.ndarray,
]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    if len(lines) < 3:
        raise RuntimeError(
            f"Malformed GRO file: {path}"
        )

    title = lines[0]
    natoms = int(
        lines[1].strip()
    )

    if len(lines) < natoms + 3:
        raise RuntimeError(
            f"Incomplete GRO file: {path}"
        )

    atoms: list[dict[str, Any]] = []

    for zero_index, line in enumerate(
        lines[2 : 2 + natoms]
    ):
        if len(line) < 44:
            raise RuntimeError(
                "Malformed atom line "
                f"{zero_index + 1} in {path}"
            )

        atoms.append(
            {
                "index": zero_index,
                "resid": int(
                    line[0:5]
                ),
                "resname": line[
                    5:10
                ].strip(),
                "atomname": line[
                    10:15
                ].strip(),
                "atomnum": int(
                    line[15:20]
                ),
                "position": np.asarray(
                    [
                        float(
                            line[20:28]
                        ),
                        float(
                            line[28:36]
                        ),
                        float(
                            line[36:44]
                        ),
                    ],
                    dtype=float,
                ),
            }
        )

    box_values = [
        float(value)
        for value in lines[
            2 + natoms
        ].split()
    ]

    box = parse_box_matrix(
        box_values
    )

    return (
        title,
        atoms,
        box,
    )


def orthorhombic_lengths(
    box: np.ndarray,
) -> np.ndarray:
    off_diagonal = (
        box
        - np.diag(
            np.diag(
                box
            )
        )
    )

    if not np.allclose(
        off_diagonal,
        0.0,
        atol=1.0e-8,
        rtol=0.0,
    ):
        raise RuntimeError(
            "The current audit requires an "
            "orthorhombic box."
        )

    lengths = np.diag(
        box
    ).copy()

    if np.any(
        lengths <= 0.0
    ):
        raise RuntimeError(
            f"Invalid box lengths: {lengths}"
        )

    return lengths


def minimum_image(
    displacement: np.ndarray,
    box_lengths: np.ndarray,
) -> np.ndarray:
    return (
        displacement
        - box_lengths
        * np.round(
            displacement
            / box_lengths
        )
    )


def periodic_cluster_coordinates(
    positions: np.ndarray,
    box_lengths: np.ndarray,
) -> np.ndarray:
    fractional = (
        positions
        / box_lengths
    )

    angles = (
        2.0
        * np.pi
        * fractional
    )

    circular_mean = np.arctan2(
        np.mean(
            np.sin(
                angles
            ),
            axis=0,
        ),
        np.mean(
            np.cos(
                angles
            ),
            axis=0,
        ),
    )

    center_fractional = (
        circular_mean
        / (
            2.0
            * np.pi
        )
    ) % 1.0

    delta_fractional = (
        fractional
        - center_fractional
    )

    delta_fractional -= np.round(
        delta_fractional
    )

    return (
        delta_fractional
        * box_lengths
    )


def deterministic_axis(
    coordinates: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    center = np.mean(
        coordinates,
        axis=0,
    )

    centered = (
        coordinates
        - center
    )

    covariance = np.cov(
        centered,
        rowvar=False,
        bias=True,
    )

    eigenvalues, eigenvectors = (
        np.linalg.eigh(
            covariance
        )
    )

    order = np.argsort(
        eigenvalues
    )[::-1]

    eigenvalues = eigenvalues[
        order
    ]

    eigenvectors = eigenvectors[
        :,
        order,
    ]

    axis = eigenvectors[
        :,
        0
    ]

    dominant_component = int(
        np.argmax(
            np.abs(
                axis
            )
        )
    )

    if axis[
        dominant_component
    ] < 0.0:
        axis = -axis

    axis /= np.linalg.norm(
        axis
    )

    return (
        center,
        axis,
        eigenvalues,
    )


def axial_radial_coordinates(
    coordinates: np.ndarray,
    center: np.ndarray,
    axis: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:
    relative_coordinates = (
        coordinates
        - center
    )

    axial = (
        relative_coordinates
        @ axis
    )

    perpendicular = (
        relative_coordinates
        - np.outer(
            axial,
            axis,
        )
    )

    radial = np.linalg.norm(
        perpendicular,
        axis=1,
    )

    return (
        axial,
        radial,
    )


def read_csv_rows(
    path: Path,
) -> list[dict[str, str]]:
    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        return list(
            csv.DictReader(
                handle
            )
        )


def discover_authoritative_timeseries() -> Path:
    preferred = (
        STAGE08_ANALYSIS_ROOT
        / "mobile_vs_frozen_water/"
        "mobile_frozen_water_timeseries.csv"
    )

    if preferred.exists():
        return preferred

    candidates = sorted(
        path
        for path in STAGE08_ANALYSIS_ROOT.rglob(
            "*water_timeseries.csv"
        )
        if (
            path.is_file()
            and "matched_mobile_vs_frozen_water"
            not in str(path)
        )
    )

    valid: list[Path] = []

    for path in candidates:
        try:
            rows = read_csv_rows(
                path
            )
        except Exception:
            continue

        frozen_rows = [
            row
            for row in rows
            if row.get(
                "dataset",
                "",
            ).strip().lower()
            == "frozen"
        ]

        if len(frozen_rows) != 201:
            continue

        try:
            times = sorted(
                float(
                    row[
                        "relative_time_ps"
                    ]
                )
                for row in frozen_rows
            )
        except (
            KeyError,
            ValueError,
        ):
            continue

        if (
            math.isclose(
                times[0],
                0.0,
                abs_tol=1.0e-9,
            )
            and math.isclose(
                times[-1],
                100.0,
                abs_tol=1.0e-9,
            )
        ):
            valid.append(
                path
            )

    if len(valid) != 1:
        raise RuntimeError(
            "Could not uniquely resolve the original "
            "accepted frozen-water time series.\n"
            + "\n".join(
                str(path)
                for path in valid
            )
        )

    return valid[0]


def extract_t0(
    gmx: str,
    xtc: Path,
    tpr: Path,
) -> None:
    if T0_GRO.exists():
        T0_GRO.unlink()

    completed = subprocess.run(
        [
            gmx,
            "trjconv",
            "-f",
            str(xtc),
            "-s",
            str(tpr),
            "-o",
            str(T0_GRO),
            "-dump",
            "0",
        ],
        cwd=OUTPUT_ROOT,
        env={
            **os.environ,
            "GMX_MAXBACKUP": "-1",
        },
        input="0\n",
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )

    TRJCONV_LOG.write_text(
        completed.stdout,
        encoding="utf-8",
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "gmx trjconv failed. "
            f"See {TRJCONV_LOG}"
        )

    require_file(
        T0_GRO
    )


def atom_positions(
    atoms: list[dict[str, Any]],
) -> np.ndarray:
    return np.asarray(
        [
            atom[
                "position"
            ]
            for atom in atoms
        ],
        dtype=float,
    )


def ordered_water_oxygens(
    water_atoms: list[dict[str, Any]],
) -> tuple[
    np.ndarray,
    bool,
    Counter[str],
]:
    if (
        len(
            water_atoms
        )
        % WATER_SITES
        != 0
    ):
        raise RuntimeError(
            "Water atom count is not divisible "
            f"by {WATER_SITES}."
        )

    oxygen_positions = []
    chunk_consistent = True
    first_site_names: Counter[str] = Counter()

    for start in range(
        0,
        len(
            water_atoms
        ),
        WATER_SITES,
    ):
        chunk = water_atoms[
            start : start
            + WATER_SITES
        ]

        residue_identity = {
            (
                atom[
                    "resid"
                ],
                atom[
                    "resname"
                ],
            )
            for atom in chunk
        }

        if len(
            residue_identity
        ) != 1:
            chunk_consistent = False

        oxygen_candidates = [
            atom
            for atom in chunk
            if (
                atom[
                    "atomname"
                ].upper().startswith(
                    "O"
                )
                or atom[
                    "atomname"
                ].upper()
                in {
                    "OW",
                    "OH2",
                }
            )
        ]

        if len(
            oxygen_candidates
        ) == 1:
            oxygen = oxygen_candidates[
                0
            ]
        else:
            oxygen = chunk[
                0
            ]

        first_site_names[
            chunk[
                0
            ][
                "atomname"
            ]
        ] += 1

        oxygen_positions.append(
            oxygen[
                "position"
            ]
        )

    return (
        np.asarray(
            oxygen_positions,
            dtype=float,
        ),
        chunk_consistent,
        first_site_names,
    )


def nearest_distances_to_subset(
    query: np.ndarray,
    reference: np.ndarray,
    box_lengths: np.ndarray,
    query_chunk: int = 512,
) -> np.ndarray:
    if (
        len(
            query
        )
        == 0
        or len(
            reference
        )
        == 0
    ):
        return np.asarray(
            [],
            dtype=float,
        )

    minima = np.full(
        len(
            query
        ),
        np.inf,
        dtype=float,
    )

    for start in range(
        0,
        len(
            query
        ),
        query_chunk,
    ):
        stop = min(
            start
            + query_chunk,
            len(
                query
            ),
        )

        difference = (
            query[
                start:stop,
                None,
                :,
            ]
            - reference[
                None,
                :,
                :,
            ]
        )

        difference = minimum_image(
            difference,
            box_lengths,
        )

        distance_squared = np.sum(
            difference
            * difference,
            axis=2,
        )

        minima[
            start:stop
        ] = np.sqrt(
            np.min(
                distance_squared,
                axis=1,
            )
        )

    return minima


def write_single_row_csv(
    path: Path,
    row: dict[str, Any],
) -> None:
    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                row.keys()
            ),
        )

        writer.writeheader()
        writer.writerow(
            row
        )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    accepted_xtc = discover_unique(
        ACCEPTED_R0,
        ".xtc",
    )

    accepted_tpr = discover_unique(
        ACCEPTED_R0,
        ".tpr",
    )

    accepted_gro = discover_unique(
        ACCEPTED_R0,
        ".gro",
    )

    authoritative_timeseries = (
        discover_authoritative_timeseries()
    )

    gmx = locate_gmx()

    extract_t0(
        gmx,
        accepted_xtc,
        accepted_tpr,
    )

    (
        t0_title,
        t0_atoms,
        t0_box,
    ) = read_gro(
        T0_GRO
    )

    (
        accepted_gro_title,
        accepted_gro_atoms,
        accepted_gro_box,
    ) = read_gro(
        accepted_gro
    )

    if len(
        t0_atoms
    ) != EXPECTED_ATOMS:
        raise RuntimeError(
            "Unexpected atom count in extracted "
            f"t=0 GRO: {len(t0_atoms)}/"
            f"{EXPECTED_ATOMS}"
        )

    if len(
        accepted_gro_atoms
    ) != EXPECTED_ATOMS:
        raise RuntimeError(
            "Unexpected atom count in accepted GRO: "
            f"{len(accepted_gro_atoms)}/"
            f"{EXPECTED_ATOMS}"
        )

    t0_time = parse_time_from_title(
        t0_title
    )

    accepted_gro_time = (
        parse_time_from_title(
            accepted_gro_title
        )
    )

    if (
        t0_time is not None
        and not math.isclose(
            t0_time,
            0.0,
            abs_tol=TIME_TOLERANCE_PS,
        )
    ):
        raise RuntimeError(
            "Extracted frame title does not correspond "
            f"to t=0 ps: {t0_time}"
        )

    t0_positions = atom_positions(
        t0_atoms
    )

    accepted_gro_positions = (
        atom_positions(
            accepted_gro_atoms
        )
    )

    box_lengths = orthorhombic_lengths(
        t0_box
    )

    accepted_box_lengths = (
        orthorhombic_lengths(
            accepted_gro_box
        )
    )

    box_difference = float(
        np.max(
            np.abs(
                box_lengths
                - accepted_box_lengths
            )
        )
    )

    solute_displacement = (
        t0_positions[
            :SOLUTE_ATOMS
        ]
        - accepted_gro_positions[
            :SOLUTE_ATOMS
        ]
    )

    solute_displacement = minimum_image(
        solute_displacement,
        box_lengths,
    )

    solute_distances = np.linalg.norm(
        solute_displacement,
        axis=1,
    )

    solute_rms_difference = float(
        np.sqrt(
            np.mean(
                solute_distances
                * solute_distances
            )
        )
    )

    solute_max_difference = float(
        np.max(
            solute_distances
        )
    )

    hbn_atoms = t0_atoms[
        :HBN_ATOMS
    ]

    pyr_atoms = t0_atoms[
        HBN_ATOMS:
        SOLUTE_ATOMS
    ]

    water_atoms = t0_atoms[
        SOLUTE_ATOMS:
    ]

    (
        water_oxygen_positions,
        water_chunks_consistent,
        first_water_site_names,
    ) = ordered_water_oxygens(
        water_atoms
    )

    if len(
        water_oxygen_positions
    ) != EXPECTED_WATERS:
        raise RuntimeError(
            "Unexpected water count: "
            f"{len(water_oxygen_positions)}/"
            f"{EXPECTED_WATERS}"
        )

    hbn_positions_wrapped = atom_positions(
        hbn_atoms
    )

    hbn_coordinates = (
        periodic_cluster_coordinates(
            hbn_positions_wrapped,
            box_lengths,
        )
    )

    (
        hbn_center_local,
        tube_axis,
        eigenvalues,
    ) = deterministic_axis(
        hbn_coordinates
    )

    (
        hbn_axial,
        hbn_radial,
    ) = axial_radial_coordinates(
        hbn_coordinates,
        hbn_center_local,
        tube_axis,
    )

    axial_low = float(
        np.quantile(
            hbn_axial,
            0.01,
        )
    )

    axial_high = float(
        np.quantile(
            hbn_axial,
            0.99,
        )
    )

    tube_length_p98 = (
        axial_high
        - axial_low
    )

    wall_radius_mean = float(
        np.mean(
            hbn_radial
        )
    )

    wall_radius_median = float(
        np.median(
            hbn_radial
        )
    )

    wall_radius_q01 = float(
        np.quantile(
            hbn_radial,
            0.01,
        )
    )

    wall_radius_q99 = float(
        np.quantile(
            hbn_radial,
            0.99,
        )
    )

    accessible_radius_provisional = max(
        0.0,
        wall_radius_median
        - PROVISIONAL_WALL_EXCLUSION_NM,
    )

    accessible_volume_provisional = (
        math.pi
        * accessible_radius_provisional
        * accessible_radius_provisional
        * tube_length_p98
    )

    eigenvalue_ratio = float(
        eigenvalues[
            0
        ]
        / max(
            eigenvalues[
                1
            ],
            1.0e-15,
        )
    )

    # Transform water oxygen positions into the same
    # periodic local coordinate system as HBN.
    hbn_center_wrapped = np.mod(
        (
            hbn_positions_wrapped[
                0
            ]
            - hbn_coordinates[
                0
            ]
            + hbn_center_local
        ),
        box_lengths,
    )

    water_relative = minimum_image(
        water_oxygen_positions
        - hbn_center_wrapped,
        box_lengths,
    )

    water_axial = (
        water_relative
        @ tube_axis
    )

    water_perpendicular = (
        water_relative
        - np.outer(
            water_axial,
            tube_axis,
        )
    )

    water_radial = np.linalg.norm(
        water_perpendicular,
        axis=1,
    )

    independent_lumen_mask = (
        (
            water_axial
            >= axial_low
        )
        & (
            water_axial
            <= axial_high
        )
        & (
            water_radial
            <= accessible_radius_provisional
        )
    )

    independent_lumen_count = int(
        np.count_nonzero(
            independent_lumen_mask
        )
    )

    lower_end_ring_mask = (
        np.abs(
            hbn_axial
            - axial_low
        )
        <= END_RING_WIDTH_NM
    )

    upper_end_ring_mask = (
        np.abs(
            hbn_axial
            - axial_high
        )
        <= END_RING_WIDTH_NM
    )

    lower_end_ring_positions = (
        hbn_positions_wrapped[
            lower_end_ring_mask
        ]
    )

    upper_end_ring_positions = (
        hbn_positions_wrapped[
            upper_end_ring_mask
        ]
    )

    lower_near_axial = (
        np.abs(
            water_axial
            - axial_low
        )
        <= END_ZONE_WIDTH_NM
    )

    upper_near_axial = (
        np.abs(
            water_axial
            - axial_high
        )
        <= END_ZONE_WIDTH_NM
    )

    near_tube_radial = (
        water_radial
        <= (
            wall_radius_q99
            + 0.30
        )
    )

    lower_inside = (
        (
            water_axial
            >= axial_low
        )
        & (
            water_axial
            <= axial_low
            + END_ZONE_WIDTH_NM
        )
        & near_tube_radial
    )

    lower_outside = (
        (
            water_axial
            < axial_low
        )
        & (
            water_axial
            >= axial_low
            - END_ZONE_WIDTH_NM
        )
        & near_tube_radial
    )

    upper_inside = (
        (
            water_axial
            <= axial_high
        )
        & (
            water_axial
            >= axial_high
            - END_ZONE_WIDTH_NM
        )
        & near_tube_radial
    )

    upper_outside = (
        (
            water_axial
            > axial_high
        )
        & (
            water_axial
            <= axial_high
            + END_ZONE_WIDTH_NM
        )
        & near_tube_radial
    )

    lower_near_water = (
        water_oxygen_positions[
            lower_near_axial
            & near_tube_radial
        ]
    )

    upper_near_water = (
        water_oxygen_positions[
            upper_near_axial
            & near_tube_radial
        ]
    )

    lower_nearest = nearest_distances_to_subset(
        lower_near_water,
        lower_end_ring_positions,
        box_lengths,
    )

    upper_nearest = nearest_distances_to_subset(
        upper_near_water,
        upper_end_ring_positions,
        box_lengths,
    )

    timeseries_rows = read_csv_rows(
        authoritative_timeseries
    )

    frozen_rows = [
        row
        for row in timeseries_rows
        if row.get(
            "dataset",
            "",
        ).strip().lower()
        == "frozen"
    ]

    if len(
        frozen_rows
    ) != 201:
        raise RuntimeError(
            "Expected 201 frozen rows in the "
            "authoritative water time series; found "
            f"{len(frozen_rows)}"
        )

    frozen_rows.sort(
        key=lambda row: float(
            row[
                "relative_time_ps"
            ]
        )
    )

    t0_row = frozen_rows[
        0
    ]

    endpoint_row = frozen_rows[
        -1
    ]

    if not math.isclose(
        float(
            t0_row[
                "relative_time_ps"
            ]
        ),
        0.0,
        abs_tol=1.0e-9,
    ):
        raise RuntimeError(
            "The first frozen time-series row is "
            "not t=0 ps."
        )

    authoritative_occupancies = np.asarray(
        [
            float(
                row[
                    "lumen_water_count"
                ]
            )
            for row in frozen_rows
        ],
        dtype=float,
    )

    authoritative_t0_occupancy = float(
        t0_row[
            "lumen_water_count"
        ]
    )

    authoritative_endpoint_occupancy = float(
        endpoint_row[
            "lumen_water_count"
        ]
    )

    authoritative_max_occupancy = float(
        np.max(
            authoritative_occupancies
        )
    )

    t0_fraction_of_maximum = (
        authoritative_t0_occupancy
        / authoritative_max_occupancy
        if authoritative_max_occupancy > 0.0
        else math.nan
    )

    authoritative_density = float(
        t0_row[
            "lumen_number_density_nm-3"
        ]
    )

    authoritative_wall_radius = float(
        t0_row[
            "wall_radius_nm"
        ]
    )

    authoritative_lumen_length = float(
        t0_row[
            "lumen_length_nm"
        ]
    )

    composition_rows = []

    for segment_name, segment_atoms in (
        (
            "HBN_ordered_segment",
            hbn_atoms,
        ),
        (
            "PYR_ordered_segment",
            pyr_atoms,
        ),
        (
            "SOL_ordered_segment",
            water_atoms,
        ),
    ):
        residue_counts = Counter(
            atom[
                "resname"
            ]
            for atom in segment_atoms
        )

        atom_name_counts = Counter(
            atom[
                "atomname"
            ]
            for atom in segment_atoms
        )

        composition_rows.append(
            {
                "segment": segment_name,
                "atom_count": len(
                    segment_atoms
                ),
                "residue_names": " | ".join(
                    f"{name}:{count}"
                    for name, count
                    in sorted(
                        residue_counts.items()
                    )
                ),
                "atom_names": " | ".join(
                    f"{name}:{count}"
                    for name, count
                    in sorted(
                        atom_name_counts.items()
                    )
                ),
            }
        )

    with COMPOSITION_CSV.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                composition_rows[
                    0
                ].keys()
            ),
        )

        writer.writeheader()
        writer.writerows(
            composition_rows
        )

    end_rows = [
        {
            "end": "lower",
            "robust_plane_coordinate_nm": axial_low,
            "end_ring_atom_count": int(
                np.count_nonzero(
                    lower_end_ring_mask
                )
            ),
            "waterO_inside_0p30nm": int(
                np.count_nonzero(
                    lower_inside
                )
            ),
            "waterO_outside_0p30nm": int(
                np.count_nonzero(
                    lower_outside
                )
            ),
            "waterO_total_near_plane_0p30nm": int(
                np.count_nonzero(
                    lower_near_axial
                    & near_tube_radial
                )
            ),
            "minimum_waterO_to_end_ring_nm": (
                float(
                    np.min(
                        lower_nearest
                    )
                )
                if len(
                    lower_nearest
                )
                else math.nan
            ),
        },
        {
            "end": "upper",
            "robust_plane_coordinate_nm": axial_high,
            "end_ring_atom_count": int(
                np.count_nonzero(
                    upper_end_ring_mask
                )
            ),
            "waterO_inside_0p30nm": int(
                np.count_nonzero(
                    upper_inside
                )
            ),
            "waterO_outside_0p30nm": int(
                np.count_nonzero(
                    upper_outside
                )
            ),
            "waterO_total_near_plane_0p30nm": int(
                np.count_nonzero(
                    upper_near_axial
                    & near_tube_radial
                )
            ),
            "minimum_waterO_to_end_ring_nm": (
                float(
                    np.min(
                        upper_nearest
                    )
                )
                if len(
                    upper_nearest
                )
                else math.nan
            ),
        },
    ]

    with END_ZONE_CSV.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                end_rows[
                    0
                ].keys()
            ),
        )

        writer.writeheader()
        writer.writerows(
            end_rows
        )

    gates = {
        "extracted_frame_has_68320_atoms": (
            len(
                t0_atoms
            )
            == EXPECTED_ATOMS
        ),
        "ordered_HBN_count_is_1680": (
            len(
                hbn_atoms
            )
            == HBN_ATOMS
        ),
        "ordered_PYR_count_is_104": (
            len(
                pyr_atoms
            )
            == PYR_ATOMS
        ),
        "ordered_water_atom_count_is_66536": (
            len(
                water_atoms
            )
            == WATER_ATOMS
        ),
        "TIP4P_water_count_is_16634": (
            len(
                water_oxygen_positions
            )
            == EXPECTED_WATERS
        ),
        "water_chunks_are_residue_consistent": (
            water_chunks_consistent
        ),
        "solute_matches_accepted_frozen_GRO": (
            solute_max_difference
            <= SOLUTE_MATCH_TOLERANCE_NM
        ),
        "box_matches_accepted_frozen_GRO": (
            box_difference
            <= 1.0e-6
        ),
        "authoritative_t0_occupancy_is_positive": (
            authoritative_t0_occupancy
            > 0.0
        ),
        "authoritative_t0_occupancy_at_least_100": (
            authoritative_t0_occupancy
            >= 100.0
        ),
        "authoritative_t0_is_at_least_90pct_of_max": (
            t0_fraction_of_maximum
            >= 0.90
        ),
        "tube_length_is_positive": (
            tube_length_p98
            > 0.0
        ),
        "wall_radius_is_positive": (
            wall_radius_median
            > 0.0
        ),
    }

    failed_gates = [
        name
        for name, passed
        in gates.items()
        if not passed
    ]

    authoritative_start_accepted = (
        len(
            failed_gates
        )
        == 0
    )

    decision = (
        "PASS"
        if authoritative_start_accepted
        else "REVIEW"
    )

    required_next_step = (
        "DEFINE_AND_GENERATE_R1_STERIC_CAP_PROTOTYPE"
        if authoritative_start_accepted
        else
        "RESOLVE_R0_T0_REFERENCE_GATE_FAILURES"
    )

    summary = {
        "decision": decision,
        "authoritative_R1_start_state_accepted": (
            authoritative_start_accepted
        ),
        "accepted_xtc": relative(
            accepted_xtc
        ),
        "accepted_tpr": relative(
            accepted_tpr
        ),
        "accepted_gro": relative(
            accepted_gro
        ),
        "authoritative_water_timeseries": (
            relative(
                authoritative_timeseries
            )
        ),
        "extracted_t0_gro": relative(
            T0_GRO
        ),
        "extracted_t0_sha256": sha256_file(
            T0_GRO
        ),
        "extracted_title": t0_title,
        "extracted_time_ps": (
            t0_time
            if t0_time is not None
            else 0.0
        ),
        "accepted_gro_title": (
            accepted_gro_title
        ),
        "accepted_gro_time_ps": (
            accepted_gro_time
            if accepted_gro_time is not None
            else ""
        ),
        "atom_count": len(
            t0_atoms
        ),
        "HBN_atom_count": len(
            hbn_atoms
        ),
        "PYR_atom_count": len(
            pyr_atoms
        ),
        "water_atom_count": len(
            water_atoms
        ),
        "water_molecule_count": len(
            water_oxygen_positions
        ),
        "water_chunk_consistency": (
            water_chunks_consistent
        ),
        "first_water_site_names": (
            " | ".join(
                f"{name}:{count}"
                for name, count
                in sorted(
                    first_water_site_names.items()
                )
            )
        ),
        "box_x_nm": box_lengths[
            0
        ],
        "box_y_nm": box_lengths[
            1
        ],
        "box_z_nm": box_lengths[
            2
        ],
        "solute_RMS_difference_vs_accepted_GRO_nm": (
            solute_rms_difference
        ),
        "solute_max_difference_vs_accepted_GRO_nm": (
            solute_max_difference
        ),
        "box_max_difference_vs_accepted_GRO_nm": (
            box_difference
        ),
        "tube_axis_x": tube_axis[
            0
        ],
        "tube_axis_y": tube_axis[
            1
        ],
        "tube_axis_z": tube_axis[
            2
        ],
        "tube_PCA_eigenvalue_ratio_1_to_2": (
            eigenvalue_ratio
        ),
        "robust_axial_low_nm": (
            axial_low
        ),
        "robust_axial_high_nm": (
            axial_high
        ),
        "tube_length_p98_nm": (
            tube_length_p98
        ),
        "wall_radius_mean_nm": (
            wall_radius_mean
        ),
        "wall_radius_median_nm": (
            wall_radius_median
        ),
        "wall_radius_q01_nm": (
            wall_radius_q01
        ),
        "wall_radius_q99_nm": (
            wall_radius_q99
        ),
        "provisional_accessible_radius_nm": (
            accessible_radius_provisional
        ),
        "provisional_accessible_volume_nm3": (
            accessible_volume_provisional
        ),
        "independent_geometric_lumen_water_count": (
            independent_lumen_count
        ),
        "authoritative_t0_lumen_water_count": (
            authoritative_t0_occupancy
        ),
        "authoritative_max_lumen_water_count": (
            authoritative_max_occupancy
        ),
        "authoritative_endpoint_lumen_water_count": (
            authoritative_endpoint_occupancy
        ),
        "authoritative_t0_fraction_of_maximum": (
            t0_fraction_of_maximum
        ),
        "authoritative_t0_density_nm-3": (
            authoritative_density
        ),
        "authoritative_t0_wall_radius_nm": (
            authoritative_wall_radius
        ),
        "authoritative_t0_lumen_length_nm": (
            authoritative_lumen_length
        ),
        "lower_end_ring_atom_count": int(
            np.count_nonzero(
                lower_end_ring_mask
            )
        ),
        "upper_end_ring_atom_count": int(
            np.count_nonzero(
                upper_end_ring_mask
            )
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

    write_single_row_csv(
        SUMMARY_CSV,
        summary,
    )

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed
        in gates.items()
    )

    REPORT_MD.write_text(
        f"""# R0 t=0 Hydrated Reference Audit

## Authoritative sources

- Accepted trajectory:
  `{relative(accepted_xtc)}`
- Accepted run input:
  `{relative(accepted_tpr)}`
- Accepted GRO:
  `{relative(accepted_gro)}`
- Validated hydration time series:
  `{relative(authoritative_timeseries)}`

## Extracted R1 starting state

- Extracted GRO:
  `{relative(T0_GRO)}`
- SHA256:
  `{summary['extracted_t0_sha256']}`
- Parsed time:
  **{summary['extracted_time_ps']} ps**
- Atoms:
  **{summary['atom_count']}**

## Composition

- HBN atoms: **{summary['HBN_atom_count']}**
- PYR atoms: **{summary['PYR_atom_count']}**
- Water atoms: **{summary['water_atom_count']}**
- TIP4P/2005 waters: **{summary['water_molecule_count']}**
- Water-residue chunk consistency:
  **{summary['water_chunk_consistency']}**
- First water-site names:
  `{summary['first_water_site_names']}`

## Consistency with accepted R0

- Solute RMS difference:
  **{solute_rms_difference:.12f} nm**
- Solute maximum difference:
  **{solute_max_difference:.12f} nm**
- Box maximum difference:
  **{box_difference:.12f} nm**

The accepted GRO may represent the final frozen trajectory state.
Solute identity is still expected because HBN and PYR were frozen.

## Tube geometry

- PCA axis:
  **({tube_axis[0]:.8f},
  {tube_axis[1]:.8f},
  {tube_axis[2]:.8f})**
- Robust axial planes:
  **{axial_low:.6f}/{axial_high:.6f} nm**
- p98 tube length:
  **{tube_length_p98:.6f} nm**
- Wall radius mean/median:
  **{wall_radius_mean:.6f}/{wall_radius_median:.6f} nm**
- Wall-radius q01/q99:
  **{wall_radius_q01:.6f}/{wall_radius_q99:.6f} nm**
- Provisional accessible radius:
  **{accessible_radius_provisional:.6f} nm**
- Provisional accessible volume:
  **{accessible_volume_provisional:.6f} nm³**

## Authoritative hydration state

- t=0 lumen occupancy:
  **{authoritative_t0_occupancy:.0f} waters**
- Maximum accepted-trajectory occupancy:
  **{authoritative_max_occupancy:.0f} waters**
- t=0 fraction of maximum:
  **{t0_fraction_of_maximum:.6f}**
- Accepted endpoint occupancy:
  **{authoritative_endpoint_occupancy:.0f} waters**
- t=0 lumen density:
  **{authoritative_density:.6f} nm^-3**
- Independent provisional geometric count:
  **{independent_lumen_count} waters**

The validated time-series occupancy remains authoritative.
The independent geometric value is retained only as a cross-check.

## End-zone audit

### Lower end

- Ring atoms: **{end_rows[0]['end_ring_atom_count']}**
- Water oxygens inside 0.30 nm:
  **{end_rows[0]['waterO_inside_0p30nm']}**
- Water oxygens outside 0.30 nm:
  **{end_rows[0]['waterO_outside_0p30nm']}**
- Minimum water-O/end-ring distance:
  **{end_rows[0]['minimum_waterO_to_end_ring_nm']:.6f} nm**

### Upper end

- Ring atoms: **{end_rows[1]['end_ring_atom_count']}**
- Water oxygens inside 0.30 nm:
  **{end_rows[1]['waterO_inside_0p30nm']}**
- Water oxygens outside 0.30 nm:
  **{end_rows[1]['waterO_outside_0p30nm']}**
- Minimum water-O/end-ring distance:
  **{end_rows[1]['minimum_waterO_to_end_ring_nm']:.6f} nm**

## Gates

{gate_lines}

## Decision

- Audit decision: **{decision}**
- Authoritative R1 start accepted:
  **{'YES' if authoritative_start_accepted else 'NO'}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Required next step:
  `{required_next_step}`

No MD, topology modification, cap generation, or QM calculation was
performed by this audit.
""",
        encoding="utf-8",
    )

    CONTRACT_UPDATE_MD.write_text(
        f"""# R1 Authoritative Starting State

The authoritative starting state for the R1 fully capped positive
control is:

`{relative(T0_GRO)}`

SHA256:

`{summary['extracted_t0_sha256']}`

Source trajectory:

`{relative(accepted_xtc)}`

Source time:

**0.0 ps**

Accepted starting hydration:

- lumen occupancy:
  **{authoritative_t0_occupancy:.0f} waters**
- maximum occupancy observed in the accepted R0 trajectory:
  **{authoritative_max_occupancy:.0f} waters**
- t=0 fraction of trajectory maximum:
  **{t0_fraction_of_maximum:.6f}**
- endpoint occupancy:
  **{authoritative_endpoint_occupancy:.0f} waters**

This state preserves the accepted R0 atom ordering:

- HBN: atoms 1-{HBN_ATOMS}
- PYR: atoms {HBN_ATOMS + 1}-{SOLUTE_ATOMS}
- TIP4P/2005 water:
  atoms {SOLUTE_ATOMS + 1}-{EXPECTED_ATOMS}

R1 cap construction must not modify this file in place. A derived
structure must be created in a new R1 design directory and must retain
a provenance link to this SHA256.
""",
        encoding="utf-8",
    )

    print(
        "Day023 R0 t=0 hydrated-state audit completed."
    )

    print(
        "Accepted R0 XTC/TPR/GRO: "
        f"{accepted_xtc.name} / "
        f"{accepted_tpr.name} / "
        f"{accepted_gro.name}"
    )

    print(
        "Extracted t=0 atoms / time: "
        f"{len(t0_atoms)}/{EXPECTED_ATOMS} / "
        f"{0.0 if t0_time is None else t0_time:.3f} ps"
    )

    print(
        "HBN/PYR/water atoms/waters: "
        f"{len(hbn_atoms)}/"
        f"{len(pyr_atoms)}/"
        f"{len(water_atoms)}/"
        f"{len(water_oxygen_positions)}"
    )

    print(
        "Water residue-chunk consistency: "
        f"{'PASS' if water_chunks_consistent else 'FAIL'}"
    )

    print(
        "Solute RMS/max difference vs accepted GRO: "
        f"{solute_rms_difference:.12f}/"
        f"{solute_max_difference:.12f} nm"
    )

    print(
        "Tube axis: "
        f"{tube_axis[0]:.8f} "
        f"{tube_axis[1]:.8f} "
        f"{tube_axis[2]:.8f}"
    )

    print(
        "Tube robust axial planes / p98 length: "
        f"{axial_low:.6f}/"
        f"{axial_high:.6f} / "
        f"{tube_length_p98:.6f} nm"
    )

    print(
        "Wall radius mean/median/q01/q99: "
        f"{wall_radius_mean:.6f}/"
        f"{wall_radius_median:.6f}/"
        f"{wall_radius_q01:.6f}/"
        f"{wall_radius_q99:.6f} nm"
    )

    print(
        "Authoritative t0/max/endpoint lumen occupancy: "
        f"{authoritative_t0_occupancy:.0f}/"
        f"{authoritative_max_occupancy:.0f}/"
        f"{authoritative_endpoint_occupancy:.0f}"
    )

    print(
        "Authoritative t0 fraction of maximum: "
        f"{t0_fraction_of_maximum:.6f}"
    )

    print(
        "Independent provisional lumen count: "
        f"{independent_lumen_count}"
    )

    print(
        "Lower end ring / inside water / outside water / "
        "minimum distance: "
        f"{end_rows[0]['end_ring_atom_count']} / "
        f"{end_rows[0]['waterO_inside_0p30nm']} / "
        f"{end_rows[0]['waterO_outside_0p30nm']} / "
        f"{end_rows[0]['minimum_waterO_to_end_ring_nm']:.6f} nm"
    )

    print(
        "Upper end ring / inside water / outside water / "
        "minimum distance: "
        f"{end_rows[1]['end_ring_atom_count']} / "
        f"{end_rows[1]['waterO_inside_0p30nm']} / "
        f"{end_rows[1]['waterO_outside_0p30nm']} / "
        f"{end_rows[1]['minimum_waterO_to_end_ring_nm']:.6f} nm"
    )

    print(
        f"Audit decision: {decision}"
    )

    print(
        "Authoritative R1 start accepted: "
        f"{'YES' if authoritative_start_accepted else 'NO'}"
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
        f"Required next step: {required_next_step}"
    )

    print(
        f"Wrote: {relative(T0_GRO)}"
    )

    print(
        f"Wrote: {relative(SUMMARY_CSV)}"
    )

    print(
        f"Wrote: {relative(END_ZONE_CSV)}"
    )

    print(
        f"Wrote: {relative(COMPOSITION_CSV)}"
    )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )

    print(
        f"Wrote: {relative(CONTRACT_UPDATE_MD)}"
    )

    if failed_gates:
        raise RuntimeError(
            "The R0 t=0 reference requires review "
            "before R1 cap generation."
        )


if __name__ == "__main__":
    main()
