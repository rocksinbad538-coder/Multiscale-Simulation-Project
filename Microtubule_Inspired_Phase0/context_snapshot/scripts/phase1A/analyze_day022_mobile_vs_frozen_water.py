#!/usr/bin/env python3

from __future__ import annotations

import csv
import importlib.util
import math
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

PROTOCOL = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

MOBILE_STAGE = "08_nvt_mobile_100ps"

MOBILE_RUN = (
    PROTOCOL
    / "execution"
    / MOBILE_STAGE
)

MOBILE_FRAME_ROOT = (
    MOBILE_RUN
    / "time_resolved_stability"
    / "solute_waterO_frames"
)

FROZEN_RUN = (
    ROOT
    / "runs/phase1A/accepted/"
    "hybrid_hydrated_gap45_pyr5shift_clean032_"
    "nvt_100ps_frozenSolute"
)

HBN_ITP = (
    PROTOCOL
    / "protocol_inputs/topology/"
    "hbn_bonded_mobile_release.itp"
)

PYR_ITP = (
    PROTOCOL
    / "protocol_inputs/topology/"
    "pyrene_mobile_release.itp"
)

INTEGRATED_MODULE_PATH = (
    ROOT
    / "scripts/phase1A/"
    "validate_day022_stage08_mobile_100ps.py"
)

ANALYSIS_ROOT = (
    MOBILE_RUN
    / "mobile_vs_frozen_water"
)

FROZEN_FRAME_ROOT = (
    ANALYSIS_ROOT
    / "frozen_solute_waterO_frames"
)

INDEX_FILE = (
    ANALYSIS_ROOT
    / "hbn_pyr_waterO.ndx"
)

FROZEN_EXTRACTION_LOG = (
    ANALYSIS_ROOT
    / "frozen_water_frame_extraction.log"
)

TIMESERIES_CSV = (
    ANALYSIS_ROOT
    / "mobile_frozen_water_timeseries.csv"
)

RADIAL_CSV = (
    ANALYSIS_ROOT
    / "mobile_frozen_radial_density.csv"
)

AXIAL_CSV = (
    ANALYSIS_ROOT
    / "mobile_frozen_axial_density.csv"
)

PYR_HYDRATION_CSV = (
    ANALYSIS_ROOT
    / "mobile_frozen_pyrene_hydration.csv"
)

SNAPSHOT_CSV = (
    ANALYSIS_ROOT
    / "mobile_frozen_representative_snapshot_pairs.csv"
)

SUMMARY_CSV = (
    ANALYSIS_ROOT
    / "mobile_frozen_water_comparison_summary.csv"
)

REPORT_MD = (
    ANALYSIS_ROOT
    / "MOBILE_VS_FROZEN_WATER_COMPARISON_DAY022.md"
)

HBN_COUNT = 1680
PYR_COUNT = 4
PYR_ATOMS = 26
SOLUTE_COUNT = HBN_COUNT + PYR_COUNT * PYR_ATOMS

TOTAL_ATOMS = 68320
FIRST_WATER_ATOM = SOLUTE_COUNT + 1
WATER_ATOMS_PER_MOLECULE = 4
WATER_COUNT = (
    TOTAL_ATOMS - SOLUTE_COUNT
) // WATER_ATOMS_PER_MOLECULE

SELECTED_ATOM_COUNT = (
    SOLUTE_COUNT + WATER_COUNT
)

EXPECTED_FRAMES = 201
EXPECTED_DURATION_PS = 100.0

RADIAL_BIN_COUNT = 24
AXIAL_BIN_COUNT = 24

HYDRATION_CUTOFFS_NM = (
    0.35,
    0.50,
)

REPRESENTATIVE_PAIR_COUNT = 21


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


def require_directory(path: Path) -> None:
    if (
        not path.exists()
        or not path.is_dir()
    ):
        raise RuntimeError(
            f"Missing directory: {path}"
        )


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
        "Could not locate GROMACS"
    )


def load_module(
    path: Path,
    name: str,
):
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


def locate_unique_file(
    directory: Path,
    preferred_name: str,
    suffix: str,
) -> Path:
    preferred = directory / preferred_name

    if (
        preferred.exists()
        and preferred.stat().st_size > 0
    ):
        return preferred

    matches = sorted(
        path
        for path in directory.glob(
            f"*{suffix}"
        )
        if (
            path.is_file()
            and path.stat().st_size > 0
        )
    )

    if len(matches) != 1:
        raise RuntimeError(
            f"Could not identify one {suffix} file "
            f"in {directory}. Found: "
            + ", ".join(
                path.name
                for path in matches
            )
        )

    return matches[0]


def natural_key(path: Path) -> int:
    match = re.search(
        r"(\d+)(?=\.gro$)",
        path.name,
    )

    if match is None:
        return -1

    return int(
        match.group(1)
    )


def list_mobile_frames() -> list[Path]:
    frames = sorted(
        MOBILE_FRAME_ROOT.glob(
            "stage08_frame*.gro"
        ),
        key=natural_key,
    )

    if len(frames) != EXPECTED_FRAMES:
        raise RuntimeError(
            "Expected 201 mobile frames; "
            f"found {len(frames)}"
        )

    return frames


def write_index() -> None:
    solute_indices = list(
        range(
            1,
            SOLUTE_COUNT + 1,
        )
    )

    water_oxygen_indices = list(
        range(
            FIRST_WATER_ATOM,
            TOTAL_ATOMS + 1,
            WATER_ATOMS_PER_MOLECULE,
        )
    )

    if len(
        water_oxygen_indices
    ) != WATER_COUNT:
        raise RuntimeError(
            "Unexpected water-oxygen index count"
        )

    selected = (
        solute_indices
        + water_oxygen_indices
    )

    with INDEX_FILE.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "[ HBN_PYR_WATERO ]\n"
        )

        for start in range(
            0,
            len(selected),
            15,
        ):
            handle.write(
                " ".join(
                    str(value)
                    for value in selected[
                        start : start + 15
                    ]
                )
                + "\n"
            )


def extract_frozen_frames(
    gmx: str,
    xtc: Path,
    tpr: Path,
) -> list[Path]:
    FROZEN_FRAME_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for old_frame in FROZEN_FRAME_ROOT.glob(
        "frozen_frame*.gro"
    ):
        old_frame.unlink()

    output_template = (
        FROZEN_FRAME_ROOT
        / "frozen_frame.gro"
    )

    command = [
        gmx,
        "trjconv",
        "-f",
        str(xtc),
        "-s",
        str(tpr),
        "-n",
        str(INDEX_FILE),
        "-o",
        str(output_template),
        "-sep",
        "-pbc",
        "atom",
    ]

    completed = subprocess.run(
        command,
        input="HBN_PYR_WATERO\n",
        cwd=FROZEN_RUN,
        env={
            **os.environ,
            "GMX_MAXBACKUP": "-1",
        },
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )

    FROZEN_EXTRACTION_LOG.write_text(
        completed.stdout,
        encoding="utf-8",
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "Frozen frame extraction failed. "
            f"See {FROZEN_EXTRACTION_LOG}"
        )

    frames = sorted(
        FROZEN_FRAME_ROOT.glob(
            "frozen_frame*.gro"
        ),
        key=natural_key,
    )

    if len(frames) != EXPECTED_FRAMES:
        raise RuntimeError(
            "Expected 201 frozen frames; "
            f"found {len(frames)}"
        )

    return frames


def read_gro(
    path: Path,
) -> tuple[
    str,
    list[dict[str, object]],
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

    try:
        natoms = int(
            lines[1].strip()
        )
    except ValueError as error:
        raise RuntimeError(
            f"Invalid GRO atom count: {path}"
        ) from error

    if natoms != SELECTED_ATOM_COUNT:
        raise RuntimeError(
            f"Expected {SELECTED_ATOM_COUNT} atoms "
            f"in {path}; found {natoms}"
        )

    atoms = []

    for local_index, line in enumerate(
        lines[
            2 : 2 + natoms
        ],
        start=1,
    ):
        atoms.append(
            {
                "local_index": local_index,
                "residue_name": (
                    line[5:10].strip()
                ),
                "atom_name": (
                    line[10:15].strip()
                ),
                "position": np.array(
                    [
                        float(line[20:28]),
                        float(line[28:36]),
                        float(line[36:44]),
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

    if len(box_values) not in {
        3,
        9,
    }:
        raise RuntimeError(
            f"Unexpected GRO box format: {path}"
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
                "This workflow currently requires "
                "an orthorhombic simulation box"
            )

    box = np.array(
        box_values[:3],
        dtype=float,
    )

    if np.any(
        box <= 0.0
    ):
        raise RuntimeError(
            f"Invalid box dimensions in {path}"
        )

    return (
        title,
        atoms,
        box,
    )


def parse_time(
    title: str,
    fallback: float,
) -> float:
    match = re.search(
        r"\bt\s*=\s*([-+0-9.eE]+)",
        title,
    )

    if match is None:
        return fallback

    return float(
        match.group(1)
    )


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


def kabsch_transform(
    reference: np.ndarray,
    target: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    reference_center = (
        reference.mean(axis=0)
    )

    target_center = (
        target.mean(axis=0)
    )

    reference_centered = (
        reference
        - reference_center
    )

    target_centered = (
        target
        - target_center
    )

    covariance = (
        target_centered.T
        @ reference_centered
    )

    u_matrix, _, vt_matrix = (
        np.linalg.svd(
            covariance
        )
    )

    correction = np.eye(3)

    correction[2, 2] = np.sign(
        np.linalg.det(
            u_matrix
            @ vt_matrix
        )
    )

    rotation = (
        u_matrix
        @ correction
        @ vt_matrix
    )

    aligned = (
        target_centered
        @ rotation
        + reference_center
    )

    return (
        aligned,
        rotation,
        target_center,
        reference_center,
    )


def apply_transform(
    coordinates: np.ndarray,
    rotation: np.ndarray,
    target_center: np.ndarray,
    reference_center: np.ndarray,
) -> np.ndarray:
    return (
        (
            coordinates
            - target_center
        )
        @ rotation
        + reference_center
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


def transverse_basis(
    axis: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:
    trial = np.array(
        [1.0, 0.0, 0.0],
        dtype=float,
    )

    if abs(
        float(
            np.dot(
                trial,
                axis,
            )
        )
    ) > 0.90:
        trial = np.array(
            [0.0, 1.0, 0.0],
            dtype=float,
        )

    basis_1 = (
        trial
        - np.dot(
            trial,
            axis,
        )
        * axis
    )

    basis_1 = (
        basis_1
        / np.linalg.norm(
            basis_1
        )
    )

    basis_2 = np.cross(
        axis,
        basis_1,
    )

    basis_2 = (
        basis_2
        / np.linalg.norm(
            basis_2
        )
    )

    return (
        basis_1,
        basis_2,
    )


def quantile(
    values: np.ndarray,
    probability: float,
) -> float:
    return float(
        np.quantile(
            values,
            probability,
        )
    )


def linear_slope(
    times: np.ndarray,
    values: np.ndarray,
) -> float:
    if (
        len(times) < 2
        or np.allclose(
            times,
            times[0],
        )
    ):
        return 0.0

    return float(
        np.polyfit(
            times,
            values,
            1,
        )[0]
    )


def percent_difference(
    mobile: float,
    frozen: float,
) -> float:
    if abs(frozen) < 1.0e-15:
        return math.nan

    return (
        100.0
        * (
            mobile
            - frozen
        )
        / frozen
    )


def js_divergence(
    first: np.ndarray,
    second: np.ndarray,
) -> float:
    first = np.asarray(
        first,
        dtype=float,
    )

    second = np.asarray(
        second,
        dtype=float,
    )

    if (
        first.sum() <= 0.0
        or second.sum() <= 0.0
    ):
        return math.nan

    first = first / first.sum()
    second = second / second.sum()

    midpoint = 0.5 * (
        first + second
    )

    first_term = np.where(
        first > 0.0,
        first
        * np.log(
            first / midpoint
        ),
        0.0,
    )

    second_term = np.where(
        second > 0.0,
        second
        * np.log(
            second / midpoint
        ),
        0.0,
    )

    return float(
        0.5
        * (
            first_term.sum()
            + second_term.sum()
        )
    )


def hydration_count(
    water_positions: np.ndarray,
    pyrene_positions: np.ndarray,
    cutoff_nm: float,
) -> int:
    delta = (
        water_positions[:, None, :]
        - pyrene_positions[None, :, :]
    )

    minimum_squared_distance = np.min(
        np.sum(
            delta ** 2,
            axis=2,
        ),
        axis=1,
    )

    return int(
        np.count_nonzero(
            minimum_squared_distance
            <= cutoff_nm ** 2
        )
    )


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
        for key in row:
            if key in seen:
                continue

            seen.add(key)
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


def process_dataset(
    label: str,
    frame_paths: list[Path],
    integrated,
    hbn_topology,
    pyr_topology,
    reference_hbn: np.ndarray,
    reference_axis: np.ndarray,
    basis_1: np.ndarray,
    basis_2: np.ndarray,
) -> dict[str, object]:
    radial_edges = np.linspace(
        0.0,
        1.0,
        RADIAL_BIN_COUNT + 1,
    )

    axial_edges = np.linspace(
        0.0,
        1.0,
        AXIAL_BIN_COUNT + 1,
    )

    radial_density_frames = []
    axial_density_frames = []

    radial_count_total = np.zeros(
        RADIAL_BIN_COUNT,
        dtype=float,
    )

    axial_count_total = np.zeros(
        AXIAL_BIN_COUNT,
        dtype=float,
    )

    rows = []
    hydration_rows = []

    first_metadata = None

    for frame_index, frame_path in enumerate(
        frame_paths
    ):
        title, atoms, box = read_gro(
            frame_path
        )

        if first_metadata is None:
            first_metadata = atoms

        time_ps = parse_time(
            title,
            fallback=0.5 * frame_index,
        )

        positions = np.array(
            [
                atom["position"]
                for atom in atoms
            ],
            dtype=float,
        )

        current_hbn_wrapped = (
            positions[:HBN_COUNT]
        )

        current_hbn = (
            integrated.unwrap_by_bonds(
                current_hbn_wrapped,
                box,
                hbn_topology["bonds"],
            )
        )

        (
            aligned_hbn,
            rotation,
            target_center,
            reference_center,
        ) = kabsch_transform(
            reference_hbn,
            current_hbn,
        )

        hbn_center_current = (
            current_hbn.mean(axis=0)
        )

        hbn_center_wrapped = np.mod(
            hbn_center_current,
            box,
        )

        water_wrapped = positions[
            SOLUTE_COUNT:
        ]

        water_near_hbn = (
            hbn_center_current
            + minimum_image(
                water_wrapped
                - hbn_center_wrapped,
                box,
            )
        )

        aligned_water = apply_transform(
            water_near_hbn,
            rotation,
            target_center,
            reference_center,
        )

        aligned_center = (
            aligned_hbn.mean(axis=0)
        )

        hbn_relative = (
            aligned_hbn
            - aligned_center
        )

        hbn_axial = (
            hbn_relative
            @ reference_axis
        )

        hbn_transverse_1 = (
            hbn_relative
            @ basis_1
        )

        hbn_transverse_2 = (
            hbn_relative
            @ basis_2
        )

        hbn_radial = np.sqrt(
            hbn_transverse_1 ** 2
            + hbn_transverse_2 ** 2
        )

        axial_lower = quantile(
            hbn_axial,
            0.01,
        )

        axial_upper = quantile(
            hbn_axial,
            0.99,
        )

        central_lower = quantile(
            hbn_axial,
            0.10,
        )

        central_upper = quantile(
            hbn_axial,
            0.90,
        )

        central_mask = (
            (hbn_axial >= central_lower)
            & (hbn_axial <= central_upper)
        )

        wall_radius = float(
            np.median(
                hbn_radial[
                    central_mask
                ]
            )
        )

        lumen_length = (
            axial_upper
            - axial_lower
        )

        if (
            wall_radius <= 0.0
            or lumen_length <= 0.0
        ):
            raise RuntimeError(
                f"Invalid lumen geometry in {frame_path}"
            )

        water_relative = (
            aligned_water
            - aligned_center
        )

        water_axial = (
            water_relative
            @ reference_axis
        )

        water_transverse_1 = (
            water_relative
            @ basis_1
        )

        water_transverse_2 = (
            water_relative
            @ basis_2
        )

        water_radial = np.sqrt(
            water_transverse_1 ** 2
            + water_transverse_2 ** 2
        )

        axial_inside = (
            (water_axial >= axial_lower)
            & (water_axial <= axial_upper)
        )

        lumen_mask = (
            axial_inside
            & (water_radial <= wall_radius)
        )

        occupancy = int(
            np.count_nonzero(
                lumen_mask
            )
        )

        lumen_volume = (
            math.pi
            * wall_radius ** 2
            * lumen_length
        )

        lumen_density = (
            occupancy
            / lumen_volume
        )

        normalized_radial = (
            water_radial[
                lumen_mask
            ]
            / wall_radius
        )

        normalized_axial = (
            (
                water_axial[
                    lumen_mask
                ]
                - axial_lower
            )
            / lumen_length
        )

        radial_counts, _ = np.histogram(
            normalized_radial,
            bins=radial_edges,
        )

        axial_counts, _ = np.histogram(
            normalized_axial,
            bins=axial_edges,
        )

        radial_count_total += (
            radial_counts
        )

        axial_count_total += (
            axial_counts
        )

        radial_volumes = (
            math.pi
            * wall_radius ** 2
            * (
                radial_edges[1:] ** 2
                - radial_edges[:-1] ** 2
            )
            * lumen_length
        )

        axial_volumes = (
            math.pi
            * wall_radius ** 2
            * lumen_length
            * np.diff(
                axial_edges
            )
        )

        radial_density_frames.append(
            radial_counts
            / radial_volumes
        )

        axial_density_frames.append(
            axial_counts
            / axial_volumes
        )

        pyrene_heavy_positions = []

        hydration_values: dict[
            tuple[int, float],
            int,
        ] = {}

        for molecule_index in range(
            PYR_COUNT
        ):
            start = (
                HBN_COUNT
                + molecule_index
                * PYR_ATOMS
            )

            stop = (
                start
                + PYR_ATOMS
            )

            current_pyr_wrapped = positions[
                start:stop
            ]

            current_pyr = (
                integrated.unwrap_by_bonds(
                    current_pyr_wrapped,
                    box,
                    pyr_topology["bonds"],
                )
            )

            pyr_center = (
                current_pyr.mean(axis=0)
            )

            pyr_center_wrapped = np.mod(
                pyr_center,
                box,
            )

            target_pyr_center = (
                hbn_center_current
                + minimum_image(
                    pyr_center_wrapped
                    - hbn_center_wrapped,
                    box,
                )
            )

            current_pyr = (
                current_pyr
                + target_pyr_center
                - pyr_center
            )

            heavy_local_indices = [
                local_index
                for local_index in range(
                    PYR_ATOMS
                )
                if not str(
                    atoms[
                        start + local_index
                    ]["atom_name"]
                ).upper().startswith("H")
            ]

            heavy_positions = current_pyr[
                heavy_local_indices
            ]

            pyrene_heavy_positions.append(
                heavy_positions
            )

            for cutoff_nm in (
                HYDRATION_CUTOFFS_NM
            ):
                hydration_values[
                    (
                        molecule_index,
                        cutoff_nm,
                    )
                ] = hydration_count(
                    water_near_hbn,
                    heavy_positions,
                    cutoff_nm,
                )

        row: dict[str, object] = {
            "dataset": label,
            "frame": frame_index,
            "time_ps": time_ps,
            "wall_radius_nm": (
                wall_radius
            ),
            "lumen_length_nm": (
                lumen_length
            ),
            "lumen_volume_nm3": (
                lumen_volume
            ),
            "lumen_water_count": (
                occupancy
            ),
            "lumen_number_density_nm-3": (
                lumen_density
            ),
            "mean_normalized_radial_position": (
                float(
                    normalized_radial.mean()
                )
                if occupancy > 0
                else math.nan
            ),
            "mean_normalized_axial_position": (
                float(
                    normalized_axial.mean()
                )
                if occupancy > 0
                else math.nan
            ),
        }

        for molecule_index in range(
            PYR_COUNT
        ):
            for cutoff_nm in (
                HYDRATION_CUTOFFS_NM
            ):
                cutoff_label = int(
                    round(
                        cutoff_nm * 100
                    )
                )

                value = hydration_values[
                    (
                        molecule_index,
                        cutoff_nm,
                    )
                ]

                row[
                    f"PYR{molecule_index + 1}_"
                    f"waterO_within_{cutoff_label}pm"
                ] = value

                hydration_rows.append(
                    {
                        "dataset": label,
                        "frame": frame_index,
                        "time_ps": time_ps,
                        "pyrene": (
                            molecule_index + 1
                        ),
                        "cutoff_nm": cutoff_nm,
                        "waterO_count": value,
                    }
                )

        rows.append(row)

    times = np.array(
        [
            float(
                row["time_ps"]
            )
            for row in rows
        ],
        dtype=float,
    )

    times -= times[0]

    if not math.isclose(
        float(
            times[-1]
        ),
        EXPECTED_DURATION_PS,
        rel_tol=0.0,
        abs_tol=1.0e-6,
    ):
        raise RuntimeError(
            f"{label} duration is not 100 ps: "
            f"{times[-1]:.6f}"
        )

    for row, relative_time in zip(
        rows,
        times,
    ):
        row["relative_time_ps"] = (
            float(relative_time)
        )

    return {
        "rows": rows,
        "hydration_rows": hydration_rows,
        "times": times,
        "radial_edges": radial_edges,
        "axial_edges": axial_edges,
        "radial_density_frames": np.array(
            radial_density_frames,
            dtype=float,
        ),
        "axial_density_frames": np.array(
            axial_density_frames,
            dtype=float,
        ),
        "radial_count_total": (
            radial_count_total
        ),
        "axial_count_total": (
            axial_count_total
        ),
    }


def representative_indices(
    rows: list[dict[str, object]],
    count: int,
) -> list[int]:
    feature_names = [
        "lumen_water_count",
        "lumen_number_density_nm-3",
        "mean_normalized_radial_position",
        "wall_radius_nm",
        "lumen_length_nm",
    ]

    for molecule_index in range(
        1,
        PYR_COUNT + 1,
    ):
        feature_names.append(
            f"PYR{molecule_index}_"
            "waterO_within_35pm"
        )

    features = np.array(
        [
            [
                float(
                    row[name]
                )
                for name in feature_names
            ]
            for row in rows
        ],
        dtype=float,
    )

    column_means = np.nanmean(
        features,
        axis=0,
    )

    nonfinite = ~np.isfinite(
        features
    )

    if np.any(nonfinite):
        features[
            nonfinite
        ] = np.take(
            column_means,
            np.where(
                nonfinite
            )[1],
        )

    means = features.mean(
        axis=0
    )

    standard_deviations = features.std(
        axis=0
    )

    standard_deviations[
        standard_deviations < 1.0e-12
    ] = 1.0

    standardized = (
        features - means
    ) / standard_deviations

    center_index = int(
        np.argmin(
            np.linalg.norm(
                standardized,
                axis=1,
            )
        )
    )

    selected = [
        0,
        len(rows) - 1,
        center_index,
    ]

    selected = list(
        dict.fromkeys(
            selected
        )
    )

    while len(selected) < count:
        selected_array = standardized[
            selected
        ]

        distances = np.linalg.norm(
            standardized[:, None, :]
            - selected_array[None, :, :],
            axis=2,
        )

        minimum_distances = distances.min(
            axis=1
        )

        minimum_distances[
            selected
        ] = -1.0

        next_index = int(
            np.argmax(
                minimum_distances
            )
        )

        selected.append(
            next_index
        )

    return sorted(
        selected[:count]
    )


def main() -> None:
    require_directory(
        MOBILE_FRAME_ROOT
    )

    require_directory(
        FROZEN_RUN
    )

    require_file(
        HBN_ITP
    )

    require_file(
        PYR_ITP
    )

    require_file(
        INTEGRATED_MODULE_PATH
    )

    ANALYSIS_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    integrated = load_module(
        INTEGRATED_MODULE_PATH,
        "stage08_water_integrated_module",
    )

    gmx = locate_gmx()

    mobile_frames = list_mobile_frames()

    frozen_xtc = locate_unique_file(
        FROZEN_RUN,
        "nvt_100ps_frozenSolute.xtc",
        ".xtc",
    )

    frozen_tpr = locate_unique_file(
        FROZEN_RUN,
        "nvt_100ps_frozenSolute.tpr",
        ".tpr",
    )

    write_index()

    frozen_frames = extract_frozen_frames(
        gmx,
        frozen_xtc,
        frozen_tpr,
    )

    hbn_topology = integrated.parse_itp(
        HBN_ITP
    )

    pyr_topology = integrated.parse_itp(
        PYR_ITP
    )

    (
        _,
        frozen_reference_atoms,
        frozen_reference_box,
    ) = read_gro(
        frozen_frames[0]
    )

    frozen_reference_positions = np.array(
        [
            atom["position"]
            for atom in frozen_reference_atoms
        ],
        dtype=float,
    )

    reference_hbn = (
        integrated.unwrap_by_bonds(
            frozen_reference_positions[
                :HBN_COUNT
            ],
            frozen_reference_box,
            hbn_topology["bonds"],
        )
    )

    reference_axis = tube_axis(
        reference_hbn
    )

    basis_1, basis_2 = transverse_basis(
        reference_axis
    )

    frozen = process_dataset(
        "frozen",
        frozen_frames,
        integrated,
        hbn_topology,
        pyr_topology,
        reference_hbn,
        reference_axis,
        basis_1,
        basis_2,
    )

    mobile = process_dataset(
        "mobile",
        mobile_frames,
        integrated,
        hbn_topology,
        pyr_topology,
        reference_hbn,
        reference_axis,
        basis_1,
        basis_2,
    )

    if not np.allclose(
        mobile["times"],
        frozen["times"],
        atol=1.0e-6,
        rtol=0.0,
    ):
        raise RuntimeError(
            "Mobile and frozen frame times do not match"
        )

    timeseries_rows = (
        frozen["rows"]
        + mobile["rows"]
    )

    write_csv(
        TIMESERIES_CSV,
        timeseries_rows,
    )

    write_csv(
        PYR_HYDRATION_CSV,
        (
            frozen["hydration_rows"]
            + mobile["hydration_rows"]
        ),
    )

    radial_rows = []

    for dataset_label, dataset in (
        ("frozen", frozen),
        ("mobile", mobile),
    ):
        density_frames = dataset[
            "radial_density_frames"
        ]

        edges = dataset[
            "radial_edges"
        ]

        for bin_index in range(
            RADIAL_BIN_COUNT
        ):
            radial_rows.append(
                {
                    "dataset": dataset_label,
                    "bin": bin_index,
                    "normalized_radius_lower": (
                        edges[bin_index]
                    ),
                    "normalized_radius_upper": (
                        edges[
                            bin_index + 1
                        ]
                    ),
                    "normalized_radius_center": (
                        0.5
                        * (
                            edges[bin_index]
                            + edges[
                                bin_index + 1
                            ]
                        )
                    ),
                    "mean_number_density_nm-3": (
                        float(
                            density_frames[
                                :,
                                bin_index
                            ].mean()
                        )
                    ),
                    "std_number_density_nm-3": (
                        float(
                            density_frames[
                                :,
                                bin_index
                            ].std()
                        )
                    ),
                    "total_waterO_observations": (
                        int(
                            dataset[
                                "radial_count_total"
                            ][bin_index]
                        )
                    ),
                }
            )

    write_csv(
        RADIAL_CSV,
        radial_rows,
    )

    axial_rows = []

    for dataset_label, dataset in (
        ("frozen", frozen),
        ("mobile", mobile),
    ):
        density_frames = dataset[
            "axial_density_frames"
        ]

        edges = dataset[
            "axial_edges"
        ]

        for bin_index in range(
            AXIAL_BIN_COUNT
        ):
            axial_rows.append(
                {
                    "dataset": dataset_label,
                    "bin": bin_index,
                    "normalized_axial_lower": (
                        edges[bin_index]
                    ),
                    "normalized_axial_upper": (
                        edges[
                            bin_index + 1
                        ]
                    ),
                    "normalized_axial_center": (
                        0.5
                        * (
                            edges[bin_index]
                            + edges[
                                bin_index + 1
                            ]
                        )
                    ),
                    "mean_number_density_nm-3": (
                        float(
                            density_frames[
                                :,
                                bin_index
                            ].mean()
                        )
                    ),
                    "std_number_density_nm-3": (
                        float(
                            density_frames[
                                :,
                                bin_index
                            ].std()
                        )
                    ),
                    "total_waterO_observations": (
                        int(
                            dataset[
                                "axial_count_total"
                            ][bin_index]
                        )
                    ),
                }
            )

    write_csv(
        AXIAL_CSV,
        axial_rows,
    )

    mobile_occupancy = np.array(
        [
            float(
                row[
                    "lumen_water_count"
                ]
            )
            for row in mobile["rows"]
        ],
        dtype=float,
    )

    frozen_occupancy = np.array(
        [
            float(
                row[
                    "lumen_water_count"
                ]
            )
            for row in frozen["rows"]
        ],
        dtype=float,
    )

    mobile_density = np.array(
        [
            float(
                row[
                    "lumen_number_density_nm-3"
                ]
            )
            for row in mobile["rows"]
        ],
        dtype=float,
    )

    frozen_density = np.array(
        [
            float(
                row[
                    "lumen_number_density_nm-3"
                ]
            )
            for row in frozen["rows"]
        ],
        dtype=float,
    )

    mobile_radial_mean = np.array(
        [
            float(
                row[
                    "mean_normalized_radial_position"
                ]
            )
            for row in mobile["rows"]
        ],
        dtype=float,
    )

    frozen_radial_mean = np.array(
        [
            float(
                row[
                    "mean_normalized_radial_position"
                ]
            )
            for row in frozen["rows"]
        ],
        dtype=float,
    )

    radial_js = js_divergence(
        mobile[
            "radial_count_total"
        ],
        frozen[
            "radial_count_total"
        ],
    )

    axial_js = js_divergence(
        mobile[
            "axial_count_total"
        ],
        frozen[
            "axial_count_total"
        ],
    )

    pyrene_summary: list[
        dict[str, object]
    ] = []

    maximum_hydration_difference_pct = 0.0

    for molecule_index in range(
        1,
        PYR_COUNT + 1,
    ):
        for cutoff_nm in (
            HYDRATION_CUTOFFS_NM
        ):
            cutoff_label = int(
                round(
                    cutoff_nm * 100
                )
            )

            field = (
                f"PYR{molecule_index}_"
                f"waterO_within_{cutoff_label}pm"
            )

            mobile_values = np.array(
                [
                    float(row[field])
                    for row in mobile[
                        "rows"
                    ]
                ],
                dtype=float,
            )

            frozen_values = np.array(
                [
                    float(row[field])
                    for row in frozen[
                        "rows"
                    ]
                ],
                dtype=float,
            )

            difference_pct = (
                percent_difference(
                    float(
                        mobile_values.mean()
                    ),
                    float(
                        frozen_values.mean()
                    ),
                )
            )

            if math.isfinite(
                difference_pct
            ):
                maximum_hydration_difference_pct = max(
                    maximum_hydration_difference_pct,
                    abs(
                        difference_pct
                    ),
                )

            pyrene_summary.append(
                {
                    "pyrene": molecule_index,
                    "cutoff_nm": cutoff_nm,
                    "mobile_mean_waterO_count": (
                        float(
                            mobile_values.mean()
                        )
                    ),
                    "mobile_std_waterO_count": (
                        float(
                            mobile_values.std()
                        )
                    ),
                    "frozen_mean_waterO_count": (
                        float(
                            frozen_values.mean()
                        )
                    ),
                    "frozen_std_waterO_count": (
                        float(
                            frozen_values.std()
                        )
                    ),
                    "mobile_vs_frozen_difference_pct": (
                        difference_pct
                    ),
                }
            )

    representative = representative_indices(
        mobile["rows"],
        REPRESENTATIVE_PAIR_COUNT,
    )

    snapshot_rows = []

    for rank, frame_index in enumerate(
        representative,
        start=1,
    ):
        mobile_row = mobile[
            "rows"
        ][frame_index]

        snapshot_rows.append(
            {
                "selection_rank": rank,
                "frame": frame_index,
                "time_ps": (
                    mobile_row["relative_time_ps"]
                ),
                "mobile_frame": relative(
                    mobile_frames[
                        frame_index
                    ]
                ),
                "matched_frozen_frame": relative(
                    frozen_frames[
                        frame_index
                    ]
                ),
                "mobile_lumen_water_count": (
                    mobile_row[
                        "lumen_water_count"
                    ]
                ),
                "mobile_lumen_number_density_nm-3": (
                    mobile_row[
                        "lumen_number_density_nm-3"
                    ]
                ),
                "mobile_mean_normalized_radial_position": (
                    mobile_row[
                        "mean_normalized_radial_position"
                    ]
                ),
                "selection_status": (
                    "CANDIDATE_NOT_YET_AUTHORIZED_FOR_QM"
                ),
            }
        )

    write_csv(
        SNAPSHOT_CSV,
        snapshot_rows,
    )

    occupancy_difference_pct = (
        percent_difference(
            float(
                mobile_occupancy.mean()
            ),
            float(
                frozen_occupancy.mean()
            ),
        )
    )

    density_difference_pct = (
        percent_difference(
            float(
                mobile_density.mean()
            ),
            float(
                frozen_density.mean()
            ),
        )
    )

    mobility_effect_flags = []

    if abs(
        occupancy_difference_pct
    ) > 5.0:
        mobility_effect_flags.append(
            "mean lumen occupancy differs by more than 5%"
        )

    if radial_js > 0.02:
        mobility_effect_flags.append(
            "radial Jensen-Shannon divergence exceeds 0.02"
        )

    if axial_js > 0.02:
        mobility_effect_flags.append(
            "axial Jensen-Shannon divergence exceeds 0.02"
        )

    if (
        maximum_hydration_difference_pct
        > 10.0
    ):
        mobility_effect_flags.append(
            "local PYR hydration differs by more than 10%"
        )

    if mobility_effect_flags:
        interpretation = (
            "MOBILITY_DEPENDENT_REORGANIZATION_CANDIDATE"
        )
    else:
        interpretation = (
            "SIMILAR_WITHIN_CURRENT_SAMPLING"
        )

    summary: dict[str, object] = {
        "mobile_frames": len(
            mobile_frames
        ),
        "frozen_frames": len(
            frozen_frames
        ),
        "duration_ps": (
            EXPECTED_DURATION_PS
        ),
        "water_molecules": (
            WATER_COUNT
        ),
        "mobile_lumen_occupancy_mean": (
            float(
                mobile_occupancy.mean()
            )
        ),
        "mobile_lumen_occupancy_std": (
            float(
                mobile_occupancy.std()
            )
        ),
        "frozen_lumen_occupancy_mean": (
            float(
                frozen_occupancy.mean()
            )
        ),
        "frozen_lumen_occupancy_std": (
            float(
                frozen_occupancy.std()
            )
        ),
        "lumen_occupancy_difference_pct": (
            occupancy_difference_pct
        ),
        "mobile_lumen_density_mean_nm-3": (
            float(
                mobile_density.mean()
            )
        ),
        "mobile_lumen_density_std_nm-3": (
            float(
                mobile_density.std()
            )
        ),
        "frozen_lumen_density_mean_nm-3": (
            float(
                frozen_density.mean()
            )
        ),
        "frozen_lumen_density_std_nm-3": (
            float(
                frozen_density.std()
            )
        ),
        "lumen_density_difference_pct": (
            density_difference_pct
        ),
        "mobile_mean_normalized_radial_position": (
            float(
                np.nanmean(
                    mobile_radial_mean
                )
            )
        ),
        "frozen_mean_normalized_radial_position": (
            float(
                np.nanmean(
                    frozen_radial_mean
                )
            )
        ),
        "radial_JS_divergence": (
            radial_js
        ),
        "axial_JS_divergence": (
            axial_js
        ),
        "mobile_occupancy_slope_water_per_ps": (
            linear_slope(
                mobile["times"],
                mobile_occupancy,
            )
        ),
        "frozen_occupancy_slope_water_per_ps": (
            linear_slope(
                frozen["times"],
                frozen_occupancy,
            )
        ),
        "maximum_absolute_PYR_hydration_difference_pct": (
            maximum_hydration_difference_pct
        ),
        "representative_snapshot_pairs": (
            len(snapshot_rows)
        ),
        "analysis_status": "COMPLETE",
        "screening_interpretation": (
            interpretation
        ),
        "screening_flags": (
            " | ".join(
                mobility_effect_flags
            )
        ),
        "authorized_next_step": (
            "REVIEW_WATER_COMPARISON_AND_"
            "PREPARE_MOBILE_ELECTRONIC_SNAPSHOT_SET"
        ),
        "electronic_recalculation_authorized": (
            False
        ),
        "longer_mobile_production_authorized": (
            False
        ),
    }

    for row in pyrene_summary:
        molecule_index = int(
            row["pyrene"]
        )

        cutoff_label = int(
            round(
                float(
                    row["cutoff_nm"]
                )
                * 100
            )
        )

        prefix = (
            f"PYR{molecule_index}_"
            f"{cutoff_label}pm"
        )

        summary[
            f"{prefix}_mobile_mean"
        ] = row[
            "mobile_mean_waterO_count"
        ]

        summary[
            f"{prefix}_frozen_mean"
        ] = row[
            "frozen_mean_waterO_count"
        ]

        summary[
            f"{prefix}_difference_pct"
        ] = row[
            "mobile_vs_frozen_difference_pct"
        ]

    write_csv(
        SUMMARY_CSV,
        [summary],
    )

    REPORT_MD.write_text(
        f"""# Mobile versus Frozen Water Comparison

## Scope

- Mobile trajectory: **201 frames / 100 ps**
- Frozen-solute trajectory: **201 frames / 100 ps**
- Water molecules: **{WATER_COUNT}**
- Geometric lumen definition: water oxygen inside the instantaneous
  HBN wall radius and between the 1st and 99th percentiles of the
  HBN axial coordinates.

## Lumen water

- Mobile occupancy mean/std:
  {summary['mobile_lumen_occupancy_mean']:.4f}/
  {summary['mobile_lumen_occupancy_std']:.4f}
- Frozen occupancy mean/std:
  {summary['frozen_lumen_occupancy_mean']:.4f}/
  {summary['frozen_lumen_occupancy_std']:.4f}
- Mobile versus frozen occupancy difference:
  {occupancy_difference_pct:.4f} %

- Mobile density mean/std:
  {summary['mobile_lumen_density_mean_nm-3']:.6f}/
  {summary['mobile_lumen_density_std_nm-3']:.6f} nm^-3
- Frozen density mean/std:
  {summary['frozen_lumen_density_mean_nm-3']:.6f}/
  {summary['frozen_lumen_density_std_nm-3']:.6f} nm^-3
- Mobile versus frozen density difference:
  {density_difference_pct:.4f} %

## Spatial distributions

- Radial Jensen-Shannon divergence:
  {radial_js:.8f}
- Axial Jensen-Shannon divergence:
  {axial_js:.8f}
- Mobile/frozen mean normalized radial position:
  {summary['mobile_mean_normalized_radial_position']:.6f}/
  {summary['frozen_mean_normalized_radial_position']:.6f}

## Temporal behavior

- Mobile occupancy slope:
  {summary['mobile_occupancy_slope_water_per_ps']:.8f}
  water molecules ps^-1
- Frozen occupancy slope:
  {summary['frozen_occupancy_slope_water_per_ps']:.8f}
  water molecules ps^-1

## PYR hydration

- Maximum absolute difference among all PYR/cutoff comparisons:
  {maximum_hydration_difference_pct:.4f} %

## Screening interpretation

- Analysis status: **COMPLETE**
- Interpretation: **{interpretation}**
- Screening flags:
  {summary['screening_flags'] or 'NONE'}

These screening thresholds are comparative diagnostics and are not
force-field acceptance criteria.

## Snapshot candidates

- Representative matched mobile/frozen snapshot pairs:
  **{len(snapshot_rows)}**
- Electronic recalculation authorized by this analysis:
  **NO**

The candidate set can be reviewed before recalculating electronic
properties and comparing against the existing time-dependent
solvent-induced site energies under frozen-solute conditions.
""",
        encoding="utf-8",
    )

    print(
        "Day022 mobile-vs-frozen water comparison completed."
    )

    print(
        "Mobile / frozen frames: "
        f"{len(mobile_frames)} / "
        f"{len(frozen_frames)}"
    )

    print(
        "Mobile lumen occupancy mean/std: "
        f"{summary['mobile_lumen_occupancy_mean']:.4f}/"
        f"{summary['mobile_lumen_occupancy_std']:.4f}"
    )

    print(
        "Frozen lumen occupancy mean/std: "
        f"{summary['frozen_lumen_occupancy_mean']:.4f}/"
        f"{summary['frozen_lumen_occupancy_std']:.4f}"
    )

    print(
        "Occupancy difference mobile vs frozen: "
        f"{occupancy_difference_pct:.4f} %"
    )

    print(
        "Mobile/frozen lumen density mean: "
        f"{summary['mobile_lumen_density_mean_nm-3']:.6f}/"
        f"{summary['frozen_lumen_density_mean_nm-3']:.6f} nm^-3"
    )

    print(
        "Radial / axial JS divergence: "
        f"{radial_js:.8f}/"
        f"{axial_js:.8f}"
    )

    print(
        "Mobile/frozen occupancy slopes: "
        f"{summary['mobile_occupancy_slope_water_per_ps']:.8f}/"
        f"{summary['frozen_occupancy_slope_water_per_ps']:.8f} "
        "water/ps"
    )

    print(
        "Maximum absolute PYR hydration difference: "
        f"{maximum_hydration_difference_pct:.4f} %"
    )

    print(
        f"Screening interpretation: {interpretation}"
    )

    if mobility_effect_flags:
        print(
            "Screening flags: "
            + " | ".join(
                mobility_effect_flags
            )
        )

    print(
        "Representative matched snapshot pairs: "
        f"{len(snapshot_rows)}"
    )

    print(
        "Electronic recalculation authorized: NO"
    )

    print(
        "Longer mobile production authorized: NO"
    )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )


if __name__ == "__main__":
    main()
