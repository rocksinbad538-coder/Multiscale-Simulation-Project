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
    / "rmsf_tube_frozen_comparison"
)

FROZEN_FRAME_ROOT = (
    ANALYSIS_ROOT
    / "frozen_solute_frames"
)

INDEX_FILE = (
    ANALYSIS_ROOT
    / "solute_index.ndx"
)

FROZEN_EXTRACTION_LOG = (
    ANALYSIS_ROOT
    / "frozen_solute_frame_extraction.log"
)

TUBE_CSV = (
    ANALYSIS_ROOT
    / "stage08_mobile_frozen_tube_geometry_timeseries.csv"
)

HBN_RMSF_CSV = (
    ANALYSIS_ROOT
    / "stage08_mobile_frozen_hbn_rmsf.csv"
)

PYR_SUMMARY_CSV = (
    ANALYSIS_ROOT
    / "stage08_mobile_frozen_pyrene_rmsf_summary.csv"
)

SUMMARY_CSV = (
    ANALYSIS_ROOT
    / "stage08_rmsf_tube_frozen_comparison_summary.csv"
)

REPORT_MD = (
    ANALYSIS_ROOT
    / "STAGE08_RMSF_TUBE_FROZEN_COMPARISON_DAY022.md"
)

HBN_COUNT = 1680
PYR_COUNT = 4
PYR_ATOMS = 26
SOLUTE_COUNT = HBN_COUNT + PYR_COUNT * PYR_ATOMS

EXPECTED_FRAMES = 201
EXPECTED_DURATION_PS = 100.0


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
    )


def require_file(path: Path) -> None:
    if (
        not path.exists()
        or path.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing or empty file: {path}"
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


def locate_gmx() -> str:
    executable = shutil.which("gmx")

    if executable:
        return executable

    default = Path(
        "/usr/local/gromacs/bin/gmx"
    )

    if default.exists():
        return str(default)

    raise RuntimeError(
        "Could not locate GROMACS"
    )


def locate_trajectory_file(
    directory: Path,
    preferred_name: str,
    suffix: str,
) -> Path:
    preferred = directory / preferred_name

    if preferred.exists():
        return preferred

    matches = sorted(
        path
        for path in directory.glob(
            f"*{suffix}"
        )
        if path.is_file()
        and path.stat().st_size > 0
    )

    if len(matches) != 1:
        raise RuntimeError(
            f"Could not identify a unique {suffix} file "
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
            "Expected exactly 201 Stage08 frames; "
            f"found {len(frames)}"
        )

    return frames


def write_index() -> None:
    with INDEX_FILE.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "[ SOLUTE ]\n"
        )

        indices = list(
            range(
                1,
                SOLUTE_COUNT + 1,
            )
        )

        for start in range(
            0,
            len(indices),
            15,
        ):
            handle.write(
                " ".join(
                    str(value)
                    for value in indices[
                        start : start + 15
                    ]
                )
                + "\n"
            )


def extract_frozen_frames(
    gmx: str,
    frozen_xtc: Path,
    frozen_tpr: Path,
) -> list[Path]:
    FROZEN_FRAME_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for old_frame in FROZEN_FRAME_ROOT.glob(
        "frozen_solute*.gro"
    ):
        old_frame.unlink()

    output_template = (
        FROZEN_FRAME_ROOT
        / "frozen_solute.gro"
    )

    command = [
        gmx,
        "trjconv",
        "-f",
        str(frozen_xtc),
        "-s",
        str(frozen_tpr),
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
        input="SOLUTE\n",
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
            "Frozen trajectory extraction failed. "
            f"See {FROZEN_EXTRACTION_LOG}"
        )

    frames = sorted(
        FROZEN_FRAME_ROOT.glob(
            "frozen_solute*.gro"
        ),
        key=natural_key,
    )

    if len(frames) != EXPECTED_FRAMES:
        raise RuntimeError(
            "Expected exactly 201 frozen frames; "
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
    natoms = int(
        lines[1].strip()
    )

    if natoms < SOLUTE_COUNT:
        raise RuntimeError(
            f"Too few atoms in {path}: {natoms}"
        )

    atoms = []

    for local_index, line in enumerate(
        lines[
            2 : 2 + SOLUTE_COUNT
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

    if len(box_values) < 3:
        raise RuntimeError(
            f"Invalid simulation box in {path}"
        )

    return (
        title,
        atoms,
        np.array(
            box_values[:3],
            dtype=float,
        ),
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
    vector: np.ndarray,
    box: np.ndarray,
) -> np.ndarray:
    return (
        vector
        - box
        * np.round(
            vector / box
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


def rmsf(
    coordinate_stack: np.ndarray,
) -> np.ndarray:
    mean_coordinates = np.mean(
        coordinate_stack,
        axis=0,
    )

    fluctuations = (
        coordinate_stack
        - mean_coordinates
    )

    return np.sqrt(
        np.mean(
            np.sum(
                fluctuations ** 2,
                axis=2,
            ),
            axis=0,
        )
    )


def rms_to_reference(
    coordinates: np.ndarray,
    reference: np.ndarray,
) -> float:
    return float(
        np.sqrt(
            np.mean(
                np.sum(
                    (
                        coordinates
                        - reference
                    )
                    ** 2,
                    axis=1,
                )
            )
        )
    )


def linear_slope(
    times: np.ndarray,
    values: np.ndarray,
) -> float:
    if len(times) < 2:
        return 0.0

    return float(
        np.polyfit(
            times,
            values,
            1,
        )[0]
    )


def q(
    values: np.ndarray,
    probability: float,
) -> float:
    return float(
        np.quantile(
            values,
            probability,
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

    axis = (
        axis
        / np.linalg.norm(axis)
    )

    return axis


def transverse_basis(
    axis: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
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


def tube_metrics(
    coordinates: np.ndarray,
    axis: np.ndarray,
    basis_1: np.ndarray,
    basis_2: np.ndarray,
) -> dict[str, float]:
    center = coordinates.mean(
        axis=0
    )

    centered = (
        coordinates
        - center
    )

    axial = (
        centered
        @ axis
    )

    transverse_1 = (
        centered
        @ basis_1
    )

    transverse_2 = (
        centered
        @ basis_2
    )

    lower_1 = q(
        axial,
        0.01,
    )

    upper_1 = q(
        axial,
        0.99,
    )

    lower_10 = q(
        axial,
        0.10,
    )

    upper_10 = q(
        axial,
        0.90,
    )

    central_mask = (
        (axial >= lower_10)
        & (axial <= upper_10)
    )

    transverse = np.column_stack(
        (
            transverse_1[
                central_mask
            ],
            transverse_2[
                central_mask
            ],
        )
    )

    radial = np.sqrt(
        np.sum(
            transverse ** 2,
            axis=1,
        )
    )

    covariance = np.cov(
        transverse,
        rowvar=False,
    )

    eigenvalues = np.linalg.eigvalsh(
        covariance
    )

    eigenvalues = np.sort(
        np.maximum(
            eigenvalues,
            0.0,
        )
    )[::-1]

    semi_axis_major = math.sqrt(
        2.0 * eigenvalues[0]
    )

    semi_axis_minor = math.sqrt(
        2.0 * eigenvalues[1]
    )

    ellipticity = (
        semi_axis_major
        / semi_axis_minor
        if semi_axis_minor > 0.0
        else math.inf
    )

    return {
        "length_p98_nm": (
            upper_1 - lower_1
        ),
        "length_full_nm": float(
            axial.max()
            - axial.min()
        ),
        "radius_mean_nm": float(
            radial.mean()
        ),
        "radius_std_nm": float(
            radial.std()
        ),
        "radius_min_nm": float(
            radial.min()
        ),
        "radius_max_nm": float(
            radial.max()
        ),
        "semi_axis_major_nm": (
            semi_axis_major
        ),
        "semi_axis_minor_nm": (
            semi_axis_minor
        ),
        "ellipticity": ellipticity,
    }


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


def process_trajectory(
    label: str,
    frame_paths: list[Path],
    integrated,
    hbn_topology,
    pyr_topology,
    reference_hbn: np.ndarray,
    reference_pyrenes: list[np.ndarray],
    axis: np.ndarray,
    basis_1: np.ndarray,
    basis_2: np.ndarray,
) -> dict[str, object]:
    hbn_stack = []
    solute_stack = []

    pyr_internal_stacks = [
        []
        for _ in range(
            PYR_COUNT
        )
    ]

    pyr_centroids = [
        []
        for _ in range(
            PYR_COUNT
        )
    ]

    geometry_rows = []
    first_atom_metadata = None
    times = []

    for frame_index, frame_path in enumerate(
        frame_paths
    ):
        title, atoms, box = read_gro(
            frame_path
        )

        if first_atom_metadata is None:
            first_atom_metadata = atoms

        time_ps = parse_time(
            title,
            fallback=0.5 * frame_index,
        )

        times.append(time_ps)

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

        aligned_solute = np.zeros(
            (
                SOLUTE_COUNT,
                3,
            ),
            dtype=float,
        )

        aligned_solute[
            :HBN_COUNT
        ] = aligned_hbn

        hbn_center = current_hbn.mean(
            axis=0
        )

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

            current_pyr_wrapped = (
                positions[start:stop]
            )

            current_pyr = (
                integrated.unwrap_by_bonds(
                    current_pyr_wrapped,
                    box,
                    pyr_topology["bonds"],
                )
            )

            pyr_center = current_pyr.mean(
                axis=0
            )

            nearest_delta = minimum_image(
                pyr_center
                - hbn_center,
                box,
            )

            target_pyr_center = (
                hbn_center
                + nearest_delta
            )

            current_pyr = (
                current_pyr
                + target_pyr_center
                - pyr_center
            )

            aligned_relative = apply_transform(
                current_pyr,
                rotation,
                target_center,
                reference_center,
            )

            aligned_solute[
                start:stop
            ] = aligned_relative

            pyr_centroids[
                molecule_index
            ].append(
                aligned_relative.mean(
                    axis=0
                )
            )

            (
                internally_aligned,
                _,
                _,
                _,
            ) = kabsch_transform(
                reference_pyrenes[
                    molecule_index
                ],
                current_pyr,
            )

            pyr_internal_stacks[
                molecule_index
            ].append(
                internally_aligned
            )

        metrics = tube_metrics(
            aligned_hbn,
            axis,
            basis_1,
            basis_2,
        )

        metrics.update(
            {
                "dataset": label,
                "frame": frame_index,
                "time_ps": time_ps,
                "HBN_rms_to_frozen_reference_nm": (
                    rms_to_reference(
                        aligned_hbn,
                        reference_hbn,
                    )
                ),
            }
        )

        geometry_rows.append(
            metrics
        )

        hbn_stack.append(
            aligned_hbn
        )

        solute_stack.append(
            aligned_solute
        )

    times_array = np.array(
        times,
        dtype=float,
    )

    times_array = (
        times_array
        - times_array[0]
    )

    if not math.isclose(
        float(
            times_array[-1]
        ),
        EXPECTED_DURATION_PS,
        rel_tol=0.0,
        abs_tol=1.0e-6,
    ):
        raise RuntimeError(
            f"{label} duration is not 100 ps: "
            f"{times_array[-1]}"
        )

    hbn_array = np.array(
        hbn_stack,
        dtype=float,
    )

    solute_array = np.array(
        solute_stack,
        dtype=float,
    )

    hbn_rmsf = rmsf(
        hbn_array
    )

    solute_rmsf = rmsf(
        solute_array
    )

    pyr_internal_results = []
    pyr_relative_results = []

    for molecule_index in range(
        PYR_COUNT
    ):
        start = (
            HBN_COUNT
            + molecule_index
            * PYR_ATOMS
        )

        stop = start + PYR_ATOMS

        internal_array = np.array(
            pyr_internal_stacks[
                molecule_index
            ],
            dtype=float,
        )

        internal_rmsf = rmsf(
            internal_array
        )

        relative_atom_rmsf = (
            solute_rmsf[
                start:stop
            ]
        )

        centroid_array = np.array(
            pyr_centroids[
                molecule_index
            ],
            dtype=float,
        )

        centroid_mean = centroid_array.mean(
            axis=0
        )

        centroid_rmsf = float(
            np.sqrt(
                np.mean(
                    np.sum(
                        (
                            centroid_array
                            - centroid_mean
                        )
                        ** 2,
                        axis=1,
                    )
                )
            )
        )

        pyr_internal_results.append(
            internal_rmsf
        )

        pyr_relative_results.append(
            {
                "dataset": label,
                "pyrene": molecule_index + 1,
                "centroid_rmsf_relative_to_HBN_nm": (
                    centroid_rmsf
                ),
                "relative_atom_rmsf_mean_nm": float(
                    relative_atom_rmsf.mean()
                ),
                "relative_atom_rmsf_q95_nm": q(
                    relative_atom_rmsf,
                    0.95,
                ),
                "relative_atom_rmsf_max_nm": float(
                    relative_atom_rmsf.max()
                ),
                "internal_atom_rmsf_mean_nm": float(
                    internal_rmsf.mean()
                ),
                "internal_atom_rmsf_q95_nm": q(
                    internal_rmsf,
                    0.95,
                ),
                "internal_atom_rmsf_max_nm": float(
                    internal_rmsf.max()
                ),
            }
        )

    geometry_by_field = {}

    for field in (
        "length_p98_nm",
        "length_full_nm",
        "radius_mean_nm",
        "radius_std_nm",
        "semi_axis_major_nm",
        "semi_axis_minor_nm",
        "ellipticity",
        "HBN_rms_to_frozen_reference_nm",
    ):
        geometry_by_field[field] = np.array(
            [
                float(row[field])
                for row in geometry_rows
            ],
            dtype=float,
        )

    return {
        "times": times_array,
        "geometry_rows": geometry_rows,
        "geometry": geometry_by_field,
        "hbn_coordinates": hbn_array,
        "hbn_rmsf": hbn_rmsf,
        "solute_rmsf": solute_rmsf,
        "pyr_internal_rmsf": (
            pyr_internal_results
        ),
        "pyr_summary_rows": (
            pyr_relative_results
        ),
        "atom_metadata": first_atom_metadata,
    }


def mean_structure_rmsd(
    coordinates: np.ndarray,
    reference: np.ndarray,
) -> float:
    mean_coordinates = np.mean(
        coordinates,
        axis=0,
    )

    return rms_to_reference(
        mean_coordinates,
        reference,
    )


def percent_difference(
    mobile: float,
    frozen: float,
) -> float:
    if frozen == 0.0:
        return math.nan

    return (
        100.0
        * (
            mobile
            - frozen
        )
        / frozen
    )


def main() -> None:
    require_file(
        INTEGRATED_MODULE_PATH
    )

    require_file(
        HBN_ITP
    )

    require_file(
        PYR_ITP
    )

    require_file(
        FROZEN_RUN
    )

    ANALYSIS_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    integrated = load_module(
        INTEGRATED_MODULE_PATH,
        "stage08_integrated_module",
    )

    gmx = locate_gmx()

    mobile_frames = list_mobile_frames()

    frozen_xtc = locate_trajectory_file(
        FROZEN_RUN,
        "nvt_100ps_frozenSolute.xtc",
        ".xtc",
    )

    frozen_tpr = locate_trajectory_file(
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

    reference_pyrenes = []

    for molecule_index in range(
        PYR_COUNT
    ):
        start = (
            HBN_COUNT
            + molecule_index
            * PYR_ATOMS
        )

        stop = start + PYR_ATOMS

        reference_pyrenes.append(
            integrated.unwrap_by_bonds(
                frozen_reference_positions[
                    start:stop
                ],
                frozen_reference_box,
                pyr_topology["bonds"],
            )
        )

    axis = tube_axis(
        reference_hbn
    )

    basis_1, basis_2 = transverse_basis(
        axis
    )

    frozen = process_trajectory(
        "frozen",
        frozen_frames,
        integrated,
        hbn_topology,
        pyr_topology,
        reference_hbn,
        reference_pyrenes,
        axis,
        basis_1,
        basis_2,
    )

    mobile = process_trajectory(
        "mobile",
        mobile_frames,
        integrated,
        hbn_topology,
        pyr_topology,
        reference_hbn,
        reference_pyrenes,
        axis,
        basis_1,
        basis_2,
    )

    tube_rows = (
        frozen["geometry_rows"]
        + mobile["geometry_rows"]
    )

    write_csv(
        TUBE_CSV,
        tube_rows,
    )

    hbn_rmsf_rows = []

    atom_metadata = mobile[
        "atom_metadata"
    ]

    for atom_index in range(
        HBN_COUNT
    ):
        atom = atom_metadata[
            atom_index
        ]

        hbn_rmsf_rows.append(
            {
                "atom_index": (
                    atom_index + 1
                ),
                "residue_name": (
                    atom["residue_name"]
                ),
                "atom_name": (
                    atom["atom_name"]
                ),
                "frozen_rmsf_nm": float(
                    frozen[
                        "hbn_rmsf"
                    ][atom_index]
                ),
                "mobile_rmsf_nm": float(
                    mobile[
                        "hbn_rmsf"
                    ][atom_index]
                ),
            }
        )

    write_csv(
        HBN_RMSF_CSV,
        hbn_rmsf_rows,
    )

    pyr_summary_rows = (
        frozen["pyr_summary_rows"]
        + mobile["pyr_summary_rows"]
    )

    write_csv(
        PYR_SUMMARY_CSV,
        pyr_summary_rows,
    )

    frozen_geometry = frozen[
        "geometry"
    ]

    mobile_geometry = mobile[
        "geometry"
    ]

    frozen_radius = float(
        frozen_geometry[
            "radius_mean_nm"
        ].mean()
    )

    mobile_radius = float(
        mobile_geometry[
            "radius_mean_nm"
        ].mean()
    )

    frozen_length = float(
        frozen_geometry[
            "length_p98_nm"
        ].mean()
    )

    mobile_length = float(
        mobile_geometry[
            "length_p98_nm"
        ].mean()
    )

    frozen_ellipticity = float(
        frozen_geometry[
            "ellipticity"
        ].mean()
    )

    mobile_ellipticity = float(
        mobile_geometry[
            "ellipticity"
        ].mean()
    )

    mobile_hbn_rmsf = mobile[
        "hbn_rmsf"
    ]

    frozen_hbn_rmsf = frozen[
        "hbn_rmsf"
    ]

    maximum_mobile_internal_pyr_rmsf = max(
        float(
            row[
                "internal_atom_rmsf_max_nm"
            ]
        )
        for row in mobile[
            "pyr_summary_rows"
        ]
    )

    blocked_reasons: list[str] = []
    review_reasons: list[str] = []

    numeric_values = np.concatenate(
        [
            mobile_hbn_rmsf,
            frozen_hbn_rmsf,
            mobile_geometry[
                "radius_mean_nm"
            ],
            mobile_geometry[
                "length_p98_nm"
            ],
            mobile_geometry[
                "ellipticity"
            ],
        ]
    )

    if not np.all(
        np.isfinite(
            numeric_values
        )
    ):
        blocked_reasons.append(
            "non-finite RMSF or tube-geometry metric"
        )

    radius_difference_pct = (
        percent_difference(
            mobile_radius,
            frozen_radius,
        )
    )

    length_difference_pct = (
        percent_difference(
            mobile_length,
            frozen_length,
        )
    )

    severe_checks = {
        "mobile HBN RMSF maximum exceeds 0.50 nm": (
            float(
                mobile_hbn_rmsf.max()
            )
            > 0.50
        ),
        "mean radius differs from frozen by more than 20%": (
            abs(
                radius_difference_pct
            )
            > 20.0
        ),
        "mean length differs from frozen by more than 20%": (
            abs(
                length_difference_pct
            )
            > 20.0
        ),
        "tube ellipticity exceeds 1.60": (
            float(
                mobile_geometry[
                    "ellipticity"
                ].max()
            )
            > 1.60
        ),
        "PYR internal RMSF exceeds 0.15 nm": (
            maximum_mobile_internal_pyr_rmsf
            > 0.15
        ),
    }

    blocked_reasons.extend(
        label
        for label, failed
        in severe_checks.items()
        if failed
    )

    if q(
        mobile_hbn_rmsf,
        0.95,
    ) > 0.15:
        review_reasons.append(
            "mobile HBN RMSF q95 exceeds 0.15 nm"
        )

    if float(
        mobile_hbn_rmsf.max()
    ) > 0.35:
        review_reasons.append(
            "mobile HBN RMSF maximum exceeds 0.35 nm"
        )

    if abs(
        radius_difference_pct
    ) > 10.0:
        review_reasons.append(
            "mean radius differs from frozen by more than 10%"
        )

    if abs(
        length_difference_pct
    ) > 10.0:
        review_reasons.append(
            "mean length differs from frozen by more than 10%"
        )

    if float(
        mobile_geometry[
            "ellipticity"
        ].mean()
    ) > 1.20:
        review_reasons.append(
            "mean tube ellipticity exceeds 1.20"
        )

    if float(
        mobile_geometry[
            "ellipticity"
        ].max()
    ) > 1.35:
        review_reasons.append(
            "maximum tube ellipticity exceeds 1.35"
        )

    if maximum_mobile_internal_pyr_rmsf > 0.08:
        review_reasons.append(
            "PYR internal RMSF exceeds 0.08 nm"
        )

    if float(
        frozen_hbn_rmsf.max()
    ) > 0.005:
        review_reasons.append(
            "frozen HBN control RMSF exceeds 0.005 nm"
        )

    radius_slope = linear_slope(
        mobile["times"],
        mobile_geometry[
            "radius_mean_nm"
        ],
    )

    length_slope = linear_slope(
        mobile["times"],
        mobile_geometry[
            "length_p98_nm"
        ],
    )

    ellipticity_slope = linear_slope(
        mobile["times"],
        mobile_geometry[
            "ellipticity"
        ],
    )

    if (
        abs(radius_slope) > 0.001
    ):
        review_reasons.append(
            "mobile tube radius slope exceeds 0.001 nm/ps"
        )

    if (
        abs(length_slope) > 0.005
    ):
        review_reasons.append(
            "mobile tube length slope exceeds 0.005 nm/ps"
        )

    if blocked_reasons:
        decision = "BLOCKED"
        next_step = (
            "REVIEW_HBN_FORCE_FIELD_OR_STRUCTURE"
        )
    elif review_reasons:
        decision = "REVIEW"
        next_step = (
            "TARGETED_RMSF_OR_TUBE_GEOMETRY_REVIEW"
        )
    else:
        decision = "PASS_WITH_MONITORING"
        next_step = (
            "MOBILE_VS_FROZEN_WATER_AND_"
            "SOLVENT_DISORDER_COMPARISON"
        )

    summary = {
        "mobile_frames": len(
            mobile_frames
        ),
        "frozen_frames": len(
            frozen_frames
        ),
        "duration_ps": (
            EXPECTED_DURATION_PS
        ),
        "mobile_HBN_rmsf_mean_nm": float(
            mobile_hbn_rmsf.mean()
        ),
        "mobile_HBN_rmsf_q95_nm": q(
            mobile_hbn_rmsf,
            0.95,
        ),
        "mobile_HBN_rmsf_max_nm": float(
            mobile_hbn_rmsf.max()
        ),
        "frozen_HBN_rmsf_max_nm": float(
            frozen_hbn_rmsf.max()
        ),
        "mobile_HBN_mean_structure_rmsd_to_frozen_nm": (
            mean_structure_rmsd(
                mobile[
                    "hbn_coordinates"
                ],
                reference_hbn,
            )
        ),
        "mobile_tube_radius_mean_nm": (
            mobile_radius
        ),
        "mobile_tube_radius_std_over_frames_nm": float(
            mobile_geometry[
                "radius_mean_nm"
            ].std()
        ),
        "frozen_tube_radius_mean_nm": (
            frozen_radius
        ),
        "tube_radius_mobile_vs_frozen_difference_pct": (
            radius_difference_pct
        ),
        "mobile_tube_length_p98_mean_nm": (
            mobile_length
        ),
        "mobile_tube_length_p98_std_nm": float(
            mobile_geometry[
                "length_p98_nm"
            ].std()
        ),
        "frozen_tube_length_p98_mean_nm": (
            frozen_length
        ),
        "tube_length_mobile_vs_frozen_difference_pct": (
            length_difference_pct
        ),
        "mobile_tube_ellipticity_mean": (
            mobile_ellipticity
        ),
        "mobile_tube_ellipticity_max": float(
            mobile_geometry[
                "ellipticity"
            ].max()
        ),
        "frozen_tube_ellipticity_mean": (
            frozen_ellipticity
        ),
        "mobile_radius_slope_nm_per_ps": (
            radius_slope
        ),
        "mobile_length_slope_nm_per_ps": (
            length_slope
        ),
        "mobile_ellipticity_slope_per_ps": (
            ellipticity_slope
        ),
        "mobile_PYR_max_internal_rmsf_nm": (
            maximum_mobile_internal_pyr_rmsf
        ),
        "comparison_decision": decision,
        "authorized_next_step": (
            next_step
        ),
        "longer_mobile_production_authorized": (
            False
        ),
        "review_reasons": (
            " | ".join(
                review_reasons
            )
        ),
        "blocked_reasons": (
            " | ".join(
                blocked_reasons
            )
        ),
    }

    write_csv(
        SUMMARY_CSV,
        [summary],
    )

    REPORT_MD.write_text(
        f"""# Stage08 RMSF, Tube Geometry, and Frozen-Control Comparison

- Mobile trajectory: **201 frames / 100 ps**
- Frozen-control trajectory: **201 frames / 100 ps**
- Decision: **{decision}**
- Authorized next step: `{next_step}`
- Longer mobile production authorized: **NO**

## HBN RMSF

- Mobile mean/q95/max RMSF: {summary['mobile_HBN_rmsf_mean_nm']:.6f}/{summary['mobile_HBN_rmsf_q95_nm']:.6f}/{summary['mobile_HBN_rmsf_max_nm']:.6f} nm
- Frozen-control maximum RMSF: {summary['frozen_HBN_rmsf_max_nm']:.6f} nm
- Mobile mean-structure RMSD to frozen reference: {summary['mobile_HBN_mean_structure_rmsd_to_frozen_nm']:.6f} nm

## Tube geometry

- Mobile/frozen mean radius: {mobile_radius:.6f}/{frozen_radius:.6f} nm
- Radius difference: {radius_difference_pct:.4f} %
- Mobile/frozen p98 length: {mobile_length:.6f}/{frozen_length:.6f} nm
- Length difference: {length_difference_pct:.4f} %
- Mobile mean/max ellipticity: {mobile_ellipticity:.6f}/{summary['mobile_tube_ellipticity_max']:.6f}
- Frozen mean ellipticity: {frozen_ellipticity:.6f}
- Radius slope: {radius_slope:.8f} nm ps^-1
- Length slope: {length_slope:.8f} nm ps^-1

## Pyrenes

- Maximum mobile internal PYR RMSF: {maximum_mobile_internal_pyr_rmsf:.6f} nm

## Scope note

The frozen-solute trajectory is used as a geometric and numerical
control. It cannot provide evidence of solute stability because HBN
and PYR coordinates were frozen by construction.
""",
        encoding="utf-8",
    )

    print(
        "Day022 Stage08 RMSF/tube/frozen comparison completed."
    )

    print(
        "Mobile / frozen frames: "
        f"{len(mobile_frames)} / "
        f"{len(frozen_frames)}"
    )

    print(
        "Mobile HBN RMSF mean/q95/max: "
        f"{summary['mobile_HBN_rmsf_mean_nm']:.6f}/"
        f"{summary['mobile_HBN_rmsf_q95_nm']:.6f}/"
        f"{summary['mobile_HBN_rmsf_max_nm']:.6f} nm"
    )

    print(
        "Frozen HBN RMSF maximum: "
        f"{summary['frozen_HBN_rmsf_max_nm']:.6f} nm"
    )

    print(
        "Mobile/frozen tube radius: "
        f"{mobile_radius:.6f}/"
        f"{frozen_radius:.6f} nm"
    )

    print(
        "Tube-radius difference: "
        f"{radius_difference_pct:.4f} %"
    )

    print(
        "Mobile/frozen p98 tube length: "
        f"{mobile_length:.6f}/"
        f"{frozen_length:.6f} nm"
    )

    print(
        "Tube-length difference: "
        f"{length_difference_pct:.4f} %"
    )

    print(
        "Mobile tube ellipticity mean/max: "
        f"{mobile_ellipticity:.6f}/"
        f"{summary['mobile_tube_ellipticity_max']:.6f}"
    )

    print(
        "Mobile radius/length slopes: "
        f"{radius_slope:.8f}/"
        f"{length_slope:.8f} nm/ps"
    )

    print(
        "Maximum mobile PYR internal RMSF: "
        f"{maximum_mobile_internal_pyr_rmsf:.6f} nm"
    )

    print(
        f"Comparison decision: {decision}"
    )

    print(
        f"Authorized next step: {next_step}"
    )

    print(
        "Longer mobile production authorized: NO"
    )

    if review_reasons:
        print(
            "Review reasons: "
            + " | ".join(
                review_reasons
            )
        )

    if blocked_reasons:
        print(
            "Blocking reasons: "
            + " | ".join(
                blocked_reasons
            )
        )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )


if __name__ == "__main__":
    main()
