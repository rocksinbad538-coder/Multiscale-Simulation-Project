#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

DAY023_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design"
)

PREPARATION_ROOT = (
    DAY023_ROOT
    / "07_r1_frozen_solute_nvt_20ps_preparation"
)

PROTOTYPE_ROOT = (
    DAY023_ROOT
    / "02_r1_steric_cap_prototype"
)

OUTPUT_ROOT = (
    DAY023_ROOT
    / "08_r1_frozen_solute_nvt_20ps"
)

SOURCE_TPR = (
    PREPARATION_ROOT
    / "r1_frozen_solute_nvt_20ps.tpr"
)

SOURCE_GRO = (
    PREPARATION_ROOT
    / "r1_frozen_solute_nvt_20ps_input.gro"
)

SOURCE_INDEX = (
    PREPARATION_ROOT
    / "r1_frozen_solute_nvt_20ps.ndx"
)

PREPARATION_SUMMARY = (
    PREPARATION_ROOT
    / "r1_frozen_solute_nvt_20ps_preparation_summary.csv"
)

WARNING_SUMMARY = (
    PREPARATION_ROOT
    / "r1_nvt_grompp_warning_authorization_summary.csv"
)

PROTOTYPE_JSON = (
    PROTOTYPE_ROOT
    / "r1_selected_steric_cap_definition.json"
)

LOCAL_TPR = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps.tpr"
)

LOCAL_INDEX = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps.ndx"
)

DEFFNM = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps"
)

FINAL_GRO = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps.gro"
)

XTC = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps.xtc"
)

EDR = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps.edr"
)

MDLOG = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps.log"
)

CPT = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps.cpt"
)

MDRUN_CONSOLE = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_mdrun_console.log"
)

TRAJECTORY_CHECK_LOG = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_trajectory_check.log"
)

TRJCONV_LOG = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_trjconv.log"
)

ENERGY_MENU_LOG = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_energy_menu.txt"
)

TEMPERATURE_XVG = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_temperature.xvg"
)

POTENTIAL_XVG = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_potential.xvg"
)

TOTAL_ENERGY_XVG = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_total_energy.xvg"
)

CAP_SOL_LJ_XVG = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_cap_sol_lj.xvg"
)

OCCUPANCY_CSV = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_lumen_occupancy.csv"
)

ENERGY_CSV = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_energy_series.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_summary.csv"
)

GATE_CSV = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_gates.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R1_FROZEN_SOLUTE_NVT_20PS_VALIDATION_DAY023.md"
)

EXPECTED_ATOMS = 68314

HBN_ATOMS = 1680
PYR_ATOMS = 104
SOLUTE_ATOMS = HBN_ATOMS + PYR_ATOMS

EXPECTED_WATERS = 16551
WATER_SITES = 4
WATER_ATOMS = EXPECTED_WATERS * WATER_SITES

CAP_ATOMS = 326
CAP_START = SOLUTE_ATOMS + WATER_ATOMS
CAP_STOP = CAP_START + CAP_ATOMS

INITIAL_LUMEN_WATERS = 428

EXPECTED_DURATION_PS = 20.0
EXPECTED_FRAME_INTERVAL_PS = 0.5
EXPECTED_FRAMES = 41

FROZEN_FINAL_TOLERANCE_NM = 1.0e-6
FROZEN_TRAJECTORY_TOLERANCE_NM = 1.5e-3

MIN_ACCEPTABLE_TEMPERATURE_K = 295.0
MAX_ACCEPTABLE_TEMPERATURE_K = 305.0
MAX_TEMPERATURE_STANDARD_DEVIATION_K = 10.0

MIN_OCCUPANCY_FRACTION = 0.90
MIN_INITIAL_WATER_RETENTION_FRACTION = 0.90

MIN_CAP_WATER_DISTANCE_NM = 0.15

INSTABILITY_PATTERNS = (
    r"\bnan\b",
    r"fatal error",
    r"segmentation fault",
    r"lincs warning",
    r"shake did not converge",
    r"constraint warning",
    r"water molecule cannot be settled",
    r"water molecule can not be settled",
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
            f"Missing or empty file: {path}"
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
        "Could not locate GROMACS."
    )


def run_command(
    command: list[str],
    *,
    cwd: Path,
    input_text: str | None = None,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=cwd,
        env={
            **os.environ,
            "GMX_MAXBACKUP": "-1",
        },
        input=input_text,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def run_live(
    command: list[str],
    *,
    cwd: Path,
    log_path: Path,
) -> int:
    process = subprocess.Popen(
        command,
        cwd=cwd,
        env={
            **os.environ,
            "GMX_MAXBACKUP": "-1",
        },
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,
    )

    if process.stdout is None:
        raise RuntimeError(
            "Could not capture mdrun output."
        )

    with log_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        for line in process.stdout:
            print(
                line,
                end="",
                flush=True,
            )

            handle.write(line)
            handle.flush()

    return process.wait()


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


def read_single_csv_row(
    path: Path,
) -> dict[str, str]:
    require_file(path)

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(
            csv.DictReader(handle)
        )

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected one row in {path}; "
            f"found {len(rows)}"
        )

    return rows[0]


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
                    key: row.get(
                        key,
                        "",
                    )
                    for key in fields
                }
            )


def parse_box(
    values: list[float],
) -> np.ndarray:
    if len(values) != 3:
        raise RuntimeError(
            "Validation requires an orthorhombic box."
        )

    box = np.asarray(
        values,
        dtype=float,
    )

    if (
        not np.all(
            np.isfinite(box)
        )
        or np.any(box <= 0.0)
    ):
        raise RuntimeError(
            "Invalid GRO box."
        )

    return box


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
    atom_count = int(
        lines[1].strip()
    )

    if len(lines) < atom_count + 3:
        raise RuntimeError(
            f"Incomplete GRO file: {path}"
        )

    atoms = []

    for index, line in enumerate(
        lines[
            2:
            2 + atom_count
        ]
    ):
        if len(line) < 44:
            raise RuntimeError(
                f"Malformed GRO atom line "
                f"{index + 1} in {path}"
            )

        atoms.append(
            {
                "index": index,
                "resid": int(
                    line[0:5]
                ),
                "resname": line[
                    5:10
                ].strip(),
                "atomname": line[
                    10:15
                ].strip(),
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

    box = parse_box(
        [
            float(value)
            for value in lines[
                2 + atom_count
            ].split()
        ]
    )

    return (
        title,
        atoms,
        box,
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


def displacement_metrics(
    reference: np.ndarray,
    current: np.ndarray,
    box: np.ndarray,
) -> tuple[
    float,
    float,
]:
    displacement = minimum_image(
        current - reference,
        box,
    )

    norms = np.linalg.norm(
        displacement,
        axis=1,
    )

    return (
        float(
            np.sqrt(
                np.mean(
                    norms * norms
                )
            )
        ),
        float(
            np.max(norms)
        ),
    )


def water_oxygen_indices(
    atoms: list[dict[str, Any]],
) -> np.ndarray:
    indices = []

    for water_index in range(
        EXPECTED_WATERS
    ):
        start = (
            SOLUTE_ATOMS
            + water_index
            * WATER_SITES
        )

        chunk = atoms[
            start:
            start + WATER_SITES
        ]

        matches = [
            start + local_index
            for local_index, atom
            in enumerate(chunk)
            if atom[
                "atomname"
            ]
            == "OW"
        ]

        if len(matches) != 1:
            raise RuntimeError(
                "Could not identify exactly one OW "
                f"in water {water_index}: "
                f"{[atom['atomname'] for atom in chunk]}"
            )

        indices.append(
            matches[0]
        )

    return np.asarray(
        indices,
        dtype=int,
    )


def lumen_mask(
    water_positions: np.ndarray,
    box: np.ndarray,
    prototype: dict[str, Any],
) -> np.ndarray:
    center = np.asarray(
        prototype[
            "tube_center_wrapped_nm"
        ],
        dtype=float,
    )

    axis = np.asarray(
        prototype[
            "tube_axis"
        ],
        dtype=float,
    )

    axis_norm = float(
        np.linalg.norm(axis)
    )

    if axis_norm <= 0.0:
        raise RuntimeError(
            "Invalid tube axis."
        )

    axis /= axis_norm

    axial_low = float(
        prototype[
            "axial_low_nm"
        ]
    )

    axial_high = float(
        prototype[
            "axial_high_nm"
        ]
    )

    accessible_radius = float(
        prototype[
            "accessible_radius_nm"
        ]
    )

    relative_coordinates = minimum_image(
        water_positions - center,
        box,
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
        (axial >= axial_low)
        & (axial <= axial_high)
        & (radial <= accessible_radius)
    )


def minimum_cap_water_distance(
    water_positions: np.ndarray,
    cap_positions: np.ndarray,
    box: np.ndarray,
) -> float:
    minimum_distance = math.inf
    chunk_size = 512

    for start in range(
        0,
        len(water_positions),
        chunk_size,
    ):
        stop = min(
            start + chunk_size,
            len(water_positions),
        )

        displacement = (
            water_positions[
                start:stop,
                None,
                :,
            ]
            - cap_positions[
                None,
                :,
                :,
            ]
        )

        displacement = minimum_image(
            displacement,
            box,
        )

        distances = np.linalg.norm(
            displacement,
            axis=2,
        )

        minimum_distance = min(
            minimum_distance,
            float(
                np.min(distances)
            ),
        )

    return minimum_distance


def parse_index_group_number(
    path: Path,
    target: str,
) -> int:
    group_number = -1

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        match = re.match(
            r"^\s*\[\s*([^\]]+?)\s*\]\s*$",
            raw_line,
        )

        if match is None:
            continue

        group_number += 1

        if (
            match.group(1).strip()
            == target
        ):
            return group_number

    raise RuntimeError(
        f"Index group {target!r} was not found."
    )


def frame_sort_key(
    path: Path,
) -> tuple[int, str]:
    match = re.search(
        r"(\d+)(?=\.gro$)",
        path.name,
    )

    if match is None:
        return (
            10**9,
            path.name,
        )

    return (
        int(
            match.group(1)
        ),
        path.name,
    )


def extract_frames(
    gmx: str,
) -> list[Path]:
    system_group = parse_index_group_number(
        LOCAL_INDEX,
        "System",
    )

    temporary_directory = tempfile.mkdtemp(
        prefix="r1_nvt20_frames_",
        dir="/tmp",
    )

    temporary_root = Path(
        temporary_directory
    )

    output_pattern = (
        temporary_root
        / "frame.gro"
    )

    completed = run_command(
        [
            gmx,
            "trjconv",
            "-s",
            str(LOCAL_TPR),
            "-f",
            str(XTC),
            "-n",
            str(LOCAL_INDEX),
            "-o",
            str(output_pattern),
            "-sep",
            "-pbc",
            "atom",
        ],
        cwd=temporary_root,
        input_text=(
            f"{system_group}\n"
        ),
    )

    TRJCONV_LOG.write_text(
        completed.stdout,
        encoding="utf-8",
    )

    if completed.returncode != 0:
        shutil.rmtree(
            temporary_root,
            ignore_errors=True,
        )

        raise RuntimeError(
            "gmx trjconv failed. "
            f"See {TRJCONV_LOG}"
        )

    frames = sorted(
        temporary_root.glob(
            "frame*.gro"
        ),
        key=frame_sort_key,
    )

    if not frames:
        shutil.rmtree(
            temporary_root,
            ignore_errors=True,
        )

        raise RuntimeError(
            "No trajectory frames were extracted."
        )

    return frames


def parse_energy_menu(
    output: str,
) -> dict[str, int]:
    terms: dict[str, int] = {}

    pattern = re.compile(
        r"(?<!\S)(\d+)\s+"
        r"([A-Za-z][A-Za-z0-9_.:+\-()]*)"
    )

    for line in output.splitlines():
        for match in pattern.finditer(line):
            terms[
                match.group(2)
            ] = int(
                match.group(1)
            )

    return terms


def probe_energy_menu(
    gmx: str,
) -> dict[str, int]:
    probe_output = (
        OUTPUT_ROOT
        / "energy_menu_probe.xvg"
    )

    completed = run_command(
        [
            gmx,
            "energy",
            "-f",
            str(EDR),
            "-o",
            str(probe_output),
        ],
        cwd=OUTPUT_ROOT,
        input_text="0\n",
    )

    ENERGY_MENU_LOG.write_text(
        completed.stdout,
        encoding="utf-8",
    )

    terms = parse_energy_menu(
        completed.stdout
    )

    if not terms:
        raise RuntimeError(
            "Could not parse the GROMACS energy menu."
        )

    return terms


def resolve_energy_term(
    terms: dict[str, int],
    *,
    exact: tuple[str, ...],
    required_tokens: tuple[str, ...],
) -> tuple[str, int]:
    for name in exact:
        if name in terms:
            return (
                name,
                terms[
                    name
                ],
            )

    for name, number in terms.items():
        if all(
            token in name
            for token in required_tokens
        ):
            return (
                name,
                number,
            )

    raise RuntimeError(
        "Could not resolve energy term. "
        f"Exact={exact}; tokens={required_tokens}"
    )


def extract_energy_series(
    gmx: str,
    term_number: int,
    output_path: Path,
) -> np.ndarray:
    completed = run_command(
        [
            gmx,
            "energy",
            "-f",
            str(EDR),
            "-o",
            str(output_path),
            "-xvg",
            "none",
            "-dp",
        ],
        cwd=OUTPUT_ROOT,
        input_text=(
            f"{term_number}\n"
            "0\n"
        ),
    )

    if (
        completed.returncode != 0
        or not output_path.exists()
        or output_path.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Could not extract energy term "
            f"{term_number}."
        )

    rows = []

    for raw_line in output_path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        line = raw_line.strip()

        if (
            not line
            or line.startswith("#")
            or line.startswith("@")
        ):
            continue

        fields = line.split()

        if len(fields) < 2:
            continue

        rows.append(
            (
                float(
                    fields[0]
                ),
                float(
                    fields[1]
                ),
            )
        )

    array = np.asarray(
        rows,
        dtype=float,
    )

    if (
        array.ndim != 2
        or array.shape[1] != 2
        or len(array) == 0
        or not np.all(
            np.isfinite(array)
        )
    ):
        raise RuntimeError(
            f"Invalid energy series in {output_path}"
        )

    return array


def align_energy_series(
    series: dict[str, np.ndarray],
) -> list[dict[str, Any]]:
    reference_name = next(
        iter(series)
    )

    reference_times = series[
        reference_name
    ][
        :,
        0
    ]

    for name, array in series.items():
        if (
            len(array)
            != len(reference_times)
            or not np.allclose(
                array[
                    :,
                    0
                ],
                reference_times,
                rtol=0.0,
                atol=1.0e-8,
            )
        ):
            raise RuntimeError(
                "Energy time grids do not match: "
                f"{reference_name} vs {name}"
            )

    rows = []

    for index, time_ps in enumerate(
        reference_times
    ):
        row: dict[str, Any] = {
            "time_ps": float(
                time_ps
            )
        }

        for name, array in series.items():
            row[
                name
            ] = float(
                array[
                    index,
                    1
                ]
            )

        rows.append(row)

    return rows


def linear_slope(
    times: np.ndarray,
    values: np.ndarray,
) -> float:
    if len(times) < 2:
        return math.nan

    coefficients = np.polyfit(
        times,
        values,
        1,
    )

    return float(
        coefficients[0]
    )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        SOURCE_TPR,
        SOURCE_GRO,
        SOURCE_INDEX,
        PREPARATION_SUMMARY,
        WARNING_SUMMARY,
        PROTOTYPE_JSON,
    ):
        require_file(required)

    preparation = read_single_csv_row(
        PREPARATION_SUMMARY
    )

    warning = read_single_csv_row(
        WARNING_SUMMARY
    )

    prototype = json.loads(
        PROTOTYPE_JSON.read_text(
            encoding="utf-8"
        )
    )

    if (
        preparation.get(
            "decision"
        )
        !=
        "R1_FROZEN_SOLUTE_NVT_20PS_PREPARED"
    ):
        raise RuntimeError(
            "The 20 ps R1 NVT preparation is not accepted."
        )

    if not parse_bool(
        preparation.get(
            "NVT_execution_authorized",
            "false",
        )
    ):
        raise RuntimeError(
            "The preparation did not authorize NVT execution."
        )

    if (
        warning.get(
            "decision"
        )
        !=
        "R1_NVT_GROMPP_WARNING_AUTHORIZED"
    ):
        raise RuntimeError(
            "The controlled grompp warning is not authorized."
        )

    if not parse_bool(
        warning.get(
            "NVT_execution_authorized",
            "false",
        )
    ):
        raise RuntimeError(
            "The warning audit did not authorize execution."
        )

    if CAP_STOP != EXPECTED_ATOMS:
        raise RuntimeError(
            "Internal atom accounting failed."
        )

    gmx = locate_gmx()

    shutil.copy2(
        SOURCE_TPR,
        LOCAL_TPR,
    )

    shutil.copy2(
        SOURCE_INDEX,
        LOCAL_INDEX,
    )

    for path in (
        FINAL_GRO,
        XTC,
        EDR,
        MDLOG,
        CPT,
        MDRUN_CONSOLE,
    ):
        if path.exists():
            path.unlink()

    print(
        "Starting R1 frozen-solute NVT 20 ps."
    )

    print(
        "HBN, PYR, and both caps remain frozen; "
        "only TIP4P/2005 water is mobile."
    )

    print(
        "Expected steps / duration / frames: "
        "40000 / 20.0 ps / 41"
    )

    mdrun_return_code = run_live(
        [
            gmx,
            "mdrun",
            "-s",
            str(LOCAL_TPR),
            "-deffnm",
            str(DEFFNM),
            "-ntmpi",
            "1",
            "-ntomp",
            "4",
        ],
        cwd=OUTPUT_ROOT,
        log_path=MDRUN_CONSOLE,
    )

    for required in (
        FINAL_GRO,
        XTC,
        EDR,
        MDLOG,
        CPT,
    ):
        require_file(required)

    trajectory_check = run_command(
        [
            gmx,
            "check",
            "-f",
            str(XTC),
        ],
        cwd=OUTPUT_ROOT,
    )

    TRAJECTORY_CHECK_LOG.write_text(
        trajectory_check.stdout,
        encoding="utf-8",
    )

    initial_title, initial_atoms, initial_box = (
        read_gro(
            SOURCE_GRO
        )
    )

    final_title, final_atoms, final_box = (
        read_gro(
            FINAL_GRO
        )
    )

    if len(initial_atoms) != EXPECTED_ATOMS:
        raise RuntimeError(
            "Unexpected initial atom count: "
            f"{len(initial_atoms)}/{EXPECTED_ATOMS}"
        )

    if len(final_atoms) != EXPECTED_ATOMS:
        raise RuntimeError(
            "Unexpected final atom count: "
            f"{len(final_atoms)}/{EXPECTED_ATOMS}"
        )

    initial_positions = atom_positions(
        initial_atoms
    )

    final_positions = atom_positions(
        final_atoms
    )

    oxygen_indices = water_oxygen_indices(
        initial_atoms
    )

    initial_water_positions = (
        initial_positions[
            oxygen_indices
        ]
    )

    final_water_positions = (
        final_positions[
            oxygen_indices
        ]
    )

    initial_cap_positions = (
        initial_positions[
            CAP_START:
            CAP_STOP
        ]
    )

    final_cap_positions = (
        final_positions[
            CAP_START:
            CAP_STOP
        ]
    )

    hbn_final_rms, hbn_final_max = (
        displacement_metrics(
            initial_positions[
                :HBN_ATOMS
            ],
            final_positions[
                :HBN_ATOMS
            ],
            initial_box,
        )
    )

    pyr_final_rms, pyr_final_max = (
        displacement_metrics(
            initial_positions[
                HBN_ATOMS:
                SOLUTE_ATOMS
            ],
            final_positions[
                HBN_ATOMS:
                SOLUTE_ATOMS
            ],
            initial_box,
        )
    )

    cap_final_rms, cap_final_max = (
        displacement_metrics(
            initial_cap_positions,
            final_cap_positions,
            initial_box,
        )
    )

    water_final_rms, water_final_max = (
        displacement_metrics(
            initial_water_positions,
            final_water_positions,
            initial_box,
        )
    )

    box_max_difference = float(
        np.max(
            np.abs(
                final_box
                - initial_box
            )
        )
    )

    initial_lumen_mask = lumen_mask(
        initial_water_positions,
        initial_box,
        prototype,
    )

    initial_lumen_count = int(
        np.count_nonzero(
            initial_lumen_mask
        )
    )

    frames = extract_frames(
        gmx
    )

    temporary_frame_root = (
        frames[
            0
        ].parent
    )

    occupancy_rows = []

    maximum_hbn_trajectory_displacement = 0.0
    maximum_pyr_trajectory_displacement = 0.0
    maximum_cap_trajectory_displacement = 0.0
    minimum_distance_across_trajectory = math.inf

    try:
        for frame_index, frame_path in enumerate(
            frames
        ):
            _, frame_atoms, frame_box = read_gro(
                frame_path
            )

            if len(frame_atoms) != EXPECTED_ATOMS:
                raise RuntimeError(
                    "Unexpected frame atom count: "
                    f"{frame_path.name} has "
                    f"{len(frame_atoms)} atoms."
                )

            frame_positions = atom_positions(
                frame_atoms
            )

            frame_water_positions = (
                frame_positions[
                    oxygen_indices
                ]
            )

            frame_cap_positions = (
                frame_positions[
                    CAP_START:
                    CAP_STOP
                ]
            )

            current_mask = lumen_mask(
                frame_water_positions,
                frame_box,
                prototype,
            )

            occupancy = int(
                np.count_nonzero(
                    current_mask
                )
            )

            retained_initial = int(
                np.count_nonzero(
                    current_mask
                    & initial_lumen_mask
                )
            )

            entered_from_outside = int(
                np.count_nonzero(
                    current_mask
                    & ~initial_lumen_mask
                )
            )

            lost_initial = (
                initial_lumen_count
                - retained_initial
            )

            cap_water_distance = (
                minimum_cap_water_distance(
                    frame_water_positions,
                    frame_cap_positions,
                    frame_box,
                )
            )

            _, hbn_frame_max = (
                displacement_metrics(
                    initial_positions[
                        :HBN_ATOMS
                    ],
                    frame_positions[
                        :HBN_ATOMS
                    ],
                    initial_box,
                )
            )

            _, pyr_frame_max = (
                displacement_metrics(
                    initial_positions[
                        HBN_ATOMS:
                        SOLUTE_ATOMS
                    ],
                    frame_positions[
                        HBN_ATOMS:
                        SOLUTE_ATOMS
                    ],
                    initial_box,
                )
            )

            _, cap_frame_max = (
                displacement_metrics(
                    initial_cap_positions,
                    frame_cap_positions,
                    initial_box,
                )
            )

            maximum_hbn_trajectory_displacement = max(
                maximum_hbn_trajectory_displacement,
                hbn_frame_max,
            )

            maximum_pyr_trajectory_displacement = max(
                maximum_pyr_trajectory_displacement,
                pyr_frame_max,
            )

            maximum_cap_trajectory_displacement = max(
                maximum_cap_trajectory_displacement,
                cap_frame_max,
            )

            minimum_distance_across_trajectory = min(
                minimum_distance_across_trajectory,
                cap_water_distance,
            )

            time_ps = (
                frame_index
                * EXPECTED_FRAME_INTERVAL_PS
            )

            occupancy_rows.append(
                {
                    "frame": frame_index,
                    "time_ps": time_ps,
                    "lumen_occupancy": occupancy,
                    "initial_lumen_waters_retained": (
                        retained_initial
                    ),
                    "initial_lumen_waters_lost": (
                        lost_initial
                    ),
                    "waters_entered_from_initial_outside": (
                        entered_from_outside
                    ),
                    "initial_lumen_retention_fraction": (
                        retained_initial
                        / initial_lumen_count
                    ),
                    "minimum_CAP_OW_distance_nm": (
                        cap_water_distance
                    ),
                    "HBN_max_displacement_nm": (
                        hbn_frame_max
                    ),
                    "PYR_max_displacement_nm": (
                        pyr_frame_max
                    ),
                    "CAPS_max_displacement_nm": (
                        cap_frame_max
                    ),
                }
            )

    finally:
        shutil.rmtree(
            temporary_frame_root,
            ignore_errors=True,
        )

    write_csv(
        OCCUPANCY_CSV,
        occupancy_rows,
    )

    frame_count = len(
        occupancy_rows
    )

    trajectory_times = np.asarray(
        [
            row[
                "time_ps"
            ]
            for row in occupancy_rows
        ],
        dtype=float,
    )

    occupancies = np.asarray(
        [
            row[
                "lumen_occupancy"
            ]
            for row in occupancy_rows
        ],
        dtype=float,
    )

    retained_initial_series = np.asarray(
        [
            row[
                "initial_lumen_waters_retained"
            ]
            for row in occupancy_rows
        ],
        dtype=float,
    )

    occupancy_mean = float(
        np.mean(
            occupancies
        )
    )

    occupancy_std = float(
        np.std(
            occupancies,
            ddof=1,
        )
        if len(occupancies) > 1
        else 0.0
    )

    occupancy_min = int(
        np.min(
            occupancies
        )
    )

    occupancy_max = int(
        np.max(
            occupancies
        )
    )

    occupancy_endpoint = int(
        occupancies[
            -1
        ]
    )

    zero_occupancy_fraction = float(
        np.mean(
            occupancies == 0.0
        )
    )

    endpoint_initial_retention = int(
        retained_initial_series[
            -1
        ]
    )

    endpoint_initial_retention_fraction = (
        endpoint_initial_retention
        / initial_lumen_count
    )

    half_start = len(
        occupancies
    ) // 2

    second_half_occupancy_mean = float(
        np.mean(
            occupancies[
                half_start:
            ]
        )
    )

    second_half_occupancy_slope = (
        linear_slope(
            trajectory_times[
                half_start:
            ],
            occupancies[
                half_start:
            ],
        )
    )

    terms = probe_energy_menu(
        gmx
    )

    temperature_name, temperature_number = (
        resolve_energy_term(
            terms,
            exact=(
                "Temperature",
            ),
            required_tokens=(
                "Temperature",
            ),
        )
    )

    potential_name, potential_number = (
        resolve_energy_term(
            terms,
            exact=(
                "Potential",
            ),
            required_tokens=(
                "Potential",
            ),
        )
    )

    total_energy_name, total_energy_number = (
        resolve_energy_term(
            terms,
            exact=(
                "Total-Energy",
                "Total_Energy",
            ),
            required_tokens=(
                "Total",
                "Energy",
            ),
        )
    )

    cap_sol_name, cap_sol_number = (
        resolve_energy_term(
            terms,
            exact=(
                "LJ-SR:CAPS-SOL",
                "LJ-SR:SOL-CAPS",
            ),
            required_tokens=(
                "LJ-SR",
                "CAPS",
                "SOL",
            ),
        )
    )

    temperature_data = extract_energy_series(
        gmx,
        temperature_number,
        TEMPERATURE_XVG,
    )

    potential_data = extract_energy_series(
        gmx,
        potential_number,
        POTENTIAL_XVG,
    )

    total_energy_data = extract_energy_series(
        gmx,
        total_energy_number,
        TOTAL_ENERGY_XVG,
    )

    cap_sol_data = extract_energy_series(
        gmx,
        cap_sol_number,
        CAP_SOL_LJ_XVG,
    )

    energy_rows = align_energy_series(
        {
            "temperature_K": (
                temperature_data
            ),
            "potential_kJ_mol": (
                potential_data
            ),
            "total_energy_kJ_mol": (
                total_energy_data
            ),
            "CAP_SOL_LJ_kJ_mol": (
                cap_sol_data
            ),
        }
    )

    write_csv(
        ENERGY_CSV,
        energy_rows,
    )

    temperatures = temperature_data[
        :,
        1
    ]

    potentials = potential_data[
        :,
        1
    ]

    total_energies = total_energy_data[
        :,
        1
    ]

    cap_sol_energies = cap_sol_data[
        :,
        1
    ]

    temperature_mean = float(
        np.mean(
            temperatures
        )
    )

    temperature_std = float(
        np.std(
            temperatures,
            ddof=1,
        )
        if len(temperatures) > 1
        else 0.0
    )

    temperature_min = float(
        np.min(
            temperatures
        )
    )

    temperature_max = float(
        np.max(
            temperatures
        )
    )

    potential_initial = float(
        potentials[
            0
        ]
    )

    potential_final = float(
        potentials[
            -1
        ]
    )

    potential_change = (
        potential_final
        - potential_initial
    )

    cap_sol_initial = float(
        cap_sol_energies[
            0
        ]
    )

    cap_sol_final = float(
        cap_sol_energies[
            -1
        ]
    )

    cap_sol_maximum = float(
        np.max(
            cap_sol_energies
        )
    )

    total_energy_initial = float(
        total_energies[
            0
        ]
    )

    total_energy_final = float(
        total_energies[
            -1
        ]
    )

    combined_log = (
        MDLOG.read_text(
            encoding="utf-8",
            errors="replace",
        )
        + "\n"
        + MDRUN_CONSOLE.read_text(
            encoding="utf-8",
            errors="replace",
        )
    )

    instability_hits = []

    for pattern in INSTABILITY_PATTERNS:
        if re.search(
            pattern,
            combined_log,
            flags=re.IGNORECASE,
        ):
            instability_hits.append(
                pattern
            )

    finished = bool(
        re.search(
            r"Finished mdrun",
            combined_log,
            flags=re.IGNORECASE,
        )
    )

    minimum_occupancy_required = int(
        math.ceil(
            MIN_OCCUPANCY_FRACTION
            * initial_lumen_count
        )
    )

    minimum_initial_retention_required = int(
        math.ceil(
            MIN_INITIAL_WATER_RETENTION_FRACTION
            * initial_lumen_count
        )
    )

    gates = {
        "mdrun_return_code_zero": (
            mdrun_return_code == 0
        ),
        "mdrun_finished": finished,
        "no_instability_signatures": (
            len(
                instability_hits
            )
            == 0
        ),
        "trajectory_check_return_code_zero": (
            trajectory_check.returncode == 0
        ),
        "final_atom_count_is_68314": (
            len(
                final_atoms
            )
            == EXPECTED_ATOMS
        ),
        "trajectory_has_41_frames": (
            frame_count
            == EXPECTED_FRAMES
        ),
        "trajectory_starts_at_0ps": (
            math.isclose(
                float(
                    trajectory_times[
                        0
                    ]
                ),
                0.0,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
        ),
        "trajectory_ends_at_20ps": (
            math.isclose(
                float(
                    trajectory_times[
                        -1
                    ]
                ),
                EXPECTED_DURATION_PS,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
        ),
        "box_is_unchanged": (
            box_max_difference
            <= 1.0e-6
        ),
        "HBN_final_coordinates_are_frozen": (
            hbn_final_max
            <= FROZEN_FINAL_TOLERANCE_NM
        ),
        "PYR_final_coordinates_are_frozen": (
            pyr_final_max
            <= FROZEN_FINAL_TOLERANCE_NM
        ),
        "CAPS_final_coordinates_are_frozen": (
            cap_final_max
            <= FROZEN_FINAL_TOLERANCE_NM
        ),
        "HBN_trajectory_is_frozen": (
            maximum_hbn_trajectory_displacement
            <= FROZEN_TRAJECTORY_TOLERANCE_NM
        ),
        "PYR_trajectory_is_frozen": (
            maximum_pyr_trajectory_displacement
            <= FROZEN_TRAJECTORY_TOLERANCE_NM
        ),
        "CAPS_trajectory_is_frozen": (
            maximum_cap_trajectory_displacement
            <= FROZEN_TRAJECTORY_TOLERANCE_NM
        ),
        "temperature_series_is_finite": (
            np.all(
                np.isfinite(
                    temperatures
                )
            )
        ),
        "temperature_mean_is_295_to_305K": (
            MIN_ACCEPTABLE_TEMPERATURE_K
            <= temperature_mean
            <= MAX_ACCEPTABLE_TEMPERATURE_K
        ),
        "temperature_standard_deviation_is_acceptable": (
            temperature_std
            <= MAX_TEMPERATURE_STANDARD_DEVIATION_K
        ),
        "potential_series_is_finite": (
            np.all(
                np.isfinite(
                    potentials
                )
            )
        ),
        "total_energy_series_is_finite": (
            np.all(
                np.isfinite(
                    total_energies
                )
            )
        ),
        "CAP_SOL_energy_series_is_finite": (
            np.all(
                np.isfinite(
                    cap_sol_energies
                )
            )
        ),
        "initial_lumen_occupancy_is_428": (
            initial_lumen_count
            == INITIAL_LUMEN_WATERS
        ),
        "no_complete_lumen_drying": (
            zero_occupancy_fraction
            == 0.0
        ),
        "minimum_lumen_occupancy_is_at_least_90_percent": (
            occupancy_min
            >= minimum_occupancy_required
        ),
        "endpoint_lumen_occupancy_is_at_least_90_percent": (
            occupancy_endpoint
            >= minimum_occupancy_required
        ),
        "endpoint_initial_lumen_retention_is_at_least_90_percent": (
            endpoint_initial_retention
            >= minimum_initial_retention_required
        ),
        "CAP_OW_distance_remains_above_0p15nm": (
            minimum_distance_across_trajectory
            >= MIN_CAP_WATER_DISTANCE_NM
        ),
        "water_coordinates_are_mobile": (
            water_final_rms
            > 0.01
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
        "R1_FROZEN_SOLUTE_NVT_20PS_VALIDATED"
        if accepted
        else
        "R1_FROZEN_SOLUTE_NVT_20PS_REQUIRES_REVIEW"
    )

    required_next_step = (
        "ASSESS_R1_EXTENSION_TO_50PS_FROZEN_SOLUTE_SCREENING"
        if accepted
        else
        "REVIEW_R1_NVT_20PS_GATE_FAILURES"
    )

    summary = {
        "decision": decision,
        "mdrun_return_code": (
            mdrun_return_code
        ),
        "mdrun_finished": finished,
        "trajectory_check_return_code": (
            trajectory_check.returncode
        ),
        "trajectory_frames": frame_count,
        "trajectory_start_ps": float(
            trajectory_times[
                0
            ]
        ),
        "trajectory_end_ps": float(
            trajectory_times[
                -1
            ]
        ),
        "temperature_term": (
            temperature_name
        ),
        "temperature_mean_K": (
            temperature_mean
        ),
        "temperature_std_K": (
            temperature_std
        ),
        "temperature_min_K": (
            temperature_min
        ),
        "temperature_max_K": (
            temperature_max
        ),
        "potential_term": (
            potential_name
        ),
        "potential_initial_kJ_mol": (
            potential_initial
        ),
        "potential_final_kJ_mol": (
            potential_final
        ),
        "potential_change_kJ_mol": (
            potential_change
        ),
        "total_energy_term": (
            total_energy_name
        ),
        "total_energy_initial_kJ_mol": (
            total_energy_initial
        ),
        "total_energy_final_kJ_mol": (
            total_energy_final
        ),
        "CAP_SOL_LJ_term": (
            cap_sol_name
        ),
        "CAP_SOL_LJ_initial_kJ_mol": (
            cap_sol_initial
        ),
        "CAP_SOL_LJ_final_kJ_mol": (
            cap_sol_final
        ),
        "CAP_SOL_LJ_maximum_kJ_mol": (
            cap_sol_maximum
        ),
        "HBN_final_RMS_displacement_nm": (
            hbn_final_rms
        ),
        "HBN_final_max_displacement_nm": (
            hbn_final_max
        ),
        "PYR_final_RMS_displacement_nm": (
            pyr_final_rms
        ),
        "PYR_final_max_displacement_nm": (
            pyr_final_max
        ),
        "CAPS_final_RMS_displacement_nm": (
            cap_final_rms
        ),
        "CAPS_final_max_displacement_nm": (
            cap_final_max
        ),
        "HBN_max_trajectory_displacement_nm": (
            maximum_hbn_trajectory_displacement
        ),
        "PYR_max_trajectory_displacement_nm": (
            maximum_pyr_trajectory_displacement
        ),
        "CAPS_max_trajectory_displacement_nm": (
            maximum_cap_trajectory_displacement
        ),
        "waterO_final_RMS_displacement_nm": (
            water_final_rms
        ),
        "waterO_final_max_displacement_nm": (
            water_final_max
        ),
        "box_max_difference_nm": (
            box_max_difference
        ),
        "initial_lumen_occupancy": (
            initial_lumen_count
        ),
        "mean_lumen_occupancy": (
            occupancy_mean
        ),
        "std_lumen_occupancy": (
            occupancy_std
        ),
        "minimum_lumen_occupancy": (
            occupancy_min
        ),
        "maximum_lumen_occupancy": (
            occupancy_max
        ),
        "endpoint_lumen_occupancy": (
            occupancy_endpoint
        ),
        "zero_occupancy_fraction": (
            zero_occupancy_fraction
        ),
        "endpoint_initial_lumen_waters_retained": (
            endpoint_initial_retention
        ),
        "endpoint_initial_lumen_retention_fraction": (
            endpoint_initial_retention_fraction
        ),
        "second_half_mean_lumen_occupancy": (
            second_half_occupancy_mean
        ),
        "second_half_occupancy_slope_waters_per_ps": (
            second_half_occupancy_slope
        ),
        "minimum_CAP_OW_distance_nm": (
            minimum_distance_across_trajectory
        ),
        "instability_signature_count": (
            len(
                instability_hits
            )
        ),
        "instability_signatures": (
            " | ".join(
                instability_hits
            )
        ),
        "failed_gates": (
            " | ".join(
                failed_gates
            )
        ),
        "extension_preparation_authorized": (
            accepted
        ),
        "long_mobile_MD_authorized": False,
        "multitemperature_MD_authorized": False,
        "QM_recalculation_authorized": False,
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

    gate_rows = [
        {
            "gate": name,
            "pass": passed,
        }
        for name, passed
        in gates.items()
    ]

    write_csv(
        GATE_CSV,
        gate_rows,
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
        f"""# R1 Frozen-Solute NVT 20 ps Validation

## Scope

The R1 fully capped steric positive control was simulated for
20 ps at 300 K.

Frozen throughout the trajectory:

- HBN;
- all pyrene chromophores;
- both steric caps.

Mobile and thermostatted:

- {EXPECTED_WATERS} TIP4P/2005 water molecules.

This is a confinement-methodology screening, not a chemically
realizable final device model.

## Execution

- Mdrun return code:
  **{mdrun_return_code}**
- Finished mdrun:
  **{'YES' if finished else 'NO'}**
- Frames:
  **{frame_count}/{EXPECTED_FRAMES}**
- Time range:
  **{trajectory_times[0]:.1f}–{trajectory_times[-1]:.1f} ps**
- Instability signatures:
  **{'NONE' if not instability_hits else ' | '.join(instability_hits)}**

## Temperature and energy

- Temperature mean ± standard deviation:
  **{temperature_mean:.4f} ± {temperature_std:.4f} K**
- Temperature minimum/maximum:
  **{temperature_min:.4f}/{temperature_max:.4f} K**
- Potential initial/final/change:
  **{potential_initial:.6f}/
  {potential_final:.6f}/
  {potential_change:.6f} kJ/mol**
- Total energy initial/final:
  **{total_energy_initial:.6f}/
  {total_energy_final:.6f} kJ/mol**
- CAP–SOL LJ initial/final/maximum:
  **{cap_sol_initial:.6f}/
  {cap_sol_final:.6f}/
  {cap_sol_maximum:.6f} kJ/mol**

## Frozen-group integrity

Final RMS/max displacement:

- HBN:
  **{hbn_final_rms:.12f}/{hbn_final_max:.12f} nm**
- PYR:
  **{pyr_final_rms:.12f}/{pyr_final_max:.12f} nm**
- CAPS:
  **{cap_final_rms:.12f}/{cap_final_max:.12f} nm**

Maximum displacement observed in compressed trajectory:

- HBN:
  **{maximum_hbn_trajectory_displacement:.6f} nm**
- PYR:
  **{maximum_pyr_trajectory_displacement:.6f} nm**
- CAPS:
  **{maximum_cap_trajectory_displacement:.6f} nm**

## Water motion and confinement

- Final water-O RMS/max displacement:
  **{water_final_rms:.6f}/{water_final_max:.6f} nm**
- Initial lumen occupancy:
  **{initial_lumen_count} waters**
- Mean ± standard deviation:
  **{occupancy_mean:.4f} ± {occupancy_std:.4f} waters**
- Minimum/maximum occupancy:
  **{occupancy_min}/{occupancy_max} waters**
- Endpoint occupancy:
  **{occupancy_endpoint} waters**
- Zero-occupancy fraction:
  **{zero_occupancy_fraction:.6f}**
- Initial lumen waters retained at endpoint:
  **{endpoint_initial_retention}/{initial_lumen_count}
  ({endpoint_initial_retention_fraction:.6f})**
- Second-half mean occupancy:
  **{second_half_occupancy_mean:.4f} waters**
- Second-half occupancy slope:
  **{second_half_occupancy_slope:.6f} waters/ps**
- Minimum CAP–OW distance:
  **{minimum_distance_across_trajectory:.6f} nm**

## Validation gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Preparation of a 50 ps frozen-solute extension authorized:
  **{'YES' if accepted else 'NO'}**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `{required_next_step}`

The R1 caps remain a neutral frozen steric positive control. A positive
20 ps result demonstrates that the retention-analysis methodology can
detect a closed-boundary confinement condition; it does not establish
chemical realizability or justify a final device architecture.
""",
        encoding="utf-8",
    )

    print()
    print(
        "Day023 R1 frozen-solute NVT 20 ps "
        "execution and validation completed."
    )

    print(
        "Mdrun / trajectory-check return codes: "
        f"{mdrun_return_code}/"
        f"{trajectory_check.returncode}"
    )

    print(
        "Finished / instability signatures: "
        f"{'YES' if finished else 'NO'} / "
        + (
            "NONE"
            if not instability_hits
            else " | ".join(
                instability_hits
            )
        )
    )

    print(
        "Frames / time range: "
        f"{frame_count}/"
        f"{trajectory_times[0]:.1f}-"
        f"{trajectory_times[-1]:.1f} ps"
    )

    print(
        "Temperature mean/std/min/max: "
        f"{temperature_mean:.4f}/"
        f"{temperature_std:.4f}/"
        f"{temperature_min:.4f}/"
        f"{temperature_max:.4f} K"
    )

    print(
        "Potential initial/final/change: "
        f"{potential_initial:.6f}/"
        f"{potential_final:.6f}/"
        f"{potential_change:.6f} kJ/mol"
    )

    print(
        "CAP-SOL LJ initial/final/maximum: "
        f"{cap_sol_initial:.6f}/"
        f"{cap_sol_final:.6f}/"
        f"{cap_sol_maximum:.6f} kJ/mol"
    )

    print(
        "HBN final RMS/max displacement: "
        f"{hbn_final_rms:.12f}/"
        f"{hbn_final_max:.12f} nm"
    )

    print(
        "PYR final RMS/max displacement: "
        f"{pyr_final_rms:.12f}/"
        f"{pyr_final_max:.12f} nm"
    )

    print(
        "CAPS final RMS/max displacement: "
        f"{cap_final_rms:.12f}/"
        f"{cap_final_max:.12f} nm"
    )

    print(
        "Maximum trajectory displacement "
        "HBN/PYR/CAPS: "
        f"{maximum_hbn_trajectory_displacement:.6f}/"
        f"{maximum_pyr_trajectory_displacement:.6f}/"
        f"{maximum_cap_trajectory_displacement:.6f} nm"
    )

    print(
        "Water-O final RMS/max displacement: "
        f"{water_final_rms:.6f}/"
        f"{water_final_max:.6f} nm"
    )

    print(
        "Lumen occupancy "
        "initial/mean/min/max/endpoint: "
        f"{initial_lumen_count}/"
        f"{occupancy_mean:.4f}/"
        f"{occupancy_min}/"
        f"{occupancy_max}/"
        f"{occupancy_endpoint}"
    )

    print(
        "Endpoint initial-lumen retention: "
        f"{endpoint_initial_retention}/"
        f"{initial_lumen_count} "
        f"({endpoint_initial_retention_fraction:.6f})"
    )

    print(
        "Second-half occupancy mean/slope: "
        f"{second_half_occupancy_mean:.4f}/"
        f"{second_half_occupancy_slope:.6f} "
        "waters/ps"
    )

    print(
        "Minimum CAP-OW distance: "
        f"{minimum_distance_across_trajectory:.6f} nm"
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
        "50 ps frozen-solute extension preparation "
        "authorized: "
        f"{'YES' if accepted else 'NO'}"
    )

    print(
        "Long mobile MD authorized: NO"
    )

    print(
        "Multitemperature MD authorized: NO"
    )

    print(
        "QM recalculation authorized: NO"
    )

    print(
        f"Required next step: {required_next_step}"
    )

    print(
        f"Wrote: {relative(FINAL_GRO)}"
    )

    print(
        f"Wrote: {relative(XTC)}"
    )

    print(
        f"Wrote: {relative(EDR)}"
    )

    print(
        f"Wrote: {relative(OCCUPANCY_CSV)}"
    )

    print(
        f"Wrote: {relative(ENERGY_CSV)}"
    )

    print(
        f"Wrote: {relative(SUMMARY_CSV)}"
    )

    print(
        f"Wrote: {relative(GATE_CSV)}"
    )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )

    if not accepted:
        raise RuntimeError(
            "R1 frozen-solute 20 ps NVT "
            "requires review."
        )


if __name__ == "__main__":
    main()
