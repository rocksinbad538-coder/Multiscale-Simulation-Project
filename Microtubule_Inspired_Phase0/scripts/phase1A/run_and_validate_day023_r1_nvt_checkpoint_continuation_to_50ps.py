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

PREPARATION20_ROOT = (
    DAY023_ROOT
    / "07_r1_frozen_solute_nvt_20ps_preparation"
)

NVT20_ROOT = (
    DAY023_ROOT
    / "08_r1_frozen_solute_nvt_20ps"
)

THERMAL_REVIEW_ROOT = (
    DAY023_ROOT
    / "09_r1_nvt_20ps_thermal_review"
)

CONTINUATION_PREPARATION_ROOT = (
    DAY023_ROOT
    / "10_r1_nvt_50ps_checkpoint_continuation_preparation"
)

PROTOTYPE_ROOT = (
    DAY023_ROOT
    / "02_r1_steric_cap_prototype"
)

OUTPUT_ROOT = (
    DAY023_ROOT
    / "11_r1_frozen_solute_nvt_20_to_50ps"
)

EXTENDED_TPR = (
    CONTINUATION_PREPARATION_ROOT
    / "r1_frozen_solute_nvt_to_50ps.tpr"
)

SOURCE_CHECKPOINT = (
    CONTINUATION_PREPARATION_ROOT
    / "r1_20ps_source_checkpoint.cpt"
)

CONTINUATION_CONTRACT = (
    CONTINUATION_PREPARATION_ROOT
    / "r1_50ps_continuation_run_contract.json"
)

CONTINUATION_PREPARATION_SUMMARY = (
    CONTINUATION_PREPARATION_ROOT
    / "r1_50ps_continuation_preparation_summary.csv"
)

THERMAL_REVIEW_SUMMARY = (
    THERMAL_REVIEW_ROOT
    / "r1_nvt_20ps_thermal_review_summary.csv"
)

SOURCE_INITIAL_GRO = (
    PREPARATION20_ROOT
    / "r1_frozen_solute_nvt_20ps_input.gro"
)

SOURCE_INDEX = (
    PREPARATION20_ROOT
    / "r1_frozen_solute_nvt_20ps.ndx"
)

SOURCE_20PS_FINAL_GRO = (
    NVT20_ROOT
    / "r1_frozen_solute_nvt_20ps.gro"
)

SOURCE_20PS_XTC = (
    NVT20_ROOT
    / "r1_frozen_solute_nvt_20ps.xtc"
)

SOURCE_20PS_OCCUPANCY = (
    NVT20_ROOT
    / "r1_frozen_solute_nvt_20ps_lumen_occupancy.csv"
)

PROTOTYPE_JSON = (
    PROTOTYPE_ROOT
    / "r1_selected_steric_cap_definition.json"
)

LOCAL_INDEX = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20_to_50ps.ndx"
)

DEFFNM = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20_to_50ps"
)

MDRUN_CONSOLE = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20_to_50ps_mdrun_console.log"
)

TRAJECTORY_CHECK_LOG = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20_to_50ps_trajectory_check.log"
)

TRJCONV_LOG = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20_to_50ps_trjconv.log"
)

ENERGY_MENU_LOG = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20_to_50ps_energy_menu.txt"
)

TEMPERATURE_XVG = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20_to_50ps_temperature.xvg"
)

POTENTIAL_XVG = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20_to_50ps_potential.xvg"
)

TOTAL_ENERGY_XVG = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20_to_50ps_total_energy.xvg"
)

CAP_SOL_LJ_XVG = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20_to_50ps_cap_sol_lj.xvg"
)

CONTINUATION_OCCUPANCY_CSV = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20_to_50ps_lumen_occupancy.csv"
)

COMBINED_OCCUPANCY_CSV = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_0_to_50ps_combined_occupancy.csv"
)

CONTINUATION_ENERGY_CSV = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20_to_50ps_energy_series.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_50ps_summary.csv"
)

GATE_CSV = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_50ps_gates.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R1_FROZEN_SOLUTE_NVT_50PS_VALIDATION_DAY023.md"
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

CONTINUATION_START_PS = 20.0
CONTINUATION_END_PS = 50.0
FRAME_INTERVAL_PS = 0.5

ALLOWED_CONTINUATION_FRAME_COUNTS = {
    60,
    61,
}

EXPECTED_COMBINED_UNIQUE_FRAMES = 101

FROZEN_FINAL_TOLERANCE_NM = 1.0e-6
FROZEN_TRAJECTORY_TOLERANCE_NM = 1.5e-3

MIN_TEMPERATURE_MEAN_K = 295.0
MAX_TEMPERATURE_MEAN_K = 305.0
MAX_TEMPERATURE_STD_K = 5.0

MIN_INSTANTANEOUS_TEMPERATURE_K = 280.0
MAX_INSTANTANEOUS_TEMPERATURE_K = 320.0

MAX_TEMPERATURE_SLOPE_K_PER_PS = 0.20

MIN_OCCUPANCY_FRACTION = 0.90
MIN_CAP_WATER_DISTANCE_NM = 0.15
MAX_SECOND_HALF_OCCUPANCY_SLOPE = 0.25

MAX_CAP_SOL_LJ_KJ_MOL = 500.0

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


def read_csv_rows(
    path: Path,
) -> list[dict[str, str]]:
    require_file(path)

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(
            csv.DictReader(handle)
        )

    if not rows:
        raise RuntimeError(
            f"No data rows in {path}"
        )

    return rows


def read_single_csv_row(
    path: Path,
) -> dict[str, str]:
    rows = read_csv_rows(path)

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
            "Invalid simulation box."
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
                f"Malformed atom line "
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


def positions(
    atoms: list[dict[str, Any]],
) -> np.ndarray:
    return np.asarray(
        [
            atom["position"]
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
            if atom["atomname"] == "OW"
        ]

        if len(matches) != 1:
            raise RuntimeError(
                "Could not identify exactly one OW "
                f"in water {water_index}."
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

    axis /= np.linalg.norm(axis)

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


def numeric_frame_key(
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
    trajectory: Path,
) -> tuple[
    list[Path],
    Path,
]:
    system_group = parse_index_group_number(
        LOCAL_INDEX,
        "System",
    )

    temporary_root = Path(
        tempfile.mkdtemp(
            prefix="r1_continuation_frames_",
            dir="/tmp",
        )
    )

    completed = run_command(
        [
            gmx,
            "trjconv",
            "-s",
            str(EXTENDED_TPR),
            "-f",
            str(trajectory),
            "-n",
            str(LOCAL_INDEX),
            "-o",
            str(
                temporary_root
                / "frame.gro"
            ),
            "-sep",
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
        key=numeric_frame_key,
    )

    if not frames:
        shutil.rmtree(
            temporary_root,
            ignore_errors=True,
        )

        raise RuntimeError(
            "No continuation frames were extracted."
        )

    return (
        frames,
        temporary_root,
    )


def parse_check_times(
    output: str,
) -> tuple[
    float | None,
    float | None,
]:
    matches = re.findall(
        r"(?:Reading frame|Last frame)\s+\d+\s+time\s+"
        r"([-+0-9.eE]+)",
        output,
    )

    if not matches:
        return (
            None,
            None,
        )

    return (
        float(matches[0]),
        float(matches[-1]),
    )


def find_product(
    suffix: str,
    *,
    exclude: tuple[Path, ...] = (),
) -> Path:
    prefix = DEFFNM.name

    candidates = [
        path
        for path in OUTPUT_ROOT.glob(
            f"{prefix}*{suffix}"
        )
        if path not in exclude
        and "_prev" not in path.name
        and path.stat().st_size > 0
    ]

    if len(candidates) != 1:
        raise RuntimeError(
            f"Expected exactly one continuation {suffix} file; "
            f"found {[path.name for path in candidates]}"
        )

    return candidates[0]


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
    edr: Path,
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
            str(edr),
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
            "Could not parse the energy menu."
        )

    return terms


def resolve_energy_term(
    terms: dict[str, int],
    exact: tuple[str, ...],
    tokens: tuple[str, ...],
) -> tuple[
    str,
    int,
]:
    for name in exact:
        if name in terms:
            return (
                name,
                terms[name],
            )

    for name, number in terms.items():
        if all(
            token in name
            for token in tokens
        ):
            return (
                name,
                number,
            )

    raise RuntimeError(
        "Could not resolve energy term. "
        f"Exact={exact}; tokens={tokens}"
    )


def extract_energy_series(
    gmx: str,
    edr: Path,
    term_number: int,
    output_path: Path,
) -> np.ndarray:
    completed = run_command(
        [
            gmx,
            "energy",
            "-f",
            str(edr),
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
            f"Could not extract energy term {term_number}."
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
                float(fields[0]),
                float(fields[1]),
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
            len(array) != len(reference_times)
            or not np.allclose(
                array[:, 0],
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
            row[name] = float(
                array[
                    index,
                    1,
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

    return float(
        np.polyfit(
            times,
            values,
            1,
        )[0]
    )


def statistics(
    times: np.ndarray,
    values: np.ndarray,
) -> dict[str, float]:
    return {
        "mean": float(
            np.mean(values)
        ),
        "std": float(
            np.std(
                values,
                ddof=1,
            )
        ),
        "min": float(
            np.min(values)
        ),
        "max": float(
            np.max(values)
        ),
        "slope": linear_slope(
            times,
            values,
        ),
    }


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    required_files = (
        EXTENDED_TPR,
        SOURCE_CHECKPOINT,
        CONTINUATION_CONTRACT,
        CONTINUATION_PREPARATION_SUMMARY,
        THERMAL_REVIEW_SUMMARY,
        SOURCE_INITIAL_GRO,
        SOURCE_INDEX,
        SOURCE_20PS_FINAL_GRO,
        SOURCE_20PS_XTC,
        SOURCE_20PS_OCCUPANCY,
        PROTOTYPE_JSON,
    )

    for required in required_files:
        require_file(required)

    contract = json.loads(
        CONTINUATION_CONTRACT.read_text(
            encoding="utf-8"
        )
    )

    preparation = read_single_csv_row(
        CONTINUATION_PREPARATION_SUMMARY
    )

    thermal_review = read_single_csv_row(
        THERMAL_REVIEW_SUMMARY
    )

    prototype = json.loads(
        PROTOTYPE_JSON.read_text(
            encoding="utf-8"
        )
    )

    if (
        preparation.get("decision")
        !=
        "R1_CHECKPOINT_CONTINUATION_TO_50PS_PREPARED"
    ):
        raise RuntimeError(
            "The continuation preparation is not accepted."
        )

    if not parse_bool(
        preparation.get(
            "continuation_execution_authorized",
            "false",
        )
    ):
        raise RuntimeError(
            "Continuation execution is not authorized."
        )

    if (
        thermal_review.get("decision")
        !=
        "R1_INITIAL_THERMALIZATION_TRANSIENT_CONFIRMED"
    ):
        raise RuntimeError(
            "The thermal-transient review is not accepted."
        )

    if not bool(
        contract.get(
            "execution_authorized",
            False,
        )
    ):
        raise RuntimeError(
            "The run contract does not authorize execution."
        )

    if bool(
        contract.get(
            "velocity_regeneration",
            True,
        )
    ):
        raise RuntimeError(
            "The contract unexpectedly allows velocity regeneration."
        )

    if int(
        contract.get(
            "expected_remaining_steps",
            -1,
        )
    ) != 60000:
        raise RuntimeError(
            "Unexpected remaining-step contract."
        )

    if not math.isclose(
        float(
            contract.get(
                "target_total_time_ps",
                math.nan,
            )
        ),
        CONTINUATION_END_PS,
        rel_tol=0.0,
        abs_tol=1.0e-10,
    ):
        raise RuntimeError(
            "Unexpected target-time contract."
        )

    if CAP_STOP != EXPECTED_ATOMS:
        raise RuntimeError(
            "Internal atom accounting failed."
        )

    existing_products = [
        path
        for path in OUTPUT_ROOT.glob(
            f"{DEFFNM.name}*"
        )
        if path != MDRUN_CONSOLE
    ]

    reuse_completed_run = (
        os.environ.get(
            "R1_REUSE_COMPLETED_RUN",
            "0",
        )
        == "1"
    )

    if existing_products and not reuse_completed_run:
        raise RuntimeError(
            "Continuation output products already exist. "
            "Set R1_REUSE_COMPLETED_RUN=1 only after confirming "
            "that the checkpoint continuation finished successfully: "
            + ", ".join(
                path.name
                for path in existing_products
            )
        )

    if reuse_completed_run and not existing_products:
        raise RuntimeError(
            "R1_REUSE_COMPLETED_RUN=1 was requested, "
            "but no completed continuation products were found."
        )

    shutil.copy2(
        SOURCE_INDEX,
        LOCAL_INDEX,
    )

    gmx = locate_gmx()

    if reuse_completed_run:
        require_file(
            MDRUN_CONSOLE
        )

        console_text = MDRUN_CONSOLE.read_text(
            encoding="utf-8",
            errors="replace",
        )

        if (
            "Writing final coordinates"
            not in console_text
        ):
            raise RuntimeError(
                "The existing run does not confirm final-coordinate "
                "output and cannot be reused."
            )

        if re.search(
            r"fatal error|segmentation fault|lincs warning|"
            r"constraint warning|\\bnan\\b",
            console_text,
            flags=re.IGNORECASE,
        ):
            raise RuntimeError(
                "The existing run contains an instability signature "
                "and cannot be reused."
            )

        print(
            "Reusing the completed R1 checkpoint continuation."
        )

        print(
            "No molecular-dynamics steps will be repeated."
        )

        print(
            "Proceeding directly to trajectory, energy, "
            "and confinement validation."
        )

        mdrun_return_code = 0

    else:
        print(
            "Starting R1 checkpoint continuation from "
            "20 to 50 ps."
        )

        print(
            "The exact 20 ps checkpoint will be used; "
            "velocities will not be regenerated."
        )

        print(
            "Expected remaining steps / duration: "
            "60000 / 30.0 ps"
        )

        mdrun_return_code = run_live(
            [
                gmx,
                "mdrun",
                "-s",
                str(EXTENDED_TPR),
                "-cpi",
                str(SOURCE_CHECKPOINT),
                "-deffnm",
                str(DEFFNM),
                "-noappend",
                "-ntmpi",
                "1",
                "-ntomp",
                "4",
            ],
            cwd=OUTPUT_ROOT,
            log_path=MDRUN_CONSOLE,
        )

    continuation_xtc = find_product(
        ".xtc"
    )

    continuation_edr = find_product(
        ".edr"
    )

    continuation_log = find_product(
        ".log",
        exclude=(
            MDRUN_CONSOLE,
        ),
    )

    continuation_gro = find_product(
        ".gro"
    )

    continuation_cpt = find_product(
        ".cpt"
    )

    trajectory_check = run_command(
        [
            gmx,
            "check",
            "-f",
            str(continuation_xtc),
        ],
        cwd=OUTPUT_ROOT,
    )

    TRAJECTORY_CHECK_LOG.write_text(
        trajectory_check.stdout,
        encoding="utf-8",
    )

    check_first_time, check_last_time = (
        parse_check_times(
            trajectory_check.stdout
        )
    )

    (
        initial_title,
        initial_atoms,
        initial_box,
    ) = read_gro(
        SOURCE_INITIAL_GRO
    )

    (
        final_title,
        final_atoms,
        final_box,
    ) = read_gro(
        continuation_gro
    )

    if len(initial_atoms) != EXPECTED_ATOMS:
        raise RuntimeError(
            "Unexpected initial atom count."
        )

    if len(final_atoms) != EXPECTED_ATOMS:
        raise RuntimeError(
            "Unexpected continuation final atom count."
        )

    initial_positions = positions(
        initial_atoms
    )

    final_positions = positions(
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

    caps_final_rms, caps_final_max = (
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

    frames, temporary_frame_root = extract_frames(
        gmx,
        continuation_xtc,
    )

    frame_count = len(frames)

    if frame_count == 61:
        inferred_first_time = 20.0

    elif frame_count == 60:
        inferred_first_time = 20.5

    else:
        inferred_first_time = math.nan

    continuation_first_time = (
        check_first_time
        if check_first_time is not None
        else inferred_first_time
    )

    continuation_last_time = (
        check_last_time
        if check_last_time is not None
        else (
            continuation_first_time
            + (
                frame_count - 1
            )
            * FRAME_INTERVAL_PS
        )
    )

    occupancy_rows = []

    maximum_hbn_trajectory_displacement = 0.0
    maximum_pyr_trajectory_displacement = 0.0
    maximum_caps_trajectory_displacement = 0.0

    minimum_cap_water_distance_seen = math.inf

    try:
        for frame_index, frame_path in enumerate(
            frames
        ):
            _, frame_atoms, frame_box = read_gro(
                frame_path
            )

            if len(frame_atoms) != EXPECTED_ATOMS:
                raise RuntimeError(
                    "Unexpected continuation frame atom count."
                )

            frame_positions = positions(
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

            current_lumen_mask = lumen_mask(
                frame_water_positions,
                frame_box,
                prototype,
            )

            occupancy = int(
                np.count_nonzero(
                    current_lumen_mask
                )
            )

            retained_initial = int(
                np.count_nonzero(
                    current_lumen_mask
                    & initial_lumen_mask
                )
            )

            entered_from_outside = int(
                np.count_nonzero(
                    current_lumen_mask
                    & ~initial_lumen_mask
                )
            )

            cap_distance = (
                minimum_cap_water_distance(
                    frame_water_positions,
                    frame_cap_positions,
                    frame_box,
                )
            )

            _, hbn_frame_max = displacement_metrics(
                initial_positions[
                    :HBN_ATOMS
                ],
                frame_positions[
                    :HBN_ATOMS
                ],
                initial_box,
            )

            _, pyr_frame_max = displacement_metrics(
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

            _, caps_frame_max = displacement_metrics(
                initial_cap_positions,
                frame_cap_positions,
                initial_box,
            )

            maximum_hbn_trajectory_displacement = max(
                maximum_hbn_trajectory_displacement,
                hbn_frame_max,
            )

            maximum_pyr_trajectory_displacement = max(
                maximum_pyr_trajectory_displacement,
                pyr_frame_max,
            )

            maximum_caps_trajectory_displacement = max(
                maximum_caps_trajectory_displacement,
                caps_frame_max,
            )

            minimum_cap_water_distance_seen = min(
                minimum_cap_water_distance_seen,
                cap_distance,
            )

            time_ps = (
                continuation_first_time
                + frame_index
                * FRAME_INTERVAL_PS
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
                        initial_lumen_count
                        - retained_initial
                    ),
                    "waters_entered_from_initial_outside": (
                        entered_from_outside
                    ),
                    "initial_lumen_retention_fraction": (
                        retained_initial
                        / initial_lumen_count
                    ),
                    "minimum_CAP_OW_distance_nm": (
                        cap_distance
                    ),
                    "HBN_max_displacement_nm": (
                        hbn_frame_max
                    ),
                    "PYR_max_displacement_nm": (
                        pyr_frame_max
                    ),
                    "CAPS_max_displacement_nm": (
                        caps_frame_max
                    ),
                }
            )

    finally:
        shutil.rmtree(
            temporary_frame_root,
            ignore_errors=True,
        )

    write_csv(
        CONTINUATION_OCCUPANCY_CSV,
        occupancy_rows,
    )

    source_occupancy_rows = read_csv_rows(
        SOURCE_20PS_OCCUPANCY
    )

    combined_by_time: dict[
        float,
        dict[str, Any],
    ] = {}

    for row in source_occupancy_rows:
        time_ps = round(
            float(
                row["time_ps"]
            ),
            6,
        )

        combined_by_time[
            time_ps
        ] = {
            "time_ps": time_ps,
            "segment": "source_0_to_20ps",
            "lumen_occupancy": int(
                float(
                    row[
                        "lumen_occupancy"
                    ]
                )
            ),
            "initial_lumen_waters_retained": int(
                float(
                    row[
                        "initial_lumen_waters_retained"
                    ]
                )
            ),
            "initial_lumen_waters_lost": int(
                float(
                    row[
                        "initial_lumen_waters_lost"
                    ]
                )
            ),
            "waters_entered_from_initial_outside": int(
                float(
                    row[
                        "waters_entered_from_initial_outside"
                    ]
                )
            ),
            "minimum_CAP_OW_distance_nm": float(
                row[
                    "minimum_CAP_OW_distance_nm"
                ]
            ),
        }

    for row in occupancy_rows:
        time_ps = round(
            float(
                row["time_ps"]
            ),
            6,
        )

        combined_by_time[
            time_ps
        ] = {
            "time_ps": time_ps,
            "segment": "continuation_20_to_50ps",
            "lumen_occupancy": int(
                row[
                    "lumen_occupancy"
                ]
            ),
            "initial_lumen_waters_retained": int(
                row[
                    "initial_lumen_waters_retained"
                ]
            ),
            "initial_lumen_waters_lost": int(
                row[
                    "initial_lumen_waters_lost"
                ]
            ),
            "waters_entered_from_initial_outside": int(
                row[
                    "waters_entered_from_initial_outside"
                ]
            ),
            "minimum_CAP_OW_distance_nm": float(
                row[
                    "minimum_CAP_OW_distance_nm"
                ]
            ),
        }

    combined_rows = [
        combined_by_time[key]
        for key in sorted(
            combined_by_time
        )
    ]

    write_csv(
        COMBINED_OCCUPANCY_CSV,
        combined_rows,
    )

    combined_times = np.asarray(
        [
            float(
                row["time_ps"]
            )
            for row in combined_rows
        ],
        dtype=float,
    )

    combined_occupancy = np.asarray(
        [
            float(
                row[
                    "lumen_occupancy"
                ]
            )
            for row in combined_rows
        ],
        dtype=float,
    )

    combined_retained = np.asarray(
        [
            float(
                row[
                    "initial_lumen_waters_retained"
                ]
            )
            for row in combined_rows
        ],
        dtype=float,
    )

    combined_cap_distance = np.asarray(
        [
            float(
                row[
                    "minimum_CAP_OW_distance_nm"
                ]
            )
            for row in combined_rows
        ],
        dtype=float,
    )

    combined_second_half_mask = (
        combined_times
        >= 25.0 - 1.0e-12
    )

    combined_second_half_slope = (
        linear_slope(
            combined_times[
                combined_second_half_mask
            ],
            combined_occupancy[
                combined_second_half_mask
            ],
        )
    )

    terms = probe_energy_menu(
        gmx,
        continuation_edr,
    )

    temperature_name, temperature_number = (
        resolve_energy_term(
            terms,
            (
                "Temperature",
            ),
            (
                "Temperature",
            ),
        )
    )

    potential_name, potential_number = (
        resolve_energy_term(
            terms,
            (
                "Potential",
            ),
            (
                "Potential",
            ),
        )
    )

    total_energy_name, total_energy_number = (
        resolve_energy_term(
            terms,
            (
                "Total-Energy",
                "Total_Energy",
            ),
            (
                "Total",
                "Energy",
            ),
        )
    )

    cap_sol_name, cap_sol_number = (
        resolve_energy_term(
            terms,
            (
                "LJ-SR:CAPS-SOL",
                "LJ-SR:SOL-CAPS",
            ),
            (
                "LJ-SR",
                "CAPS",
                "SOL",
            ),
        )
    )

    temperature_data = extract_energy_series(
        gmx,
        continuation_edr,
        temperature_number,
        TEMPERATURE_XVG,
    )

    potential_data = extract_energy_series(
        gmx,
        continuation_edr,
        potential_number,
        POTENTIAL_XVG,
    )

    total_energy_data = extract_energy_series(
        gmx,
        continuation_edr,
        total_energy_number,
        TOTAL_ENERGY_XVG,
    )

    cap_sol_data = extract_energy_series(
        gmx,
        continuation_edr,
        cap_sol_number,
        CAP_SOL_LJ_XVG,
    )

    continuation_energy_rows = align_energy_series(
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
        CONTINUATION_ENERGY_CSV,
        continuation_energy_rows,
    )

    energy_times = temperature_data[
        :,
        0
    ]

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

    continuation_temperature_stats = statistics(
        energy_times,
        temperatures,
    )

    last15_mask = (
        energy_times
        >= 35.0 - 1.0e-12
    )

    last15_temperature_stats = statistics(
        energy_times[
            last15_mask
        ],
        temperatures[
            last15_mask
        ],
    )

    continuation_log_text = (
        continuation_log.read_text(
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
            continuation_log_text,
            flags=re.IGNORECASE,
        ):
            instability_hits.append(
                pattern
            )

    finished = bool(
        re.search(
            r"Finished mdrun",
            continuation_log_text,
            flags=re.IGNORECASE,
        )
    )

    checkpoint_reference_present = bool(
        re.search(
            r"checkpoint",
            continuation_log_text,
            flags=re.IGNORECASE,
        )
    )

    minimum_occupancy_required = int(
        math.ceil(
            MIN_OCCUPANCY_FRACTION
            * initial_lumen_count
        )
    )

    continuation_occupancies = np.asarray(
        [
            float(
                row[
                    "lumen_occupancy"
                ]
            )
            for row in occupancy_rows
        ],
        dtype=float,
    )

    continuation_retained = np.asarray(
        [
            float(
                row[
                    "initial_lumen_waters_retained"
                ]
            )
            for row in occupancy_rows
        ],
        dtype=float,
    )

    final_unique_time_spacing = (
        np.diff(
            combined_times
        )
    )

    gates = {
        "mdrun_return_code_zero": (
            mdrun_return_code == 0
        ),
        "mdrun_finished": finished,
        "checkpoint_continuation_was_reported": (
            checkpoint_reference_present
        ),
        "no_instability_signatures": (
            len(instability_hits) == 0
        ),
        "trajectory_check_return_code_zero": (
            trajectory_check.returncode == 0
        ),
        "continuation_frame_count_is_60_or_61": (
            frame_count
            in ALLOWED_CONTINUATION_FRAME_COUNTS
        ),
        "continuation_last_time_is_50ps": (
            math.isclose(
                continuation_last_time,
                CONTINUATION_END_PS,
                rel_tol=0.0,
                abs_tol=1.0e-6,
            )
        ),
        "combined_trajectory_has_101_unique_frames": (
            len(combined_rows)
            == EXPECTED_COMBINED_UNIQUE_FRAMES
        ),
        "combined_trajectory_starts_at_0ps": (
            math.isclose(
                float(
                    combined_times[0]
                ),
                0.0,
                rel_tol=0.0,
                abs_tol=1.0e-10,
            )
        ),
        "combined_trajectory_ends_at_50ps": (
            math.isclose(
                float(
                    combined_times[-1]
                ),
                CONTINUATION_END_PS,
                rel_tol=0.0,
                abs_tol=1.0e-10,
            )
        ),
        "combined_frame_spacing_is_0p5ps": (
            np.allclose(
                final_unique_time_spacing,
                FRAME_INTERVAL_PS,
                rtol=0.0,
                atol=1.0e-8,
            )
        ),
        "final_atom_count_is_68314": (
            len(final_atoms)
            == EXPECTED_ATOMS
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
            caps_final_max
            <= FROZEN_FINAL_TOLERANCE_NM
        ),
        "HBN_continuation_trajectory_is_frozen": (
            maximum_hbn_trajectory_displacement
            <= FROZEN_TRAJECTORY_TOLERANCE_NM
        ),
        "PYR_continuation_trajectory_is_frozen": (
            maximum_pyr_trajectory_displacement
            <= FROZEN_TRAJECTORY_TOLERANCE_NM
        ),
        "CAPS_continuation_trajectory_is_frozen": (
            maximum_caps_trajectory_displacement
            <= FROZEN_TRAJECTORY_TOLERANCE_NM
        ),
        "water_coordinates_are_mobile": (
            water_final_rms > 0.01
        ),
        "continuation_temperature_mean_is_295_to_305K": (
            MIN_TEMPERATURE_MEAN_K
            <= continuation_temperature_stats[
                "mean"
            ]
            <= MAX_TEMPERATURE_MEAN_K
        ),
        "continuation_temperature_std_is_at_most_5K": (
            continuation_temperature_stats[
                "std"
            ]
            <= MAX_TEMPERATURE_STD_K
        ),
        "continuation_temperature_range_is_280_to_320K": (
            continuation_temperature_stats[
                "min"
            ]
            >= MIN_INSTANTANEOUS_TEMPERATURE_K
            and continuation_temperature_stats[
                "max"
            ]
            <= MAX_INSTANTANEOUS_TEMPERATURE_K
        ),
        "continuation_temperature_slope_is_small": (
            abs(
                continuation_temperature_stats[
                    "slope"
                ]
            )
            <= MAX_TEMPERATURE_SLOPE_K_PER_PS
        ),
        "last15ps_temperature_mean_is_295_to_305K": (
            MIN_TEMPERATURE_MEAN_K
            <= last15_temperature_stats[
                "mean"
            ]
            <= MAX_TEMPERATURE_MEAN_K
        ),
        "last15ps_temperature_std_is_at_most_5K": (
            last15_temperature_stats[
                "std"
            ]
            <= MAX_TEMPERATURE_STD_K
        ),
        "last15ps_temperature_slope_is_small": (
            abs(
                last15_temperature_stats[
                    "slope"
                ]
            )
            <= MAX_TEMPERATURE_SLOPE_K_PER_PS
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
        "CAP_SOL_energy_remains_below_500kJmol": (
            float(
                np.max(
                    cap_sol_energies
                )
            )
            <= MAX_CAP_SOL_LJ_KJ_MOL
        ),
        "initial_lumen_occupancy_is_428": (
            initial_lumen_count
            == INITIAL_LUMEN_WATERS
        ),
        "no_complete_lumen_drying_over_0_to_50ps": (
            np.min(
                combined_occupancy
            )
            > 0.0
        ),
        "combined_minimum_occupancy_is_at_least_90_percent": (
            int(
                np.min(
                    combined_occupancy
                )
            )
            >= minimum_occupancy_required
        ),
        "continuation_endpoint_occupancy_is_at_least_90_percent": (
            int(
                continuation_occupancies[-1]
            )
            >= minimum_occupancy_required
        ),
        "continuation_endpoint_initial_retention_is_at_least_90_percent": (
            int(
                continuation_retained[-1]
            )
            >= minimum_occupancy_required
        ),
        "combined_second_half_occupancy_slope_is_small": (
            abs(
                combined_second_half_slope
            )
            <= MAX_SECOND_HALF_OCCUPANCY_SLOPE
        ),
        "CAP_OW_distance_remains_above_0p15nm": (
            float(
                np.min(
                    combined_cap_distance
                )
            )
            >= MIN_CAP_WATER_DISTANCE_NM
        ),
    }

    failed_gates = [
        name
        for name, passed
        in gates.items()
        if not passed
    ]

    accepted = (
        len(failed_gates) == 0
    )

    decision = (
        "R1_FROZEN_SOLUTE_50PS_POSITIVE_CONTROL_VALIDATED"
        if accepted
        else
        "R1_FROZEN_SOLUTE_50PS_REQUIRES_REVIEW"
    )

    required_next_step = (
        "BEGIN_R2_PARTIAL_CAP_DESIGN_AND_STATIC_GATE"
        if accepted
        else
        "REVIEW_R1_50PS_CONTINUATION_GATE_FAILURES"
    )

    summary = {
        "decision": decision,
        "mdrun_return_code": (
            mdrun_return_code
        ),
        "mdrun_finished": finished,
        "checkpoint_reference_present": (
            checkpoint_reference_present
        ),
        "continuation_XTC": relative(
            continuation_xtc
        ),
        "continuation_EDR": relative(
            continuation_edr
        ),
        "continuation_LOG": relative(
            continuation_log
        ),
        "continuation_final_GRO": relative(
            continuation_gro
        ),
        "continuation_CPT": relative(
            continuation_cpt
        ),
        "continuation_frames": frame_count,
        "continuation_first_time_ps": (
            continuation_first_time
        ),
        "continuation_last_time_ps": (
            continuation_last_time
        ),
        "combined_unique_frames": (
            len(combined_rows)
        ),
        "combined_start_time_ps": (
            float(
                combined_times[0]
            )
        ),
        "combined_end_time_ps": (
            float(
                combined_times[-1]
            )
        ),
        "temperature_term": (
            temperature_name
        ),
        "continuation_temperature_mean_K": (
            continuation_temperature_stats[
                "mean"
            ]
        ),
        "continuation_temperature_std_K": (
            continuation_temperature_stats[
                "std"
            ]
        ),
        "continuation_temperature_min_K": (
            continuation_temperature_stats[
                "min"
            ]
        ),
        "continuation_temperature_max_K": (
            continuation_temperature_stats[
                "max"
            ]
        ),
        "continuation_temperature_slope_K_per_ps": (
            continuation_temperature_stats[
                "slope"
            ]
        ),
        "last15_temperature_mean_K": (
            last15_temperature_stats[
                "mean"
            ]
        ),
        "last15_temperature_std_K": (
            last15_temperature_stats[
                "std"
            ]
        ),
        "last15_temperature_min_K": (
            last15_temperature_stats[
                "min"
            ]
        ),
        "last15_temperature_max_K": (
            last15_temperature_stats[
                "max"
            ]
        ),
        "last15_temperature_slope_K_per_ps": (
            last15_temperature_stats[
                "slope"
            ]
        ),
        "potential_term": potential_name,
        "potential_initial_kJ_mol": float(
            potentials[0]
        ),
        "potential_final_kJ_mol": float(
            potentials[-1]
        ),
        "potential_change_kJ_mol": float(
            potentials[-1]
            - potentials[0]
        ),
        "total_energy_term": (
            total_energy_name
        ),
        "total_energy_initial_kJ_mol": float(
            total_energies[0]
        ),
        "total_energy_final_kJ_mol": float(
            total_energies[-1]
        ),
        "CAP_SOL_LJ_term": (
            cap_sol_name
        ),
        "CAP_SOL_LJ_initial_kJ_mol": float(
            cap_sol_energies[0]
        ),
        "CAP_SOL_LJ_final_kJ_mol": float(
            cap_sol_energies[-1]
        ),
        "CAP_SOL_LJ_maximum_kJ_mol": float(
            np.max(
                cap_sol_energies
            )
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
            caps_final_rms
        ),
        "CAPS_final_max_displacement_nm": (
            caps_final_max
        ),
        "HBN_max_continuation_displacement_nm": (
            maximum_hbn_trajectory_displacement
        ),
        "PYR_max_continuation_displacement_nm": (
            maximum_pyr_trajectory_displacement
        ),
        "CAPS_max_continuation_displacement_nm": (
            maximum_caps_trajectory_displacement
        ),
        "waterO_final_RMS_displacement_nm": (
            water_final_rms
        ),
        "waterO_final_max_displacement_nm": (
            water_final_max
        ),
        "initial_lumen_occupancy": (
            initial_lumen_count
        ),
        "combined_mean_lumen_occupancy": float(
            np.mean(
                combined_occupancy
            )
        ),
        "combined_std_lumen_occupancy": float(
            np.std(
                combined_occupancy,
                ddof=1,
            )
        ),
        "combined_minimum_lumen_occupancy": int(
            np.min(
                combined_occupancy
            )
        ),
        "combined_maximum_lumen_occupancy": int(
            np.max(
                combined_occupancy
            )
        ),
        "endpoint_lumen_occupancy": int(
            continuation_occupancies[-1]
        ),
        "endpoint_initial_lumen_waters_retained": int(
            continuation_retained[-1]
        ),
        "endpoint_initial_lumen_retention_fraction": float(
            continuation_retained[-1]
            / initial_lumen_count
        ),
        "combined_zero_occupancy_fraction": float(
            np.mean(
                combined_occupancy == 0.0
            )
        ),
        "combined_second_half_occupancy_slope_waters_per_ps": (
            combined_second_half_slope
        ),
        "minimum_CAP_OW_distance_0_to_50ps_nm": float(
            np.min(
                combined_cap_distance
            )
        ),
        "instability_signature_count": (
            len(instability_hits)
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
        "R2_static_design_authorized": (
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
        [summary],
    )

    gate_rows = [
        {
            "gate": name,
            "pass": passed,
        }
        for name, passed in gates.items()
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
        for name, passed in gates.items()
    )

    REPORT_MD.write_text(
        f"""# R1 Frozen-Solute NVT 50 ps Validation

## Scope

The R1 fully capped neutral steric positive control was continued from
20 to 50 ps using the exact 20 ps checkpoint.

No velocities were regenerated. HBN, all pyrenes, and both caps
remained frozen. Only TIP4P/2005 water was mobile and thermostatted.

## Execution

- Mdrun return code:
  **{mdrun_return_code}**
- Finished mdrun:
  **{'YES' if finished else 'NO'}**
- Checkpoint continuation reported:
  **{'YES' if checkpoint_reference_present else 'NO'}**
- Continuation frames:
  **{frame_count}**
- Continuation time range:
  **{continuation_first_time:.3f}–{continuation_last_time:.3f} ps**
- Combined unique frames:
  **{len(combined_rows)}**
- Combined time range:
  **{combined_times[0]:.3f}–{combined_times[-1]:.3f} ps**
- Instability signatures:
  **{'NONE' if not instability_hits else ' | '.join(instability_hits)}**

## Temperature

Continuation 20–50 ps:

- Mean ± standard deviation:
  **{continuation_temperature_stats['mean']:.4f}
  ± {continuation_temperature_stats['std']:.4f} K**
- Minimum/maximum:
  **{continuation_temperature_stats['min']:.4f}/
  {continuation_temperature_stats['max']:.4f} K**
- Linear slope:
  **{continuation_temperature_stats['slope']:.6f} K/ps**

Final 15 ps:

- Mean ± standard deviation:
  **{last15_temperature_stats['mean']:.4f}
  ± {last15_temperature_stats['std']:.4f} K**
- Minimum/maximum:
  **{last15_temperature_stats['min']:.4f}/
  {last15_temperature_stats['max']:.4f} K**
- Linear slope:
  **{last15_temperature_stats['slope']:.6f} K/ps**

## Frozen-group integrity

Final RMS/max displacement:

- HBN:
  **{hbn_final_rms:.12f}/{hbn_final_max:.12f} nm**
- PYR:
  **{pyr_final_rms:.12f}/{pyr_final_max:.12f} nm**
- CAPS:
  **{caps_final_rms:.12f}/{caps_final_max:.12f} nm**

Maximum displacement in the continuation trajectory:

- HBN:
  **{maximum_hbn_trajectory_displacement:.6f} nm**
- PYR:
  **{maximum_pyr_trajectory_displacement:.6f} nm**
- CAPS:
  **{maximum_caps_trajectory_displacement:.6f} nm**

## Confinement over 0–50 ps

- Initial lumen occupancy:
  **{initial_lumen_count} waters**
- Mean ± standard deviation:
  **{np.mean(combined_occupancy):.4f}
  ± {np.std(combined_occupancy, ddof=1):.4f} waters**
- Minimum/maximum:
  **{int(np.min(combined_occupancy))}/
  {int(np.max(combined_occupancy))} waters**
- Endpoint occupancy:
  **{int(continuation_occupancies[-1])} waters**
- Endpoint initially luminal waters retained:
  **{int(continuation_retained[-1])}/{initial_lumen_count}
  ({continuation_retained[-1] / initial_lumen_count:.6f})**
- Zero-occupancy fraction:
  **{np.mean(combined_occupancy == 0.0):.6f}**
- Occupancy slope over 25–50 ps:
  **{combined_second_half_slope:.6f} waters/ps**
- Minimum CAP–OW distance:
  **{np.min(combined_cap_distance):.6f} nm**

## Energy

- Potential initial/final/change:
  **{potentials[0]:.6f}/
  {potentials[-1]:.6f}/
  {potentials[-1] - potentials[0]:.6f} kJ/mol**
- Total energy initial/final:
  **{total_energies[0]:.6f}/
  {total_energies[-1]:.6f} kJ/mol**
- CAP–SOL LJ initial/final/maximum:
  **{cap_sol_energies[0]:.6f}/
  {cap_sol_energies[-1]:.6f}/
  {np.max(cap_sol_energies):.6f} kJ/mol**

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- R2 partial-cap static design authorized:
  **{'YES' if accepted else 'NO'}**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `{required_next_step}`

A passing result validates R1 only as a frozen neutral steric positive
control for the confinement-analysis methodology. It does not establish
that R1 is chemically realizable or appropriate as the final device.
""",
        encoding="utf-8",
    )

    print()
    print(
        "Day023 R1 checkpoint continuation to "
        "50 ps execution and validation completed."
    )

    print(
        "Mdrun / trajectory-check return codes: "
        f"{mdrun_return_code}/"
        f"{trajectory_check.returncode}"
    )

    print(
        "Finished / checkpoint reference / instabilities: "
        f"{'YES' if finished else 'NO'} / "
        f"{'YES' if checkpoint_reference_present else 'NO'} / "
        + (
            "NONE"
            if not instability_hits
            else " | ".join(
                instability_hits
            )
        )
    )

    print(
        "Continuation frames / time range: "
        f"{frame_count} / "
        f"{continuation_first_time:.3f}-"
        f"{continuation_last_time:.3f} ps"
    )

    print(
        "Combined unique frames / time range: "
        f"{len(combined_rows)} / "
        f"{combined_times[0]:.3f}-"
        f"{combined_times[-1]:.3f} ps"
    )

    print(
        "Continuation temperature mean/std/min/max/slope: "
        f"{continuation_temperature_stats['mean']:.4f}/"
        f"{continuation_temperature_stats['std']:.4f}/"
        f"{continuation_temperature_stats['min']:.4f}/"
        f"{continuation_temperature_stats['max']:.4f}/"
        f"{continuation_temperature_stats['slope']:.6f} "
        "K / K / K / K / K ps^-1"
    )

    print(
        "Last-15ps temperature mean/std/min/max/slope: "
        f"{last15_temperature_stats['mean']:.4f}/"
        f"{last15_temperature_stats['std']:.4f}/"
        f"{last15_temperature_stats['min']:.4f}/"
        f"{last15_temperature_stats['max']:.4f}/"
        f"{last15_temperature_stats['slope']:.6f} "
        "K / K / K / K / K ps^-1"
    )

    print(
        "Potential initial/final/change: "
        f"{potentials[0]:.6f}/"
        f"{potentials[-1]:.6f}/"
        f"{potentials[-1] - potentials[0]:.6f} kJ/mol"
    )

    print(
        "CAP-SOL LJ initial/final/maximum: "
        f"{cap_sol_energies[0]:.6f}/"
        f"{cap_sol_energies[-1]:.6f}/"
        f"{np.max(cap_sol_energies):.6f} kJ/mol"
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
        f"{caps_final_rms:.12f}/"
        f"{caps_final_max:.12f} nm"
    )

    print(
        "Maximum continuation displacement HBN/PYR/CAPS: "
        f"{maximum_hbn_trajectory_displacement:.6f}/"
        f"{maximum_pyr_trajectory_displacement:.6f}/"
        f"{maximum_caps_trajectory_displacement:.6f} nm"
    )

    print(
        "Water-O final RMS/max displacement: "
        f"{water_final_rms:.6f}/"
        f"{water_final_max:.6f} nm"
    )

    print(
        "Combined lumen occupancy "
        "initial/mean/min/max/endpoint: "
        f"{initial_lumen_count}/"
        f"{np.mean(combined_occupancy):.4f}/"
        f"{int(np.min(combined_occupancy))}/"
        f"{int(np.max(combined_occupancy))}/"
        f"{int(continuation_occupancies[-1])}"
    )

    print(
        "Endpoint initial-lumen retention: "
        f"{int(continuation_retained[-1])}/"
        f"{initial_lumen_count} "
        f"({continuation_retained[-1] / initial_lumen_count:.6f})"
    )

    print(
        "Combined second-half occupancy slope: "
        f"{combined_second_half_slope:.6f} waters/ps"
    )

    print(
        "Minimum CAP-OW distance over 0-50 ps: "
        f"{np.min(combined_cap_distance):.6f} nm"
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
        "R2 partial-cap static design authorized: "
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
        f"Wrote: {relative(continuation_gro)}"
    )

    print(
        f"Wrote: {relative(continuation_xtc)}"
    )

    print(
        f"Wrote: {relative(continuation_edr)}"
    )

    print(
        f"Wrote: {relative(continuation_cpt)}"
    )

    print(
        f"Wrote: {relative(CONTINUATION_OCCUPANCY_CSV)}"
    )

    print(
        f"Wrote: {relative(COMBINED_OCCUPANCY_CSV)}"
    )

    print(
        f"Wrote: {relative(CONTINUATION_ENERGY_CSV)}"
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
            "R1 frozen-solute 50 ps validation "
            "requires review."
        )


if __name__ == "__main__":
    main()
