#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Iterator

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

DAY023_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design"
)

PROTOTYPE_ROOT = (
    DAY023_ROOT
    / "02_r1_steric_cap_prototype"
)

PREPARATION_20PS_ROOT = (
    DAY023_ROOT
    / "15_r2_frozen_solute_nvt_20ps_preparation"
)

RUN_20PS_ROOT = (
    DAY023_ROOT
    / "16_r2_frozen_solute_nvt_20ps"
)

PREPARATION_50PS_ROOT = (
    DAY023_ROOT
    / "18_r2_nvt_50ps_checkpoint_continuation_preparation"
)

AUTHORIZATION_ROOT = (
    PREPARATION_50PS_ROOT
    / "difference_classification"
)

OUTPUT_ROOT = (
    DAY023_ROOT
    / "19_r2_frozen_solute_nvt_20_to_50ps"
)

INPUT_GRO = (
    PREPARATION_20PS_ROOT
    / "r2_frozen_solute_nvt_20ps_input.gro"
)

SOURCE_XTC = (
    RUN_20PS_ROOT
    / "r2_frozen_solute_nvt_20ps.xtc"
)

SOURCE_EDR = (
    RUN_20PS_ROOT
    / "r2_frozen_solute_nvt_20ps.edr"
)

SOURCE_CHECKPOINT = (
    PREPARATION_50PS_ROOT
    / "r2_20ps_source_checkpoint.cpt"
)

EXTENDED_TPR = (
    PREPARATION_50PS_ROOT
    / "r2_frozen_solute_nvt_to_50ps.tpr"
)

AUTHORIZED_CONTRACT = (
    AUTHORIZATION_ROOT
    / "r2_50ps_continuation_authorized_run_contract.json"
)

PROTOTYPE_JSON = (
    PROTOTYPE_ROOT
    / "r1_selected_steric_cap_definition.json"
)

DEFFNM = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20_to_50ps"
)

MDRUN_CONSOLE = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20_to_50ps_mdrun_console.log"
)

TRAJECTORY_CHECK_LOG = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20_to_50ps_trajectory_check.log"
)

FINAL_CHECKPOINT_DUMP = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_50ps_checkpoint_dump.txt"
)

FINAL_CHECKPOINT_DUMP_STDERR = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_50ps_checkpoint_dump.stderr.log"
)

ENERGY_MENU_LOG = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20_to_50ps_energy_menu.txt"
)

CONTINUATION_TEMPERATURE_XVG = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20_to_50ps_temperature.xvg"
)

CONTINUATION_POTENTIAL_XVG = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20_to_50ps_potential.xvg"
)

CONTINUATION_TOTAL_ENERGY_XVG = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20_to_50ps_total_energy.xvg"
)

CONTINUATION_CAP_SOL_LJ_XVG = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20_to_50ps_cap_sol_lj.xvg"
)

COMBINED_OCCUPANCY_CSV = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_combined_0_to_50ps_occupancy.csv"
)

OCCUPANCY_WINDOWS_CSV = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_50ps_occupancy_windows.csv"
)

OCCUPANCY_BLOCKS_CSV = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_50ps_occupancy_blocks.csv"
)

CONTINUATION_ENERGY_CSV = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20_to_50ps_energy_series.csv"
)

OUTPUT_MANIFEST_JSON = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20_to_50ps_output_manifest.json"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_50ps_validation_summary.csv"
)

GATES_CSV = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_50ps_validation_gates.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_FROZEN_SOLUTE_NVT_50PS_VALIDATION_DAY023.md"
)

EXPECTED_AUTHORIZATION_DECISION = (
    "R2_CHECKPOINT_CONTINUATION_TO_50PS_AUTHORIZED"
)

EXPECTED_ATOMS = 68332

HBN_ATOMS = 1680
PYR_ATOMS = 104
SOLUTE_ATOMS = HBN_ATOMS + PYR_ATOMS

WATERS = 16565
WATER_SITES = 4
WATER_ATOMS = WATERS * WATER_SITES

TOTAL_CAPS = 288
CAP_START = SOLUTE_ATOMS + WATER_ATOMS
CAP_STOP = CAP_START + TOTAL_CAPS

INITIAL_LUMEN_OCCUPANCY = 428

SOURCE_TIME_PS = 20.0
TARGET_TIME_PS = 50.0

SOURCE_STEP = 40000
TARGET_STEP = 100000

FRAME_INTERVAL_PS = 0.5
EXPECTED_COMBINED_UNIQUE_FRAMES = 101

FROZEN_TOLERANCE_NM = 1.0e-6
MINIMUM_CAP_OW_DISTANCE_NM = 0.15

TEMPERATURE_TARGET_K = 300.0
TEMPERATURE_MEAN_TOLERANCE_K = 5.0
MAXIMUM_TEMPERATURE_SLOPE_K_PS = 0.10

SOL_DEGREES_OF_FREEDOM = 99387.0

THEORETICAL_TEMPERATURE_STD_K = (
    TEMPERATURE_TARGET_K
    * math.sqrt(
        2.0
        / SOL_DEGREES_OF_FREEDOM
    )
)

MAXIMUM_TEMPERATURE_STD_K = (
    3.0
    * THEORETICAL_TEMPERATURE_STD_K
)

MINIMUM_OCCUPANCY_FRACTION = 0.80
MINIMUM_FINAL_WINDOW_MEAN_FRACTION = 0.90
MINIMUM_ENDPOINT_OCCUPANCY_FRACTION = 0.90
MINIMUM_ENDPOINT_IDENTITY_FRACTION = 0.50

OCCUPANCY_SLOPE_LIMIT_WATER_PS = 0.50
MAXIMUM_LAST10_NET_CHANGE_WATERS = 5

MAXIMUM_CAP_SOL_LJ_KJ_MOL = 100.0

WINDOW_STARTS_PS = (
    0.0,
    10.0,
    20.0,
    25.0,
    30.0,
    35.0,
    40.0,
)

BLOCKS_PS = (
    (0.0, 10.0),
    (10.0, 20.0),
    (20.0, 30.0),
    (30.0, 40.0),
    (40.0, 50.0),
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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

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
) -> tuple[float, float]:
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


def minimum_pair_distance(
    first: np.ndarray,
    second: np.ndarray,
    box: np.ndarray,
    *,
    chunk_size: int = 768,
) -> float:
    minimum_distance = math.inf

    for start in range(
        0,
        len(first),
        chunk_size,
    ):
        stop = min(
            start + chunk_size,
            len(first),
        )

        displacement = (
            first[
                start:stop,
                None,
                :,
            ]
            - second[
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


def linear_slope(
    time: np.ndarray,
    values: np.ndarray,
) -> float:
    if len(time) < 2:
        return math.nan

    return float(
        np.polyfit(
            time,
            values,
            1,
        )[0]
    )


def describe_window(
    time: np.ndarray,
    occupancy: np.ndarray,
    retained: np.ndarray,
    start_ps: float,
) -> dict[str, Any]:
    mask = (
        time >= start_ps - 1.0e-8
    )

    window_time = time[mask]
    window_occupancy = occupancy[mask]
    window_retained = retained[mask]

    if len(window_time) < 2:
        raise RuntimeError(
            f"Insufficient records for {start_ps}-50 ps."
        )

    return {
        "window_start_ps": start_ps,
        "window_end_ps": float(
            window_time[-1]
        ),
        "points": len(window_time),
        "occupancy_start": int(
            window_occupancy[0]
        ),
        "occupancy_end": int(
            window_occupancy[-1]
        ),
        "occupancy_change": int(
            window_occupancy[-1]
            - window_occupancy[0]
        ),
        "occupancy_mean": float(
            np.mean(
                window_occupancy
            )
        ),
        "occupancy_minimum": int(
            np.min(
                window_occupancy
            )
        ),
        "occupancy_maximum": int(
            np.max(
                window_occupancy
            )
        ),
        "occupancy_slope_water_ps": (
            linear_slope(
                window_time,
                window_occupancy,
            )
        ),
        "retained_initial_start": int(
            window_retained[0]
        ),
        "retained_initial_end": int(
            window_retained[-1]
        ),
        "retained_initial_change": int(
            window_retained[-1]
            - window_retained[0]
        ),
        "retained_initial_slope_water_ps": (
            linear_slope(
                window_time,
                window_retained,
            )
        ),
    }


def describe_block(
    time: np.ndarray,
    occupancy: np.ndarray,
    retained: np.ndarray,
    start_ps: float,
    end_ps: float,
    *,
    include_end: bool,
) -> dict[str, Any]:
    if include_end:
        mask = (
            (time >= start_ps - 1.0e-8)
            & (time <= end_ps + 1.0e-8)
        )
    else:
        mask = (
            (time >= start_ps - 1.0e-8)
            & (time < end_ps - 1.0e-8)
        )

    block_time = time[mask]
    block_occupancy = occupancy[mask]
    block_retained = retained[mask]

    if len(block_time) < 2:
        raise RuntimeError(
            f"Insufficient records for block "
            f"{start_ps}-{end_ps} ps."
        )

    return {
        "block_start_ps": start_ps,
        "block_end_ps": end_ps,
        "points": len(block_time),
        "occupancy_start": int(
            block_occupancy[0]
        ),
        "occupancy_end": int(
            block_occupancy[-1]
        ),
        "occupancy_change": int(
            block_occupancy[-1]
            - block_occupancy[0]
        ),
        "occupancy_mean": float(
            np.mean(
                block_occupancy
            )
        ),
        "occupancy_minimum": int(
            np.min(
                block_occupancy
            )
        ),
        "occupancy_maximum": int(
            np.max(
                block_occupancy
            )
        ),
        "occupancy_slope_water_ps": (
            linear_slope(
                block_time,
                block_occupancy,
            )
        ),
        "retained_initial_start": int(
            block_retained[0]
        ),
        "retained_initial_end": int(
            block_retained[-1]
        ),
    }


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


def resolve_energy_term(
    terms: dict[str, int],
    exact: tuple[str, ...],
    required_tokens: tuple[str, ...],
) -> tuple[str, int]:
    for candidate in exact:
        if candidate in terms:
            return (
                candidate,
                terms[candidate],
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

        if len(fields) >= 2:
            rows.append(
                (
                    float(fields[0]),
                    float(fields[1]),
                )
            )

    data = np.asarray(
        rows,
        dtype=float,
    )

    if (
        data.ndim != 2
        or data.shape[1] != 2
        or len(data) == 0
        or not np.all(
            np.isfinite(data)
        )
    ):
        raise RuntimeError(
            f"Invalid energy series: {output_path}"
        )

    return data


def parse_checkpoint_step(
    text: str,
) -> int:
    for pattern in (
        r"^\s*step\s*=\s*(\d+)",
        r"^\s*step\s+(\d+)",
    ):
        match = re.search(
            pattern,
            text,
            flags=re.MULTILINE,
        )

        if match is not None:
            return int(
                match.group(1)
            )

    raise RuntimeError(
        "Could not parse checkpoint step."
    )


def parse_checkpoint_time(
    text: str,
) -> float:
    for pattern in (
        r"^\s*t\s*=\s*([-+0-9.eE]+)",
        r"^\s*time\s*=\s*([-+0-9.eE]+)",
    ):
        match = re.search(
            pattern,
            text,
            flags=re.MULTILINE,
        )

        if match is not None:
            return float(
                match.group(1)
            )

    raise RuntimeError(
        "Could not parse checkpoint time."
    )


def find_output_product(
    extension: str,
) -> Path:
    candidates = []

    for path in OUTPUT_ROOT.glob(
        f"{DEFFNM.name}*.{extension}"
    ):
        if (
            path == MDRUN_CONSOLE
            or "_prev." in path.name
            or not path.is_file()
            or path.stat().st_size == 0
        ):
            continue

        candidates.append(path)

    if not candidates:
        raise RuntimeError(
            f"No nonempty continuation .{extension} "
            "product was found."
        )

    candidates.sort(
        key=lambda path: (
            path.stat().st_mtime_ns,
            path.stat().st_size,
        ),
        reverse=True,
    )

    return candidates[0]


def trajectory_frames(
    topology: Path,
    trajectory: Path,
) -> Iterator[
    tuple[
        float,
        np.ndarray,
        np.ndarray,
        np.ndarray,
    ]
]:
    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise RuntimeError(
            "MDAnalysis is required for trajectory analysis."
        ) from exc

    universe = mda.Universe(
        str(topology),
        str(trajectory),
    )

    if universe.atoms.n_atoms != EXPECTED_ATOMS:
        raise RuntimeError(
            "Unexpected trajectory atom count: "
            f"{universe.atoms.n_atoms}/{EXPECTED_ATOMS}"
        )

    names = np.asarray(
        universe.atoms.names,
        dtype=object,
    )

    oxygen_indices = np.arange(
        SOLUTE_ATOMS,
        CAP_START,
        WATER_SITES,
        dtype=int,
    )

    if (
        len(oxygen_indices) != WATERS
        or not np.all(
            names[oxygen_indices] == "OW"
        )
    ):
        raise RuntimeError(
            "Water-oxygen atom accounting failed."
        )

    for timestep in universe.trajectory:
        positions = (
            universe.atoms.positions.astype(
                float,
                copy=True,
            )
            / 10.0
        )

        box = (
            np.asarray(
                timestep.dimensions[:3],
                dtype=float,
            )
            / 10.0
        )

        if (
            len(box) != 3
            or not np.all(
                np.isfinite(box)
            )
            or np.any(box <= 0.0)
        ):
            raise RuntimeError(
                "Invalid trajectory box."
            )

        yield (
            float(timestep.time),
            positions,
            box,
            oxygen_indices,
        )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        INPUT_GRO,
        SOURCE_XTC,
        SOURCE_EDR,
        SOURCE_CHECKPOINT,
        EXTENDED_TPR,
        AUTHORIZED_CONTRACT,
        PROTOTYPE_JSON,
    ):
        require_file(required)

    contract = json.loads(
        AUTHORIZED_CONTRACT.read_text(
            encoding="utf-8"
        )
    )

    prototype = json.loads(
        PROTOTYPE_JSON.read_text(
            encoding="utf-8"
        )
    )

    if (
        contract.get("decision")
        != EXPECTED_AUTHORIZATION_DECISION
    ):
        raise RuntimeError(
            "The continuation contract is not authorized."
        )

    if not bool(
        contract.get(
            "execution_authorized",
            False,
        )
    ):
        raise RuntimeError(
            "Execution is not authorized by the contract."
        )

    if bool(
        contract.get(
            "velocity_regeneration",
            True,
        )
    ):
        raise RuntimeError(
            "The contract unexpectedly authorizes "
            "velocity regeneration."
        )

    if bool(
        contract.get(
            "thermostat_state_regeneration",
            True,
        )
    ):
        raise RuntimeError(
            "The contract unexpectedly authorizes "
            "thermostat-state regeneration."
        )

    if bool(
        contract.get(
            "source_0_to_20ps_rerun",
            True,
        )
    ):
        raise RuntimeError(
            "The contract unexpectedly authorizes "
            "rerunning 0-20 ps."
        )

    expected_checkpoint_hash = str(
        contract.get(
            "copied_checkpoint_sha256",
            "",
        )
    )

    observed_checkpoint_hash = sha256(
        SOURCE_CHECKPOINT
    )

    if (
        not expected_checkpoint_hash
        or observed_checkpoint_hash
        != expected_checkpoint_hash
    ):
        raise RuntimeError(
            "The 20 ps checkpoint hash does not match "
            "the authorized contract."
        )

    if int(
        contract.get(
            "source_step",
            -1,
        )
    ) != SOURCE_STEP:
        raise RuntimeError(
            "Unexpected source step in contract."
        )

    if int(
        contract.get(
            "target_step",
            -1,
        )
    ) != TARGET_STEP:
        raise RuntimeError(
            "Unexpected target step in contract."
        )

    reuse_completed_run = (
        os.environ.get(
            "R2_CONTINUATION_REUSE_COMPLETED_RUN",
            "0",
        )
        == "1"
    )

    existing_dynamic_products = [
        path
        for path in OUTPUT_ROOT.glob(
            f"{DEFFNM.name}*"
        )
        if (
            path.is_file()
            and path.stat().st_size > 0
            and path != MDRUN_CONSOLE
            and not path.name.endswith(
                "_output_manifest.json"
            )
        )
    ]

    gmx = locate_gmx()

    if reuse_completed_run:
        require_file(MDRUN_CONSOLE)

        if not existing_dynamic_products:
            raise RuntimeError(
                "Reuse was requested, but no continuation "
                "products were found."
            )

        prior_text = MDRUN_CONSOLE.read_text(
            encoding="utf-8",
            errors="replace",
        )

        if not (
            "Writing final coordinates"
            in prior_text
            or re.search(
                r"\bFinished\s+mdrun\b",
                prior_text,
                flags=re.IGNORECASE,
            )
            is not None
        ):
            raise RuntimeError(
                "Existing console output does not confirm "
                "a completed continuation."
            )

        print(
            "Reusing the completed R2 20-to-50 ps "
            "checkpoint continuation."
        )

        print(
            "No molecular-dynamics steps will be repeated."
        )

        mdrun_return_code = 0

    else:
        if existing_dynamic_products:
            raise RuntimeError(
                "Continuation products already exist. "
                "Use R2_CONTINUATION_REUSE_COMPLETED_RUN=1 "
                "only after confirming that mdrun completed."
            )

        print(
            "Starting exact R2 checkpoint continuation "
            "from 20 to 50 ps."
        )

        print(
            "Source checkpoint: step 40000 / 20.0 ps."
        )

        print(
            "Target: step 100000 / 50.0 ps."
        )

        print(
            "Only the remaining 60000 steps will be run."
        )

        print(
            "Velocity and thermostat-state regeneration: DISABLED."
        )

        print(
            "Source 0-to-20 ps rerun: DISABLED."
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
                "12",
            ],
            cwd=OUTPUT_ROOT,
            log_path=MDRUN_CONSOLE,
        )

    continuation_xtc = find_output_product(
        "xtc"
    )

    continuation_edr = find_output_product(
        "edr"
    )

    continuation_log = find_output_product(
        "log"
    )

    continuation_gro = find_output_product(
        "gro"
    )

    continuation_cpt = find_output_product(
        "cpt"
    )

    combined_log = (
        MDRUN_CONSOLE.read_text(
            encoding="utf-8",
            errors="replace",
        )
        + "\n"
        + continuation_log.read_text(
            encoding="utf-8",
            errors="replace",
        )
    )

    completion_confirmed = (
        "Writing final coordinates"
        in combined_log
        or re.search(
            r"\bFinished\s+mdrun\b",
            combined_log,
            flags=re.IGNORECASE,
        )
        is not None
    )

    checkpoint_continuation_confirmed = any(
        pattern.search(combined_log)
        is not None
        for pattern in (
            re.compile(
                r"reading\s+checkpoint",
                flags=re.IGNORECASE,
            ),
            re.compile(
                r"continuing\s+from\s+checkpoint",
                flags=re.IGNORECASE,
            ),
            re.compile(
                r"starting\s+from\s+checkpoint",
                flags=re.IGNORECASE,
            ),
        )
    )

    instability_patterns = {
        "standalone_nan": re.compile(
            r"\bnan\b",
            flags=re.IGNORECASE,
        ),
        "fatal_error": re.compile(
            r"\bfatal\s+error\b",
            flags=re.IGNORECASE,
        ),
        "segmentation_fault": re.compile(
            r"\bsegmentation\s+fault\b",
            flags=re.IGNORECASE,
        ),
        "lincs_warning": re.compile(
            r"\blincs\s+warning\b",
            flags=re.IGNORECASE,
        ),
        "constraint_warning": re.compile(
            r"\bconstraint\s+warning\b",
            flags=re.IGNORECASE,
        ),
        "settle_failure": re.compile(
            r"\bwater\s+molecule\s+ca(?:n\s*not|nnot)\s+"
            r"be\s+settled\b",
            flags=re.IGNORECASE,
        ),
    }

    instability_hits = [
        name
        for name, pattern
        in instability_patterns.items()
        if pattern.search(
            combined_log
        )
    ]

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

    checkpoint_dump = subprocess.run(
        [
            gmx,
            "dump",
            "-cp",
            str(continuation_cpt),
        ],
        cwd=OUTPUT_ROOT,
        env={
            **os.environ,
            "GMX_MAXBACKUP": "-1",
        },
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )

    FINAL_CHECKPOINT_DUMP.write_text(
        checkpoint_dump.stdout,
        encoding="utf-8",
    )

    FINAL_CHECKPOINT_DUMP_STDERR.write_text(
        checkpoint_dump.stderr,
        encoding="utf-8",
    )

    if checkpoint_dump.returncode != 0:
        raise RuntimeError(
            "Could not dump the final continuation checkpoint."
        )

    final_checkpoint_step = parse_checkpoint_step(
        checkpoint_dump.stdout
    )

    final_checkpoint_time = parse_checkpoint_time(
        checkpoint_dump.stdout
    )

    source_raw_times = []
    continuation_raw_times = []

    reference_positions: np.ndarray | None = None
    reference_box: np.ndarray | None = None
    initial_lumen_mask: np.ndarray | None = None
    initial_lumen_count: int | None = None

    maximum_hbn_displacement = 0.0
    maximum_pyr_displacement = 0.0
    maximum_caps_displacement = 0.0

    minimum_cap_ow_distance_seen = math.inf

    final_water_rms = math.nan
    final_water_max = math.nan

    combined_rows: list[dict[str, Any]] = []
    seen_times: set[float] = set()

    trajectory_sources = (
        (
            "source_0_to_20ps",
            SOURCE_XTC,
            source_raw_times,
        ),
        (
            "continuation_20_to_50ps",
            continuation_xtc,
            continuation_raw_times,
        ),
    )

    for (
        source_label,
        trajectory_path,
        raw_times,
    ) in trajectory_sources:
        for (
            time_ps,
            current_positions,
            current_box,
            oxygen_indices,
        ) in trajectory_frames(
            INPUT_GRO,
            trajectory_path,
        ):
            raw_times.append(
                time_ps
            )

            time_key = round(
                time_ps,
                6,
            )

            if time_key in seen_times:
                continue

            seen_times.add(
                time_key
            )

            water_positions = (
                current_positions[
                    oxygen_indices
                ]
            )

            cap_positions = (
                current_positions[
                    CAP_START:
                    CAP_STOP
                ]
            )

            current_lumen_mask = lumen_mask(
                water_positions,
                current_box,
                prototype,
            )

            current_occupancy = int(
                np.count_nonzero(
                    current_lumen_mask
                )
            )

            current_cap_ow_distance = minimum_pair_distance(
                water_positions,
                cap_positions,
                current_box,
            )

            minimum_cap_ow_distance_seen = min(
                minimum_cap_ow_distance_seen,
                current_cap_ow_distance,
            )

            if reference_positions is None:
                reference_positions = (
                    current_positions.copy()
                )

                reference_box = (
                    current_box.copy()
                )

                initial_lumen_mask = (
                    current_lumen_mask.copy()
                )

                initial_lumen_count = int(
                    np.count_nonzero(
                        initial_lumen_mask
                    )
                )

            assert reference_positions is not None
            assert reference_box is not None
            assert initial_lumen_mask is not None
            assert initial_lumen_count is not None

            _, hbn_max = displacement_metrics(
                reference_positions[
                    :HBN_ATOMS
                ],
                current_positions[
                    :HBN_ATOMS
                ],
                current_box,
            )

            _, pyr_max = displacement_metrics(
                reference_positions[
                    HBN_ATOMS:
                    SOLUTE_ATOMS
                ],
                current_positions[
                    HBN_ATOMS:
                    SOLUTE_ATOMS
                ],
                current_box,
            )

            _, caps_max = displacement_metrics(
                reference_positions[
                    CAP_START:
                    CAP_STOP
                ],
                current_positions[
                    CAP_START:
                    CAP_STOP
                ],
                current_box,
            )

            maximum_hbn_displacement = max(
                maximum_hbn_displacement,
                hbn_max,
            )

            maximum_pyr_displacement = max(
                maximum_pyr_displacement,
                pyr_max,
            )

            maximum_caps_displacement = max(
                maximum_caps_displacement,
                caps_max,
            )

            retained_initial = int(
                np.count_nonzero(
                    initial_lumen_mask
                    & current_lumen_mask
                )
            )

            combined_rows.append(
                {
                    "source_segment": source_label,
                    "time_ps": time_ps,
                    "lumen_occupancy": (
                        current_occupancy
                    ),
                    "initial_lumen_waters_retained": (
                        retained_initial
                    ),
                    "initial_lumen_identity_retention_fraction": (
                        retained_initial
                        / initial_lumen_count
                    ),
                    "noninitial_lumen_waters": (
                        current_occupancy
                        - retained_initial
                    ),
                    "minimum_CAP_OW_distance_nm": (
                        current_cap_ow_distance
                    ),
                }
            )

            if math.isclose(
                time_ps,
                TARGET_TIME_PS,
                rel_tol=0.0,
                abs_tol=1.0e-3,
            ):
                (
                    final_water_rms,
                    final_water_max,
                ) = displacement_metrics(
                    reference_positions[
                        oxygen_indices
                    ],
                    water_positions,
                    current_box,
                )

    combined_rows.sort(
        key=lambda row: float(
            row["time_ps"]
        )
    )

    write_csv(
        COMBINED_OCCUPANCY_CSV,
        combined_rows,
    )

    source_raw_times_array = np.asarray(
        source_raw_times,
        dtype=float,
    )

    continuation_raw_times_array = np.asarray(
        continuation_raw_times,
        dtype=float,
    )

    combined_time = np.asarray(
        [
            float(row["time_ps"])
            for row in combined_rows
        ],
        dtype=float,
    )

    combined_occupancy = np.asarray(
        [
            int(row["lumen_occupancy"])
            for row in combined_rows
        ],
        dtype=int,
    )

    combined_retained = np.asarray(
        [
            int(
                row[
                    "initial_lumen_waters_retained"
                ]
            )
            for row in combined_rows
        ],
        dtype=int,
    )

    combined_intervals = np.diff(
        combined_time
    )

    window_rows = [
        describe_window(
            combined_time,
            combined_occupancy,
            combined_retained,
            start_ps,
        )
        for start_ps in WINDOW_STARTS_PS
    ]

    block_rows = []

    for block_index, (
        start_ps,
        end_ps,
    ) in enumerate(BLOCKS_PS):
        block_rows.append(
            describe_block(
                combined_time,
                combined_occupancy,
                combined_retained,
                start_ps,
                end_ps,
                include_end=(
                    block_index
                    == len(BLOCKS_PS) - 1
                ),
            )
        )

    write_csv(
        OCCUPANCY_WINDOWS_CSV,
        window_rows,
    )

    write_csv(
        OCCUPANCY_BLOCKS_CSV,
        block_rows,
    )

    windows_by_start = {
        float(
            row["window_start_ps"]
        ): row
        for row in window_rows
    }

    final20 = windows_by_start[
        30.0
    ]

    final15 = windows_by_start[
        35.0
    ]

    final10 = windows_by_start[
        40.0
    ]

    menu_probe = run_command(
        [
            gmx,
            "energy",
            "-f",
            str(continuation_edr),
            "-o",
            str(
                OUTPUT_ROOT
                / "r2_frozen_solute_nvt_20_to_50ps_energy_menu_probe.xvg"
            ),
        ],
        cwd=OUTPUT_ROOT,
        input_text="0\n",
    )

    ENERGY_MENU_LOG.write_text(
        menu_probe.stdout,
        encoding="utf-8",
    )

    terms = parse_energy_menu(
        menu_probe.stdout
    )

    energy_requests = {
        "Temperature": (
            ("Temperature",),
            ("Temperature",),
            CONTINUATION_TEMPERATURE_XVG,
        ),
        "Potential": (
            ("Potential",),
            ("Potential",),
            CONTINUATION_POTENTIAL_XVG,
        ),
        "Total-Energy": (
            ("Total-Energy",),
            ("Total-Energy",),
            CONTINUATION_TOTAL_ENERGY_XVG,
        ),
        "CAP_SOL_LJ": (
            (
                "LJ-SR:CAPS-SOL",
                "LJ-SR:SOL-CAPS",
            ),
            (
                "LJ-SR",
                "CAPS",
                "SOL",
            ),
            CONTINUATION_CAP_SOL_LJ_XVG,
        ),
    }

    extracted: dict[str, np.ndarray] = {}
    resolved_names: dict[str, str] = {}

    for label, (
        exact,
        required_tokens,
        output_path,
    ) in energy_requests.items():
        term_name, term_number = resolve_energy_term(
            terms,
            exact,
            required_tokens,
        )

        resolved_names[label] = term_name

        extracted[label] = extract_energy_series(
            gmx,
            continuation_edr,
            term_number,
            output_path,
        )

    energy_time = extracted[
        "Temperature"
    ][
        :,
        0,
    ]

    for label, data in extracted.items():
        if not np.allclose(
            data[:, 0],
            energy_time,
            rtol=0.0,
            atol=1.0e-8,
        ):
            raise RuntimeError(
                f"Energy time grid mismatch for {label}."
            )

    temperature = extracted[
        "Temperature"
    ][
        :,
        1,
    ]

    potential = extracted[
        "Potential"
    ][
        :,
        1,
    ]

    total_energy = extracted[
        "Total-Energy"
    ][
        :,
        1,
    ]

    cap_sol_lj = extracted[
        "CAP_SOL_LJ"
    ][
        :,
        1,
    ]

    energy_rows = []

    for index in range(
        len(energy_time)
    ):
        energy_rows.append(
            {
                "time_ps": float(
                    energy_time[index]
                ),
                "temperature_K": float(
                    temperature[index]
                ),
                "potential_kJ_mol": float(
                    potential[index]
                ),
                "total_energy_kJ_mol": float(
                    total_energy[index]
                ),
                "CAP_SOL_LJ_kJ_mol": float(
                    cap_sol_lj[index]
                ),
            }
        )

    write_csv(
        CONTINUATION_ENERGY_CSV,
        energy_rows,
    )

    final15_energy_mask = (
        energy_time
        >= 35.0 - 1.0e-8
    )

    continuation_temperature_mean = float(
        np.mean(
            temperature
        )
    )

    continuation_temperature_std = float(
        np.std(
            temperature,
            ddof=1,
        )
    )

    continuation_temperature_min = float(
        np.min(
            temperature
        )
    )

    continuation_temperature_max = float(
        np.max(
            temperature
        )
    )

    continuation_temperature_slope = linear_slope(
        energy_time,
        temperature,
    )

    final15_temperature = (
        temperature[
            final15_energy_mask
        ]
    )

    final15_energy_time = (
        energy_time[
            final15_energy_mask
        ]
    )

    final15_temperature_mean = float(
        np.mean(
            final15_temperature
        )
    )

    final15_temperature_std = float(
        np.std(
            final15_temperature,
            ddof=1,
        )
    )

    final15_temperature_slope = linear_slope(
        final15_energy_time,
        final15_temperature,
    )

    occupancy_initial = int(
        combined_occupancy[0]
    )

    occupancy_mean = float(
        np.mean(
            combined_occupancy
        )
    )

    occupancy_minimum = int(
        np.min(
            combined_occupancy
        )
    )

    occupancy_maximum = int(
        np.max(
            combined_occupancy
        )
    )

    occupancy_endpoint = int(
        combined_occupancy[-1]
    )

    endpoint_retained_initial = int(
        combined_retained[-1]
    )

    endpoint_identity_fraction = (
        endpoint_retained_initial
        / occupancy_initial
    )

    minimum_allowed_occupancy = (
        MINIMUM_OCCUPANCY_FRACTION
        * INITIAL_LUMEN_OCCUPANCY
    )

    minimum_final_window_mean = (
        MINIMUM_FINAL_WINDOW_MEAN_FRACTION
        * INITIAL_LUMEN_OCCUPANCY
    )

    minimum_endpoint_occupancy = (
        MINIMUM_ENDPOINT_OCCUPANCY_FRACTION
        * INITIAL_LUMEN_OCCUPANCY
    )

    output_manifest = {
        "continuation_XTC": (
            relative(
                continuation_xtc
            )
        ),
        "continuation_EDR": (
            relative(
                continuation_edr
            )
        ),
        "continuation_LOG": (
            relative(
                continuation_log
            )
        ),
        "continuation_final_GRO": (
            relative(
                continuation_gro
            )
        ),
        "continuation_final_checkpoint": (
            relative(
                continuation_cpt
            )
        ),
        "continuation_XTC_sha256": (
            sha256(
                continuation_xtc
            )
        ),
        "continuation_EDR_sha256": (
            sha256(
                continuation_edr
            )
        ),
        "continuation_checkpoint_sha256": (
            sha256(
                continuation_cpt
            )
        ),
    }

    OUTPUT_MANIFEST_JSON.write_text(
        json.dumps(
            output_manifest,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    continuation_frame_count_acceptable = (
        len(
            continuation_raw_times_array
        )
        in {
            60,
            61,
        }
    )

    continuation_start_acceptable = (
        len(
            continuation_raw_times_array
        )
        > 0
        and (
            math.isclose(
                float(
                    continuation_raw_times_array[
                        0
                    ]
                ),
                20.0,
                rel_tol=0.0,
                abs_tol=1.0e-3,
            )
            or math.isclose(
                float(
                    continuation_raw_times_array[
                        0
                    ]
                ),
                20.5,
                rel_tol=0.0,
                abs_tol=1.0e-3,
            )
        )
    )

    gates = {
        "authorized_contract_is_valid": (
            contract.get("decision")
            == EXPECTED_AUTHORIZATION_DECISION
            and bool(
                contract.get(
                    "execution_authorized",
                    False,
                )
            )
        ),
        "source_checkpoint_hash_matches_contract": (
            observed_checkpoint_hash
            == expected_checkpoint_hash
        ),
        "velocity_regeneration_remained_disabled": (
            not bool(
                contract.get(
                    "velocity_regeneration",
                    True,
                )
            )
        ),
        "thermostat_state_regeneration_remained_disabled": (
            not bool(
                contract.get(
                    "thermostat_state_regeneration",
                    True,
                )
            )
        ),
        "source_0_to_20ps_was_not_rerun": (
            not bool(
                contract.get(
                    "source_0_to_20ps_rerun",
                    True,
                )
            )
        ),
        "mdrun_return_code_zero": (
            mdrun_return_code == 0
        ),
        "mdrun_completion_confirmed": (
            completion_confirmed
        ),
        "checkpoint_continuation_confirmed_in_logs": (
            checkpoint_continuation_confirmed
        ),
        "no_instability_signatures": (
            len(
                instability_hits
            )
            == 0
        ),
        "continuation_trajectory_check_return_code_zero": (
            trajectory_check.returncode
            == 0
        ),
        "final_checkpoint_dump_return_code_zero": (
            checkpoint_dump.returncode
            == 0
        ),
        "final_checkpoint_step_is_100000": (
            final_checkpoint_step
            == TARGET_STEP
        ),
        "final_checkpoint_time_is_50ps": (
            math.isclose(
                final_checkpoint_time,
                TARGET_TIME_PS,
                rel_tol=0.0,
                abs_tol=1.0e-6,
            )
        ),
        "source_trajectory_has_41_frames": (
            len(
                source_raw_times_array
            )
            == 41
        ),
        "continuation_trajectory_has_60_or_61_frames": (
            continuation_frame_count_acceptable
        ),
        "continuation_starts_at_20_or_20p5ps": (
            continuation_start_acceptable
        ),
        "continuation_ends_at_50ps": (
            len(
                continuation_raw_times_array
            )
            > 0
            and math.isclose(
                float(
                    continuation_raw_times_array[
                        -1
                    ]
                ),
                TARGET_TIME_PS,
                rel_tol=0.0,
                abs_tol=1.0e-3,
            )
        ),
        "combined_trajectory_has_101_unique_frames": (
            len(
                combined_time
            )
            == EXPECTED_COMBINED_UNIQUE_FRAMES
        ),
        "combined_trajectory_starts_at_0ps": (
            math.isclose(
                float(
                    combined_time[0]
                ),
                0.0,
                rel_tol=0.0,
                abs_tol=1.0e-3,
            )
        ),
        "combined_trajectory_ends_at_50ps": (
            math.isclose(
                float(
                    combined_time[-1]
                ),
                TARGET_TIME_PS,
                rel_tol=0.0,
                abs_tol=1.0e-3,
            )
        ),
        "combined_trajectory_interval_is_0p5ps": (
            len(
                combined_intervals
            )
            > 0
            and np.allclose(
                combined_intervals,
                FRAME_INTERVAL_PS,
                rtol=0.0,
                atol=1.0e-3,
            )
        ),
        "HBN_remained_frozen": (
            maximum_hbn_displacement
            <= FROZEN_TOLERANCE_NM
        ),
        "PYR_remained_frozen": (
            maximum_pyr_displacement
            <= FROZEN_TOLERANCE_NM
        ),
        "CAPS_remained_frozen": (
            maximum_caps_displacement
            <= FROZEN_TOLERANCE_NM
        ),
        "water_displacement_is_finite": (
            math.isfinite(
                final_water_rms
            )
            and math.isfinite(
                final_water_max
            )
        ),
        "minimum_CAP_OW_distance_is_safe": (
            minimum_cap_ow_distance_seen
            >= MINIMUM_CAP_OW_DISTANCE_NM
        ),
        "continuation_temperature_mean_is_295_to_305K": (
            abs(
                continuation_temperature_mean
                - TEMPERATURE_TARGET_K
            )
            <= TEMPERATURE_MEAN_TOLERANCE_K
        ),
        "continuation_temperature_std_is_acceptable": (
            continuation_temperature_std
            <= MAXIMUM_TEMPERATURE_STD_K
        ),
        "continuation_temperature_slope_is_acceptable": (
            abs(
                continuation_temperature_slope
            )
            <= MAXIMUM_TEMPERATURE_SLOPE_K_PS
        ),
        "final15_temperature_mean_is_295_to_305K": (
            abs(
                final15_temperature_mean
                - TEMPERATURE_TARGET_K
            )
            <= TEMPERATURE_MEAN_TOLERANCE_K
        ),
        "final15_temperature_std_is_acceptable": (
            final15_temperature_std
            <= MAXIMUM_TEMPERATURE_STD_K
        ),
        "final15_temperature_slope_is_acceptable": (
            abs(
                final15_temperature_slope
            )
            <= MAXIMUM_TEMPERATURE_SLOPE_K_PS
        ),
        "initial_lumen_occupancy_is_428": (
            occupancy_initial
            == INITIAL_LUMEN_OCCUPANCY
        ),
        "combined_minimum_occupancy_remains_above_80_percent": (
            occupancy_minimum
            >= minimum_allowed_occupancy
        ),
        "final20_mean_occupancy_is_at_least_90_percent": (
            float(
                final20[
                    "occupancy_mean"
                ]
            )
            >= minimum_final_window_mean
        ),
        "endpoint_occupancy_is_at_least_90_percent": (
            occupancy_endpoint
            >= minimum_endpoint_occupancy
        ),
        "endpoint_initial_identity_retention_is_at_least_50_percent": (
            endpoint_identity_fraction
            >= MINIMUM_ENDPOINT_IDENTITY_FRACTION
        ),
        "final20_occupancy_slope_is_stationary": (
            abs(
                float(
                    final20[
                        "occupancy_slope_water_ps"
                    ]
                )
            )
            <= OCCUPANCY_SLOPE_LIMIT_WATER_PS
        ),
        "final15_occupancy_slope_is_stationary": (
            abs(
                float(
                    final15[
                        "occupancy_slope_water_ps"
                    ]
                )
            )
            <= OCCUPANCY_SLOPE_LIMIT_WATER_PS
        ),
        "final10_occupancy_slope_is_stationary": (
            abs(
                float(
                    final10[
                        "occupancy_slope_water_ps"
                    ]
                )
            )
            <= OCCUPANCY_SLOPE_LIMIT_WATER_PS
        ),
        "final10_net_occupancy_change_is_at_most_5_waters": (
            abs(
                int(
                    final10[
                        "occupancy_change"
                    ]
                )
            )
            <= MAXIMUM_LAST10_NET_CHANGE_WATERS
        ),
        "temperature_series_is_finite": (
            np.all(
                np.isfinite(
                    temperature
                )
            )
        ),
        "potential_series_is_finite": (
            np.all(
                np.isfinite(
                    potential
                )
            )
        ),
        "total_energy_series_is_finite": (
            np.all(
                np.isfinite(
                    total_energy
                )
            )
        ),
        "CAP_SOL_LJ_series_is_finite": (
            np.all(
                np.isfinite(
                    cap_sol_lj
                )
            )
        ),
        "CAP_SOL_LJ_remains_below_100kJmol": (
            float(
                np.max(
                    cap_sol_lj
                )
            )
            <= MAXIMUM_CAP_SOL_LJ_KJ_MOL
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
        "R2_FROZEN_SOLUTE_NVT_50PS_VALIDATED"
        if accepted
        else
        "R2_FROZEN_SOLUTE_NVT_50PS_REQUIRES_REVIEW"
    )

    required_next_step = (
        "COMPARE_R2_WITH_R1_AND_DECIDE_NEXT_ARCHITECTURE_GATE"
        if accepted
        else
        "REVIEW_R2_50PS_STATIONARITY_AND_CONFINEMENT_FAILURES"
    )

    summary = {
        "decision": decision,
        "mdrun_return_code": (
            mdrun_return_code
        ),
        "completion_confirmed": (
            completion_confirmed
        ),
        "checkpoint_continuation_confirmed": (
            checkpoint_continuation_confirmed
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
        "continuation_trajectory_check_return_code": (
            trajectory_check.returncode
        ),
        "final_checkpoint_step": (
            final_checkpoint_step
        ),
        "final_checkpoint_time_ps": (
            final_checkpoint_time
        ),
        "source_raw_frames": (
            len(
                source_raw_times_array
            )
        ),
        "continuation_raw_frames": (
            len(
                continuation_raw_times_array
            )
        ),
        "continuation_start_ps": float(
            continuation_raw_times_array[
                0
            ]
        ),
        "continuation_end_ps": float(
            continuation_raw_times_array[
                -1
            ]
        ),
        "combined_unique_frames": (
            len(
                combined_time
            )
        ),
        "combined_start_ps": float(
            combined_time[0]
        ),
        "combined_end_ps": float(
            combined_time[-1]
        ),
        "combined_frame_interval_ps": float(
            np.mean(
                combined_intervals
            )
        ),
        "continuation_temperature_mean_K": (
            continuation_temperature_mean
        ),
        "continuation_temperature_std_K": (
            continuation_temperature_std
        ),
        "continuation_temperature_min_K": (
            continuation_temperature_min
        ),
        "continuation_temperature_max_K": (
            continuation_temperature_max
        ),
        "continuation_temperature_slope_K_ps": (
            continuation_temperature_slope
        ),
        "final15_temperature_mean_K": (
            final15_temperature_mean
        ),
        "final15_temperature_std_K": (
            final15_temperature_std
        ),
        "final15_temperature_slope_K_ps": (
            final15_temperature_slope
        ),
        "theoretical_temperature_std_K": (
            THEORETICAL_TEMPERATURE_STD_K
        ),
        "HBN_maximum_displacement_nm": (
            maximum_hbn_displacement
        ),
        "PYR_maximum_displacement_nm": (
            maximum_pyr_displacement
        ),
        "CAPS_maximum_displacement_nm": (
            maximum_caps_displacement
        ),
        "waterO_final_RMS_displacement_nm": (
            final_water_rms
        ),
        "waterO_final_max_displacement_nm": (
            final_water_max
        ),
        "lumen_occupancy_initial": (
            occupancy_initial
        ),
        "lumen_occupancy_mean": (
            occupancy_mean
        ),
        "lumen_occupancy_minimum": (
            occupancy_minimum
        ),
        "lumen_occupancy_maximum": (
            occupancy_maximum
        ),
        "lumen_occupancy_endpoint": (
            occupancy_endpoint
        ),
        "endpoint_initial_lumen_waters_retained": (
            endpoint_retained_initial
        ),
        "endpoint_initial_identity_retention_fraction": (
            endpoint_identity_fraction
        ),
        "final20_occupancy_mean": (
            final20[
                "occupancy_mean"
            ]
        ),
        "final20_occupancy_slope_water_ps": (
            final20[
                "occupancy_slope_water_ps"
            ]
        ),
        "final15_occupancy_mean": (
            final15[
                "occupancy_mean"
            ]
        ),
        "final15_occupancy_slope_water_ps": (
            final15[
                "occupancy_slope_water_ps"
            ]
        ),
        "final10_occupancy_mean": (
            final10[
                "occupancy_mean"
            ]
        ),
        "final10_occupancy_change": (
            final10[
                "occupancy_change"
            ]
        ),
        "final10_occupancy_slope_water_ps": (
            final10[
                "occupancy_slope_water_ps"
            ]
        ),
        "minimum_CAP_OW_distance_nm": (
            minimum_cap_ow_distance_seen
        ),
        "continuation_potential_initial_kJ_mol": float(
            potential[0]
        ),
        "continuation_potential_final_kJ_mol": float(
            potential[-1]
        ),
        "continuation_potential_change_kJ_mol": float(
            potential[-1]
            - potential[0]
        ),
        "continuation_CAP_SOL_LJ_initial_kJ_mol": float(
            cap_sol_lj[0]
        ),
        "continuation_CAP_SOL_LJ_final_kJ_mol": float(
            cap_sol_lj[-1]
        ),
        "continuation_CAP_SOL_LJ_maximum_kJ_mol": float(
            np.max(
                cap_sol_lj
            )
        ),
        "failed_gates": (
            " | ".join(
                failed_gates
            )
        ),
        "R2_vs_R1_comparison_authorized": (
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

    window_lines = "\n".join(
        (
            f"- {row['window_start_ps']:.1f}–"
            f"{row['window_end_ps']:.1f} ps: "
            f"mean={row['occupancy_mean']:.4f}; "
            f"min/max={row['occupancy_minimum']}/"
            f"{row['occupancy_maximum']}; "
            f"change={row['occupancy_change']:+d}; "
            f"slope="
            f"{row['occupancy_slope_water_ps']:.6f} waters/ps"
        )
        for row in window_rows
    )

    block_lines = "\n".join(
        (
            f"- {row['block_start_ps']:.1f}–"
            f"{row['block_end_ps']:.1f} ps: "
            f"mean={row['occupancy_mean']:.4f}; "
            f"min/max={row['occupancy_minimum']}/"
            f"{row['occupancy_maximum']}; "
            f"change={row['occupancy_change']:+d}; "
            f"slope="
            f"{row['occupancy_slope_water_ps']:.6f} waters/ps"
        )
        for row in block_rows
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
        f"""# R2 Frozen-Solute NVT — Combined 0–50 ps Validation

## Execution scope

The accepted 0–20 ps trajectory was not repeated. The simulation was
continued from the exact 20 ps checkpoint to the 50 ps target using the
extended TPR.

- Source checkpoint:
  **step {SOURCE_STEP}, {SOURCE_TIME_PS:.1f} ps**
- Final checkpoint:
  **step {final_checkpoint_step}, {final_checkpoint_time:.6f} ps**
- Velocity regeneration:
  **NO**
- Thermostat-state regeneration:
  **NO**
- Source 0–20 ps rerun:
  **NO**
- Continuation instability signatures:
  **{'NONE' if not instability_hits else ' | '.join(instability_hits)}**

## Trajectory integrity

- Source raw frames:
  **{len(source_raw_times_array)}**
- Continuation raw frames:
  **{len(continuation_raw_times_array)}**
- Continuation time range:
  **{continuation_raw_times_array[0]:.3f}–
  {continuation_raw_times_array[-1]:.3f} ps**
- Combined unique frames:
  **{len(combined_time)}**
- Combined time range:
  **{combined_time[0]:.3f}–{combined_time[-1]:.3f} ps**
- Combined frame interval:
  **{np.mean(combined_intervals):.6f} ps**

## Thermal behavior

- Continuation mean/std/min/max:
  **{continuation_temperature_mean:.4f}/
  {continuation_temperature_std:.4f}/
  {continuation_temperature_min:.4f}/
  {continuation_temperature_max:.4f} K**
- Continuation slope:
  **{continuation_temperature_slope:.6f} K/ps**
- Final 15 ps mean/std/slope:
  **{final15_temperature_mean:.4f}/
  {final15_temperature_std:.4f}/
  {final15_temperature_slope:.6f} K, K, K/ps**
- Theoretical canonical temperature standard deviation:
  **{THEORETICAL_TEMPERATURE_STD_K:.6f} K**

## Frozen-group integrity

Maximum displacement over the combined trajectory:

- HBN:
  **{maximum_hbn_displacement:.12f} nm**
- PYR:
  **{maximum_pyr_displacement:.12f} nm**
- CAPS:
  **{maximum_caps_displacement:.12f} nm**

Water-O final RMS/max displacement:

- **{final_water_rms:.6f}/{final_water_max:.6f} nm**

## Lumen-water behavior

- Initial/mean/minimum/maximum/endpoint occupancy:
  **{occupancy_initial}/
  {occupancy_mean:.4f}/
  {occupancy_minimum}/
  {occupancy_maximum}/
  {occupancy_endpoint}**
- Endpoint retained initial identities:
  **{endpoint_retained_initial}/{occupancy_initial}**
- Endpoint identity-retention fraction:
  **{endpoint_identity_fraction:.6f}**
- Minimum CAP–OW distance:
  **{minimum_cap_ow_distance_seen:.6f} nm**

### Cumulative windows

{window_lines}

### Ten-picosecond blocks

{block_lines}

## Continuation energetics

- Potential initial/final/change:
  **{potential[0]:.6f}/
  {potential[-1]:.6f}/
  {potential[-1] - potential[0]:.6f} kJ/mol**
- CAP–SOL LJ initial/final/maximum:
  **{cap_sol_lj[0]:.6f}/
  {cap_sol_lj[-1]:.6f}/
  {np.max(cap_sol_lj):.6f} kJ/mol**

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- R2-versus-R1 comparison authorized:
  **{'YES' if accepted else 'NO'}**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `{required_next_step}`

R2 remains a neutral frozen steric screening architecture. Even if the
50 ps gate passes, this result does not establish chemical realizability
or authorize long mobile, multitemperature, or quantum calculations.
""",
        encoding="utf-8",
    )

    print()
    print(
        "Day023 R2 checkpoint continuation and "
        "combined 50 ps validation completed."
    )

    print(
        "Mdrun / trajectory-check / checkpoint-dump "
        "return codes: "
        f"{mdrun_return_code}/"
        f"{trajectory_check.returncode}/"
        f"{checkpoint_dump.returncode}"
    )

    print(
        "Completion / checkpoint continuation / "
        "instabilities: "
        f"{'YES' if completion_confirmed else 'NO'} / "
        f"{'YES' if checkpoint_continuation_confirmed else 'NO'} / "
        + (
            "NONE"
            if not instability_hits
            else " | ".join(
                instability_hits
            )
        )
    )

    print(
        "Final checkpoint step / time: "
        f"{final_checkpoint_step}/"
        f"{final_checkpoint_time:.6f} ps"
    )

    print(
        "Source / continuation / combined-unique frames: "
        f"{len(source_raw_times_array)}/"
        f"{len(continuation_raw_times_array)}/"
        f"{len(combined_time)}"
    )

    print(
        "Continuation time range: "
        f"{continuation_raw_times_array[0]:.3f}-"
        f"{continuation_raw_times_array[-1]:.3f} ps"
    )

    print(
        "Combined time range / interval: "
        f"{combined_time[0]:.3f}-"
        f"{combined_time[-1]:.3f} ps / "
        f"{np.mean(combined_intervals):.6f} ps"
    )

    print(
        "Continuation temperature mean/std/min/max/slope: "
        f"{continuation_temperature_mean:.4f}/"
        f"{continuation_temperature_std:.4f}/"
        f"{continuation_temperature_min:.4f}/"
        f"{continuation_temperature_max:.4f}/"
        f"{continuation_temperature_slope:.6f} "
        "K / K / K / K / K ps^-1"
    )

    print(
        "Final-15ps temperature mean/std/slope: "
        f"{final15_temperature_mean:.4f}/"
        f"{final15_temperature_std:.4f}/"
        f"{final15_temperature_slope:.6f} "
        "K / K / K ps^-1"
    )

    print(
        "Maximum HBN/PYR/CAPS displacement: "
        f"{maximum_hbn_displacement:.12f}/"
        f"{maximum_pyr_displacement:.12f}/"
        f"{maximum_caps_displacement:.12f} nm"
    )

    print(
        "Water-O final RMS/max displacement: "
        f"{final_water_rms:.6f}/"
        f"{final_water_max:.6f} nm"
    )

    print(
        "Combined lumen occupancy "
        "initial/mean/min/max/endpoint: "
        f"{occupancy_initial}/"
        f"{occupancy_mean:.4f}/"
        f"{occupancy_minimum}/"
        f"{occupancy_maximum}/"
        f"{occupancy_endpoint}"
    )

    print(
        "Endpoint initial-lumen identity retention: "
        f"{endpoint_retained_initial}/"
        f"{occupancy_initial} "
        f"({endpoint_identity_fraction:.6f})"
    )

    print(
        "Final-20ps occupancy mean/slope: "
        f"{float(final20['occupancy_mean']):.4f}/"
        f"{float(final20['occupancy_slope_water_ps']):.6f} "
        "waters / waters ps^-1"
    )

    print(
        "Final-15ps occupancy mean/slope: "
        f"{float(final15['occupancy_mean']):.4f}/"
        f"{float(final15['occupancy_slope_water_ps']):.6f} "
        "waters / waters ps^-1"
    )

    print(
        "Final-10ps occupancy mean/change/slope: "
        f"{float(final10['occupancy_mean']):.4f}/"
        f"{int(final10['occupancy_change']):+d}/"
        f"{float(final10['occupancy_slope_water_ps']):.6f} "
        "waters / waters / waters ps^-1"
    )

    print(
        "Minimum CAP-OW distance: "
        f"{minimum_cap_ow_distance_seen:.6f} nm"
    )

    print(
        "Continuation potential initial/final/change: "
        f"{potential[0]:.6f}/"
        f"{potential[-1]:.6f}/"
        f"{potential[-1] - potential[0]:.6f} kJ/mol"
    )

    print(
        "Continuation CAP-SOL LJ initial/final/maximum: "
        f"{cap_sol_lj[0]:.6f}/"
        f"{cap_sol_lj[-1]:.6f}/"
        f"{np.max(cap_sol_lj):.6f} kJ/mol"
    )

    for row in block_rows:
        print(
            "Block "
            f"{row['block_start_ps']:.1f}-"
            f"{row['block_end_ps']:.1f} ps | "
            f"mean/min/max="
            f"{row['occupancy_mean']:.4f}/"
            f"{row['occupancy_minimum']}/"
            f"{row['occupancy_maximum']} | "
            f"change="
            f"{row['occupancy_change']:+d} | "
            f"slope="
            f"{row['occupancy_slope_water_ps']:.6f} waters/ps"
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
        "R2-versus-R1 comparison authorized: "
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

    for path in (
        continuation_gro,
        continuation_xtc,
        continuation_edr,
        continuation_cpt,
        COMBINED_OCCUPANCY_CSV,
        OCCUPANCY_WINDOWS_CSV,
        OCCUPANCY_BLOCKS_CSV,
        CONTINUATION_ENERGY_CSV,
        OUTPUT_MANIFEST_JSON,
        SUMMARY_CSV,
        GATES_CSV,
        REPORT_MD,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not accepted:
        raise RuntimeError(
            "R2 combined 50 ps validation requires review."
        )


if __name__ == "__main__":
    main()
