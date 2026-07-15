#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

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

PREPARATION_ROOT = (
    DAY023_ROOT
    / "15_r2_frozen_solute_nvt_20ps_preparation"
)

OUTPUT_ROOT = (
    DAY023_ROOT
    / "16_r2_frozen_solute_nvt_20ps"
)

INPUT_GRO = (
    PREPARATION_ROOT
    / "r2_frozen_solute_nvt_20ps_input.gro"
)

INPUT_TPR = (
    PREPARATION_ROOT
    / "r2_frozen_solute_nvt_20ps.tpr"
)

INPUT_INDEX = (
    PREPARATION_ROOT
    / "r2_frozen_solute_nvt_20ps.ndx"
)

PREPARATION_SUMMARY = (
    PREPARATION_ROOT
    / "r2_frozen_solute_nvt_20ps_preparation_summary.csv"
)

RUN_CONTRACT = (
    PREPARATION_ROOT
    / "r2_frozen_solute_nvt_20ps_run_contract.json"
)

PROTOTYPE_JSON = (
    PROTOTYPE_ROOT
    / "r1_selected_steric_cap_definition.json"
)

DEFFNM = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps"
)

MDRUN_CONSOLE = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_mdrun_console.log"
)

TRAJECTORY_CHECK_LOG = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_trajectory_check.log"
)

ENERGY_MENU_LOG = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_energy_menu.txt"
)

TEMPERATURE_XVG = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_temperature.xvg"
)

POTENTIAL_XVG = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_potential.xvg"
)

TOTAL_ENERGY_XVG = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_total_energy.xvg"
)

CAP_SOL_LJ_XVG = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_cap_sol_lj.xvg"
)

OCCUPANCY_CSV = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_lumen_occupancy.csv"
)

ENERGY_SERIES_CSV = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_energy_series.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_summary.csv"
)

GATES_CSV = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_gates.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_FROZEN_SOLUTE_NVT_20PS_VALIDATION_DAY023.md"
)

EXPECTED_PREPARATION_DECISION = (
    "R2_FROZEN_SOLUTE_NVT_20PS_PREPARED"
)

EXPECTED_ATOMS = 68332

HBN_ATOMS = 1680
PYR_ATOMS = 104
SOLUTE_ATOMS = HBN_ATOMS + PYR_ATOMS

WATERS = 16565
WATER_SITES = 4
WATER_ATOMS = WATERS * WATER_SITES

CAPS_PER_END = 144
TOTAL_CAPS = 288

CAP_START = SOLUTE_ATOMS + WATER_ATOMS
CAP_STOP = CAP_START + TOTAL_CAPS

DT_PS = 0.0005
NSTEPS = 40000
TOTAL_TIME_PS = 20.0

EXPECTED_FRAMES = 41
FRAME_INTERVAL_PS = 0.5

EXPECTED_INITIAL_LUMEN_WATERS = 428

FROZEN_TOLERANCE_NM = 1.0e-6
MIN_WATER_RMS_DISPLACEMENT_NM = 0.05

TEMPERATURE_ANALYSIS_START_PS = 5.0
TEMPERATURE_TARGET_K = 300.0
TEMPERATURE_MEAN_TOLERANCE_K = 5.0

SOL_DEGREES_OF_FREEDOM = 99387.0

THEORETICAL_TEMPERATURE_STD_K = (
    TEMPERATURE_TARGET_K
    * math.sqrt(
        2.0
        / SOL_DEGREES_OF_FREEDOM
    )
)

MAX_POST5_TEMPERATURE_STD_K = (
    3.0
    * THEORETICAL_TEMPERATURE_STD_K
)

MAX_POST5_TEMPERATURE_SLOPE_K_PS = 0.10

MIN_OCCUPANCY_FRACTION = 0.80
MIN_SECOND_HALF_MEAN_OCCUPANCY_FRACTION = 0.90
MIN_ENDPOINT_OCCUPANCY_FRACTION = 0.90
MIN_ENDPOINT_IDENTITY_RETENTION_FRACTION = 0.50

MAX_SECOND_HALF_OCCUPANCY_SLOPE_WATER_PS = 0.50

MIN_CAP_OW_DISTANCE_NM = 0.15
MAX_CAP_SOL_LJ_KJ_MOL = 100.0


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
            f"No rows found in {path}"
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


def minimum_pair_distance(
    first: np.ndarray,
    second: np.ndarray,
    box: np.ndarray,
    *,
    chunk_size: int = 512,
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
            return candidate, terms[candidate]

    for name, number in terms.items():
        if all(
            token in name
            for token in required_tokens
        ):
            return name, number

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


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        INPUT_GRO,
        INPUT_TPR,
        INPUT_INDEX,
        PREPARATION_SUMMARY,
        RUN_CONTRACT,
        PROTOTYPE_JSON,
    ):
        require_file(required)

    preparation = read_single_csv_row(
        PREPARATION_SUMMARY
    )

    contract = json.loads(
        RUN_CONTRACT.read_text(
            encoding="utf-8"
        )
    )

    prototype = json.loads(
        PROTOTYPE_JSON.read_text(
            encoding="utf-8"
        )
    )

    if (
        preparation.get("decision")
        != EXPECTED_PREPARATION_DECISION
    ):
        raise RuntimeError(
            "R2 NVT preparation is not accepted."
        )

    if not parse_bool(
        preparation.get(
            "NVT_execution_authorized",
            "false",
        )
    ):
        raise RuntimeError(
            "R2 NVT execution is not authorized."
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
            "velocity_regeneration_during_mdrun",
            True,
        )
    ):
        raise RuntimeError(
            "Run contract unexpectedly authorizes "
            "velocity regeneration."
        )

    if int(
        contract.get(
            "nsteps",
            -1,
        )
    ) != NSTEPS:
        raise RuntimeError(
            "Run-contract step count is inconsistent."
        )

    reuse_completed_run = (
        os.environ.get(
            "R2_REUSE_COMPLETED_RUN",
            "0",
        )
        == "1"
    )

    existing_products = [
        path
        for path in OUTPUT_ROOT.glob(
            f"{DEFFNM.name}.*"
        )
        if (
            path.is_file()
            and path.stat().st_size > 0
            and path != MDRUN_CONSOLE
        )
    ]

    if existing_products and not reuse_completed_run:
        raise RuntimeError(
            "R2 NVT products already exist. "
            "Set R2_REUSE_COMPLETED_RUN=1 only after "
            "confirming that mdrun completed successfully: "
            + ", ".join(
                path.name
                for path in existing_products
            )
        )

    if reuse_completed_run and not existing_products:
        raise RuntimeError(
            "Reuse mode was requested, but no completed "
            "trajectory products were found."
        )

    gmx = locate_gmx()

    if reuse_completed_run:
        require_file(MDRUN_CONSOLE)

        prior_console = MDRUN_CONSOLE.read_text(
            encoding="utf-8",
            errors="replace",
        )

        if not (
            "Writing final coordinates"
            in prior_console
            or "Finished mdrun"
            in prior_console
        ):
            raise RuntimeError(
                "Existing mdrun output does not confirm completion."
            )

        print(
            "Reusing the completed R2 frozen-solute NVT run."
        )

        print(
            "No molecular-dynamics steps will be repeated."
        )

        print(
            "Proceeding directly to trajectory, energy, "
            "and water-confinement validation."
        )

        mdrun_return_code = 0

    else:
        print(
            "Starting R2 frozen-solute NVT screening run."
        )

        print(
            "HBN, all PYR molecules, and both partial "
            "caps are frozen."
        )

        print(
            "Only TIP4P/2005 water is mobile and "
            "thermostatted."
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
                str(INPUT_TPR),
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

    xtc = Path(
        str(DEFFNM)
        + ".xtc"
    )

    edr = Path(
        str(DEFFNM)
        + ".edr"
    )

    log = Path(
        str(DEFFNM)
        + ".log"
    )

    final_gro = Path(
        str(DEFFNM)
        + ".gro"
    )

    checkpoint = Path(
        str(DEFFNM)
        + ".cpt"
    )

    for product in (
        xtc,
        edr,
        log,
        final_gro,
        checkpoint,
        MDRUN_CONSOLE,
    ):
        require_file(product)

    combined_log = (
        MDRUN_CONSOLE.read_text(
            encoding="utf-8",
            errors="replace",
        )
        + "\n"
        + log.read_text(
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

    checkpoint_reference_confirmed = (
        checkpoint.exists()
        and checkpoint.stat().st_size > 0
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
            str(xtc),
        ],
        cwd=OUTPUT_ROOT,
    )

    TRAJECTORY_CHECK_LOG.write_text(
        trajectory_check.stdout,
        encoding="utf-8",
    )

    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise RuntimeError(
            "MDAnalysis is required for trajectory validation."
        ) from exc

    universe = mda.Universe(
        str(INPUT_GRO),
        str(xtc),
    )

    if universe.atoms.n_atoms != EXPECTED_ATOMS:
        raise RuntimeError(
            "Unexpected trajectory atom count: "
            f"{universe.atoms.n_atoms}/"
            f"{EXPECTED_ATOMS}"
        )

    atom_names = np.asarray(
        universe.atoms.names,
        dtype=object,
    )

    oxygen_indices = np.arange(
        SOLUTE_ATOMS,
        CAP_START,
        WATER_SITES,
        dtype=int,
    )

    if len(oxygen_indices) != WATERS:
        raise RuntimeError(
            "Water-oxygen index accounting failed."
        )

    if not np.all(
        atom_names[
            oxygen_indices
        ]
        == "OW"
    ):
        raise RuntimeError(
            "The expected water-oxygen positions "
            "are not named OW."
        )

    times = []
    occupancy_rows = []

    reference_positions = None
    reference_box = None
    initial_lumen_mask = None
    initial_lumen_count = None

    maximum_hbn_displacement = 0.0
    maximum_pyr_displacement = 0.0
    maximum_caps_displacement = 0.0

    minimum_cap_ow_distance_seen = math.inf

    final_water_rms = math.nan
    final_water_max = math.nan

    for frame_index, timestep in enumerate(
        universe.trajectory
    ):
        current_positions = (
            universe.atoms.positions.astype(
                float,
                copy=True,
            )
            / 10.0
        )

        current_box = (
            np.asarray(
                timestep.dimensions[:3],
                dtype=float,
            )
            / 10.0
        )

        if (
            len(current_box) != 3
            or not np.all(
                np.isfinite(current_box)
            )
            or np.any(current_box <= 0.0)
        ):
            raise RuntimeError(
                f"Invalid box in frame {frame_index}."
            )

        time_ps = float(
            timestep.time
        )

        times.append(time_ps)

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

        current_cap_ow_distance = (
            minimum_pair_distance(
                water_positions,
                cap_positions,
                current_box,
            )
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

        retained_initial_lumen = int(
            np.count_nonzero(
                initial_lumen_mask
                & current_lumen_mask
            )
        )

        occupancy_rows.append(
            {
                "frame": frame_index,
                "time_ps": time_ps,
                "lumen_occupancy": (
                    current_occupancy
                ),
                "initial_lumen_waters_retained": (
                    retained_initial_lumen
                ),
                "initial_lumen_identity_retention_fraction": (
                    retained_initial_lumen
                    / initial_lumen_count
                ),
                "minimum_CAP_OW_distance_nm": (
                    current_cap_ow_distance
                ),
            }
        )

        if frame_index == len(
            universe.trajectory
        ) - 1:
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

    write_csv(
        OCCUPANCY_CSV,
        occupancy_rows,
    )

    times_array = np.asarray(
        times,
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

    retained_identities = np.asarray(
        [
            row[
                "initial_lumen_waters_retained"
            ]
            for row in occupancy_rows
        ],
        dtype=float,
    )

    frame_intervals = np.diff(
        times_array
    )

    second_half_mask = (
        times_array
        >= TOTAL_TIME_PS / 2.0
    )

    second_half_times = (
        times_array[
            second_half_mask
        ]
    )

    second_half_occupancies = (
        occupancies[
            second_half_mask
        ]
    )

    second_half_occupancy_slope = (
        linear_slope(
            second_half_times,
            second_half_occupancies,
        )
    )

    gmx_energy_probe = run_command(
        [
            gmx,
            "energy",
            "-f",
            str(edr),
            "-o",
            str(
                OUTPUT_ROOT
                / "r2_frozen_solute_nvt_20ps_energy_menu_probe.xvg"
            ),
        ],
        cwd=OUTPUT_ROOT,
        input_text="0\n",
    )

    ENERGY_MENU_LOG.write_text(
        gmx_energy_probe.stdout,
        encoding="utf-8",
    )

    terms = parse_energy_menu(
        gmx_energy_probe.stdout
    )

    energy_requests = {
        "Temperature": (
            ("Temperature",),
            ("Temperature",),
            TEMPERATURE_XVG,
        ),
        "Potential": (
            ("Potential",),
            ("Potential",),
            POTENTIAL_XVG,
        ),
        "Total-Energy": (
            ("Total-Energy",),
            ("Total-Energy",),
            TOTAL_ENERGY_XVG,
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
            CAP_SOL_LJ_XVG,
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

        resolved_names[label] = (
            term_name
        )

        extracted[label] = extract_energy_series(
            gmx,
            edr,
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
        ENERGY_SERIES_CSV,
        energy_rows,
    )

    post5_mask = (
        energy_time
        >= TEMPERATURE_ANALYSIS_START_PS
    )

    post5_time = (
        energy_time[
            post5_mask
        ]
    )

    post5_temperature = (
        temperature[
            post5_mask
        ]
    )

    last10_mask = (
        energy_time >= 10.0
    )

    last10_time = (
        energy_time[
            last10_mask
        ]
    )

    last10_temperature = (
        temperature[
            last10_mask
        ]
    )

    temperature_mean = float(
        np.mean(temperature)
    )

    temperature_std = float(
        np.std(
            temperature,
            ddof=1,
        )
    )

    temperature_min = float(
        np.min(temperature)
    )

    temperature_max = float(
        np.max(temperature)
    )

    post5_temperature_mean = float(
        np.mean(
            post5_temperature
        )
    )

    post5_temperature_std = float(
        np.std(
            post5_temperature,
            ddof=1,
        )
    )

    post5_temperature_min = float(
        np.min(
            post5_temperature
        )
    )

    post5_temperature_max = float(
        np.max(
            post5_temperature
        )
    )

    post5_temperature_slope = linear_slope(
        post5_time,
        post5_temperature,
    )

    last10_temperature_mean = float(
        np.mean(
            last10_temperature
        )
    )

    last10_temperature_std = float(
        np.std(
            last10_temperature,
            ddof=1,
        )
    )

    last10_temperature_slope = linear_slope(
        last10_time,
        last10_temperature,
    )

    occupancy_initial = int(
        occupancies[0]
    )

    occupancy_mean = float(
        np.mean(
            occupancies
        )
    )

    occupancy_minimum = int(
        np.min(
            occupancies
        )
    )

    occupancy_maximum = int(
        np.max(
            occupancies
        )
    )

    occupancy_endpoint = int(
        occupancies[-1]
    )

    second_half_occupancy_mean = float(
        np.mean(
            second_half_occupancies
        )
    )

    endpoint_retained_initial = int(
        retained_identities[-1]
    )

    minimum_allowed_occupancy = int(
        math.floor(
            MIN_OCCUPANCY_FRACTION
            * EXPECTED_INITIAL_LUMEN_WATERS
        )
    )

    minimum_second_half_mean = (
        MIN_SECOND_HALF_MEAN_OCCUPANCY_FRACTION
        * EXPECTED_INITIAL_LUMEN_WATERS
    )

    minimum_endpoint_occupancy = int(
        math.floor(
            MIN_ENDPOINT_OCCUPANCY_FRACTION
            * EXPECTED_INITIAL_LUMEN_WATERS
        )
    )

    minimum_endpoint_identity_retention = int(
        math.floor(
            MIN_ENDPOINT_IDENTITY_RETENTION_FRACTION
            * EXPECTED_INITIAL_LUMEN_WATERS
        )
    )

    gates = {
        "preparation_decision_is_valid": (
            preparation.get("decision")
            == EXPECTED_PREPARATION_DECISION
        ),
        "run_contract_authorized_execution": (
            bool(
                contract.get(
                    "execution_authorized",
                    False,
                )
            )
        ),
        "velocity_regeneration_remained_disabled": (
            not bool(
                contract.get(
                    "velocity_regeneration_during_mdrun",
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
        "checkpoint_was_written": (
            checkpoint_reference_confirmed
        ),
        "no_instability_signatures": (
            len(
                instability_hits
            )
            == 0
        ),
        "trajectory_check_return_code_zero": (
            trajectory_check.returncode
            == 0
        ),
        "trajectory_has_41_frames": (
            len(times_array)
            == EXPECTED_FRAMES
        ),
        "trajectory_starts_at_0ps": (
            math.isclose(
                float(times_array[0]),
                0.0,
                rel_tol=0.0,
                abs_tol=1.0e-4,
            )
        ),
        "trajectory_ends_at_20ps": (
            math.isclose(
                float(times_array[-1]),
                TOTAL_TIME_PS,
                rel_tol=0.0,
                abs_tol=1.0e-3,
            )
        ),
        "trajectory_interval_is_0p5ps": (
            len(frame_intervals) > 0
            and np.allclose(
                frame_intervals,
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
        "water_is_mobile": (
            math.isfinite(
                final_water_rms
            )
            and final_water_rms
            >= MIN_WATER_RMS_DISPLACEMENT_NM
        ),
        "water_displacement_is_finite": (
            math.isfinite(
                final_water_rms
            )
            and math.isfinite(
                final_water_max
            )
        ),
        "post5_temperature_mean_is_295_to_305K": (
            abs(
                post5_temperature_mean
                - TEMPERATURE_TARGET_K
            )
            <= TEMPERATURE_MEAN_TOLERANCE_K
        ),
        "post5_temperature_std_is_acceptable": (
            post5_temperature_std
            <= MAX_POST5_TEMPERATURE_STD_K
        ),
        "post5_temperature_slope_is_acceptable": (
            abs(
                post5_temperature_slope
            )
            <= MAX_POST5_TEMPERATURE_SLOPE_K_PS
        ),
        "initial_lumen_occupancy_is_428": (
            occupancy_initial
            == EXPECTED_INITIAL_LUMEN_WATERS
        ),
        "trajectory_does_not_approach_complete_drying": (
            occupancy_minimum
            >= minimum_allowed_occupancy
        ),
        "second_half_mean_occupancy_is_at_least_90_percent": (
            second_half_occupancy_mean
            >= minimum_second_half_mean
        ),
        "endpoint_occupancy_is_at_least_90_percent": (
            occupancy_endpoint
            >= minimum_endpoint_occupancy
        ),
        "endpoint_identity_retention_is_at_least_50_percent": (
            endpoint_retained_initial
            >= minimum_endpoint_identity_retention
        ),
        "second_half_occupancy_slope_is_acceptable": (
            abs(
                second_half_occupancy_slope
            )
            <= MAX_SECOND_HALF_OCCUPANCY_SLOPE_WATER_PS
        ),
        "minimum_CAP_OW_distance_is_safe": (
            minimum_cap_ow_distance_seen
            >= MIN_CAP_OW_DISTANCE_NM
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
            <= MAX_CAP_SOL_LJ_KJ_MOL
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
        "R2_FROZEN_SOLUTE_NVT_20PS_VALIDATED"
        if accepted
        else
        "R2_FROZEN_SOLUTE_NVT_20PS_REQUIRES_REVIEW"
    )

    required_next_step = (
        "PREPARE_R2_30PS_CHECKPOINT_CONTINUATION_TO_50PS"
        if accepted
        else
        "REVIEW_R2_FROZEN_SOLUTE_NVT_20PS_GATE_FAILURES"
    )

    summary = {
        "decision": decision,
        "mdrun_return_code": (
            mdrun_return_code
        ),
        "trajectory_check_return_code": (
            trajectory_check.returncode
        ),
        "completion_confirmed": (
            completion_confirmed
        ),
        "checkpoint_written": (
            checkpoint_reference_confirmed
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
        "trajectory_frames": (
            len(times_array)
        ),
        "trajectory_start_ps": float(
            times_array[0]
        ),
        "trajectory_end_ps": float(
            times_array[-1]
        ),
        "trajectory_interval_ps": float(
            np.mean(
                frame_intervals
            )
        ),
        "temperature_full_mean_K": (
            temperature_mean
        ),
        "temperature_full_std_K": (
            temperature_std
        ),
        "temperature_full_min_K": (
            temperature_min
        ),
        "temperature_full_max_K": (
            temperature_max
        ),
        "temperature_post5_mean_K": (
            post5_temperature_mean
        ),
        "temperature_post5_std_K": (
            post5_temperature_std
        ),
        "temperature_post5_min_K": (
            post5_temperature_min
        ),
        "temperature_post5_max_K": (
            post5_temperature_max
        ),
        "temperature_post5_slope_K_ps": (
            post5_temperature_slope
        ),
        "temperature_last10_mean_K": (
            last10_temperature_mean
        ),
        "temperature_last10_std_K": (
            last10_temperature_std
        ),
        "temperature_last10_slope_K_ps": (
            last10_temperature_slope
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
        "second_half_occupancy_mean": (
            second_half_occupancy_mean
        ),
        "second_half_occupancy_slope_water_ps": (
            second_half_occupancy_slope
        ),
        "endpoint_initial_lumen_waters_retained": (
            endpoint_retained_initial
        ),
        "endpoint_initial_lumen_identity_retention_fraction": (
            endpoint_retained_initial
            / occupancy_initial
        ),
        "minimum_CAP_OW_distance_nm": (
            minimum_cap_ow_distance_seen
        ),
        "potential_initial_kJ_mol": float(
            potential[0]
        ),
        "potential_final_kJ_mol": float(
            potential[-1]
        ),
        "potential_change_kJ_mol": float(
            potential[-1]
            - potential[0]
        ),
        "CAP_SOL_LJ_initial_kJ_mol": float(
            cap_sol_lj[0]
        ),
        "CAP_SOL_LJ_final_kJ_mol": float(
            cap_sol_lj[-1]
        ),
        "CAP_SOL_LJ_maximum_kJ_mol": float(
            np.max(
                cap_sol_lj
            )
        ),
        "failed_gates": (
            " | ".join(
                failed_gates
            )
        ),
        "checkpoint_continuation_preparation_authorized": (
            accepted
        ),
        "checkpoint_continuation_execution_authorized": False,
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

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed in gates.items()
    )

    REPORT_MD.write_text(
        f"""# R2 Frozen-Solute NVT — 20 ps Validation

## Scope

R2 was screened for 20 ps with HBN, all four pyrenes, and both
partial-cap assemblies frozen. Only TIP4P/2005 water was mobile and
thermostatted.

## Execution

- Mdrun return code:
  **{mdrun_return_code}**
- Trajectory-check return code:
  **{trajectory_check.returncode}**
- Completion confirmed:
  **{'YES' if completion_confirmed else 'NO'}**
- Checkpoint written:
  **{'YES' if checkpoint_reference_confirmed else 'NO'}**
- Instability signatures:
  **{'NONE' if not instability_hits else ' | '.join(instability_hits)}**
- Frames:
  **{len(times_array)}**
- Time range:
  **{times_array[0]:.3f}–{times_array[-1]:.3f} ps**
- Mean frame interval:
  **{np.mean(frame_intervals):.6f} ps**

## Thermal behavior

- Full mean/std/min/max:
  **{temperature_mean:.4f}/
  {temperature_std:.4f}/
  {temperature_min:.4f}/
  {temperature_max:.4f} K**
- 5–20 ps mean/std/min/max:
  **{post5_temperature_mean:.4f}/
  {post5_temperature_std:.4f}/
  {post5_temperature_min:.4f}/
  {post5_temperature_max:.4f} K**
- 5–20 ps slope:
  **{post5_temperature_slope:.6f} K/ps**
- Theoretical canonical standard deviation:
  **{THEORETICAL_TEMPERATURE_STD_K:.6f} K**
- 10–20 ps mean/std/slope:
  **{last10_temperature_mean:.4f}/
  {last10_temperature_std:.4f}/
  {last10_temperature_slope:.6f} K, K, K/ps**

## Frozen-group integrity

Maximum displacement over all frames:

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
- Second-half mean:
  **{second_half_occupancy_mean:.4f} waters**
- Second-half slope:
  **{second_half_occupancy_slope:.6f} waters/ps**
- Endpoint retention of initially luminal waters:
  **{endpoint_retained_initial}/{occupancy_initial}**
- Endpoint initial-identity retention fraction:
  **{endpoint_retained_initial / occupancy_initial:.6f}**
- Minimum CAP–OW distance:
  **{minimum_cap_ow_distance_seen:.6f} nm**

## Energetics

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
- Checkpoint-continuation preparation authorized:
  **{'YES' if accepted else 'NO'}**
- Checkpoint-continuation execution authorized:
  **NO**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `{required_next_step}`

R2 remains a neutral frozen steric screening architecture. Passing
this gate would justify only preparation of the matched checkpoint
continuation from 20 to 50 ps.
""",
        encoding="utf-8",
    )

    print()
    print(
        "Day023 R2 frozen-solute NVT 20 ps "
        "execution and validation completed."
    )

    print(
        "Mdrun / trajectory-check return codes: "
        f"{mdrun_return_code}/"
        f"{trajectory_check.returncode}"
    )

    print(
        "Finished / checkpoint / instabilities: "
        f"{'YES' if completion_confirmed else 'NO'} / "
        f"{'YES' if checkpoint_reference_confirmed else 'NO'} / "
        + (
            "NONE"
            if not instability_hits
            else " | ".join(
                instability_hits
            )
        )
    )

    print(
        "Trajectory frames / time range / interval: "
        f"{len(times_array)} / "
        f"{times_array[0]:.3f}-"
        f"{times_array[-1]:.3f} ps / "
        f"{np.mean(frame_intervals):.6f} ps"
    )

    print(
        "Full temperature mean/std/min/max: "
        f"{temperature_mean:.4f}/"
        f"{temperature_std:.4f}/"
        f"{temperature_min:.4f}/"
        f"{temperature_max:.4f} K"
    )

    print(
        "Post-5ps temperature mean/std/min/max/slope: "
        f"{post5_temperature_mean:.4f}/"
        f"{post5_temperature_std:.4f}/"
        f"{post5_temperature_min:.4f}/"
        f"{post5_temperature_max:.4f}/"
        f"{post5_temperature_slope:.6f} "
        "K / K / K / K / K ps^-1"
    )

    print(
        "Theoretical canonical temperature std: "
        f"{THEORETICAL_TEMPERATURE_STD_K:.6f} K"
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
        "Lumen occupancy initial/mean/min/max/endpoint: "
        f"{occupancy_initial}/"
        f"{occupancy_mean:.4f}/"
        f"{occupancy_minimum}/"
        f"{occupancy_maximum}/"
        f"{occupancy_endpoint}"
    )

    print(
        "Second-half occupancy mean/slope: "
        f"{second_half_occupancy_mean:.4f}/"
        f"{second_half_occupancy_slope:.6f} "
        "waters / waters ps^-1"
    )

    print(
        "Endpoint initial-lumen identity retention: "
        f"{endpoint_retained_initial}/"
        f"{occupancy_initial} "
        f"({endpoint_retained_initial / occupancy_initial:.6f})"
    )

    print(
        "Minimum CAP-OW distance: "
        f"{minimum_cap_ow_distance_seen:.6f} nm"
    )

    print(
        "Potential initial/final/change: "
        f"{potential[0]:.6f}/"
        f"{potential[-1]:.6f}/"
        f"{potential[-1] - potential[0]:.6f} kJ/mol"
    )

    print(
        "CAP-SOL LJ initial/final/maximum: "
        f"{cap_sol_lj[0]:.6f}/"
        f"{cap_sol_lj[-1]:.6f}/"
        f"{np.max(cap_sol_lj):.6f} kJ/mol"
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
        "Checkpoint-continuation preparation authorized: "
        f"{'YES' if accepted else 'NO'}"
    )

    print(
        "Checkpoint-continuation execution authorized: NO"
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
        final_gro,
        xtc,
        edr,
        checkpoint,
        OCCUPANCY_CSV,
        ENERGY_SERIES_CSV,
        SUMMARY_CSV,
        GATES_CSV,
        REPORT_MD,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not accepted:
        raise RuntimeError(
            "R2 frozen-solute NVT 20 ps "
            "requires review."
        )


if __name__ == "__main__":
    main()
