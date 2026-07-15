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

STAGE06 = "06_nvt_unrestrained_10ps"
STAGE07 = "07_nvt_unrestrained_25ps"

RUN06 = PROTOCOL / "execution" / STAGE06
RUN07 = PROTOCOL / "execution" / STAGE07

REFERENCE_GRO = RUN06 / f"{STAGE06}.gro"

TRAJECTORY = RUN07 / f"{STAGE07}.xtc"
TPR = RUN07 / f"{STAGE07}.tpr"
EDR = RUN07 / f"{STAGE07}.edr"

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
    "validate_day022_stage07_mobile_25ps.py"
)

IMPROPER_MODULE_PATH = (
    ROOT
    / "scripts/phase1A/"
    "diagnose_day022_stage07_hbn_improper_phase.py"
)

ANALYSIS_ROOT = (
    RUN07
    / "time_resolved_stability"
)

FRAME_ROOT = (
    ANALYSIS_ROOT
    / "solute_waterO_frames"
)

INDEX_FILE = (
    ANALYSIS_ROOT
    / "hbn_pyr_waterO.ndx"
)

TRJCONV_LOG = (
    ANALYSIS_ROOT
    / "stage07_frame_extraction.log"
)

PER_FRAME_CSV = (
    ANALYSIS_ROOT
    / "stage07_time_resolved_metrics.csv"
)

BLOCK_CSV = (
    ANALYSIS_ROOT
    / "stage07_5ps_block_metrics.csv"
)

SUMMARY_CSV = (
    ANALYSIS_ROOT
    / "stage07_time_resolved_stability_summary.csv"
)

REPORT_MD = (
    ANALYSIS_ROOT
    / "STAGE07_TIME_RESOLVED_STABILITY_DAY022.md"
)

HBN_COUNT = 1680
PYR_COUNT = 4
PYR_ATOMS = 26
SOLUTE_COUNT = HBN_COUNT + PYR_COUNT * PYR_ATOMS

TOTAL_ATOMS = 68320
FIRST_WATER_ATOM = SOLUTE_COUNT + 1
WATER_ATOMS_PER_MOLECULE = 4
BLOCK_COUNT = 5


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
    )


def load_module(
    module_path: Path,
    module_name: str,
):
    specification = (
        importlib.util.spec_from_file_location(
            module_name,
            module_path,
        )
    )

    if (
        specification is None
        or specification.loader is None
    ):
        raise RuntimeError(
            f"Could not load module: {module_path}"
        )

    module = importlib.util.module_from_spec(
        specification
    )

    specification.loader.exec_module(module)

    return module


def require_inputs() -> None:
    required = (
        REFERENCE_GRO,
        TRAJECTORY,
        TPR,
        EDR,
        HBN_ITP,
        PYR_ITP,
        INTEGRATED_MODULE_PATH,
        IMPROPER_MODULE_PATH,
    )

    missing = [
        path
        for path in required
        if (
            not path.exists()
            or path.stat().st_size == 0
        )
    ]

    if missing:
        raise RuntimeError(
            "Missing or empty required files:\n"
            + "\n".join(
                str(path)
                for path in missing
            )
        )


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
        "Could not locate the GROMACS executable"
    )


def write_index() -> None:
    solute = list(
        range(
            1,
            SOLUTE_COUNT + 1,
        )
    )

    water_oxygen = list(
        range(
            FIRST_WATER_ATOM,
            TOTAL_ATOMS + 1,
            WATER_ATOMS_PER_MOLECULE,
        )
    )

    indices = solute + water_oxygen

    with INDEX_FILE.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "[ HBN_PYR_WATERO ]\n"
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


def natural_frame_key(path: Path) -> int:
    match = re.search(
        r"(\d+)(?=\.gro$)",
        path.name,
    )

    if match is None:
        return -1

    return int(
        match.group(1)
    )


def extract_frames(
    gmx: str,
) -> list[Path]:
    FRAME_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for old_frame in FRAME_ROOT.glob(
        "stage07_frame*.gro"
    ):
        old_frame.unlink()

    output_template = (
        FRAME_ROOT
        / "stage07_frame.gro"
    )

    command = [
        gmx,
        "trjconv",
        "-f",
        str(TRAJECTORY),
        "-s",
        str(TPR),
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
        cwd=RUN07,
        env={
            **os.environ,
            "GMX_MAXBACKUP": "-1",
        },
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
            "GROMACS frame extraction failed. "
            f"See {TRJCONV_LOG}"
        )

    frames = sorted(
        FRAME_ROOT.glob(
            "stage07_frame*.gro"
        ),
        key=natural_frame_key,
    )

    if len(frames) < 40:
        raise RuntimeError(
            "Too few frames were extracted: "
            f"{len(frames)}"
        )

    return frames


def read_gro(
    path: Path,
    parse_count: int | None = None,
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

    count = (
        natoms
        if parse_count is None
        else min(
            natoms,
            parse_count,
        )
    )

    atoms = []

    for local_index, line in enumerate(
        lines[2 : 2 + count],
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
            f"Invalid box in {path}"
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
        r"\bt\s*=\s*"
        r"([-+0-9.eE]+)",
        title,
    )

    if match is None:
        return fallback

    return float(
        match.group(1)
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


def extract_energy_term(
    gmx: str,
    term: str,
) -> tuple[np.ndarray, np.ndarray]:
    safe_name = (
        term.lower()
        .replace(
            "-",
            "_",
        )
        .replace(
            " ",
            "_",
        )
    )

    output = (
        ANALYSIS_ROOT
        / f"stage07_{safe_name}.xvg"
    )

    command = [
        gmx,
        "energy",
        "-f",
        str(EDR),
        "-o",
        str(output),
    ]

    completed = subprocess.run(
        command,
        input=f"{term}\n0\n",
        cwd=RUN07,
        env={
            **os.environ,
            "GMX_MAXBACKUP": "-1",
        },
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )

    if (
        completed.returncode != 0
        or not output.exists()
    ):
        raise RuntimeError(
            f"Could not extract energy term: {term}"
        )

    rows = []

    for raw_line in output.read_text(
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

    if not rows:
        raise RuntimeError(
            f"No numeric values in {output}"
        )

    data = np.array(
        rows,
        dtype=float,
    )

    times = (
        data[:, 0]
        - data[0, 0]
    )

    return (
        times,
        data[:, 1],
    )


def block_mean(
    times: np.ndarray,
    values: np.ndarray,
    lower: float,
    upper: float,
    last_block: bool,
) -> float:
    if last_block:
        mask = (
            (times >= lower)
            & (times <= upper)
        )
    else:
        mask = (
            (times >= lower)
            & (times < upper)
        )

    if not np.any(mask):
        return math.nan

    return float(
        np.mean(
            values[mask]
        )
    )


def main() -> None:
    require_inputs()

    ANALYSIS_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    integrated = load_module(
        INTEGRATED_MODULE_PATH,
        "stage07_integrated_module",
    )

    improper_module = load_module(
        IMPROPER_MODULE_PATH,
        "stage07_improper_module",
    )

    gmx = locate_gmx()

    write_index()

    frames = extract_frames(
        gmx
    )

    (
        _,
        reference_atoms,
        reference_box,
    ) = read_gro(
        REFERENCE_GRO,
        parse_count=SOLUTE_COUNT,
    )

    reference_positions = np.array(
        [
            atom["position"]
            for atom in reference_atoms
        ],
        dtype=float,
    )

    hbn_topology = integrated.parse_itp(
        HBN_ITP
    )

    pyr_topology = integrated.parse_itp(
        PYR_ITP
    )

    impropers = improper_module.parse_impropers(
        HBN_ITP
    )

    hbn_reference_wrapped = (
        reference_positions[:HBN_COUNT]
    )

    hbn_reference_unwrapped = (
        integrated.unwrap_by_bonds(
            hbn_reference_wrapped,
            reference_box,
            hbn_topology["bonds"],
        )
    )

    pyr_reference_unwrapped = []

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

        pyr_reference_unwrapped.append(
            integrated.unwrap_by_bonds(
                reference_positions[
                    start:stop
                ],
                reference_box,
                pyr_topology["bonds"],
            )
        )

    hbn_bond_equilibrium = (
        integrated.finite_equilibria(
            hbn_topology["bonds"]
        )
    )

    hbn_angle_equilibrium = (
        integrated.finite_equilibria(
            hbn_topology["angles"]
        )
    )

    pyr_bond_equilibrium = (
        integrated.finite_equilibria(
            pyr_topology["bonds"]
        )
    )

    pyr_angle_equilibrium = (
        integrated.finite_equilibria(
            pyr_topology["angles"]
        )
    )

    per_frame_rows = []
    planarity_matrix = []
    equilibrium_matrix = []

    top_count = max(
        1,
        math.ceil(
            0.01
            * len(impropers)
        ),
    )

    top_sets = []

    for frame_index, frame_path in enumerate(
        frames
    ):
        title, atoms, box = read_gro(
            frame_path
        )

        expected_selected_atoms = (
            SOLUTE_COUNT
            + (
                TOTAL_ATOMS
                - FIRST_WATER_ATOM
            )
            // WATER_ATOMS_PER_MOLECULE
            + 1
        )

        if len(atoms) != expected_selected_atoms:
            raise RuntimeError(
                f"Unexpected selected atom count "
                f"in {frame_path}: {len(atoms)}"
            )

        positions = np.array(
            [
                atom["position"]
                for atom in atoms
            ],
            dtype=float,
        )

        absolute_time = parse_time(
            title,
            fallback=0.5 * frame_index,
        )

        hbn_wrapped = positions[
            :HBN_COUNT
        ]

        hbn_unwrapped = (
            integrated.unwrap_by_bonds(
                hbn_wrapped,
                box,
                hbn_topology["bonds"],
            )
        )

        hbn_residuals, _ = (
            integrated.kabsch_metrics(
                hbn_reference_unwrapped,
                hbn_unwrapped,
            )
        )

        hbn_aligned_rms = float(
            np.sqrt(
                np.mean(
                    hbn_residuals ** 2
                )
            )
        )

        hbn_bonds = (
            integrated.bond_lengths(
                hbn_unwrapped,
                hbn_topology["bonds"],
            )
        )

        hbn_bond_deviation = np.abs(
            hbn_bonds
            - hbn_bond_equilibrium
        )

        hbn_angles = (
            integrated.angle_values(
                hbn_unwrapped,
                hbn_topology["angles"],
            )
        )

        hbn_angle_deviation = np.abs(
            hbn_angles
            - hbn_angle_equilibrium
        )

        improper_angles = np.array(
            [
                improper_module.dihedral_angle(
                    hbn_wrapped,
                    box,
                    improper,
                )
                for improper in impropers
            ],
            dtype=float,
        )

        planarity = (
            improper_module.planarity_deviation(
                improper_angles
            )
        )

        transformed_angles = (
            improper_module.wrap_degrees(
                -improper_angles
                + 180.0
            )
        )

        equilibrium_residual = (
            improper_module.equilibrium_residual(
                transformed_angles,
                impropers,
            )
        )

        planarity_matrix.append(
            planarity
        )

        equilibrium_matrix.append(
            equilibrium_residual
        )

        top_sets.append(
            set(
                np.argsort(
                    planarity
                )[::-1][
                    :top_count
                ].tolist()
            )
        )

        max_pyr_aligned_rms = 0.0
        max_pyr_bond_deviation = 0.0
        max_pyr_angle_deviation = 0.0

        pyrene_heavy_positions = []
        pyrene_molecule_heavy = []

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

            current_wrapped = positions[
                start:stop
            ]

            current_unwrapped = (
                integrated.unwrap_by_bonds(
                    current_wrapped,
                    box,
                    pyr_topology["bonds"],
                )
            )

            residuals, _ = (
                integrated.kabsch_metrics(
                    pyr_reference_unwrapped[
                        molecule_index
                    ],
                    current_unwrapped,
                )
            )

            aligned_rms = float(
                np.sqrt(
                    np.mean(
                        residuals ** 2
                    )
                )
            )

            max_pyr_aligned_rms = max(
                max_pyr_aligned_rms,
                aligned_rms,
            )

            current_bonds = (
                integrated.bond_lengths(
                    current_unwrapped,
                    pyr_topology["bonds"],
                )
            )

            current_angles = (
                integrated.angle_values(
                    current_unwrapped,
                    pyr_topology["angles"],
                )
            )

            max_pyr_bond_deviation = max(
                max_pyr_bond_deviation,
                float(
                    np.max(
                        np.abs(
                            current_bonds
                            - pyr_bond_equilibrium
                        )
                    )
                ),
            )

            max_pyr_angle_deviation = max(
                max_pyr_angle_deviation,
                float(
                    np.max(
                        np.abs(
                            current_angles
                            - pyr_angle_equilibrium
                        )
                    )
                ),
            )

            heavy_indices = [
                local_index
                for local_index in range(
                    start,
                    stop,
                )
                if not str(
                    atoms[local_index][
                        "atom_name"
                    ]
                ).upper().startswith("H")
            ]

            molecule_heavy = positions[
                heavy_indices
            ]

            pyrene_molecule_heavy.append(
                molecule_heavy
            )

            pyrene_heavy_positions.append(
                molecule_heavy
            )

        pyrene_heavy = np.concatenate(
            pyrene_heavy_positions,
            axis=0,
        )

        water_oxygen = positions[
            SOLUTE_COUNT:
        ]

        hbn_pyr_contact = (
            integrated.pair_contact_metrics(
                hbn_wrapped,
                pyrene_heavy,
                box,
            )["minimum_nm"]
        )

        pyr_water_contact = (
            integrated.pair_contact_metrics(
                pyrene_heavy,
                water_oxygen,
                box,
            )["minimum_nm"]
        )

        pyr_pyr_contact = math.inf

        for first in range(
            PYR_COUNT
        ):
            for second in range(
                first + 1,
                PYR_COUNT,
            ):
                contact = (
                    integrated.pair_contact_metrics(
                        pyrene_molecule_heavy[
                            first
                        ],
                        pyrene_molecule_heavy[
                            second
                        ],
                        box,
                    )["minimum_nm"]
                )

                pyr_pyr_contact = min(
                    pyr_pyr_contact,
                    float(contact),
                )

        per_frame_rows.append(
            {
                "frame": frame_index,
                "absolute_time_ps": (
                    absolute_time
                ),
                "relative_time_ps": 0.0,
                "HBN_aligned_rms_nm": (
                    hbn_aligned_rms
                ),
                "HBN_aligned_max_nm": (
                    float(
                        hbn_residuals.max()
                    )
                ),
                "HBN_bond_deviation_q99_nm": (
                    q(
                        hbn_bond_deviation,
                        0.99,
                    )
                ),
                "HBN_bond_deviation_max_nm": (
                    float(
                        hbn_bond_deviation.max()
                    )
                ),
                "HBN_angle_deviation_q99_deg": (
                    q(
                        hbn_angle_deviation,
                        0.99,
                    )
                ),
                "HBN_angle_deviation_max_deg": (
                    float(
                        hbn_angle_deviation.max()
                    )
                ),
                "HBN_planarity_q95_deg": (
                    q(
                        planarity,
                        0.95,
                    )
                ),
                "HBN_planarity_q99_deg": (
                    q(
                        planarity,
                        0.99,
                    )
                ),
                "HBN_planarity_max_deg": (
                    float(
                        planarity.max()
                    )
                ),
                "HBN_equilibrium_q99_deg": (
                    q(
                        equilibrium_residual,
                        0.99,
                    )
                ),
                "HBN_equilibrium_max_deg": (
                    float(
                        equilibrium_residual.max()
                    )
                ),
                "PYR_max_aligned_rms_nm": (
                    max_pyr_aligned_rms
                ),
                "PYR_max_bond_deviation_nm": (
                    max_pyr_bond_deviation
                ),
                "PYR_max_angle_deviation_deg": (
                    max_pyr_angle_deviation
                ),
                "HBN_PYR_minimum_contact_nm": (
                    float(
                        hbn_pyr_contact
                    )
                ),
                "PYR_PYR_minimum_contact_nm": (
                    float(
                        pyr_pyr_contact
                    )
                ),
                "PYR_waterO_minimum_contact_nm": (
                    float(
                        pyr_water_contact
                    )
                ),
            }
        )

    initial_time = float(
        per_frame_rows[0][
            "absolute_time_ps"
        ]
    )

    for row in per_frame_rows:
        row["relative_time_ps"] = (
            float(
                row["absolute_time_ps"]
            )
            - initial_time
        )

    times = np.array(
        [
            float(
                row["relative_time_ps"]
            )
            for row in per_frame_rows
        ],
        dtype=float,
    )

    planarity_array = np.array(
        planarity_matrix,
        dtype=float,
    )

    equilibrium_array = np.array(
        equilibrium_matrix,
        dtype=float,
    )

    temperature_time, temperature = (
        extract_energy_term(
            gmx,
            "Temperature",
        )
    )

    potential_time, potential = (
        extract_energy_term(
            gmx,
            "Potential",
        )
    )

    pressure_time, pressure = (
        extract_energy_term(
            gmx,
            "Pressure",
        )
    )

    duration = float(
        times.max()
    )

    edges = np.linspace(
        0.0,
        duration,
        BLOCK_COUNT + 1,
    )

    block_rows = []
    block_top_sets = []

    metric_names = [
        "HBN_aligned_rms_nm",
        "HBN_bond_deviation_q99_nm",
        "HBN_bond_deviation_max_nm",
        "HBN_angle_deviation_q99_deg",
        "HBN_angle_deviation_max_deg",
        "HBN_planarity_q99_deg",
        "HBN_planarity_max_deg",
        "HBN_equilibrium_q99_deg",
        "HBN_equilibrium_max_deg",
        "PYR_max_aligned_rms_nm",
        "PYR_max_bond_deviation_nm",
        "PYR_max_angle_deviation_deg",
        "HBN_PYR_minimum_contact_nm",
        "PYR_PYR_minimum_contact_nm",
        "PYR_waterO_minimum_contact_nm",
    ]

    for block_index in range(
        BLOCK_COUNT
    ):
        lower = float(
            edges[block_index]
        )

        upper = float(
            edges[block_index + 1]
        )

        last_block = (
            block_index
            == BLOCK_COUNT - 1
        )

        if last_block:
            frame_mask = (
                (times >= lower)
                & (times <= upper)
            )
        else:
            frame_mask = (
                (times >= lower)
                & (times < upper)
            )

        indices = np.flatnonzero(
            frame_mask
        )

        if len(indices) == 0:
            raise RuntimeError(
                f"No frames in block {block_index + 1}"
            )

        row = {
            "block": block_index + 1,
            "start_ps": lower,
            "end_ps": upper,
            "frame_count": len(indices),
            "temperature_mean_K": block_mean(
                temperature_time,
                temperature,
                lower,
                upper,
                last_block,
            ),
            "potential_mean_kJ_mol": block_mean(
                potential_time,
                potential,
                lower,
                upper,
                last_block,
            ),
            "pressure_mean_bar": block_mean(
                pressure_time,
                pressure,
                lower,
                upper,
                last_block,
            ),
        }

        for metric_name in metric_names:
            values = np.array(
                [
                    float(
                        per_frame_rows[index][
                            metric_name
                        ]
                    )
                    for index in indices
                ],
                dtype=float,
            )

            row[
                f"{metric_name}_mean"
            ] = float(
                np.mean(values)
            )

            if "minimum_contact" in metric_name:
                row[
                    f"{metric_name}_extreme"
                ] = float(
                    np.min(values)
                )
            else:
                row[
                    f"{metric_name}_extreme"
                ] = float(
                    np.max(values)
                )

        block_mean_planarity = np.mean(
            planarity_array[
                indices
            ],
            axis=0,
        )

        block_top_sets.append(
            set(
                np.argsort(
                    block_mean_planarity
                )[::-1][
                    :top_count
                ].tolist()
            )
        )

        block_rows.append(row)

    above20_fraction = np.mean(
        planarity_array > 20.0,
        axis=0,
    )

    persistent_50 = int(
        np.count_nonzero(
            above20_fraction >= 0.50
        )
    )

    persistent_80 = int(
        np.count_nonzero(
            above20_fraction >= 0.80
        )
    )

    persistent_80_indices = np.flatnonzero(
        above20_fraction >= 0.80
    )

    persistent_80_max_planarity = (
        float(
            np.max(
                planarity_array[
                    :,
                    persistent_80_indices
                ]
            )
        )
        if len(
            persistent_80_indices
        )
        else 0.0
    )

    consecutive_frame_overlaps = [
        len(
            top_sets[index]
            & top_sets[index + 1]
        )
        for index in range(
            len(top_sets) - 1
        )
    ]

    consecutive_block_overlaps = [
        len(
            block_top_sets[index]
            & block_top_sets[index + 1]
        )
        for index in range(
            len(block_top_sets) - 1
        )
    ]

    planarity_q99 = np.array(
        [
            float(
                row[
                    "HBN_planarity_q99_deg"
                ]
            )
            for row in per_frame_rows
        ],
        dtype=float,
    )

    equilibrium_q99 = np.array(
        [
            float(
                row[
                    "HBN_equilibrium_q99_deg"
                ]
            )
            for row in per_frame_rows
        ],
        dtype=float,
    )

    hbn_rms = np.array(
        [
            float(
                row[
                    "HBN_aligned_rms_nm"
                ]
            )
            for row in per_frame_rows
        ],
        dtype=float,
    )

    first_block = block_rows[0]
    last_block = block_rows[-1]

    planarity_block_change = (
        float(
            last_block[
                "HBN_planarity_q99_deg_mean"
            ]
        )
        - float(
            first_block[
                "HBN_planarity_q99_deg_mean"
            ]
        )
    )

    equilibrium_block_change = (
        float(
            last_block[
                "HBN_equilibrium_q99_deg_mean"
            ]
        )
        - float(
            first_block[
                "HBN_equilibrium_q99_deg_mean"
            ]
        )
    )

    hbn_rms_block_change = (
        float(
            last_block[
                "HBN_aligned_rms_nm_mean"
            ]
        )
        - float(
            first_block[
                "HBN_aligned_rms_nm_mean"
            ]
        )
    )

    blocked_reasons = []
    review_reasons = []

    all_numeric_values = []

    for row in per_frame_rows:
        for key, value in row.items():
            if isinstance(
                value,
                (int, float),
            ):
                all_numeric_values.append(
                    float(value)
                )

    if not np.all(
        np.isfinite(
            np.array(
                all_numeric_values,
                dtype=float,
            )
        )
    ):
        blocked_reasons.append(
            "non-finite time-resolved metric"
        )

    severe_checks = {
        "HBN bond deviation exceeds 0.03 nm": (
            max(
                row[
                    "HBN_bond_deviation_max_nm"
                ]
                for row in per_frame_rows
            )
            > 0.03
        ),
        "HBN angle deviation exceeds 25 degrees": (
            max(
                row[
                    "HBN_angle_deviation_max_deg"
                ]
                for row in per_frame_rows
            )
            > 25.0
        ),
        "HBN planarity exceeds 60 degrees": (
            float(
                planarity_array.max()
            )
            > 60.0
        ),
        "calibrated equilibrium deviation exceeds 60 degrees": (
            float(
                equilibrium_array.max()
            )
            > 60.0
        ),
        "PYR bond deviation exceeds 0.03 nm": (
            max(
                row[
                    "PYR_max_bond_deviation_nm"
                ]
                for row in per_frame_rows
            )
            > 0.03
        ),
        "PYR angle deviation exceeds 25 degrees": (
            max(
                row[
                    "PYR_max_angle_deviation_deg"
                ]
                for row in per_frame_rows
            )
            > 25.0
        ),
        "intermolecular contact below 0.10 nm": (
            min(
                min(
                    row[
                        "HBN_PYR_minimum_contact_nm"
                    ],
                    row[
                        "PYR_PYR_minimum_contact_nm"
                    ],
                    row[
                        "PYR_waterO_minimum_contact_nm"
                    ],
                )
                for row in per_frame_rows
            )
            < 0.10
        ),
        "temperature outside 270-330 K": (
            float(
                temperature.min()
            )
            < 270.0
            or float(
                temperature.max()
            )
            > 330.0
        ),
    }

    blocked_reasons.extend(
        label
        for label, failed
        in severe_checks.items()
        if failed
    )

    planarity_slope = linear_slope(
        times,
        planarity_q99,
    )

    equilibrium_slope = linear_slope(
        times,
        equilibrium_q99,
    )

    hbn_rms_slope = linear_slope(
        times,
        hbn_rms,
    )

    if (
        planarity_slope > 0.15
        and planarity_block_change > 2.5
    ):
        review_reasons.append(
            "progressive HBN planarity increase"
        )

    if (
        equilibrium_slope > 0.15
        and equilibrium_block_change > 2.5
    ):
        review_reasons.append(
            "progressive calibrated-equilibrium deviation"
        )

    if (
        hbn_rms_slope > 0.004
        and hbn_rms_block_change > 0.08
    ):
        review_reasons.append(
            "progressive aligned HBN RMS increase"
        )

    if (
        persistent_80 > 0
        and persistent_80_max_planarity
        > 35.0
    ):
        review_reasons.append(
            "persistent localized improper deformation"
        )

    maximum_pyr_rms = max(
        row[
            "PYR_max_aligned_rms_nm"
        ]
        for row in per_frame_rows
    )

    if maximum_pyr_rms > 0.08:
        review_reasons.append(
            "PYR internal aligned RMS exceeds 0.08 nm"
        )

    minimum_contact = min(
        min(
            row[
                "HBN_PYR_minimum_contact_nm"
            ],
            row[
                "PYR_PYR_minimum_contact_nm"
            ],
            row[
                "PYR_waterO_minimum_contact_nm"
            ],
        )
        for row in per_frame_rows
    )

    if minimum_contact < 0.14:
        review_reasons.append(
            "intermolecular contact below 0.14 nm"
        )

    temperature_block_change = abs(
        float(
            last_block[
                "temperature_mean_K"
            ]
        )
        - float(
            first_block[
                "temperature_mean_K"
            ]
        )
    )

    if temperature_block_change > 5.0:
        review_reasons.append(
            "temperature block drift exceeds 5 K"
        )

    if blocked_reasons:
        decision = "BLOCKED"
        next_step = "REVIEW_FORCE_FIELD_OR_STRUCTURE"
    elif review_reasons:
        decision = "REVIEW"
        next_step = "TARGETED_STABILITY_ANALYSIS"
    else:
        decision = "PASS_WITH_MONITORING"
        next_step = (
            "100PS_MOBILE_PRODUCTION_CANDIDATE"
        )

    summary = {
        "stage": STAGE07,
        "frame_count": len(frames),
        "duration_ps": duration,
        "temperature_mean_K": float(
            np.mean(temperature)
        ),
        "temperature_min_K": float(
            np.min(temperature)
        ),
        "temperature_max_K": float(
            np.max(temperature)
        ),
        "potential_first_block_mean_kJ_mol": (
            first_block[
                "potential_mean_kJ_mol"
            ]
        ),
        "potential_last_block_mean_kJ_mol": (
            last_block[
                "potential_mean_kJ_mol"
            ]
        ),
        "HBN_planarity_q99_slope_deg_per_ps": (
            planarity_slope
        ),
        "HBN_planarity_q99_first_last_block_change_deg": (
            planarity_block_change
        ),
        "HBN_equilibrium_q99_slope_deg_per_ps": (
            equilibrium_slope
        ),
        "HBN_equilibrium_q99_first_last_block_change_deg": (
            equilibrium_block_change
        ),
        "HBN_aligned_rms_slope_nm_per_ps": (
            hbn_rms_slope
        ),
        "HBN_aligned_rms_first_last_block_change_nm": (
            hbn_rms_block_change
        ),
        "HBN_max_bond_deviation_nm": max(
            row[
                "HBN_bond_deviation_max_nm"
            ]
            for row in per_frame_rows
        ),
        "HBN_max_angle_deviation_deg": max(
            row[
                "HBN_angle_deviation_max_deg"
            ]
            for row in per_frame_rows
        ),
        "HBN_max_planarity_deg": float(
            planarity_array.max()
        ),
        "HBN_max_equilibrium_deviation_deg": float(
            equilibrium_array.max()
        ),
        "HBN_impropers_above20_persistent50_count": (
            persistent_50
        ),
        "HBN_impropers_above20_persistent80_count": (
            persistent_80
        ),
        "HBN_persistent80_max_planarity_deg": (
            persistent_80_max_planarity
        ),
        "top1pct_consecutive_frame_overlap_mean": (
            float(
                np.mean(
                    consecutive_frame_overlaps
                )
            )
            if consecutive_frame_overlaps
            else 0.0
        ),
        "top1pct_consecutive_block_overlap_mean": (
            float(
                np.mean(
                    consecutive_block_overlaps
                )
            )
            if consecutive_block_overlaps
            else 0.0
        ),
        "PYR_max_aligned_rms_nm": (
            maximum_pyr_rms
        ),
        "minimum_intergroup_contact_nm": (
            minimum_contact
        ),
        "time_resolved_decision": (
            decision
        ),
        "authorized_next_step": (
            next_step
        ),
        "long_mobile_production_authorized": (
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
        PER_FRAME_CSV,
        per_frame_rows,
    )

    write_csv(
        BLOCK_CSV,
        block_rows,
    )

    write_csv(
        SUMMARY_CSV,
        [summary],
    )

    REPORT_MD.write_text(
        f"""# Day022 Stage07 Time-Resolved Stability

- Stage: `{STAGE07}`
- Frames analyzed: {len(frames)}
- Duration: {duration:.3f} ps
- Time-resolved decision: **{decision}**
- Authorized next step: `{next_step}`
- Long mobile production authorized by this script: **NO**

## Thermodynamic behavior

- Temperature mean/min/max: {summary['temperature_mean_K']:.4f}/{summary['temperature_min_K']:.4f}/{summary['temperature_max_K']:.4f} K
- First/last block potential means: {summary['potential_first_block_mean_kJ_mol']:.4f}/{summary['potential_last_block_mean_kJ_mol']:.4f} kJ mol^-1

## HBN trends

- Planarity q99 slope: {planarity_slope:.6f} deg ps^-1
- Planarity q99 first-to-last block change: {planarity_block_change:.4f} deg
- Calibrated-equilibrium q99 slope: {equilibrium_slope:.6f} deg ps^-1
- Calibrated-equilibrium first-to-last block change: {equilibrium_block_change:.4f} deg
- Aligned RMS slope: {hbn_rms_slope:.8f} nm ps^-1
- Aligned RMS first-to-last block change: {hbn_rms_block_change:.6f} nm
- Maximum bond deviation: {summary['HBN_max_bond_deviation_nm']:.6f} nm
- Maximum angle deviation: {summary['HBN_max_angle_deviation_deg']:.4f} deg
- Persistent impropers above 20 deg in >=80% frames: {persistent_80}

## PYR and contacts

- Maximum PYR aligned RMS: {maximum_pyr_rms:.6f} nm
- Minimum intergroup contact: {minimum_contact:.6f} nm
""",
        encoding="utf-8",
    )

    print(
        "Day022 Stage07 time-resolved stability analysis completed."
    )

    print(
        f"Frames / duration: "
        f"{len(frames)} / {duration:.3f} ps"
    )

    print(
        "Temperature mean/min/max: "
        f"{summary['temperature_mean_K']:.4f}/"
        f"{summary['temperature_min_K']:.4f}/"
        f"{summary['temperature_max_K']:.4f} K"
    )

    print(
        "HBN planarity q99 slope / "
        "first-last block change: "
        f"{planarity_slope:.6f} deg/ps / "
        f"{planarity_block_change:.4f} deg"
    )

    print(
        "HBN calibrated-equilibrium q99 slope / "
        "first-last block change: "
        f"{equilibrium_slope:.6f} deg/ps / "
        f"{equilibrium_block_change:.4f} deg"
    )

    print(
        "HBN aligned RMS slope / "
        "first-last block change: "
        f"{hbn_rms_slope:.8f} nm/ps / "
        f"{hbn_rms_block_change:.6f} nm"
    )

    print(
        "Persistent impropers above 20 deg "
        "in >=50% / >=80% frames: "
        f"{persistent_50}/{persistent_80}"
    )

    print(
        "Top 1% overlap mean "
        "consecutive frames / blocks: "
        f"{summary['top1pct_consecutive_frame_overlap_mean']:.3f}/"
        f"{summary['top1pct_consecutive_block_overlap_mean']:.3f}"
    )

    print(
        "Maximum PYR aligned RMS: "
        f"{maximum_pyr_rms:.6f} nm"
    )

    print(
        "Minimum intergroup contact: "
        f"{minimum_contact:.6f} nm"
    )

    print(
        f"Time-resolved decision: {decision}"
    )

    print(
        f"Authorized next step: {next_step}"
    )

    print(
        "Long mobile production authorized: NO"
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
