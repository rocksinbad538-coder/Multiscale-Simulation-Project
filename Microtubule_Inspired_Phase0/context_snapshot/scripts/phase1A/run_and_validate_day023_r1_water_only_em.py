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

PROTOCOL_ROOT = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

PROTOTYPE_ROOT = (
    DAY023_ROOT
    / "02_r1_steric_cap_prototype"
)

STATIC_ROOT = (
    DAY023_ROOT
    / "05_r1_full_static_scan"
)

OUTPUT_ROOT = (
    DAY023_ROOT
    / "06_r1_water_only_em"
)

INPUT_GRO = (
    PROTOTYPE_ROOT
    / "r1_t0_hydrated_with_steric_caps_geometry_only.gro"
)

PROTOTYPE_JSON = (
    PROTOTYPE_ROOT
    / "r1_selected_steric_cap_definition.json"
)

MODEL_JSON = (
    STATIC_ROOT
    / "r1_selected_cap_model.json"
)

INPUT_TOP = (
    STATIC_ROOT
    / "selected/r1_selected_cap_model.top"
)

INPUT_INDEX = (
    STATIC_ROOT
    / "selected/r1_selected_groups.ndx"
)

TEMPLATE_CANDIDATES = (
    PROTOCOL_ROOT
    / "protocol_inputs/mdp/01_em_k10000.mdp",
    PROTOCOL_ROOT
    / "protocol_inputs/mdp/00_em_k100000.mdp",
)

EM_MDP = (
    OUTPUT_ROOT
    / "r1_water_only_em.mdp"
)

LOCAL_GRO = (
    OUTPUT_ROOT
    / "r1_water_only_em_input.gro"
)

LOCAL_TOP = (
    OUTPUT_ROOT
    / "r1_water_only_em.top"
)

LOCAL_INDEX = (
    OUTPUT_ROOT
    / "r1_water_only_em.ndx"
)

TPR = (
    OUTPUT_ROOT
    / "r1_water_only_em.tpr"
)

PROCESSED_TOP = (
    OUTPUT_ROOT
    / "r1_water_only_em_processed.top"
)

MDOUT = (
    OUTPUT_ROOT
    / "r1_water_only_em_mdout.mdp"
)

GROMPP_LOG = (
    OUTPUT_ROOT
    / "r1_water_only_em_grompp.log"
)

MDRUN_CONSOLE = (
    OUTPUT_ROOT
    / "r1_water_only_em_mdrun_console.log"
)

DEFFNM = (
    OUTPUT_ROOT
    / "r1_water_only_em"
)

FINAL_GRO = (
    OUTPUT_ROOT
    / "r1_water_only_em.gro"
)

EDR = (
    OUTPUT_ROOT
    / "r1_water_only_em.edr"
)

MDLOG = (
    OUTPUT_ROOT
    / "r1_water_only_em.log"
)

POTENTIAL_XVG = (
    OUTPUT_ROOT
    / "r1_water_only_em_potential.xvg"
)

CAP_SOL_XVG = (
    OUTPUT_ROOT
    / "r1_water_only_em_cap_sol_lj.xvg"
)

ENERGY_MENU = (
    OUTPUT_ROOT
    / "r1_water_only_em_energy_menu.txt"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r1_water_only_em_summary.csv"
)

GATE_CSV = (
    OUTPUT_ROOT
    / "r1_water_only_em_gates.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R1_WATER_ONLY_EM_VALIDATION_DAY023.md"
)

EXPECTED_ATOMS = 68314
HBN_ATOMS = 1680
PYR_ATOMS = 104
SOLUTE_ATOMS = HBN_ATOMS + PYR_ATOMS
WATER_SITES = 4
EXPECTED_WATERS = 16551
CAP_BEADS_PER_END = 163
TOTAL_CAP_BEADS = 326
INITIAL_LUMEN_WATERS = 428

EMTOL_KJ_MOL_NM = 500.0
EMSTEP_NM = 0.001
EM_NSTEPS = 50000

FROZEN_RMS_TOLERANCE_NM = 1.0e-6
FROZEN_MAX_TOLERANCE_NM = 1.0e-6
MIN_FINAL_CAP_WATER_DISTANCE_NM = 0.17
MIN_FINAL_LUMEN_WATERS = 407
MAX_WATER_O_DISPLACEMENT_NM = 1.0

INSTABILITY_PATTERNS = (
    r"\bnan\b",
    r"fatal error",
    r"segmentation fault",
    r"constraint warning",
    r"shake did not converge",
    r"water molecule cannot be settled",
    r"water molecule can not be settled",
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


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


def select_template() -> Path:
    existing = [
        path
        for path in TEMPLATE_CANDIDATES
        if path.exists() and path.is_file()
    ]

    if not existing:
        raise RuntimeError(
            "No validated Phase 1A EM template was found."
        )

    return existing[0]


def mdp_key(line: str) -> str | None:
    active = line.split(
        ";",
        1,
    )[0].strip()

    if (
        not active
        or active.startswith("#")
        or "=" not in active
    ):
        return None

    return (
        active.split(
            "=",
            1,
        )[0]
        .strip()
        .lower()
        .replace("-", "_")
    )


def build_em_mdp(
    template: Path,
    destination: Path,
) -> None:
    overridden = {
        "integrator",
        "emtol",
        "emstep",
        "nsteps",
        "continuation",
        "define",
        "freezegrps",
        "freezedim",
        "gen_vel",
        "tcoupl",
        "pcoupl",
        "comm_mode",
        "energygrps",
        "nstcalcenergy",
        "nstenergy",
        "nstlog",
        "nstxout",
        "nstvout",
        "nstfout",
        "nstxout_compressed",
    }

    retained_lines = []

    for line in template.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        key = mdp_key(line)

        if key in overridden:
            continue

        retained_lines.append(line)

    retained_lines.extend(
        [
            "",
            "; Day023 R1 water-only minimization overrides",
            "integrator               = steep",
            f"emtol                    = {EMTOL_KJ_MOL_NM:.1f}",
            f"emstep                   = {EMSTEP_NM:.6f}",
            f"nsteps                   = {EM_NSTEPS}",
            "continuation             = no",
            "define                   =",
            "freezegrps               = HBN_PYR CAPS",
            "freezedim                = Y Y Y Y Y Y",
            "gen_vel                  = no",
            "tcoupl                   = no",
            "pcoupl                   = no",
            "comm_mode                = none",
            "energygrps               = HBN_PYR SOL CAPS",
            "nstcalcenergy            = 10",
            "nstenergy                = 10",
            "nstlog                   = 10",
            "nstxout                  = 0",
            "nstvout                  = 0",
            "nstfout                  = 0",
            "nstxout_compressed       = 0",
            "",
        ]
    )

    destination.write_text(
        "\n".join(retained_lines),
        encoding="utf-8",
    )


def parse_box(
    values: list[float],
) -> np.ndarray:
    if len(values) == 3:
        return np.asarray(
            values,
            dtype=float,
        )

    raise RuntimeError(
        "This validation requires an orthorhombic box."
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
    atom_count = int(lines[1].strip())

    if len(lines) < atom_count + 3:
        raise RuntimeError(
            f"Incomplete GRO file: {path}"
        )

    atoms = []

    for index, line in enumerate(
        lines[2 : 2 + atom_count]
    ):
        if len(line) < 44:
            raise RuntimeError(
                f"Malformed atom line {index + 1} in {path}"
            )

        atoms.append(
            {
                "index": index,
                "resid": int(line[0:5]),
                "resname": line[5:10].strip(),
                "atomname": line[10:15].strip(),
                "position": np.asarray(
                    [
                        float(line[20:28]),
                        float(line[28:36]),
                        float(line[36:44]),
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

    return title, atoms, box


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


def displacement_metrics(
    initial: np.ndarray,
    final: np.ndarray,
    box: np.ndarray,
) -> tuple[float, float]:
    difference = minimum_image(
        final - initial,
        box,
    )

    norms = np.linalg.norm(
        difference,
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
        float(np.max(norms)),
    )


def water_oxygen_indices(
    atoms: list[dict[str, Any]],
    water_start: int,
    water_count: int,
    oxygen_name: str,
) -> list[int]:
    indices = []

    for water_index in range(water_count):
        start = (
            water_start
            + water_index * WATER_SITES
        )

        chunk = atoms[
            start:
            start + WATER_SITES
        ]

        candidates = [
            start + local_index
            for local_index, atom in enumerate(chunk)
            if atom["atomname"] == oxygen_name
        ]

        if len(candidates) != 1:
            raise RuntimeError(
                "Could not identify exactly one oxygen "
                f"in water {water_index}: "
                f"{[atom['atomname'] for atom in chunk]}"
            )

        indices.append(candidates[0])

    return indices


def nearest_cap_water_distance(
    water_oxygen_positions: np.ndarray,
    cap_positions: np.ndarray,
    box: np.ndarray,
) -> float:
    minimum = math.inf
    chunk_size = 512

    for start in range(
        0,
        len(water_oxygen_positions),
        chunk_size,
    ):
        stop = min(
            start + chunk_size,
            len(water_oxygen_positions),
        )

        displacement = (
            water_oxygen_positions[
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

        minimum = min(
            minimum,
            float(np.min(distances)),
        )

    return minimum


def lumen_occupancy(
    water_positions: np.ndarray,
    box: np.ndarray,
    prototype: dict[str, Any],
) -> int:
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
        prototype["axial_low_nm"]
    )

    axial_high = float(
        prototype["axial_high_nm"]
    )

    accessible_radius = float(
        prototype["accessible_radius_nm"]
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

    mask = (
        (axial >= axial_low)
        & (axial <= axial_high)
        & (radial <= accessible_radius)
    )

    return int(np.count_nonzero(mask))


def parse_energy_menu(
    output: str,
) -> dict[str, int]:
    terms: dict[str, int] = {}

    pattern = re.compile(
        r"(?<!\S)(\d+)\s+"
        r"([A-Za-z][A-Za-z0-9_.:+\-]*)"
    )

    for line in output.splitlines():
        for match in pattern.finditer(line):
            terms[
                match.group(2)
            ] = int(match.group(1))

    return terms


def probe_energy_menu(
    gmx: str,
    edr: Path,
) -> dict[str, int]:
    probe = run_command(
        [
            gmx,
            "energy",
            "-f",
            str(edr),
            "-o",
            str(
                OUTPUT_ROOT
                / "energy_menu_probe.xvg"
            ),
        ],
        cwd=OUTPUT_ROOT,
        input_text="0\n",
    )

    ENERGY_MENU.write_text(
        probe.stdout,
        encoding="utf-8",
    )

    terms = parse_energy_menu(
        probe.stdout
    )

    if not terms:
        raise RuntimeError(
            "Could not parse the GROMACS energy menu."
        )

    return terms


def resolve_term(
    terms: dict[str, int],
    exact_names: tuple[str, ...],
    required_tokens: tuple[str, ...] = (),
) -> tuple[str, int]:
    for name in exact_names:
        if name in terms:
            return name, terms[name]

    for name, number in terms.items():
        if all(
            token in name
            for token in required_tokens
        ):
            return name, number

    raise RuntimeError(
        "Could not resolve energy term. "
        f"Exact={exact_names}; tokens={required_tokens}"
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

    array = np.asarray(
        rows,
        dtype=float,
    )

    if (
        array.ndim != 2
        or array.shape[1] != 2
        or len(array) == 0
    ):
        raise RuntimeError(
            f"No valid numeric data in {output_path}"
        )

    return array


def parse_final_force(
    log_text: str,
) -> tuple[
    float | None,
    float | None,
    bool,
]:
    maximum_force_matches = re.findall(
        r"Maximum force\s*=\s*"
        r"([-+0-9.eE]+)",
        log_text,
        flags=re.IGNORECASE,
    )

    norm_force_matches = re.findall(
        r"Norm of force\s*=\s*"
        r"([-+0-9.eE]+)",
        log_text,
        flags=re.IGNORECASE,
    )

    maximum_force = (
        float(maximum_force_matches[-1])
        if maximum_force_matches
        else None
    )

    norm_force = (
        float(norm_force_matches[-1])
        if norm_force_matches
        else None
    )

    converged = bool(
        re.search(
            r"Steepest Descents converged",
            log_text,
            flags=re.IGNORECASE,
        )
    )

    return (
        maximum_force,
        norm_force,
        converged,
    )


def write_csv(
    path: Path,
    rows: list[dict[str, Any]],
) -> None:
    if not rows:
        raise RuntimeError(
            f"No rows available for {path}"
        )

    fields = []

    for row in rows:
        for field in row:
            if field not in fields:
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
                    field: row.get(field, "")
                    for field in fields
                }
            )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        INPUT_GRO,
        PROTOTYPE_JSON,
        MODEL_JSON,
        INPUT_TOP,
        INPUT_INDEX,
    ):
        require_file(required)

    prototype = json.loads(
        PROTOTYPE_JSON.read_text(
            encoding="utf-8"
        )
    )

    model = json.loads(
        MODEL_JSON.read_text(
            encoding="utf-8"
        )
    )

    if (
        model.get("decision")
        != "FULL_R1_STATIC_CAP_WATER_MODEL_VALIDATED"
    ):
        raise RuntimeError(
            "The selected R1 cap model is not validated."
        )

    if not bool(
        model.get(
            "energy_minimization_authorized"
        )
    ):
        raise RuntimeError(
            "Energy minimization is not authorized."
        )

    if bool(
        model.get(
            "MD_execution_authorized"
        )
    ):
        raise RuntimeError(
            "Unexpected MD authorization state."
        )

    retained_waters = int(
        model["retained_water_molecules"]
    )

    cap_beads_per_end = int(
        model["cap_beads_per_end"]
    )

    if retained_waters != EXPECTED_WATERS:
        raise RuntimeError(
            "Unexpected retained-water count."
        )

    if cap_beads_per_end != CAP_BEADS_PER_END:
        raise RuntimeError(
            "Unexpected cap size."
        )

    template = select_template()
    gmx = locate_gmx()

    shutil.copy2(
        INPUT_GRO,
        LOCAL_GRO,
    )

    shutil.copy2(
        INPUT_TOP,
        LOCAL_TOP,
    )

    shutil.copy2(
        INPUT_INDEX,
        LOCAL_INDEX,
    )

    build_em_mdp(
        template,
        EM_MDP,
    )

    grompp = run_command(
        [
            gmx,
            "grompp",
            "-f",
            str(EM_MDP),
            "-c",
            str(LOCAL_GRO),
            "-r",
            str(LOCAL_GRO),
            "-p",
            str(LOCAL_TOP),
            "-n",
            str(LOCAL_INDEX),
            "-o",
            str(TPR),
            "-po",
            str(MDOUT),
            "-pp",
            str(PROCESSED_TOP),
            "-maxwarn",
            "0",
        ],
        cwd=OUTPUT_ROOT,
    )

    GROMPP_LOG.write_text(
        grompp.stdout,
        encoding="utf-8",
    )

    if (
        grompp.returncode != 0
        or not TPR.exists()
    ):
        raise RuntimeError(
            "R1 water-only EM grompp failed. "
            f"See {GROMPP_LOG}"
        )

    mdrun = run_command(
        [
            gmx,
            "mdrun",
            "-s",
            str(TPR),
            "-deffnm",
            str(DEFFNM),
            "-ntmpi",
            "1",
            "-ntomp",
            "4",
        ],
        cwd=OUTPUT_ROOT,
    )

    MDRUN_CONSOLE.write_text(
        mdrun.stdout,
        encoding="utf-8",
    )

    for required in (
        FINAL_GRO,
        EDR,
        MDLOG,
    ):
        require_file(required)

    initial_title, initial_atoms, initial_box = (
        read_gro(LOCAL_GRO)
    )

    final_title, final_atoms, final_box = (
        read_gro(FINAL_GRO)
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

    box_difference = float(
        np.max(
            np.abs(
                final_box - initial_box
            )
        )
    )

    initial_positions = positions(
        initial_atoms
    )

    final_positions = positions(
        final_atoms
    )

    water_start = SOLUTE_ATOMS
    water_stop = (
        water_start
        + retained_waters
        * WATER_SITES
    )

    cap_start = water_stop
    cap_stop = (
        cap_start
        + TOTAL_CAP_BEADS
    )

    if cap_stop != EXPECTED_ATOMS:
        raise RuntimeError(
            "R1 segment accounting failed."
        )

    oxygen_name = "OW"

    water_oxygen_atom_indices = (
        water_oxygen_indices(
            initial_atoms,
            water_start,
            retained_waters,
            oxygen_name,
        )
    )

    initial_water_oxygen_positions = (
        initial_positions[
            water_oxygen_atom_indices
        ]
    )

    final_water_oxygen_positions = (
        final_positions[
            water_oxygen_atom_indices
        ]
    )

    initial_cap_positions = (
        initial_positions[
            cap_start:
            cap_stop
        ]
    )

    final_cap_positions = (
        final_positions[
            cap_start:
            cap_stop
        ]
    )

    hbn_rms, hbn_max = displacement_metrics(
        initial_positions[:HBN_ATOMS],
        final_positions[:HBN_ATOMS],
        initial_box,
    )

    pyr_rms, pyr_max = displacement_metrics(
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

    cap_rms, cap_max = displacement_metrics(
        initial_cap_positions,
        final_cap_positions,
        initial_box,
    )

    water_o_rms, water_o_max = (
        displacement_metrics(
            initial_water_oxygen_positions,
            final_water_oxygen_positions,
            initial_box,
        )
    )

    initial_cap_water_distance = (
        nearest_cap_water_distance(
            initial_water_oxygen_positions,
            initial_cap_positions,
            initial_box,
        )
    )

    final_cap_water_distance = (
        nearest_cap_water_distance(
            final_water_oxygen_positions,
            final_cap_positions,
            final_box,
        )
    )

    initial_lumen_occupancy = lumen_occupancy(
        initial_water_oxygen_positions,
        initial_box,
        prototype,
    )

    final_lumen_occupancy = lumen_occupancy(
        final_water_oxygen_positions,
        final_box,
        prototype,
    )

    terms = probe_energy_menu(
        gmx,
        EDR,
    )

    potential_name, potential_number = (
        resolve_term(
            terms,
            ("Potential",),
            ("Potential",),
        )
    )

    cap_sol_name, cap_sol_number = (
        resolve_term(
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

    potential_data = extract_energy_series(
        gmx,
        EDR,
        potential_number,
        POTENTIAL_XVG,
    )

    cap_sol_data = extract_energy_series(
        gmx,
        EDR,
        cap_sol_number,
        CAP_SOL_XVG,
    )

    potential_initial = float(
        potential_data[0, 1]
    )

    potential_final = float(
        potential_data[-1, 1]
    )

    potential_change = (
        potential_final
        - potential_initial
    )

    cap_sol_initial = float(
        cap_sol_data[0, 1]
    )

    cap_sol_final = float(
        cap_sol_data[-1, 1]
    )

    mdlog_text = MDLOG.read_text(
        encoding="utf-8",
        errors="replace",
    )

    console_text = MDRUN_CONSOLE.read_text(
        encoding="utf-8",
        errors="replace",
    )

    combined_log = (
        mdlog_text
        + "\n"
        + console_text
    )

    maximum_force, norm_force, converged = (
        parse_final_force(
            combined_log
        )
    )

    instability_hits = []

    for pattern in INSTABILITY_PATTERNS:
        if re.search(
            pattern,
            combined_log,
            flags=re.IGNORECASE,
        ):
            instability_hits.append(pattern)

    finished = bool(
        re.search(
            r"Finished mdrun",
            combined_log,
            flags=re.IGNORECASE,
        )
    )

    force_gate = bool(
        converged
        or (
            maximum_force is not None
            and maximum_force
            <= EMTOL_KJ_MOL_NM * 1.01
        )
    )

    gates = {
        "grompp_return_code_zero": (
            grompp.returncode == 0
        ),
        "mdrun_return_code_zero": (
            mdrun.returncode == 0
        ),
        "mdrun_finished": finished,
        "no_instability_signatures": (
            len(instability_hits) == 0
        ),
        "final_atom_count_is_68314": (
            len(final_atoms)
            == EXPECTED_ATOMS
        ),
        "box_is_unchanged": (
            box_difference <= 1.0e-6
        ),
        "HBN_is_frozen": (
            hbn_rms
            <= FROZEN_RMS_TOLERANCE_NM
            and hbn_max
            <= FROZEN_MAX_TOLERANCE_NM
        ),
        "PYR_is_frozen": (
            pyr_rms
            <= FROZEN_RMS_TOLERANCE_NM
            and pyr_max
            <= FROZEN_MAX_TOLERANCE_NM
        ),
        "caps_are_frozen": (
            cap_rms
            <= FROZEN_RMS_TOLERANCE_NM
            and cap_max
            <= FROZEN_MAX_TOLERANCE_NM
        ),
        "energy_is_finite": (
            math.isfinite(potential_initial)
            and math.isfinite(potential_final)
        ),
        "potential_did_not_increase": (
            potential_final
            <= potential_initial + 1.0e-3
        ),
        "force_convergence": force_gate,
        "final_cap_water_distance": (
            final_cap_water_distance
            >= MIN_FINAL_CAP_WATER_DISTANCE_NM
        ),
        "lumen_water_retention": (
            final_lumen_occupancy
            >= MIN_FINAL_LUMEN_WATERS
        ),
        "water_displacement_is_finite": (
            math.isfinite(water_o_rms)
            and math.isfinite(water_o_max)
        ),
        "water_displacement_is_local": (
            water_o_max
            <= MAX_WATER_O_DISPLACEMENT_NM
        ),
        "cap_water_energy_did_not_increase": (
            cap_sol_final
            <= cap_sol_initial + 1.0e-3
        ),
    }

    failed_gates = [
        name
        for name, passed in gates.items()
        if not passed
    ]

    accepted = (
        len(failed_gates) == 0
    )

    decision = (
        "R1_WATER_ONLY_EM_VALIDATED"
        if accepted
        else
        "R1_WATER_ONLY_EM_REQUIRES_REVIEW"
    )

    required_next_step = (
        "PREPARE_R1_FROZEN_SOLUTE_SHORT_NVT_SCREENING"
        if accepted
        else
        "REVIEW_R1_WATER_ONLY_EM_GATE_FAILURES"
    )

    summary = {
        "decision": decision,
        "template_mdp": relative(template),
        "grompp_return_code": grompp.returncode,
        "mdrun_return_code": mdrun.returncode,
        "mdrun_finished": finished,
        "em_converged": converged,
        "emtol_kJ_mol_nm": EMTOL_KJ_MOL_NM,
        "final_maximum_force_kJ_mol_nm": (
            maximum_force
            if maximum_force is not None
            else ""
        ),
        "final_norm_force_kJ_mol_nm": (
            norm_force
            if norm_force is not None
            else ""
        ),
        "potential_term": potential_name,
        "potential_initial_kJ_mol": (
            potential_initial
        ),
        "potential_final_kJ_mol": (
            potential_final
        ),
        "potential_change_kJ_mol": (
            potential_change
        ),
        "cap_SOL_LJ_term": cap_sol_name,
        "cap_SOL_LJ_initial_kJ_mol": (
            cap_sol_initial
        ),
        "cap_SOL_LJ_final_kJ_mol": (
            cap_sol_final
        ),
        "HBN_RMS_displacement_nm": hbn_rms,
        "HBN_max_displacement_nm": hbn_max,
        "PYR_RMS_displacement_nm": pyr_rms,
        "PYR_max_displacement_nm": pyr_max,
        "caps_RMS_displacement_nm": cap_rms,
        "caps_max_displacement_nm": cap_max,
        "waterO_RMS_displacement_nm": water_o_rms,
        "waterO_max_displacement_nm": water_o_max,
        "initial_cap_water_distance_nm": (
            initial_cap_water_distance
        ),
        "final_cap_water_distance_nm": (
            final_cap_water_distance
        ),
        "initial_lumen_water_count": (
            initial_lumen_occupancy
        ),
        "final_lumen_water_count": (
            final_lumen_occupancy
        ),
        "lumen_water_change": (
            final_lumen_occupancy
            - initial_lumen_occupancy
        ),
        "box_max_difference_nm": box_difference,
        "instability_signature_count": (
            len(instability_hits)
        ),
        "instability_signatures": (
            " | ".join(instability_hits)
        ),
        "failed_gates": (
            " | ".join(failed_gates)
        ),
        "short_NVT_preparation_authorized": (
            accepted
        ),
        "MD_execution_authorized": False,
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
        f"""# R1 Water-Only Energy Minimization

## Scope

The validated R1 positive-control geometry was minimized with:

- HBN frozen;
- all pyrene atoms frozen;
- both steric caps frozen;
- only the {retained_waters} TIP4P/2005 water molecules mobile.

No molecular dynamics was performed.

## Inputs

- Coordinates:
  `{relative(LOCAL_GRO)}`
- Topology:
  `{relative(LOCAL_TOP)}`
- Index:
  `{relative(LOCAL_INDEX)}`
- EM template:
  `{relative(template)}`
- EM tolerance:
  **{EMTOL_KJ_MOL_NM:.1f} kJ mol^-1 nm^-1**
- Maximum steps:
  **{EM_NSTEPS}**

## Energy minimization result

- Grompp return code:
  **{grompp.returncode}**
- Mdrun return code:
  **{mdrun.returncode}**
- Finished mdrun:
  **{'YES' if finished else 'NO'}**
- Explicit convergence message:
  **{'YES' if converged else 'NO'}**
- Final maximum force:
  **{maximum_force if maximum_force is not None else 'UNAVAILABLE'}
  kJ mol^-1 nm^-1**
- Final force norm:
  **{norm_force if norm_force is not None else 'UNAVAILABLE'}
  kJ mol^-1 nm^-1**

## Energy

- Initial potential:
  **{potential_initial:.6f} kJ/mol**
- Final potential:
  **{potential_final:.6f} kJ/mol**
- Potential change:
  **{potential_change:.6f} kJ/mol**
- Initial CAP–SOL LJ:
  **{cap_sol_initial:.6f} kJ/mol**
- Final CAP–SOL LJ:
  **{cap_sol_final:.6f} kJ/mol**

## Frozen-group integrity

- HBN RMS/max displacement:
  **{hbn_rms:.12f}/{hbn_max:.12f} nm**
- PYR RMS/max displacement:
  **{pyr_rms:.12f}/{pyr_max:.12f} nm**
- Cap RMS/max displacement:
  **{cap_rms:.12f}/{cap_max:.12f} nm**

## Water relaxation

- Water-O RMS displacement:
  **{water_o_rms:.6f} nm**
- Water-O maximum displacement:
  **{water_o_max:.6f} nm**
- Initial/final minimum CAP–OW distance:
  **{initial_cap_water_distance:.6f}/
  {final_cap_water_distance:.6f} nm**
- Initial/final lumen occupancy:
  **{initial_lumen_occupancy}/
  {final_lumen_occupancy} waters**
- Lumen occupancy change:
  **{final_lumen_occupancy - initial_lumen_occupancy} waters**

## Validation gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Short frozen-solute NVT preparation authorized:
  **{'YES' if accepted else 'NO'}**
- Molecular dynamics execution authorized:
  **NO**
- Required next step:
  `{required_next_step}`

The next stage must prepare and statically validate a short
frozen-solute NVT screening before any dynamics are executed.
""",
        encoding="utf-8",
    )

    print(
        "Day023 R1 water-only energy minimization "
        "and validation completed."
    )

    print(
        "Grompp / mdrun return codes: "
        f"{grompp.returncode}/{mdrun.returncode}"
    )

    print(
        "Finished / explicit convergence: "
        f"{'YES' if finished else 'NO'} / "
        f"{'YES' if converged else 'NO'}"
    )

    print(
        "Final maximum force / norm: "
        f"{maximum_force if maximum_force is not None else 'UNAVAILABLE'} / "
        f"{norm_force if norm_force is not None else 'UNAVAILABLE'} "
        "kJ mol^-1 nm^-1"
    )

    print(
        "Potential initial / final / change: "
        f"{potential_initial:.6f}/"
        f"{potential_final:.6f}/"
        f"{potential_change:.6f} kJ/mol"
    )

    print(
        "CAP-SOL LJ initial / final: "
        f"{cap_sol_initial:.6f}/"
        f"{cap_sol_final:.6f} kJ/mol"
    )

    print(
        "HBN RMS/max displacement: "
        f"{hbn_rms:.12f}/"
        f"{hbn_max:.12f} nm"
    )

    print(
        "PYR RMS/max displacement: "
        f"{pyr_rms:.12f}/"
        f"{pyr_max:.12f} nm"
    )

    print(
        "Caps RMS/max displacement: "
        f"{cap_rms:.12f}/"
        f"{cap_max:.12f} nm"
    )

    print(
        "Water-O RMS/max displacement: "
        f"{water_o_rms:.6f}/"
        f"{water_o_max:.6f} nm"
    )

    print(
        "Minimum CAP-OW initial/final: "
        f"{initial_cap_water_distance:.6f}/"
        f"{final_cap_water_distance:.6f} nm"
    )

    print(
        "Lumen occupancy initial/final/change: "
        f"{initial_lumen_occupancy}/"
        f"{final_lumen_occupancy}/"
        f"{final_lumen_occupancy - initial_lumen_occupancy}"
    )

    print(
        "Instability signatures: "
        + (
            "NONE"
            if not instability_hits
            else " | ".join(instability_hits)
        )
    )

    print(
        f"Decision: {decision}"
    )

    print(
        "Failed gates: "
        + (
            "NONE"
            if not failed_gates
            else " | ".join(failed_gates)
        )
    )

    print(
        "Short frozen-solute NVT preparation authorized: "
        f"{'YES' if accepted else 'NO'}"
    )

    print(
        "MD execution authorized: NO"
    )

    print(
        f"Required next step: {required_next_step}"
    )

    print(
        f"Wrote: {relative(EM_MDP)}"
    )

    print(
        f"Wrote: {relative(TPR)}"
    )

    print(
        f"Wrote: {relative(FINAL_GRO)}"
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
            "R1 water-only minimization requires review."
        )


if __name__ == "__main__":
    main()
