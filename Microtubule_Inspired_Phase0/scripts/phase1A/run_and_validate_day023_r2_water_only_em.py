#!/usr/bin/env python3

from __future__ import annotations

import csv
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

R1_PROTOTYPE_ROOT = (
    DAY023_ROOT
    / "02_r1_steric_cap_prototype"
)

R2_GEOMETRY_ROOT = (
    DAY023_ROOT
    / "12_r2_partial_cap_geometry_design"
)

R2_STATIC_ROOT = (
    DAY023_ROOT
    / "13_r2_topology_static_scan"
)

R2_STATIC_SELECTED_ROOT = (
    R2_STATIC_ROOT
    / "selected"
)

R2_INSTABILITY_AUDIT_ROOT = (
    R2_STATIC_ROOT
    / "postrun_instability_audit"
)

OUTPUT_ROOT = (
    DAY023_ROOT
    / "14_r2_water_only_em"
)

SELECTED_ROOT = (
    OUTPUT_ROOT
    / "selected"
)

SOURCE_GRO = (
    R2_GEOMETRY_ROOT
    / "selected"
    / "r2_selected_partial_cap_geometry_only.gro"
)

SOURCE_TOPOLOGY = (
    R2_STATIC_SELECTED_ROOT
    / "r2_selected_cap_model.top"
)

SOURCE_INDEX = (
    R2_STATIC_SELECTED_ROOT
    / "r2_static_energy_groups.ndx"
)

R1_PROTOTYPE_JSON = (
    R1_PROTOTYPE_ROOT
    / "r1_selected_steric_cap_definition.json"
)

R2_STATIC_SUMMARY = (
    R2_STATIC_ROOT
    / "r2_topology_static_scan_summary.csv"
)

R2_INSTABILITY_SUMMARY = (
    R2_INSTABILITY_AUDIT_ROOT
    / "r2_static_instability_audit_summary.csv"
)

LOCAL_GRO = (
    SELECTED_ROOT
    / "r2_water_only_em_input.gro"
)

LOCAL_TOPOLOGY = (
    SELECTED_ROOT
    / "r2_water_only_em.top"
)

LOCAL_INDEX = (
    SELECTED_ROOT
    / "r2_water_only_em.ndx"
)

EM_MDP = (
    SELECTED_ROOT
    / "r2_water_only_em.mdp"
)

EM_TPR = (
    SELECTED_ROOT
    / "r2_water_only_em.tpr"
)

GROMPP_LOG = (
    OUTPUT_ROOT
    / "r2_water_only_em_grompp.log"
)

MDRUN_CONSOLE = (
    OUTPUT_ROOT
    / "r2_water_only_em_mdrun_console.log"
)

DEFFNM = (
    OUTPUT_ROOT
    / "r2_water_only_em"
)

ENERGY_MENU_LOG = (
    OUTPUT_ROOT
    / "r2_water_only_em_energy_menu.txt"
)

POTENTIAL_XVG = (
    OUTPUT_ROOT
    / "r2_water_only_em_potential.xvg"
)

CAP_SOL_LJ_XVG = (
    OUTPUT_ROOT
    / "r2_water_only_em_cap_sol_lj.xvg"
)

CAP_SOL_COUL_XVG = (
    OUTPUT_ROOT
    / "r2_water_only_em_cap_sol_coulomb.xvg"
)

ENERGY_SERIES_CSV = (
    OUTPUT_ROOT
    / "r2_water_only_em_energy_series.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_water_only_em_summary.csv"
)

GATE_CSV = (
    OUTPUT_ROOT
    / "r2_water_only_em_gates.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_WATER_ONLY_ENERGY_MINIMIZATION_DAY023.md"
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

EXPECTED_INITIAL_LUMEN_WATERS = 428

EMTOL_KJ_MOL_NM = 500.0
EMSTEP_NM = 0.001
MAX_EM_STEPS = 5000

FROZEN_COORDINATE_TOLERANCE_NM = 1.0e-6
MAX_WATER_O_DISPLACEMENT_NM = 0.10

MIN_LUMEN_RETENTION_FRACTION = 0.98
MIN_CAP_OW_DISTANCE_NM = 0.15

MAX_CAP_SOL_LJ_KJ_MOL = 100.0
ZERO_ENERGY_TOLERANCE_KJ_MOL = 1.0e-4

EXPECTED_STATIC_DECISION = (
    "R2_TOPOLOGY_AND_STATIC_CAP_WATER_MODEL_VALIDATED"
)

EXPECTED_AUDIT_DECISION = (
    "R2_STATIC_INSTABILITY_SCAN_CONFIRMED_CLEAN"
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
            "Could not capture process output."
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
            "An orthorhombic box is required."
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
        WATERS
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
                f"for water {water_index}."
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

    axis /= np.linalg.norm(
        axis
    )

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
    if (
        len(first) == 0
        or len(second) == 0
    ):
        raise RuntimeError(
            "Cannot calculate a distance for an empty group."
        )

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
) -> tuple[
    str,
    int,
]:
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

        if len(fields) < 2:
            continue

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
            f"Invalid energy series in {output_path}"
        )

    return data


def parse_em_metrics(
    text: str,
) -> dict[str, Any]:
    converged_match = re.search(
        r"Steepest\s+Descents\s+converged\s+to\s+"
        r"Fmax\s*<\s*([-+0-9.eE]+)\s+in\s+"
        r"(\d+)\s+steps",
        text,
        flags=re.IGNORECASE,
    )

    maximum_force_match = re.search(
        r"Maximum\s+force\s*=\s*"
        r"([-+0-9.eE]+)"
        r"(?:\s+on\s+atom\s+(\d+))?",
        text,
        flags=re.IGNORECASE,
    )

    norm_force_match = re.search(
        r"Norm\s+of\s+force\s*=\s*"
        r"([-+0-9.eE]+)",
        text,
        flags=re.IGNORECASE,
    )

    potential_matches = re.findall(
        r"Potential\s+Energy\s*=\s*"
        r"([-+0-9.eE]+)",
        text,
        flags=re.IGNORECASE,
    )

    return {
        "converged": (
            converged_match is not None
            and re.search(
                r"did\s+not\s+converge",
                text,
                flags=re.IGNORECASE,
            )
            is None
        ),
        "convergence_target_kJ_mol_nm": (
            float(
                converged_match.group(1)
            )
            if converged_match
            else math.nan
        ),
        "steps": (
            int(
                converged_match.group(2)
            )
            if converged_match
            else -1
        ),
        "maximum_force_kJ_mol_nm": (
            float(
                maximum_force_match.group(1)
            )
            if maximum_force_match
            else math.nan
        ),
        "maximum_force_atom": (
            int(
                maximum_force_match.group(2)
            )
            if (
                maximum_force_match
                and maximum_force_match.group(2)
            )
            else -1
        ),
        "norm_force_kJ_mol_nm": (
            float(
                norm_force_match.group(1)
            )
            if norm_force_match
            else math.nan
        ),
        "logged_potential_initial_kJ_mol": (
            float(
                potential_matches[0]
            )
            if potential_matches
            else math.nan
        ),
        "logged_potential_final_kJ_mol": (
            float(
                potential_matches[-1]
            )
            if potential_matches
            else math.nan
        ),
    }


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    SELECTED_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        SOURCE_GRO,
        SOURCE_TOPOLOGY,
        SOURCE_INDEX,
        R1_PROTOTYPE_JSON,
        R2_STATIC_SUMMARY,
        R2_INSTABILITY_SUMMARY,
    ):
        require_file(required)

    static_summary = read_single_csv_row(
        R2_STATIC_SUMMARY
    )

    instability_summary = read_single_csv_row(
        R2_INSTABILITY_SUMMARY
    )

    if (
        static_summary.get(
            "decision"
        )
        != EXPECTED_STATIC_DECISION
    ):
        raise RuntimeError(
            "R2 static topology gate is not accepted."
        )

    if not parse_bool(
        static_summary.get(
            "water_only_energy_minimization_authorized",
            "false",
        )
    ):
        raise RuntimeError(
            "The static topology gate did not authorize EM."
        )

    if (
        instability_summary.get(
            "decision"
        )
        != EXPECTED_AUDIT_DECISION
    ):
        raise RuntimeError(
            "The strict instability audit is not accepted."
        )

    if not parse_bool(
        instability_summary.get(
            "water_only_energy_minimization_authorized",
            "false",
        )
    ):
        raise RuntimeError(
            "The strict instability audit did not authorize EM."
        )

    existing_products = [
        path
        for path in OUTPUT_ROOT.glob(
            f"{DEFFNM.name}.*"
        )
        if path.stat().st_size > 0
    ]

    if existing_products:
        raise RuntimeError(
            "R2 water-only EM products already exist: "
            + ", ".join(
                path.name
                for path in existing_products
            )
        )

    shutil.copy2(
        SOURCE_GRO,
        LOCAL_GRO,
    )

    shutil.copy2(
        SOURCE_TOPOLOGY,
        LOCAL_TOPOLOGY,
    )

    shutil.copy2(
        SOURCE_INDEX,
        LOCAL_INDEX,
    )

    EM_MDP.write_text(
        f"""integrator               = steep
emtol                    = {EMTOL_KJ_MOL_NM:.1f}
emstep                   = {EMSTEP_NM:.6f}
nsteps                   = {MAX_EM_STEPS}

cutoff-scheme            = Verlet
nstlist                  = 20
rlist                    = 1.20

coulombtype              = PME
rcoulomb                 = 1.20
fourierspacing           = 0.12
pme-order                = 4
ewald-rtol               = 1e-5

vdwtype                  = Cut-off
rvdw                     = 1.20
vdw-modifier             = Potential-shift
DispCorr                 = EnerPres

constraints              = h-bonds
constraint-algorithm     = lincs
lincs-iter               = 1
lincs-order              = 4

pbc                      = xyz

freezegrps               = HBN_PYR CAPS
freezedim                = Y Y Y Y Y Y

nstlog                   = 1
nstenergy                = 1
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstxout-compressed       = 0

energygrps               = HBN_PYR CAPS SOL
""",
        encoding="utf-8",
    )

    gmx = locate_gmx()

    grompp = run_command(
        [
            gmx,
            "grompp",
            "-f",
            str(EM_MDP),
            "-c",
            str(LOCAL_GRO),
            "-p",
            str(LOCAL_TOPOLOGY),
            "-n",
            str(LOCAL_INDEX),
            "-o",
            str(EM_TPR),
            "-po",
            str(
                OUTPUT_ROOT
                / "r2_water_only_em_processed.mdp"
            ),
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
        or not EM_TPR.exists()
        or EM_TPR.stat().st_size == 0
    ):
        raise RuntimeError(
            "R2 water-only EM grompp failed. "
            f"See {GROMPP_LOG}"
        )

    warning_count = len(
        re.findall(
            r"^\s*WARNING\s+\d+",
            grompp.stdout,
            flags=re.MULTILINE,
        )
    )

    print(
        "Starting R2 water-only energy minimization."
    )

    print(
        "HBN, all PYR molecules, and both partial caps "
        "are frozen in all dimensions."
    )

    print(
        "Only the 16565 TIP4P/2005 water molecules "
        "are allowed to minimize."
    )

    mdrun_return_code = run_live(
        [
            gmx,
            "mdrun",
            "-s",
            str(EM_TPR),
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

    final_gro = Path(
        str(DEFFNM)
        + ".gro"
    )

    edr = Path(
        str(DEFFNM)
        + ".edr"
    )

    log = Path(
        str(DEFFNM)
        + ".log"
    )

    for product in (
        final_gro,
        edr,
        log,
    ):
        require_file(product)

    mdrun_text = (
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

    em_metrics = parse_em_metrics(
        mdrun_text
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
            mdrun_text
        )
    ]

    prototype = __import__(
        "json"
    ).loads(
        R1_PROTOTYPE_JSON.read_text(
            encoding="utf-8"
        )
    )

    (
        initial_title,
        initial_atoms,
        initial_box,
    ) = read_gro(
        LOCAL_GRO
    )

    (
        final_title,
        final_atoms,
        final_box,
    ) = read_gro(
        final_gro
    )

    if len(initial_atoms) != EXPECTED_ATOMS:
        raise RuntimeError(
            "Unexpected initial atom count."
        )

    if len(final_atoms) != EXPECTED_ATOMS:
        raise RuntimeError(
            "Unexpected final atom count."
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

    final_lumen_mask = lumen_mask(
        final_water_positions,
        final_box,
        prototype,
    )

    initial_lumen_count = int(
        np.count_nonzero(
            initial_lumen_mask
        )
    )

    final_lumen_count = int(
        np.count_nonzero(
            final_lumen_mask
        )
    )

    retained_initial_lumen_count = int(
        np.count_nonzero(
            initial_lumen_mask
            & final_lumen_mask
        )
    )

    minimum_retained_lumen_waters = int(
        math.ceil(
            MIN_LUMEN_RETENTION_FRACTION
            * initial_lumen_count
        )
    )

    hbn_rms, hbn_max = displacement_metrics(
        initial_positions[
            :HBN_ATOMS
        ],
        final_positions[
            :HBN_ATOMS
        ],
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

    caps_rms, caps_max = displacement_metrics(
        initial_cap_positions,
        final_cap_positions,
        initial_box,
    )

    water_rms, water_max = displacement_metrics(
        initial_water_positions,
        final_water_positions,
        initial_box,
    )

    initial_cap_ow_distance = minimum_pair_distance(
        initial_cap_positions,
        initial_water_positions,
        initial_box,
    )

    final_cap_ow_distance = minimum_pair_distance(
        final_cap_positions,
        final_water_positions,
        final_box,
    )

    box_max_difference = float(
        np.max(
            np.abs(
                final_box
                - initial_box
            )
        )
    )

    menu_probe = run_command(
        [
            gmx,
            "energy",
            "-f",
            str(edr),
            "-o",
            str(
                OUTPUT_ROOT
                / "r2_water_only_em_energy_menu_probe.xvg"
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

    potential_name, potential_number = resolve_energy_term(
        terms,
        (
            "Potential",
        ),
        (
            "Potential",
        ),
    )

    cap_sol_lj_name, cap_sol_lj_number = resolve_energy_term(
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

    cap_sol_coul_name, cap_sol_coul_number = resolve_energy_term(
        terms,
        (
            "Coul-SR:CAPS-SOL",
            "Coul-SR:SOL-CAPS",
        ),
        (
            "Coul-SR",
            "CAPS",
            "SOL",
        ),
    )

    potential_data = extract_energy_series(
        gmx,
        edr,
        potential_number,
        POTENTIAL_XVG,
    )

    cap_sol_lj_data = extract_energy_series(
        gmx,
        edr,
        cap_sol_lj_number,
        CAP_SOL_LJ_XVG,
    )

    cap_sol_coul_data = extract_energy_series(
        gmx,
        edr,
        cap_sol_coul_number,
        CAP_SOL_COUL_XVG,
    )

    if not (
        len(potential_data)
        == len(cap_sol_lj_data)
        == len(cap_sol_coul_data)
    ):
        raise RuntimeError(
            "Energy-series lengths do not match."
        )

    if not (
        np.allclose(
            potential_data[:, 0],
            cap_sol_lj_data[:, 0],
            rtol=0.0,
            atol=1.0e-10,
        )
        and np.allclose(
            potential_data[:, 0],
            cap_sol_coul_data[:, 0],
            rtol=0.0,
            atol=1.0e-10,
        )
    ):
        raise RuntimeError(
            "Energy time/step grids do not match."
        )

    energy_rows = []

    for index in range(
        len(potential_data)
    ):
        energy_rows.append(
            {
                "energy_record": index,
                "step_or_time": float(
                    potential_data[
                        index,
                        0,
                    ]
                ),
                "potential_kJ_mol": float(
                    potential_data[
                        index,
                        1,
                    ]
                ),
                "CAP_SOL_LJ_kJ_mol": float(
                    cap_sol_lj_data[
                        index,
                        1,
                    ]
                ),
                "CAP_SOL_Coulomb_kJ_mol": float(
                    cap_sol_coul_data[
                        index,
                        1,
                    ]
                ),
            }
        )

    write_csv(
        ENERGY_SERIES_CSV,
        energy_rows,
    )

    potentials = potential_data[
        :,
        1
    ]

    cap_sol_lj = cap_sol_lj_data[
        :,
        1
    ]

    cap_sol_coul = cap_sol_coul_data[
        :,
        1
    ]

    maximum_force = float(
        em_metrics[
            "maximum_force_kJ_mol_nm"
        ]
    )

    norm_force = float(
        em_metrics[
            "norm_force_kJ_mol_nm"
        ]
    )

    gates = {
        "R2_static_gate_is_validated": (
            static_summary.get(
                "decision"
            )
            == EXPECTED_STATIC_DECISION
        ),
        "R2_strict_instability_audit_is_clean": (
            instability_summary.get(
                "decision"
            )
            == EXPECTED_AUDIT_DECISION
        ),
        "grompp_return_code_zero": (
            grompp.returncode == 0
        ),
        "grompp_warning_count_zero": (
            warning_count == 0
        ),
        "mdrun_return_code_zero": (
            mdrun_return_code == 0
        ),
        "no_instability_signatures": (
            len(
                instability_hits
            )
            == 0
        ),
        "steepest_descents_converged": (
            bool(
                em_metrics[
                    "converged"
                ]
            )
        ),
        "maximum_force_is_finite": (
            math.isfinite(
                maximum_force
            )
        ),
        "maximum_force_is_at_most_emtol": (
            math.isfinite(
                maximum_force
            )
            and maximum_force
            <= EMTOL_KJ_MOL_NM
            * 1.001
        ),
        "norm_force_is_finite": (
            math.isfinite(
                norm_force
            )
        ),
        "initial_atom_count_is_68332": (
            len(initial_atoms)
            == EXPECTED_ATOMS
        ),
        "final_atom_count_is_68332": (
            len(final_atoms)
            == EXPECTED_ATOMS
        ),
        "box_is_unchanged": (
            box_max_difference
            <= 1.0e-6
        ),
        "HBN_is_exactly_frozen": (
            hbn_max
            <= FROZEN_COORDINATE_TOLERANCE_NM
        ),
        "PYR_is_exactly_frozen": (
            pyr_max
            <= FROZEN_COORDINATE_TOLERANCE_NM
        ),
        "CAPS_are_exactly_frozen": (
            caps_max
            <= FROZEN_COORDINATE_TOLERANCE_NM
        ),
        "water_displacement_is_finite": (
            math.isfinite(
                water_rms
            )
            and math.isfinite(
                water_max
            )
        ),
        "water_displacement_is_local": (
            water_max
            <= MAX_WATER_O_DISPLACEMENT_NM
        ),
        "initial_lumen_occupancy_is_428": (
            initial_lumen_count
            == EXPECTED_INITIAL_LUMEN_WATERS
        ),
        "final_lumen_occupancy_retains_at_least_98_percent": (
            final_lumen_count
            >= minimum_retained_lumen_waters
        ),
        "initial_lumen_identity_retains_at_least_98_percent": (
            retained_initial_lumen_count
            >= minimum_retained_lumen_waters
        ),
        "initial_CAP_OW_distance_is_safe": (
            initial_cap_ow_distance
            >= MIN_CAP_OW_DISTANCE_NM
        ),
        "final_CAP_OW_distance_is_safe": (
            final_cap_ow_distance
            >= MIN_CAP_OW_DISTANCE_NM
        ),
        "potential_series_is_finite": (
            np.all(
                np.isfinite(
                    potentials
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
        "CAP_SOL_Coulomb_is_zero": (
            float(
                np.max(
                    np.abs(
                        cap_sol_coul
                    )
                )
            )
            <= ZERO_ENERGY_TOLERANCE_KJ_MOL
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
        "R2_WATER_ONLY_ENERGY_MINIMIZATION_VALIDATED"
        if accepted
        else
        "R2_WATER_ONLY_ENERGY_MINIMIZATION_REQUIRES_REVIEW"
    )

    required_next_step = (
        "PREPARE_R2_FROZEN_SOLUTE_NVT_20PS"
        if accepted
        else
        "REVIEW_R2_WATER_ONLY_EM_GATE_FAILURES"
    )

    summary = {
        "decision": decision,
        "grompp_return_code": (
            grompp.returncode
        ),
        "grompp_warning_count": (
            warning_count
        ),
        "mdrun_return_code": (
            mdrun_return_code
        ),
        "EM_converged": (
            em_metrics[
                "converged"
            ]
        ),
        "EM_steps": (
            em_metrics[
                "steps"
            ]
        ),
        "EM_target_kJ_mol_nm": (
            em_metrics[
                "convergence_target_kJ_mol_nm"
            ]
        ),
        "maximum_force_kJ_mol_nm": (
            maximum_force
        ),
        "maximum_force_atom": (
            em_metrics[
                "maximum_force_atom"
            ]
        ),
        "norm_force_kJ_mol_nm": (
            norm_force
        ),
        "energy_record_count": (
            len(
                potential_data
            )
        ),
        "potential_term": (
            potential_name
        ),
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
        "CAP_SOL_LJ_term": (
            cap_sol_lj_name
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
        "CAP_SOL_Coulomb_term": (
            cap_sol_coul_name
        ),
        "CAP_SOL_Coulomb_maximum_absolute_kJ_mol": float(
            np.max(
                np.abs(
                    cap_sol_coul
                )
            )
        ),
        "HBN_RMS_displacement_nm": (
            hbn_rms
        ),
        "HBN_max_displacement_nm": (
            hbn_max
        ),
        "PYR_RMS_displacement_nm": (
            pyr_rms
        ),
        "PYR_max_displacement_nm": (
            pyr_max
        ),
        "CAPS_RMS_displacement_nm": (
            caps_rms
        ),
        "CAPS_max_displacement_nm": (
            caps_max
        ),
        "waterO_RMS_displacement_nm": (
            water_rms
        ),
        "waterO_max_displacement_nm": (
            water_max
        ),
        "initial_lumen_occupancy": (
            initial_lumen_count
        ),
        "final_lumen_occupancy": (
            final_lumen_count
        ),
        "initial_lumen_waters_retained": (
            retained_initial_lumen_count
        ),
        "initial_lumen_retention_fraction": (
            retained_initial_lumen_count
            / initial_lumen_count
        ),
        "initial_CAP_OW_distance_nm": (
            initial_cap_ow_distance
        ),
        "final_CAP_OW_distance_nm": (
            final_cap_ow_distance
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
        "short_frozen_solute_NVT_preparation_authorized": (
            accepted
        ),
        "short_frozen_solute_NVT_execution_authorized": False,
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
        GATE_CSV,
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
        f"""# R2 Water-Only Energy Minimization

## Scope

The R2 symmetric partial-cap system was minimized with HBN, all four
pyrenes, and both cap assemblies frozen in all dimensions.

Only the {WATERS} TIP4P/2005 water molecules were allowed to move.

## Execution

- Grompp return code:
  **{grompp.returncode}**
- Grompp warnings:
  **{warning_count}**
- Mdrun return code:
  **{mdrun_return_code}**
- Steepest-descents convergence:
  **{'YES' if em_metrics['converged'] else 'NO'}**
- EM steps:
  **{em_metrics['steps']}**
- Maximum force:
  **{maximum_force:.6f} kJ mol^-1 nm^-1**
- Norm of force:
  **{norm_force:.6f} kJ mol^-1 nm^-1**

## Frozen-group integrity

RMS/max displacement:

- HBN:
  **{hbn_rms:.12f}/{hbn_max:.12f} nm**
- PYR:
  **{pyr_rms:.12f}/{pyr_max:.12f} nm**
- CAPS:
  **{caps_rms:.12f}/{caps_max:.12f} nm**
- Water oxygen:
  **{water_rms:.9f}/{water_max:.9f} nm**

## Water confinement

- Initial lumen occupancy:
  **{initial_lumen_count}**
- Final lumen occupancy:
  **{final_lumen_count}**
- Initially luminal waters retained:
  **{retained_initial_lumen_count}/{initial_lumen_count}**
- Initial retention fraction:
  **{retained_initial_lumen_count / initial_lumen_count:.6f}**
- Initial CAP–OW distance:
  **{initial_cap_ow_distance:.6f} nm**
- Final CAP–OW distance:
  **{final_cap_ow_distance:.6f} nm**

## Energy

- Energy records:
  **{len(potential_data)}**
- Potential initial/final/change:
  **{potentials[0]:.6f}/
  {potentials[-1]:.6f}/
  {potentials[-1] - potentials[0]:.6f} kJ/mol**
- CAP–SOL LJ initial/final/maximum:
  **{cap_sol_lj[0]:.6f}/
  {cap_sol_lj[-1]:.6f}/
  {np.max(cap_sol_lj):.6f} kJ/mol**
- Maximum absolute CAP–SOL Coulomb:
  **{np.max(np.abs(cap_sol_coul)):.9f} kJ/mol**

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Short frozen-solute NVT preparation authorized:
  **{'YES' if accepted else 'NO'}**
- Short frozen-solute NVT execution authorized:
  **NO**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `{required_next_step}`

This result applies only to the neutral frozen steric R2 screening
model. It does not establish chemical realizability or long-time water
retention.
""",
        encoding="utf-8",
    )

    print()
    print(
        "Day023 R2 water-only energy minimization "
        "and validation completed."
    )

    print(
        "Grompp / mdrun return codes: "
        f"{grompp.returncode}/"
        f"{mdrun_return_code}"
    )

    print(
        "Grompp warnings / instability signatures: "
        f"{warning_count}/"
        + (
            "NONE"
            if not instability_hits
            else " | ".join(
                instability_hits
            )
        )
    )

    print(
        "EM converged / steps: "
        f"{'YES' if em_metrics['converged'] else 'NO'} / "
        f"{em_metrics['steps']}"
    )

    print(
        "Maximum force / norm: "
        f"{maximum_force:.6f}/"
        f"{norm_force:.6f} kJ mol^-1 nm^-1"
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
        "CAPS RMS/max displacement: "
        f"{caps_rms:.12f}/"
        f"{caps_max:.12f} nm"
    )

    print(
        "Water-O RMS/max displacement: "
        f"{water_rms:.9f}/"
        f"{water_max:.9f} nm"
    )

    print(
        "Lumen occupancy initial/final/retained: "
        f"{initial_lumen_count}/"
        f"{final_lumen_count}/"
        f"{retained_initial_lumen_count}"
    )

    print(
        "Initial/final minimum CAP-OW distance: "
        f"{initial_cap_ow_distance:.6f}/"
        f"{final_cap_ow_distance:.6f} nm"
    )

    print(
        "Potential initial/final/change: "
        f"{potentials[0]:.6f}/"
        f"{potentials[-1]:.6f}/"
        f"{potentials[-1] - potentials[0]:.6f} kJ/mol"
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
        "Short frozen-solute NVT preparation authorized: "
        f"{'YES' if accepted else 'NO'}"
    )

    print(
        "Short frozen-solute NVT execution authorized: NO"
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
        f"Wrote: {relative(final_gro)}"
    )

    print(
        f"Wrote: {relative(edr)}"
    )

    print(
        f"Wrote: {relative(log)}"
    )

    print(
        f"Wrote: {relative(ENERGY_SERIES_CSV)}"
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
            "R2 water-only energy minimization "
            "requires review."
        )


if __name__ == "__main__":
    main()
