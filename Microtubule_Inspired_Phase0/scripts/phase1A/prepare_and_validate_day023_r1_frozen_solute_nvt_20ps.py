#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
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

DAY023_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design"
)

PROTOCOL_ROOT = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

EM_ROOT = (
    DAY023_ROOT
    / "06_r1_water_only_em"
)

STATIC_ROOT = (
    DAY023_ROOT
    / "05_r1_full_static_scan"
)

PROTOTYPE_ROOT = (
    DAY023_ROOT
    / "02_r1_steric_cap_prototype"
)

OUTPUT_ROOT = (
    DAY023_ROOT
    / "07_r1_frozen_solute_nvt_20ps_preparation"
)

INPUT_GRO = (
    EM_ROOT
    / "r1_water_only_em.gro"
)

INPUT_TOP = (
    STATIC_ROOT
    / "selected/r1_selected_cap_model.top"
)

INPUT_INDEX = (
    STATIC_ROOT
    / "selected/r1_selected_groups.ndx"
)

SELECTED_MODEL_JSON = (
    STATIC_ROOT
    / "r1_selected_cap_model.json"
)

PROTOTYPE_JSON = (
    PROTOTYPE_ROOT
    / "r1_selected_steric_cap_definition.json"
)

EM_SUMMARY_CSV = (
    EM_ROOT
    / "r1_water_only_em_summary.csv"
)

EM_POTENTIAL_XVG = (
    EM_ROOT
    / "r1_water_only_em_potential.xvg"
)

TEMPLATE_CANDIDATES = (
    PROTOCOL_ROOT
    / "protocol_inputs/mdp/02_nvt_k10000_1ps.mdp",
    PROTOCOL_ROOT
    / "protocol_inputs/mdp/03_nvt_k1000_2ps.mdp",
)

LOCAL_GRO = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_input.gro"
)

LOCAL_TOP = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps.top"
)

LOCAL_INDEX = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps.ndx"
)

NVT_MDP = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps.mdp"
)

TPR = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps.tpr"
)

MDOUT = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_mdout.mdp"
)

PROCESSED_TOP = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_processed.top"
)

GROMPP_LOG = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_grompp.log"
)

TPR_DUMP_STDOUT = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_tpr_dump_stdout.txt"
)

TPR_DUMP_STDERR = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_tpr_dump_stderr.txt"
)

VELOCITY_CSV = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_velocity_audit.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_preparation_summary.csv"
)

GATE_CSV = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_20ps_preparation_gates.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R1_FROZEN_SOLUTE_NVT_20PS_PREPARATION_DAY023.md"
)

EXPECTED_ATOMS = 68314
HBN_ATOMS = 1680
PYR_ATOMS = 104
SOLUTE_ATOMS = HBN_ATOMS + PYR_ATOMS
EXPECTED_WATERS = 16551
WATER_SITES = 4
WATER_ATOMS = EXPECTED_WATERS * WATER_SITES
CAP_ATOMS = 326

DT_PS = 0.0005
DURATION_PS = 20.0
NSTEPS = int(round(DURATION_PS / DT_PS))

TRAJECTORY_INTERVAL_PS = 0.5
NSTXTC = int(
    round(
        TRAJECTORY_INTERVAL_PS / DT_PS
    )
)

EXPECTED_FRAMES = (
    NSTEPS // NSTXTC
    + 1
)

TEMPERATURE_K = 300.0
TAU_T_PS = 0.1
GEN_SEED = 20260708

MIN_WATER_NONZERO_VELOCITY_FRACTION = 0.70
MAX_WATER_NONZERO_VELOCITY_FRACTION = 0.80

CAP_SIGMA_EXPECTED_NM = -0.17
CAP_EPSILON_EXPECTED_KJ_MOL = 3.1179234818

VELOCITY_PATTERN = re.compile(
    r"^\s*v\[\s*(\d+)\s*\]\s*=\s*"
    r"\{\s*([-+0-9.eE]+)\s*,\s*"
    r"([-+0-9.eE]+)\s*,\s*"
    r"([-+0-9.eE]+)\s*\}"
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


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


def select_template() -> Path:
    candidates = [
        path
        for path in TEMPLATE_CANDIDATES
        if path.exists() and path.is_file()
    ]

    if not candidates:
        raise RuntimeError(
            "No validated NVT template was found."
        )

    return candidates[0]


def normalized_mdp_key(
    line: str,
) -> str | None:
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


def build_nvt_mdp(
    template: Path,
    destination: Path,
) -> None:
    overridden = {
        "integrator",
        "dt",
        "nsteps",
        "tinit",
        "init_step",
        "continuation",
        "define",
        "freezegrps",
        "freezedim",
        "tcoupl",
        "tc_grps",
        "tau_t",
        "ref_t",
        "pcoupl",
        "gen_vel",
        "gen_temp",
        "gen_seed",
        "comm_mode",
        "comm_grps",
        "nstcomm",
        "energygrps",
        "nstcalcenergy",
        "nstenergy",
        "nstlog",
        "nstxout",
        "nstvout",
        "nstfout",
        "nstxout_compressed",
        "compressed_x_precision",
    }

    retained = []

    for line in template.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        key = normalized_mdp_key(
            line
        )

        if key in overridden:
            continue

        retained.append(line)

    retained.extend(
        [
            "",
            "; Day023 R1 frozen-solute 20 ps screening",
            "integrator               = md",
            f"dt                       = {DT_PS:.7f}",
            f"nsteps                   = {NSTEPS}",
            "tinit                    = 0.0",
            "init-step                = 0",
            "continuation             = no",
            "define                   =",
            "",
            "freezegrps               = HBN_PYR CAPS",
            "freezedim                = Y Y Y Y Y Y",
            "",
            "tcoupl                   = V-rescale",
            "tc-grps                  = SOL HBN_PYR CAPS",
            f"tau-t                    = {TAU_T_PS:.3f} -1 -1",
            f"ref-t                    = {TEMPERATURE_K:.1f} {TEMPERATURE_K:.1f} {TEMPERATURE_K:.1f}",
            "pcoupl                   = no",
            "",
            "gen-vel                  = yes",
            f"gen-temp                 = {TEMPERATURE_K:.1f}",
            f"gen-seed                 = {GEN_SEED}",
            "",
            "comm-mode                = Linear",
            "comm-grps                = SOL",
            "nstcomm                  = 100",
            "",
            "energygrps               = HBN_PYR SOL CAPS",
            "nstcalcenergy            = 100",
            "nstenergy                = 100",
            "nstlog                   = 100",
            "nstxout                  = 0",
            "nstvout                  = 0",
            "nstfout                  = 0",
            f"nstxout-compressed       = {NSTXTC}",
            "compressed-x-precision   = 1000",
            "",
        ]
    )

    destination.write_text(
        "\n".join(retained),
        encoding="utf-8",
    )


def parse_mdp(
    path: Path,
) -> dict[str, str]:
    values: dict[str, str] = {}

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        active = raw_line.split(
            ";",
            1,
        )[0].strip()

        if (
            not active
            or active.startswith("#")
            or "=" not in active
        ):
            continue

        key, value = active.split(
            "=",
            1,
        )

        normalized = (
            key.strip()
            .lower()
            .replace("-", "_")
        )

        values[normalized] = (
            " ".join(
                value.strip().split()
            )
        )

    return values


def read_gro_atom_count(
    path: Path,
) -> int:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    if len(lines) < 2:
        raise RuntimeError(
            f"Malformed GRO file: {path}"
        )

    return int(
        lines[1].strip()
    )


def count_numeric_xvg_rows(
    path: Path,
) -> int:
    if not path.exists():
        return 0

    count = 0

    for raw_line in path.read_text(
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
            try:
                float(fields[0])
                float(fields[1])
            except ValueError:
                continue

            count += 1

    return count


def dump_tpr(
    gmx: str,
) -> int:
    with TPR_DUMP_STDOUT.open(
        "w",
        encoding="utf-8",
    ) as stdout_handle, TPR_DUMP_STDERR.open(
        "w",
        encoding="utf-8",
    ) as stderr_handle:
        completed = subprocess.run(
            [
                gmx,
                "dump",
                "-s",
                str(TPR),
            ],
            cwd=OUTPUT_ROOT,
            env={
                **os.environ,
                "GMX_MAXBACKUP": "-1",
            },
            text=True,
            stdout=stdout_handle,
            stderr=stderr_handle,
            check=False,
        )

    return completed.returncode


def parse_tpr_dump(
    path: Path,
) -> tuple[
    int | None,
    np.ndarray,
]:
    natoms: int | None = None

    velocities: dict[
        int,
        tuple[
            float,
            float,
            float,
        ],
    ] = {}

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        for line in handle:
            if natoms is None:
                match = re.search(
                    r"\bnatoms\s*=\s*(\d+)",
                    line,
                )

                if match is not None:
                    natoms = int(
                        match.group(1)
                    )

            velocity_match = (
                VELOCITY_PATTERN.match(
                    line
                )
            )

            if velocity_match is None:
                continue

            index = int(
                velocity_match.group(1)
            )

            velocities[index] = (
                float(
                    velocity_match.group(2)
                ),
                float(
                    velocity_match.group(3)
                ),
                float(
                    velocity_match.group(4)
                ),
            )

    if velocities:
        maximum_index = max(
            velocities
        )

        array = np.full(
            (
                maximum_index + 1,
                3,
            ),
            np.nan,
            dtype=float,
        )

        for index, vector in velocities.items():
            array[index] = vector
    else:
        array = np.empty(
            (0, 3),
            dtype=float,
        )

    return natoms, array


def velocity_metrics(
    velocities: np.ndarray,
    first: int,
    last: int,
) -> dict[str, float]:
    subset = velocities[
        first:last
    ]

    if len(subset) == 0:
        raise RuntimeError(
            "Empty velocity group."
        )

    if not np.all(
        np.isfinite(subset)
    ):
        raise RuntimeError(
            "Non-finite or missing velocities "
            f"in range {first}:{last}."
        )

    norms = np.linalg.norm(
        subset,
        axis=1,
    )

    nonzero = (
        norms > 1.0e-12
    )

    return {
        "atom_count": float(
            len(subset)
        ),
        "nonzero_count": float(
            np.count_nonzero(
                nonzero
            )
        ),
        "nonzero_fraction": float(
            np.mean(
                nonzero
            )
        ),
        "rms_speed_nm_ps": float(
            np.sqrt(
                np.mean(
                    norms * norms
                )
            )
        ),
        "maximum_speed_nm_ps": float(
            np.max(norms)
        ),
    }


def count_active_sections(
    path: Path,
    target: str,
) -> int:
    pattern = re.compile(
        r"^\s*\[\s*"
        + re.escape(target)
        + r"\s*\]\s*$",
        flags=re.IGNORECASE,
    )

    count = 0

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        line = raw_line.split(
            ";",
            1,
        )[0].strip()

        if pattern.match(line):
            count += 1

    return count


def parse_processed_topology(
    path: Path,
) -> dict[str, Any]:
    current_section = ""
    molecule_counts: dict[str, int] = {}
    cap_ow_overrides = []

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        line = raw_line.split(
            ";",
            1,
        )[0].strip()

        if not line:
            continue

        section_match = re.match(
            r"^\[\s*([^\]]+?)\s*\]$",
            line,
        )

        if section_match is not None:
            current_section = (
                section_match.group(1)
                .strip()
                .lower()
            )

            continue

        if line.startswith("#"):
            continue

        fields = line.split()

        if (
            current_section
            == "nonbond_params"
            and len(fields) >= 5
        ):
            type_i = fields[0]
            type_j = fields[1]

            if {
                type_i.upper(),
                type_j.upper(),
            } == {
                "CAP",
                "OW",
            }:
                cap_ow_overrides.append(
                    {
                        "type_i": type_i,
                        "type_j": type_j,
                        "function": int(
                            fields[2]
                        ),
                        "sigma_nm": float(
                            fields[3]
                        ),
                        "epsilon_kJ_mol": float(
                            fields[4]
                        ),
                    }
                )

        elif (
            current_section
            == "molecules"
            and len(fields) >= 2
        ):
            try:
                molecule_counts[
                    fields[0]
                ] = int(
                    fields[1]
                )
            except ValueError:
                continue

    return {
        "molecule_counts": (
            molecule_counts
        ),
        "cap_ow_overrides": (
            cap_ow_overrides
        ),
    }


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


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        INPUT_GRO,
        INPUT_TOP,
        INPUT_INDEX,
        SELECTED_MODEL_JSON,
        PROTOTYPE_JSON,
        EM_SUMMARY_CSV,
    ):
        require_file(required)

    model = json.loads(
        SELECTED_MODEL_JSON.read_text(
            encoding="utf-8"
        )
    )

    prototype = json.loads(
        PROTOTYPE_JSON.read_text(
            encoding="utf-8"
        )
    )

    em_summary = read_single_csv_row(
        EM_SUMMARY_CSV
    )

    if (
        model.get("decision")
        !=
        "FULL_R1_STATIC_CAP_WATER_MODEL_VALIDATED"
    ):
        raise RuntimeError(
            "The selected R1 cap model is not valid."
        )

    if (
        em_summary.get("decision")
        != "R1_WATER_ONLY_EM_VALIDATED"
    ):
        raise RuntimeError(
            "The R1 water-only EM is not accepted."
        )

    if not parse_bool(
        em_summary.get(
            "short_NVT_preparation_authorized",
            "false",
        )
    ):
        raise RuntimeError(
            "Short NVT preparation is not authorized."
        )

    if int(
        model[
            "R1_atom_count"
        ]
    ) != EXPECTED_ATOMS:
        raise RuntimeError(
            "Unexpected selected-model atom count."
        )

    if int(
        prototype[
            "retained_water_molecules"
        ]
    ) != EXPECTED_WATERS:
        raise RuntimeError(
            "Unexpected retained-water count."
        )

    if (
        HBN_ATOMS
        + PYR_ATOMS
        + WATER_ATOMS
        + CAP_ATOMS
        != EXPECTED_ATOMS
    ):
        raise RuntimeError(
            "Internal atom accounting failed."
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

    build_nvt_mdp(
        template,
        NVT_MDP,
    )

    grompp = run_command(
        [
            gmx,
            "grompp",
            "-f",
            str(NVT_MDP),
            "-c",
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
            "1",
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
            "R1 20 ps NVT grompp failed. "
            f"See {GROMPP_LOG}"
        )

    dump_return_code = dump_tpr(
        gmx
    )

    if dump_return_code != 0:
        raise RuntimeError(
            "gmx dump failed for the R1 NVT TPR."
        )

    mdp_values = parse_mdp(
        MDOUT
    )

    tpr_natoms, velocities = parse_tpr_dump(
        TPR_DUMP_STDOUT
    )

    if velocities.shape != (
        EXPECTED_ATOMS,
        3,
    ):
        raise RuntimeError(
            "Unexpected velocity array shape: "
            f"{velocities.shape}; expected "
            f"({EXPECTED_ATOMS}, 3)"
        )

    ranges = {
        "HBN": (
            0,
            HBN_ATOMS,
        ),
        "PYR": (
            HBN_ATOMS,
            SOLUTE_ATOMS,
        ),
        "SOL": (
            SOLUTE_ATOMS,
            SOLUTE_ATOMS
            + WATER_ATOMS,
        ),
        "CAPS": (
            SOLUTE_ATOMS
            + WATER_ATOMS,
            EXPECTED_ATOMS,
        ),
    }

    velocity_rows = []

    for group, (
        first,
        last,
    ) in ranges.items():
        metrics = velocity_metrics(
            velocities,
            first,
            last,
        )

        velocity_rows.append(
            {
                "group": group,
                "first_atom_one_based": (
                    first + 1
                ),
                "last_atom_one_based": last,
                **metrics,
            }
        )

    write_csv(
        VELOCITY_CSV,
        velocity_rows,
    )

    velocity_by_group = {
        row["group"]: row
        for row in velocity_rows
    }

    processed = parse_processed_topology(
        PROCESSED_TOP
    )

    molecule_counts = processed[
        "molecule_counts"
    ]

    cap_ow_overrides = processed[
        "cap_ow_overrides"
    ]

    cap_override_valid = bool(
        len(cap_ow_overrides) == 1
        and math.isclose(
            cap_ow_overrides[0][
                "sigma_nm"
            ],
            CAP_SIGMA_EXPECTED_NM,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
        and math.isclose(
            cap_ow_overrides[0][
                "epsilon_kJ_mol"
            ],
            CAP_EPSILON_EXPECTED_KJ_MOL,
            rel_tol=1.0e-9,
            abs_tol=1.0e-9,
        )
    )

    position_restraint_sections = (
        count_active_sections(
            PROCESSED_TOP,
            "position_restraints",
        )
    )

    gro_atoms = read_gro_atom_count(
        LOCAL_GRO
    )

    em_energy_points = (
        count_numeric_xvg_rows(
            EM_POTENTIAL_XVG
        )
    )

    em_energy_change_interpretable = (
        em_energy_points >= 2
    )

    expected_freezegrps = (
        mdp_values.get(
            "freezegrps",
            "",
        )
    )

    expected_freezedim = (
        mdp_values.get(
            "freezedim",
            "",
        )
    )

    water_nonzero_fraction = float(
        velocity_by_group[
            "SOL"
        ][
            "nonzero_fraction"
        ]
    )

    all_velocity_values_finite = bool(
        np.all(
            np.isfinite(
                velocities
            )
        )
    )

    tc_group_values = (
        mdp_values.get(
            "tc_grps",
            "",
        ).split()
    )

    tau_t_values = [
        float(value)
        for value in mdp_values.get(
            "tau_t",
            "",
        ).split()
    ]

    ref_t_values = [
        float(value)
        for value in mdp_values.get(
            "ref_t",
            "",
        ).split()
    ]

    gates = {
        "grompp_return_code_zero": (
            grompp.returncode == 0
        ),
        "TPR_dump_return_code_zero": (
            dump_return_code == 0
        ),
        "input_GRO_atom_count_is_68314": (
            gro_atoms == EXPECTED_ATOMS
        ),
        "TPR_atom_count_is_68314": (
            tpr_natoms == EXPECTED_ATOMS
        ),
        "all_TPR_velocities_parsed": (
            velocities.shape
            == (
                EXPECTED_ATOMS,
                3,
            )
        ),
        "all_TPR_velocities_are_finite": (
            all_velocity_values_finite
        ),
        "water_velocity_fraction_is_TIP4P_consistent": (
            MIN_WATER_NONZERO_VELOCITY_FRACTION
            <= water_nonzero_fraction
            <= MAX_WATER_NONZERO_VELOCITY_FRACTION
        ),
        "integrator_is_md": (
            mdp_values.get(
                "integrator",
                "",
            ).lower()
            == "md"
        ),
        "dt_is_0p0005_ps": (
            math.isclose(
                float(
                    mdp_values.get(
                        "dt",
                        "nan",
                    )
                ),
                DT_PS,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
        ),
        "nsteps_is_40000": (
            int(
                mdp_values.get(
                    "nsteps",
                    "-1",
                )
            )
            == NSTEPS
        ),
        "continuation_is_no": (
            mdp_values.get(
                "continuation",
                "",
            ).lower()
            == "no"
        ),
        "generation_temperature_is_300K": (
            math.isclose(
                float(
                    mdp_values.get(
                        "gen_temp",
                        "nan",
                    )
                ),
                TEMPERATURE_K,
                rel_tol=0.0,
                abs_tol=1.0e-8,
            )
        ),
        "generation_seed_is_fixed": (
            int(
                mdp_values.get(
                    "gen_seed",
                    "-1",
                )
            )
            == GEN_SEED
        ),
        "freeze_groups_are_HBN_PYR_and_CAPS": (
            expected_freezegrps
            == "HBN_PYR CAPS"
        ),
        "all_freeze_dimensions_are_enabled": (
            expected_freezedim
            == "Y Y Y Y Y Y"
        ),
        "temperature_groups_partition_mobile_and_frozen_atoms": (
            tc_group_values
            == [
                "SOL",
                "HBN_PYR",
                "CAPS",
            ]
        ),
        "mobile_group_tau_t_is_0p1_ps": (
            len(
                tau_t_values
            )
            == 3
            and math.isclose(
                tau_t_values[
                    0
                ],
                TAU_T_PS,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
        ),
        "frozen_groups_are_not_thermostatted": (
            len(
                tau_t_values
            )
            == 3
            and tau_t_values[
                1:
            ]
            == [
                -1.0,
                -1.0,
            ]
        ),
        "reference_temperatures_have_three_300K_entries": (
            len(
                ref_t_values
            )
            == 3
            and all(
                math.isclose(
                    value,
                    TEMPERATURE_K,
                    rel_tol=0.0,
                    abs_tol=1.0e-8,
                )
                for value in ref_t_values
            )
        ),
        "pressure_coupling_is_disabled": (
            mdp_values.get(
                "pcoupl",
                "",
            ).lower()
            == "no"
        ),
        "trajectory_interval_is_0p5ps": (
            int(
                mdp_values.get(
                    "nstxout_compressed",
                    "-1",
                )
            )
            == NSTXTC
        ),
        "no_active_position_restraints": (
            position_restraint_sections == 0
        ),
        "SOL_count_is_16551": (
            molecule_counts.get(
                "SOL"
            )
            == EXPECTED_WATERS
        ),
        "CAPL_count_is_one": (
            molecule_counts.get(
                "CAPL"
            )
            == 1
        ),
        "CAPU_count_is_one": (
            molecule_counts.get(
                "CAPU"
            )
            == 1
        ),
        "selected_CAP_OW_override_is_preserved": (
            cap_override_valid
        ),
    }

    failed_gates = [
        name
        for name, passed in gates.items()
        if not passed
    ]

    prepared = (
        len(failed_gates) == 0
    )

    decision = (
        "R1_FROZEN_SOLUTE_NVT_20PS_PREPARED"
        if prepared
        else
        "R1_FROZEN_SOLUTE_NVT_20PS_PREPARATION_REQUIRES_REVIEW"
    )

    required_next_step = (
        "RUN_R1_FROZEN_SOLUTE_NVT_20PS"
        if prepared
        else
        "RESOLVE_R1_NVT_PREPARATION_GATE_FAILURES"
    )

    summary = {
        "decision": decision,
        "template_mdp": relative(
            template
        ),
        "grompp_return_code": (
            grompp.returncode
        ),
        "TPR_dump_return_code": (
            dump_return_code
        ),
        "GRO_atom_count": gro_atoms,
        "TPR_atom_count": tpr_natoms,
        "dt_ps": DT_PS,
        "nsteps": NSTEPS,
        "duration_ps": DURATION_PS,
        "trajectory_interval_ps": (
            TRAJECTORY_INTERVAL_PS
        ),
        "expected_trajectory_frames": (
            EXPECTED_FRAMES
        ),
        "temperature_K": TEMPERATURE_K,
        "tau_t_ps": TAU_T_PS,
        "generation_seed": GEN_SEED,
        "freeze_groups": (
            expected_freezegrps
        ),
        "freeze_dimensions": (
            expected_freezedim
        ),
        "temperature_groups": (
            mdp_values.get(
                "tc_grps",
                "",
            )
        ),
        "temperature_coupling_time_constants_ps": (
            mdp_values.get(
                "tau_t",
                "",
            )
        ),
        "reference_temperatures_K": (
            mdp_values.get(
                "ref_t",
                "",
            )
        ),
        "water_nonzero_velocity_fraction": (
            water_nonzero_fraction
        ),
        "water_RMS_speed_nm_ps": (
            velocity_by_group[
                "SOL"
            ][
                "rms_speed_nm_ps"
            ]
        ),
        "HBN_nonzero_velocity_fraction": (
            velocity_by_group[
                "HBN"
            ][
                "nonzero_fraction"
            ]
        ),
        "PYR_nonzero_velocity_fraction": (
            velocity_by_group[
                "PYR"
            ][
                "nonzero_fraction"
            ]
        ),
        "CAPS_nonzero_velocity_fraction": (
            velocity_by_group[
                "CAPS"
            ][
                "nonzero_fraction"
            ]
        ),
        "active_position_restraint_sections": (
            position_restraint_sections
        ),
        "processed_SOL_count": (
            molecule_counts.get(
                "SOL",
                "",
            )
        ),
        "processed_CAPL_count": (
            molecule_counts.get(
                "CAPL",
                "",
            )
        ),
        "processed_CAPU_count": (
            molecule_counts.get(
                "CAPU",
                "",
            )
        ),
        "CAP_OW_override_count": (
            len(cap_ow_overrides)
        ),
        "CAP_OW_override_valid": (
            cap_override_valid
        ),
        "EM_potential_energy_points": (
            em_energy_points
        ),
        "EM_initial_final_energy_change_interpretable": (
            em_energy_change_interpretable
        ),
        "failed_gates": (
            " | ".join(
                failed_gates
            )
        ),
        "NVT_execution_authorized": (
            prepared
        ),
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

    velocity_lines = "\n".join(
        (
            f"- `{row['group']}`: "
            f"atoms={int(row['atom_count'])}; "
            f"nonzero fraction="
            f"{row['nonzero_fraction']:.6f}; "
            f"RMS speed="
            f"{row['rms_speed_nm_ps']:.6f} nm/ps"
        )
        for row in velocity_rows
    )

    REPORT_MD.write_text(
        f"""# R1 Frozen-Solute NVT 20 ps Preparation

## Scope

This gate prepared and statically validated the first short R1
confinement screening trajectory.

No molecular dynamics was executed.

## Inputs

- Minimized coordinates:
  `{relative(LOCAL_GRO)}`
- Selected topology:
  `{relative(LOCAL_TOP)}`
- Index:
  `{relative(LOCAL_INDEX)}`
- NVT template:
  `{relative(template)}`

## Protocol

- Integrator: **MD**
- Temperature: **{TEMPERATURE_K:.1f} K**
- Time step: **{DT_PS:.7f} ps**
- Steps: **{NSTEPS}**
- Duration: **{DURATION_PS:.1f} ps**
- Trajectory interval:
  **{TRAJECTORY_INTERVAL_PS:.1f} ps**
- Expected frames:
  **{EXPECTED_FRAMES}**
- Velocity seed:
  **{GEN_SEED}**
- Frozen groups:
  **HBN_PYR CAPS**
- T-coupling partition:
  **SOL HBN_PYR CAPS**
- Mobile and thermostatted group:
  **SOL, tau-t = 0.100 ps**
- Frozen and unthermostatted groups:
  **HBN_PYR and CAPS, tau-t = -1 ps**
- Pressure coupling:
  **disabled**

## TPR and topology audit

- GRO atoms:
  **{gro_atoms}**
- TPR atoms:
  **{tpr_natoms}**
- Active position-restraint sections:
  **{position_restraint_sections}**
- Processed SOL/CAPL/CAPU counts:
  **{molecule_counts.get('SOL')}/
  {molecule_counts.get('CAPL')}/
  {molecule_counts.get('CAPU')}**
- CAP–OW override preserved:
  **{'YES' if cap_override_valid else 'NO'}**

## Velocity audit

{velocity_lines}

The expected nonzero fraction for four-site TIP4P/2005 water is
approximately 0.75 because the virtual M site has no independently
generated velocity.

Frozen-group velocities are recorded for provenance. Their coordinates
remain fixed by the six active freeze dimensions.

## Energy-minimization reporting note

- Potential-energy records in the EM XVG:
  **{em_energy_points}**
- Initial/final EM energy change independently resolvable:
  **{'YES' if em_energy_change_interpretable else 'NO'}**

If only one energy record is present, the previously reported zero
energy change is a sampling limitation, not evidence that the
minimization failed. EM acceptance remains based on explicit
convergence, final force, structural integrity, and preserved hydration.

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- NVT execution authorized:
  **{'YES' if prepared else 'NO'}**
- Required next step:
  `{required_next_step}`

The 20 ps run is an initial positive-control screening. Extension to a
longer frozen-solute trajectory will require validation of temperature,
cap integrity, lumen occupancy, zero-occupancy fraction, and axial
water retention.
""",
        encoding="utf-8",
    )

    print(
        "Day023 R1 frozen-solute NVT 20 ps "
        "preparation completed."
    )

    print(
        "Grompp / TPR dump return codes: "
        f"{grompp.returncode}/"
        f"{dump_return_code}"
    )

    print(
        "GRO / TPR atoms: "
        f"{gro_atoms}/"
        f"{tpr_natoms}"
    )

    print(
        "Protocol dt / steps / duration: "
        f"{DT_PS:.7f} ps / "
        f"{NSTEPS} / "
        f"{DURATION_PS:.1f} ps"
    )

    print(
        "Trajectory interval / expected frames: "
        f"{TRAJECTORY_INTERVAL_PS:.1f} ps / "
        f"{EXPECTED_FRAMES}"
    )

    print(
        "Freeze groups / dimensions: "
        f"{expected_freezegrps} / "
        f"{expected_freezedim}"
    )

    print(
        "T-coupling groups / reference temperatures / seed: "
        f"{mdp_values.get('tc_grps', '')} / "
        f"{mdp_values.get('ref_t', '')} K / "
        f"{mdp_values.get('gen_seed', '')}"
    )

    for row in velocity_rows:
        print(
            "Velocity group "
            f"{row['group']}: "
            f"nonzero={row['nonzero_fraction']:.6f}; "
            f"RMS={row['rms_speed_nm_ps']:.6f} nm/ps; "
            f"max={row['maximum_speed_nm_ps']:.6f} nm/ps"
        )

    print(
        "Active position-restraint sections: "
        f"{position_restraint_sections}"
    )

    print(
        "Processed SOL/CAPL/CAPU counts: "
        f"{molecule_counts.get('SOL')}/"
        f"{molecule_counts.get('CAPL')}/"
        f"{molecule_counts.get('CAPU')}"
    )

    print(
        "CAP-OW override count / valid: "
        f"{len(cap_ow_overrides)} / "
        f"{'YES' if cap_override_valid else 'NO'}"
    )

    print(
        "EM potential-energy points / "
        "energy-change interpretable: "
        f"{em_energy_points} / "
        f"{'YES' if em_energy_change_interpretable else 'NO'}"
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
        "NVT execution authorized: "
        f"{'YES' if prepared else 'NO'}"
    )

    print(
        f"Required next step: {required_next_step}"
    )

    print(
        f"Wrote: {relative(NVT_MDP)}"
    )

    print(
        f"Wrote: {relative(TPR)}"
    )

    print(
        f"Wrote: {relative(PROCESSED_TOP)}"
    )

    print(
        f"Wrote: {relative(VELOCITY_CSV)}"
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

    if not prepared:
        raise RuntimeError(
            "R1 20 ps NVT preparation requires review."
        )


if __name__ == "__main__":
    main()
