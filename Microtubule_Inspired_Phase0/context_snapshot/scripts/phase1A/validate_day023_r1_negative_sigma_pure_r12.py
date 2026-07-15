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

TOPOLOGY_MODEL_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design/"
    "03_r1_topology_model"
)

OUTPUT_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design/"
    "04_r1_negative_sigma_validation"
)

ATOMTYPE_CSV = (
    TOPOLOGY_MODEL_ROOT
    / "r1_active_atomtype_inventory.csv"
)

PREVIOUS_CONTRACT = (
    TOPOLOGY_MODEL_ROOT
    / "r1_cap_nonbonded_model_contract.json"
)

VALIDATION_CSV = (
    OUTPUT_ROOT
    / "r1_negative_sigma_pure_r12_validation.csv"
)

CANDIDATE_SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r1_negative_sigma_candidate_summary.csv"
)

CORRECTED_CONTRACT_JSON = (
    OUTPUT_ROOT
    / "r1_cap_nonbonded_model_contract_corrected.json"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R1_NEGATIVE_SIGMA_PURE_R12_VALIDATION_DAY023.md"
)

TEMPERATURE_K = 300.0
GAS_CONSTANT_KJ_MOL_K = 8.31446261815324e-3
KBT_KJ_MOL = (
    GAS_CONSTANT_KJ_MOL_K
    * TEMPERATURE_K
)

SIGMA_ABS_NM = 0.17

TARGET_LEVELS_KBT = (
    5.0,
    10.0,
    20.0,
    40.0,
)

DISTANCES_NM = (
    0.150,
    0.170,
    0.200,
    0.220,
    0.250,
    0.300,
    0.500,
)

BOX_LENGTH_NM = 3.0
VDW_CUTOFF_NM = 1.0

MAX_RELATIVE_ERROR = 2.0e-3
MAX_ABSOLUTE_ERROR_KJ_MOL = 2.0e-4

EXPECTED_COMB_RULE = 2
EXPECTED_NBFUNC = 1
EXPECTED_WATER_TYPE = "OW"


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


def read_csv_rows(
    path: Path,
) -> list[dict[str, str]]:
    require_file(path)

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        return list(
            csv.DictReader(handle)
        )


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


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


def resolve_water_oxygen_parameters() -> tuple[
    float,
    float,
]:
    rows = read_csv_rows(
        ATOMTYPE_CSV
    )

    candidates = [
        row
        for row in rows
        if (
            row.get(
                "atomtype",
                "",
            ).strip()
            == EXPECTED_WATER_TYPE
            or parse_bool(
                row.get(
                    "used_by_water_oxygen",
                    "false",
                )
            )
        )
    ]

    unique = {
        (
            row[
                "atomtype"
            ].strip(),
            float(
                row[
                    "parameter_v"
                ]
            ),
            float(
                row[
                    "parameter_w"
                ]
            ),
        )
        for row in candidates
    }

    if len(unique) != 1:
        raise RuntimeError(
            "Could not uniquely resolve OW sigma/epsilon: "
            + repr(
                sorted(unique)
            )
        )

    atomtype, sigma, epsilon = (
        next(
            iter(unique)
        )
    )

    if atomtype != EXPECTED_WATER_TYPE:
        raise RuntimeError(
            f"Expected water type OW; found {atomtype}"
        )

    if (
        sigma <= 0.0
        or epsilon <= 0.0
    ):
        raise RuntimeError(
            "The baseline OW sigma/epsilon are not "
            "positive."
        )

    return (
        sigma,
        epsilon,
    )


def write_topology(
    path: Path,
    ow_sigma_nm: float,
    ow_epsilon_kj_mol: float,
    pair_epsilon_kj_mol: float,
) -> None:
    path.write_text(
        f"""[ defaults ]
; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ
1         2          no         1.0      1.0

[ atomtypes ]
; name  at.num  mass       charge  ptype  sigma       epsilon
CAP     0       12.011000  0.0000  A      0.00000000  0.00000000
OWT     8       15.999400  0.0000  A      {ow_sigma_nm:.12e}  {ow_epsilon_kj_mol:.12e}

[ nonbond_params ]
; i    j    funct  sigma                    epsilon
CAP    OWT  1      {-SIGMA_ABS_NM:.12e}     {pair_epsilon_kj_mol:.12e}

[ moleculetype ]
; name  nrexcl
CAPM    0

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1     CAP   1      CAP      CAP   1     0.0     12.011

[ moleculetype ]
; name  nrexcl
OWM     0

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge  mass
1     OWT   1      OWM      OW    1     0.0     15.9994

[ system ]
R1 negative-sigma pure-r12 CAP-OW validation

[ molecules ]
CAPM    1
OWM     1
""",
        encoding="utf-8",
    )


def gro_frame(
    distance_nm: float,
    time_ps: float,
) -> str:
    first = np.asarray(
        [
            1.0,
            1.0,
            1.0,
        ],
        dtype=float,
    )

    second = np.asarray(
        [
            1.0 + distance_nm,
            1.0,
            1.0,
        ],
        dtype=float,
    )

    return (
        f"CAP-OW distance scan t={time_ps:.3f}\n"
        "2\n"
        f"{1:5d}{'CAP':<5s}{'CAP':>5s}{1:5d}"
        f"{first[0]:8.3f}{first[1]:8.3f}{first[2]:8.3f}\n"
        f"{2:5d}{'OWM':<5s}{'OW':>5s}{2:5d}"
        f"{second[0]:8.3f}{second[1]:8.3f}{second[2]:8.3f}\n"
        f"{BOX_LENGTH_NM:10.5f}"
        f"{BOX_LENGTH_NM:10.5f}"
        f"{BOX_LENGTH_NM:10.5f}\n"
    )


def write_trajectory(
    path: Path,
) -> None:
    frames = [
        gro_frame(
            distance_nm,
            float(index),
        )
        for index, distance_nm in enumerate(
            DISTANCES_NM
        )
    ]

    path.write_text(
        "".join(frames),
        encoding="utf-8",
    )


def write_index(
    path: Path,
) -> None:
    path.write_text(
        """[ System ]
1 2

[ CAP ]
1

[ OW ]
2
""",
        encoding="utf-8",
    )


def write_mdp(
    path: Path,
) -> None:
    path.write_text(
        f"""integrator               = md
dt                       = 0.001
nsteps                   = 1
continuation             = no

cutoff-scheme            = Verlet
nstlist                  = 1
rlist                    = {VDW_CUTOFF_NM:.3f}

coulombtype              = Cut-off
coulomb-modifier         = None
rcoulomb                 = {VDW_CUTOFF_NM:.3f}

vdwtype                  = Cut-off
vdw-modifier             = None
rvdw                     = {VDW_CUTOFF_NM:.3f}
DispCorr                 = no

constraints              = none
pbc                      = xyz
periodic-molecules       = no

tcoupl                   = no
pcoupl                   = no
gen-vel                  = no
comm-mode                = none

nstcalcenergy            = 1
nstenergy                = 1
nstlog                   = 1
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstxout-compressed       = 0

energygrps               = CAP OW
""",
        encoding="utf-8",
    )


def parse_numeric_xvg(
    path: Path,
) -> np.ndarray:
    require_file(path)

    rows = []

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
        or len(array) != len(
            DISTANCES_NM
        )
    ):
        raise RuntimeError(
            f"Unexpected energy XVG shape for {path}: "
            f"{array.shape}; expected "
            f"({len(DISTANCES_NM)}, 2)"
        )

    if not np.all(
        np.isfinite(array)
    ):
        raise RuntimeError(
            f"Non-finite energy values in {path}"
        )

    return array


def extract_pair_energy(
    gmx: str,
    working_directory: Path,
    edr: Path,
    output_xvg: Path,
) -> tuple[
    np.ndarray,
    str,
    str,
]:
    selectors = (
        "LJ-SR:CAP-OW",
        "LJ-SR:OW-CAP",
    )

    attempts = []

    for selector in selectors:
        if output_xvg.exists():
            output_xvg.unlink()

        completed = run_command(
            [
                gmx,
                "energy",
                "-f",
                str(edr),
                "-o",
                str(output_xvg),
                "-xvg",
                "none",
                "-dp",
            ],
            cwd=working_directory,
            input_text=(
                f"{selector}\n"
                "0\n"
            ),
        )

        attempts.append(
            (
                selector,
                completed.returncode,
                completed.stdout,
            )
        )

        if (
            completed.returncode == 0
            and output_xvg.exists()
            and output_xvg.stat().st_size > 0
        ):
            try:
                values = parse_numeric_xvg(
                    output_xvg
                )
            except Exception:
                continue

            return (
                values,
                selector,
                completed.stdout,
            )

    diagnostic = "\n\n".join(
        (
            f"SELECTOR: {selector}\n"
            f"RETURN CODE: {return_code}\n"
            f"{output}"
        )
        for selector, return_code, output
        in attempts
    )

    raise RuntimeError(
        "Could not extract CAP-OW LJ-SR energy.\n"
        + diagnostic
    )


def parse_dump_nbfp(
    dump_text: str,
) -> tuple[
    bool,
    bool,
]:
    lower = dump_text.lower()

    contains_negative_sigma = (
        f"{-SIGMA_ABS_NM:.5e}"
        in lower
        or f"{-SIGMA_ABS_NM:.6f}"
        in lower
        or "-0.17"
        in lower
    )

    contains_nonbond_data = any(
        token in lower
        for token in (
            "nbfp",
            "nonbond",
            "lj",
        )
    )

    return (
        contains_negative_sigma,
        contains_nonbond_data,
    )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        ATOMTYPE_CSV,
        PREVIOUS_CONTRACT,
    ):
        require_file(required)

    previous_contract = json.loads(
        PREVIOUS_CONTRACT.read_text(
            encoding="utf-8"
        )
    )

    if int(
        previous_contract[
            "nbfunc"
        ]
    ) != EXPECTED_NBFUNC:
        raise RuntimeError(
            "Unexpected baseline nbfunc."
        )

    if int(
        previous_contract[
            "combination_rule"
        ]
    ) != EXPECTED_COMB_RULE:
        raise RuntimeError(
            "Unexpected baseline combination rule."
        )

    if (
        previous_contract[
            "water_oxygen_atom_type"
        ]
        != EXPECTED_WATER_TYPE
    ):
        raise RuntimeError(
            "Unexpected water oxygen atom type."
        )

    ow_sigma_nm, ow_epsilon_kj_mol = (
        resolve_water_oxygen_parameters()
    )

    gmx = locate_gmx()

    validation_rows: list[
        dict[str, Any]
    ] = []

    summary_rows: list[
        dict[str, Any]
    ] = []

    all_candidates_pass = True

    for target_level_kbt in TARGET_LEVELS_KBT:
        candidate_name = (
            f"target_{target_level_kbt:.0f}kBT"
        )

        candidate_root = (
            OUTPUT_ROOT
            / candidate_name
        )

        candidate_root.mkdir(
            parents=True,
            exist_ok=True,
        )

        topology = (
            candidate_root
            / "cap_ow_pure_r12.top"
        )

        mdp = (
            candidate_root
            / "rerun.mdp"
        )

        index = (
            candidate_root
            / "cap_ow.ndx"
        )

        first_gro = (
            candidate_root
            / "first_frame.gro"
        )

        trajectory = (
            candidate_root
            / "distance_scan.gro"
        )

        tpr = (
            candidate_root
            / "cap_ow_pure_r12.tpr"
        )

        mdout = (
            candidate_root
            / "mdout.mdp"
        )

        processed_top = (
            candidate_root
            / "processed.top"
        )

        grompp_log = (
            candidate_root
            / "grompp.log"
        )

        dump_stdout = (
            candidate_root
            / "tpr_dump_stdout.txt"
        )

        dump_stderr = (
            candidate_root
            / "tpr_dump_stderr.txt"
        )

        rerun_log = (
            candidate_root
            / "rerun_console.log"
        )

        energy_xvg = (
            candidate_root
            / "cap_ow_lj_energy.xvg"
        )

        extraction_log = (
            candidate_root
            / "energy_extraction.log"
        )

        target_energy_kj_mol = (
            target_level_kbt
            * KBT_KJ_MOL
        )

        # With negative sigma under comb-rule 2,
        # GROMACS sets C6=0 and evaluates:
        # C12 = 4 epsilon |sigma|^12.
        pair_epsilon_kj_mol = (
            target_energy_kj_mol
            / 4.0
        )

        c12_kj_mol_nm12 = (
            4.0
            * pair_epsilon_kj_mol
            * SIGMA_ABS_NM
            ** 12
        )

        write_topology(
            topology,
            ow_sigma_nm,
            ow_epsilon_kj_mol,
            pair_epsilon_kj_mol,
        )

        write_mdp(
            mdp
        )

        write_index(
            index
        )

        first_gro.write_text(
            gro_frame(
                DISTANCES_NM[0],
                0.0,
            ),
            encoding="utf-8",
        )

        write_trajectory(
            trajectory
        )

        grompp = run_command(
            [
                gmx,
                "grompp",
                "-f",
                str(mdp),
                "-c",
                str(first_gro),
                "-p",
                str(topology),
                "-n",
                str(index),
                "-o",
                str(tpr),
                "-po",
                str(mdout),
                "-pp",
                str(processed_top),
                "-maxwarn",
                "0",
            ],
            cwd=candidate_root,
        )

        grompp_log.write_text(
            grompp.stdout,
            encoding="utf-8",
        )

        if (
            grompp.returncode != 0
            or not tpr.exists()
        ):
            raise RuntimeError(
                f"grompp failed for {candidate_name}. "
                f"See {grompp_log}"
            )

        dump = subprocess.run(
            [
                gmx,
                "dump",
                "-s",
                str(tpr),
            ],
            cwd=candidate_root,
            env={
                **os.environ,
                "GMX_MAXBACKUP": "-1",
            },
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

        dump_stdout.write_text(
            dump.stdout,
            encoding="utf-8",
        )

        dump_stderr.write_text(
            dump.stderr,
            encoding="utf-8",
        )

        if dump.returncode != 0:
            raise RuntimeError(
                f"gmx dump failed for {candidate_name}"
            )

        (
            dump_contains_negative_sigma,
            dump_contains_nonbond_data,
        ) = parse_dump_nbfp(
            dump.stdout
        )

        rerun = run_command(
            [
                gmx,
                "mdrun",
                "-s",
                str(tpr),
                "-rerun",
                str(trajectory),
                "-deffnm",
                "rerun",
                "-ntmpi",
                "1",
                "-ntomp",
                "1",
            ],
            cwd=candidate_root,
        )

        rerun_log.write_text(
            rerun.stdout,
            encoding="utf-8",
        )

        edr = (
            candidate_root
            / "rerun.edr"
        )

        if (
            rerun.returncode != 0
            or not edr.exists()
        ):
            raise RuntimeError(
                f"Rerun failed for {candidate_name}. "
                f"See {rerun_log}"
            )

        (
            energy_array,
            selected_energy_term,
            energy_output,
        ) = extract_pair_energy(
            gmx,
            candidate_root,
            edr,
            energy_xvg,
        )

        extraction_log.write_text(
            energy_output,
            encoding="utf-8",
        )

        measured_energies = (
            energy_array[:, 1]
        )

        expected_energies = np.asarray(
            [
                c12_kj_mol_nm12
                / distance_nm
                ** 12
                for distance_nm in DISTANCES_NM
            ],
            dtype=float,
        )

        absolute_errors = np.abs(
            measured_energies
            - expected_energies
        )

        relative_errors = np.divide(
            absolute_errors,
            np.abs(
                expected_energies
            ),
            out=np.zeros_like(
                absolute_errors
            ),
            where=(
                np.abs(
                    expected_energies
                )
                > 1.0e-12
            ),
        )

        row_passes = []

        for frame_index, distance_nm in enumerate(
            DISTANCES_NM
        ):
            expected = float(
                expected_energies[
                    frame_index
                ]
            )

            measured = float(
                measured_energies[
                    frame_index
                ]
            )

            absolute_error = float(
                absolute_errors[
                    frame_index
                ]
            )

            relative_error = float(
                relative_errors[
                    frame_index
                ]
            )

            if abs(expected) >= 0.1:
                point_pass = (
                    relative_error
                    <= MAX_RELATIVE_ERROR
                )
            else:
                point_pass = (
                    absolute_error
                    <= MAX_ABSOLUTE_ERROR_KJ_MOL
                )

            row_passes.append(
                point_pass
            )

            validation_rows.append(
                {
                    "candidate": candidate_name,
                    "target_energy_kBT_at_0p17nm": (
                        target_level_kbt
                    ),
                    "pair_sigma_nm": (
                        -SIGMA_ABS_NM
                    ),
                    "pair_epsilon_kJ_mol": (
                        pair_epsilon_kj_mol
                    ),
                    "C6_kJ_mol_nm6": 0.0,
                    "C12_kJ_mol_nm12": (
                        c12_kj_mol_nm12
                    ),
                    "distance_nm": (
                        distance_nm
                    ),
                    "expected_energy_kJ_mol": (
                        expected
                    ),
                    "measured_energy_kJ_mol": (
                        measured
                    ),
                    "expected_energy_kBT": (
                        expected
                        / KBT_KJ_MOL
                    ),
                    "measured_energy_kBT": (
                        measured
                        / KBT_KJ_MOL
                    ),
                    "absolute_error_kJ_mol": (
                        absolute_error
                    ),
                    "relative_error": (
                        relative_error
                    ),
                    "point_pass": (
                        point_pass
                    ),
                }
            )

        candidate_pass = bool(
            grompp.returncode == 0
            and rerun.returncode == 0
            and all(
                row_passes
            )
            and np.all(
                measured_energies > 0.0
            )
        )

        all_candidates_pass = (
            all_candidates_pass
            and candidate_pass
        )

        summary_rows.append(
            {
                "candidate": candidate_name,
                "target_energy_kBT_at_0p17nm": (
                    target_level_kbt
                ),
                "pair_sigma_nm": (
                    -SIGMA_ABS_NM
                ),
                "pair_epsilon_kJ_mol": (
                    pair_epsilon_kj_mol
                ),
                "C6_kJ_mol_nm6": 0.0,
                "C12_kJ_mol_nm12": (
                    c12_kj_mol_nm12
                ),
                "grompp_return_code": (
                    grompp.returncode
                ),
                "rerun_return_code": (
                    rerun.returncode
                ),
                "selected_energy_term": (
                    selected_energy_term
                ),
                "all_measured_energies_positive": (
                    bool(
                        np.all(
                            measured_energies > 0.0
                        )
                    )
                ),
                "maximum_absolute_error_kJ_mol": (
                    float(
                        np.max(
                            absolute_errors
                        )
                    )
                ),
                "maximum_relative_error_for_energy_ge_0p1": (
                    float(
                        np.max(
                            relative_errors[
                                expected_energies
                                >= 0.1
                            ]
                        )
                    )
                    if np.any(
                        expected_energies
                        >= 0.1
                    )
                    else 0.0
                ),
                "dump_contains_negative_sigma_text": (
                    dump_contains_negative_sigma
                ),
                "dump_contains_nonbond_data": (
                    dump_contains_nonbond_data
                ),
                "candidate_pass": (
                    candidate_pass
                ),
            }
        )

    write_csv(
        VALIDATION_CSV,
        validation_rows,
    )

    write_csv(
        CANDIDATE_SUMMARY_CSV,
        summary_rows,
    )

    decision = (
        "PURE_R12_NEGATIVE_SIGMA_OVERRIDE_VALIDATED"
        if all_candidates_pass
        else
        "NEGATIVE_SIGMA_OVERRIDE_REQUIRES_REVIEW"
    )

    required_next_step = (
        "BUILD_FULL_R1_TOPOLOGY_CANDIDATES_AND_RUN_"
        "STATIC_CAP_WATER_ENERGY_SCAN"
        if all_candidates_pass
        else
        "REVIEW_NEGATIVE_SIGMA_MICROVALIDATION"
    )

    corrected_contract = {
        **previous_contract,
        "previous_decision": (
            previous_contract[
                "decision"
            ]
        ),
        "decision": decision,
        "parameter_semantics": (
            "SIGMA_EPSILON_WITH_NEGATIVE_SIGMA_"
            "SPECIAL_CASE"
        ),
        "negative_sigma_special_case": True,
        "proposed_cap_water_O_interaction": (
            "PURE_R12_USING_NEGATIVE_SIGMA_OVERRIDE"
        ),
        "pair_sigma_nm": (
            -SIGMA_ABS_NM
        ),
        "candidate_target_levels_kBT_at_0p17nm": list(
            TARGET_LEVELS_KBT
        ),
        "all_microvalidation_candidates_pass": (
            all_candidates_pass
        ),
        "pure_r12_standard_nonbond_override_feasible": (
            all_candidates_pass
        ),
        "tabulated_interaction_required": False,
        "full_R1_topology_built": False,
        "energy_minimization_authorized": False,
        "MD_execution_authorized": False,
        "required_next_step": (
            required_next_step
        ),
        "supersedes_contract": relative(
            PREVIOUS_CONTRACT
        ),
    }

    CORRECTED_CONTRACT_JSON.write_text(
        json.dumps(
            corrected_contract,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    candidate_lines = "\n".join(
        (
            f"- {row['candidate']}: "
            f"epsilon={row['pair_epsilon_kJ_mol']:.8f} "
            f"kJ mol^-1; "
            f"C12={row['C12_kJ_mol_nm12']:.12e} "
            f"kJ mol^-1 nm^12; "
            f"max relative error="
            f"{row['maximum_relative_error_for_energy_ge_0p1']:.6e}; "
            f"**{'PASS' if row['candidate_pass'] else 'FAIL'}**"
        )
        for row in summary_rows
    )

    REPORT_MD.write_text(
        f"""# R1 Negative-Sigma Pure-r12 Validation

## Correction to the previous topology audit

The previous audit treated `comb-rule = 2` as incompatible with a
pure C12 interaction. GROMACS provides a documented special case:
when sigma is negative, C6 is set to zero and C12 is calculated from
the absolute sigma value and epsilon.

Therefore, a CAP–OW pure repulsive interaction can be represented as:

- sigma = **{-SIGMA_ABS_NM:.6f} nm**
- epsilon > 0
- C6 = **0**
- C12 = **4 epsilon |sigma|^12**

The previous contract is retained for provenance but is superseded by:

`{relative(CORRECTED_CONTRACT_JSON)}`

## Microvalidation setup

- GROMACS executable: `{gmx}`
- Temperature reference: **{TEMPERATURE_K:.1f} K**
- kBT: **{KBT_KJ_MOL:.8f} kJ mol^-1**
- Baseline OW sigma:
  **{ow_sigma_nm:.10f} nm**
- Baseline OW epsilon:
  **{ow_epsilon_kj_mol:.10f} kJ mol^-1**
- Negative pair sigma:
  **{-SIGMA_ABS_NM:.6f} nm**
- Distances:
  **{', '.join(f'{value:.3f}' for value in DISTANCES_NM)} nm**

Each candidate was evaluated using `grompp`, `mdrun -rerun`, and the
CAP–OW short-range Lennard-Jones energy extracted from the resulting
energy file.

## Candidate results

{candidate_lines}

## Decision

- Decision: **{decision}**
- All candidates passed:
  **{'YES' if all_candidates_pass else 'NO'}**
- Pure-r12 standard nonbonded override feasible:
  **{'YES' if all_candidates_pass else 'NO'}**
- Tabulated interaction required: **NO**
- Full R1 topology built: **NO**
- Energy minimization authorized: **NO**
- MD execution authorized: **NO**
- Required next step:
  `{required_next_step}`

The next gate must evaluate the candidate interaction strengths in the
complete 68,314-atom R1 system before selecting a production value.
""",
        encoding="utf-8",
    )

    print(
        "Day023 R1 negative-sigma pure-r12 validation completed."
    )

    print(
        "Baseline OW sigma / epsilon: "
        f"{ow_sigma_nm:.10f} nm / "
        f"{ow_epsilon_kj_mol:.10f} kJ/mol"
    )

    print(
        "Negative pair sigma: "
        f"{-SIGMA_ABS_NM:.6f} nm"
    )

    print(
        "Candidates tested / passed: "
        f"{len(summary_rows)}/"
        f"{sum(bool(row['candidate_pass']) for row in summary_rows)}"
    )

    for row in summary_rows:
        print(
            "Candidate "
            f"{row['target_energy_kBT_at_0p17nm']:.0f} kBT: "
            f"epsilon={row['pair_epsilon_kJ_mol']:.8f} kJ/mol; "
            f"C12={row['C12_kJ_mol_nm12']:.12e}; "
            "max relative error="
            f"{row['maximum_relative_error_for_energy_ge_0p1']:.6e}; "
            f"{'PASS' if row['candidate_pass'] else 'FAIL'}"
        )

    print(
        f"Decision: {decision}"
    )

    print(
        "Pure r^-12 standard override feasible: "
        f"{'YES' if all_candidates_pass else 'NO'}"
    )

    print(
        "Tabulated interaction required: NO"
    )

    print(
        "Full R1 topology built: NO"
    )

    print(
        "Energy minimization authorized: NO"
    )

    print(
        "MD execution authorized: NO"
    )

    print(
        f"Required next step: {required_next_step}"
    )

    print(
        f"Wrote: {relative(VALIDATION_CSV)}"
    )

    print(
        f"Wrote: {relative(CANDIDATE_SUMMARY_CSV)}"
    )

    print(
        f"Wrote: {relative(CORRECTED_CONTRACT_JSON)}"
    )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )

    if not all_candidates_pass:
        raise RuntimeError(
            "Negative-sigma pure-r12 validation "
            "requires review."
        )


if __name__ == "__main__":
    main()
