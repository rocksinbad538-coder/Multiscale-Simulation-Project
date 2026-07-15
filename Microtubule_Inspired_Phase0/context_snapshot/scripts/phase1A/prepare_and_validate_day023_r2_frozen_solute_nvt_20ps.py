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
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

DAY023_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design"
)

EM_ROOT = (
    DAY023_ROOT
    / "14_r2_water_only_em"
)

EM_SELECTED_ROOT = (
    EM_ROOT
    / "selected"
)

OUTPUT_ROOT = (
    DAY023_ROOT
    / "15_r2_frozen_solute_nvt_20ps_preparation"
)

SOURCE_GRO = (
    EM_ROOT
    / "r2_water_only_em.gro"
)

SOURCE_TOPOLOGY = (
    EM_SELECTED_ROOT
    / "r2_water_only_em.top"
)

SOURCE_INDEX = (
    EM_SELECTED_ROOT
    / "r2_water_only_em.ndx"
)

SOURCE_SUMMARY = (
    EM_ROOT
    / "r2_water_only_em_summary.csv"
)

LOCAL_GRO = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_input.gro"
)

LOCAL_TOPOLOGY = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps.top"
)

LOCAL_INDEX = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps.ndx"
)

NVT_MDP = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps.mdp"
)

PROCESSED_MDP = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_processed.mdp"
)

PROBE_TPR = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_probe.tpr"
)

FINAL_TPR = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps.tpr"
)

GROMPP_PROBE_LOG = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_grompp_probe.log"
)

GROMPP_AUTHORIZED_LOG = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_grompp_authorized.log"
)

TPR_DUMP = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_tpr_dump.txt"
)

TPR_DUMP_STDERR = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_tpr_dump.stderr.log"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_preparation_summary.csv"
)

GATES_CSV = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_preparation_gates.csv"
)

CHECKSUMS_CSV = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_checksums.csv"
)

RUN_CONTRACT_JSON = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20ps_run_contract.json"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_FROZEN_SOLUTE_NVT_20PS_PREPARATION_DAY023.md"
)

EXPECTED_EM_DECISION = (
    "R2_WATER_ONLY_ENERGY_MINIMIZATION_VALIDATED"
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
TOTAL_TIME_PS = DT_PS * NSTEPS

NSTXOUT_COMPRESSED = 1000
NSTENERGY = 100
NSTLOG = 100

TRAJECTORY_INTERVAL_PS = (
    DT_PS
    * NSTXOUT_COMPRESSED
)

EXPECTED_XTC_FRAMES = (
    NSTEPS
    // NSTXOUT_COMPRESSED
    + 1
)

TEMPERATURE_K = 300.0
TAU_T_PS = 0.1
GEN_SEED = 20260708

EXPECTED_SOL_DOF = (
    6 * WATERS
    - 3
)

EXPECTED_WARNING_MAXIMUM = 1

EXPECTED_WARNING_REQUIRED_WORDS = (
    "verlet",
    "frozen",
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
                    key: row.get(key, "")
                    for key in fields
                }
            )


def read_gro(
    path: Path,
) -> tuple[
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
            2 + atom_count
        ].split()
    ]

    if len(box_values) != 3:
        raise RuntimeError(
            "An orthorhombic box is required."
        )

    box = np.asarray(
        box_values,
        dtype=float,
    )

    return atoms, box


def parse_index_groups(
    path: Path,
) -> dict[str, list[int]]:
    require_file(path)

    groups: dict[str, list[int]] = {}
    current: str | None = None

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        match = re.match(
            r"^\s*\[\s*([^\]]+?)\s*\]\s*$",
            raw_line,
        )

        if match is not None:
            current = match.group(1).strip()
            groups[current] = []
            continue

        if current is None:
            continue

        stripped = raw_line.strip()

        if not stripped:
            continue

        groups[current].extend(
            int(token)
            for token in stripped.split()
        )

    return groups


def parse_topology_molecule_counts(
    path: Path,
) -> dict[str, int]:
    lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    in_molecules = False
    counts: dict[str, int] = {}

    for raw_line in lines:
        section_match = re.match(
            r"^\s*\[\s*([^\]]+?)\s*\]\s*$",
            raw_line,
        )

        if section_match is not None:
            in_molecules = (
                section_match.group(1)
                .strip()
                .lower()
                == "molecules"
            )
            continue

        if not in_molecules:
            continue

        stripped = (
            raw_line
            .split(";", 1)[0]
            .strip()
        )

        if not stripped:
            continue

        fields = stripped.split()

        if len(fields) < 2:
            continue

        try:
            counts[fields[0]] = int(
                fields[1]
            )
        except ValueError:
            continue

    return counts


def warning_count(text: str) -> int:
    return len(
        re.findall(
            r"^\s*WARNING\s+\d+",
            text,
            flags=re.MULTILINE,
        )
    )


def expected_freeze_warning(
    text: str,
) -> bool:
    lowered = text.lower()

    return all(
        word in lowered
        for word in EXPECTED_WARNING_REQUIRED_WORDS
    )


def parse_scalar(
    dump_text: str,
    name: str,
) -> str:
    match = re.search(
        rf"^\s*{re.escape(name)}\s*=\s*(\S+)",
        dump_text,
        flags=re.MULTILINE,
    )

    if match is None:
        raise RuntimeError(
            f"Could not parse {name!r} from TPR dump."
        )

    return match.group(1)


def parse_tpr_velocities(
    dump_text: str,
    atom_count: int,
) -> np.ndarray:
    velocities = np.full(
        (
            atom_count,
            3,
        ),
        np.nan,
        dtype=float,
    )

    pattern = re.compile(
        r"v\[\s*(\d+)\s*\]\s*=\s*"
        r"\{\s*([-+0-9.eE]+)\s*,\s*"
        r"([-+0-9.eE]+)\s*,\s*"
        r"([-+0-9.eE]+)\s*\}"
    )

    for match in pattern.finditer(
        dump_text
    ):
        index = int(
            match.group(1)
        )

        if (
            index < 0
            or index >= atom_count
        ):
            raise RuntimeError(
                "Velocity index outside atom range."
            )

        velocities[index] = [
            float(match.group(2)),
            float(match.group(3)),
            float(match.group(4)),
        ]

    represented = int(
        np.count_nonzero(
            np.all(
                np.isfinite(velocities),
                axis=1,
            )
        )
    )

    if represented != atom_count:
        raise RuntimeError(
            "TPR velocity table is incomplete: "
            f"{represented}/{atom_count}"
        )

    return velocities


def parse_dof(
    text: str,
) -> dict[str, float]:
    result: dict[str, float] = {}

    pattern = re.compile(
        r"Number of degrees of freedom in "
        r"T-Coupling group\s+(\S+)\s+is\s+"
        r"([-+0-9.eE]+)",
        flags=re.IGNORECASE,
    )

    for match in pattern.finditer(text):
        result[
            match.group(1)
        ] = float(
            match.group(2)
        )

    return result


def parse_processed_mdp(
    path: Path,
) -> dict[str, str]:
    require_file(path)

    values: dict[str, str] = {}

    for raw_line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        line = (
            raw_line
            .split(";", 1)[0]
            .strip()
        )

        if (
            not line
            or "=" not in line
        ):
            continue

        key, value = line.split(
            "=",
            1,
        )

        values[
            key.strip().lower()
        ] = value.strip()

    return values


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        SOURCE_GRO,
        SOURCE_TOPOLOGY,
        SOURCE_INDEX,
        SOURCE_SUMMARY,
    ):
        require_file(required)

    summary = read_single_csv_row(
        SOURCE_SUMMARY
    )

    if (
        summary.get("decision")
        != EXPECTED_EM_DECISION
    ):
        raise RuntimeError(
            "R2 water-only minimization is not accepted."
        )

    if not parse_bool(
        summary.get(
            "short_frozen_solute_NVT_preparation_authorized",
            "false",
        )
    ):
        raise RuntimeError(
            "R2 water-only minimization did not "
            "authorize NVT preparation."
        )

    existing_outputs = [
        path
        for path in (
            FINAL_TPR,
            PROBE_TPR,
            SUMMARY_CSV,
            GATES_CSV,
            RUN_CONTRACT_JSON,
        )
        if path.exists()
    ]

    if existing_outputs:
        raise RuntimeError(
            "R2 NVT preparation products already exist: "
            + ", ".join(
                path.name
                for path in existing_outputs
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

    atoms, box = read_gro(
        LOCAL_GRO
    )

    if len(atoms) != EXPECTED_ATOMS:
        raise RuntimeError(
            "Unexpected R2 atom count: "
            f"{len(atoms)}/{EXPECTED_ATOMS}"
        )

    if CAP_STOP != EXPECTED_ATOMS:
        raise RuntimeError(
            "Internal atom accounting failed."
        )

    molecule_counts = (
        parse_topology_molecule_counts(
            LOCAL_TOPOLOGY
        )
    )

    index_groups = parse_index_groups(
        LOCAL_INDEX
    )

    expected_group_counts = {
        "System": EXPECTED_ATOMS,
        "HBN": HBN_ATOMS,
        "PYR": PYR_ATOMS,
        "HBN_PYR": SOLUTE_ATOMS,
        "SOL": WATER_ATOMS,
        "CAPL": CAPS_PER_END,
        "CAPU": CAPS_PER_END,
        "CAPS": TOTAL_CAPS,
    }

    NVT_MDP.write_text(
        f"""integrator               = md
dt                       = {DT_PS:.7f}
nsteps                   = {NSTEPS}
continuation             = no

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

tcoupl                   = V-rescale
tc-grps                  = SOL HBN_PYR CAPS
tau-t                    = {TAU_T_PS:.3f} -1 -1
ref-t                    = {TEMPERATURE_K:.1f} {TEMPERATURE_K:.1f} {TEMPERATURE_K:.1f}

pcoupl                   = no

comm-mode                = Linear
nstcomm                  = 100
comm-grps                = SOL

freezegrps               = HBN_PYR CAPS
freezedim                = Y Y Y Y Y Y

gen-vel                  = yes
gen-temp                 = {TEMPERATURE_K:.1f}
gen-seed                 = {GEN_SEED}

nstlog                   = {NSTLOG}
nstenergy                = {NSTENERGY}
nstxout                   = 0
nstvout                  = 0
nstfout                  = 0
nstxout-compressed       = {NSTXOUT_COMPRESSED}
compressed-x-precision   = 1000

energygrps               = HBN_PYR CAPS SOL
""",
        encoding="utf-8",
    )

    gmx = locate_gmx()

    probe = run_command(
        [
            gmx,
            "grompp",
            "-f",
            str(NVT_MDP),
            "-c",
            str(LOCAL_GRO),
            "-p",
            str(LOCAL_TOPOLOGY),
            "-n",
            str(LOCAL_INDEX),
            "-o",
            str(PROBE_TPR),
            "-po",
            str(PROCESSED_MDP),
            "-maxwarn",
            "0",
        ],
        cwd=OUTPUT_ROOT,
    )

    GROMPP_PROBE_LOG.write_text(
        probe.stdout,
        encoding="utf-8",
    )

    probe_warning_count = warning_count(
        probe.stdout
    )

    warning_is_authorized = (
        probe_warning_count == 1
        and expected_freeze_warning(
            probe.stdout
        )
    )

    controlled_maxwarn_used = False

    if probe.returncode == 0:
        if probe_warning_count != 0:
            raise RuntimeError(
                "Grompp succeeded but warning accounting "
                "was inconsistent."
            )

        require_file(PROBE_TPR)

        shutil.move(
            str(PROBE_TPR),
            str(FINAL_TPR),
        )

        GROMPP_AUTHORIZED_LOG.write_text(
            probe.stdout,
            encoding="utf-8",
        )

        final_grompp_return_code = 0
        final_warning_count = 0

    else:
        if not warning_is_authorized:
            raise RuntimeError(
                "R2 NVT grompp failed for a reason other "
                "than the single expected frozen-particle "
                "Verlet-buffer warning. See "
                f"{GROMPP_PROBE_LOG}"
            )

        authorized = run_command(
            [
                gmx,
                "grompp",
                "-f",
                str(NVT_MDP),
                "-c",
                str(LOCAL_GRO),
                "-p",
                str(LOCAL_TOPOLOGY),
                "-n",
                str(LOCAL_INDEX),
                "-o",
                str(FINAL_TPR),
                "-po",
                str(PROCESSED_MDP),
                "-maxwarn",
                "1",
            ],
            cwd=OUTPUT_ROOT,
        )

        GROMPP_AUTHORIZED_LOG.write_text(
            authorized.stdout,
            encoding="utf-8",
        )

        final_grompp_return_code = (
            authorized.returncode
        )

        final_warning_count = warning_count(
            authorized.stdout
        )

        controlled_maxwarn_used = True

        if (
            authorized.returncode != 0
            or final_warning_count != 1
            or not expected_freeze_warning(
                authorized.stdout
            )
        ):
            raise RuntimeError(
                "Controlled grompp did not reproduce "
                "exactly the authorized warning."
            )

    require_file(FINAL_TPR)
    require_file(PROCESSED_MDP)

    dump = subprocess.run(
        [
            gmx,
            "dump",
            "-s",
            str(FINAL_TPR),
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

    TPR_DUMP.write_text(
        dump.stdout,
        encoding="utf-8",
    )

    TPR_DUMP_STDERR.write_text(
        dump.stderr,
        encoding="utf-8",
    )

    if dump.returncode != 0:
        raise RuntimeError(
            "Could not dump the prepared R2 TPR."
        )

    dump_text = dump.stdout

    tpr_atoms = int(
        parse_scalar(
            dump_text,
            "natoms",
        )
    )

    tpr_integrator = parse_scalar(
        dump_text,
        "integrator",
    )

    tpr_dt = float(
        parse_scalar(
            dump_text,
            "dt",
        )
    )

    tpr_nsteps = int(
        parse_scalar(
            dump_text,
            "nsteps",
        )
    )

    tpr_nstlog = int(
        parse_scalar(
            dump_text,
            "nstlog",
        )
    )

    tpr_nstenergy = int(
        parse_scalar(
            dump_text,
            "nstenergy",
        )
    )

    tpr_nstxout_compressed = int(
        parse_scalar(
            dump_text,
            "nstxout-compressed",
        )
    )

    tpr_continuation = parse_scalar(
        dump_text,
        "continuation",
    ).lower()

    velocities = parse_tpr_velocities(
        dump_text,
        EXPECTED_ATOMS,
    )

    oxygen_indices = np.asarray(
        [
            atom["index"]
            for atom in atoms
            if (
                atom["resname"] == "SOL"
                and atom["atomname"] == "OW"
            )
        ],
        dtype=int,
    )

    if len(oxygen_indices) != WATERS:
        raise RuntimeError(
            "Unexpected number of water oxygens: "
            f"{len(oxygen_indices)}/{WATERS}"
        )

    water_ow_speeds = np.linalg.norm(
        velocities[
            oxygen_indices
        ],
        axis=1,
    )

    hbn_pyr_speeds = np.linalg.norm(
        velocities[
            :SOLUTE_ATOMS
        ],
        axis=1,
    )

    caps_speeds = np.linalg.norm(
        velocities[
            CAP_START:
            CAP_STOP
        ],
        axis=1,
    )

    processed = parse_processed_mdp(
        PROCESSED_MDP
    )

    grompp_text = (
        GROMPP_AUTHORIZED_LOG.read_text(
            encoding="utf-8",
            errors="replace",
        )
    )

    dof = parse_dof(
        grompp_text
    )

    sol_dof = dof.get(
        "SOL",
        math.nan,
    )

    hbn_pyr_dof = dof.get(
        "HBN_PYR",
        math.nan,
    )

    caps_dof = dof.get(
        "CAPS",
        math.nan,
    )

    index_count_gates = {
        name: (
            name in index_groups
            and len(
                index_groups[name]
            )
            == expected_count
        )
        for name, expected_count
        in expected_group_counts.items()
    }

    molecule_count_gates = {
        "HBN": (
            molecule_counts.get("HBN")
            == 1
        ),
        "PYR": (
            molecule_counts.get("PYR")
            == 4
        ),
        "SOL": (
            molecule_counts.get("SOL")
            == WATERS
        ),
        "CAPL": (
            molecule_counts.get("CAPL")
            == 1
        ),
        "CAPU": (
            molecule_counts.get("CAPU")
            == 1
        ),
    }

    processed_mdp_gates = {
        "processed_integrator_is_md": (
            processed.get(
                "integrator",
                "",
            ).lower()
            == "md"
        ),
        "processed_nsteps_is_40000": (
            int(
                processed.get(
                    "nsteps",
                    "-1",
                )
            )
            == NSTEPS
        ),
        "processed_continuation_is_no": (
            processed.get(
                "continuation",
                "",
            ).lower()
            in {
                "no",
                "false",
            }
        ),
        "processed_gen_vel_is_yes": (
            processed.get(
                "gen-vel",
                "",
            ).lower()
            in {
                "yes",
                "true",
            }
        ),
        "processed_gen_seed_is_20260708": (
            int(
                processed.get(
                    "gen-seed",
                    "-1",
                )
            )
            == GEN_SEED
        ),
        "processed_freezegrps_are_correct": (
            processed.get(
                "freezegrps",
                "",
            ).split()
            == [
                "HBN_PYR",
                "CAPS",
            ]
        ),
        "processed_freezedim_is_correct": (
            processed.get(
                "freezedim",
                "",
            ).split()
            == [
                "Y",
                "Y",
                "Y",
                "Y",
                "Y",
                "Y",
            ]
        ),
        "processed_tc_grps_are_correct": (
            processed.get(
                "tc-grps",
                "",
            ).split()
            == [
                "SOL",
                "HBN_PYR",
                "CAPS",
            ]
        ),
    }

    gates: dict[str, bool] = {
        "R2_water_only_EM_is_validated": (
            summary.get("decision")
            == EXPECTED_EM_DECISION
        ),
        "R2_EM_authorized_NVT_preparation": (
            parse_bool(
                summary.get(
                    "short_frozen_solute_NVT_preparation_authorized",
                    "false",
                )
            )
        ),
        "source_GRO_has_68332_atoms": (
            len(atoms)
            == EXPECTED_ATOMS
        ),
        "source_topology_molecule_counts_are_correct": (
            all(
                molecule_count_gates.values()
            )
        ),
        "source_index_group_counts_are_correct": (
            all(
                index_count_gates.values()
            )
        ),
        "grompp_return_code_zero": (
            final_grompp_return_code == 0
        ),
        "grompp_warning_policy_is_valid": (
            (
                final_warning_count == 0
                and not controlled_maxwarn_used
            )
            or (
                final_warning_count == 1
                and controlled_maxwarn_used
                and expected_freeze_warning(
                    grompp_text
                )
            )
        ),
        "TPR_dump_return_code_zero": (
            dump.returncode == 0
        ),
        "TPR_has_68332_atoms": (
            tpr_atoms
            == EXPECTED_ATOMS
        ),
        "TPR_integrator_is_md": (
            tpr_integrator.lower()
            == "md"
        ),
        "TPR_dt_is_0p0005ps": (
            math.isclose(
                tpr_dt,
                DT_PS,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
        ),
        "TPR_nsteps_is_40000": (
            tpr_nsteps
            == NSTEPS
        ),
        "TPR_total_time_is_20ps": (
            math.isclose(
                tpr_dt
                * tpr_nsteps,
                TOTAL_TIME_PS,
                rel_tol=0.0,
                abs_tol=1.0e-10,
            )
        ),
        "TPR_continuation_is_false": (
            tpr_continuation
            in {
                "false",
                "no",
            }
        ),
        "TPR_nstlog_is_100": (
            tpr_nstlog
            == NSTLOG
        ),
        "TPR_nstenergy_is_100": (
            tpr_nstenergy
            == NSTENERGY
        ),
        "TPR_XTC_stride_is_1000": (
            tpr_nstxout_compressed
            == NSTXOUT_COMPRESSED
        ),
        "TPR_has_all_velocity_entries": (
            np.all(
                np.isfinite(
                    velocities
                )
            )
        ),
        "water_oxygen_velocities_are_nonzero": (
            float(
                np.sqrt(
                    np.mean(
                        water_ow_speeds
                        * water_ow_speeds
                    )
                )
            )
            > 0.01
        ),
        "SOL_degrees_of_freedom_are_expected": (
            math.isfinite(sol_dof)
            and math.isclose(
                sol_dof,
                EXPECTED_SOL_DOF,
                rel_tol=0.0,
                abs_tol=0.01,
            )
        ),
        "HBN_PYR_degrees_of_freedom_are_zero": (
            math.isfinite(hbn_pyr_dof)
            and math.isclose(
                hbn_pyr_dof,
                0.0,
                rel_tol=0.0,
                abs_tol=0.01,
            )
        ),
        "CAPS_degrees_of_freedom_are_zero": (
            math.isfinite(caps_dof)
            and math.isclose(
                caps_dof,
                0.0,
                rel_tol=0.0,
                abs_tol=0.01,
            )
        ),
        "processed_MDP_contract_is_correct": (
            all(
                processed_mdp_gates.values()
            )
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
        "R2_FROZEN_SOLUTE_NVT_20PS_PREPARED"
        if accepted
        else
        "R2_FROZEN_SOLUTE_NVT_20PS_PREPARATION_REQUIRES_REVIEW"
    )

    required_next_step = (
        "RUN_R2_FROZEN_SOLUTE_NVT_20PS"
        if accepted
        else
        "REVIEW_R2_FROZEN_SOLUTE_NVT_PREPARATION_FAILURES"
    )

    input_files = (
        LOCAL_GRO,
        LOCAL_TOPOLOGY,
        LOCAL_INDEX,
        NVT_MDP,
        FINAL_TPR,
    )

    checksum_rows = [
        {
            "file": relative(path),
            "sha256": sha256(path),
            "size_bytes": path.stat().st_size,
        }
        for path in input_files
    ]

    write_csv(
        CHECKSUMS_CSV,
        checksum_rows,
    )

    summary_row = {
        "decision": decision,
        "atoms": EXPECTED_ATOMS,
        "waters": WATERS,
        "CAPL_atoms": CAPS_PER_END,
        "CAPU_atoms": CAPS_PER_END,
        "dt_ps": DT_PS,
        "nsteps": NSTEPS,
        "total_time_ps": TOTAL_TIME_PS,
        "nstxout_compressed": (
            NSTXOUT_COMPRESSED
        ),
        "trajectory_interval_ps": (
            TRAJECTORY_INTERVAL_PS
        ),
        "expected_XTC_frames": (
            EXPECTED_XTC_FRAMES
        ),
        "nstenergy": NSTENERGY,
        "nstlog": NSTLOG,
        "temperature_K": (
            TEMPERATURE_K
        ),
        "tau_t_SOL_ps": (
            TAU_T_PS
        ),
        "generation_seed": (
            GEN_SEED
        ),
        "grompp_probe_return_code": (
            probe.returncode
        ),
        "grompp_final_return_code": (
            final_grompp_return_code
        ),
        "grompp_warning_count": (
            final_warning_count
        ),
        "controlled_maxwarn_used": (
            controlled_maxwarn_used
        ),
        "SOL_degrees_of_freedom": (
            sol_dof
        ),
        "expected_SOL_degrees_of_freedom": (
            EXPECTED_SOL_DOF
        ),
        "HBN_PYR_degrees_of_freedom": (
            hbn_pyr_dof
        ),
        "CAPS_degrees_of_freedom": (
            caps_dof
        ),
        "water_O_RMS_speed_nm_ps": float(
            np.sqrt(
                np.mean(
                    water_ow_speeds
                    * water_ow_speeds
                )
            )
        ),
        "water_O_max_speed_nm_ps": float(
            np.max(
                water_ow_speeds
            )
        ),
        "HBN_PYR_RMS_input_speed_nm_ps": float(
            np.sqrt(
                np.mean(
                    hbn_pyr_speeds
                    * hbn_pyr_speeds
                )
            )
        ),
        "CAPS_RMS_input_speed_nm_ps": float(
            np.sqrt(
                np.mean(
                    caps_speeds
                    * caps_speeds
                )
            )
        ),
        "failed_gates": (
            " | ".join(
                failed_gates
            )
        ),
        "NVT_execution_authorized": (
            accepted
        ),
        "velocity_regeneration_during_mdrun_authorized": False,
        "long_mobile_MD_authorized": False,
        "multitemperature_MD_authorized": False,
        "QM_recalculation_authorized": False,
        "required_next_step": (
            required_next_step
        ),
    }

    write_csv(
        SUMMARY_CSV,
        [summary_row],
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

    contract = {
        "decision": decision,
        "execution_authorized": (
            accepted
        ),
        "gmx_command": [
            "gmx",
            "mdrun",
            "-s",
            relative(FINAL_TPR),
            "-deffnm",
            (
                "runs/phase1A/"
                "day023_confinement_design/"
                "16_r2_frozen_solute_nvt_20ps/"
                "r2_frozen_solute_nvt_20ps"
            ),
            "-ntmpi",
            "1",
            "-ntomp",
            "4",
        ],
        "input_TPR": relative(
            FINAL_TPR
        ),
        "atoms": EXPECTED_ATOMS,
        "waters": WATERS,
        "dt_ps": DT_PS,
        "nsteps": NSTEPS,
        "target_time_ps": TOTAL_TIME_PS,
        "expected_XTC_frames": (
            EXPECTED_XTC_FRAMES
        ),
        "expected_XTC_interval_ps": (
            TRAJECTORY_INTERVAL_PS
        ),
        "temperature_K": TEMPERATURE_K,
        "frozen_groups": [
            "HBN_PYR",
            "CAPS",
        ],
        "thermostatted_group": "SOL",
        "velocity_generation_performed_by_grompp": True,
        "velocity_generation_seed": GEN_SEED,
        "velocity_regeneration_during_mdrun": False,
        "checkpoint_continuation": False,
        "long_mobile_MD_authorized": False,
        "multitemperature_MD_authorized": False,
        "QM_recalculation_authorized": False,
        "required_next_step": (
            required_next_step
        ),
    }

    RUN_CONTRACT_JSON.write_text(
        json.dumps(
            contract,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed in gates.items()
    )

    REPORT_MD.write_text(
        f"""# R2 Frozen-Solute NVT 20 ps Preparation

## Scope

This stage prepares, but does not execute, a 20 ps frozen-solute NVT
screen for R2.

The input state is the validated R2 water-only minimized structure.
HBN, all four pyrenes, and both partial-cap assemblies will remain
frozen. Only TIP4P/2005 water has nonzero dynamical degrees of freedom.

## Protocol

- Atoms:
  **{EXPECTED_ATOMS}**
- Waters:
  **{WATERS}**
- CAPL/CAPU:
  **{CAPS_PER_END}/{CAPS_PER_END} beads**
- Integrator:
  **md**
- Time step:
  **{DT_PS:.7f} ps**
- Steps:
  **{NSTEPS}**
- Duration:
  **{TOTAL_TIME_PS:.3f} ps**
- Temperature:
  **{TEMPERATURE_K:.1f} K**
- Thermostat:
  **V-rescale**
- Effective thermostat group:
  **SOL**
- Frozen groups:
  **HBN_PYR, CAPS**
- Velocity seed:
  **{GEN_SEED}**
- XTC interval:
  **{TRAJECTORY_INTERVAL_PS:.3f} ps**
- Expected XTC frames:
  **{EXPECTED_XTC_FRAMES}**

## Grompp audit

- Probe return code:
  **{probe.returncode}**
- Final return code:
  **{final_grompp_return_code}**
- Warning count:
  **{final_warning_count}**
- Controlled `-maxwarn 1` used:
  **{'YES' if controlled_maxwarn_used else 'NO'}**
- SOL degrees of freedom:
  **{sol_dof:.3f}**
- Expected SOL degrees of freedom:
  **{EXPECTED_SOL_DOF}**
- HBN_PYR degrees of freedom:
  **{hbn_pyr_dof:.3f}**
- CAPS degrees of freedom:
  **{caps_dof:.3f}**

Any accepted warning is restricted to the known GROMACS
Verlet-buffer/frozen-particle warning. No unrelated warning is
authorized.

## Generated velocities

- OW RMS speed:
  **{np.sqrt(np.mean(water_ow_speeds * water_ow_speeds)):.6f} nm/ps**
- OW maximum speed:
  **{np.max(water_ow_speeds):.6f} nm/ps**

Velocity generation has already been encoded in the TPR. Mdrun must
not regenerate velocities.

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- NVT execution authorized:
  **{'YES' if accepted else 'NO'}**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `{required_next_step}`

R2 remains a frozen neutral steric screening design. Preparation of
this short trajectory does not establish chemical realizability,
long-time retention, or a final device architecture.
""",
        encoding="utf-8",
    )

    print(
        "Day023 R2 frozen-solute NVT 20 ps "
        "preparation completed."
    )

    print(
        "Atoms / waters / CAPL / CAPU: "
        f"{EXPECTED_ATOMS}/"
        f"{WATERS}/"
        f"{CAPS_PER_END}/"
        f"{CAPS_PER_END}"
    )

    print(
        "Grompp probe / final return codes: "
        f"{probe.returncode}/"
        f"{final_grompp_return_code}"
    )

    print(
        "Grompp warnings / controlled maxwarn used: "
        f"{final_warning_count}/"
        f"{'YES' if controlled_maxwarn_used else 'NO'}"
    )

    print(
        "TPR dt / nsteps / duration: "
        f"{tpr_dt:.7f} ps / "
        f"{tpr_nsteps} / "
        f"{tpr_dt * tpr_nsteps:.3f} ps"
    )

    print(
        "XTC stride / interval / expected frames: "
        f"{tpr_nstxout_compressed}/"
        f"{TRAJECTORY_INTERVAL_PS:.3f} ps/"
        f"{EXPECTED_XTC_FRAMES}"
    )

    print(
        "SOL / HBN_PYR / CAPS degrees of freedom: "
        f"{sol_dof:.3f}/"
        f"{hbn_pyr_dof:.3f}/"
        f"{caps_dof:.3f}"
    )

    print(
        "OW RMS / maximum generated speed: "
        f"{np.sqrt(np.mean(water_ow_speeds * water_ow_speeds)):.6f}/"
        f"{np.max(water_ow_speeds):.6f} nm/ps"
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
        f"{'YES' if accepted else 'NO'}"
    )

    print(
        "Velocity regeneration during mdrun authorized: NO"
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
        LOCAL_GRO,
        LOCAL_TOPOLOGY,
        LOCAL_INDEX,
        NVT_MDP,
        PROCESSED_MDP,
        FINAL_TPR,
        TPR_DUMP,
        CHECKSUMS_CSV,
        SUMMARY_CSV,
        GATES_CSV,
        RUN_CONTRACT_JSON,
        REPORT_MD,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not accepted:
        raise RuntimeError(
            "R2 frozen-solute NVT preparation "
            "requires review."
        )


if __name__ == "__main__":
    main()
