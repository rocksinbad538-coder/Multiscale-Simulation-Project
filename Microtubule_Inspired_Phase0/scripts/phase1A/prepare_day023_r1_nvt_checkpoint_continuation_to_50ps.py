#!/usr/bin/env python3

from __future__ import annotations

import csv
import difflib
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

DAY023_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design"
)

PREPARATION_ROOT = (
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

OUTPUT_ROOT = (
    DAY023_ROOT
    / "10_r1_nvt_50ps_checkpoint_continuation_preparation"
)

SOURCE_TPR = (
    NVT20_ROOT
    / "r1_frozen_solute_nvt_20ps.tpr"
)

SOURCE_CPT = (
    NVT20_ROOT
    / "r1_frozen_solute_nvt_20ps.cpt"
)

SOURCE_FINAL_GRO = (
    NVT20_ROOT
    / "r1_frozen_solute_nvt_20ps.gro"
)

SOURCE_XTC = (
    NVT20_ROOT
    / "r1_frozen_solute_nvt_20ps.xtc"
)

SOURCE_EDR = (
    NVT20_ROOT
    / "r1_frozen_solute_nvt_20ps.edr"
)

SOURCE_NVT_SUMMARY = (
    NVT20_ROOT
    / "r1_frozen_solute_nvt_20ps_summary.csv"
)

THERMAL_REVIEW_SUMMARY = (
    THERMAL_REVIEW_ROOT
    / "r1_nvt_20ps_thermal_review_summary.csv"
)

WARNING_AUTHORIZATION_SUMMARY = (
    PREPARATION_ROOT
    / "r1_nvt_grompp_warning_authorization_summary.csv"
)

LOCAL_SOURCE_TPR = (
    OUTPUT_ROOT
    / "r1_20ps_source.tpr"
)

LOCAL_SOURCE_CPT = (
    OUTPUT_ROOT
    / "r1_20ps_source_checkpoint.cpt"
)

LOCAL_SOURCE_FINAL_GRO = (
    OUTPUT_ROOT
    / "r1_20ps_source_final.gro"
)

EXTENDED_TPR = (
    OUTPUT_ROOT
    / "r1_frozen_solute_nvt_to_50ps.tpr"
)

SOURCE_TPR_DUMP = (
    OUTPUT_ROOT
    / "r1_20ps_source_tpr_dump.txt"
)

EXTENDED_TPR_DUMP = (
    OUTPUT_ROOT
    / "r1_50ps_extended_tpr_dump.txt"
)

CHECKPOINT_DUMP = (
    OUTPUT_ROOT
    / "r1_20ps_checkpoint_dump.txt"
)

TPR_DIFF = (
    OUTPUT_ROOT
    / "r1_source_vs_extended_tpr_dump.diff"
)

CHECKSUM_CSV = (
    OUTPUT_ROOT
    / "r1_continuation_source_checksums.csv"
)

GATE_CSV = (
    OUTPUT_ROOT
    / "r1_50ps_continuation_preparation_gates.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r1_50ps_continuation_preparation_summary.csv"
)

RUN_CONTRACT_JSON = (
    OUTPUT_ROOT
    / "r1_50ps_continuation_run_contract.json"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R1_NVT_50PS_CHECKPOINT_CONTINUATION_PREPARATION_DAY023.md"
)

EXPECTED_ATOMS = 68314

DT_PS = 0.0005

SOURCE_NSTEPS = 40000
CHECKPOINT_STEP = 40000
CHECKPOINT_TIME_PS = 20.0

EXTENSION_PS = 30.0
EXTENSION_STEPS = 60000

TOTAL_TIME_PS = 50.0
TOTAL_NSTEPS = 100000

NSTXOUT_COMPRESSED = 1000
EXPECTED_CONTINUATION_FRAMES = 61


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
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=cwd,
        env={
            **os.environ,
            "GMX_MAXBACKUP": "-1",
        },
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def run_dump(
    gmx: str,
    arguments: list[str],
    destination: Path,
) -> int:
    """
    Run gmx dump while keeping serialized data and diagnostic
    output in separate streams.

    GROMACS banners, file-reading messages, and the randomized
    closing quotation are written to stderr. Merging stderr into
    stdout can insert those messages inside coordinate or velocity
    records and create false TPR differences.
    """
    completed = subprocess.run(
        [
            gmx,
            "dump",
            *arguments,
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

    destination.write_text(
        completed.stdout,
        encoding="utf-8",
    )

    stderr_path = destination.with_name(
        destination.name
        + ".stderr.log"
    )

    stderr_path.write_text(
        completed.stderr,
        encoding="utf-8",
    )

    return completed.returncode


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


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            chunk = handle.read(
                1024 * 1024
            )

            if not chunk:
                break

            digest.update(chunk)

    return digest.hexdigest()


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


def read_gro_atom_count(path: Path) -> int:
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


def parse_integer(
    text: str,
    patterns: tuple[str, ...],
    label: str,
) -> int:
    for pattern in patterns:
        match = re.search(
            pattern,
            text,
            flags=re.MULTILINE
            | re.IGNORECASE,
        )

        if match is not None:
            return int(
                match.group(1)
            )

    raise RuntimeError(
        f"Could not parse integer field {label}."
    )


def parse_float(
    text: str,
    patterns: tuple[str, ...],
    label: str,
) -> float:
    for pattern in patterns:
        match = re.search(
            pattern,
            text,
            flags=re.MULTILINE
            | re.IGNORECASE,
        )

        if match is not None:
            value = float(
                match.group(1)
            )

            if not math.isfinite(value):
                raise RuntimeError(
                    f"Non-finite field {label}."
                )

            return value

    raise RuntimeError(
        f"Could not parse floating field {label}."
    )


def parse_tpr_dump(
    path: Path,
) -> dict[str, Any]:
    text = path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    return {
        "natoms": parse_integer(
            text,
            (
                r"^\s*natoms\s*=\s*(\d+)",
            ),
            "natoms",
        ),
        "nsteps": parse_integer(
            text,
            (
                r"^\s*nsteps\s*=\s*(-?\d+)",
            ),
            "nsteps",
        ),
        "init_step": parse_integer(
            text,
            (
                r"^\s*init[-_]step\s*=\s*(-?\d+)",
            ),
            "init-step",
        ),
        "dt_ps": parse_float(
            text,
            (
                r"^\s*dt\s*=\s*([-+0-9.eE]+)",
            ),
            "dt",
        ),
        "nstxout_compressed": parse_integer(
            text,
            (
                r"^\s*nstxout[-_]compressed\s*=\s*(\d+)",
            ),
            "nstxout-compressed",
        ),
    }


def parse_checkpoint_dump(
    path: Path,
) -> dict[str, Any]:
    text = path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    step = parse_integer(
        text,
        (
            r"^\s*step\s*=\s*(\d+)",
            r"^\s*step\s+(\d+)",
        ),
        "checkpoint step",
    )

    time_ps = parse_float(
        text,
        (
            r"^\s*t\s*=\s*([-+0-9.eE]+)",
            r"^\s*time\s*=\s*([-+0-9.eE]+)",
        ),
        "checkpoint time",
    )

    velocity_entries = len(
        re.findall(
            r"^\s*v\[\s*\d+\s*\]\s*=",
            text,
            flags=re.MULTILINE,
        )
    )

    coordinate_entries = len(
        re.findall(
            r"^\s*x\[\s*\d+\s*\]\s*=",
            text,
            flags=re.MULTILINE,
        )
    )

    return {
        "step": step,
        "time_ps": time_ps,
        "velocity_entries": (
            velocity_entries
        ),
        "coordinate_entries": (
            coordinate_entries
        ),
    }


def normalize_tpr_dump(
    path: Path,
) -> list[str]:
    """
    Retain only serialized physical TPR content beginning at
    the inputrec block.

    The sole permitted difference introduced by
    gmx convert-tpr -extend is the nsteps value. That field is
    normalized without a regular expression to avoid escaping
    ambiguity.
    """
    raw_lines = path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    inputrec_index = None

    for index, raw_line in enumerate(
        raw_lines
    ):
        if raw_line.strip() == "inputrec:":
            inputrec_index = index
            break

    if inputrec_index is None:
        raise RuntimeError(
            f"Could not locate inputrec block in {path}"
        )

    result = []
    normalized_nsteps_count = 0

    for raw_line in raw_lines[
        inputrec_index:
    ]:
        line = raw_line.rstrip()

        left, separator, right = line.partition(
            "="
        )

        if (
            separator
            and left.strip().lower() == "nsteps"
        ):
            line = (
                f"{left}= <NSTEPS>"
            )

            normalized_nsteps_count += 1

        result.append(line)

    if normalized_nsteps_count != 1:
        raise RuntimeError(
            "Expected exactly one nsteps field in "
            f"{path}; found "
            f"{normalized_nsteps_count}"
        )

    if not result:
        raise RuntimeError(
            f"No physical TPR content was retained from {path}"
        )

    return result


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    required_files = (
        SOURCE_TPR,
        SOURCE_CPT,
        SOURCE_FINAL_GRO,
        SOURCE_XTC,
        SOURCE_EDR,
        SOURCE_NVT_SUMMARY,
        THERMAL_REVIEW_SUMMARY,
        WARNING_AUTHORIZATION_SUMMARY,
    )

    for required in required_files:
        require_file(required)

    nvt_summary = read_single_csv_row(
        SOURCE_NVT_SUMMARY
    )

    thermal_review = read_single_csv_row(
        THERMAL_REVIEW_SUMMARY
    )

    warning_summary = read_single_csv_row(
        WARNING_AUTHORIZATION_SUMMARY
    )

    if (
        thermal_review.get(
            "decision"
        )
        !=
        "R1_INITIAL_THERMALIZATION_TRANSIENT_CONFIRMED"
    ):
        raise RuntimeError(
            "The thermal-transient review did not "
            "authorize checkpoint continuation."
        )

    if not parse_bool(
        thermal_review.get(
            "checkpoint_continuation_preparation_authorized",
            "false",
        )
    ):
        raise RuntimeError(
            "Checkpoint-continuation preparation "
            "is not authorized."
        )

    if (
        warning_summary.get(
            "decision"
        )
        !=
        "R1_NVT_GROMPP_WARNING_AUTHORIZED"
    ):
        raise RuntimeError(
            "The original TPR warning authorization "
            "is not valid."
        )

    source_failed_gates = [
        item.strip()
        for item in nvt_summary.get(
            "failed_gates",
            "",
        ).split("|")
        if item.strip()
    ]

    if source_failed_gates != [
        "temperature_standard_deviation_is_acceptable"
    ]:
        raise RuntimeError(
            "The source NVT has unexpected failed gates: "
            + " | ".join(
                source_failed_gates
            )
        )

    shutil.copy2(
        SOURCE_TPR,
        LOCAL_SOURCE_TPR,
    )

    shutil.copy2(
        SOURCE_CPT,
        LOCAL_SOURCE_CPT,
    )

    shutil.copy2(
        SOURCE_FINAL_GRO,
        LOCAL_SOURCE_FINAL_GRO,
    )

    if EXTENDED_TPR.exists():
        EXTENDED_TPR.unlink()

    gmx = locate_gmx()

    convert = run_command(
        [
            gmx,
            "convert-tpr",
            "-s",
            str(LOCAL_SOURCE_TPR),
            "-extend",
            f"{EXTENSION_PS:.6f}",
            "-o",
            str(EXTENDED_TPR),
        ],
        cwd=OUTPUT_ROOT,
    )

    convert_log = (
        OUTPUT_ROOT
        / "r1_50ps_convert_tpr.log"
    )

    convert_log.write_text(
        convert.stdout,
        encoding="utf-8",
    )

    if (
        convert.returncode != 0
        or not EXTENDED_TPR.exists()
        or EXTENDED_TPR.stat().st_size == 0
    ):
        raise RuntimeError(
            "gmx convert-tpr failed. "
            f"See {convert_log}"
        )

    source_dump_rc = run_dump(
        gmx,
        [
            "-s",
            str(LOCAL_SOURCE_TPR),
        ],
        SOURCE_TPR_DUMP,
    )

    extended_dump_rc = run_dump(
        gmx,
        [
            "-s",
            str(EXTENDED_TPR),
        ],
        EXTENDED_TPR_DUMP,
    )

    checkpoint_dump_rc = run_dump(
        gmx,
        [
            "-cp",
            str(LOCAL_SOURCE_CPT),
        ],
        CHECKPOINT_DUMP,
    )

    if source_dump_rc != 0:
        raise RuntimeError(
            "Could not dump the source TPR."
        )

    if extended_dump_rc != 0:
        raise RuntimeError(
            "Could not dump the extended TPR."
        )

    if checkpoint_dump_rc != 0:
        raise RuntimeError(
            "Could not dump the source checkpoint."
        )

    source_tpr = parse_tpr_dump(
        SOURCE_TPR_DUMP
    )

    extended_tpr = parse_tpr_dump(
        EXTENDED_TPR_DUMP
    )

    checkpoint = parse_checkpoint_dump(
        CHECKPOINT_DUMP
    )

    source_normalized = normalize_tpr_dump(
        SOURCE_TPR_DUMP
    )

    extended_normalized = normalize_tpr_dump(
        EXTENDED_TPR_DUMP
    )

    differences = list(
        difflib.unified_diff(
            source_normalized,
            extended_normalized,
            fromfile="source_20ps_tpr",
            tofile="extended_50ps_tpr",
            lineterm="",
        )
    )

    TPR_DIFF.write_text(
        "\n".join(differences)
        + (
            "\n"
            if differences
            else ""
        ),
        encoding="utf-8",
    )

    gro_atoms = read_gro_atom_count(
        LOCAL_SOURCE_FINAL_GRO
    )

    source_total_time = (
        source_tpr[
            "nsteps"
        ]
        * source_tpr[
            "dt_ps"
        ]
    )

    extended_total_time = (
        extended_tpr[
            "nsteps"
        ]
        * extended_tpr[
            "dt_ps"
        ]
    )

    remaining_steps = (
        extended_tpr[
            "nsteps"
        ]
        - checkpoint[
            "step"
        ]
    )

    remaining_time_ps = (
        extended_total_time
        - checkpoint[
            "time_ps"
        ]
    )

    expected_continuation_frames = (
        remaining_steps
        // extended_tpr[
            "nstxout_compressed"
        ]
        + 1
    )

    checksum_rows = []

    for label, path in (
        (
            "source_TPR",
            SOURCE_TPR,
        ),
        (
            "source_checkpoint",
            SOURCE_CPT,
        ),
        (
            "source_final_GRO",
            SOURCE_FINAL_GRO,
        ),
        (
            "source_XTC",
            SOURCE_XTC,
        ),
        (
            "source_EDR",
            SOURCE_EDR,
        ),
        (
            "local_checkpoint_copy",
            LOCAL_SOURCE_CPT,
        ),
        (
            "extended_TPR",
            EXTENDED_TPR,
        ),
    ):
        checksum_rows.append(
            {
                "label": label,
                "path": relative(path),
                "size_bytes": (
                    path.stat().st_size
                ),
                "sha256": sha256(path),
            }
        )

    write_csv(
        CHECKSUM_CSV,
        checksum_rows,
    )

    source_checkpoint_hash = sha256(
        SOURCE_CPT
    )

    local_checkpoint_hash = sha256(
        LOCAL_SOURCE_CPT
    )

    gates = {
        "convert_tpr_return_code_zero": (
            convert.returncode == 0
        ),
        "source_TPR_dump_return_code_zero": (
            source_dump_rc == 0
        ),
        "extended_TPR_dump_return_code_zero": (
            extended_dump_rc == 0
        ),
        "checkpoint_dump_return_code_zero": (
            checkpoint_dump_rc == 0
        ),
        "source_TPR_has_68314_atoms": (
            source_tpr[
                "natoms"
            ]
            == EXPECTED_ATOMS
        ),
        "extended_TPR_has_68314_atoms": (
            extended_tpr[
                "natoms"
            ]
            == EXPECTED_ATOMS
        ),
        "source_final_GRO_has_68314_atoms": (
            gro_atoms
            == EXPECTED_ATOMS
        ),
        "source_TPR_has_40000_steps": (
            source_tpr[
                "nsteps"
            ]
            == SOURCE_NSTEPS
        ),
        "extended_TPR_has_100000_steps": (
            extended_tpr[
                "nsteps"
            ]
            == TOTAL_NSTEPS
        ),
        "source_dt_is_0p0005ps": (
            math.isclose(
                source_tpr[
                    "dt_ps"
                ],
                DT_PS,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
        ),
        "extended_dt_is_0p0005ps": (
            math.isclose(
                extended_tpr[
                    "dt_ps"
                ],
                DT_PS,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
        ),
        "source_total_time_is_20ps": (
            math.isclose(
                source_total_time,
                CHECKPOINT_TIME_PS,
                rel_tol=0.0,
                abs_tol=1.0e-9,
            )
        ),
        "extended_total_time_is_50ps": (
            math.isclose(
                extended_total_time,
                TOTAL_TIME_PS,
                rel_tol=0.0,
                abs_tol=1.0e-9,
            )
        ),
        "checkpoint_step_is_40000": (
            checkpoint[
                "step"
            ]
            == CHECKPOINT_STEP
        ),
        "checkpoint_time_is_20ps": (
            math.isclose(
                checkpoint[
                    "time_ps"
                ],
                CHECKPOINT_TIME_PS,
                rel_tol=0.0,
                abs_tol=1.0e-8,
            )
        ),
        "remaining_steps_are_60000": (
            remaining_steps
            == EXTENSION_STEPS
        ),
        "remaining_time_is_30ps": (
            math.isclose(
                remaining_time_ps,
                EXTENSION_PS,
                rel_tol=0.0,
                abs_tol=1.0e-8,
            )
        ),
        "trajectory_stride_is_preserved": (
            source_tpr[
                "nstxout_compressed"
            ]
            == NSTXOUT_COMPRESSED
            and extended_tpr[
                "nstxout_compressed"
            ]
            == NSTXOUT_COMPRESSED
        ),
        "expected_continuation_frames_are_61": (
            expected_continuation_frames
            == EXPECTED_CONTINUATION_FRAMES
        ),
        "no_TPR_changes_beyond_nsteps": (
            len(differences) == 0
        ),
        "checkpoint_copy_is_bitwise_identical": (
            source_checkpoint_hash
            == local_checkpoint_hash
        ),
        "thermal_review_authorized_continuation": (
            parse_bool(
                thermal_review[
                    "checkpoint_continuation_preparation_authorized"
                ]
            )
        ),
        "no_unreviewed_source_gate_failures": (
            source_failed_gates
            == [
                "temperature_standard_deviation_is_acceptable"
            ]
        ),
    }

    failed_gates = [
        name
        for name, passed
        in gates.items()
        if not passed
    ]

    prepared = (
        len(failed_gates) == 0
    )

    decision = (
        "R1_CHECKPOINT_CONTINUATION_TO_50PS_PREPARED"
        if prepared
        else
        "R1_CHECKPOINT_CONTINUATION_PREPARATION_REQUIRES_REVIEW"
    )

    required_next_step = (
        "RUN_R1_20_TO_50PS_CHECKPOINT_CONTINUATION"
        if prepared
        else
        "RESOLVE_R1_CONTINUATION_PREPARATION_GATE_FAILURES"
    )

    continuation_output_root = (
        DAY023_ROOT
        / "11_r1_frozen_solute_nvt_20_to_50ps"
    )

    future_deffnm = (
        continuation_output_root
        / "r1_frozen_solute_nvt_20_to_50ps"
    )

    run_contract = {
        "decision": decision,
        "source_trajectory_end_ps": (
            CHECKPOINT_TIME_PS
        ),
        "target_total_time_ps": (
            TOTAL_TIME_PS
        ),
        "continuation_duration_ps": (
            EXTENSION_PS
        ),
        "source_checkpoint_step": (
            checkpoint[
                "step"
            ]
        ),
        "source_checkpoint_time_ps": (
            checkpoint[
                "time_ps"
            ]
        ),
        "source_checkpoint": relative(
            LOCAL_SOURCE_CPT
        ),
        "extended_TPR": relative(
            EXTENDED_TPR
        ),
        "expected_remaining_steps": (
            remaining_steps
        ),
        "expected_continuation_frames": (
            expected_continuation_frames
        ),
        "velocity_regeneration": False,
        "new_independent_trajectory": False,
        "checkpoint_state_must_be_used": True,
        "append_to_original_files": False,
        "future_output_root": relative(
            continuation_output_root
        ),
        "future_deffnm": relative(
            future_deffnm
        ),
        "future_mdrun_arguments": [
            gmx,
            "mdrun",
            "-s",
            str(
                EXTENDED_TPR.resolve()
            ),
            "-cpi",
            str(
                LOCAL_SOURCE_CPT.resolve()
            ),
            "-deffnm",
            str(
                future_deffnm.resolve()
            ),
            "-noappend",
            "-ntmpi",
            "1",
            "-ntomp",
            "4",
        ],
        "execution_authorized": (
            prepared
        ),
        "long_mobile_MD_authorized": False,
        "multitemperature_MD_authorized": False,
        "QM_recalculation_authorized": False,
        "required_next_step": (
            required_next_step
        ),
    }

    RUN_CONTRACT_JSON.write_text(
        json.dumps(
            run_contract,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
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

    summary = {
        "decision": decision,
        "convert_tpr_return_code": (
            convert.returncode
        ),
        "source_TPR_atoms": (
            source_tpr[
                "natoms"
            ]
        ),
        "extended_TPR_atoms": (
            extended_tpr[
                "natoms"
            ]
        ),
        "source_TPR_nsteps": (
            source_tpr[
                "nsteps"
            ]
        ),
        "extended_TPR_nsteps": (
            extended_tpr[
                "nsteps"
            ]
        ),
        "dt_ps": (
            extended_tpr[
                "dt_ps"
            ]
        ),
        "source_total_time_ps": (
            source_total_time
        ),
        "target_total_time_ps": (
            extended_total_time
        ),
        "checkpoint_step": (
            checkpoint[
                "step"
            ]
        ),
        "checkpoint_time_ps": (
            checkpoint[
                "time_ps"
            ]
        ),
        "remaining_steps": (
            remaining_steps
        ),
        "remaining_time_ps": (
            remaining_time_ps
        ),
        "nstxout_compressed": (
            extended_tpr[
                "nstxout_compressed"
            ]
        ),
        "expected_continuation_frames": (
            expected_continuation_frames
        ),
        "checkpoint_coordinate_entries_in_dump": (
            checkpoint[
                "coordinate_entries"
            ]
        ),
        "checkpoint_velocity_entries_in_dump": (
            checkpoint[
                "velocity_entries"
            ]
        ),
        "TPR_dump_differences_beyond_nsteps": (
            len(differences)
        ),
        "checkpoint_copy_sha256": (
            local_checkpoint_hash
        ),
        "failed_gates": (
            " | ".join(
                failed_gates
            )
        ),
        "continuation_execution_authorized": (
            prepared
        ),
        "velocity_regeneration_authorized": False,
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

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed in gates.items()
    )

    REPORT_MD.write_text(
        f"""# R1 NVT Checkpoint Continuation to 50 ps

## Purpose

The 20 ps R1 frozen-solute trajectory showed a short initial
thermalization transient followed by stable behavior:

- 5–20 ps temperature:
  **{float(thermal_review['post_5ps_temperature_mean_K']):.4f}
  ± {float(thermal_review['post_5ps_temperature_std_K']):.4f} K**
- 10–20 ps temperature:
  **{float(thermal_review['last_10ps_temperature_mean_K']):.4f}
  ± {float(thermal_review['last_10ps_temperature_std_K']):.4f} K**
- Endpoint lumen occupancy:
  **428 waters**
- Endpoint initially luminal waters retained:
  **428 waters**

A 30 ps checkpoint continuation was therefore prepared to extend the
same trajectory from 20 to 50 ps.

## Continuation method

The original TPR was extended using `gmx convert-tpr`.

No new `grompp` operation was performed. No velocities were generated.

The future run must use:

- the extended TPR;
- the exact 20 ps checkpoint;
- `mdrun -cpi`;
- `-noappend`, preserving the original 20 ps result.

This retains the checkpoint coordinates, velocities, thermostat state,
random state, and current integration step.

## Source state

- Source atoms:
  **{source_tpr['natoms']}**
- Source steps:
  **{source_tpr['nsteps']}**
- Source nominal duration:
  **{source_total_time:.6f} ps**
- Checkpoint step:
  **{checkpoint['step']}**
- Checkpoint time:
  **{checkpoint['time_ps']:.6f} ps**
- Final GRO atoms:
  **{gro_atoms}**
- Checkpoint SHA256:
  `{local_checkpoint_hash}`

## Extended TPR

- Extended atoms:
  **{extended_tpr['natoms']}**
- Extended total steps:
  **{extended_tpr['nsteps']}**
- Extended target time:
  **{extended_total_time:.6f} ps**
- Remaining steps:
  **{remaining_steps}**
- Remaining duration:
  **{remaining_time_ps:.6f} ps**
- Compressed trajectory stride:
  **{extended_tpr['nstxout_compressed']} steps**
- Expected continuation frames:
  **{expected_continuation_frames}**
- TPR differences beyond `nsteps`:
  **{len(differences)}**

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Checkpoint-continuation execution authorized:
  **{'YES' if prepared else 'NO'}**
- Velocity regeneration authorized:
  **NO**
- New independent trajectory authorized:
  **NO**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `{required_next_step}`

The continuation remains part of the R1 frozen steric positive-control
screening. It does not establish chemical realizability.
""",
        encoding="utf-8",
    )

    print(
        "Day023 R1 checkpoint continuation to "
        "50 ps preparation completed."
    )

    print(
        "convert-tpr / source dump / extended dump / "
        "checkpoint dump return codes: "
        f"{convert.returncode}/"
        f"{source_dump_rc}/"
        f"{extended_dump_rc}/"
        f"{checkpoint_dump_rc}"
    )

    print(
        "Source / extended atoms: "
        f"{source_tpr['natoms']}/"
        f"{extended_tpr['natoms']}"
    )

    print(
        "Source / extended nsteps: "
        f"{source_tpr['nsteps']}/"
        f"{extended_tpr['nsteps']}"
    )

    print(
        "Source / target total time: "
        f"{source_total_time:.6f}/"
        f"{extended_total_time:.6f} ps"
    )

    print(
        "Checkpoint step / time: "
        f"{checkpoint['step']} / "
        f"{checkpoint['time_ps']:.6f} ps"
    )

    print(
        "Remaining steps / duration: "
        f"{remaining_steps} / "
        f"{remaining_time_ps:.6f} ps"
    )

    print(
        "Trajectory stride / expected continuation frames: "
        f"{extended_tpr['nstxout_compressed']} / "
        f"{expected_continuation_frames}"
    )

    print(
        "Checkpoint coordinates / velocities represented "
        "in dump: "
        f"{checkpoint['coordinate_entries']}/"
        f"{checkpoint['velocity_entries']}"
    )

    print(
        "TPR differences beyond nsteps: "
        f"{len(differences)}"
    )

    print(
        "Checkpoint copy bitwise identical: "
        f"{'YES' if source_checkpoint_hash == local_checkpoint_hash else 'NO'}"
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
        "Checkpoint continuation execution authorized: "
        f"{'YES' if prepared else 'NO'}"
    )

    print(
        "Velocity regeneration authorized: NO"
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
        f"Wrote: {relative(EXTENDED_TPR)}"
    )

    print(
        f"Wrote: {relative(LOCAL_SOURCE_CPT)}"
    )

    print(
        f"Wrote: {relative(SOURCE_TPR_DUMP)}"
    )

    print(
        f"Wrote: {relative(EXTENDED_TPR_DUMP)}"
    )

    print(
        f"Wrote: {relative(CHECKPOINT_DUMP)}"
    )

    print(
        f"Wrote: {relative(TPR_DIFF)}"
    )

    print(
        f"Wrote: {relative(CHECKSUM_CSV)}"
    )

    print(
        f"Wrote: {relative(GATE_CSV)}"
    )

    print(
        f"Wrote: {relative(SUMMARY_CSV)}"
    )

    print(
        f"Wrote: {relative(RUN_CONTRACT_JSON)}"
    )

    print(
        f"Wrote: {relative(REPORT_MD)}"
    )

    if not prepared:
        raise RuntimeError(
            "R1 checkpoint continuation preparation "
            "requires review."
        )


if __name__ == "__main__":
    main()
