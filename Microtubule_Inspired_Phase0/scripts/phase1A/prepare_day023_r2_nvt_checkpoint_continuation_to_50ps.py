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

PREPARATION_20PS_ROOT = (
    DAY023_ROOT
    / "15_r2_frozen_solute_nvt_20ps_preparation"
)

RUN_20PS_ROOT = (
    DAY023_ROOT
    / "16_r2_frozen_solute_nvt_20ps"
)

TRANSIENT_AUDIT_ROOT = (
    DAY023_ROOT
    / "17_r2_nvt_20ps_occupancy_transient_audit"
)

OUTPUT_ROOT = (
    DAY023_ROOT
    / "18_r2_nvt_50ps_checkpoint_continuation_preparation"
)

SOURCE_TPR = (
    PREPARATION_20PS_ROOT
    / "r2_frozen_solute_nvt_20ps.tpr"
)

SOURCE_PREPARATION_SUMMARY = (
    PREPARATION_20PS_ROOT
    / "r2_frozen_solute_nvt_20ps_preparation_summary.csv"
)

SOURCE_CHECKPOINT = (
    RUN_20PS_ROOT
    / "r2_frozen_solute_nvt_20ps.cpt"
)

SOURCE_FINAL_GRO = (
    RUN_20PS_ROOT
    / "r2_frozen_solute_nvt_20ps.gro"
)

SOURCE_RUN_SUMMARY = (
    RUN_20PS_ROOT
    / "r2_frozen_solute_nvt_20ps_summary.csv"
)

SOURCE_MDRUN_CONSOLE = (
    RUN_20PS_ROOT
    / "r2_frozen_solute_nvt_20ps_mdrun_console.log"
)

SOURCE_MDRUN_LOG = (
    RUN_20PS_ROOT
    / "r2_frozen_solute_nvt_20ps.log"
)

TRANSIENT_AUDIT_SUMMARY = (
    TRANSIENT_AUDIT_ROOT
    / "r2_nvt_20ps_occupancy_transient_audit_summary.csv"
)

EXTENDED_TPR = (
    OUTPUT_ROOT
    / "r2_frozen_solute_nvt_to_50ps.tpr"
)

COPIED_CHECKPOINT = (
    OUTPUT_ROOT
    / "r2_20ps_source_checkpoint.cpt"
)

SOURCE_TPR_DUMP = (
    OUTPUT_ROOT
    / "r2_20ps_source_tpr_dump.txt"
)

SOURCE_TPR_DUMP_STDERR = (
    OUTPUT_ROOT
    / "r2_20ps_source_tpr_dump.stderr.log"
)

EXTENDED_TPR_DUMP = (
    OUTPUT_ROOT
    / "r2_50ps_extended_tpr_dump.txt"
)

EXTENDED_TPR_DUMP_STDERR = (
    OUTPUT_ROOT
    / "r2_50ps_extended_tpr_dump.stderr.log"
)

CHECKPOINT_DUMP = (
    OUTPUT_ROOT
    / "r2_20ps_checkpoint_dump.txt"
)

CHECKPOINT_DUMP_STDERR = (
    OUTPUT_ROOT
    / "r2_20ps_checkpoint_dump.stderr.log"
)

CONVERT_TPR_LOG = (
    OUTPUT_ROOT
    / "r2_convert_tpr_extend_30ps.log"
)

PHYSICAL_DIFF = (
    OUTPUT_ROOT
    / "r2_source_vs_extended_tpr_physical_diff.txt"
)

CHECKSUMS_CSV = (
    OUTPUT_ROOT
    / "r2_50ps_continuation_checksums.csv"
)

SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_50ps_continuation_preparation_summary.csv"
)

GATES_CSV = (
    OUTPUT_ROOT
    / "r2_50ps_continuation_preparation_gates.csv"
)

RUN_CONTRACT_JSON = (
    OUTPUT_ROOT
    / "r2_50ps_continuation_run_contract.json"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_CHECKPOINT_CONTINUATION_TO_50PS_PREPARATION_DAY023.md"
)

EXPECTED_PREPARATION_DECISION = (
    "R2_FROZEN_SOLUTE_NVT_20PS_PREPARED"
)

EXPECTED_RUN_DECISION = (
    "R2_FROZEN_SOLUTE_NVT_20PS_REQUIRES_REVIEW"
)

EXPECTED_AUDIT_DECISION = (
    "R2_OCCUPANCY_TRANSIENT_CHECKPOINT_EXTENSION_JUSTIFIED"
)

EXPECTED_FAILED_SOURCE_GATE = (
    "second_half_occupancy_slope_is_acceptable"
)

EXPECTED_ATOMS = 68332
EXPECTED_COORDINATE_COMPONENTS = EXPECTED_ATOMS * 3

DT_PS = 0.0005

SOURCE_NSTEPS = 40000
TARGET_NSTEPS = 100000
REMAINING_STEPS = TARGET_NSTEPS - SOURCE_NSTEPS

SOURCE_TIME_PS = SOURCE_NSTEPS * DT_PS
TARGET_TIME_PS = TARGET_NSTEPS * DT_PS
EXTENSION_TIME_PS = TARGET_TIME_PS - SOURCE_TIME_PS

XTC_STRIDE = 1000
XTC_INTERVAL_PS = XTC_STRIDE * DT_PS

EXPECTED_CONTINUATION_FRAMES = (
    REMAINING_STEPS
    // XTC_STRIDE
    + 1
)

EXECUTION_OUTPUT_ROOT = (
    DAY023_ROOT
    / "19_r2_frozen_solute_nvt_20_to_50ps"
)

EXECUTION_DEFFNM = (
    EXECUTION_OUTPUT_ROOT
    / "r2_frozen_solute_nvt_20_to_50ps"
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
            f"Expected one row in {path}; found {len(rows)}"
        )

    return rows[0]


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


def as_float(
    row: dict[str, str],
    key: str,
) -> float:
    try:
        value = float(row[key])
    except (
        KeyError,
        TypeError,
        ValueError,
    ) as exc:
        raise RuntimeError(
            f"Could not parse numeric field {key!r}"
        ) from exc

    if not math.isfinite(value):
        raise RuntimeError(
            f"Non-finite field {key!r}"
        )

    return value


def as_int(
    row: dict[str, str],
    key: str,
) -> int:
    return int(
        round(
            as_float(row, key)
        )
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
        stderr=subprocess.PIPE,
        check=False,
    )


def parse_scalar(
    text: str,
    key: str,
) -> str:
    match = re.search(
        rf"^\s*{re.escape(key)}\s*=\s*(\S+)",
        text,
        flags=re.MULTILINE,
    )

    if match is None:
        raise RuntimeError(
            f"Could not parse {key!r}"
        )

    return match.group(1)


def parse_checkpoint_step(
    text: str,
) -> int:
    patterns = (
        r"^\s*step\s*=\s*(\d+)",
        r"^\s*step\s+(\d+)",
    )

    for pattern in patterns:
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
    patterns = (
        r"^\s*t\s*=\s*([-+0-9.eE]+)",
        r"^\s*time\s*=\s*([-+0-9.eE]+)",
    )

    for pattern in patterns:
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


def count_state_entries(
    text: str,
    symbol: str,
) -> int:
    return len(
        re.findall(
            rf"(?m)^\s*{re.escape(symbol)}"
            r"\[\s*\d+\s*\]",
            text,
        )
    )


def normalize_tpr_dump(
    text: str,
) -> list[str]:
    normalized = []

    for raw_line in text.splitlines():
        line = raw_line.rstrip()

        if re.match(
            r"^\s*nsteps\s*=",
            line,
        ):
            line = re.sub(
                r"(nsteps\s*=\s*)\S+",
                r"\1<NORMALIZED_NSTEPS>",
                line,
                count=1,
            )

        normalized.append(line)

    return normalized


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    required_files = (
        SOURCE_TPR,
        SOURCE_PREPARATION_SUMMARY,
        SOURCE_CHECKPOINT,
        SOURCE_FINAL_GRO,
        SOURCE_RUN_SUMMARY,
        SOURCE_MDRUN_CONSOLE,
        SOURCE_MDRUN_LOG,
        TRANSIENT_AUDIT_SUMMARY,
    )

    for required in required_files:
        require_file(required)

    for product in (
        EXTENDED_TPR,
        COPIED_CHECKPOINT,
        SUMMARY_CSV,
        GATES_CSV,
        RUN_CONTRACT_JSON,
    ):
        if product.exists():
            raise RuntimeError(
                f"Preparation product already exists: {product}"
            )

    preparation = read_single_csv_row(
        SOURCE_PREPARATION_SUMMARY
    )

    run_summary = read_single_csv_row(
        SOURCE_RUN_SUMMARY
    )

    transient_audit = read_single_csv_row(
        TRANSIENT_AUDIT_SUMMARY
    )

    if (
        preparation.get("decision")
        != EXPECTED_PREPARATION_DECISION
    ):
        raise RuntimeError(
            "The original 20 ps preparation is not accepted."
        )

    if (
        run_summary.get("decision")
        != EXPECTED_RUN_DECISION
    ):
        raise RuntimeError(
            "The 20 ps source run is not in the expected "
            "nonstationary review state."
        )

    if (
        transient_audit.get("decision")
        != EXPECTED_AUDIT_DECISION
    ):
        raise RuntimeError(
            "The transient audit did not authorize "
            "checkpoint-extension preparation."
        )

    if not parse_bool(
        transient_audit.get(
            "checkpoint_continuation_preparation_authorized",
            "false",
        )
    ):
        raise RuntimeError(
            "Checkpoint-continuation preparation is not authorized."
        )

    source_failed_gates = (
        run_summary.get(
            "failed_gates",
            "",
        )
        .strip()
    )

    if source_failed_gates != EXPECTED_FAILED_SOURCE_GATE:
        raise RuntimeError(
            "The source run failed an unexpected gate: "
            f"{source_failed_gates!r}"
        )

    if as_int(
        run_summary,
        "mdrun_return_code",
    ) != 0:
        raise RuntimeError(
            "The source mdrun did not complete successfully."
        )

    if as_int(
        run_summary,
        "trajectory_check_return_code",
    ) != 0:
        raise RuntimeError(
            "The source trajectory check did not pass."
        )

    if not parse_bool(
        run_summary.get(
            "completion_confirmed",
            "false",
        )
    ):
        raise RuntimeError(
            "Source-run completion was not confirmed."
        )

    if not parse_bool(
        run_summary.get(
            "checkpoint_written",
            "false",
        )
    ):
        raise RuntimeError(
            "The source checkpoint was not confirmed."
        )

    if as_int(
        run_summary,
        "instability_signature_count",
    ) != 0:
        raise RuntimeError(
            "The source run contains instability signatures."
        )

    if as_int(
        run_summary,
        "trajectory_frames",
    ) != 41:
        raise RuntimeError(
            "The source trajectory does not contain 41 frames."
        )

    if not math.isclose(
        as_float(
            run_summary,
            "trajectory_end_ps",
        ),
        SOURCE_TIME_PS,
        rel_tol=0.0,
        abs_tol=1.0e-6,
    ):
        raise RuntimeError(
            "The source trajectory does not end at 20 ps."
        )

    completion_text = (
        SOURCE_MDRUN_CONSOLE.read_text(
            encoding="utf-8",
            errors="replace",
        )
        + "\n"
        + SOURCE_MDRUN_LOG.read_text(
            encoding="utf-8",
            errors="replace",
        )
    )

    if not (
        "Writing final coordinates"
        in completion_text
        or re.search(
            r"\bFinished\s+mdrun\b",
            completion_text,
            flags=re.IGNORECASE,
        )
        is not None
    ):
        raise RuntimeError(
            "The source logs do not confirm completion."
        )

    instability_patterns = (
        re.compile(
            r"\bnan\b",
            flags=re.IGNORECASE,
        ),
        re.compile(
            r"\bfatal\s+error\b",
            flags=re.IGNORECASE,
        ),
        re.compile(
            r"\blincs\s+warning\b",
            flags=re.IGNORECASE,
        ),
        re.compile(
            r"\bconstraint\s+warning\b",
            flags=re.IGNORECASE,
        ),
        re.compile(
            r"\bsegmentation\s+fault\b",
            flags=re.IGNORECASE,
        ),
    )

    if any(
        pattern.search(completion_text)
        for pattern in instability_patterns
    ):
        raise RuntimeError(
            "A strict instability signature was found "
            "in the source logs."
        )

    gmx = locate_gmx()

    convert = run_command(
        [
            gmx,
            "convert-tpr",
            "-s",
            str(SOURCE_TPR),
            "-o",
            str(EXTENDED_TPR),
            "-extend",
            f"{EXTENSION_TIME_PS:.6f}",
        ],
        cwd=OUTPUT_ROOT,
    )

    CONVERT_TPR_LOG.write_text(
        convert.stdout
        + "\n"
        + convert.stderr,
        encoding="utf-8",
    )

    if (
        convert.returncode != 0
        or not EXTENDED_TPR.exists()
        or EXTENDED_TPR.stat().st_size == 0
    ):
        raise RuntimeError(
            "gmx convert-tpr failed. "
            f"See {CONVERT_TPR_LOG}"
        )

    source_dump = run_command(
        [
            gmx,
            "dump",
            "-s",
            str(SOURCE_TPR),
        ],
        cwd=OUTPUT_ROOT,
    )

    SOURCE_TPR_DUMP.write_text(
        source_dump.stdout,
        encoding="utf-8",
    )

    SOURCE_TPR_DUMP_STDERR.write_text(
        source_dump.stderr,
        encoding="utf-8",
    )

    extended_dump = run_command(
        [
            gmx,
            "dump",
            "-s",
            str(EXTENDED_TPR),
        ],
        cwd=OUTPUT_ROOT,
    )

    EXTENDED_TPR_DUMP.write_text(
        extended_dump.stdout,
        encoding="utf-8",
    )

    EXTENDED_TPR_DUMP_STDERR.write_text(
        extended_dump.stderr,
        encoding="utf-8",
    )

    checkpoint_dump = run_command(
        [
            gmx,
            "dump",
            "-cp",
            str(SOURCE_CHECKPOINT),
        ],
        cwd=OUTPUT_ROOT,
    )

    CHECKPOINT_DUMP.write_text(
        checkpoint_dump.stdout,
        encoding="utf-8",
    )

    CHECKPOINT_DUMP_STDERR.write_text(
        checkpoint_dump.stderr,
        encoding="utf-8",
    )

    if (
        source_dump.returncode != 0
        or extended_dump.returncode != 0
        or checkpoint_dump.returncode != 0
    ):
        raise RuntimeError(
            "One or more GROMACS dump operations failed."
        )

    source_text = source_dump.stdout
    extended_text = extended_dump.stdout
    checkpoint_text = checkpoint_dump.stdout

    source_atoms = int(
        parse_scalar(
            source_text,
            "natoms",
        )
    )

    extended_atoms = int(
        parse_scalar(
            extended_text,
            "natoms",
        )
    )

    source_nsteps = int(
        parse_scalar(
            source_text,
            "nsteps",
        )
    )

    extended_nsteps = int(
        parse_scalar(
            extended_text,
            "nsteps",
        )
    )

    source_dt = float(
        parse_scalar(
            source_text,
            "dt",
        )
    )

    extended_dt = float(
        parse_scalar(
            extended_text,
            "dt",
        )
    )

    source_integrator = parse_scalar(
        source_text,
        "integrator",
    ).lower()

    extended_integrator = parse_scalar(
        extended_text,
        "integrator",
    ).lower()

    source_continuation = parse_scalar(
        source_text,
        "continuation",
    ).lower()

    extended_continuation = parse_scalar(
        extended_text,
        "continuation",
    ).lower()

    source_stride = int(
        parse_scalar(
            source_text,
            "nstxout-compressed",
        )
    )

    extended_stride = int(
        parse_scalar(
            extended_text,
            "nstxout-compressed",
        )
    )

    normalized_source = normalize_tpr_dump(
        source_text
    )

    normalized_extended = normalize_tpr_dump(
        extended_text
    )

    diff_lines = list(
        difflib.unified_diff(
            normalized_source,
            normalized_extended,
            fromfile="source_20ps_tpr",
            tofile="extended_50ps_tpr",
            lineterm="",
        )
    )

    PHYSICAL_DIFF.write_text(
        (
            "\n".join(diff_lines)
            + ("\n" if diff_lines else "")
        ),
        encoding="utf-8",
    )

    checkpoint_step = parse_checkpoint_step(
        checkpoint_text
    )

    checkpoint_time = parse_checkpoint_time(
        checkpoint_text
    )

    checkpoint_coordinate_entries = (
        count_state_entries(
            checkpoint_text,
            "x",
        )
    )

    checkpoint_velocity_entries = (
        count_state_entries(
            checkpoint_text,
            "v",
        )
    )

    allowed_state_entry_counts = {
        EXPECTED_ATOMS,
        EXPECTED_COORDINATE_COMPONENTS,
    }

    shutil.copy2(
        SOURCE_CHECKPOINT,
        COPIED_CHECKPOINT,
    )

    source_checkpoint_hash = sha256(
        SOURCE_CHECKPOINT
    )

    copied_checkpoint_hash = sha256(
        COPIED_CHECKPOINT
    )

    checkpoint_bitwise_copy = (
        source_checkpoint_hash
        == copied_checkpoint_hash
        and SOURCE_CHECKPOINT.stat().st_size
        == COPIED_CHECKPOINT.stat().st_size
    )

    gates = {
        "source_preparation_is_valid": (
            preparation.get("decision")
            == EXPECTED_PREPARATION_DECISION
        ),
        "source_run_is_nonstationary_review_state": (
            run_summary.get("decision")
            == EXPECTED_RUN_DECISION
        ),
        "transient_audit_authorized_extension_preparation": (
            transient_audit.get("decision")
            == EXPECTED_AUDIT_DECISION
            and parse_bool(
                transient_audit.get(
                    "checkpoint_continuation_preparation_authorized",
                    "false",
                )
            )
        ),
        "source_failed_only_stationarity_gate": (
            source_failed_gates
            == EXPECTED_FAILED_SOURCE_GATE
        ),
        "source_run_completed_without_instability": (
            as_int(
                run_summary,
                "mdrun_return_code",
            )
            == 0
            and as_int(
                run_summary,
                "instability_signature_count",
            )
            == 0
        ),
        "convert_tpr_return_code_zero": (
            convert.returncode == 0
        ),
        "source_TPR_dump_return_code_zero": (
            source_dump.returncode == 0
        ),
        "extended_TPR_dump_return_code_zero": (
            extended_dump.returncode == 0
        ),
        "checkpoint_dump_return_code_zero": (
            checkpoint_dump.returncode == 0
        ),
        "source_TPR_has_68332_atoms": (
            source_atoms
            == EXPECTED_ATOMS
        ),
        "extended_TPR_has_68332_atoms": (
            extended_atoms
            == EXPECTED_ATOMS
        ),
        "source_TPR_nsteps_is_40000": (
            source_nsteps
            == SOURCE_NSTEPS
        ),
        "extended_TPR_nsteps_is_100000": (
            extended_nsteps
            == TARGET_NSTEPS
        ),
        "source_and_extended_dt_are_0p0005ps": (
            math.isclose(
                source_dt,
                DT_PS,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            and math.isclose(
                extended_dt,
                DT_PS,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
        ),
        "source_and_extended_integrators_are_md": (
            source_integrator == "md"
            and extended_integrator == "md"
        ),
        "continuation_flag_is_preserved": (
            source_continuation
            == extended_continuation
        ),
        "XTC_stride_is_preserved_at_1000": (
            source_stride
            == XTC_STRIDE
            and extended_stride
            == XTC_STRIDE
        ),
        "TPRs_differ_only_by_nsteps": (
            len(diff_lines) == 0
        ),
        "checkpoint_step_is_40000": (
            checkpoint_step
            == SOURCE_NSTEPS
        ),
        "checkpoint_time_is_20ps": (
            math.isclose(
                checkpoint_time,
                SOURCE_TIME_PS,
                rel_tol=0.0,
                abs_tol=1.0e-6,
            )
        ),
        "checkpoint_coordinate_entries_are_complete": (
            checkpoint_coordinate_entries
            in allowed_state_entry_counts
        ),
        "checkpoint_velocity_entries_are_complete": (
            checkpoint_velocity_entries
            in allowed_state_entry_counts
        ),
        "checkpoint_coordinate_and_velocity_counts_match": (
            checkpoint_coordinate_entries
            == checkpoint_velocity_entries
        ),
        "checkpoint_copy_is_bitwise_identical": (
            checkpoint_bitwise_copy
        ),
        "remaining_steps_are_60000": (
            REMAINING_STEPS
            == 60000
        ),
        "remaining_time_is_30ps": (
            math.isclose(
                EXTENSION_TIME_PS,
                30.0,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
        ),
        "expected_continuation_frames_are_61": (
            EXPECTED_CONTINUATION_FRAMES
            == 61
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
        "R2_CHECKPOINT_CONTINUATION_TO_50PS_PREPARED"
        if accepted
        else
        "R2_CHECKPOINT_CONTINUATION_PREPARATION_REQUIRES_REVIEW"
    )

    required_next_step = (
        "RUN_R2_20_TO_50PS_CHECKPOINT_CONTINUATION"
        if accepted
        else
        "REVIEW_R2_CHECKPOINT_CONTINUATION_PREPARATION_FAILURES"
    )

    checksum_paths = (
        SOURCE_TPR,
        EXTENDED_TPR,
        SOURCE_CHECKPOINT,
        COPIED_CHECKPOINT,
        SOURCE_FINAL_GRO,
    )

    write_csv(
        CHECKSUMS_CSV,
        [
            {
                "file": relative(path),
                "sha256": sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for path in checksum_paths
        ],
    )

    summary = {
        "decision": decision,
        "convert_tpr_return_code": (
            convert.returncode
        ),
        "source_TPR_dump_return_code": (
            source_dump.returncode
        ),
        "extended_TPR_dump_return_code": (
            extended_dump.returncode
        ),
        "checkpoint_dump_return_code": (
            checkpoint_dump.returncode
        ),
        "source_atoms": (
            source_atoms
        ),
        "extended_atoms": (
            extended_atoms
        ),
        "source_nsteps": (
            source_nsteps
        ),
        "extended_nsteps": (
            extended_nsteps
        ),
        "source_dt_ps": (
            source_dt
        ),
        "extended_dt_ps": (
            extended_dt
        ),
        "source_time_ps": (
            source_nsteps
            * source_dt
        ),
        "target_time_ps": (
            extended_nsteps
            * extended_dt
        ),
        "remaining_steps": (
            REMAINING_STEPS
        ),
        "remaining_time_ps": (
            EXTENSION_TIME_PS
        ),
        "XTC_stride": (
            extended_stride
        ),
        "XTC_interval_ps": (
            XTC_INTERVAL_PS
        ),
        "expected_continuation_frames": (
            EXPECTED_CONTINUATION_FRAMES
        ),
        "checkpoint_step": (
            checkpoint_step
        ),
        "checkpoint_time_ps": (
            checkpoint_time
        ),
        "checkpoint_coordinate_entries": (
            checkpoint_coordinate_entries
        ),
        "checkpoint_velocity_entries": (
            checkpoint_velocity_entries
        ),
        "TPR_differences_beyond_nsteps": (
            len(diff_lines)
        ),
        "checkpoint_bitwise_copy": (
            checkpoint_bitwise_copy
        ),
        "source_checkpoint_sha256": (
            source_checkpoint_hash
        ),
        "copied_checkpoint_sha256": (
            copied_checkpoint_hash
        ),
        "failed_gates": (
            " | ".join(
                failed_gates
            )
        ),
        "checkpoint_continuation_execution_authorized": (
            accepted
        ),
        "velocity_regeneration_authorized": False,
        "source_0_to_20ps_rerun_authorized": False,
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

    run_contract = {
        "decision": decision,
        "execution_authorized": (
            accepted
        ),
        "purpose": (
            "Test whether R2 lumen occupancy approaches "
            "a plateau or continues to decline."
        ),
        "input_extended_TPR": relative(
            EXTENDED_TPR
        ),
        "input_exact_20ps_checkpoint": relative(
            COPIED_CHECKPOINT
        ),
        "source_checkpoint_sha256": (
            source_checkpoint_hash
        ),
        "copied_checkpoint_sha256": (
            copied_checkpoint_hash
        ),
        "output_directory": relative(
            EXECUTION_OUTPUT_ROOT
        ),
        "output_deffnm": relative(
            EXECUTION_DEFFNM
        ),
        "command": [
            "gmx",
            "mdrun",
            "-s",
            relative(EXTENDED_TPR),
            "-cpi",
            relative(COPIED_CHECKPOINT),
            "-deffnm",
            relative(EXECUTION_DEFFNM),
            "-noappend",
            "-ntmpi",
            "1",
            "-ntomp",
            "4",
        ],
        "source_step": SOURCE_NSTEPS,
        "source_time_ps": SOURCE_TIME_PS,
        "target_step": TARGET_NSTEPS,
        "target_time_ps": TARGET_TIME_PS,
        "remaining_steps": REMAINING_STEPS,
        "remaining_time_ps": EXTENSION_TIME_PS,
        "expected_continuation_frames": (
            EXPECTED_CONTINUATION_FRAMES
        ),
        "expected_continuation_time_range_ps": [
            SOURCE_TIME_PS,
            TARGET_TIME_PS,
        ],
        "checkpoint_continuation": True,
        "velocity_regeneration": False,
        "thermostat_state_regeneration": False,
        "source_0_to_20ps_rerun": False,
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

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed
        in gates.items()
    )

    REPORT_MD.write_text(
        f"""# R2 Checkpoint Continuation to 50 ps — Preparation

## Purpose

The R2 20 ps frozen-solute screen remains formally nonstationary.
The final trajectory windows show flattening, so an exact checkpoint
continuation is prepared to test plateau formation versus continued
depletion.

The completed 0–20 ps trajectory will not be repeated.

## Source state

- Source atoms:
  **{source_atoms}**
- Source steps:
  **{source_nsteps}**
- Source time:
  **{source_nsteps * source_dt:.6f} ps**
- Source checkpoint step/time:
  **{checkpoint_step}/{checkpoint_time:.6f} ps**
- Source failed gate:
  **{source_failed_gates}**

## Extended TPR

- Extended atoms:
  **{extended_atoms}**
- Extended steps:
  **{extended_nsteps}**
- Target time:
  **{extended_nsteps * extended_dt:.6f} ps**
- Remaining steps:
  **{REMAINING_STEPS}**
- Remaining time:
  **{EXTENSION_TIME_PS:.6f} ps**
- XTC stride:
  **{extended_stride}**
- Expected continuation frames:
  **{EXPECTED_CONTINUATION_FRAMES}**
- TPR differences beyond normalized `nsteps`:
  **{len(diff_lines)}**

## Checkpoint integrity

- Coordinate entries:
  **{checkpoint_coordinate_entries}**
- Velocity entries:
  **{checkpoint_velocity_entries}**
- Source checkpoint SHA256:
  `{source_checkpoint_hash}`
- Copied checkpoint SHA256:
  `{copied_checkpoint_hash}`
- Bitwise-identical copy:
  **{'YES' if checkpoint_bitwise_copy else 'NO'}**

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Checkpoint-continuation execution authorized:
  **{'YES' if accepted else 'NO'}**
- Velocity regeneration authorized:
  **NO**
- Source 0–20 ps rerun authorized:
  **NO**
- Long mobile MD authorized:
  **NO**
- Multitemperature MD authorized:
  **NO**
- QM recalculation authorized:
  **NO**
- Required next step:
  `{required_next_step}`
""",
        encoding="utf-8",
    )

    print(
        "Day023 R2 checkpoint continuation to 50 ps "
        "preparation completed."
    )

    print(
        "Convert-tpr / source dump / extended dump / "
        "checkpoint dump return codes: "
        f"{convert.returncode}/"
        f"{source_dump.returncode}/"
        f"{extended_dump.returncode}/"
        f"{checkpoint_dump.returncode}"
    )

    print(
        "Source / extended atoms: "
        f"{source_atoms}/{extended_atoms}"
    )

    print(
        "Source / extended nsteps: "
        f"{source_nsteps}/{extended_nsteps}"
    )

    print(
        "Source / target time: "
        f"{source_nsteps * source_dt:.6f}/"
        f"{extended_nsteps * extended_dt:.6f} ps"
    )

    print(
        "Checkpoint step / time: "
        f"{checkpoint_step}/"
        f"{checkpoint_time:.6f} ps"
    )

    print(
        "Remaining steps / time: "
        f"{REMAINING_STEPS}/"
        f"{EXTENSION_TIME_PS:.6f} ps"
    )

    print(
        "XTC stride / expected continuation frames: "
        f"{extended_stride}/"
        f"{EXPECTED_CONTINUATION_FRAMES}"
    )

    print(
        "Checkpoint coordinate / velocity entries: "
        f"{checkpoint_coordinate_entries}/"
        f"{checkpoint_velocity_entries}"
    )

    print(
        "TPR differences beyond nsteps: "
        f"{len(diff_lines)}"
    )

    print(
        "Checkpoint bitwise copy: "
        f"{'YES' if checkpoint_bitwise_copy else 'NO'}"
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
        "Checkpoint-continuation execution authorized: "
        f"{'YES' if accepted else 'NO'}"
    )

    print(
        "Velocity regeneration authorized: NO"
    )

    print(
        "Source 0-to-20 ps rerun authorized: NO"
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
        EXTENDED_TPR,
        COPIED_CHECKPOINT,
        SOURCE_TPR_DUMP,
        EXTENDED_TPR_DUMP,
        CHECKPOINT_DUMP,
        PHYSICAL_DIFF,
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
            "R2 checkpoint-continuation preparation "
            "requires review."
        )


if __name__ == "__main__":
    main()
