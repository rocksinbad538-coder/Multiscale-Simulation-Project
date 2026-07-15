#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

STAGE = (
    ROOT
    / "runs/phase1A/day023_confinement_design/"
    "18_r2_nvt_50ps_checkpoint_continuation_preparation"
)

OUTPUT_ROOT = (
    STAGE
    / "difference_classification"
)

SOURCE_DUMP = (
    STAGE
    / "r2_20ps_source_tpr_dump.txt"
)

EXTENDED_DUMP = (
    STAGE
    / "r2_50ps_extended_tpr_dump.txt"
)

MISMATCH_CSV = (
    STAGE
    / "r2_source_vs_extended_tpr_mismatch_records.csv"
)

AUDIT_JSON = (
    STAGE
    / "r2_source_vs_extended_tpr_difference_audit.json"
)

PREPARATION_SUMMARY = (
    STAGE
    / "r2_50ps_continuation_preparation_summary.csv"
)

PREPARATION_GATES = (
    STAGE
    / "r2_50ps_continuation_preparation_gates.csv"
)

ORIGINAL_RUN_CONTRACT = (
    STAGE
    / "r2_50ps_continuation_run_contract.json"
)

CLASSIFICATION_SUMMARY = (
    OUTPUT_ROOT
    / "r2_tpr_extension_difference_classification_summary.csv"
)

CLASSIFICATION_GATES = (
    OUTPUT_ROOT
    / "r2_tpr_extension_difference_classification_gates.csv"
)

AUTHORIZED_RUN_CONTRACT = (
    OUTPUT_ROOT
    / "r2_50ps_continuation_authorized_run_contract.json"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_TPR_EXTENSION_DIFFERENCE_CLASSIFICATION_DAY023.md"
)

EXPECTED_ORIGINAL_DECISION = (
    "R2_CHECKPOINT_CONTINUATION_PREPARATION_REQUIRES_REVIEW"
)

EXPECTED_ORIGINAL_FAILED_GATE = (
    "TPRs_differ_only_by_nsteps"
)

AUTHORIZED_DECISION = (
    "R2_CHECKPOINT_CONTINUATION_TO_50PS_AUTHORIZED"
)

EXPECTED_SOURCE_NSTEPS = 40000
EXPECTED_EXTENDED_NSTEPS = 100000
EXPECTED_ATOMS = 68332


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
            f"Expected one row in {path}; found {len(rows)}"
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
            f"Could not parse {key!r} from TPR dump."
        )

    return match.group(1)


def is_dump_filename_header(line: str) -> bool:
    stripped = line.strip()

    return bool(
        re.fullmatch(
            r".+\.tpr:",
            stripped,
            flags=re.IGNORECASE,
        )
    )


def normalize_expected_run_length(
    line: str,
) -> str:
    if re.match(
        r"^\s*nsteps\s*=",
        line,
    ):
        return re.sub(
            r"(nsteps\s*=\s*)\S+",
            r"\1<EXPECTED_RUN_LENGTH>",
            line,
            count=1,
        )

    return line.rstrip()


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        SOURCE_DUMP,
        EXTENDED_DUMP,
        MISMATCH_CSV,
        AUDIT_JSON,
        PREPARATION_SUMMARY,
        PREPARATION_GATES,
        ORIGINAL_RUN_CONTRACT,
    ):
        require_file(required)

    source_text = SOURCE_DUMP.read_text(
        encoding="utf-8",
        errors="replace",
    )

    extended_text = EXTENDED_DUMP.read_text(
        encoding="utf-8",
        errors="replace",
    )

    source_lines = source_text.splitlines()
    extended_lines = extended_text.splitlines()

    if (
        not source_lines
        or not extended_lines
    ):
        raise RuntimeError(
            "One or both TPR dumps are empty."
        )

    mismatch_rows = read_csv_rows(
        MISMATCH_CSV
    )

    audit = json.loads(
        AUDIT_JSON.read_text(
            encoding="utf-8"
        )
    )

    preparation_summary = read_single_csv_row(
        PREPARATION_SUMMARY
    )

    preparation_gate_rows = read_csv_rows(
        PREPARATION_GATES
    )

    original_contract = json.loads(
        ORIGINAL_RUN_CONTRACT.read_text(
            encoding="utf-8"
        )
    )

    failed_preparation_gates = [
        row.get("gate", "")
        for row in preparation_gate_rows
        if not parse_bool(
            row.get(
                "pass",
                "false",
            )
        )
    ]

    source_header = source_lines[0]
    extended_header = extended_lines[0]

    source_body = [
        normalize_expected_run_length(line)
        for line in source_lines[1:]
    ]

    extended_body = [
        normalize_expected_run_length(line)
        for line in extended_lines[1:]
    ]

    residual_body_mismatches = []

    shared_length = min(
        len(source_body),
        len(extended_body),
    )

    for index in range(shared_length):
        if source_body[index] != extended_body[index]:
            residual_body_mismatches.append(
                {
                    "body_line": index + 2,
                    "source": source_lines[
                        index + 1
                    ],
                    "extended": extended_lines[
                        index + 1
                    ],
                }
            )

    if len(source_body) != len(extended_body):
        residual_body_mismatches.append(
            {
                "body_line": "length",
                "source": len(source_body),
                "extended": len(extended_body),
            }
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

    single_mismatch = (
        len(mismatch_rows) == 1
    )

    mismatch = (
        mismatch_rows[0]
        if single_mismatch
        else {}
    )

    mismatch_is_filename_header = (
        single_mismatch
        and mismatch.get(
            "line_number",
            "",
        )
        == "1"
        and mismatch.get(
            "family",
            "",
        )
        == "other"
        and not mismatch.get(
            "parameter",
            "",
        ).strip()
        and is_dump_filename_header(
            mismatch.get(
                "source_line",
                "",
            )
        )
        and is_dump_filename_header(
            mismatch.get(
                "extended_line",
                "",
            )
        )
    )

    audit_changed_parameters = (
        audit.get(
            "changed_parameters",
            [],
        )
    )

    audit_family_counts = (
        audit.get(
            "family_counts",
            {},
        )
    )

    gates = {
        "original_preparation_is_expected_review_state": (
            preparation_summary.get(
                "decision"
            )
            == EXPECTED_ORIGINAL_DECISION
        ),
        "original_preparation_failed_exactly_one_gate": (
            failed_preparation_gates
            == [
                EXPECTED_ORIGINAL_FAILED_GATE
            ]
        ),
        "audit_found_exactly_one_semantic_mismatch": (
            int(
                audit.get(
                    "semantic_mismatch_records_after_nsteps_normalization",
                    -1,
                )
            )
            == 1
        ),
        "audit_found_zero_coordinate_mismatches": (
            int(
                audit.get(
                    "coordinate_mismatches",
                    -1,
                )
            )
            == 0
        ),
        "audit_found_zero_velocity_mismatches": (
            int(
                audit.get(
                    "velocity_mismatches",
                    -1,
                )
            )
            == 0
        ),
        "audit_found_zero_box_state_mismatches": (
            int(
                audit.get(
                    "box_state_mismatches",
                    -1,
                )
            )
            == 0
        ),
        "audit_found_zero_scalar_parameter_mismatches": (
            int(
                audit.get(
                    "scalar_parameter_mismatches",
                    -1,
                )
            )
            == 0
        ),
        "audit_found_no_changed_parameters": (
            audit_changed_parameters == []
        ),
        "audit_mismatch_family_is_only_other": (
            audit_family_counts == {
                "other": 1
            }
        ),
        "mismatch_record_is_only_dump_filename_header": (
            mismatch_is_filename_header
        ),
        "source_first_line_is_dump_filename_header": (
            is_dump_filename_header(
                source_header
            )
        ),
        "extended_first_line_is_dump_filename_header": (
            is_dump_filename_header(
                extended_header
            )
        ),
        "normalized_TPR_bodies_are_exactly_equal": (
            len(
                residual_body_mismatches
            )
            == 0
        ),
        "source_nsteps_is_40000": (
            source_nsteps
            == EXPECTED_SOURCE_NSTEPS
        ),
        "extended_nsteps_is_100000": (
            extended_nsteps
            == EXPECTED_EXTENDED_NSTEPS
        ),
        "source_and_extended_atom_counts_are_68332": (
            source_atoms
            == EXPECTED_ATOMS
            and extended_atoms
            == EXPECTED_ATOMS
        ),
        "checkpoint_copy_is_bitwise_identical": (
            parse_bool(
                preparation_summary.get(
                    "checkpoint_bitwise_copy",
                    "false",
                )
            )
        ),
        "checkpoint_step_is_40000": (
            int(
                float(
                    preparation_summary.get(
                        "checkpoint_step",
                        "-1",
                    )
                )
            )
            == EXPECTED_SOURCE_NSTEPS
        ),
        "checkpoint_time_is_20ps": (
            abs(
                float(
                    preparation_summary.get(
                        "checkpoint_time_ps",
                        "nan",
                    )
                )
                - 20.0
            )
            <= 1.0e-9
        ),
        "original_contract_did_not_authorize_execution": (
            not bool(
                original_contract.get(
                    "execution_authorized",
                    True,
                )
            )
        ),
        "original_contract_forbids_velocity_regeneration": (
            not bool(
                original_contract.get(
                    "velocity_regeneration",
                    True,
                )
            )
        ),
        "original_contract_forbids_source_rerun": (
            not bool(
                original_contract.get(
                    "source_0_to_20ps_rerun",
                    True,
                )
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
        AUTHORIZED_DECISION
        if accepted
        else
        "R2_TPR_EXTENSION_DIFFERENCE_CLASSIFICATION_REQUIRES_REVIEW"
    )

    required_next_step = (
        "RUN_R2_20_TO_50PS_CHECKPOINT_CONTINUATION"
        if accepted
        else
        "REVIEW_R2_TPR_EXTENSION_CLASSIFICATION_FAILURES"
    )

    classification_summary = {
        "decision": decision,
        "original_preparation_decision": (
            preparation_summary.get(
                "decision",
                "",
            )
        ),
        "original_failed_gates": (
            " | ".join(
                failed_preparation_gates
            )
        ),
        "semantic_mismatch_count": (
            len(mismatch_rows)
        ),
        "classified_benign_header_mismatch_count": (
            1
            if mismatch_is_filename_header
            else 0
        ),
        "residual_TPR_body_mismatch_count_after_header_and_nsteps_normalization": (
            len(
                residual_body_mismatches
            )
        ),
        "source_dump_header": (
            source_header
        ),
        "extended_dump_header": (
            extended_header
        ),
        "source_nsteps": (
            source_nsteps
        ),
        "extended_nsteps": (
            extended_nsteps
        ),
        "source_atoms": (
            source_atoms
        ),
        "extended_atoms": (
            extended_atoms
        ),
        "difference_classification": (
            "BENIGN_GMX_DUMP_FILENAME_HEADER"
            if mismatch_is_filename_header
            else "UNRESOLVED"
        ),
        "physical_TPR_differences_beyond_nsteps": (
            len(
                residual_body_mismatches
            )
        ),
        "failed_classification_gates": (
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
        CLASSIFICATION_SUMMARY,
        [
            classification_summary
        ],
    )

    write_csv(
        CLASSIFICATION_GATES,
        [
            {
                "gate": name,
                "pass": passed,
            }
            for name, passed
            in gates.items()
        ],
    )

    authorized_contract = dict(
        original_contract
    )

    authorized_contract.update(
        {
            "decision": decision,
            "execution_authorized": accepted,
            "authorization_basis": (
                "The only residual mismatch after normalizing "
                "nsteps is the nonphysical gmx dump filename header."
            ),
            "difference_classification": (
                "BENIGN_GMX_DUMP_FILENAME_HEADER"
                if accepted
                else "UNRESOLVED"
            ),
            "physical_TPR_differences_beyond_nsteps": (
                len(
                    residual_body_mismatches
                )
            ),
            "source_dump_header": (
                source_header
            ),
            "extended_dump_header": (
                extended_header
            ),
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
    )

    AUTHORIZED_RUN_CONTRACT.write_text(
        json.dumps(
            authorized_contract,
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
        f"""# R2 TPR-Extension Difference Classification

## Finding

The source and extended `gmx dump` outputs contain one residual
difference after normalization of `nsteps`.

That difference is the first-line filename header:

- Source:
  `{source_header}`
- Extended:
  `{extended_header}`

This line identifies the file supplied to `gmx dump`; it is not part of
the simulation input record, coordinates, velocities, box state, force
field, thermostat state, or integration state.

## State comparison

- Coordinate mismatches:
  **{audit.get('coordinate_mismatches', 'UNKNOWN')}**
- Velocity mismatches:
  **{audit.get('velocity_mismatches', 'UNKNOWN')}**
- Box-state mismatches:
  **{audit.get('box_state_mismatches', 'UNKNOWN')}**
- Scalar-parameter mismatches:
  **{audit.get('scalar_parameter_mismatches', 'UNKNOWN')}**
- Changed parameters:
  **{'NONE' if not audit_changed_parameters else ' | '.join(audit_changed_parameters)}**
- Residual body mismatches after removing the filename header and
  normalizing `nsteps`:
  **{len(residual_body_mismatches)}**

## Run-length change

- Source `nsteps`:
  **{source_nsteps}**
- Extended `nsteps`:
  **{extended_nsteps}**
- Source atoms:
  **{source_atoms}**
- Extended atoms:
  **{extended_atoms}**

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Difference classification:
  **{'BENIGN_GMX_DUMP_FILENAME_HEADER' if accepted else 'UNRESOLVED'}**
- Physical TPR differences beyond `nsteps`:
  **{len(residual_body_mismatches)}**
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
        "Day023 R2 TPR-extension difference "
        "classification completed."
    )

    print(
        "Original preparation decision: "
        f"{preparation_summary.get('decision', '')}"
    )

    print(
        "Original failed gates: "
        + (
            "NONE"
            if not failed_preparation_gates
            else " | ".join(
                failed_preparation_gates
            )
        )
    )

    print(
        "Semantic mismatch count: "
        f"{len(mismatch_rows)}"
    )

    print(
        "Difference classification: "
        + (
            "BENIGN_GMX_DUMP_FILENAME_HEADER"
            if mismatch_is_filename_header
            else "UNRESOLVED"
        )
    )

    print(
        "Residual TPR-body mismatches after header "
        "and nsteps normalization: "
        f"{len(residual_body_mismatches)}"
    )

    print(
        "Source / extended nsteps: "
        f"{source_nsteps}/"
        f"{extended_nsteps}"
    )

    print(
        "Source / extended atoms: "
        f"{source_atoms}/"
        f"{extended_atoms}"
    )

    print(
        f"Decision: {decision}"
    )

    print(
        "Failed classification gates: "
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
        CLASSIFICATION_SUMMARY,
        CLASSIFICATION_GATES,
        AUTHORIZED_RUN_CONTRACT,
        REPORT_MD,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not accepted:
        raise RuntimeError(
            "R2 TPR-extension difference classification "
            "requires review."
        )


if __name__ == "__main__":
    main()
