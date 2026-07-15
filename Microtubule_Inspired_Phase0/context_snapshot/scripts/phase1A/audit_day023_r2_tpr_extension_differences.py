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

SOURCE_DUMP = (
    STAGE
    / "r2_20ps_source_tpr_dump.txt"
)

EXTENDED_DUMP = (
    STAGE
    / "r2_50ps_extended_tpr_dump.txt"
)

EXISTING_DIFF = (
    STAGE
    / "r2_source_vs_extended_tpr_physical_diff.txt"
)

MISMATCH_CSV = (
    STAGE
    / "r2_source_vs_extended_tpr_mismatch_records.csv"
)

SUMMARY_JSON = (
    STAGE
    / "r2_source_vs_extended_tpr_difference_audit.json"
)

REPORT_MD = (
    STAGE
    / "R2_SOURCE_VS_EXTENDED_TPR_DIFFERENCE_AUDIT_DAY023.md"
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


def normalize_expected_nsteps(
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

    return line


def record_family(
    line: str,
) -> str:
    stripped = line.strip()

    if re.match(
        r"^x\[\s*\d+\s*\]",
        stripped,
    ):
        return "coordinates"

    if re.match(
        r"^v\[\s*\d+\s*\]",
        stripped,
    ):
        return "velocities"

    if re.match(
        r"^(box|box-rel|boxv)\[",
        stripped,
    ):
        return "box_state"

    if re.match(
        r"^[A-Za-z0-9_.-]+\s*=",
        stripped,
    ):
        return "scalar_parameter"

    return "other"


def parameter_name(
    line: str,
) -> str:
    match = re.match(
        r"^\s*([A-Za-z0-9_.-]+)\s*=",
        line,
    )

    if match is None:
        return ""

    return match.group(1)


def write_csv(
    path: Path,
    rows: list[dict[str, Any]],
) -> None:
    fields = [
        "line_number",
        "family",
        "parameter",
        "source_line",
        "extended_line",
    ]

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
            writer.writerow(row)


def main() -> None:
    for path in (
        SOURCE_DUMP,
        EXTENDED_DUMP,
        EXISTING_DIFF,
    ):
        require_file(path)

    source_lines = SOURCE_DUMP.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    extended_lines = EXTENDED_DUMP.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    existing_diff_lines = EXISTING_DIFF.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()

    source_normalized = [
        normalize_expected_nsteps(line)
        for line in source_lines
    ]

    extended_normalized = [
        normalize_expected_nsteps(line)
        for line in extended_lines
    ]

    mismatch_rows: list[dict[str, Any]] = []

    shared_length = min(
        len(source_normalized),
        len(extended_normalized),
    )

    for index in range(shared_length):
        source_line = source_normalized[index]
        extended_line = extended_normalized[index]

        if source_line == extended_line:
            continue

        source_family = record_family(
            source_line
        )

        extended_family = record_family(
            extended_line
        )

        family = (
            source_family
            if source_family == extended_family
            else (
                source_family
                + " -> "
                + extended_family
            )
        )

        source_parameter = parameter_name(
            source_line
        )

        extended_parameter = parameter_name(
            extended_line
        )

        parameter = (
            source_parameter
            if source_parameter == extended_parameter
            else (
                source_parameter
                + " -> "
                + extended_parameter
            ).strip(" ->")
        )

        mismatch_rows.append(
            {
                "line_number": index + 1,
                "family": family,
                "parameter": parameter,
                "source_line": source_lines[index],
                "extended_line": extended_lines[index],
            }
        )

    if len(source_normalized) != len(extended_normalized):
        longer = (
            source_lines
            if len(source_lines) > len(extended_lines)
            else extended_lines
        )

        side = (
            "source_only"
            if len(source_lines) > len(extended_lines)
            else "extended_only"
        )

        for index in range(
            shared_length,
            len(longer),
        ):
            mismatch_rows.append(
                {
                    "line_number": index + 1,
                    "family": side,
                    "parameter": "",
                    "source_line": (
                        source_lines[index]
                        if side == "source_only"
                        else ""
                    ),
                    "extended_line": (
                        extended_lines[index]
                        if side == "extended_only"
                        else ""
                    ),
                }
            )

    write_csv(
        MISMATCH_CSV,
        mismatch_rows,
    )

    unified_removed = [
        line
        for line in existing_diff_lines
        if (
            line.startswith("-")
            and not line.startswith("---")
        )
    ]

    unified_added = [
        line
        for line in existing_diff_lines
        if (
            line.startswith("+")
            and not line.startswith("+++")
        )
    ]

    family_counts: dict[str, int] = {}

    for row in mismatch_rows:
        family = str(row["family"])
        family_counts[family] = (
            family_counts.get(family, 0)
            + 1
        )

    changed_parameters = sorted(
        {
            str(row["parameter"])
            for row in mismatch_rows
            if str(row["parameter"])
        }
    )

    coordinate_mismatches = sum(
        row["family"] == "coordinates"
        for row in mismatch_rows
    )

    velocity_mismatches = sum(
        row["family"] == "velocities"
        for row in mismatch_rows
    )

    box_state_mismatches = sum(
        row["family"] == "box_state"
        for row in mismatch_rows
    )

    scalar_mismatches = sum(
        row["family"] == "scalar_parameter"
        for row in mismatch_rows
    )

    summary = {
        "source_dump_lines": len(source_lines),
        "extended_dump_lines": len(extended_lines),
        "unified_diff_total_lines": len(
            existing_diff_lines
        ),
        "unified_diff_removed_records": len(
            unified_removed
        ),
        "unified_diff_added_records": len(
            unified_added
        ),
        "semantic_mismatch_records_after_nsteps_normalization": len(
            mismatch_rows
        ),
        "coordinate_mismatches": (
            coordinate_mismatches
        ),
        "velocity_mismatches": (
            velocity_mismatches
        ),
        "box_state_mismatches": (
            box_state_mismatches
        ),
        "scalar_parameter_mismatches": (
            scalar_mismatches
        ),
        "changed_parameters": (
            changed_parameters
        ),
        "family_counts": (
            family_counts
        ),
        "continuation_execution_authorized_by_this_audit": False,
        "required_next_step": (
            "CLASSIFY_TPR_EXTENSION_RESIDUAL_DIFFERENCES"
        ),
    }

    SUMMARY_JSON.write_text(
        json.dumps(
            summary,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    mismatch_lines = []

    for row in mismatch_rows:
        mismatch_lines.append(
            f"- Line {row['line_number']} "
            f"[{row['family']}] "
            f"`{row['parameter']}`\n"
            f"  - source: `{row['source_line']}`\n"
            f"  - extended: `{row['extended_line']}`"
        )

    REPORT_MD.write_text(
        f"""# R2 Source-versus-Extended TPR Difference Audit

## Scope

This audit compares the complete `gmx dump` representations of the
20 ps source TPR and the 50 ps extended TPR.

The expected `nsteps` change is normalized before comparison. No TPR,
checkpoint, trajectory, or molecular-dynamics result is modified.

## Diff accounting

- Unified-diff total lines:
  **{len(existing_diff_lines)}**
- Actual removed records:
  **{len(unified_removed)}**
- Actual added records:
  **{len(unified_added)}**
- Semantic mismatch records after normalizing `nsteps`:
  **{len(mismatch_rows)}**

## State integrity

- Coordinate mismatches:
  **{coordinate_mismatches}**
- Velocity mismatches:
  **{velocity_mismatches}**
- Box-state mismatches:
  **{box_state_mismatches}**
- Scalar-parameter mismatches:
  **{scalar_mismatches}**

## Changed parameters

{chr(10).join(f'- `{name}`' for name in changed_parameters) if changed_parameters else '- NONE'}

## Exact mismatch records

{chr(10).join(mismatch_lines) if mismatch_lines else '- NONE'}

## Current authorization

- Checkpoint-continuation execution authorized:
  **NO**
- Required next step:
  `CLASSIFY_TPR_EXTENSION_RESIDUAL_DIFFERENCES`
""",
        encoding="utf-8",
    )

    print(
        "Day023 R2 TPR-extension difference audit completed."
    )

    print(
        "Unified diff total / removed / added lines: "
        f"{len(existing_diff_lines)}/"
        f"{len(unified_removed)}/"
        f"{len(unified_added)}"
    )

    print(
        "Semantic mismatches after nsteps normalization: "
        f"{len(mismatch_rows)}"
    )

    print(
        "Coordinate / velocity / box-state mismatches: "
        f"{coordinate_mismatches}/"
        f"{velocity_mismatches}/"
        f"{box_state_mismatches}"
    )

    print(
        "Scalar-parameter mismatches: "
        f"{scalar_mismatches}"
    )

    print(
        "Changed parameters: "
        + (
            "NONE"
            if not changed_parameters
            else " | ".join(
                changed_parameters
            )
        )
    )

    for row in mismatch_rows:
        print(
            "MISMATCH "
            f"line={row['line_number']} "
            f"family={row['family']} "
            f"parameter={row['parameter']}"
        )

        print(
            "  SOURCE:   "
            f"{row['source_line']}"
        )

        print(
            "  EXTENDED: "
            f"{row['extended_line']}"
        )

    print(
        "Checkpoint-continuation execution authorized: NO"
    )

    print(
        "Required next step: "
        "CLASSIFY_TPR_EXTENSION_RESIDUAL_DIFFERENCES"
    )

    print(
        f"Wrote: {MISMATCH_CSV.relative_to(ROOT)}"
    )

    print(
        f"Wrote: {SUMMARY_JSON.relative_to(ROOT)}"
    )

    print(
        f"Wrote: {REPORT_MD.relative_to(ROOT)}"
    )


if __name__ == "__main__":
    main()
