#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import math
import re
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

RUN_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol/"
    "execution/02_nvt_k10000_1ps"
)

STAGE = "02_nvt_k10000_1ps"

SUMMARY = (
    RUN_ROOT
    / f"{STAGE}_summary.csv"
)

REPORT = (
    RUN_ROOT
    / f"{STAGE.upper()}_DAY021.md"
)

MDRUN_LOG = (
    RUN_ROOT
    / f"{STAGE}.log"
)

CONSOLE_LOG = (
    RUN_ROOT
    / f"{STAGE}_mdrun_console.log"
)

RECONCILIATION_REPORT = (
    RUN_ROOT
    / "STAGE02_FALSE_POSITIVE_RECONCILIATION_DAY021.md"
)

OUTPUT_MANIFEST = (
    RUN_ROOT
    / f"{STAGE}_output_manifest.csv"
)

NONFINITE_PATTERN = re.compile(
    r"(?<![A-Za-z])"
    r"(?:nan|[-+]?inf(?:inity)?)"
    r"(?![A-Za-z])",
    re.IGNORECASE,
)

HARMLESS_EPSILON_RF = re.compile(
    r"^\s*epsilon-rf\s*=\s*"
    r"(?:inf|infinity)\s*$",
    re.IGNORECASE,
)

SERIOUS_PATTERN = re.compile(
    r"LINCS WARNING|"
    r"Fatal error|"
    r"SETTLE.*(?:error|failed|cannot)|"
    r"SHAKE.*(?:failed|converge)|"
    r"constraint.*(?:error|failed)",
    re.IGNORECASE,
)


def relative(path: Path) -> str:
    try:
        return str(
            path.resolve().relative_to(
                PROJECT_ROOT
            )
        )
    except ValueError:
        return str(path.resolve())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            block = handle.read(
                1024 * 1024
            )

            if not block:
                break

            digest.update(block)

    return digest.hexdigest()


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {
        "true",
        "1",
        "yes",
        "pass",
    }


def require_files() -> None:
    required = (
        SUMMARY,
        REPORT,
        MDRUN_LOG,
        CONSOLE_LOG,
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
            "Missing or empty Stage02 files:\n"
            + "\n".join(
                str(path)
                for path in missing
            )
        )


def audit_text_logs() -> tuple[
    list[tuple[str, int, str]],
    list[tuple[str, int, str]],
]:
    nonfinite_hits = []
    serious_hits = []

    for path in (
        MDRUN_LOG,
        CONSOLE_LOG,
    ):
        lines = path.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines()

        for line_number, line in enumerate(
            lines,
            start=1,
        ):
            if SERIOUS_PATTERN.search(line):
                serious_hits.append(
                    (
                        path.name,
                        line_number,
                        line.strip(),
                    )
                )

            if NONFINITE_PATTERN.search(line):
                nonfinite_hits.append(
                    (
                        path.name,
                        line_number,
                        line.strip(),
                    )
                )

    return nonfinite_hits, serious_hits


def audit_xvg_files() -> tuple[int, int]:
    files = sorted(
        RUN_ROOT.glob("*.xvg")
    )

    if not files:
        raise RuntimeError(
            "No Stage02 XVG files were found"
        )

    numeric_rows = 0
    nonfinite_values = 0

    for path in files:
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

            numeric_rows += 1

            for token in line.split():
                try:
                    value = float(token)
                except ValueError:
                    continue

                if not math.isfinite(value):
                    nonfinite_values += 1

    return numeric_rows, nonfinite_values


def read_summary() -> tuple[
    dict[str, str],
    list[str],
]:
    with SUMMARY.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = list(
            reader.fieldnames or []
        )

    if len(rows) != 1:
        raise RuntimeError(
            "Stage02 summary must contain "
            "exactly one data row"
        )

    return rows[0], fieldnames


def write_summary(
    row: dict[str, str],
    fieldnames: list[str],
) -> None:
    additions = (
        "reconciliation_status",
        "reconciliation_reason",
        "rerun_required",
    )

    for field in additions:
        if field not in fieldnames:
            fieldnames.append(field)

    with SUMMARY.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerow(
            {
                field: row.get(field, "")
                for field in fieldnames
            }
        )


def update_report() -> None:
    text = REPORT.read_text(
        encoding="utf-8",
        errors="replace",
    )

    text = text.replace(
        "- Decision: **BLOCKED**",
        "- Decision: **PASS**",
        1,
    )

    text = re.sub(
        r"- Instability signatures:.*",
        "- Instability signatures: none",
        text,
        count=1,
    )

    reconciliation_section = """
## False-positive reconciliation

The original automated detector matched the normal GROMACS
parameter record `epsilon-rf = inf`. With PME electrostatics this
is metadata, not a non-finite runtime observable. All extracted
temperature, potential-energy, total-energy, and pressure values
were finite, and no LINCS, SETTLE, SHAKE, constraint, or fatal
errors were detected.

- Scientific decision: **PASS**
- Rerun required: **NO**
"""

    if (
        "## False-positive reconciliation"
        not in text
    ):
        text = (
            text.rstrip()
            + "\n\n"
            + reconciliation_section.strip()
            + "\n"
        )

    REPORT.write_text(
        text,
        encoding="utf-8",
    )


def write_reconciliation_report(
    nonfinite_hits,
    serious_hits,
    numeric_rows: int,
    numeric_nonfinite: int,
) -> None:
    with RECONCILIATION_REPORT.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day021 Stage02 False-Positive Reconciliation\n\n"
        )

        handle.write(
            f"- Stage: `{STAGE}`\n"
        )

        handle.write(
            f"- Textual non-finite matches: "
            f"{len(nonfinite_hits)}\n"
        )

        for filename, line_number, line in (
            nonfinite_hits
        ):
            handle.write(
                f"  - `{filename}:{line_number}`: "
                f"`{line}`\n"
            )

        handle.write(
            f"- Serious instability matches: "
            f"{len(serious_hits)}\n"
        )

        handle.write(
            f"- Numeric XVG rows audited: "
            f"{numeric_rows}\n"
        )

        handle.write(
            f"- Non-finite XVG values: "
            f"{numeric_nonfinite}\n\n"
        )

        handle.write(
            "The only textual match was the normal "
            "`epsilon-rf = inf` parameter record. "
            "No runtime observable was non-finite.\n\n"
        )

        handle.write(
            "- Reconciled decision: **PASS**\n"
        )

        handle.write(
            "- Rerun required: **NO**\n"
        )


def rewrite_manifest() -> None:
    files = sorted(
        path
        for path in RUN_ROOT.iterdir()
        if (
            path.is_file()
            and path != OUTPUT_MANIFEST
            and not path.name.endswith(
                "_wrapper.log"
            )
        )
    )

    rows = [
        {
            "path": relative(path),
            "size_bytes": (
                path.stat().st_size
            ),
            "sha256": sha256(path),
        }
        for path in files
    ]

    with OUTPUT_MANIFEST.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "path",
                "size_bytes",
                "sha256",
            ],
        )

        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    require_files()

    row, fieldnames = read_summary()

    required_pass_fields = (
        "atom_count_pass",
        "finite_coordinate_box_pass",
        "temperature_pass",
        "incremental_displacement_pass",
        "cumulative_displacement_pass",
        "velocity_pass",
    )

    failed_fields = [
        field
        for field in required_pass_fields
        if not parse_bool(
            row.get(field, "")
        )
    ]

    if failed_fields:
        raise RuntimeError(
            "Stage02 cannot be reconciled because "
            "scientific acceptance checks failed:\n"
            + "\n".join(failed_fields)
        )

    nonfinite_hits, serious_hits = (
        audit_text_logs()
    )

    numeric_rows, numeric_nonfinite = (
        audit_xvg_files()
    )

    harmful_nonfinite_hits = [
        hit
        for hit in nonfinite_hits
        if not HARMLESS_EPSILON_RF.match(
            hit[2]
        )
    ]

    if serious_hits:
        raise RuntimeError(
            "Serious instability signatures remain; "
            "Stage02 was not reclassified"
        )

    if harmful_nonfinite_hits:
        raise RuntimeError(
            "A non-finite value other than the known "
            "epsilon-rf metadata record remains"
        )

    if numeric_nonfinite != 0:
        raise RuntimeError(
            "Non-finite numeric XVG values remain"
        )

    if len(nonfinite_hits) != 1:
        raise RuntimeError(
            "Expected exactly one harmless textual "
            "non-finite record"
        )

    row["decision"] = "PASS"
    row["instability_signatures"] = ""
    row["blocked_reasons"] = ""
    row["reconciliation_status"] = "PASS"
    row["reconciliation_reason"] = (
        "False positive from normal GROMACS "
        "parameter record epsilon-rf = inf"
    )
    row["rerun_required"] = "False"

    write_summary(
        row,
        fieldnames,
    )

    update_report()

    write_reconciliation_report(
        nonfinite_hits,
        serious_hits,
        numeric_rows,
        numeric_nonfinite,
    )

    rewrite_manifest()

    print(
        "Day021 Stage02 false-positive "
        "reconciliation completed."
    )

    print(
        f"Textual non-finite matches: "
        f"{len(nonfinite_hits)}"
    )

    print(
        "Harmless epsilon-rf metadata matches: "
        f"{len(nonfinite_hits)}"
    )

    print(
        f"Serious instability signatures: "
        f"{len(serious_hits)}"
    )

    print(
        f"Numeric XVG rows audited: "
        f"{numeric_rows}"
    )

    print(
        f"Non-finite XVG values: "
        f"{numeric_nonfinite}"
    )

    print(
        "Stage02 reconciled decision: PASS"
    )

    print(
        "Stage02 rerun required: NO"
    )

    print(
        f"Wrote: {relative(RECONCILIATION_REPORT)}"
    )


if __name__ == "__main__":
    main()
