#!/usr/bin/env python3

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

STAGE_ROOT = (
    ROOT
    / "runs/phase1A/day023_confinement_design/"
    "13_r2_topology_static_scan"
)

OUTPUT_ROOT = (
    STAGE_ROOT
    / "postrun_instability_audit"
)

CONSOLE_LOG = (
    STAGE_ROOT
    / "r2_static_single_point_mdrun_console.log"
)

MDRUN_LOG = (
    STAGE_ROOT
    / "r2_static_single_point.log"
)

SUMMARY_CSV = (
    STAGE_ROOT
    / "r2_topology_static_scan_summary.csv"
)

GATES_CSV = (
    STAGE_ROOT
    / "r2_topology_static_scan_gates.csv"
)

MATCHES_CSV = (
    OUTPUT_ROOT
    / "r2_static_instability_pattern_matches.csv"
)

AUDIT_SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_static_instability_audit_summary.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_STATIC_INSTABILITY_SIGNATURE_AUDIT_DAY023.md"
)

EXPECTED_DECISION = (
    "R2_TOPOLOGY_AND_STATIC_CAP_WATER_MODEL_VALIDATED"
)

BROAD_LEGACY_PATTERN = re.compile(
    r"Fatal error|Segmentation fault|LINCS WARNING|"
    r"constraint warning|nan",
    flags=re.IGNORECASE,
)

STRICT_PATTERNS = {
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
    "standalone_nan": re.compile(
        r"\bnan\b",
        flags=re.IGNORECASE,
    ),
}


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
            f"No CSV rows found in {path}"
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
    fields: list[str],
) -> None:
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
        CONSOLE_LOG,
        MDRUN_LOG,
        SUMMARY_CSV,
        GATES_CSV,
    ):
        require_file(required)

    summary = read_single_csv_row(
        SUMMARY_CSV
    )

    gate_rows = read_csv_rows(
        GATES_CSV
    )

    scan_rows: list[dict[str, Any]] = []

    broad_match_count = 0
    strict_match_count = 0
    broad_only_count = 0

    for source_path in (
        CONSOLE_LOG,
        MDRUN_LOG,
    ):
        lines = source_path.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines()

        for line_number, line in enumerate(
            lines,
            start=1,
        ):
            broad_match = (
                BROAD_LEGACY_PATTERN.search(line)
            )

            strict_names = [
                name
                for name, pattern
                in STRICT_PATTERNS.items()
                if pattern.search(line)
            ]

            if broad_match is not None:
                broad_match_count += 1

            if strict_names:
                strict_match_count += len(
                    strict_names
                )

            if (
                broad_match is not None
                and not strict_names
            ):
                broad_only_count += 1

            if (
                broad_match is None
                and not strict_names
            ):
                continue

            classification = (
                "STRICT_INSTABILITY"
                if strict_names
                else "BROAD_PATTERN_FALSE_POSITIVE"
            )

            scan_rows.append(
                {
                    "source_file": (
                        source_path.name
                    ),
                    "line_number": (
                        line_number
                    ),
                    "classification": (
                        classification
                    ),
                    "broad_pattern_match": (
                        broad_match.group(0)
                        if broad_match is not None
                        else ""
                    ),
                    "strict_patterns": (
                        " | ".join(
                            strict_names
                        )
                    ),
                    "line_text": line,
                }
            )

    write_csv(
        MATCHES_CSV,
        scan_rows,
        [
            "source_file",
            "line_number",
            "classification",
            "broad_pattern_match",
            "strict_patterns",
            "line_text",
        ],
    )

    decision = summary.get(
        "decision",
        "",
    )

    all_core_gates_passed = all(
        parse_bool(
            row.get(
                "pass",
                "false",
            )
        )
        for row in gate_rows
    )

    water_only_em_authorized = parse_bool(
        summary.get(
            "water_only_energy_minimization_authorized",
            "false",
        )
    )

    core_decision_valid = (
        decision
        == EXPECTED_DECISION
    )

    strict_scan_clean = (
        strict_match_count == 0
    )

    accepted = (
        core_decision_valid
        and all_core_gates_passed
        and strict_scan_clean
        and water_only_em_authorized
    )

    audit_decision = (
        "R2_STATIC_INSTABILITY_SCAN_CONFIRMED_CLEAN"
        if accepted
        else
        "R2_STATIC_INSTABILITY_AUDIT_REQUIRES_REVIEW"
    )

    required_next_step = (
        "RUN_R2_WATER_ONLY_ENERGY_MINIMIZATION"
        if accepted
        else
        "REVIEW_R2_STATIC_INSTABILITY_AUDIT"
    )

    audit_summary = {
        "decision": audit_decision,
        "core_gate_decision": decision,
        "core_gate_decision_valid": (
            core_decision_valid
        ),
        "core_gate_count": (
            len(gate_rows)
        ),
        "all_core_gates_passed": (
            all_core_gates_passed
        ),
        "legacy_broad_pattern_matches": (
            broad_match_count
        ),
        "strict_instability_matches": (
            strict_match_count
        ),
        "broad_pattern_false_positives": (
            broad_only_count
        ),
        "water_only_energy_minimization_authorized": (
            accepted
        ),
        "short_frozen_solute_NVT_authorized": False,
        "long_mobile_MD_authorized": False,
        "multitemperature_MD_authorized": False,
        "QM_recalculation_authorized": False,
        "MD_or_static_rerun_required": False,
        "required_next_step": (
            required_next_step
        ),
    }

    write_csv(
        AUDIT_SUMMARY_CSV,
        [audit_summary],
        list(
            audit_summary.keys()
        ),
    )

    broad_only_lines = [
        (
            f"- `{row['source_file']}:{row['line_number']}` — "
            f"`{row['broad_pattern_match']}` in: "
            f"`{row['line_text']}`"
        )
        for row in scan_rows
        if (
            row["classification"]
            ==
            "BROAD_PATTERN_FALSE_POSITIVE"
        )
    ]

    strict_lines = [
        (
            f"- `{row['source_file']}:{row['line_number']}` — "
            f"{row['strict_patterns']}: "
            f"`{row['line_text']}`"
        )
        for row in scan_rows
        if (
            row["classification"]
            ==
            "STRICT_INSTABILITY"
        )
    ]

    REPORT_MD.write_text(
        f"""# R2 Static Instability-Signature Audit

## Scope

This audit evaluates the completed R2 zero-step single-point run. It
does not rerun `grompp`, `mdrun`, minimization, or molecular dynamics.

The earlier terminal wrapper used the unrestricted expression `nan`,
which can match ordinary text containing those three characters. The
strict audit uses the bounded numerical token `\\bnan\\b` and explicit
GROMACS failure signatures.

## Core Gate 2B state

- Decision:
  **{decision}**
- Expected decision:
  **{EXPECTED_DECISION}**
- Core gates:
  **{len(gate_rows)}**
- All core gates passed:
  **{'YES' if all_core_gates_passed else 'NO'}**

## Pattern results

- Legacy broad-pattern matches:
  **{broad_match_count}**
- Strict instability matches:
  **{strict_match_count}**
- Broad-pattern false positives:
  **{broad_only_count}**

### Broad-pattern-only matches

{chr(10).join(broad_only_lines) if broad_only_lines else '- NONE'}

### Strict instability matches

{chr(10).join(strict_lines) if strict_lines else '- NONE'}

## Decision

- Audit decision:
  **{audit_decision}**
- Static or MD rerun required:
  **NO**
- Water-only energy minimization authorized:
  **{'YES' if accepted else 'NO'}**
- Short frozen-solute NVT authorized:
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
        "Day023 R2 static instability-signature "
        "audit completed."
    )

    print(
        "Core Gate 2B decision: "
        f"{decision}"
    )

    print(
        "Core gates / all passed: "
        f"{len(gate_rows)} / "
        f"{'YES' if all_core_gates_passed else 'NO'}"
    )

    print(
        "Legacy broad-pattern matches: "
        f"{broad_match_count}"
    )

    print(
        "Strict instability matches: "
        f"{strict_match_count}"
    )

    print(
        "Broad-pattern false positives: "
        f"{broad_only_count}"
    )

    for row in scan_rows:
        print(
            "MATCH: "
            f"{row['classification']} | "
            f"{row['source_file']}:"
            f"{row['line_number']} | "
            f"{row['line_text']}"
        )

    print(
        f"Decision: {audit_decision}"
    )

    print(
        "Static or MD rerun required: NO"
    )

    print(
        "Water-only energy minimization authorized: "
        f"{'YES' if accepted else 'NO'}"
    )

    print(
        "Short frozen-solute NVT authorized: NO"
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
        "Wrote: "
        + str(
            MATCHES_CSV.relative_to(ROOT)
        )
    )

    print(
        "Wrote: "
        + str(
            AUDIT_SUMMARY_CSV.relative_to(ROOT)
        )
    )

    print(
        "Wrote: "
        + str(
            REPORT_MD.relative_to(ROOT)
        )
    )

    if not accepted:
        raise RuntimeError(
            "R2 static instability audit requires review."
        )


if __name__ == "__main__":
    main()
