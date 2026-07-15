#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

STAGE = "04_nvt_k100_2ps"

RUN_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol/"
    "execution/04_nvt_k100_2ps"
)

STRUCTURAL_SUMMARY = (
    RUN_ROOT
    / f"{STAGE}_structural_summary.csv"
)

STAGE_SUMMARY = (
    RUN_ROOT
    / f"{STAGE}_summary.csv"
)

STAGE_REPORT = (
    RUN_ROOT
    / f"{STAGE.upper()}_DAY021.md"
)

RECONCILIATION_REPORT = (
    RUN_ROOT
    / "STAGE04_STRUCTURAL_PASS_RECONCILIATION_DAY021.md"
)

OUTPUT_MANIFEST = (
    RUN_ROOT
    / f"{STAGE}_output_manifest.csv"
)


def relative(path: Path) -> str:
    try:
        return str(
            path.resolve().relative_to(PROJECT_ROOT)
        )
    except ValueError:
        return str(path.resolve())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            block = handle.read(1024 * 1024)

            if not block:
                break

            digest.update(block)

    return digest.hexdigest()


def read_single_row(
    path: Path,
) -> tuple[dict[str, str], list[str]]:
    if not path.exists() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty file: {path}"
        )

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fields = list(reader.fieldnames or [])

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected exactly one row in {path}"
        )

    return rows[0], fields


def number(
    row: dict[str, str],
    field: str,
) -> float:
    try:
        return float(row[field])
    except (KeyError, ValueError) as error:
        raise RuntimeError(
            f"Invalid field: {field}"
        ) from error


def validate_structural_result(
    row: dict[str, str],
) -> None:
    if (
        row.get("structural_screen", "")
        .strip()
        .upper()
        != "STABLE_CANDIDATE"
    ):
        raise RuntimeError(
            "Stage04 is not a STABLE_CANDIDATE"
        )

    zero_fields = (
        "harmful_nonfinite_matches",
        "serious_instability_matches",
        "numeric_xvg_nonfinite_values",
    )

    failures = [
        field
        for field in zero_fields
        if number(row, field) != 0.0
    ]

    checks = {
        "maximum HBN bond length <= 0.20 nm": (
            number(
                row,
                "HBN_bond_length_max_nm",
            )
            <= 0.20
        ),
        "q99 bond deviation <= 0.015 nm": (
            number(
                row,
                "HBN_bond_deviation_q99_nm",
            )
            <= 0.015
        ),
        "maximum bond deviation <= 0.03 nm": (
            number(
                row,
                "HBN_bond_deviation_max_nm",
            )
            <= 0.03
        ),
        "q95 angle deviation <= 20 degrees": (
            number(
                row,
                "HBN_angle_deviation_q95_deg",
            )
            <= 20.0
        ),
        "maximum angle deviation <= 60 degrees": (
            number(
                row,
                "HBN_angle_deviation_max_deg",
            )
            <= 60.0
        ),
    }

    failures.extend(
        label
        for label, passed in checks.items()
        if not passed
    )

    if failures:
        raise RuntimeError(
            "Stage04 structural acceptance failed:\n"
            + "\n".join(failures)
        )


def update_stage_summary() -> None:
    row, fields = read_single_row(
        STAGE_SUMMARY
    )

    if (
        row.get("decision", "")
        .strip()
        .upper()
        != "PASS"
    ):
        raise RuntimeError(
            "Stage04 operational decision is not PASS"
        )

    additions = {
        "structural_reconciliation_status": "PASS",
        "structural_screen": "STABLE_CANDIDATE",
        "structural_reconciliation_reason": (
            "HBN bond and angle distributions remain "
            "within accepted structural thresholds at k=100"
        ),
        "rerun_required": "False",
        "next_action": "05_nvt_unrestrained_2ps",
    }

    row.update(additions)

    for field in additions:
        if field not in fields:
            fields.append(field)

    with STAGE_SUMMARY.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
        )

        writer.writeheader()
        writer.writerow(
            {
                field: row.get(field, "")
                for field in fields
            }
        )


def update_stage_report() -> None:
    if not STAGE_REPORT.exists():
        raise RuntimeError(
            f"Missing stage report: {STAGE_REPORT}"
        )

    text = STAGE_REPORT.read_text(
        encoding="utf-8",
        errors="replace",
    )

    section = """
## Structural acceptance

The k = 100 restrained NVT stage was accepted after a
bonded-network diagnostic.

- Structural screen: **STABLE_CANDIDATE**
- HBN bond-length range: 0.13851715–0.15524175 nm
- Maximum bond-equilibrium deviation: 0.00973175 nm
- Maximum angle-equilibrium deviation: 3.7094 degrees
- Serious instability signatures: 0
- Harmful non-finite runtime values: 0

Scientific decision: **PASS**.

- Stage04 rerun required: **NO**
- Stage05 unrestrained pilot authorized: **YES**
"""

    if "## Structural acceptance" not in text:
        text = (
            text.rstrip()
            + "\n\n"
            + section.strip()
            + "\n"
        )

    STAGE_REPORT.write_text(
        text,
        encoding="utf-8",
    )


def write_reconciliation_report() -> None:
    RECONCILIATION_REPORT.write_text(
        """# Day021 Stage04 Structural Pass Reconciliation

- Operational decision: **PASS**
- Structural screen: **STABLE_CANDIDATE**
- Bonded HBN failure: **NO**
- Numerical instability: **NO**
- Rerun required: **NO**
- Authorized next stage: `05_nvt_unrestrained_2ps`

The HBN bonded geometry remained stable after reducing the
position-restraint strength to 100 kJ mol^-1 nm^-2.
""",
        encoding="utf-8",
    )


def rewrite_manifest() -> None:
    files = sorted(
        path
        for path in RUN_ROOT.iterdir()
        if (
            path.is_file()
            and path != OUTPUT_MANIFEST
            and not path.name.endswith("_wrapper.log")
        )
    )

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

        for path in files:
            writer.writerow(
                {
                    "path": relative(path),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )


def main() -> None:
    structural_row, _ = read_single_row(
        STRUCTURAL_SUMMARY
    )

    validate_structural_result(
        structural_row
    )

    update_stage_summary()
    update_stage_report()
    write_reconciliation_report()
    rewrite_manifest()

    print(
        "Day021 Stage04 structural reconciliation completed."
    )
    print(
        "Stage04 operational decision: PASS"
    )
    print(
        "Stage04 structural screen: STABLE_CANDIDATE"
    )
    print(
        "Stage04 scientific decision: PASS"
    )
    print(
        "Stage04 rerun required: NO"
    )
    print(
        "Stage05 unrestrained pilot authorized: YES"
    )
    print(
        f"Wrote: {relative(RECONCILIATION_REPORT)}"
    )


if __name__ == "__main__":
    main()
