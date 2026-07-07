#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

PROTOCOL = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

STAGE = "07_nvt_unrestrained_25ps"

STATIC_ROOT = (
    PROTOCOL
    / "static_validation"
)

STAGE_STATIC_ROOT = (
    STATIC_ROOT
    / STAGE
)

MASTER_CSV = (
    STATIC_ROOT
    / "mobile_release_static_validation.csv"
)

PREPARATION_SUMMARY = (
    STAGE_STATIC_ROOT
    / f"{STAGE}_preparation_summary.csv"
)

TPR = (
    STAGE_STATIC_ROOT
    / f"{STAGE}.tpr"
)

PROCESSED_TOP = (
    STAGE_STATIC_ROOT
    / f"{STAGE}_processed.top"
)

GROMPP_LOG = (
    STAGE_STATIC_ROOT
    / f"{STAGE}_grompp.log"
)


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
    )


def read_single_row(path: Path) -> dict[str, str]:
    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle))

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected exactly one row in {path}, found {len(rows)}"
        )

    return rows[0]


def main() -> None:
    required = (
        MASTER_CSV,
        PREPARATION_SUMMARY,
        TPR,
        PROCESSED_TOP,
        GROMPP_LOG,
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
            "Missing or empty required files:\n"
            + "\n".join(
                str(path)
                for path in missing
            )
        )

    preparation = read_single_row(
        PREPARATION_SUMMARY
    )

    if (
        preparation.get(
            "preparation_decision",
            "",
        )
        .strip()
        .upper()
        != "PASS"
    ):
        raise RuntimeError(
            "Stage07 preparation decision is not PASS"
        )

    if (
        preparation.get(
            "grompp_pass",
            "",
        )
        .strip()
        .lower()
        not in {"true", "1", "yes"}
    ):
        raise RuntimeError(
            "Stage07 static grompp is not PASS"
        )

    if int(
        preparation.get(
            "position_restraint_entries",
            "-1",
        )
    ) != 0:
        raise RuntimeError(
            "Stage07 contains active position restraints"
        )

    with MASTER_CSV.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fields = list(reader.fieldnames or [])

    if not fields:
        raise RuntimeError(
            "Static-validation master CSV has no header"
        )

    new_row = {
        "stage": STAGE,
        "kind": "nvt",
        "restraint_k": "0",
        "grompp_return_code": "0",
        "grompp_pass": "True",
        "position_restraint_entries": "0",
        "expected_position_restraint_entries": "0",
        "restraint_entry_validation": "True",
        "tpr_path": relative(TPR),
        "processed_topology_path": relative(
            PROCESSED_TOP
        ),
        "grompp_log": relative(GROMPP_LOG),
    }

    for field in new_row:
        if field not in fields:
            fields.append(field)

    rows = [
        row
        for row in rows
        if row.get(
            "stage",
            "",
        ).strip() != STAGE
    ]

    rows.append(new_row)

    temporary = MASTER_CSV.with_suffix(
        ".csv.tmp"
    )

    with temporary.open(
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

    temporary.replace(MASTER_CSV)

    with MASTER_CSV.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        registered_rows = [
            row
            for row in csv.DictReader(handle)
            if row.get(
                "stage",
                "",
            ).strip() == STAGE
        ]

    if len(registered_rows) != 1:
        raise RuntimeError(
            "Stage07 was not uniquely registered"
        )

    print(
        "Day022 Stage07 static-validation registration completed."
    )
    print(
        f"Stage: {STAGE}"
    )
    print(
        "Matching master-validation rows: 1"
    )
    print(
        "Static grompp: PASS"
    )
    print(
        "Active position restraints: 0/0"
    )
    print(
        f"TPR: {relative(TPR)}"
    )
    print(
        "Stage07 execution authorization: PASS"
    )


if __name__ == "__main__":
    main()
