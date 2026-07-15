#!/usr/bin/env python3
"""
Summarize LAMMPS thermo output for Phase 0 runs.

The script extracts thermo tables containing:
Step Atoms Temp c_twater PotEng KinEng TotEng Press c_msdwater[4]

It writes:
- one combined CSV with all thermo records
- one summary CSV with initial/final/min/max values per run
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable


THERMO_HEADER_START = "Step"
EXPECTED_COLUMNS = [
    "Step",
    "Atoms",
    "Temp",
    "c_twater",
    "PotEng",
    "KinEng",
    "TotEng",
    "Press",
    "c_msdwater[4]",
]


def is_float_like(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def parse_log(log_path: Path, run_label: str) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    in_table = False
    columns: list[str] = []

    with log_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.strip().split()
            if not parts:
                continue

            if parts[0] == THERMO_HEADER_START and "c_twater" in parts:
                in_table = True
                columns = parts
                continue

            if in_table:
                if parts[0] in {"Loop", "Performance:", "MPI", "Section", "Nlocal:", "Histogram:"}:
                    in_table = False
                    continue

                if len(parts) != len(columns):
                    continue

                if not is_float_like(parts[0]):
                    continue

                record: dict[str, float | str] = {"run": run_label, "log_file": str(log_path)}
                for col, val in zip(columns, parts):
                    record[col] = float(val)
                rows.append(record)

    return rows


def summarize(records: list[dict[str, float | str]]) -> dict[str, float | str]:
    if not records:
        return {}

    run_label = str(records[0]["run"])
    log_file = str(records[0]["log_file"])

    first = records[0]
    last = records[-1]

    numeric_cols = [
        "Step",
        "Atoms",
        "Temp",
        "c_twater",
        "PotEng",
        "KinEng",
        "TotEng",
        "Press",
        "c_msdwater[4]",
    ]

    summary: dict[str, float | str] = {
        "run": run_label,
        "log_file": log_file,
        "n_records": len(records),
        "initial_step": float(first["Step"]),
        "final_step": float(last["Step"]),
    }

    for col in numeric_cols:
        values = [float(r[col]) for r in records if col in r]
        summary[f"{col}_initial"] = float(first[col])
        summary[f"{col}_final"] = float(last[col])
        summary[f"{col}_min"] = min(values)
        summary[f"{col}_max"] = max(values)

    return summary


def write_csv(path: Path, rows: Iterable[dict[str, float | str]]) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)

    if not rows:
        path.write_text("", encoding="utf-8")
        return

    fieldnames: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)

    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--log", action="append", nargs=2, metavar=("LABEL", "PATH"), required=True)
    parser.add_argument("--records-csv", required=True)
    parser.add_argument("--summary-csv", required=True)
    args = parser.parse_args()

    all_records: list[dict[str, float | str]] = []
    summaries: list[dict[str, float | str]] = []

    for label, path_str in args.log:
        log_path = Path(path_str)
        records = parse_log(log_path, label)
        if not records:
            print(f"WARNING: no thermo records found in {log_path}")
            continue
        all_records.extend(records)
        summaries.append(summarize(records))

    write_csv(Path(args.records_csv), all_records)
    write_csv(Path(args.summary_csv), summaries)

    print("LAMMPS thermo summary complete.")
    print(f"Records CSV: {args.records_csv}")
    print(f"Summary CSV: {args.summary_csv}")
    print(f"Runs parsed: {len(summaries)}")


if __name__ == "__main__":
    main()
