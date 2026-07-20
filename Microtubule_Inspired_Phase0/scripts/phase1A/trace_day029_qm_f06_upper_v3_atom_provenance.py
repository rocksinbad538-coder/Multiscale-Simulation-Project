#!/usr/bin/env python3
"""
Trace provenance of the critical QM_F06 UPPER V3 atoms.

Targets:
- local indices 11, 14, 17
- N11, H14, H17
- known UPPER cap identifiers
- restored source identifiers used during V2/V3 construction

This script performs text-level provenance discovery only.
It does not modify geometries or authorize QM execution.
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

SEARCH_ROOTS = [
    ROOT / "scripts/phase1A",
    ROOT / "runs/phase1A/day028_qm_f06_upper_transferability",
]

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day029_qm_f06_upper_v3_provenance"
)

TARGET_PATTERNS = {
    "local_index_11": re.compile(
        r"(?<!\d)11(?!\d)|N11|atom[_\s-]*11",
        re.IGNORECASE,
    ),
    "local_index_14": re.compile(
        r"(?<!\d)14(?!\d)|H14|atom[_\s-]*14",
        re.IGNORECASE,
    ),
    "local_index_17": re.compile(
        r"(?<!\d)17(?!\d)|H17|atom[_\s-]*17",
        re.IGNORECASE,
    ),
    "upper_caps": re.compile(
        r"HCAP(?::|V2:|V3:)?UPPER",
        re.IGNORECASE,
    ),
    "restored_atom_A_UPPER_10_4": re.compile(
        r"A:UPPER:10:4",
        re.IGNORECASE,
    ),
    "restored_atom_A_UPPER_8_4": re.compile(
        r"A:UPPER:8:4",
        re.IGNORECASE,
    ),
}

TEXT_SUFFIXES = {
    ".py",
    ".json",
    ".csv",
    ".txt",
    ".log",
    ".xyz",
    ".pdb",
    ".inp",
    ".md",
    ".tex",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def iter_candidate_files():
    for root in SEARCH_ROOTS:
        if not root.exists():
            continue

        for path in root.rglob("*"):
            if (
                path.is_file()
                and path.suffix.lower() in TEXT_SUFFIXES
            ):
                yield path


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    matches = []

    for path in sorted(iter_candidate_files()):
        try:
            text = path.read_text(
                encoding="utf-8",
                errors="replace",
            )
        except OSError:
            continue

        lines = text.splitlines()

        for line_number, line in enumerate(
            lines,
            start=1,
        ):
            labels = [
                label
                for label, pattern in TARGET_PATTERNS.items()
                if pattern.search(line)
            ]

            if not labels:
                continue

            start = max(0, line_number - 3)
            end = min(len(lines), line_number + 2)

            context = "\n".join(
                f"{index + 1}: {lines[index]}"
                for index in range(start, end)
            )

            matches.append(
                {
                    "file": str(path.relative_to(ROOT)),
                    "line_number": line_number,
                    "matched_labels": labels,
                    "line": line,
                    "context": context,
                    "file_sha256": sha256(path),
                }
            )

    summary = {
        "decision": (
            "QM_F06_UPPER_V3_CRITICAL_ATOM_"
            "PROVENANCE_TRACE_GENERATED"
        ),
        "generated_utc": datetime.now(
            timezone.utc
        ).isoformat(),
        "search_roots": [
            str(path.relative_to(ROOT))
            for path in SEARCH_ROOTS
            if path.exists()
        ],
        "targets": list(TARGET_PATTERNS),
        "match_count": len(matches),
        "matches": matches,
        "authorization": {
            "boundary_reconstruction_authorized": False,
            "qm_execution_authorized": False,
            "geometry_reference_accepted": False,
            "RESP_execution_authorized": False,
        },
    }

    report_path = (
        OUTPUT_DIR
        / "QM_F06_UPPER_V3_CRITICAL_ATOM_PROVENANCE.json"
    )

    report_path.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    text_report = (
        OUTPUT_DIR
        / "QM_F06_UPPER_V3_CRITICAL_ATOM_PROVENANCE.txt"
    )

    sections = [
        "=" * 78,
        "QM_F06 UPPER V3 CRITICAL ATOM PROVENANCE",
        "=" * 78,
        f"Matches: {len(matches)}",
        "",
    ]

    for number, match in enumerate(matches, start=1):
        sections.extend(
            [
                f"[{number}] {match['file']}:{match['line_number']}",
                "Matched: "
                + ", ".join(match["matched_labels"]),
                match["context"],
                "-" * 78,
            ]
        )

    text_report.write_text(
        "\n".join(sections) + "\n",
        encoding="utf-8",
    )

    print("=" * 78)
    print("QM_F06 UPPER V3 CRITICAL ATOM PROVENANCE")
    print("=" * 78)
    print("Matches:", len(matches))
    print("JSON:", report_path)
    print("TXT:", text_report)
    print("QM execution authorized: False")


if __name__ == "__main__":
    main()
