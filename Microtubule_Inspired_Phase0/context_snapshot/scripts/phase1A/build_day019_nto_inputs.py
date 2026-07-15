#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import re
import shutil
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_ROOT = (
    PROJECT_ROOT
    / "runs/phase1A/day016_md_bath_extraction"
)

SELECTION_CSV = (
    DATA_ROOT
    / "day019_vacuum_reference_analysis"
    / "NTO_SELECTION_DAY019.csv"
)

EMBEDDED_ROOT = (
    DATA_ROOT
    / "orca_embedding_pilot_inputs"
)

VACUUM_ROOT = (
    DATA_ROOT
    / "day019_vacuum_reference_inputs"
)

OUTPUT_ROOT = (
    DATA_ROOT
    / "day019_nto_inputs"
)

MANIFEST_CSV = (
    OUTPUT_ROOT
    / "NTO_INPUT_MANIFEST_DAY019.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "NTO_INPUT_AUDIT_DAY019.md"
)

NTO_STATES = [1, 2]
NTO_THRESHOLD = "1e-4"

TDDFT_START_RE = re.compile(
    r"^\s*%tddft\s*$",
    flags=re.IGNORECASE,
)

BLOCK_END_RE = re.compile(
    r"^\s*end\s*$",
    flags=re.IGNORECASE,
)

POINTCHARGE_RE = re.compile(
    r'^\s*%pointcharges\s+"([^"]+)"\s*$',
    flags=re.IGNORECASE | re.MULTILINE,
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def locate_tddft_block(
    lines: list[str],
) -> tuple[int, int]:
    starts = [
        index
        for index, line in enumerate(lines)
        if TDDFT_START_RE.match(line)
    ]

    if len(starts) != 1:
        raise RuntimeError(
            "Expected exactly one %tddft block; "
            f"found {len(starts)}."
        )

    start = starts[0]
    stop = None

    for index in range(start + 1, len(lines)):
        if BLOCK_END_RE.match(lines[index]):
            stop = index
            break

    if stop is None:
        raise RuntimeError(
            "Could not locate end of %tddft block."
        )

    return start, stop


def validate_source_input(
    text: str,
    calculation_type: str,
    source_path: Path,
) -> dict:
    lines = text.splitlines()

    tddft_start, tddft_stop = locate_tddft_block(
        lines
    )

    tddft_text = "\n".join(
        lines[tddft_start : tddft_stop + 1]
    )

    if not re.search(
        r"(?mi)^\s*nroots\s+10\s*$",
        tddft_text,
    ):
        raise RuntimeError(
            f"Missing 'nroots 10' in {source_path}"
        )

    if not re.search(
        r"(?mi)^\s*tda\s+true\s*$",
        tddft_text,
    ):
        raise RuntimeError(
            f"Missing 'tda true' in {source_path}"
        )

    forbidden = [
        "donto",
        "ntostates",
        "ntothresh",
    ]

    for keyword in forbidden:
        if re.search(
            rf"(?mi)^\s*{keyword}\b",
            tddft_text,
        ):
            raise RuntimeError(
                f"Source already contains {keyword}: "
                f"{source_path}"
            )

    pointcharge_lines = [
        (index, match.group(1))
        for index, line in enumerate(lines)
        if (
            match := POINTCHARGE_RE.match(line)
        )
    ]

    if calculation_type == "embedded":
        if len(pointcharge_lines) != 1:
            raise RuntimeError(
                "Embedded source must contain exactly "
                f"one %pointcharges line: {source_path}"
            )
    elif calculation_type == "vacuum_reference":
        if pointcharge_lines:
            raise RuntimeError(
                "Vacuum source unexpectedly contains "
                f"%pointcharges: {source_path}"
            )
    else:
        raise RuntimeError(
            f"Unsupported calculation type: "
            f"{calculation_type}"
        )

    return {
        "lines": lines,
        "tddft_start": tddft_start,
        "tddft_stop": tddft_stop,
        "pointcharge_lines": pointcharge_lines,
    }


def build_nto_input(
    source_text: str,
) -> str:
    lines = source_text.splitlines()
    _, tddft_stop = locate_tddft_block(lines)

    nto_lines = [
        "  DoNTO true",
        "  NTOStates 1,2",
        f"  NTOThresh {NTO_THRESHOLD}",
    ]

    target_lines = (
        lines[:tddft_stop]
        + nto_lines
        + lines[tddft_stop:]
    )

    return "\n".join(target_lines).rstrip() + "\n"


def validate_target_input(
    source_text: str,
    target_text: str,
    calculation_type: str,
    target_path: Path,
) -> None:
    source_lines = source_text.splitlines()
    target_lines = target_text.splitlines()

    _, source_stop = locate_tddft_block(
        source_lines
    )
    _, target_stop = locate_tddft_block(
        target_lines
    )

    expected_inserted = [
        "  DoNTO true",
        "  NTOStates 1,2",
        f"  NTOThresh {NTO_THRESHOLD}",
    ]

    observed_inserted = (
        target_lines[source_stop:target_stop]
    )

    if observed_inserted != expected_inserted:
        raise RuntimeError(
            "Unexpected lines added to TDDFT block "
            f"in {target_path}: {observed_inserted}"
        )

    reconstructed_source = (
        target_lines[:source_stop]
        + target_lines[target_stop:]
    )

    if reconstructed_source != source_lines:
        raise RuntimeError(
            "Target differs from source by more than "
            f"the three NTO lines: {target_path}"
        )

    required_patterns = {
        "DoNTO": r"(?mi)^\s*DoNTO\s+true\s*$",
        "NTOStates": (
            r"(?mi)^\s*NTOStates\s+1\s*,\s*2\s*$"
        ),
        "NTOThresh": (
            rf"(?mi)^\s*NTOThresh\s+"
            rf"{re.escape(NTO_THRESHOLD)}\s*$"
        ),
    }

    for label, pattern in required_patterns.items():
        if len(re.findall(pattern, target_text)) != 1:
            raise RuntimeError(
                f"Invalid {label} directive in "
                f"{target_path}"
            )

    pointcharge_matches = [
        match
        for line in target_text.splitlines()
        if (match := POINTCHARGE_RE.match(line))
    ]

    if calculation_type == "embedded":
        if len(pointcharge_matches) != 1:
            raise RuntimeError(
                "Embedded NTO target must retain exactly "
                f"one %pointcharges directive: {target_path}"
            )
    else:
        if pointcharge_matches:
            raise RuntimeError(
                "Vacuum NTO target contains "
                f"%pointcharges: {target_path}"
            )


def source_paths(
    calculation_type: str,
    source_job: str,
) -> tuple[Path, Path]:
    if calculation_type == "embedded":
        source_dir = EMBEDDED_ROOT / source_job
    elif calculation_type == "vacuum_reference":
        source_dir = VACUUM_ROOT / source_job
    else:
        raise RuntimeError(
            f"Unsupported calculation type: "
            f"{calculation_type}"
        )

    source_input = (
        source_dir
        / f"{source_job}.inp"
    )

    return source_dir, source_input


def main() -> None:
    if not SELECTION_CSV.is_file():
        raise SystemExit(
            f"Missing NTO selection file: "
            f"{SELECTION_CSV}"
        )

    selection = pd.read_csv(SELECTION_CSV)

    required_columns = {
        "calculation_type",
        "frame",
        "cluster",
        "job",
        "tracked_root",
        "reason",
    }

    missing = sorted(
        required_columns - set(selection.columns)
    )

    if missing:
        raise SystemExit(
            f"NTO selection is missing columns: "
            f"{missing}"
        )

    if len(selection) != 8:
        raise SystemExit(
            "Expected exactly 8 selected NTO cases; "
            f"found {len(selection)}."
        )

    if selection["job"].duplicated().any():
        duplicates = selection.loc[
            selection["job"].duplicated(
                keep=False
            ),
            "job",
        ].tolist()

        raise SystemExit(
            f"Duplicate NTO jobs selected: "
            f"{duplicates}"
        )

    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    rows: list[dict] = []

    for record in selection.itertuples(
        index=False
    ):
        calculation_type = str(
            record.calculation_type
        )
        source_job = str(record.job)

        source_dir, source_input = source_paths(
            calculation_type,
            source_job,
        )

        if not source_input.is_file():
            raise SystemExit(
                f"Missing source input: "
                f"{source_input}"
            )

        source_text = source_input.read_text(
            encoding="utf-8"
        )

        source_info = validate_source_input(
            source_text,
            calculation_type,
            source_input,
        )

        target_job = (
            f"{source_job}_nto_s1_s2"
        )
        target_dir = (
            OUTPUT_ROOT
            / target_job
        )
        target_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        target_input = (
            target_dir
            / f"{target_job}.inp"
        )

        target_text = build_nto_input(
            source_text
        )

        target_input.write_text(
            target_text,
            encoding="utf-8",
        )

        validate_target_input(
            source_text=source_text,
            target_text=target_text,
            calculation_type=calculation_type,
            target_path=target_input,
        )

        copied_pc = ""
        copied_pc_sha256 = ""

        if calculation_type == "embedded":
            pc_filename = (
                source_info[
                    "pointcharge_lines"
                ][0][1]
            )
            source_pc = (
                source_dir
                / pc_filename
            )

            if not source_pc.is_file():
                raise SystemExit(
                    f"Missing source point-charge "
                    f"file: {source_pc}"
                )

            target_pc = (
                target_dir
                / pc_filename
            )

            shutil.copy2(
                source_pc,
                target_pc,
            )

            if (
                sha256_file(source_pc)
                != sha256_file(target_pc)
            ):
                raise RuntimeError(
                    "Point-charge copy hash mismatch: "
                    f"{target_pc}"
                )

            copied_pc = str(
                target_pc.relative_to(
                    PROJECT_ROOT
                )
            )
            copied_pc_sha256 = (
                sha256_file(target_pc)
            )

        rows.append(
            {
                "calculation_type":
                    calculation_type,
                "frame":
                    int(record.frame),
                "cluster":
                    str(record.cluster),
                "source_job":
                    source_job,
                "target_job":
                    target_job,
                "tracked_root":
                    int(record.tracked_root),
                "NTO_states":
                    "1,2",
                "NTO_threshold":
                    NTO_THRESHOLD,
                "reason":
                    str(record.reason),
                "source_input":
                    str(
                        source_input.relative_to(
                            PROJECT_ROOT
                        )
                    ),
                "target_input":
                    str(
                        target_input.relative_to(
                            PROJECT_ROOT
                        )
                    ),
                "source_input_sha256":
                    sha256_file(source_input),
                "target_input_sha256":
                    sha256_file(target_input),
                "copied_pointcharge_file":
                    copied_pc,
                "copied_pointcharge_sha256":
                    copied_pc_sha256,
            }
        )

        print(
            "Generated and validated: "
            f"{target_input.relative_to(PROJECT_ROOT)}"
        )

    with MANIFEST_CSV.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(rows)

    n_embedded = sum(
        row["calculation_type"] == "embedded"
        for row in rows
    )
    n_vacuum = sum(
        row["calculation_type"]
        == "vacuum_reference"
        for row in rows
    )
    n_pc = sum(
        bool(
            row[
                "copied_pointcharge_file"
            ]
        )
        for row in rows
    )

    with REPORT_MD.open(
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write(
            "# Day019 NTO input audit\n\n"
        )
        handle.write(
            "## Scope\n\n"
        )
        handle.write(
            "Eight representative calculations were "
            "prepared for natural transition orbital "
            "analysis. Both S1 and S2 are requested for "
            "every case to support comparison of the two "
            "low-lying local states.\n\n"
        )
        handle.write(
            "## NTO settings\n\n"
        )
        handle.write(
            "- `DoNTO true`\n"
        )
        handle.write(
            "- `NTOStates 1,2`\n"
        )
        handle.write(
            f"- `NTOThresh {NTO_THRESHOLD}`\n\n"
        )
        handle.write(
            "## Validation\n\n"
        )
        handle.write(
            f"- Inputs generated: "
            f"{len(rows)}/8\n"
        )
        handle.write(
            f"- Embedded cases: "
            f"{n_embedded}\n"
        )
        handle.write(
            f"- Vacuum-reference cases: "
            f"{n_vacuum}\n"
        )
        handle.write(
            f"- Embedded point-charge files "
            f"copied and hash-verified: "
            f"{n_pc}/{n_embedded}\n"
        )
        handle.write(
            "- Unexpected source-input changes: 0\n"
        )
        handle.write(
            "- States requested per case: S1 and S2\n\n"
        )
        handle.write(
            "## Cases\n\n"
        )
        handle.write(
            "| Type | Frame | Site | "
            "Tracked root | Reason | Target job |\n"
        )
        handle.write(
            "|---|---:|---|---:|---|---|\n"
        )

        for row in rows:
            handle.write(
                f"| {row['calculation_type']} "
                f"| {row['frame']} "
                f"| {row['cluster']} "
                f"| S{row['tracked_root']} "
                f"| {row['reason']} "
                f"| `{row['target_job']}` |\n"
            )

    print(
        f"Wrote manifest: "
        f"{MANIFEST_CSV.relative_to(PROJECT_ROOT)}"
    )
    print(
        f"Wrote report: "
        f"{REPORT_MD.relative_to(PROJECT_ROOT)}"
    )
    print(
        "Day019 NTO input generation: PASS"
    )


if __name__ == "__main__":
    main()
