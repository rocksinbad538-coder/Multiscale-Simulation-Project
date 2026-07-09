#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/phase1A/day024_chemical_end_rim_design"

GATE3R = (
    BASE
    / "30_r2_force_field_coverage_preflight"
)

FILE_INVENTORY = (
    GATE3R
    / "r2_local_force_field_file_inventory.csv"
)

TEXT_HITS = (
    GATE3R
    / "r2_local_force_field_bn_h_text_hits.csv"
)

PREFLIGHT_SUMMARY = (
    GATE3R
    / "r2_force_field_coverage_preflight_summary.csv"
)

OUT = (
    BASE
    / "31_r2_local_force_field_asset_audit"
)

ASSET_CLASSIFICATION = (
    OUT
    / "r2_local_force_field_asset_classification.csv"
)

SECTION_INVENTORY = (
    OUT
    / "r2_local_force_field_section_inventory.csv"
)

PARAMETER_LINES = (
    OUT
    / "r2_local_bn_h_parameter_lines.csv"
)

CANDIDATE_MODELS = (
    OUT
    / "r2_local_physical_force_field_candidates.csv"
)

REJECTED_ASSETS = (
    OUT
    / "r2_local_force_field_rejected_assets.csv"
)

SUMMARY = (
    OUT
    / "r2_local_force_field_asset_audit_summary.csv"
)

JSON_OUT = (
    OUT
    / "r2_local_force_field_asset_audit.json"
)

MANIFEST = (
    OUT
    / "r2_local_force_field_asset_audit_manifest.csv"
)

REPORT = (
    OUT
    / "R2_LOCAL_FORCE_FIELD_ASSET_AUDIT_DAY024.md"
)

EXPECTED_PREFLIGHT_DECISION = (
    "R2_FORCE_FIELD_COVERAGE_PREFLIGHT_COMPLETED"
)

FF_PARAMETER_SECTIONS = {
    "atomtypes",
    "bondtypes",
    "angletypes",
    "dihedraltypes",
    "pairtypes",
    "nonbond_params",
    "constrainttypes",
}

TOPOLOGY_INSTANCE_SECTIONS = {
    "moleculetype",
    "atoms",
    "bonds",
    "angles",
    "dihedrals",
    "pairs",
    "constraints",
    "settles",
    "exclusions",
    "system",
    "molecules",
}

DUMMY_MARKERS = (
    "not a physical h-bn force field",
    "fixed dummy scaffold",
    "topology assembly only",
    "dummy scaffold",
    "nonphysical",
    "non-physical",
    "placeholder parameter",
    "placeholder topology",
)

PHYSICAL_MARKERS = (
    "force field",
    "force-field",
    "parameterized",
    "parameterization",
    "lennard-jones",
    "lennard jones",
    "buckingham",
    "morse",
    "charges derived",
    "partial charge",
    "bond force constant",
)

REFERENCE_MARKERS = (
    "doi",
    "journal",
    "reference",
    "citation",
    "published",
    "literature",
)

BN_TYPE_PATTERN = re.compile(
    r"(^|[^A-Za-z0-9])"
    r"(B|N|H|B0|N0|BH|NH|BN|HB|HN)"
    r"([^A-Za-z0-9]|$)",
    re.IGNORECASE,
)

SECTION_PATTERN = re.compile(
    r"^\s*\[\s*([^\]]+?)\s*\]"
)


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty required file: {path}"
        )


def read_rows(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(
            f"No rows found in {path}"
        )

    return rows


def read_one(path: Path) -> dict[str, str]:
    rows = read_rows(path)

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected one row in {path}; found {len(rows)}"
        )

    return rows[0]


def write_rows(
    path: Path,
    rows: list[dict[str, Any]],
) -> None:
    if not rows:
        rows = [
            {
                "status": "NO_ROWS",
            }
        ]

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


def resolve_inventory_path(raw: str) -> Path:
    path = Path(raw)

    if path.is_absolute():
        return path

    return ROOT / path


def normalize_section(value: str) -> str:
    return re.sub(
        r"\s+",
        "_",
        value.strip().lower(),
    )


def strip_comment(line: str) -> str:
    for marker in (
        ";",
        "#",
    ):
        if marker in line:
            line = line.split(
                marker,
                1,
            )[0]

    return line.strip()


def classify_asset(
    path: Path,
    text_lower: str,
    sections: set[str],
    parameter_line_count: int,
    reference_marker_count: int,
) -> tuple[str, str]:
    suffix = path.suffix.lower()

    if any(
        marker in text_lower
        for marker in DUMMY_MARKERS
    ):
        return (
            "DUMMY_OR_TOPOLOGY_ONLY",
            "Explicit dummy/nonphysical marker",
        )

    if path.name.lower() in {
        "sha256sums.txt",
        "checksums.txt",
        "manifest.txt",
    }:
        return (
            "MANIFEST_OR_CHECKSUM",
            "Checksum or manifest file",
        )

    if suffix in {
        ".gro",
        ".pdb",
        ".xyz",
    }:
        return (
            "COORDINATE_FILE",
            "Coordinate file, not a parameter model",
        )

    if (
        sections
        and sections.issubset(
            TOPOLOGY_INSTANCE_SECTIONS
        )
        and not (
            sections
            & FF_PARAMETER_SECTIONS
        )
    ):
        return (
            "INSTANCE_TOPOLOGY_ONLY",
            "Contains instantiated molecular topology but no parameter-type sections",
        )

    if (
        sections
        & FF_PARAMETER_SECTIONS
        and parameter_line_count > 0
    ):
        if reference_marker_count > 0:
            return (
                "PHYSICAL_MODEL_CANDIDATE_WITH_REFERENCE",
                "Contains parameter sections, BN/H terms and reference markers",
            )

        return (
            "PHYSICAL_MODEL_CANDIDATE_UNREFERENCED",
            "Contains parameter sections and BN/H terms but no explicit source marker",
        )

    if parameter_line_count > 0:
        return (
            "BN_H_DATA_WITHOUT_RECOGNIZED_PARAMETER_SECTION",
            "BN/H numerical content found outside recognized parameter-type sections",
        )

    if any(
        marker in text_lower
        for marker in PHYSICAL_MARKERS
    ):
        return (
            "DESCRIPTIVE_OR_DOCUMENTATION_ASSET",
            "Mentions force-field concepts but no explicit BN/H parameter records were identified",
        )

    return (
        "NO_USABLE_BN_H_PARAMETER_CONTENT",
        "No explicit reusable BN/H parameter terms identified",
    )


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        FILE_INVENTORY,
        TEXT_HITS,
        PREFLIGHT_SUMMARY,
    ):
        require_file(required)

    preflight_summary = read_one(
        PREFLIGHT_SUMMARY
    )

    if preflight_summary.get(
        "decision"
    ) != EXPECTED_PREFLIGHT_DECISION:
        raise RuntimeError(
            "Gate 3R preflight is not accepted."
        )

    inventory_rows = read_rows(
        FILE_INVENTORY
    )

    classification_rows = []
    section_rows = []
    parameter_rows = []

    for inventory_row in inventory_rows:
        raw_file = inventory_row.get(
            "file",
            "",
        )

        if not raw_file:
            continue

        path = resolve_inventory_path(
            raw_file
        )

        if not path.is_file():
            classification_rows.append(
                {
                    "file": raw_file,
                    "classification": "MISSING_FILE",
                    "classification_reason": (
                        "Inventory path no longer exists"
                    ),
                    "recognized_sections": "",
                    "parameter_line_count": 0,
                    "reference_marker_count": 0,
                    "dummy_marker_count": 0,
                }
            )

            continue

        try:
            lines = path.read_text(
                encoding="utf-8",
                errors="ignore",
            ).splitlines()
        except OSError:
            classification_rows.append(
                {
                    "file": raw_file,
                    "classification": "UNREADABLE_FILE",
                    "classification_reason": (
                        "Could not read file"
                    ),
                    "recognized_sections": "",
                    "parameter_line_count": 0,
                    "reference_marker_count": 0,
                    "dummy_marker_count": 0,
                }
            )

            continue

        text_lower = "\n".join(
            lines
        ).lower()

        current_section = ""
        sections = set()
        file_parameter_count = 0

        reference_marker_count = sum(
            text_lower.count(marker)
            for marker in REFERENCE_MARKERS
        )

        dummy_marker_count = sum(
            text_lower.count(marker)
            for marker in DUMMY_MARKERS
        )

        for line_number, original_line in enumerate(
            lines,
            start=1,
        ):
            section_match = SECTION_PATTERN.match(
                original_line
            )

            if section_match:
                current_section = normalize_section(
                    section_match.group(1)
                )

                sections.add(
                    current_section
                )

                section_rows.append(
                    {
                        "file": raw_file,
                        "line_number": line_number,
                        "section": current_section,
                    }
                )

                continue

            data_line = strip_comment(
                original_line
            )

            if not data_line:
                continue

            if current_section not in FF_PARAMETER_SECTIONS:
                continue

            if not BN_TYPE_PATTERN.search(
                data_line
            ):
                continue

            tokens = data_line.split()

            numeric_token_count = 0

            for token in tokens:
                try:
                    float(token)
                except ValueError:
                    continue
                else:
                    numeric_token_count += 1

            if numeric_token_count == 0:
                continue

            file_parameter_count += 1

            parameter_rows.append(
                {
                    "file": raw_file,
                    "line_number": line_number,
                    "section": current_section,
                    "text": data_line[:1000],
                    "token_count": len(tokens),
                    "numeric_token_count": numeric_token_count,
                    "has_reference_marker_in_file": (
                        reference_marker_count > 0
                    ),
                    "dummy_marker_in_file": (
                        dummy_marker_count > 0
                    ),
                }
            )

        classification, reason = classify_asset(
            path,
            text_lower,
            sections,
            file_parameter_count,
            reference_marker_count,
        )

        classification_rows.append(
            {
                "file": raw_file,
                "suffix": path.suffix.lower(),
                "size_bytes": path.stat().st_size,
                "classification": classification,
                "classification_reason": reason,
                "recognized_sections": (
                    " | ".join(
                        sorted(
                            sections
                        )
                    )
                ),
                "parameter_sections": (
                    " | ".join(
                        sorted(
                            sections
                            & FF_PARAMETER_SECTIONS
                        )
                    )
                ),
                "instance_sections": (
                    " | ".join(
                        sorted(
                            sections
                            & TOPOLOGY_INSTANCE_SECTIONS
                        )
                    )
                ),
                "parameter_line_count": (
                    file_parameter_count
                ),
                "reference_marker_count": (
                    reference_marker_count
                ),
                "dummy_marker_count": (
                    dummy_marker_count
                ),
                "sha256": (
                    sha256(path)
                    if path.stat().st_size
                    <= 50 * 1024 * 1024
                    else "SKIPPED_GT_50MB"
                ),
            }
        )

    write_rows(
        ASSET_CLASSIFICATION,
        classification_rows,
    )

    write_rows(
        SECTION_INVENTORY,
        section_rows,
    )

    write_rows(
        PARAMETER_LINES,
        parameter_rows,
    )

    candidate_classes = {
        "PHYSICAL_MODEL_CANDIDATE_WITH_REFERENCE",
        "PHYSICAL_MODEL_CANDIDATE_UNREFERENCED",
        "BN_H_DATA_WITHOUT_RECOGNIZED_PARAMETER_SECTION",
    }

    candidate_rows = [
        row
        for row in classification_rows
        if row.get(
            "classification"
        )
        in candidate_classes
    ]

    rejected_rows = [
        row
        for row in classification_rows
        if row.get(
            "classification"
        )
        not in candidate_classes
    ]

    write_rows(
        CANDIDATE_MODELS,
        candidate_rows,
    )

    write_rows(
        REJECTED_ASSETS,
        rejected_rows,
    )

    class_counts = Counter(
        row.get(
            "classification",
            "UNKNOWN",
        )
        for row in classification_rows
    )

    parameter_section_counts = Counter(
        row.get(
            "section",
            "",
        )
        for row in parameter_rows
    )

    referenced_candidates = [
        row
        for row in candidate_rows
        if row.get(
            "classification"
        )
        == "PHYSICAL_MODEL_CANDIDATE_WITH_REFERENCE"
    ]

    unreferenced_candidates = [
        row
        for row in candidate_rows
        if row.get(
            "classification"
        )
        == "PHYSICAL_MODEL_CANDIDATE_UNREFERENCED"
    ]

    dummy_assets = [
        row
        for row in classification_rows
        if row.get(
            "classification"
        )
        == "DUMMY_OR_TOPOLOGY_ONLY"
    ]

    usable_parameter_lines = [
        row
        for row in parameter_rows
        if not bool(
            row.get(
                "dummy_marker_in_file"
            )
        )
    ]

    audit_outcome = (
        "LOCAL_PHYSICAL_BN_H_PARAMETER_CANDIDATES_REQUIRE_SOURCE_VALIDATION"
        if candidate_rows
        else
        "NO_LOCAL_PHYSICAL_BN_H_FORCE_FIELD_CANDIDATE_IDENTIFIED"
    )

    required_next_step = (
        "VALIDATE_LOCAL_CANDIDATE_PROVENANCE_AND_COMPARE_TERM_COVERAGE"
        if candidate_rows
        else
        "PERFORM_PRIMARY_SOURCE_BN_H_FORCE_FIELD_LITERATURE_AUDIT"
    )

    summary = {
        "decision": (
            "R2_LOCAL_FORCE_FIELD_ASSET_AUDIT_COMPLETED"
        ),
        "inventory_files_processed": len(
            classification_rows
        ),
        "dummy_or_topology_only_assets": len(
            dummy_assets
        ),
        "physical_candidates_with_reference_markers": len(
            referenced_candidates
        ),
        "physical_candidates_without_reference_markers": len(
            unreferenced_candidates
        ),
        "all_candidate_assets": len(
            candidate_rows
        ),
        "explicit_BN_H_parameter_lines": len(
            parameter_rows
        ),
        "non_dummy_explicit_BN_H_parameter_lines": len(
            usable_parameter_lines
        ),
        "parameter_section_counts": (
            " | ".join(
                f"{key}:{value}"
                for key, value
                in sorted(
                    parameter_section_counts.items()
                )
            )
        ),
        "audit_outcome": (
            audit_outcome
        ),
        "force_field_coverage_established": False,
        "topology_generation_authorized": False,
        "formal_charge_assignment_authorized": False,
        "force_field_parameterization_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "QM_calculation_authorized": False,
        "required_next_step": (
            required_next_step
        ),
    }

    write_rows(
        SUMMARY,
        [
            summary
        ],
    )

    JSON_OUT.write_text(
        json.dumps(
            {
                "summary": summary,
                "classification_counts": dict(
                    class_counts
                ),
                "parameter_section_counts": dict(
                    parameter_section_counts
                ),
                "limitations": [
                    (
                        "A file is only a candidate if explicit parameter "
                        "records are present; this still does not establish "
                        "scientific validity or transferability."
                    ),
                    (
                        "Reference-marker detection is textual and must be "
                        "verified against the original publication."
                    ),
                    (
                        "No topology, charges or parameters are assigned."
                    ),
                    (
                        "No minimization, MD or QM calculation is performed."
                    ),
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    write_rows(
        MANIFEST,
        [
            {
                "role": "Gate3R_file_inventory",
                "file": str(
                    FILE_INVENTORY.relative_to(
                        ROOT
                    )
                ),
                "sha256": sha256(
                    FILE_INVENTORY
                ),
            },
            {
                "role": "Gate3R_text_hits",
                "file": str(
                    TEXT_HITS.relative_to(
                        ROOT
                    )
                ),
                "sha256": sha256(
                    TEXT_HITS
                ),
            },
            {
                "role": "Gate3R_summary",
                "file": str(
                    PREFLIGHT_SUMMARY.relative_to(
                        ROOT
                    )
                ),
                "sha256": sha256(
                    PREFLIGHT_SUMMARY
                ),
            },
        ],
    )

    class_lines = "\n".join(
        (
            f"- {name}: **{count}**"
        )
        for name, count in sorted(
            class_counts.items()
        )
    )

    candidate_lines = "\n".join(
        (
            f"- `{row['file']}` — "
            f"{row['classification']}; "
            f"parameter lines="
            f"{row['parameter_line_count']}; "
            f"sections="
            f"{row['parameter_sections'] or 'NONE'}"
        )
        for row in candidate_rows
    )

    if not candidate_lines:
        candidate_lines = (
            "- No local physical BN/H parameter candidate was identified."
        )

    REPORT.write_text(
        f"""# R2 Local Force-Field Asset Audit

## Files processed

- Inventory files processed:
  **{len(classification_rows)}**

## Classification

{class_lines}

## Explicit parameter records

- All explicit BN/H parameter lines:
  **{len(parameter_rows)}**
- Non-dummy explicit BN/H parameter lines:
  **{len(usable_parameter_lines)}**
- Parameter sections:
  **{summary['parameter_section_counts'] or 'NONE'}**

## Candidate assets

{candidate_lines}

## Audit outcome

- Outcome:
  **{audit_outcome}**
- Force-field coverage established:
  **NO**
- Topology generation authorized:
  **NO**
- Charge assignment authorized:
  **NO**
- Parameterization authorized:
  **NO**
- Energy minimization authorized:
  **NO**
- MD authorized:
  **NO**
- QM calculation authorized:
  **NO**

## Required next step

`{required_next_step}`
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 local force-field asset audit completed."
    )

    print(
        "Files processed / dummy-topology-only assets: "
        f"{len(classification_rows)}/"
        f"{len(dummy_assets)}"
    )

    print(
        "Physical candidates referenced/unreferenced/total: "
        f"{len(referenced_candidates)}/"
        f"{len(unreferenced_candidates)}/"
        f"{len(candidate_rows)}"
    )

    print(
        "Explicit BN-H parameter lines all/non-dummy: "
        f"{len(parameter_rows)}/"
        f"{len(usable_parameter_lines)}"
    )

    print(
        "Parameter section counts: "
        f"{summary['parameter_section_counts'] or 'NONE'}"
    )

    print(
        "Classification counts: "
        + " | ".join(
            f"{key}:{value}"
            for key, value in sorted(
                class_counts.items()
            )
        )
    )

    print(
        f"Audit outcome: {audit_outcome}"
    )

    print(
        "Force-field coverage established: NO"
    )

    print(
        "Topology generation authorized: NO"
    )

    print(
        "Formal charge assignment authorized: NO"
    )

    print(
        "Force-field parameterization authorized: NO"
    )

    print(
        "Energy minimization authorized: NO"
    )

    print(
        "MD authorized: NO"
    )

    print(
        "QM calculation authorized: NO"
    )

    print(
        f"Required next step: {required_next_step}"
    )

    for path in (
        ASSET_CLASSIFICATION,
        SECTION_INVENTORY,
        PARAMETER_LINES,
        CANDIDATE_MODELS,
        REJECTED_ASSETS,
        SUMMARY,
        JSON_OUT,
        MANIFEST,
        REPORT,
    ):
        print(
            f"Wrote: {path.relative_to(ROOT)}"
        )


if __name__ == "__main__":
    main()
