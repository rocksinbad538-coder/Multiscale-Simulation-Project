#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import mimetypes
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/phase1A/day024_chemical_end_rim_design"

GATE3R3A = (
    BASE
    / "34_r2_primary_source_provenance_correction"
)

PRIORITY_LEDGER = (
    GATE3R3A
    / "r2_prioritized_source_artifact_retrieval.csv"
)

PROVENANCE_SUMMARY = (
    GATE3R3A
    / "r2_primary_source_provenance_correction_summary.csv"
)

OUT = (
    BASE
    / "35_r2_p1_source_artifact_retrieval"
)

RAW = OUT / "raw"
EXTRACTED = OUT / "extracted"

RETRIEVAL_LEDGER = (
    OUT
    / "r2_p1_artifact_retrieval_ledger.csv"
)

FILE_INVENTORY = (
    OUT
    / "r2_p1_retrieved_file_inventory.csv"
)

ARCHIVE_INVENTORY = (
    OUT
    / "r2_p1_archive_member_inventory.csv"
)

SUMMARY = (
    OUT
    / "r2_p1_artifact_retrieval_summary.csv"
)

JSON_OUT = (
    OUT
    / "r2_p1_artifact_retrieval.json"
)

REPORT = (
    OUT
    / "R2_P1_SOURCE_ARTIFACT_RETRIEVAL_DAY024.md"
)

EXPECTED_PROVENANCE_DECISION = (
    "R2_PRIMARY_SOURCE_PROVENANCE_CORRECTED_"
    "ARTIFACT_RETRIEVAL_PRIORITIZED"
)


ARTIFACTS = [
    {
        "artifact_id": "GHORAI_ARTICLE_LANDING",
        "model_name": "GHORAI_FUNCTIONALIZED_HBN",
        "DOI": "10.1063/5.0242541",
        "artifact_type": "FULL_ARTICLE_LANDING_PAGE",
        "url": "https://doi.org/10.1063/5.0242541",
        "filename": "ghorai_2025_article_landing.html",
        "expected_kind": "HTML",
        "priority": 12,
    },
    {
        "artifact_id": "LELE_ARTICLE_OPEN_COPY",
        "model_name": "LELE_GAS_PHASE_HBN_REAXFF",
        "DOI": "10.1021/acs.jpca.1c09648",
        "artifact_type": "FULL_ARTICLE_OPEN_COPY",
        "url": "https://par.nsf.gov/servlets/purl/10321378",
        "filename": "lele_2022_reaxff_article.pdf",
        "expected_kind": "PDF",
        "priority": 22,
    },
    {
        "artifact_id": "LELE_SI_PDF",
        "model_name": "LELE_GAS_PHASE_HBN_REAXFF",
        "DOI": "10.1021/acs.jpca.1c09648",
        "artifact_type": "SUPPORTING_INFORMATION",
        "url": "https://doi.org/10.1021/acs.jpca.1c09648.s001",
        "filename": "jp1c09648_si_001_download",
        "expected_kind": "PDF_OR_HTML",
        "priority": 21,
    },
    {
        "artifact_id": "LELE_REAXFF_PARAMETERS",
        "model_name": "LELE_GAS_PHASE_HBN_REAXFF",
        "DOI": "10.1021/acs.jpca.1c09648",
        "artifact_type": "PARAMETER_FILE",
        "url": "https://doi.org/10.1021/acs.jpca.1c09648.s002",
        "filename": "jp1c09648_si_002_download",
        "expected_kind": "TEXT_OR_HTML",
        "priority": 20,
    },
    {
        "artifact_id": "IFFR_BNNT_DATASET",
        "model_name": "IFF_R_BNNT",
        "DOI": "10.1021/acsanm.2c05285",
        "artifact_type": "PARAMETER_FILE_OR_REPOSITORY",
        "url": (
            "https://figshare.com/ndownloader/articles/22144074/"
            "versions/1"
        ),
        "filename": "bamane_2023_bnnt_iff_r_dataset.zip",
        "expected_kind": "ZIP",
        "priority": 30,
    },
    {
        "artifact_id": "IFFR_BNNT_ARTICLE_COPY",
        "model_name": "IFF_R_BNNT",
        "DOI": "10.1021/acsanm.2c05285",
        "artifact_type": "FULL_ARTICLE_OPEN_COPY",
        "url": (
            "https://digitalcommons.mtu.edu/cgi/viewcontent.cgi"
            "?article=36286&context=michigantech-p"
        ),
        "filename": "bamane_2023_bnnt_article.pdf",
        "expected_kind": "PDF",
        "priority": 32,
    },
    {
        "artifact_id": "RAJAN_SI",
        "model_name": "RAJAN_HBN_CLASSICAL_FF",
        "DOI": "10.1021/acs.jpclett.7b03443",
        "artifact_type": "SUPPORTING_INFORMATION",
        "url": "https://doi.org/10.1021/acs.jpclett.7b03443.s001",
        "filename": "jz7b03443_si_001_download",
        "expected_kind": "PDF_OR_HTML",
        "priority": 41,
    },
    {
        "artifact_id": "RAJAN_ARTICLE_METADATA",
        "model_name": "RAJAN_HBN_CLASSICAL_FF",
        "DOI": "10.1021/acs.jpclett.7b03443",
        "artifact_type": "ARTICLE_METADATA",
        "url": "https://doi.org/10.1021/acs.jpclett.7b03443",
        "filename": "rajan_2018_article_landing.html",
        "expected_kind": "HTML",
        "priority": 42,
    },
]


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
        rows = [{"status": "NO_ROWS"}]

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


def detect_kind(path: Path) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        return "EMPTY_OR_MISSING"

    with path.open("rb") as handle:
        head = handle.read(4096)

    if head.startswith(b"%PDF-"):
        return "PDF"

    if head.startswith(b"PK\x03\x04"):
        return "ZIP"

    lower = head.lstrip().lower()

    if (
        lower.startswith(b"<!doctype html")
        or lower.startswith(b"<html")
        or b"<html" in lower[:1000]
    ):
        return "HTML"

    if b"\x00" not in head:
        return "TEXT"

    guessed, _ = mimetypes.guess_type(path.name)

    if guessed:
        return guessed.upper()

    return "BINARY"


def download(
    url: str,
    destination: Path,
) -> tuple[int, str]:
    command = [
        "curl",
        "--location",
        "--fail-with-body",
        "--silent",
        "--show-error",
        "--retry",
        "2",
        "--connect-timeout",
        "20",
        "--max-time",
        "180",
        "--user-agent",
        (
            "Mozilla/5.0 "
            "(compatible; R2-primary-source-audit/1.0)"
        ),
        "--output",
        str(destination),
        "--write-out",
        "%{http_code}",
        url,
    ]

    process = subprocess.run(
        command,
        text=True,
        capture_output=True,
        check=False,
    )

    http_code = process.stdout.strip()

    return (
        process.returncode,
        http_code,
    )


def safe_extract_zip(
    archive: Path,
    destination: Path,
) -> list[dict[str, Any]]:
    destination.mkdir(
        parents=True,
        exist_ok=True,
    )

    member_rows = []

    with zipfile.ZipFile(
        archive,
        "r",
    ) as handle:
        destination_root = destination.resolve()

        for info in handle.infolist():
            target = (
                destination
                / info.filename
            ).resolve()

            if (
                target != destination_root
                and destination_root not in target.parents
            ):
                raise RuntimeError(
                    f"Unsafe ZIP member path: {info.filename}"
                )

            member_rows.append(
                {
                    "archive": str(
                        archive.relative_to(ROOT)
                    ),
                    "member": info.filename,
                    "uncompressed_size_bytes": (
                        info.file_size
                    ),
                    "compressed_size_bytes": (
                        info.compress_size
                    ),
                    "is_directory": (
                        info.is_dir()
                    ),
                }
            )

        handle.extractall(
            destination
        )

    return member_rows


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    RAW.mkdir(
        parents=True,
        exist_ok=True,
    )

    EXTRACTED.mkdir(
        parents=True,
        exist_ok=True,
    )

    require_file(
        PRIORITY_LEDGER
    )

    provenance = read_one(
        PROVENANCE_SUMMARY
    )

    if provenance.get(
        "decision"
    ) != EXPECTED_PROVENANCE_DECISION:
        raise RuntimeError(
            "Gate 3R.3a provenance correction is not accepted."
        )

    priority_rows = read_rows(
        PRIORITY_LEDGER
    )

    p1_dois = {
        row["DOI"]
        for row in priority_rows
        if row["priority_group"] == "P1_IMMEDIATE"
    }

    expected_p1_dois = {
        "10.1063/5.0242541",
        "10.1021/acs.jpca.1c09648",
        "10.1021/acsanm.2c05285",
        "10.1021/acs.jpclett.7b03443",
    }

    if p1_dois != expected_p1_dois:
        raise RuntimeError(
            "Priority-1 DOI set does not match the approved ledger."
        )

    retrieval_rows = []
    archive_rows = []

    for artifact in sorted(
        ARTIFACTS,
        key=lambda row: int(row["priority"]),
    ):
        destination = (
            RAW
            / artifact["filename"]
        )

        if destination.exists():
            destination.unlink()

        return_code, http_code = download(
            artifact["url"],
            destination,
        )

        actual_kind = detect_kind(
            destination
        )

        size_bytes = (
            destination.stat().st_size
            if destination.exists()
            else 0
        )

        status = (
            "RETRIEVED"
            if (
                return_code == 0
                and size_bytes > 0
                and actual_kind
                not in {
                    "EMPTY_OR_MISSING",
                }
            )
            else "FAILED"
        )

        expected = artifact[
            "expected_kind"
        ]

        if (
            status == "RETRIEVED"
            and expected == "PDF"
            and actual_kind != "PDF"
        ):
            status = "RETRIEVED_WRONG_FORMAT"

        elif (
            status == "RETRIEVED"
            and expected == "ZIP"
            and actual_kind != "ZIP"
        ):
            status = "RETRIEVED_WRONG_FORMAT"

        elif (
            status == "RETRIEVED"
            and expected == "HTML"
            and actual_kind != "HTML"
        ):
            status = "RETRIEVED_UNEXPECTED_FORMAT"

        elif (
            status == "RETRIEVED"
            and expected == "PDF_OR_HTML"
            and actual_kind not in {
                "PDF",
                "HTML",
            }
        ):
            status = "RETRIEVED_UNEXPECTED_FORMAT"

        elif (
            status == "RETRIEVED"
            and expected == "TEXT_OR_HTML"
            and actual_kind not in {
                "TEXT",
                "HTML",
            }
        ):
            status = "RETRIEVED_UNEXPECTED_FORMAT"

        if (
            actual_kind == "HTML"
            and expected
            in {
                "PDF",
                "ZIP",
                "TEXT_OR_HTML",
                "PDF_OR_HTML",
            }
        ):
            status = (
                "LANDING_OR_ACCESS_PAGE_RETRIEVED"
            )

        extracted_directory = ""

        if (
            actual_kind == "ZIP"
            and status == "RETRIEVED"
        ):
            extraction_target = (
                EXTRACTED
                / artifact["artifact_id"]
            )

            if extraction_target.exists():
                shutil.rmtree(
                    extraction_target
                )

            members = safe_extract_zip(
                destination,
                extraction_target,
            )

            archive_rows.extend(
                members
            )

            extracted_directory = str(
                extraction_target.relative_to(
                    ROOT
                )
            )

        retrieval_rows.append(
            {
                **artifact,
                "curl_return_code": return_code,
                "http_code": http_code,
                "status": status,
                "actual_kind": actual_kind,
                "size_bytes": size_bytes,
                "local_file": (
                    str(
                        destination.relative_to(
                            ROOT
                        )
                    )
                    if destination.exists()
                    else ""
                ),
                "sha256": (
                    sha256(destination)
                    if destination.exists()
                    and size_bytes > 0
                    else ""
                ),
                "extracted_directory": (
                    extracted_directory
                ),
                "parameter_adoption_authorized": False,
            }
        )

    write_rows(
        RETRIEVAL_LEDGER,
        retrieval_rows,
    )

    write_rows(
        ARCHIVE_INVENTORY,
        archive_rows,
    )

    file_rows = []

    for path in sorted(
        OUT.rglob("*")
    ):
        if not path.is_file():
            continue

        if path in {
            RETRIEVAL_LEDGER,
            FILE_INVENTORY,
            ARCHIVE_INVENTORY,
            SUMMARY,
            JSON_OUT,
            REPORT,
        }:
            continue

        file_rows.append(
            {
                "file": str(
                    path.relative_to(
                        ROOT
                    )
                ),
                "kind": detect_kind(
                    path
                ),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(
                    path
                ),
            }
        )

    write_rows(
        FILE_INVENTORY,
        file_rows,
    )

    successful_binary_artifacts = [
        row
        for row in retrieval_rows
        if (
            row["status"] == "RETRIEVED"
            and row["actual_kind"]
            in {
                "PDF",
                "ZIP",
                "TEXT",
            }
        )
    ]

    access_or_landing_rows = [
        row
        for row in retrieval_rows
        if row["status"]
        in {
            "LANDING_OR_ACCESS_PAGE_RETRIEVED",
            "RETRIEVED_WRONG_FORMAT",
            "RETRIEVED_UNEXPECTED_FORMAT",
        }
    ]

    failed_rows = [
        row
        for row in retrieval_rows
        if row["status"] == "FAILED"
    ]

    retrieved_models = sorted(
        {
            row["model_name"]
            for row in successful_binary_artifacts
        }
    )

    gates = {
        "approved_priority1_DOI_set_was_used": (
            p1_dois == expected_p1_dois
        ),
        "all_retrieved_files_have_SHA256": all(
            (
                row["sha256"] != ""
                if row["size_bytes"] > 0
                else True
            )
            for row in retrieval_rows
        ),
        "ZIP_files_are_extracted_only_after_magic_validation": True,
        "no_artifact_authorizes_parameter_adoption": all(
            row[
                "parameter_adoption_authorized"
            ]
            is False
            for row in retrieval_rows
        ),
        "no_topology_charges_parameters_MD_or_QM_generated": True,
    }

    failed_gates = [
        name
        for name, passed in gates.items()
        if not passed
    ]

    accepted = (
        len(
            failed_gates
        )
        == 0
    )

    decision = (
        "R2_P1_SOURCE_ARTIFACT_RETRIEVAL_AUDITED"
        if accepted
        else
        "R2_P1_SOURCE_ARTIFACT_RETRIEVAL_REQUIRES_REVIEW"
    )

    required_next_step = (
        "AUDIT_RETRIEVED_P1_ARTIFACT_CONTENT_AND_PARAMETER_PROVENANCE"
        if successful_binary_artifacts
        else
        "RESOLVE_P1_ARTIFACT_ACCESS_AND_DOWNLOAD_BLOCKS"
    )

    summary = {
        "decision": decision,
        "artifacts_attempted": len(
            retrieval_rows
        ),
        "binary_or_text_artifacts_retrieved": len(
            successful_binary_artifacts
        ),
        "landing_or_access_pages_retrieved": len(
            access_or_landing_rows
        ),
        "failed_retrievals": len(
            failed_rows
        ),
        "archives_extracted": sum(
            row["actual_kind"] == "ZIP"
            and row["status"] == "RETRIEVED"
            for row in retrieval_rows
        ),
        "archive_members": len(
            archive_rows
        ),
        "models_with_at_least_one_binary_artifact": (
            " | ".join(
                retrieved_models
            )
        ),
        "force_field_coverage_established": False,
        "parameter_adoption_authorized": False,
        "topology_generation_authorized": False,
        "formal_charge_assignment_authorized": False,
        "force_field_parameterization_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "QM_calculation_authorized": False,
        "failed_gates": (
            " | ".join(
                failed_gates
            )
        ),
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
                "retrievals": retrieval_rows,
                "gates": gates,
                "limitations": [
                    (
                        "A successful download does not establish "
                        "scientific validity or R2 transferability."
                    ),
                    (
                        "HTML landing or access-control pages are "
                        "retained as evidence but are not treated as "
                        "articles, SI or parameter files."
                    ),
                    (
                        "No downloaded parameter set is copied into "
                        "an accepted-parameters directory."
                    ),
                    (
                        "No topology, charge assignment, "
                        "parameterization, minimization, MD or QM "
                        "calculation is performed."
                    ),
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    status_lines = "\n".join(
        (
            f"- `{row['artifact_id']}`: "
            f"{row['status']}; kind={row['actual_kind']}; "
            f"bytes={row['size_bytes']}; HTTP={row['http_code']}"
        )
        for row in retrieval_rows
    )

    REPORT.write_text(
        f"""# R2 Priority-1 Source Artifact Retrieval

## Retrieval results

{status_lines}

## Summary

- Artifacts attempted:
  **{len(retrieval_rows)}**
- Binary/text artifacts retrieved:
  **{len(successful_binary_artifacts)}**
- Landing/access pages:
  **{len(access_or_landing_rows)}**
- Failed retrievals:
  **{len(failed_rows)}**
- Archives extracted:
  **{summary['archives_extracted']}**
- Archive members:
  **{len(archive_rows)}**

## Restrictions

- Force-field coverage established: **NO**
- Parameter adoption authorized: **NO**
- Topology generation authorized: **NO**
- Charge assignment authorized: **NO**
- Parameterization authorized: **NO**
- Energy minimization authorized: **NO**
- MD authorized: **NO**
- QM calculation authorized: **NO**

## Required next step

`{required_next_step}`
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 priority-1 source artifact "
        "retrieval completed."
    )

    print(
        "Artifacts attempted / binary-text retrieved / "
        "landing-access / failed: "
        f"{len(retrieval_rows)}/"
        f"{len(successful_binary_artifacts)}/"
        f"{len(access_or_landing_rows)}/"
        f"{len(failed_rows)}"
    )

    print(
        "Archives extracted / archive members: "
        f"{summary['archives_extracted']}/"
        f"{len(archive_rows)}"
    )

    print(
        "Models with binary artifacts: "
        + (
            summary[
                "models_with_at_least_one_binary_artifact"
            ]
            or "NONE"
        )
    )

    for row in retrieval_rows:
        print(
            f"{row['artifact_id']} "
            f"status/kind/http/bytes: "
            f"{row['status']}/"
            f"{row['actual_kind']}/"
            f"{row['http_code']}/"
            f"{row['size_bytes']}"
        )

    print(
        f"Decision: {decision}"
    )

    print(
        "Failed gates: "
        + (
            "NONE"
            if not failed_gates
            else " | ".join(
                failed_gates
            )
        )
    )

    print(
        "Force-field coverage established: NO"
    )

    print(
        "Parameter adoption authorized: NO"
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
        RETRIEVAL_LEDGER,
        FILE_INVENTORY,
        ARCHIVE_INVENTORY,
        SUMMARY,
        JSON_OUT,
        REPORT,
    ):
        print(
            f"Wrote: {path.relative_to(ROOT)}"
        )

    if not accepted:
        raise RuntimeError(
            "Priority-1 artifact retrieval audit requires review."
        )


if __name__ == "__main__":
    main()
