#!/usr/bin/env python3
"""
DAY040 / D040-A7a

Official supplementary-archive recognition and safe extraction for
Ghorai et al. 2025.

This block repairs the artifact model used by A7:

- The canonical article is a PDF.
- The official supplementary artifact is ghorai2025_SI.zip.
- The ZIP contains hBN.tar.
- A separate supplementary PDF is not required.
- The ZIP and nested TAR are extracted with path-traversal protection.
- Every extracted artifact is inventoried with SHA256.

No force-field parameter is adopted.
No accepted topology or coordinate file is modified.
No MD calculation is executed.
"""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import tarfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

A7_REPORT = (
    ROOT
    / "runs/phase1A/day040_phase1A_G_primary_literature_extraction"
    / "QM_F06_UPPER_V7A_R1_PRIMARY_LITERATURE_EXTRACTION.json"
)

CANONICAL_ARTICLE = (
    ROOT
    / "references/phase1A/primary_sources/ghorai_2025"
    / "044705_1_5.0242541.pdf"
)

SUPPLEMENTARY_ROOT = (
    ROOT
    / "references/phase1A/primary_sources/ghorai_2025"
    / "supplementary"
)

OFFICIAL_SI_ZIP = (
    SUPPLEMENTARY_ROOT
    / "ghorai2025_SI.zip"
)

EXTRACTION_ROOT = (
    SUPPLEMENTARY_ROOT
    / "official_extracted"
)

ZIP_EXTRACTION_DIR = (
    EXTRACTION_ROOT
    / "zip"
)

TAR_EXTRACTION_DIR = (
    EXTRACTION_ROOT
    / "hBN_tar"
)

OUTPUT_DIR = (
    ROOT
    / "runs/phase1A/day040_phase1A_G_supplementary_archive_recognition"
)

ARCHIVE_INVENTORY_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_GHORAI2025_SUPPLEMENTARY_ARCHIVE_INVENTORY.csv"
)

EXTRACTED_FILE_INVENTORY_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_GHORAI2025_EXTRACTED_FILE_INVENTORY.csv"
)

REPORT_JSON = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_GHORAI2025_SUPPLEMENTARY_ARCHIVE_RECOGNITION.json"
)

EXPECTED_ARTICLE_SHA256 = (
    "01bd5938fc1e66fb77e4bbd9aba599f51be22e9c4b59853b9e3bb6c9268fff8c"
)

EXPECTED_SI_ZIP_SHA256 = (
    "95cbef4d91ed6e1ac235c7dd265c6b94bca0ee0ce01313117b0cf8e82708b3dc"
)

EXPECTED_NESTED_MEMBER = "hBN.tar"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def relative(path: Path) -> str:
    return str(path.relative_to(ROOT))


def ensure_within_directory(
    destination: Path,
    candidate: Path,
) -> None:
    destination_resolved = destination.resolve()
    candidate_resolved = candidate.resolve()

    if (
        candidate_resolved != destination_resolved
        and destination_resolved
        not in candidate_resolved.parents
    ):
        raise RuntimeError(
            "Unsafe archive member path detected: "
            f"{candidate}"
        )


def safe_extract_zip(
    archive_path: Path,
    destination: Path,
) -> list[Path]:
    destination.mkdir(
        parents=True,
        exist_ok=True,
    )

    extracted = []

    with zipfile.ZipFile(
        archive_path,
        "r",
    ) as archive:
        for info in archive.infolist():
            target = destination / info.filename

            ensure_within_directory(
                destination,
                target,
            )

            if info.is_dir():
                target.mkdir(
                    parents=True,
                    exist_ok=True,
                )
                continue

            target.parent.mkdir(
                parents=True,
                exist_ok=True,
            )

            with archive.open(
                info,
                "r",
            ) as source_handle:
                with target.open(
                    "wb",
                ) as target_handle:
                    shutil.copyfileobj(
                        source_handle,
                        target_handle,
                    )

            extracted.append(target)

    return extracted


def safe_extract_tar(
    archive_path: Path,
    destination: Path,
) -> list[Path]:
    destination.mkdir(
        parents=True,
        exist_ok=True,
    )

    extracted = []

    with tarfile.open(
        archive_path,
        "r:*",
    ) as archive:
        members = archive.getmembers()

        for member in members:
            target = destination / member.name

            ensure_within_directory(
                destination,
                target,
            )

            if member.issym() or member.islnk():
                raise RuntimeError(
                    "Symbolic and hard links are not accepted "
                    f"in the supplementary archive: {member.name}"
                )

        archive.extractall(
            destination,
            members=members,
            filter="data",
        )

        for member in members:
            target = destination / member.name

            if target.is_file():
                extracted.append(target)

    return extracted


def file_record(
    path: Path,
    artifact_level: str,
    parent_archive: str,
) -> dict:
    return {
        "artifact_level": artifact_level,
        "parent_archive": parent_archive,
        "path": relative(path),
        "suffix": path.suffix.lower(),
        "bytes": int(path.stat().st_size),
        "sha256": sha256(path),
    }


print("=" * 100)
print("DAY040 / D040-A7a — SUPPLEMENTARY ARCHIVE RECOGNITION AND EXTRACTION")
print("=" * 100)


print("\n[1] SOURCE ARTIFACTS")

for path in (
    A7_REPORT,
    CANONICAL_ARTICLE,
    OFFICIAL_SI_ZIP,
):
    if not path.is_file():
        raise FileNotFoundError(path)

    print(
        f"FOUND bytes={path.stat().st_size:8d} "
        f"sha256={sha256(path)} "
        f"{path}"
    )

article_sha_gate = (
    sha256(CANONICAL_ARTICLE)
    == EXPECTED_ARTICLE_SHA256
)

si_zip_sha_gate = (
    sha256(OFFICIAL_SI_ZIP)
    == EXPECTED_SI_ZIP_SHA256
)

print(
    f"canonical_article_sha_gate="
    f"{'PASS' if article_sha_gate else 'FAIL'}"
)
print(
    f"official_SI_zip_sha_gate="
    f"{'PASS' if si_zip_sha_gate else 'FAIL'}"
)


print("\n[2] PRESERVE A7 DIAGNOSTIC")

a7_report = json.loads(
    A7_REPORT.read_text(
        encoding="utf-8"
    )
)

a7_decision = a7_report.get(
    "decision"
)

print(
    f"A7_recorded_decision={a7_decision}"
)
print(
    "A7_interpretation="
    "ARTIFACT_RECOGNITION_LOGIC_REVIEW_REQUIRED"
)


print("\n[3] ZIP CONTENT AUDIT")

with zipfile.ZipFile(
    OFFICIAL_SI_ZIP,
    "r",
) as archive:
    zip_members = archive.infolist()

print(
    f"zip_member_count={len(zip_members)}"
)

for info in zip_members:
    print(
        f"member={info.filename} "
        f"bytes={info.file_size} "
        f"compressed_bytes={info.compress_size}"
    )

nested_tar_members = [
    info
    for info in zip_members
    if (
        not info.is_dir()
        and Path(
            info.filename
        ).name
        == EXPECTED_NESTED_MEMBER
    )
]

nested_tar_unique_gate = (
    len(nested_tar_members)
    == 1
)

print(
    f"unique_hBN_tar_member_gate="
    f"{'PASS' if nested_tar_unique_gate else 'FAIL'}"
)

if not nested_tar_unique_gate:
    raise RuntimeError(
        "The official SI ZIP must contain exactly one hBN.tar"
    )


print("\n[4] CLEAN ISOLATED EXTRACTION")

if EXTRACTION_ROOT.exists():
    shutil.rmtree(
        EXTRACTION_ROOT
    )

ZIP_EXTRACTION_DIR.mkdir(
    parents=True,
    exist_ok=False,
)

zip_extracted_files = safe_extract_zip(
    OFFICIAL_SI_ZIP,
    ZIP_EXTRACTION_DIR,
)

print(
    f"zip_extracted_file_count="
    f"{len(zip_extracted_files)}"
)

for path in zip_extracted_files:
    print(
        f"ZIP_EXTRACTED bytes={path.stat().st_size} "
        f"sha256={sha256(path)} "
        f"path={path}"
    )

nested_tar_paths = [
    path
    for path in zip_extracted_files
    if path.name == EXPECTED_NESTED_MEMBER
]

if len(nested_tar_paths) != 1:
    raise RuntimeError(
        "The extracted ZIP does not contain exactly one hBN.tar"
    )

nested_tar_path = nested_tar_paths[0]


print("\n[5] NESTED TAR CONTENT AUDIT")

with tarfile.open(
    nested_tar_path,
    "r:*",
) as archive:
    tar_members = archive.getmembers()

print(
    f"tar_member_count={len(tar_members)}"
)

for member in tar_members:
    print(
        f"member={member.name} "
        f"type={'file' if member.isfile() else 'directory_or_other'} "
        f"bytes={member.size}"
    )


print("\n[6] SAFE NESTED TAR EXTRACTION")

tar_extracted_files = safe_extract_tar(
    nested_tar_path,
    TAR_EXTRACTION_DIR,
)

print(
    f"tar_extracted_file_count="
    f"{len(tar_extracted_files)}"
)

for path in tar_extracted_files:
    print(
        f"TAR_EXTRACTED bytes={path.stat().st_size} "
        f"sha256={sha256(path)} "
        f"path={path}"
    )


print("\n[7] CONTENT CLASSIFICATION")

recognized_suffixes = {
    ".in",
    ".lmp",
    ".data",
    ".txt",
    ".dat",
    ".log",
    ".restart",
    ".lammpstrj",
}

parameter_relevant_files = [
    path
    for path in tar_extracted_files
    if (
        path.suffix.lower()
        in recognized_suffixes
        or any(
            token in path.name.lower()
            for token in (
                "param",
                "force",
                "ffield",
                "input",
                "lammps",
                "hbn",
                "borazine",
            )
        )
    )
]

print(
    f"parameter_relevant_file_count="
    f"{len(parameter_relevant_files)}"
)

for path in parameter_relevant_files:
    print(
        f"PARAMETER_RELEVANT_FILE "
        f"path={relative(path)}"
    )


print("\n[8] WRITE INVENTORIES")

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

archive_rows = [
    {
        "artifact_type": "PRIMARY_ARTICLE",
        "path": relative(
            CANONICAL_ARTICLE
        ),
        "bytes": int(
            CANONICAL_ARTICLE.stat().st_size
        ),
        "sha256": sha256(
            CANONICAL_ARTICLE
        ),
        "recognition_status": (
            "CANONICAL_PRIMARY_ARTICLE"
        ),
    },
    {
        "artifact_type": (
            "OFFICIAL_SUPPLEMENTARY_ZIP"
        ),
        "path": relative(
            OFFICIAL_SI_ZIP
        ),
        "bytes": int(
            OFFICIAL_SI_ZIP.stat().st_size
        ),
        "sha256": sha256(
            OFFICIAL_SI_ZIP
        ),
        "recognition_status": (
            "OFFICIAL_SUPPLEMENTARY_ARTIFACT"
        ),
    },
    {
        "artifact_type": (
            "NESTED_PARAMETER_TAR"
        ),
        "path": relative(
            nested_tar_path
        ),
        "bytes": int(
            nested_tar_path.stat().st_size
        ),
        "sha256": sha256(
            nested_tar_path
        ),
        "recognition_status": (
            "EXTRACTED_FROM_OFFICIAL_SUPPLEMENTARY_ZIP"
        ),
    },
]

with ARCHIVE_INVENTORY_CSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(
            archive_rows[0].keys()
        ),
    )
    writer.writeheader()
    writer.writerows(
        archive_rows
    )

extracted_rows = []

for path in zip_extracted_files:
    extracted_rows.append(
        file_record(
            path,
            "ZIP_EXTRACTED",
            relative(
                OFFICIAL_SI_ZIP
            ),
        )
    )

for path in tar_extracted_files:
    extracted_rows.append(
        file_record(
            path,
            "TAR_EXTRACTED",
            relative(
                nested_tar_path
            ),
        )
    )

extracted_fieldnames = [
    "artifact_level",
    "parent_archive",
    "path",
    "suffix",
    "bytes",
    "sha256",
]

with EXTRACTED_FILE_INVENTORY_CSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=extracted_fieldnames,
    )
    writer.writeheader()
    writer.writerows(
        extracted_rows
    )


print("\n[9] GATES")

gates = {
    "A7_diagnostic_preserved_gate": (
        A7_REPORT.is_file()
    ),
    "canonical_article_present_gate": (
        CANONICAL_ARTICLE.is_file()
    ),
    "canonical_article_sha_gate": (
        article_sha_gate
    ),
    "official_SI_zip_present_gate": (
        OFFICIAL_SI_ZIP.is_file()
    ),
    "official_SI_zip_sha_gate": (
        si_zip_sha_gate
    ),
    "unique_hBN_tar_member_gate": (
        nested_tar_unique_gate
    ),
    "zip_safe_extraction_gate": (
        len(
            zip_extracted_files
        )
        > 0
    ),
    "nested_tar_safe_extraction_gate": (
        len(
            tar_extracted_files
        )
        > 0
    ),
    "parameter_relevant_content_gate": (
        len(
            parameter_relevant_files
        )
        > 0
    ),
    "archive_inventory_created_gate": (
        ARCHIVE_INVENTORY_CSV.is_file()
        and ARCHIVE_INVENTORY_CSV.stat().st_size
        > 0
    ),
    "extracted_inventory_created_gate": (
        EXTRACTED_FILE_INVENTORY_CSV.is_file()
        and EXTRACTED_FILE_INVENTORY_CSV.stat().st_size
        > 0
    ),
    "no_parameter_adopted_gate": True,
    "no_accepted_topology_modified_gate": True,
    "no_accepted_coordinates_modified_gate": True,
    "force_field_adoption_blocked_gate": True,
    "MD_execution_blocked_gate": True,
}

for gate_name, value in gates.items():
    print(
        f"{gate_name}="
        f"{'PASS' if value else 'FAIL'}"
    )

all_gates_pass = all(
    gates.values()
)


print("\n[10] WRITE REPORT")

decision = (
    "D040_A7A_SUPPLEMENTARY_ARCHIVE_RECOGNITION_PASS_"
    "FULL_PARAMETER_CONTENT_AUDIT_AUTHORIZED"
    if all_gates_pass
    else
    "D040_A7A_SUPPLEMENTARY_ARCHIVE_RECOGNITION_"
    "REVIEW_REQUIRED"
)

report = {
    "generated_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "decision": decision,
    "source_identity": {
        "canonical_article": {
            "path": relative(
                CANONICAL_ARTICLE
            ),
            "bytes": int(
                CANONICAL_ARTICLE.stat().st_size
            ),
            "sha256": sha256(
                CANONICAL_ARTICLE
            ),
        },
        "official_supplementary_zip": {
            "path": relative(
                OFFICIAL_SI_ZIP
            ),
            "bytes": int(
                OFFICIAL_SI_ZIP.stat().st_size
            ),
            "sha256": sha256(
                OFFICIAL_SI_ZIP
            ),
        },
        "nested_parameter_tar": {
            "path": relative(
                nested_tar_path
            ),
            "bytes": int(
                nested_tar_path.stat().st_size
            ),
            "sha256": sha256(
                nested_tar_path
            ),
        },
    },
    "A7_diagnostic": {
        "path": relative(
            A7_REPORT
        ),
        "decision": (
            a7_decision
        ),
        "interpretation": (
            "PRESERVED_AS_ARTIFACT_RECOGNITION_LOGIC_DIAGNOSTIC"
        ),
    },
    "artifact_model": {
        "separate_supplementary_PDF_required": False,
        "official_supplementary_artifact": (
            "ghorai2025_SI.zip"
        ),
        "nested_parameter_archive": (
            "hBN.tar"
        ),
        "supplementary_artifact_complete": (
            all_gates_pass
        ),
    },
    "content_counts": {
        "zip_member_count": len(
            zip_members
        ),
        "tar_member_count": len(
            tar_members
        ),
        "zip_extracted_file_count": len(
            zip_extracted_files
        ),
        "tar_extracted_file_count": len(
            tar_extracted_files
        ),
        "parameter_relevant_file_count": len(
            parameter_relevant_files
        ),
    },
    "parameter_relevant_files": [
        relative(path)
        for path in parameter_relevant_files
    ],
    "gates": gates,
    "authorizations": {
        "supplementary_artifact_recognition_completed": (
            all_gates_pass
        ),
        "full_parameter_content_audit_authorized": (
            all_gates_pass
        ),
        "parameter_translation_review_authorized": False,
        "parameter_comparison_authorized": False,
        "parameter_adoption_authorized": False,
        "new_atom_type_definition_authorized": False,
        "bonded_parameter_modification_authorized": False,
        "hydrogen_coordinate_insertion_authorized": False,
        "charge_to_topology_mapping_execution_authorized": False,
        "topology_modification_authorized": False,
        "force_field_adoption_authorized": False,
        "energy_execution_authorized": False,
        "minimization_execution_authorized": False,
        "validation_MD_execution_authorized": False,
        "production_MD_authorized": False,
    },
    "next_required_block": {
        "name": (
            "D040_A8_FULL_PARAMETER_CONTENT_AUDIT"
        ),
        "required_actions": [
            (
                "Parse the extracted LAMMPS and supporting files "
                "without executing them."
            ),
            (
                "Extract HN and HB atom definitions, nonbonded terms, "
                "charges, B-H and N-H bonds, angles, dihedrals, "
                "impropers, exclusions and combining rules."
            ),
            (
                "Record source units and functional forms before any "
                "GROMACS translation."
            ),
            (
                "Keep all parameter adoption and topology modification "
                "blocked."
            ),
        ],
    },
}

REPORT_JSON.write_text(
    json.dumps(
        report,
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)

print(
    f"archive_inventory_csv="
    f"{ARCHIVE_INVENTORY_CSV}"
)
print(
    f"archive_inventory_csv_sha256="
    f"{sha256(ARCHIVE_INVENTORY_CSV)}"
)
print(
    f"extracted_file_inventory_csv="
    f"{EXTRACTED_FILE_INVENTORY_CSV}"
)
print(
    f"extracted_file_inventory_csv_sha256="
    f"{sha256(EXTRACTED_FILE_INVENTORY_CSV)}"
)
print(
    f"report_json={REPORT_JSON}"
)
print(
    f"report_json_sha256="
    f"{sha256(REPORT_JSON)}"
)


print("\n[11] DECISION")

print(f"decision={decision}")
print(
    "supplementary_artifact_recognition_completed="
    f"{all_gates_pass}"
)
print(
    "full_parameter_content_audit_authorized="
    f"{all_gates_pass}"
)
print(
    "parameter_translation_review_authorized=False"
)
print(
    "parameter_adoption_authorized=False"
)
print(
    "topology_modification_authorized=False"
)
print(
    "force_field_adoption_authorized=False"
)
print(
    "validation_MD_execution_authorized=False"
)
print(
    "production_MD_authorized=False"
)
print("=" * 100)

if not all_gates_pass:
    raise SystemExit(2)
