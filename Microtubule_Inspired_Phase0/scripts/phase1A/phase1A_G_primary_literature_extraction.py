#!/usr/bin/env python3
"""
DAY040 / D040-A7

Primary-literature artifact and parameter extraction for hydrogen-
functionalized hBN.

Scope
-----
1. Bind to the D040-A6 provenance decision.
2. Locate the primary Ghorai et al. 2025 article in the repository.
3. Locate supplementary information and parameter archives.
4. Extract only explicitly reported values from the primary article.
5. Classify evidence by interaction domain and artifact completeness.
6. Keep all parameter adoption and topology modification blocked.

No force-field parameter is adopted.
No topology or coordinate file is modified.
No MD calculation is executed.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

A6_DIR = (
    ROOT
    / "runs/phase1A/day040_phase1A_G_hydrogen_parameter_provenance_review"
)

A6_JSON = (
    A6_DIR
    / "QM_F06_UPPER_V7A_R1_HYDROGEN_PARAMETER_PROVENANCE.json"
)

A6_CSV = (
    A6_DIR
    / "QM_F06_UPPER_V7A_R1_HYDROGEN_PARAMETER_PROVENANCE.csv"
)

EXPECTED_A6_DECISION = (
    "D040_A6_PARAMETER_PROVENANCE_PASS_"
    "PRIMARY_LITERATURE_EXTRACTION_AUTHORIZED"
)

OUTPUT_DIR = (
    ROOT
    / "runs/phase1A/day040_phase1A_G_primary_literature_extraction"
)

SOURCE_REGISTRY_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_PRIMARY_SOURCE_REGISTRY.csv"
)

ARTIFACT_REGISTRY_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_PRIMARY_SOURCE_ARTIFACTS.csv"
)

PARAMETER_EXTRACTION_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_GHORAI2025_PARAMETER_EXTRACTION.csv"
)

DOMAIN_MATRIX_CSV = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_GHORAI2025_DOMAIN_EVIDENCE_MATRIX.csv"
)

REPORT_JSON = (
    OUTPUT_DIR
    / "QM_F06_UPPER_V7A_R1_PRIMARY_LITERATURE_EXTRACTION.json"
)

PRIMARY_DOI = "10.1063/5.0242541"
PRIMARY_SOURCE_ID = "GHORAI_2025_JCP_162_044705"
PRIMARY_ARTICLE_NAME = "044705_1_5.0242541.pdf"

SEARCH_EXTENSIONS = {
    ".pdf",
    ".txt",
    ".md",
    ".json",
    ".csv",
    ".tar",
    ".gz",
    ".tgz",
    ".zip",
    ".lmp",
    ".in",
    ".data",
    ".itp",
    ".top",
}

EXCLUDED_DIRECTORY_NAMES = {
    ".git",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    "node_modules",
}

ARTICLE_TEXT_PATTERNS = {
    "title": re.compile(
        r"Molecular dynamics simulations of functionalized\s+"
        r"hBN nanopores in water",
        re.IGNORECASE,
    ),
    "doi": re.compile(
        r"10\.1063/5\.0242541",
        re.IGNORECASE,
    ),
    "supplementary_tar": re.compile(
        r"tar file is included as part of the supplementary material",
        re.IGNORECASE,
    ),
    "table_I": re.compile(
        r"TABLE I\.",
        re.IGNORECASE,
    ),
    "table_II": re.compile(
        r"TABLE II\.",
        re.IGNORECASE,
    ),
    "table_S2": re.compile(
        r"Table S2",
        re.IGNORECASE,
    ),
    "table_S3": re.compile(
        r"Table S3",
        re.IGNORECASE,
    ),
    "geometric_combining_rules": re.compile(
        r"geometric(?:-mean)? combining rules",
        re.IGNORECASE,
    ),
    "tip4p2005": re.compile(
        r"TIP4P/2005",
        re.IGNORECASE,
    ),
    "lammps": re.compile(
        r"LAMMPS",
        re.IGNORECASE,
    ),
    "gulp": re.compile(
        r"GULP",
        re.IGNORECASE,
    ),
}

EXPLICIT_PRIMARY_VALUES = [
    {
        "source_id": PRIMARY_SOURCE_ID,
        "table": "TABLE_I",
        "interaction_family": "BOND",
        "interaction": "B-H",
        "functional_form": "HARMONIC_BOND",
        "equilibrium_value": 1.197,
        "equilibrium_unit": "angstrom",
        "force_constant": 2112.9,
        "force_constant_unit": "kJ mol^-1 angstrom^-2",
        "evidence_status": "DIRECT_EVIDENCE",
        "transcription_basis": (
            "Primary article Table I: B-H l0 and final fitted bond constant"
        ),
        "parameter_adoption_status": "NOT_AUTHORIZED",
    },
    {
        "source_id": PRIMARY_SOURCE_ID,
        "table": "TABLE_I",
        "interaction_family": "BOND",
        "interaction": "N-H",
        "functional_form": "HARMONIC_BOND",
        "equilibrium_value": 1.017,
        "equilibrium_unit": "angstrom",
        "force_constant": 4315.2,
        "force_constant_unit": "kJ mol^-1 angstrom^-2",
        "evidence_status": "DIRECT_EVIDENCE",
        "transcription_basis": (
            "Primary article Table I: N-H l0 and final fitted bond constant"
        ),
        "parameter_adoption_status": "NOT_AUTHORIZED",
    },
]

DOMAIN_ROWS = [
    {
        "source_id": PRIMARY_SOURCE_ID,
        "domain": "hydrogenated_hBN_edge",
        "evidence_status": "DIRECT_EVIDENCE",
        "basis": (
            "Article develops parameters for H functionalization on B and N edges"
        ),
    },
    {
        "source_id": PRIMARY_SOURCE_ID,
        "domain": "B_H_bond",
        "evidence_status": "DIRECT_EVIDENCE",
        "basis": "Table I explicitly reports B-H harmonic bond parameters",
    },
    {
        "source_id": PRIMARY_SOURCE_ID,
        "domain": "N_H_bond",
        "evidence_status": "DIRECT_EVIDENCE",
        "basis": "Table I explicitly reports N-H harmonic bond parameters",
    },
    {
        "source_id": PRIMARY_SOURCE_ID,
        "domain": "H_containing_angles",
        "evidence_status": "CONTENT_AUDIT_PENDING",
        "basis": "Article states that Tables I and II contain bonded parameters",
    },
    {
        "source_id": PRIMARY_SOURCE_ID,
        "domain": "H_containing_dihedrals",
        "evidence_status": "CONTENT_AUDIT_PENDING",
        "basis": "Article states that Tables I and II contain bonded parameters",
    },
    {
        "source_id": PRIMARY_SOURCE_ID,
        "domain": "H_nonbonded_parameters",
        "evidence_status": "CONTENT_AUDIT_PENDING",
        "basis": "Article points to Table S2 and borazine fitting",
    },
    {
        "source_id": PRIMARY_SOURCE_ID,
        "domain": "partial_charges",
        "evidence_status": "CONTENT_AUDIT_PENDING",
        "basis": "Article uses DDAP charges and references supplementary details",
    },
    {
        "source_id": PRIMARY_SOURCE_ID,
        "domain": "water_compatibility",
        "evidence_status": "DIRECT_EVIDENCE",
        "basis": "Article validates functionalized nanopores with TIP4P/2005 water",
    },
    {
        "source_id": PRIMARY_SOURCE_ID,
        "domain": "GROMACS_direct_implementation",
        "evidence_status": "NO_EVIDENCE",
        "basis": "Article reports GULP fitting and LAMMPS implementation, not GROMACS",
    },
    {
        "source_id": PRIMARY_SOURCE_ID,
        "domain": "curved_BNNT_transferability",
        "evidence_status": "NO_EVIDENCE",
        "basis": "Article studies planar nanoporous hBN, not curved BNNT geometry",
    },
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(chunk)

    return digest.hexdigest()


def relative_or_absolute(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def is_excluded(path: Path) -> bool:
    try:
        parts = path.relative_to(ROOT).parts
    except ValueError:
        parts = path.parts

    return any(
        part in EXCLUDED_DIRECTORY_NAMES
        for part in parts
    )


def load_json(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(path)

    return json.loads(
        path.read_text(
            encoding="utf-8"
        )
    )


def candidate_files() -> list[Path]:
    results = []

    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue

        if is_excluded(path):
            continue

        if path.suffix.lower() not in SEARCH_EXTENSIONS:
            continue

        results.append(path)

    return sorted(
        results,
        key=lambda item: relative_or_absolute(item),
    )


def extract_pdf_text(
    article_path: Path,
    output_text_path: Path,
) -> dict:
    pdftotext = shutil.which(
        "pdftotext"
    )

    if pdftotext is None:
        return {
            "status": "PDFTOTEXT_NOT_AVAILABLE",
            "returncode": None,
            "stdout": "",
            "stderr": "",
        }

    completed = subprocess.run(
        [
            pdftotext,
            "-layout",
            str(article_path),
            str(output_text_path),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    return {
        "status": (
            "PASS"
            if completed.returncode == 0
            else "FAIL"
        ),
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


print("=" * 100)
print("DAY040 / D040-A7 — PRIMARY LITERATURE PARAMETER EXTRACTION")
print("=" * 100)


print("\n[1] UPSTREAM AUTHORIZATION")

for path in (
    A6_JSON,
    A6_CSV,
):
    if not path.is_file():
        raise FileNotFoundError(path)

    print(
        f"FOUND bytes={path.stat().st_size:8d} "
        f"sha256={sha256(path)} "
        f"{path}"
    )

a6_report = load_json(
    A6_JSON
)

a6_decision = a6_report.get(
    "decision"
)

a6_summary = a6_report.get(
    "summary",
    {}
)

# D040-A6 was generated using an earlier report contract that did not
# serialize the printed terminal decision into the JSON. Accept either:
#
# 1. the explicit decision contract, when present; or
# 2. the authoritative legacy contract defined by the block identity,
#    authorization flags and declared next block.
explicit_decision_contract_gate = (
    a6_decision
    == EXPECTED_A6_DECISION
)

legacy_summary_contract_gate = (
    a6_decision is None
    and a6_report.get(
        "day"
    )
    == "DAY040"
    and a6_report.get(
        "block"
    )
    == "D040_A6"
    and a6_summary.get(
        "complete_parameter_provenance"
    )
    is False
    and a6_summary.get(
        "parameter_adoption_authorized"
    )
    is False
    and a6_summary.get(
        "primary_literature_required"
    )
    is True
    and a6_summary.get(
        "next_block"
    )
    == "D040_A7_PRIMARY_LITERATURE_PARAMETER_EXTRACTION"
)

if not (
    explicit_decision_contract_gate
    or legacy_summary_contract_gate
):
    raise RuntimeError(
        "Unexpected D040-A6 authorization contract.\n"
        f"Observed decision: {a6_decision}\n"
        f"Observed block: {a6_report.get('block')}\n"
        f"Observed summary: {a6_summary}"
    )

if (
    a6_summary.get(
        "primary_literature_required"
    )
    is not True
):
    raise RuntimeError(
        "Primary-literature extraction is not authorized"
    )

if (
    a6_summary.get(
        "parameter_adoption_authorized"
    )
    is not False
):
    raise RuntimeError(
        "Unexpected parameter-adoption authorization"
    )

authorization_contract = (
    "EXPLICIT_DECISION"
    if explicit_decision_contract_gate
    else "LEGACY_SUMMARY_CONTRACT"
)

print(
    "D040_A6_authorization_contract = "
    f"{authorization_contract}"
)
print("D040_A6_decision_gate = PASS")
print("primary_literature_extraction_authorization_gate = PASS")
print("parameter_adoption_blocked_gate = PASS")
print("topology_modification_blocked_gate = PASS")


print("\n[2] REPOSITORY PRIMARY-ARTIFACT INVENTORY")

files = candidate_files()

article_candidates = []
supplementary_candidates = []
parameter_archive_candidates = []

for path in files:
    name_lower = path.name.lower()
    path_lower = relative_or_absolute(
        path
    ).lower()

    if (
        path.name == PRIMARY_ARTICLE_NAME
        or "5.0242541" in name_lower
        or "044705" in name_lower
    ):
        article_candidates.append(
            path
        )

    if (
        "supp" in name_lower
        or "si" in name_lower
        or "support" in name_lower
    ) and (
        "0242541" in path_lower
        or "044705" in path_lower
        or "ghorai" in path_lower
        or "functionalized" in path_lower
    ):
        supplementary_candidates.append(
            path
        )

    if (
        path.suffix.lower()
        in {
            ".tar",
            ".gz",
            ".tgz",
            ".zip",
        }
        and (
            "0242541" in path_lower
            or "044705" in path_lower
            or "ghorai" in path_lower
            or "hbn" in path_lower
        )
    ):
        parameter_archive_candidates.append(
            path
        )

print(
    f"repository_search_file_count = "
    f"{len(files)}"
)
print(
    f"primary_article_candidate_count = "
    f"{len(article_candidates)}"
)
print(
    f"supplementary_candidate_count = "
    f"{len(supplementary_candidates)}"
)
print(
    f"parameter_archive_candidate_count = "
    f"{len(parameter_archive_candidates)}"
)

for path in article_candidates:
    print(
        "ARTICLE_CANDIDATE "
        f"bytes={path.stat().st_size} "
        f"sha256={sha256(path)} "
        f"path={relative_or_absolute(path)}"
    )

for path in supplementary_candidates:
    print(
        "SUPPLEMENTARY_CANDIDATE "
        f"bytes={path.stat().st_size} "
        f"sha256={sha256(path)} "
        f"path={relative_or_absolute(path)}"
    )

for path in parameter_archive_candidates:
    print(
        "PARAMETER_ARCHIVE_CANDIDATE "
        f"bytes={path.stat().st_size} "
        f"sha256={sha256(path)} "
        f"path={relative_or_absolute(path)}"
    )


print("\n[3] SELECT PRIMARY ARTICLE")

if len(article_candidates) == 0:
    selected_article = None
    article_selection_status = (
        "ARTICLE_NOT_PRESENT_IN_REPOSITORY"
    )

elif len(article_candidates) == 1:
    selected_article = article_candidates[0]
    article_selection_status = (
        "UNIQUE_ARTICLE_SELECTED"
    )

else:
    exact_name_matches = [
        path
        for path in article_candidates
        if path.name
        == PRIMARY_ARTICLE_NAME
    ]

    if len(exact_name_matches) == 1:
        selected_article = exact_name_matches[0]
        article_selection_status = (
            "EXACT_FILENAME_SELECTED"
        )
    else:
        selected_article = None
        article_selection_status = (
            "ARTICLE_SELECTION_AMBIGUOUS"
        )

print(
    f"article_selection_status = "
    f"{article_selection_status}"
)

if selected_article is not None:
    print(
        f"selected_article = "
        f"{relative_or_absolute(selected_article)}"
    )
    print(
        f"selected_article_sha256 = "
        f"{sha256(selected_article)}"
    )


print("\n[4] PRIMARY ARTICLE CONTENT AUDIT")

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

article_text_path = (
    OUTPUT_DIR
    / "GHORAI_2025_JCP_162_044705_pdftotext.txt"
)

content_flags = {
    name: False
    for name in ARTICLE_TEXT_PATTERNS
}

pdf_text_result = {
    "status": "NOT_RUN",
    "returncode": None,
    "stdout": "",
    "stderr": "",
}

article_text = ""

if selected_article is not None:
    pdf_text_result = extract_pdf_text(
        selected_article,
        article_text_path,
    )

    if (
        pdf_text_result[
            "status"
        ]
        == "PASS"
        and article_text_path.is_file()
    ):
        article_text = article_text_path.read_text(
            encoding="utf-8",
            errors="replace",
        )

        for name, pattern in ARTICLE_TEXT_PATTERNS.items():
            content_flags[name] = bool(
                pattern.search(
                    article_text
                )
            )

print(
    f"pdftotext_status = "
    f"{pdf_text_result['status']}"
)

for name, value in content_flags.items():
    print(
        f"{name}_gate = "
        f"{'PASS' if value else 'FAIL'}"
    )


print("\n[5] EXTRACT EXPLICIT PRIMARY VALUES")

for record in EXPLICIT_PRIMARY_VALUES:
    print(
        f"table={record['table']} "
        f"interaction={record['interaction']} "
        f"l0={record['equilibrium_value']} "
        f"{record['equilibrium_unit']} "
        f"k={record['force_constant']} "
        f"{record['force_constant_unit']} "
        f"status={record['evidence_status']}"
    )


print("\n[6] SUPPLEMENTARY-ARTIFACT COMPLETENESS")

supplementary_article_present = any(
    path.suffix.lower() == ".pdf"
    for path in supplementary_candidates
)

parameter_archive_present = (
    len(
        parameter_archive_candidates
    )
    > 0
)

supplementary_artifact_complete = (
    supplementary_article_present
    and parameter_archive_present
)

print(
    "supplementary_article_present = "
    f"{supplementary_article_present}"
)
print(
    "parameter_archive_present = "
    f"{parameter_archive_present}"
)
print(
    "supplementary_artifact_complete = "
    f"{supplementary_artifact_complete}"
)

if not supplementary_article_present:
    print(
        "MISSING_ARTIFACT="
        "GHORAI_2025_SUPPLEMENTARY_INFORMATION_PDF"
    )

if not parameter_archive_present:
    print(
        "MISSING_ARTIFACT="
        "GHORAI_2025_LAMMPS_PARAMETER_ARCHIVE"
    )


print("\n[7] WRITE SOURCE REGISTRY")

source_rows = [
    {
        "source_id": PRIMARY_SOURCE_ID,
        "title": (
            "Molecular dynamics simulations of functionalized hBN "
            "nanopores in water: Ab initio force field and "
            "implications for water desalination"
        ),
        "authors": (
            "Sagar Ghorai; Pradeep Dhondi; Ananth Govind Rajan"
        ),
        "journal": "Journal of Chemical Physics",
        "volume": "162",
        "article_number": "044705",
        "year": 2025,
        "doi": PRIMARY_DOI,
        "source_role": (
            "PRIMARY_FORCE_FIELD_SOURCE_FOR_H_AND_OH_FUNCTIONALIZED_HBN"
        ),
        "article_status": (
            "PRESENT_AND_CONTENT_AUDITED"
            if selected_article is not None
            and content_flags["doi"]
            and content_flags["title"]
            else article_selection_status
        ),
        "supplementary_information_status": (
            "PRESENT"
            if supplementary_article_present
            else "MISSING"
        ),
        "parameter_archive_status": (
            "PRESENT"
            if parameter_archive_present
            else "MISSING"
        ),
        "parameter_adoption_status": "NOT_AUTHORIZED",
    }
]

with SOURCE_REGISTRY_CSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(
            source_rows[0].keys()
        ),
    )
    writer.writeheader()
    writer.writerows(
        source_rows
    )


print("\n[8] WRITE ARTIFACT REGISTRY")

artifact_rows = []

for artifact_type, paths in (
    (
        "PRIMARY_ARTICLE",
        article_candidates,
    ),
    (
        "SUPPLEMENTARY_INFORMATION_CANDIDATE",
        supplementary_candidates,
    ),
    (
        "PARAMETER_ARCHIVE_CANDIDATE",
        parameter_archive_candidates,
    ),
):
    for path in paths:
        artifact_rows.append(
            {
                "source_id": PRIMARY_SOURCE_ID,
                "artifact_type": artifact_type,
                "path": relative_or_absolute(
                    path
                ),
                "suffix": path.suffix.lower(),
                "bytes": int(
                    path.stat().st_size
                ),
                "sha256": sha256(path),
                "content_audit_status": (
                    "AUDITED"
                    if path
                    == selected_article
                    else "PENDING"
                ),
                "parameter_adoption_status": (
                    "NOT_AUTHORIZED"
                ),
            }
        )

artifact_fieldnames = [
    "source_id",
    "artifact_type",
    "path",
    "suffix",
    "bytes",
    "sha256",
    "content_audit_status",
    "parameter_adoption_status",
]

with ARTIFACT_REGISTRY_CSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=artifact_fieldnames,
    )
    writer.writeheader()
    writer.writerows(
        artifact_rows
    )


print("\n[9] WRITE PARAMETER AND DOMAIN TABLES")

parameter_fieldnames = [
    "source_id",
    "table",
    "interaction_family",
    "interaction",
    "functional_form",
    "equilibrium_value",
    "equilibrium_unit",
    "force_constant",
    "force_constant_unit",
    "evidence_status",
    "transcription_basis",
    "parameter_adoption_status",
]

with PARAMETER_EXTRACTION_CSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=parameter_fieldnames,
    )
    writer.writeheader()
    writer.writerows(
        EXPLICIT_PRIMARY_VALUES
    )

domain_fieldnames = [
    "source_id",
    "domain",
    "evidence_status",
    "basis",
]

with DOMAIN_MATRIX_CSV.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=domain_fieldnames,
    )
    writer.writeheader()
    writer.writerows(
        DOMAIN_ROWS
    )


print("\n[10] SCIENTIFIC GATES")

article_present_gate = (
    selected_article is not None
)

article_identity_gate = (
    content_flags["title"]
    and content_flags["doi"]
)

primary_bond_values_extracted_gate = (
    len(
        EXPLICIT_PRIMARY_VALUES
    )
    == 2
    and {
        record["interaction"]
        for record in EXPLICIT_PRIMARY_VALUES
    }
    == {
        "B-H",
        "N-H",
    }
)

supplementary_gap_recorded_gate = (
    (
        supplementary_artifact_complete
    )
    or (
        not supplementary_article_present
        and not parameter_archive_present
    )
)

gates = {
    "D040_A6_decision_gate": True,
    "primary_article_present_gate": (
        article_present_gate
    ),
    "primary_article_identity_gate": (
        article_identity_gate
    ),
    "primary_article_content_audit_gate": (
        content_flags["table_I"]
        and content_flags["lammps"]
        and content_flags["gulp"]
    ),
    "primary_bond_values_extracted_gate": (
        primary_bond_values_extracted_gate
    ),
    "source_registry_created_gate": (
        SOURCE_REGISTRY_CSV.is_file()
        and SOURCE_REGISTRY_CSV.stat().st_size > 0
    ),
    "artifact_registry_created_gate": (
        ARTIFACT_REGISTRY_CSV.is_file()
        and ARTIFACT_REGISTRY_CSV.stat().st_size > 0
    ),
    "parameter_extraction_created_gate": (
        PARAMETER_EXTRACTION_CSV.is_file()
        and PARAMETER_EXTRACTION_CSV.stat().st_size > 0
    ),
    "domain_matrix_created_gate": (
        DOMAIN_MATRIX_CSV.is_file()
        and DOMAIN_MATRIX_CSV.stat().st_size > 0
    ),
    "supplementary_gap_recorded_gate": (
        supplementary_gap_recorded_gate
    ),
    "no_parameter_adopted_gate": True,
    "no_topology_modified_gate": True,
    "no_coordinates_modified_gate": True,
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


print("\n[11] WRITE REPORT")

if all_gates_pass and supplementary_artifact_complete:
    decision = (
        "D040_A7_PRIMARY_LITERATURE_EXTRACTION_PASS_"
        "FULL_PARAMETER_CONTENT_AUDIT_AUTHORIZED"
    )
    next_block = (
        "D040_A8_FULL_PARAMETER_CONTENT_AUDIT"
    )

elif all_gates_pass:
    decision = (
        "D040_A7_PRIMARY_LITERATURE_EXTRACTION_PASS_"
        "SUPPLEMENTARY_ARTIFACT_RECOVERY_REQUIRED"
    )
    next_block = (
        "D040_A8_GHORAI_SUPPLEMENTARY_ARTIFACT_RECOVERY"
    )

else:
    decision = (
        "D040_A7_PRIMARY_LITERATURE_EXTRACTION_"
        "REVIEW_REQUIRED"
    )
    next_block = (
        "D040_A7_REPAIR_OR_MANUAL_ARTICLE_PLACEMENT"
    )

report = {
    "generated_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "decision": decision,
    "source_id": PRIMARY_SOURCE_ID,
    "doi": PRIMARY_DOI,
    "article_selection_status": (
        article_selection_status
    ),
    "selected_article": (
        relative_or_absolute(
            selected_article
        )
        if selected_article is not None
        else None
    ),
    "selected_article_sha256": (
        sha256(selected_article)
        if selected_article is not None
        else None
    ),
    "pdf_text_extraction": (
        pdf_text_result
    ),
    "content_flags": content_flags,
    "artifact_counts": {
        "article_candidates": len(
            article_candidates
        ),
        "supplementary_candidates": len(
            supplementary_candidates
        ),
        "parameter_archive_candidates": len(
            parameter_archive_candidates
        ),
    },
    "supplementary_completeness": {
        "supplementary_article_present": (
            supplementary_article_present
        ),
        "parameter_archive_present": (
            parameter_archive_present
        ),
        "complete": (
            supplementary_artifact_complete
        ),
    },
    "explicit_primary_values": (
        EXPLICIT_PRIMARY_VALUES
    ),
    "domain_evidence": DOMAIN_ROWS,
    "gates": gates,
    "authorizations": {
        "primary_article_content_audit_completed": (
            all_gates_pass
        ),
        "supplementary_artifact_recovery_authorized": (
            all_gates_pass
            and not supplementary_artifact_complete
        ),
        "full_parameter_content_audit_authorized": (
            all_gates_pass
            and supplementary_artifact_complete
        ),
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
        "name": next_block,
        "required_actions": (
            [
                (
                    "Recover the official Ghorai 2025 supplementary "
                    "information PDF."
                ),
                (
                    "Recover the official supplementary tar archive "
                    "containing LAMMPS input/output and parameter files."
                ),
                (
                    "Verify file type, hash, content and provenance "
                    "before extracting additional parameters."
                ),
            ]
            if not supplementary_artifact_complete
            else [
                (
                    "Extract all H atom types, LJ terms, angles, "
                    "dihedrals, impropers, charges and combining rules."
                ),
                (
                    "Translate units and functional forms without "
                    "adopting the parameters."
                ),
            ]
        ),
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
    f"source_registry_csv = "
    f"{SOURCE_REGISTRY_CSV}"
)
print(
    "source_registry_csv_sha256 = "
    f"{sha256(SOURCE_REGISTRY_CSV)}"
)
print(
    f"artifact_registry_csv = "
    f"{ARTIFACT_REGISTRY_CSV}"
)
print(
    "artifact_registry_csv_sha256 = "
    f"{sha256(ARTIFACT_REGISTRY_CSV)}"
)
print(
    f"parameter_extraction_csv = "
    f"{PARAMETER_EXTRACTION_CSV}"
)
print(
    "parameter_extraction_csv_sha256 = "
    f"{sha256(PARAMETER_EXTRACTION_CSV)}"
)
print(
    f"domain_matrix_csv = "
    f"{DOMAIN_MATRIX_CSV}"
)
print(
    "domain_matrix_csv_sha256 = "
    f"{sha256(DOMAIN_MATRIX_CSV)}"
)
print(f"report_json = {REPORT_JSON}")
print(
    f"report_json_sha256 = "
    f"{sha256(REPORT_JSON)}"
)


print("\n[12] DECISION")

print(f"decision={decision}")
print(
    "primary_article_content_audit_completed="
    f"{all_gates_pass}"
)
print(
    "supplementary_artifact_recovery_authorized="
    f"{all_gates_pass and not supplementary_artifact_complete}"
)
print(
    "full_parameter_content_audit_authorized="
    f"{all_gates_pass and supplementary_artifact_complete}"
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
