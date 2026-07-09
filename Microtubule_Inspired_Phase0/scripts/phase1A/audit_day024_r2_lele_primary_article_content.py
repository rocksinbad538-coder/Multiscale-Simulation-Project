#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/phase1A/day024_chemical_end_rim_design"

GATE3R4 = (
    BASE
    / "35_r2_p1_source_artifact_retrieval"
)

PDF = (
    GATE3R4
    / "raw"
    / "lele_2022_reaxff_article.pdf"
)

RETRIEVAL_LEDGER = (
    GATE3R4
    / "r2_p1_artifact_retrieval_ledger.csv"
)

RETRIEVAL_SUMMARY = (
    GATE3R4
    / "r2_p1_artifact_retrieval_summary.csv"
)

OUT = (
    BASE
    / "36_r2_lele_primary_article_content_audit"
)

PAGE_TEXT = (
    OUT
    / "lele_2022_article_text_by_page.txt"
)

EVIDENCE = (
    OUT
    / "lele_2022_article_evidence.csv"
)

DOMAIN_MATRIX = (
    OUT
    / "lele_2022_R2_domain_assessment.csv"
)

REQUIRED_ARTIFACTS = (
    OUT
    / "lele_2022_required_followup_artifacts.csv"
)

SUMMARY = (
    OUT
    / "lele_2022_primary_article_content_audit_summary.csv"
)

GATES = (
    OUT
    / "lele_2022_primary_article_content_audit_gates.csv"
)

JSON_OUT = (
    OUT
    / "lele_2022_primary_article_content_audit.json"
)

MANIFEST = (
    OUT
    / "lele_2022_primary_article_content_audit_manifest.csv"
)

REPORT = (
    OUT
    / "LELE_2022_PRIMARY_ARTICLE_CONTENT_AUDIT_DAY024.md"
)

EXPECTED_RETRIEVAL_DECISION = (
    "R2_P1_SOURCE_ARTIFACT_RETRIEVAL_AUDITED"
)

EXPECTED_ARTIFACT_ID = (
    "LELE_ARTICLE_OPEN_COPY"
)

EXPECTED_DOI = (
    "10.1021/acs.jpca.1c09648"
)

PASS_DECISION = (
    "R2_LELE_PRIMARY_ARTICLE_CONTENT_AUDITED_"
    "SI_AND_PARAMETER_FILES_STILL_REQUIRED"
)

REVIEW_DECISION = (
    "R2_LELE_PRIMARY_ARTICLE_CONTENT_AUDIT_REQUIRES_REVIEW"
)


PATTERNS = {
    "ARTICLE_IDENTITY": [
        r"10\.1021/acs\.jpca\.1c09648",
        r"ReaxFF Force Field Development",
        r"gas[- ]phase",
        r"hBN nanostructure synthesis",
    ],
    "REAXFF_METHOD": [
        r"\bReaxFF\b",
        r"reactive force field",
        r"bond order",
        r"charge equilibration",
        r"\bQEq\b",
    ],
    "ELEMENT_DOMAIN": [
        r"\bboron\b",
        r"\bnitrogen\b",
        r"\bhydrogen\b",
        r"\bB/N/H\b",
        r"\bB[-–]H\b",
        r"\bN[-–]H\b",
        r"\bB[-–]N\b",
    ],
    "TRAINING_SET": [
        r"training set",
        r"training data",
        r"quantum mechanical",
        r"\bDFT\b",
        r"heat of formation",
        r"bond dissociation",
        r"reaction energ",
        r"angle scan",
        r"torsion",
    ],
    "BORAZINE_AND_PRECURSORS": [
        r"borazine",
        r"aminoborane",
        r"ammonia borane",
        r"boron nitride precursor",
    ],
    "PARAMETER_AVAILABILITY": [
        r"Supporting Information",
        r"supporting information",
        r"parameter file",
        r"force field parameters",
        r"available.*supporting",
    ],
    "VALIDATION": [
        r"validation",
        r"validated",
        r"molecular dynamics",
        r"nanostructure",
        r"nucleation",
        r"growth",
        r"synthesis",
    ],
    "WATER_OR_AQUEOUS": [
        r"\bwater\b",
        r"aqueous",
        r"hydration",
        r"solvation",
        r"liquid water",
    ],
    "BNNT_OR_CURVATURE": [
        r"boron nitride nanotube",
        r"\bBNNT\b",
        r"nanotube",
        r"curvature",
        r"curved",
    ],
    "EDGE_OR_DEFECT": [
        r"\bedge\b",
        r"defect",
        r"vacancy",
        r"termination",
        r"hydrogen[- ]terminated",
        r"functionalized",
    ],
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty required file: {path}"
        )


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
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


def normalize_space(text: str) -> str:
    return re.sub(
        r"\s+",
        " ",
        text,
    ).strip()


def extract_with_pymupdf(
    path: Path,
) -> list[str] | None:
    try:
        import fitz
    except ImportError:
        return None

    document = fitz.open(
        path
    )

    pages = []

    try:
        for page in document:
            pages.append(
                page.get_text(
                    "text"
                )
            )
    finally:
        document.close()

    return pages


def extract_with_pypdf(
    path: Path,
) -> list[str] | None:
    try:
        from pypdf import PdfReader
    except ImportError:
        return None

    reader = PdfReader(
        str(path)
    )

    return [
        page.extract_text() or ""
        for page in reader.pages
    ]


def extract_with_pdftotext(
    path: Path,
) -> list[str] | None:
    executable = shutil.which(
        "pdftotext"
    )

    if executable is None:
        return None

    with tempfile.TemporaryDirectory() as temporary:
        output = (
            Path(temporary)
            / "article.txt"
        )

        process = subprocess.run(
            [
                executable,
                "-layout",
                str(path),
                str(output),
            ],
            text=True,
            capture_output=True,
            check=False,
        )

        if process.returncode != 0:
            return None

        text = output.read_text(
            encoding="utf-8",
            errors="replace",
        )

    pages = text.split(
        "\f"
    )

    while (
        pages
        and not pages[-1].strip()
    ):
        pages.pop()

    return pages


def extract_pages(
    path: Path,
) -> tuple[list[str], str]:
    for method_name, extractor in (
        (
            "PYMUPDF",
            extract_with_pymupdf,
        ),
        (
            "PYPDF",
            extract_with_pypdf,
        ),
        (
            "PDFTOTEXT",
            extract_with_pdftotext,
        ),
    ):
        pages = extractor(
            path
        )

        if (
            pages is not None
            and len(pages) > 0
            and sum(
                len(page.strip())
                for page in pages
            )
            > 1000
        ):
            return (
                pages,
                method_name,
            )

    raise RuntimeError(
        "Could not extract sufficient text from the PDF. "
        "Install PyMuPDF, pypdf or poppler/pdftotext."
    )


def snippet_around_match(
    text: str,
    start: int,
    end: int,
    radius: int = 260,
) -> str:
    left = max(
        0,
        start - radius,
    )

    right = min(
        len(text),
        end + radius,
    )

    return normalize_space(
        text[left:right]
    )


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        PDF,
        RETRIEVAL_LEDGER,
        RETRIEVAL_SUMMARY,
    ):
        require_file(
            required
        )

    retrieval_summary = read_one(
        RETRIEVAL_SUMMARY
    )

    if retrieval_summary.get(
        "decision"
    ) != EXPECTED_RETRIEVAL_DECISION:
        raise RuntimeError(
            "Gate 3R.4 retrieval audit is not accepted."
        )

    retrieval_rows = read_rows(
        RETRIEVAL_LEDGER
    )

    artifact_matches = [
        row
        for row in retrieval_rows
        if row.get(
            "artifact_id"
        )
        == EXPECTED_ARTIFACT_ID
    ]

    if len(
        artifact_matches
    ) != 1:
        raise RuntimeError(
            "Expected exactly one Lele article retrieval row; "
            f"found {len(artifact_matches)}."
        )

    artifact = artifact_matches[0]

    if artifact.get(
        "status"
    ) != "RETRIEVED":
        raise RuntimeError(
            "Lele article was not marked RETRIEVED."
        )

    if artifact.get(
        "actual_kind"
    ) != "PDF":
        raise RuntimeError(
            "Lele article artifact is not a validated PDF."
        )

    ledger_hash = artifact.get(
        "sha256",
        "",
    )

    actual_hash = sha256(
        PDF
    )

    if (
        ledger_hash
        and ledger_hash
        != actual_hash
    ):
        raise RuntimeError(
            "PDF SHA-256 does not match the retrieval ledger."
        )

    pages, extraction_method = extract_pages(
        PDF
    )

    page_blocks = []

    for page_number, page_text in enumerate(
        pages,
        start=1,
    ):
        page_blocks.append(
            (
                f"\n===== PAGE {page_number:04d} =====\n"
                f"{page_text.rstrip()}\n"
            )
        )

    PAGE_TEXT.write_text(
        "".join(
            page_blocks
        ),
        encoding="utf-8",
    )

    evidence_rows = []
    seen_evidence = set()
    category_counts = Counter()

    compiled_patterns = {
        category: [
            re.compile(
                pattern,
                re.IGNORECASE,
            )
            for pattern in patterns
        ]
        for category, patterns in PATTERNS.items()
    }

    for page_number, page_text in enumerate(
        pages,
        start=1,
    ):
        normalized_page = normalize_space(
            page_text
        )

        if not normalized_page:
            continue

        for category, patterns in compiled_patterns.items():
            for pattern in patterns:
                for match in pattern.finditer(
                    normalized_page
                ):
                    snippet = snippet_around_match(
                        normalized_page,
                        match.start(),
                        match.end(),
                    )

                    deduplication_key = (
                        category,
                        page_number,
                        snippet,
                    )

                    if deduplication_key in seen_evidence:
                        continue

                    seen_evidence.add(
                        deduplication_key
                    )

                    category_counts[
                        category
                    ] += 1

                    evidence_rows.append(
                        {
                            "evidence_id": (
                                f"EVID_{len(evidence_rows) + 1:05d}"
                            ),
                            "category": category,
                            "page": page_number,
                            "matched_text": match.group(
                                0
                            ),
                            "snippet": snippet,
                            "interpretation_status": (
                                "TEXTUAL_EVIDENCE_ONLY"
                            ),
                            "R2_coverage_established": False,
                        }
                    )

    write_rows(
        EVIDENCE,
        evidence_rows,
    )

    full_text = "\n".join(
        pages
    )

    DOI_present = (
        re.search(
            re.escape(
                EXPECTED_DOI
            ),
            full_text,
            re.IGNORECASE,
        )
        is not None
    )

    reaxff_present = (
        re.search(
            r"\bReaxFF\b",
            full_text,
            re.IGNORECASE,
        )
        is not None
    )

    supporting_information_present = (
        re.search(
            r"Supporting Information",
            full_text,
            re.IGNORECASE,
        )
        is not None
    )

    parameter_language_present = (
        re.search(
            r"(parameter file|force field parameters|"
            r"parameter set|parameters are available)",
            full_text,
            re.IGNORECASE,
        )
        is not None
    )

    water_hits = category_counts[
        "WATER_OR_AQUEOUS"
    ]

    BNNT_hits = category_counts[
        "BNNT_OR_CURVATURE"
    ]

    edge_hits = category_counts[
        "EDGE_OR_DEFECT"
    ]

    domain_rows = [
        {
            "R2_domain": "B_N_REACTIVE_CHEMISTRY",
            "article_evidence": (
                "PRESENT"
                if (
                    reaxff_present
                    and category_counts[
                        "ELEMENT_DOMAIN"
                    ]
                    > 0
                )
                else "NOT_CONFIRMED"
            ),
            "direct_R2_transferability": "NOT_ESTABLISHED",
            "parameter_file_required": True,
            "additional_validation_required": True,
            "assessment": (
                "Article establishes a B/N/H ReaxFF development "
                "context, not direct R2 coverage."
            ),
        },
        {
            "R2_domain": "B_H_AND_N_H_CHEMISTRY",
            "article_evidence": (
                "TEXTUAL_CHEMICAL_DOMAIN_ONLY"
                if category_counts[
                    "ELEMENT_DOMAIN"
                ]
                > 0
                else "NOT_CONFIRMED"
            ),
            "direct_R2_transferability": "NOT_ESTABLISHED",
            "parameter_file_required": True,
            "additional_validation_required": True,
            "assessment": (
                "Presence of B/N/H chemistry does not prove accurate "
                "edge-passivant parameters for R2."
            ),
        },
        {
            "R2_domain": "FOUR_ATOM_BN_BRIDGE",
            "article_evidence": "NO_DIRECT_VALIDATION_IDENTIFIED",
            "direct_R2_transferability": "NOT_ESTABLISHED",
            "parameter_file_required": True,
            "additional_validation_required": True,
            "assessment": (
                "Dedicated bridge and attachment validation remains required."
            ),
        },
        {
            "R2_domain": "CURVED_BNNT",
            "article_evidence": (
                "MENTION_FOUND"
                if BNNT_hits > 0
                else "NO_DIRECT_EVIDENCE_FOUND"
            ),
            "direct_R2_transferability": "NOT_ESTABLISHED",
            "parameter_file_required": True,
            "additional_validation_required": True,
            "assessment": (
                "Gas-phase synthesis validation is not equivalent "
                "to equilibrium mechanics of the selected R2 BNNT."
            ),
        },
        {
            "R2_domain": "WATER_AND_AQUEOUS_INTERFACE",
            "article_evidence": (
                "MENTION_FOUND_NOT_VALIDATION"
                if water_hits > 0
                else "NO_AQUEOUS_VALIDATION_IDENTIFIED"
            ),
            "direct_R2_transferability": "NOT_ESTABLISHED",
            "parameter_file_required": True,
            "additional_validation_required": True,
            "assessment": (
                "No water compatibility can be inferred from a "
                "gas-phase ReaxFF training domain."
            ),
        },
        {
            "R2_domain": "EDGE_AND_TERMINATION_ENVIRONMENTS",
            "article_evidence": (
                "MENTION_FOUND_NOT_R2_VALIDATION"
                if edge_hits > 0
                else "NO_DIRECT_R2_EDGE_VALIDATION_IDENTIFIED"
            ),
            "direct_R2_transferability": "NOT_ESTABLISHED",
            "parameter_file_required": True,
            "additional_validation_required": True,
            "assessment": (
                "R2 annulus, seed and passivated edge environments "
                "require explicit validation."
            ),
        },
        {
            "R2_domain": "CHARGE_MODEL",
            "article_evidence": (
                "REAXFF_CHARGE_MODEL_CONTEXT"
                if (
                    category_counts[
                        "REAXFF_METHOD"
                    ]
                    > 0
                )
                else "NOT_CONFIRMED"
            ),
            "direct_R2_transferability": "NOT_ESTABLISHED",
            "parameter_file_required": True,
            "additional_validation_required": True,
            "assessment": (
                "ReaxFF charge equilibration must not be treated as "
                "an independently validated aqueous polarization model."
            ),
        },
    ]

    write_rows(
        DOMAIN_MATRIX,
        domain_rows,
    )

    required_artifact_rows = [
        {
            "artifact_id": "LELE_2022_SUPPORTING_INFORMATION_PDF",
            "expected_content": (
                "Training-set details, parameterization procedure, "
                "validation data and supplementary figures/tables"
            ),
            "retrieval_status": "BLOCKED_OR_NOT_RETRIEVED",
            "required_for_parameter_assessment": True,
            "parameter_adoption_authorized_without_artifact": False,
        },
        {
            "artifact_id": "LELE_2022_REAXFF_PARAMETER_FILE",
            "expected_content": (
                "Complete ReaxFF B/N/H parameter set, element order, "
                "functional terms and numerical values"
            ),
            "retrieval_status": "BLOCKED_OR_NOT_RETRIEVED",
            "required_for_parameter_assessment": True,
            "parameter_adoption_authorized_without_artifact": False,
        },
        {
            "artifact_id": "LELE_2022_TRAINING_AND_VALIDATION_DATA",
            "expected_content": (
                "Reference structures, energies, forces, reactions "
                "and error metrics used in fitting and validation"
            ),
            "retrieval_status": "PARTIAL_ARTICLE_DESCRIPTION_ONLY",
            "required_for_parameter_assessment": True,
            "parameter_adoption_authorized_without_artifact": False,
        },
    ]

    write_rows(
        REQUIRED_ARTIFACTS,
        required_artifact_rows,
    )

    gates = {
        "retrieval_gate_is_accepted": (
            retrieval_summary.get(
                "decision"
            )
            == EXPECTED_RETRIEVAL_DECISION
        ),
        "artifact_is_validated_PDF": (
            artifact.get(
                "status"
            )
            == "RETRIEVED"
            and artifact.get(
                "actual_kind"
            )
            == "PDF"
        ),
        "PDF_SHA256_matches_retrieval_ledger": (
            not ledger_hash
            or ledger_hash
            == actual_hash
        ),
        "article_text_was_extracted": (
            len(
                full_text.strip()
            )
            > 1000
        ),
        "DOI_is_present_in_extracted_text": (
            DOI_present
        ),
        "ReaxFF_is_present_in_extracted_text": (
            reaxff_present
        ),
        "all_evidence_is_marked_textual_only": all(
            row[
                "interpretation_status"
            ]
            == "TEXTUAL_EVIDENCE_ONLY"
            and row[
                "R2_coverage_established"
            ]
            is False
            for row in evidence_rows
        ),
        "no_domain_is_declared_directly_transferable_to_R2": all(
            row[
                "direct_R2_transferability"
            ]
            == "NOT_ESTABLISHED"
            for row in domain_rows
        ),
        "SI_and_parameter_file_remain_required": all(
            row[
                "required_for_parameter_assessment"
            ]
            is True
            for row in required_artifact_rows
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
        PASS_DECISION
        if accepted
        else REVIEW_DECISION
    )

    required_next_step = (
        "RETRIEVE_AND_AUDIT_LELE_2022_SI_AND_REAXFF_PARAMETER_FILE"
        if accepted
        else
        "REVIEW_LELE_2022_PRIMARY_ARTICLE_CONTENT_AUDIT"
    )

    summary = {
        "decision": decision,
        "DOI": EXPECTED_DOI,
        "PDF_file": relative(
            PDF
        ),
        "PDF_sha256": actual_hash,
        "PDF_pages": len(
            pages
        ),
        "text_extraction_method": extraction_method,
        "extracted_characters": len(
            full_text
        ),
        "evidence_rows": len(
            evidence_rows
        ),
        "evidence_categories_with_hits": sum(
            count > 0
            for count in category_counts.values()
        ),
        "DOI_present": DOI_present,
        "ReaxFF_present": reaxff_present,
        "supporting_information_mentioned": (
            supporting_information_present
        ),
        "parameter_availability_language_present": (
            parameter_language_present
        ),
        "water_or_aqueous_hits": water_hits,
        "BNNT_or_curvature_hits": BNNT_hits,
        "edge_or_defect_hits": edge_hits,
        "R2_force_field_coverage_established": False,
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

    write_rows(
        GATES,
        [
            {
                "gate": name,
                "pass": passed,
            }
            for name, passed in gates.items()
        ],
    )

    JSON_OUT.write_text(
        json.dumps(
            {
                "summary": summary,
                "category_counts": dict(
                    category_counts
                ),
                "domain_assessment": domain_rows,
                "gates": gates,
                "limitations": [
                    (
                        "Textual article evidence does not establish "
                        "numerical parameter coverage."
                    ),
                    (
                        "The Supporting Information and full ReaxFF "
                        "parameter file remain unavailable locally."
                    ),
                    (
                        "The gas-phase training domain does not establish "
                        "aqueous, curved-BNNT or R2 bridge transferability."
                    ),
                    (
                        "No topology, charges, force-field assignment, "
                        "minimization, MD or QM calculation is generated."
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
                "role": "Lele_2022_primary_article_PDF",
                "file": relative(
                    PDF
                ),
                "sha256": actual_hash,
            },
            {
                "role": "Gate3R4_retrieval_ledger",
                "file": relative(
                    RETRIEVAL_LEDGER
                ),
                "sha256": sha256(
                    RETRIEVAL_LEDGER
                ),
            },
            {
                "role": "Gate3R4_retrieval_summary",
                "file": relative(
                    RETRIEVAL_SUMMARY
                ),
                "sha256": sha256(
                    RETRIEVAL_SUMMARY
                ),
            },
        ],
    )

    category_lines = "\n".join(
        (
            f"- {category}: **{category_counts[category]}** evidence rows"
        )
        for category in PATTERNS
    )

    domain_lines = "\n".join(
        (
            f"- **{row['R2_domain']}**: "
            f"{row['article_evidence']}; "
            f"R2 transferability={row['direct_R2_transferability']}."
        )
        for row in domain_rows
    )

    REPORT.write_text(
        f"""# Lele 2022 Primary Article Content Audit

## Artifact

- DOI: `{EXPECTED_DOI}`
- PDF: `{relative(PDF)}`
- SHA-256: `{actual_hash}`
- Pages: **{len(pages)}**
- Extraction method: **{extraction_method}**
- Extracted characters: **{len(full_text)}**

## Evidence inventory

{category_lines}

All evidence rows are textual evidence only. No numerical force-field
coverage is inferred from keyword occurrence.

## R2 domain assessment

{domain_lines}

## Artifact status

- Supporting Information mentioned:
  **{supporting_information_present}**
- Parameter-availability language found:
  **{parameter_language_present}**
- Supporting Information retrieved:
  **NO**
- ReaxFF parameter file retrieved:
  **NO**

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- R2 force-field coverage established:
  **NO**
- Parameter adoption authorized:
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
- Required next step:
  `{required_next_step}`
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 Lele 2022 primary-article "
        "content audit completed."
    )

    print(
        "PDF pages / extraction method / characters: "
        f"{len(pages)}/"
        f"{extraction_method}/"
        f"{len(full_text)}"
    )

    print(
        "DOI / ReaxFF / Supporting Information mentioned: "
        f"{DOI_present}/"
        f"{reaxff_present}/"
        f"{supporting_information_present}"
    )

    print(
        "Evidence rows / categories with hits: "
        f"{len(evidence_rows)}/"
        f"{sum(count > 0 for count in category_counts.values())}"
    )

    print(
        "Training / parameter / water / BNNT / edge evidence counts: "
        f"{category_counts['TRAINING_SET']}/"
        f"{category_counts['PARAMETER_AVAILABILITY']}/"
        f"{category_counts['WATER_OR_AQUEOUS']}/"
        f"{category_counts['BNNT_OR_CURVATURE']}/"
        f"{category_counts['EDGE_OR_DEFECT']}"
    )

    print(
        "Supporting Information retrieved: NO"
    )

    print(
        "ReaxFF parameter file retrieved: NO"
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
        "R2 force-field coverage established: NO"
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
        PAGE_TEXT,
        EVIDENCE,
        DOMAIN_MATRIX,
        REQUIRED_ARTIFACTS,
        SUMMARY,
        GATES,
        JSON_OUT,
        MANIFEST,
        REPORT,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not accepted:
        raise RuntimeError(
            "Lele 2022 primary-article content audit requires review."
        )


if __name__ == "__main__":
    main()
