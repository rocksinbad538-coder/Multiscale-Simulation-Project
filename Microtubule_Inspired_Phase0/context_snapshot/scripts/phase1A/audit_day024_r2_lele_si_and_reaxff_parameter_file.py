#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

import fitz


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/phase1A/day024_chemical_end_rim_design"

PREVIOUS_GATE = (
    BASE
    / "36_r2_lele_primary_article_content_audit"
)

PREVIOUS_SUMMARY = (
    PREVIOUS_GATE
    / "lele_2022_primary_article_content_audit_summary.csv"
)

OUT = (
    BASE
    / "37_r2_lele_si_parameter_file_audit"
)

RAW = OUT / "raw"

ARTICLE_PDF = RAW / "jp1c09648.pdf"
SI_PDF = RAW / "jp1c09648_si_001.pdf"
PARAMETER_FILE = RAW / "jp1c09648_si_002.txt"

SI_TEXT = (
    OUT
    / "lele_2022_supporting_information_text_by_page.txt"
)

SI_EVIDENCE = (
    OUT
    / "lele_2022_supporting_information_evidence.csv"
)

PARAMETER_SECTION_INVENTORY = (
    OUT
    / "lele_2022_reaxff_section_inventory.csv"
)

ELEMENT_INVENTORY = (
    OUT
    / "lele_2022_reaxff_element_inventory.csv"
)

BNH_TERM_INVENTORY = (
    OUT
    / "lele_2022_reaxff_bnh_term_inventory.csv"
)

DOMAIN_ASSESSMENT = (
    OUT
    / "lele_2022_R2_domain_assessment.csv"
)

SUMMARY = (
    OUT
    / "lele_2022_si_parameter_file_audit_summary.csv"
)

GATES = (
    OUT
    / "lele_2022_si_parameter_file_audit_gates.csv"
)

MANIFEST = (
    OUT
    / "lele_2022_si_parameter_file_audit_manifest.csv"
)

JSON_OUT = (
    OUT
    / "lele_2022_si_parameter_file_audit.json"
)

REPORT = (
    OUT
    / "LELE_2022_SI_PARAMETER_FILE_AUDIT_DAY024.md"
)

EXPECTED_PREVIOUS_DECISION = (
    "R2_LELE_PRIMARY_ARTICLE_CONTENT_AUDITED_"
    "SI_AND_PARAMETER_FILES_STILL_REQUIRED"
)

PASS_DECISION = (
    "R2_LELE_2022_SI_AND_REAXFF_PARAMETER_FILE_"
    "AUDITED_NOT_AUTHORIZED_FOR_R2"
)

REVIEW_DECISION = (
    "R2_LELE_2022_SI_OR_REAXFF_PARAMETER_FILE_REQUIRES_REVIEW"
)

EXPECTED_HASHES = {
    ARTICLE_PDF: (
        8390448,
        "3f997ca60f0242bfc85cd54d8d6bc5a9453ec2b2180fee603e6aaeb6defbeb45",
    ),
    SI_PDF: (
        199670,
        "6fa080473219004e8ed753788943243ac6ffa27597063d24716b7eb36195d040",
    ),
    PARAMETER_FILE: (
        21491,
        "5c5f0826e0dd7b1b7bc6d843437b1325f3d4fa6931046c015ab1ecbd904f47ae",
    ),
}

SECTION_PATTERNS = [
    (
        "GENERAL_PARAMETERS",
        re.compile(
            r"^\s*(\d+)\s+!\s*Number of general parameters",
            re.IGNORECASE,
        ),
    ),
    (
        "ATOMS",
        re.compile(
            r"^\s*(\d+)\s+!\s*Nr of atoms",
            re.IGNORECASE,
        ),
    ),
    (
        "BONDS",
        re.compile(
            r"^\s*(\d+)\s+!\s*Nr of bonds",
            re.IGNORECASE,
        ),
    ),
    (
        "OFF_DIAGONAL",
        re.compile(
            r"^\s*(\d+)\s+!\s*Nr of off-diagonal terms",
            re.IGNORECASE,
        ),
    ),
    (
        "ANGLES",
        re.compile(
            r"^\s*(\d+)\s+!\s*Nr of angles",
            re.IGNORECASE,
        ),
    ),
    (
        "TORSIONS",
        re.compile(
            r"^\s*(\d+)\s+!\s*Nr of torsions",
            re.IGNORECASE,
        ),
    ),
    (
        "HYDROGEN_BONDS",
        re.compile(
            r"^\s*(\d+)\s+!\s*Nr of hydrogen bonds",
            re.IGNORECASE,
        ),
    ),
]

EXPECTED_SECTION_COUNTS = {
    "GENERAL_PARAMETERS": 39,
    "ATOMS": 6,
    "BONDS": 21,
    "OFF_DIAGONAL": 15,
    "ANGLES": 105,
    "TORSIONS": 58,
    "HYDROGEN_BONDS": 4,
}

ELEMENT_ORDER_EXPECTED = [
    "C",
    "H",
    "O",
    "N",
    "B",
    "Al",
]

ELEMENT_TO_INDEX = {
    element: index
    for index, element in enumerate(
        ELEMENT_ORDER_EXPECTED,
        start=1,
    )
}

SI_PATTERNS = {
    "QUALITY_FACTOR": [
        r"quality factor",
        r"six-membered rings",
        r"maximum number of hexagons",
    ],
    "BORAZINE_SNAPSHOT": [
        r"borazine",
        r"2000 K",
        r"4 ns",
    ],
    "PARAMETER_TABLE": [
        r"force field parameter",
        r"bond parameter",
        r"angle parameter",
        r"torsion parameter",
    ],
    "WATER": [
        r"\bwater\b",
        r"aqueous",
        r"hydration",
        r"solvation",
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


def detect_pdf(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(5) == b"%PDF-"


def parameter_file_is_html(text: str) -> bool:
    head = text.lstrip()[:1000].lower()

    return (
        "<!doctype html" in head
        or "<html" in head
    )


def extract_pdf_pages(path: Path) -> list[str]:
    document = fitz.open(path)

    try:
        return [
            page.get_text("text")
            for page in document
        ]
    finally:
        document.close()


def integer_tokens(line: str) -> list[int]:
    tokens = []

    for token in line.split():
        try:
            tokens.append(
                int(token)
            )
        except ValueError:
            continue

    return tokens


def finite_numeric_tokens(line: str) -> list[float]:
    values = []

    for token in line.split():
        cleaned = token.strip(
            ",;()[]{}"
        )

        try:
            value = float(cleaned)
        except ValueError:
            continue

        if math.isfinite(value):
            values.append(value)

    return values


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        PREVIOUS_SUMMARY,
        ARTICLE_PDF,
        SI_PDF,
        PARAMETER_FILE,
    ):
        require_file(
            required
        )

    previous_summary = read_one(
        PREVIOUS_SUMMARY
    )

    if previous_summary.get(
        "decision"
    ) != EXPECTED_PREVIOUS_DECISION:
        raise RuntimeError(
            "Gate 3R.5 article audit is not accepted."
        )

    integrity_results = {}

    for path, (
        expected_size,
        expected_hash,
    ) in EXPECTED_HASHES.items():
        actual_size = path.stat().st_size
        actual_hash = sha256(
            path
        )

        integrity_results[
            relative(path)
        ] = {
            "expected_size": expected_size,
            "actual_size": actual_size,
            "size_match": actual_size == expected_size,
            "expected_sha256": expected_hash,
            "actual_sha256": actual_hash,
            "hash_match": actual_hash == expected_hash,
        }

    article_pdf_valid = detect_pdf(
        ARTICLE_PDF
    )

    SI_pdf_valid = detect_pdf(
        SI_PDF
    )

    parameter_text = PARAMETER_FILE.read_text(
        encoding="utf-8",
        errors="replace",
    )

    parameter_text_valid = (
        len(
            parameter_text.strip()
        )
        > 1000
        and not parameter_file_is_html(
            parameter_text
        )
    )

    article_pages = extract_pdf_pages(
        ARTICLE_PDF
    )

    SI_pages = extract_pdf_pages(
        SI_PDF
    )

    article_text = "\n".join(
        article_pages
    )

    SI_full_text = "\n".join(
        SI_pages
    )

    SI_TEXT.write_text(
        "".join(
            (
                f"\n===== PAGE {index:04d} =====\n"
                f"{page.rstrip()}\n"
            )
            for index, page in enumerate(
                SI_pages,
                start=1,
            )
        ),
        encoding="utf-8",
    )

    SI_evidence_rows = []
    SI_category_counts = {
        category: 0
        for category in SI_PATTERNS
    }

    for page_number, page_text in enumerate(
        SI_pages,
        start=1,
    ):
        normalized = normalize_space(
            page_text
        )

        for category, patterns in SI_PATTERNS.items():
            for raw_pattern in patterns:
                pattern = re.compile(
                    raw_pattern,
                    re.IGNORECASE,
                )

                for match in pattern.finditer(
                    normalized
                ):
                    SI_category_counts[
                        category
                    ] += 1

                    left = max(
                        0,
                        match.start() - 180,
                    )

                    right = min(
                        len(normalized),
                        match.end() + 180,
                    )

                    SI_evidence_rows.append(
                        {
                            "evidence_id": (
                                f"SI_EVID_{len(SI_evidence_rows) + 1:04d}"
                            ),
                            "category": category,
                            "page": page_number,
                            "matched_text": match.group(
                                0
                            ),
                            "snippet": normalize_space(
                                normalized[left:right]
                            ),
                            "R2_parameter_coverage_established": False,
                        }
                    )

    write_rows(
        SI_EVIDENCE,
        SI_evidence_rows,
    )

    lines = parameter_text.splitlines()

    section_counts_found: dict[str, int] = {}
    section_line_numbers: dict[str, int] = {}

    for line_number, line in enumerate(
        lines,
        start=1,
    ):
        for section_name, pattern in SECTION_PATTERNS:
            match = pattern.search(
                line
            )

            if match:
                section_counts_found[
                    section_name
                ] = int(
                    match.group(1)
                )

                section_line_numbers[
                    section_name
                ] = line_number

    section_rows = []

    for section_name in EXPECTED_SECTION_COUNTS:
        expected_count = EXPECTED_SECTION_COUNTS[
            section_name
        ]

        actual_count = section_counts_found.get(
            section_name
        )

        section_rows.append(
            {
                "section": section_name,
                "expected_count": expected_count,
                "actual_count": (
                    actual_count
                    if actual_count is not None
                    else ""
                ),
                "count_match": (
                    actual_count == expected_count
                ),
                "header_line": section_line_numbers.get(
                    section_name,
                    "",
                ),
            }
        )

    write_rows(
        PARAMETER_SECTION_INVENTORY,
        section_rows,
    )

    atom_header_line = section_line_numbers.get(
        "ATOMS"
    )

    found_elements: list[str] = []

    if atom_header_line is not None:
        for line in lines[
            atom_header_line:
        ]:
            match = re.match(
                r"^\s*([A-Z][a-z]?)\s+[-+0-9.]",
                line,
            )

            if not match:
                continue

            element = match.group(
                1
            )

            if element not in found_elements:
                found_elements.append(
                    element
                )

            if len(
                found_elements
            ) == section_counts_found.get(
                "ATOMS",
                0,
            ):
                break

    element_rows = []

    for expected_index, element in enumerate(
        ELEMENT_ORDER_EXPECTED,
        start=1,
    ):
        element_rows.append(
            {
                "element": element,
                "expected_index": expected_index,
                "present": element in found_elements,
                "observed_index": (
                    found_elements.index(
                        element
                    )
                    + 1
                    if element in found_elements
                    else ""
                ),
                "order_match": (
                    len(found_elements)
                    >= expected_index
                    and found_elements[
                        expected_index - 1
                    ]
                    == element
                ),
            }
        )

    write_rows(
        ELEMENT_INVENTORY,
        element_rows,
    )

    B_index = ELEMENT_TO_INDEX[
        "B"
    ]

    N_index = ELEMENT_TO_INDEX[
        "N"
    ]

    H_index = ELEMENT_TO_INDEX[
        "H"
    ]

    section_boundaries = []

    for section_name, line_number in section_line_numbers.items():
        section_boundaries.append(
            (
                line_number,
                section_name,
            )
        )

    section_boundaries.sort()

    def section_for_line(
        line_number: int,
    ) -> str:
        current = "PREAMBLE"

        for start_line, section_name in section_boundaries:
            if line_number >= start_line:
                current = section_name
            else:
                break

        return current

    BNH_term_rows = []

    for line_number, line in enumerate(
        lines,
        start=1,
    ):
        section = section_for_line(
            line_number
        )

        stripped = line.strip()

        if not stripped:
            continue

        integers = integer_tokens(
            stripped
        )

        if section == "BONDS" and len(integers) >= 2:
            atom_indices = integers[:2]

        elif section == "OFF_DIAGONAL" and len(integers) >= 2:
            atom_indices = integers[:2]

        elif section == "ANGLES" and len(integers) >= 3:
            atom_indices = integers[:3]

        elif section == "TORSIONS" and len(integers) >= 4:
            atom_indices = integers[:4]

        elif section == "HYDROGEN_BONDS" and len(integers) >= 3:
            atom_indices = integers[:3]

        else:
            continue

        relevant_indices = {
            H_index,
            N_index,
            B_index,
        }

        if not any(
            index in relevant_indices
            for index in atom_indices
        ):
            continue

        mapped_elements = [
            ELEMENT_ORDER_EXPECTED[
                index - 1
            ]
            if 1 <= index <= len(
                ELEMENT_ORDER_EXPECTED
            )
            else (
                "WILDCARD"
                if index == 0
                else f"UNKNOWN_{index}"
            )
            for index in atom_indices
        ]

        numeric_values = finite_numeric_tokens(
            stripped
        )

        BNH_term_rows.append(
            {
                "line_number": line_number,
                "section": section,
                "atom_indices": " ".join(
                    str(index)
                    for index in atom_indices
                ),
                "elements": "-".join(
                    mapped_elements
                ),
                "contains_H": (
                    H_index in atom_indices
                ),
                "contains_N": (
                    N_index in atom_indices
                ),
                "contains_B": (
                    B_index in atom_indices
                ),
                "contains_B_and_N": (
                    B_index in atom_indices
                    and N_index in atom_indices
                ),
                "contains_B_and_H": (
                    B_index in atom_indices
                    and H_index in atom_indices
                ),
                "contains_N_and_H": (
                    N_index in atom_indices
                    and H_index in atom_indices
                ),
                "numeric_token_count": len(
                    numeric_values
                ),
                "text": stripped[:1000],
                "R2_transferability_established": False,
            }
        )

    write_rows(
        BNH_TERM_INVENTORY,
        BNH_term_rows,
    )

    BNH_section_counts: dict[str, int] = {}

    for row in BNH_term_rows:
        section = row[
            "section"
        ]

        BNH_section_counts[
            section
        ] = (
            BNH_section_counts.get(
                section,
                0,
            )
            + 1
        )

    direct_article_findings = {
        "DOI_present": (
            "10.1021/acs.jpca.1c09648"
            in article_text
        ),
        "gas_phase_scope_present": (
            re.search(
                r"gas[- ]phase",
                article_text,
                re.IGNORECASE,
            )
            is not None
        ),
        "B3LYP_method_present": (
            re.search(
                r"\bB3LYP\b",
                article_text,
                re.IGNORECASE,
            )
            is not None
        ),
        "basis_6311_present": (
            re.search(
                r"6[-–]311G",
                article_text,
                re.IGNORECASE,
            )
            is not None
        ),
        "water_validation_present": (
            re.search(
                r"(aqueous|liquid water|water interface|hydration)",
                article_text,
                re.IGNORECASE,
            )
            is not None
        ),
        "R2_bridge_validation_present": False,
    }

    domain_rows = [
        {
            "domain": "B_N_H_REACTIVE_ELEMENTAL_DOMAIN",
            "evidence": (
                "CONFIRMED"
                if {
                    "B",
                    "N",
                    "H",
                }.issubset(
                    set(
                        found_elements
                    )
                )
                else "NOT_CONFIRMED"
            ),
            "R2_coverage": "NOT_ESTABLISHED",
            "interpretation": (
                "The parameter file contains B, N and H, but elemental "
                "presence does not establish transferability to R2."
            ),
        },
        {
            "domain": "B_N_BONDED_REACTIVE_CHEMISTRY",
            "evidence": (
                "PARAMETER_RECORDS_PRESENT"
                if any(
                    row[
                        "contains_B_and_N"
                    ]
                    for row in BNH_term_rows
                )
                else "NOT_FOUND"
            ),
            "R2_coverage": "NOT_ESTABLISHED",
            "interpretation": (
                "The model was trained for gas-phase reactions and "
                "nanostructure synthesis, not equilibrium R2 mechanics."
            ),
        },
        {
            "domain": "B_H_REACTIVE_CHEMISTRY",
            "evidence": (
                "PARAMETER_RECORDS_PRESENT"
                if any(
                    row[
                        "contains_B_and_H"
                    ]
                    for row in BNH_term_rows
                )
                else "NOT_FOUND"
            ),
            "R2_coverage": "NOT_ESTABLISHED",
            "interpretation": (
                "B-H chemistry is represented, but R2 annulus and "
                "bridge-edge environments were not directly validated."
            ),
        },
        {
            "domain": "N_H_REACTIVE_CHEMISTRY",
            "evidence": (
                "PARAMETER_RECORDS_PRESENT"
                if any(
                    row[
                        "contains_N_and_H"
                    ]
                    for row in BNH_term_rows
                )
                else "NOT_FOUND"
            ),
            "R2_coverage": "NOT_ESTABLISHED",
            "interpretation": (
                "N-H chemistry is represented, but R2 termination "
                "transferability remains unverified."
            ),
        },
        {
            "domain": "CURVED_BNNT_EQUILIBRIUM_MECHANICS",
            "evidence": "NOT_DIRECTLY_VALIDATED",
            "R2_coverage": "NOT_ESTABLISHED",
            "interpretation": (
                "Formation of curved or closed BN structures does not "
                "constitute validation of equilibrium mechanics at the R2 radius."
            ),
        },
        {
            "domain": "FOUR_ATOM_B_N_B_N_BRIDGE",
            "evidence": "NO_DIRECT_REFERENCE_CONFIGURATION",
            "R2_coverage": "NOT_ESTABLISHED",
            "interpretation": (
                "The selected bridge and attachment environments remain "
                "outside the demonstrated validation set."
            ),
        },
        {
            "domain": "WATER_AND_AQUEOUS_INTERFACE",
            "evidence": (
                "NO_VALIDATION_IDENTIFIED"
                if not direct_article_findings[
                    "water_validation_present"
                ]
                and SI_category_counts[
                    "WATER"
                ]
                == 0
                else "MENTION_ONLY"
            ),
            "R2_coverage": "NOT_ESTABLISHED",
            "interpretation": (
                "The force field was developed for gas-phase synthesis "
                "and cannot be assumed compatible with confined water."
            ),
        },
        {
            "domain": "ANISOTROPIC_POLARIZATION",
            "evidence": "EEM_CHARGE_EQUILIBRATION_ONLY",
            "R2_coverage": "NOT_ESTABLISHED",
            "interpretation": (
                "ReaxFF EEM charges are not equivalent to an independently "
                "validated polarizable h-BN/water model."
            ),
        },
    ]

    write_rows(
        DOMAIN_ASSESSMENT,
        domain_rows,
    )

    gates = {
        "previous_article_audit_is_accepted": (
            previous_summary.get(
                "decision"
            )
            == EXPECTED_PREVIOUS_DECISION
        ),
        "all_artifact_hashes_match": all(
            row[
                "hash_match"
            ]
            for row in integrity_results.values()
        ),
        "all_artifact_sizes_match": all(
            row[
                "size_match"
            ]
            for row in integrity_results.values()
        ),
        "article_is_valid_PDF": (
            article_pdf_valid
        ),
        "SI_is_valid_PDF": (
            SI_pdf_valid
        ),
        "SI_has_two_pages": (
            len(
                SI_pages
            )
            == 2
        ),
        "parameter_file_is_valid_text": (
            parameter_text_valid
        ),
        "all_ReaxFF_section_counts_match": all(
            row[
                "count_match"
            ]
            for row in section_rows
        ),
        "element_order_matches_C_H_O_N_B_Al": (
            found_elements
            == ELEMENT_ORDER_EXPECTED
        ),
        "BNH_terms_are_present": (
            any(
                row[
                    "contains_B_and_N"
                ]
                for row in BNH_term_rows
            )
            and any(
                row[
                    "contains_B_and_H"
                ]
                for row in BNH_term_rows
            )
            and any(
                row[
                    "contains_N_and_H"
                ]
                for row in BNH_term_rows
            )
        ),
        "all_R2_domains_remain_unvalidated": all(
            row[
                "R2_coverage"
            ]
            == "NOT_ESTABLISHED"
            for row in domain_rows
        ),
        "no_parameter_adoption_topology_MD_or_QM_generated": True,
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
        "COMPARE_LELE_REAXFF_DOMAIN_WITH_FUNCTIONALIZED_HBN_"
        "FIXED_TOPOLOGY_AND_QM_REFERENCE_ROUTES"
        if accepted
        else
        "REVIEW_LELE_2022_SI_OR_PARAMETER_FILE_AUDIT"
    )

    summary = {
        "decision": decision,
        "article_pages": len(
            article_pages
        ),
        "SI_pages": len(
            SI_pages
        ),
        "SI_extracted_characters": len(
            SI_full_text
        ),
        "SI_evidence_rows": len(
            SI_evidence_rows
        ),
        "SI_parameter_table_hits": SI_category_counts[
            "PARAMETER_TABLE"
        ],
        "SI_water_hits": SI_category_counts[
            "WATER"
        ],
        "parameter_file_lines": len(
            lines
        ),
        "general_parameters": section_counts_found.get(
            "GENERAL_PARAMETERS",
            "",
        ),
        "atoms": section_counts_found.get(
            "ATOMS",
            "",
        ),
        "bonds": section_counts_found.get(
            "BONDS",
            "",
        ),
        "off_diagonal_terms": section_counts_found.get(
            "OFF_DIAGONAL",
            "",
        ),
        "angles": section_counts_found.get(
            "ANGLES",
            "",
        ),
        "torsions": section_counts_found.get(
            "TORSIONS",
            "",
        ),
        "hydrogen_bonds": section_counts_found.get(
            "HYDROGEN_BONDS",
            "",
        ),
        "element_order": (
            " | ".join(
                found_elements
            )
        ),
        "BNH_relevant_parameter_records": len(
            BNH_term_rows
        ),
        "BNH_records_by_section": (
            " | ".join(
                f"{section}:{count}"
                for section, count in sorted(
                    BNH_section_counts.items()
                )
            )
        ),
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

    manifest_rows = []

    for path, role in (
        (
            ARTICLE_PDF,
            "Lele_2022_primary_article",
        ),
        (
            SI_PDF,
            "Lele_2022_supporting_information",
        ),
        (
            PARAMETER_FILE,
            "Lele_2022_ReaxFF_parameter_file",
        ),
        (
            PREVIOUS_SUMMARY,
            "Gate3R5_article_audit_summary",
        ),
    ):
        manifest_rows.append(
            {
                "role": role,
                "file": relative(
                    path
                ),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(
                    path
                ),
            }
        )

    write_rows(
        MANIFEST,
        manifest_rows,
    )

    JSON_OUT.write_text(
        json.dumps(
            {
                "summary": summary,
                "integrity_results": integrity_results,
                "SI_category_counts": SI_category_counts,
                "section_counts": section_counts_found,
                "elements_found": found_elements,
                "BNH_section_counts": BNH_section_counts,
                "domain_assessment": domain_rows,
                "article_findings": direct_article_findings,
                "gates": gates,
                "limitations": [
                    (
                        "The SI PDF contains quality-factor methodology "
                        "and an additional borazine snapshot, not an "
                        "expanded parameter-validation data set."
                    ),
                    (
                        "The ReaxFF file is a composite C/H/O/N/B/Al set, "
                        "not an R2-specific B/N/H parameterization."
                    ),
                    (
                        "Presence of B-N, B-H and N-H records does not "
                        "establish transferability to the R2 annulus, "
                        "bridge, water or polarization domains."
                    ),
                    (
                        "No parameter file is copied to an accepted "
                        "force-field directory."
                    ),
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    REPORT.write_text(
        f"""# Lele 2022 SI and ReaxFF Parameter-File Audit

## Verified artifacts

- Primary article:
  `{relative(ARTICLE_PDF)}`
- Supporting Information:
  `{relative(SI_PDF)}`
- ReaxFF parameter file:
  `{relative(PARAMETER_FILE)}`

All three artifacts match the expected SHA-256 hashes.

## Supporting Information

- Pages: **{len(SI_pages)}**
- Extracted characters: **{len(SI_full_text)}**
- Quality-factor evidence rows:
  **{SI_category_counts['QUALITY_FACTOR']}**
- Borazine-snapshot evidence rows:
  **{SI_category_counts['BORAZINE_SNAPSHOT']}**
- Parameter-table hits:
  **{SI_category_counts['PARAMETER_TABLE']}**
- Water hits:
  **{SI_category_counts['WATER']}**

The SI contains the quality-factor derivation and the 4 ns
borazine simulation snapshot. It does not supply a separate
R2-relevant validation set.

## ReaxFF parameter-file structure

- General parameters:
  **{section_counts_found.get('GENERAL_PARAMETERS')}**
- Elements:
  **{section_counts_found.get('ATOMS')}**
- Bonds:
  **{section_counts_found.get('BONDS')}**
- Off-diagonal terms:
  **{section_counts_found.get('OFF_DIAGONAL')}**
- Angles:
  **{section_counts_found.get('ANGLES')}**
- Torsions:
  **{section_counts_found.get('TORSIONS')}**
- Hydrogen bonds:
  **{section_counts_found.get('HYDROGEN_BONDS')}**
- Element order:
  **{' | '.join(found_elements)}**
- B/N/H-relevant parameter records:
  **{len(BNH_term_rows)}**

## Scientific conclusion

The file is a genuine composite ReaxFF parameter set containing
C, H, O, N, B and Al. It includes reactive records involving
B-N, B-H and N-H chemistry. Its demonstrated target domain is
gas-phase B/N/H chemistry and high-temperature BN nanostructure
formation.

It does not establish validated transferability to:

- equilibrium mechanics of the selected R2 BNNT;
- reconstructed annulus environments;
- four-atom B-N-B-N bridges;
- confined water;
- anisotropic scaffold-water polarization.

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
- Force-field parameterization authorized:
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
        "Day024 R2 Lele 2022 SI and ReaxFF "
        "parameter-file audit completed."
    )

    print(
        "Article pages / SI pages / SI characters: "
        f"{len(article_pages)}/"
        f"{len(SI_pages)}/"
        f"{len(SI_full_text)}"
    )

    print(
        "SI quality-factor / borazine / parameter-table / water hits: "
        f"{SI_category_counts['QUALITY_FACTOR']}/"
        f"{SI_category_counts['BORAZINE_SNAPSHOT']}/"
        f"{SI_category_counts['PARAMETER_TABLE']}/"
        f"{SI_category_counts['WATER']}"
    )

    print(
        "ReaxFF general/atoms/bonds/offdiag/angles/torsions/hbonds: "
        f"{section_counts_found.get('GENERAL_PARAMETERS')}/"
        f"{section_counts_found.get('ATOMS')}/"
        f"{section_counts_found.get('BONDS')}/"
        f"{section_counts_found.get('OFF_DIAGONAL')}/"
        f"{section_counts_found.get('ANGLES')}/"
        f"{section_counts_found.get('TORSIONS')}/"
        f"{section_counts_found.get('HYDROGEN_BONDS')}"
    )

    print(
        "Element order: "
        + " | ".join(
            found_elements
        )
    )

    print(
        "B/N/H-relevant parameter records: "
        f"{len(BNH_term_rows)}"
    )

    print(
        "B/N/H records by section: "
        + (
            " | ".join(
                f"{section}:{count}"
                for section, count in sorted(
                    BNH_section_counts.items()
                )
            )
            or "NONE"
        )
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
        SI_TEXT,
        SI_EVIDENCE,
        PARAMETER_SECTION_INVENTORY,
        ELEMENT_INVENTORY,
        BNH_TERM_INVENTORY,
        DOMAIN_ASSESSMENT,
        SUMMARY,
        GATES,
        MANIFEST,
        JSON_OUT,
        REPORT,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not accepted:
        raise RuntimeError(
            "Lele 2022 SI or parameter-file audit requires review."
        )


if __name__ == "__main__":
    main()
