#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

DAY024 = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design"
)

DAY025 = (
    ROOT
    / "runs/phase1A/day025_force_field_route_comparison"
)

GATE34 = (
    DAY024
    / "34_r2_primary_source_provenance_correction"
)

GATE35 = (
    DAY024
    / "35_r2_p1_source_artifact_retrieval"
)

GATE39 = (
    DAY025
    / "39_r2_force_field_route_coverage_and_risk_matrix"
)

INPUTS = {
    "gate39_summary": (
        GATE39
        / "r2_force_field_route_coverage_and_risk_summary.csv"
    ),
    "gate39_environment_matrix": (
        GATE39
        / "r2_environment_route_coverage_and_risk_matrix.csv"
    ),
    "gate39_term_matrix": (
        GATE39
        / "r2_bonded_term_route_coverage_and_risk_matrix.csv"
    ),
    "gate34_source_ledger": (
        GATE34
        / "r2_corrected_primary_source_ledger.csv"
    ),
    "gate34_domain_coverage": (
        GATE34
        / "r2_corrected_candidate_model_domain_coverage.csv"
    ),
    "gate34_model_decisions": (
        GATE34
        / "r2_corrected_preliminary_model_decisions.csv"
    ),
    "gate34_prioritized_artifacts": (
        GATE34
        / "r2_prioritized_source_artifact_retrieval.csv"
    ),
    "gate35_retrieval_ledger": (
        GATE35
        / "r2_p1_artifact_retrieval_ledger.csv"
    ),
    "gate35_file_inventory": (
        GATE35
        / "r2_p1_retrieved_file_inventory.csv"
    ),
}

OUT = (
    DAY025
    / "40_r2_functionalized_hbn_fixed_topology_primary_source_audit"
)

SOURCE_RECORDS = (
    OUT
    / "r2_functionalized_hbn_primary_source_records.csv"
)

SOURCE_DOMAIN_MATRIX = (
    OUT
    / "r2_functionalized_hbn_source_domain_matrix.csv"
)

ARTIFACT_GAP_MATRIX = (
    OUT
    / "r2_functionalized_hbn_source_artifact_gap_matrix.csv"
)

R2_SCOPE_MAPPING = (
    OUT
    / "r2_functionalized_hbn_source_to_R2_scope_mapping.csv"
)

ROUTE_ASSESSMENT = (
    OUT
    / "r2_functionalized_hbn_fixed_topology_route_assessment.csv"
)

SUMMARY = (
    OUT
    / "r2_functionalized_hbn_fixed_topology_primary_source_audit_summary.csv"
)

GATES = (
    OUT
    / "r2_functionalized_hbn_fixed_topology_primary_source_audit_gates.csv"
)

MANIFEST = (
    OUT
    / "r2_functionalized_hbn_fixed_topology_primary_source_audit_manifest.csv"
)

JSON_OUT = (
    OUT
    / "r2_functionalized_hbn_fixed_topology_primary_source_audit.json"
)

REPORT = (
    OUT
    / "R2_FUNCTIONALIZED_HBN_FIXED_TOPOLOGY_PRIMARY_SOURCE_AUDIT_DAY025.md"
)

EXPECTED_GATE39_DECISION = (
    "R2_FORCE_FIELD_ROUTE_COVERAGE_AND_RISK_MATRIX_BUILT_"
    "NO_ROUTE_YET_AUTHORIZED"
)

PASS_DECISION = (
    "R2_FUNCTIONALIZED_HBN_FIXED_TOPOLOGY_PRIMARY_SOURCES_"
    "AUDITED_ARTIFACTS_STILL_REQUIRED"
)

REVIEW_DECISION = (
    "R2_FUNCTIONALIZED_HBN_FIXED_TOPOLOGY_PRIMARY_SOURCE_"
    "AUDIT_REQUIRES_REVIEW"
)

SOURCES = [
    {
        "source_id": "GHORAI_2025_FUNCTIONALIZED_HBN_NANOPORE",
        "short_name": "Ghorai 2025",
        "doi": "10.1063/5.0242541",
        "route_role": (
            "FUNCTIONALIZED_HBN_AQUEOUS_NONBONDED_AND_EDGE_REFERENCE"
        ),
        "model_class": (
            "FIXED_TOPOLOGY_FUNCTIONALIZED_HBN_FORCE_FIELD"
        ),
        "priority": "P1",
    },
    {
        "source_id": "BAMANE_2023_IFF_R_BNNT",
        "short_name": "Bamane/IFF-R 2023",
        "doi": "10.1021/acsanm.2c05285",
        "route_role": (
            "CURVED_BNNT_STRUCTURE_AND_MECHANICS_REFERENCE"
        ),
        "model_class": (
            "IFF_R_FIXED_TOPOLOGY_BNNT_FORCE_FIELD"
        ),
        "priority": "P1",
    },
    {
        "source_id": "RAJAN_2018_HBN_FORCE_FIELD",
        "short_name": "Rajan 2018",
        "doi": "10.1021/acs.jpclett.7b03443",
        "route_role": (
            "BASELINE_HBN_BONDED_AND_NONBONDED_REFERENCE"
        ),
        "model_class": (
            "FIXED_TOPOLOGY_HBN_FORCE_FIELD"
        ),
        "priority": "P1",
    },
]

DOMAINS = [
    "PARENT_HBN_BULK_LIKE",
    "CURVED_BNNT_MECHANICS",
    "HYDROGENATED_OR_FUNCTIONALIZED_EDGE",
    "WATER_AND_ION_INTERFACE",
    "R2_FOUR_ATOM_BN_BRIDGE_AND_JUNCTION",
]

R2_SCOPE_KEYWORDS = {
    "PARENT_HBN_BULK_LIKE": (
        "PARENT",
        "BULK",
        "SCAFFOLD",
    ),
    "CURVED_BNNT_MECHANICS": (
        "PARENT",
        "CURVED",
        "BNNT",
        "SCAFFOLD",
    ),
    "HYDROGENATED_OR_FUNCTIONALIZED_EDGE": (
        "ANNULUS",
        "EDGE",
        "RIM",
        "PASSIV",
        "TERMIN",
        "HYDROGEN",
    ),
    "WATER_AND_ION_INTERFACE": (
        "WATER",
        "AQUEOUS",
        "INTERFACE",
        "ION",
        "SOLVENT",
    ),
    "R2_FOUR_ATOM_BN_BRIDGE_AND_JUNCTION": (
        "FOUR_ATOM",
        "BRIDGE",
        "JUNCTION",
    ),
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


def read_csv(
    path: Path,
) -> tuple[list[str], list[dict[str, str]]]:
    require_file(path)

    with path.open(
        "r",
        newline="",
        encoding="utf-8",
    ) as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fields = list(reader.fieldnames or [])

    if not fields:
        raise RuntimeError(
            f"No CSV header found: {path}"
        )

    return fields, rows


def read_one(path: Path) -> dict[str, str]:
    _, rows = read_csv(path)

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected one row in {path}; found {len(rows)}"
        )

    return rows[0]


def write_csv(
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


def normalize(text: str) -> str:
    return re.sub(
        r"\s+",
        " ",
        text,
    ).strip()


def row_text(row: dict[str, str]) -> str:
    return normalize(
        " | ".join(
            f"{key}={value}"
            for key, value in row.items()
            if value and value.strip()
        )
    )


def rows_matching_doi(
    rows: list[dict[str, str]],
    doi: str,
) -> list[dict[str, str]]:
    doi_lower = doi.lower()

    return [
        row
        for row in rows
        if doi_lower in row_text(row).lower()
    ]


def rows_matching_source(
    rows: list[dict[str, str]],
    source: dict[str, str],
) -> list[dict[str, str]]:
    doi_matches = rows_matching_doi(
        rows,
        source["doi"],
    )

    if doi_matches:
        return doi_matches

    name_tokens = [
        token.lower()
        for token in re.split(
            r"[^A-Za-z0-9]+",
            source["short_name"],
        )
        if len(token) >= 5
    ]

    matches = []

    for row in rows:
        text = row_text(row).lower()

        if any(
            token in text
            for token in name_tokens
        ):
            matches.append(row)

    return matches


def compact_rows(
    rows: list[dict[str, str]],
    maximum_rows: int = 5,
    maximum_chars: int = 500,
) -> str:
    texts = [
        row_text(row)
        for row in rows[:maximum_rows]
    ]

    joined = " || ".join(texts)

    if len(joined) > maximum_chars:
        return joined[:maximum_chars - 3] + "..."

    return joined


def boolean_text(value: str) -> bool:
    return value.strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
        "passed",
        "retrieved",
        "available",
    }


def infer_artifact_status(
    source: dict[str, str],
    prioritized_rows: list[dict[str, str]],
    retrieval_rows: list[dict[str, str]],
    file_rows: list[dict[str, str]],
) -> tuple[str, str, str]:
    all_text = " ".join(
        row_text(row)
        for row in (
            prioritized_rows
            + retrieval_rows
            + file_rows
        )
    ).lower()

    retrieved_terms = (
        "retrieved",
        "binary",
        "pdf",
        "text",
        "archive",
        "parameter file",
        "dataset",
    )

    failure_terms = (
        "failed",
        "403",
        "202",
        "empty",
        "missing",
        "landing",
        "access page",
    )

    retrieved_signal = any(
        term in all_text
        for term in retrieved_terms
    )

    failure_signal = any(
        term in all_text
        for term in failure_terms
    )

    real_file_signal = any(
        row
        for row in file_rows
        if any(
            key.lower() in {
                "size_bytes",
                "bytes",
                "sha256",
                "actual_kind",
                "file",
                "path",
                "filename",
            }
            and value.strip()
            for key, value in row.items()
        )
    )

    if real_file_signal and retrieved_signal and not failure_signal:
        return (
            "LOCAL_ARTIFACT_PRESENT_REQUIRES_CONTENT_AUDIT",
            "YES",
            "A local candidate artifact is recorded, but scientific content must still be audited.",
        )

    if real_file_signal and retrieved_signal:
        return (
            "PARTIAL_OR_AMBIGUOUS_LOCAL_ARTIFACT",
            "PARTIAL",
            "At least one candidate file is recorded, but retrieval status is mixed or ambiguous.",
        )

    if failure_signal:
        return (
            "NOT_RETRIEVED_OR_ACCESS_UNRESOLVED",
            "NO",
            "The prior retrieval attempt did not yield a verified scientific artifact.",
        )

    return (
        "RETRIEVAL_STATUS_UNRESOLVED",
        "NO",
        "No verified local primary artifact could be established from the existing ledgers.",
    )


def source_domain_assessment(
    source_id: str,
    domain: str,
) -> tuple[str, str, str]:
    if source_id == "GHORAI_2025_FUNCTIONALIZED_HBN_NANOPORE":
        mapping = {
            "PARENT_HBN_BULK_LIKE": (
                "SECONDARY_SUPPORT_ONLY",
                "NOT_ESTABLISHED",
                "The source is prioritized for functionalized nanopore environments rather than general BNNT bulk mechanics.",
            ),
            "CURVED_BNNT_MECHANICS": (
                "NO_DIRECT_EVIDENCE_YET",
                "NOT_ESTABLISHED",
                "Curved BNNT equilibrium mechanics must not be inferred from nanopore functionalization.",
            ),
            "HYDROGENATED_OR_FUNCTIONALIZED_EDGE": (
                "PRIMARY_RELEVANCE",
                "SOURCE_ARTIFACT_AUDIT_REQUIRED",
                "This is the most relevant candidate for H/OH-functionalized h-BN local environments.",
            ),
            "WATER_AND_ION_INTERFACE": (
                "PRIMARY_RELEVANCE",
                "SOURCE_ARTIFACT_AUDIT_REQUIRED",
                "This source is prioritized for aqueous and ion-interaction parameters.",
            ),
            "R2_FOUR_ATOM_BN_BRIDGE_AND_JUNCTION": (
                "NO_DIRECT_REFERENCE_CONFIGURATION",
                "NOT_ESTABLISHED",
                "The R2 B-N-B-N bridge and its attachment junctions are not assumed to be covered.",
            ),
        }

    elif source_id == "BAMANE_2023_IFF_R_BNNT":
        mapping = {
            "PARENT_HBN_BULK_LIKE": (
                "PRIMARY_RELEVANCE",
                "SOURCE_ARTIFACT_AUDIT_REQUIRED",
                "The source is prioritized as a transferable BNNT structural-mechanics candidate.",
            ),
            "CURVED_BNNT_MECHANICS": (
                "PRIMARY_RELEVANCE",
                "SOURCE_ARTIFACT_AUDIT_REQUIRED",
                "This is the most relevant candidate for curved BNNT equilibrium mechanics.",
            ),
            "HYDROGENATED_OR_FUNCTIONALIZED_EDGE": (
                "POSSIBLE_PARTIAL_RELEVANCE",
                "NOT_ESTABLISHED",
                "Exact edge termination coverage must be demonstrated from the parameter artifacts.",
            ),
            "WATER_AND_ION_INTERFACE": (
                "POSSIBLE_NONBONDED_RELEVANCE",
                "NOT_ESTABLISHED",
                "Water compatibility cannot be inferred without explicit source validation.",
            ),
            "R2_FOUR_ATOM_BN_BRIDGE_AND_JUNCTION": (
                "NO_DIRECT_REFERENCE_CONFIGURATION",
                "NOT_ESTABLISHED",
                "Standard BNNT parameters are not assumed to cover the reconstructed bridge.",
            ),
        }

    elif source_id == "RAJAN_2018_HBN_FORCE_FIELD":
        mapping = {
            "PARENT_HBN_BULK_LIKE": (
                "PRIMARY_RELEVANCE",
                "SOURCE_ARTIFACT_AUDIT_REQUIRED",
                "This source is prioritized as a baseline fixed-topology h-BN model.",
            ),
            "CURVED_BNNT_MECHANICS": (
                "POSSIBLE_TRANSFERABILITY",
                "SOURCE_ARTIFACT_AUDIT_REQUIRED",
                "Curvature applicability must be verified from the article and validation set.",
            ),
            "HYDROGENATED_OR_FUNCTIONALIZED_EDGE": (
                "NO_DIRECT_EVIDENCE_YET",
                "NOT_ESTABLISHED",
                "Hydrogenated edge terms must not be inferred without explicit parameter evidence.",
            ),
            "WATER_AND_ION_INTERFACE": (
                "NO_DIRECT_EVIDENCE_YET",
                "NOT_ESTABLISHED",
                "Aqueous compatibility has not yet been established.",
            ),
            "R2_FOUR_ATOM_BN_BRIDGE_AND_JUNCTION": (
                "NO_DIRECT_REFERENCE_CONFIGURATION",
                "NOT_ESTABLISHED",
                "No direct R2 bridge analogue is assumed.",
            ),
        }

    else:
        raise ValueError(
            source_id
        )

    return mapping[domain]


def environment_matches_domain(
    row: dict[str, str],
    domain: str,
) -> bool:
    text = row_text(row).upper()

    return any(
        keyword in text
        for keyword in R2_SCOPE_KEYWORDS[
            domain
        ]
    )


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    data: dict[str, dict[str, Any]] = {}

    for name, path in INPUTS.items():
        fields, rows = read_csv(path)

        data[name] = {
            "path": path,
            "fields": fields,
            "rows": rows,
        }

    gate39_summary = data[
        "gate39_summary"
    ]["rows"][0]

    if gate39_summary.get(
        "decision"
    ) != EXPECTED_GATE39_DECISION:
        raise RuntimeError(
            "Gate 3S.1 decision is not accepted."
        )

    environment_rows = data[
        "gate39_environment_matrix"
    ]["rows"]

    term_rows = data[
        "gate39_term_matrix"
    ]["rows"]

    source_records = []
    domain_rows = []
    artifact_rows = []
    scope_mapping_rows = []

    source_ledger = data[
        "gate34_source_ledger"
    ]["rows"]

    domain_coverage = data[
        "gate34_domain_coverage"
    ]["rows"]

    model_decisions = data[
        "gate34_model_decisions"
    ]["rows"]

    prioritized_artifacts = data[
        "gate34_prioritized_artifacts"
    ]["rows"]

    retrieval_ledger = data[
        "gate35_retrieval_ledger"
    ]["rows"]

    file_inventory = data[
        "gate35_file_inventory"
    ]["rows"]

    for source in SOURCES:
        ledger_matches = rows_matching_source(
            source_ledger,
            source,
        )

        coverage_matches = rows_matching_source(
            domain_coverage,
            source,
        )

        decision_matches = rows_matching_source(
            model_decisions,
            source,
        )

        priority_matches = rows_matching_source(
            prioritized_artifacts,
            source,
        )

        retrieval_matches = rows_matching_source(
            retrieval_ledger,
            source,
        )

        file_matches = rows_matching_source(
            file_inventory,
            source,
        )

        (
            artifact_status,
            verified_artifact_available,
            artifact_interpretation,
        ) = infer_artifact_status(
            source,
            priority_matches,
            retrieval_matches,
            file_matches,
        )

        source_found = bool(
            ledger_matches
            or coverage_matches
            or decision_matches
            or priority_matches
            or retrieval_matches
            or file_matches
        )

        source_records.append(
            {
                **source,
                "source_found_in_curated_ledgers": source_found,
                "source_ledger_matches": len(
                    ledger_matches
                ),
                "domain_coverage_matches": len(
                    coverage_matches
                ),
                "model_decision_matches": len(
                    decision_matches
                ),
                "prioritized_artifact_matches": len(
                    priority_matches
                ),
                "retrieval_ledger_matches": len(
                    retrieval_matches
                ),
                "retrieved_file_inventory_matches": len(
                    file_matches
                ),
                "artifact_status": artifact_status,
                "verified_local_scientific_artifact_available": (
                    verified_artifact_available
                ),
                "artifact_interpretation": (
                    artifact_interpretation
                ),
                "curated_source_evidence": compact_rows(
                    ledger_matches
                    + coverage_matches
                    + decision_matches
                ),
                "retrieval_evidence": compact_rows(
                    priority_matches
                    + retrieval_matches
                    + file_matches
                ),
                "parameter_adoption_authorized": False,
            }
        )

        required_artifacts = [
            {
                "artifact_class": "PRIMARY_ARTICLE",
                "required_for": (
                    "Scientific domain, validation protocol and limitations"
                ),
            },
            {
                "artifact_class": "SUPPORTING_INFORMATION",
                "required_for": (
                    "Detailed force-field form, fitting data and validation"
                ),
            },
            {
                "artifact_class": "PARAMETER_OR_DATASET_FILES",
                "required_for": (
                    "Exact atom types, bonded terms, charges and nonbonded parameters"
                ),
            },
        ]

        for artifact in required_artifacts:
            artifact_rows.append(
                {
                    "source_id": source[
                        "source_id"
                    ],
                    "short_name": source[
                        "short_name"
                    ],
                    "doi": source[
                        "doi"
                    ],
                    **artifact,
                    "current_overall_artifact_status": artifact_status,
                    "verified_local_scientific_artifact_available": (
                        verified_artifact_available
                    ),
                    "content_audit_complete": False,
                    "artifact_requirement_closed": False,
                }
            )

        for domain in DOMAINS:
            (
                evidence_class,
                coverage_status,
                interpretation,
            ) = source_domain_assessment(
                source[
                    "source_id"
                ],
                domain,
            )

            matched_environments = [
                row
                for row in environment_rows
                if environment_matches_domain(
                    row,
                    domain,
                )
            ]

            matched_term_ids = sorted(
                {
                    row.get(
                        "term_id",
                        "",
                    )
                    for row in term_rows
                    if any(
                        keyword in row_text(
                            row
                        ).upper()
                        for keyword in R2_SCOPE_KEYWORDS[
                            domain
                        ]
                    )
                }
                - {
                    "",
                }
            )

            domain_rows.append(
                {
                    "source_id": source[
                        "source_id"
                    ],
                    "short_name": source[
                        "short_name"
                    ],
                    "doi": source[
                        "doi"
                    ],
                    "domain": domain,
                    "evidence_class": evidence_class,
                    "R2_coverage_status": coverage_status,
                    "interpretation": interpretation,
                    "matched_R2_environment_classes": len(
                        matched_environments
                    ),
                    "matched_R2_bonded_term_classes": len(
                        matched_term_ids
                    ),
                    "direct_R2_validation_established": False,
                    "parameter_assignment_authorized": False,
                }
            )

            for environment in matched_environments:
                scope_mapping_rows.append(
                    {
                        "source_id": source[
                            "source_id"
                        ],
                        "short_name": source[
                            "short_name"
                        ],
                        "domain": domain,
                        "environment_id": environment.get(
                            "environment_id",
                            "",
                        ),
                        "element": environment.get(
                            "element",
                            "",
                        ),
                        "node_type": environment.get(
                            "node_type",
                            "",
                        ),
                        "region": environment.get(
                            "region",
                            "",
                        ),
                        "environment_count": environment.get(
                            "count",
                            "",
                        ),
                        "current_scientific_risk": environment.get(
                            "scientific_risk",
                            "",
                        ),
                        "source_coverage_status": coverage_status,
                        "source_transferability_established": False,
                    }
                )

    write_csv(
        SOURCE_RECORDS,
        source_records,
    )

    write_csv(
        SOURCE_DOMAIN_MATRIX,
        domain_rows,
    )

    write_csv(
        ARTIFACT_GAP_MATRIX,
        artifact_rows,
    )

    write_csv(
        R2_SCOPE_MAPPING,
        scope_mapping_rows,
    )

    source_found_count = sum(
        1
        for row in source_records
        if row[
            "source_found_in_curated_ledgers"
        ]
    )

    verified_artifact_source_count = sum(
        1
        for row in source_records
        if row[
            "verified_local_scientific_artifact_available"
        ]
        == "YES"
    )

    partial_artifact_source_count = sum(
        1
        for row in source_records
        if row[
            "verified_local_scientific_artifact_available"
        ]
        == "PARTIAL"
    )

    directly_validated_domain_count = sum(
        1
        for row in domain_rows
        if row[
            "direct_R2_validation_established"
        ]
    )

    primary_relevance_rows = [
        row
        for row in domain_rows
        if row[
            "evidence_class"
        ]
        == "PRIMARY_RELEVANCE"
    ]

    route_rows = [
        {
            "route": "FUNCTIONALIZED_HBN_FIXED_TOPOLOGY",
            "sources_audited": len(
                SOURCES
            ),
            "sources_found_in_curated_ledgers": (
                source_found_count
            ),
            "sources_with_verified_local_artifacts": (
                verified_artifact_source_count
            ),
            "sources_with_partial_local_artifacts": (
                partial_artifact_source_count
            ),
            "primary_relevance_domain_rows": len(
                primary_relevance_rows
            ),
            "directly_validated_R2_domains": (
                directly_validated_domain_count
            ),
            "bulk_or_curved_hBN_candidate": (
                "BAMANE_IFF_R_AND_RAJAN"
            ),
            "functionalized_edge_and_aqueous_candidate": (
                "GHORAI_2025"
            ),
            "R2_bridge_candidate": (
                "NONE_DIRECTLY_IDENTIFIED"
            ),
            "route_status": (
                "PROMISING_BUT_INCOMPLETE_PRIMARY_ARTIFACT_AUDIT"
            ),
            "force_field_route_selected": False,
            "parameter_adoption_authorized": False,
            "required_next_step": (
                "RETRIEVE_AND_CONTENT_AUDIT_GHORAI_IFF_R_AND_RAJAN_"
                "PRIMARY_ARTIFACTS"
            ),
        }
    ]

    write_csv(
        ROUTE_ASSESSMENT,
        route_rows,
    )

    gates = {
        "gate39_decision_is_accepted": (
            gate39_summary.get(
                "decision"
            )
            == EXPECTED_GATE39_DECISION
        ),
        "three_priority_sources_defined": (
            len(
                SOURCES
            )
            == 3
        ),
        "all_priority_sources_found_in_curated_ledgers": (
            source_found_count
            == 3
        ),
        "five_domains_assessed_for_each_source": (
            len(
                domain_rows
            )
            == 15
        ),
        "three_required_artifact_classes_defined_per_source": (
            len(
                artifact_rows
            )
            == 9
        ),
        "no_direct_R2_validation_inferred": (
            directly_validated_domain_count
            == 0
        ),
        "no_fixed_topology_route_selected": True,
        "no_parameter_adoption_or_simulation_performed": True,
    }

    failed_gates = [
        name
        for name, passed in gates.items()
        if not passed
    ]

    decision = (
        PASS_DECISION
        if not failed_gates
        else REVIEW_DECISION
    )

    required_next_step = (
        "RETRIEVE_AND_CONTENT_AUDIT_GHORAI_IFF_R_AND_RAJAN_"
        "PRIMARY_ARTIFACTS"
        if not failed_gates
        else
        "REVIEW_FIXED_TOPOLOGY_SOURCE_PROVENANCE_AND_LEDGER_MATCHING"
    )

    summary = {
        "decision": decision,
        "priority_sources": len(
            SOURCES
        ),
        "sources_found_in_curated_ledgers": (
            source_found_count
        ),
        "sources_with_verified_local_artifacts": (
            verified_artifact_source_count
        ),
        "sources_with_partial_local_artifacts": (
            partial_artifact_source_count
        ),
        "source_domain_rows": len(
            domain_rows
        ),
        "primary_relevance_domain_rows": len(
            primary_relevance_rows
        ),
        "directly_validated_R2_domains": (
            directly_validated_domain_count
        ),
        "artifact_requirement_rows": len(
            artifact_rows
        ),
        "R2_scope_mapping_rows": len(
            scope_mapping_rows
        ),
        "fixed_topology_route_selected": False,
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
        "required_next_step": required_next_step,
    }

    write_csv(
        SUMMARY,
        [
            summary
        ],
    )

    write_csv(
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

    for name, record in INPUTS.items():
        path = record

        manifest_rows.append(
            {
                "role": name,
                "file": relative(
                    path
                ),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(
                    path
                ),
            }
        )

    write_csv(
        MANIFEST,
        manifest_rows,
    )

    JSON_OUT.write_text(
        json.dumps(
            {
                "summary": summary,
                "source_records": source_records,
                "domain_matrix": domain_rows,
                "artifact_gaps": artifact_rows,
                "route_assessment": route_rows,
                "gates": gates,
                "limitations": [
                    (
                        "The audit uses curated provenance and retrieval "
                        "ledgers already present in the repository."
                    ),
                    (
                        "Source relevance is not equivalent to exact "
                        "parameter coverage."
                    ),
                    (
                        "No direct R2 coverage is inferred without article, "
                        "SI and parameter-file content validation."
                    ),
                    (
                        "The R2 four-atom bridge remains without a direct "
                        "fixed-topology literature analogue."
                    ),
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    report_lines = [
        "# R2 Functionalized h-BN Fixed-Topology Primary-Source Audit — Day025",
        "",
        "## Priority sources",
        "",
    ]

    for row in source_records:
        report_lines.extend(
            [
                f"### {row['short_name']}",
                "",
                f"- DOI: `{row['doi']}`",
                f"- Role: **{row['route_role']}**",
                f"- Found in curated ledgers: **{row['source_found_in_curated_ledgers']}**",
                f"- Artifact status: **{row['artifact_status']}**",
                f"- Verified local scientific artifact: "
                f"**{row['verified_local_scientific_artifact_available']}**",
                "",
            ]
        )

    report_lines.extend(
        [
            "## Route interpretation",
            "",
            "- Ghorai 2025 is the primary candidate for functionalized "
            "edge and aqueous-interface evidence.",
            "- Bamane/IFF-R 2023 is the primary curved-BNNT structural "
            "and mechanics candidate.",
            "- Rajan 2018 is the baseline fixed-topology h-BN candidate.",
            "- No source is assumed to cover the R2 four-atom bridge.",
            "- No source currently establishes direct R2 validation.",
            "",
            "## Decision",
            "",
            f"- Decision: **{decision}**",
            (
                "- Failed gates: **NONE**"
                if not failed_gates
                else
                f"- Failed gates: **{' | '.join(failed_gates)}**"
            ),
            "- Fixed-topology route selected: **NO**",
            "- Parameter adoption authorized: **NO**",
            "- Topology/minimization/MD/QM: **NOT AUTHORIZED**",
            f"- Required next step: `{required_next_step}`",
            "",
        ]
    )

    REPORT.write_text(
        "\n".join(
            report_lines
        ),
        encoding="utf-8",
    )

    print(
        "Day025 R2 functionalized h-BN fixed-topology "
        "primary-source audit completed."
    )

    print(
        "Priority sources / found in curated ledgers: "
        f"{len(SOURCES)}/"
        f"{source_found_count}"
    )

    print(
        "Sources with verified / partial local artifacts: "
        f"{verified_artifact_source_count}/"
        f"{partial_artifact_source_count}"
    )

    print(
        "Domain rows / primary-relevance rows / "
        "directly validated R2 domains: "
        f"{len(domain_rows)}/"
        f"{len(primary_relevance_rows)}/"
        f"{directly_validated_domain_count}"
    )

    print(
        "Artifact requirements / R2 mapping rows: "
        f"{len(artifact_rows)}/"
        f"{len(scope_mapping_rows)}"
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
        "Fixed-topology route selected: NO"
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

    print()
    print("SOURCE ASSESSMENT")

    for row in source_records:
        print(
            f"- {row['short_name']}: "
            f"found={row['source_found_in_curated_ledgers']} "
            f"artifact={row['artifact_status']} "
            f"verified_local="
            f"{row['verified_local_scientific_artifact_available']}"
        )

    for path in (
        SOURCE_RECORDS,
        SOURCE_DOMAIN_MATRIX,
        ARTIFACT_GAP_MATRIX,
        R2_SCOPE_MAPPING,
        ROUTE_ASSESSMENT,
        SUMMARY,
        GATES,
        MANIFEST,
        JSON_OUT,
        REPORT,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if failed_gates:
        raise RuntimeError(
            "Functionalized h-BN primary-source audit requires review."
        )


if __name__ == "__main__":
    main()
