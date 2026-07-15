#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
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

PREFLIGHT = (
    DAY025
    / "38_r2_force_field_route_comparison_preflight"
)

OUT = (
    DAY025
    / "39_r2_force_field_route_coverage_and_risk_matrix"
)

INPUTS = {
    "preflight_summary": (
        PREFLIGHT
        / "r2_force_field_route_comparison_preflight_summary.csv"
    ),
    "environment_inventory": (
        DAY024
        / "29_r2_chemical_realizability_and_parameterization_scope"
        / "r2_local_chemical_environment_inventory.csv"
    ),
    "critical_centers": (
        DAY024
        / "29_r2_chemical_realizability_and_parameterization_scope"
        / "r2_parameterization_critical_centers.csv"
    ),
    "scope": (
        DAY024
        / "29_r2_chemical_realizability_and_parameterization_scope"
        / "r2_parameterization_scope.csv"
    ),
    "bonds": (
        DAY024
        / "30_r2_force_field_coverage_preflight"
        / "r2_required_bond_terms.csv"
    ),
    "angles": (
        DAY024
        / "30_r2_force_field_coverage_preflight"
        / "r2_required_angle_terms.csv"
    ),
    "torsions": (
        DAY024
        / "30_r2_force_field_coverage_preflight"
        / "r2_required_torsion_terms.csv"
    ),
    "impropers": (
        DAY024
        / "30_r2_force_field_coverage_preflight"
        / "r2_required_improper_centers.csv"
    ),
    "qm_fragments": (
        DAY024
        / "30_r2_force_field_coverage_preflight"
        / "r2_preliminary_qm_fragment_classes.csv"
    ),
    "lele_terms": (
        DAY024
        / "37_r2_lele_si_parameter_file_audit"
        / "lele_2022_reaxff_bnh_term_inventory.csv"
    ),
    "lele_domains": (
        DAY024
        / "37_r2_lele_si_parameter_file_audit"
        / "lele_2022_R2_domain_assessment.csv"
    ),
}

ENV_MATRIX = OUT / "r2_environment_route_coverage_and_risk_matrix.csv"
TERM_MATRIX = OUT / "r2_bonded_term_route_coverage_and_risk_matrix.csv"
IMPROPER_MATRIX = OUT / "r2_improper_region_route_coverage_matrix.csv"
ROUTE_MATRIX = OUT / "r2_force_field_route_comparison_matrix.csv"
QM_MATRIX = OUT / "r2_qm_reference_requirement_matrix.csv"
CRITICAL_REGION_SUMMARY = OUT / "r2_critical_center_region_summary.csv"
SUMMARY = OUT / "r2_force_field_route_coverage_and_risk_summary.csv"
GATES = OUT / "r2_force_field_route_coverage_and_risk_gates.csv"
JSON_OUT = OUT / "r2_force_field_route_coverage_and_risk.json"
REPORT = OUT / "R2_FORCE_FIELD_ROUTE_COVERAGE_AND_RISK_MATRIX_DAY025.md"

EXPECTED_PREFLIGHT_DECISION = (
    "R2_FORCE_FIELD_ROUTE_COMPARISON_DATA_CONTRACT_VALIDATED"
)

PASS_DECISION = (
    "R2_FORCE_FIELD_ROUTE_COVERAGE_AND_RISK_MATRIX_BUILT_"
    "NO_ROUTE_YET_AUTHORIZED"
)

REVIEW_DECISION = (
    "R2_FORCE_FIELD_ROUTE_COVERAGE_AND_RISK_MATRIX_REQUIRES_REVIEW"
)


def require_file(path: Path) -> None:
    if not path.is_file():
        raise RuntimeError(f"Missing required input: {path}")

    if path.stat().st_size == 0:
        raise RuntimeError(f"Empty required input: {path}")


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(
        "r",
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(f"No rows found in {path}")

    return rows


def write_csv(
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


def as_int(value: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def atom_type_element(atom_type: str) -> str:
    return atom_type.split("|", 1)[0].strip()


def atom_type_label(atom_type: str) -> str:
    parts = atom_type.split("|", 1)

    if len(parts) == 1:
        return parts[0].strip()

    return parts[1].strip()


def classify_region_risk(
    region: str,
    node_type: str,
) -> tuple[str, str]:
    text = f"{region} {node_type}".upper()

    if "FOUR_ATOM_BRIDGE" in text:
        return (
            "HIGH",
            "Novel four-atom B-N-B-N bridge and attachment environment.",
        )

    if (
        "ANNULUS_EDGE" in text
        or "OUTER_BOUNDARY" in text
        or "INNER_BOUNDARY" in text
        or "RIM" in text
    ):
        return (
            "HIGH",
            "Reconstructed polar edge/rim environment with termination sensitivity.",
        )

    if "ANNULUS" in text:
        return (
            "MEDIUM_HIGH",
            "Reconstructed annulus environment outside ordinary bulk h-BN.",
        )

    if (
        "PASSIVANT" in text
        or "TERMINAL" in text
        or "HYDROGEN" in text
    ):
        return (
            "MEDIUM_HIGH",
            "Hydrogen termination chemistry requires source-specific validation.",
        )

    if (
        "PARENT" in text
        or "BULK" in text
        or "SCAFFOLD" in text
    ):
        return (
            "MEDIUM",
            "Bulk-like curved h-BN may be transferable but curvature must be validated.",
        )

    return (
        "HIGH",
        "Environment does not map unambiguously to a validated literature domain.",
    )


def lele_environment_status(
    element: str,
    region: str,
    node_type: str,
) -> tuple[str, str]:
    text = f"{region} {node_type}".upper()

    if element not in {"B", "N", "H"}:
        return (
            "OUTSIDE_PRIMARY_ELEMENTAL_SCOPE",
            "Element is outside the B/N/H subset relevant to R2.",
        )

    if "FOUR_ATOM_BRIDGE" in text:
        return (
            "ELEMENTAL_RECORDS_PRESENT_TRANSFERABILITY_NOT_ESTABLISHED",
            "Lele contains B/N/H reactive records but no direct R2 bridge validation.",
        )

    if (
        "ANNULUS" in text
        or "RIM" in text
        or "BOUNDARY" in text
    ):
        return (
            "ELEMENTAL_RECORDS_PRESENT_TRANSFERABILITY_NOT_ESTABLISHED",
            "Lele was trained for gas-phase synthesis, not reconstructed equilibrium rims.",
        )

    if (
        "PASSIVANT" in text
        or "HYDROGEN" in text
        or "TERMINAL" in text
    ):
        return (
            "B_H_OR_N_H_REACTIVE_RECORDS_PRESENT_TRANSFERABILITY_NOT_ESTABLISHED",
            "Reactive B-H/N-H chemistry exists, but equilibrium termination behavior is unvalidated.",
        )

    if (
        "PARENT" in text
        or "BULK" in text
        or "SCAFFOLD" in text
    ):
        return (
            "BN_REACTIVE_RECORDS_PRESENT_EQUILIBRIUM_CURVATURE_NOT_VALIDATED",
            "BN formation capability is not equivalent to equilibrium BNNT mechanics.",
        )

    return (
        "ELEMENTAL_RECORDS_PRESENT_TRANSFERABILITY_NOT_ESTABLISHED",
        "No direct validation match identified.",
    )


def fixed_topology_environment_status(
    region: str,
    node_type: str,
) -> tuple[str, str]:
    text = f"{region} {node_type}".upper()

    if (
        "PARENT" in text
        or "BULK" in text
        or "SCAFFOLD" in text
    ):
        return (
            "POTENTIALLY_TRANSFERABLE_SOURCE_EVIDENCE_PENDING",
            "Bulk or curved h-BN fixed-topology models may cover this domain.",
        )

    if (
        "PASSIVANT" in text
        or "HYDROGEN" in text
        or "TERMINAL" in text
        or "ANNULUS" in text
        or "RIM" in text
        or "BRIDGE" in text
    ):
        return (
            "FUNCTIONALIZED_HBN_SOURCE_EVIDENCE_PENDING",
            "Requires direct audit of H/OH-functionalized h-BN source parameters.",
        )

    return (
        "SOURCE_EVIDENCE_PENDING",
        "No primary-source fixed-topology evidence has yet been mapped.",
    )


def qm_environment_status(
    risk: str,
) -> tuple[str, str]:
    if risk in {"HIGH", "MEDIUM_HIGH"}:
        return (
            "QM_REFERENCE_REQUIRED_HIGH_PRIORITY",
            "Dedicated fragment reference data are required before parameter assignment.",
        )

    return (
        "QM_REFERENCE_RECOMMENDED",
        "Reference calculations are recommended to test transferability.",
    )


def term_elements(row: dict[str, str], term_kind: str) -> list[str]:
    if term_kind == "BOND":
        fields = [
            "atom_type_1",
            "atom_type_2",
        ]

    elif term_kind == "ANGLE":
        fields = [
            "atom_type_1",
            "center_type",
            "atom_type_3",
        ]

    elif term_kind == "TORSION":
        fields = [
            "atom_type_1",
            "atom_type_2",
            "atom_type_3",
            "atom_type_4",
        ]

    else:
        raise ValueError(term_kind)

    return [
        atom_type_element(row[field])
        for field in fields
    ]


def term_labels(row: dict[str, str], term_kind: str) -> list[str]:
    if term_kind == "BOND":
        fields = [
            "atom_type_1",
            "atom_type_2",
        ]

    elif term_kind == "ANGLE":
        fields = [
            "atom_type_1",
            "center_type",
            "atom_type_3",
        ]

    elif term_kind == "TORSION":
        fields = [
            "atom_type_1",
            "atom_type_2",
            "atom_type_3",
            "atom_type_4",
        ]

    else:
        raise ValueError(term_kind)

    return [
        atom_type_label(row[field])
        for field in fields
    ]


def classify_term_risk(labels: list[str]) -> tuple[str, str]:
    text = " ".join(labels).upper()

    if "FOUR_ATOM_BRIDGE" in text:
        return (
            "HIGH",
            "Term directly involves the novel four-atom bridge.",
        )

    if (
        "ANNULUS_OUTER_BOUNDARY" in text
        or "ANNULUS_INNER_BOUNDARY" in text
        or "ANNULUS_EDGE" in text
        or "RIM" in text
    ):
        return (
            "HIGH",
            "Term involves reconstructed edge or boundary chemistry.",
        )

    if "ANNULUS" in text:
        return (
            "MEDIUM_HIGH",
            "Term involves reconstructed annulus chemistry.",
        )

    if (
        "PASSIVANT" in text
        or "TERMINAL" in text
    ):
        return (
            "MEDIUM_HIGH",
            "Term involves B-H or N-H termination.",
        )

    if (
        "PARENT" in text
        or "BULK" in text
        or "SCAFFOLD" in text
    ):
        return (
            "MEDIUM",
            "Bulk-like h-BN term requires curvature transferability validation.",
        )

    return (
        "HIGH",
        "Term has no unambiguous validated literature analogue.",
    )


def lele_term_status(
    term_kind: str,
    elements: list[str],
) -> tuple[str, str]:
    element_set = set(elements)

    if not element_set.issubset({"B", "N", "H"}):
        return (
            "OUTSIDE_R2_BNH_SUBSET",
            "The term contains elements outside the R2 B/N/H chemistry.",
        )

    if term_kind == "BOND":
        if element_set == {"B", "N"}:
            return (
                "BN_REACTIVE_RECORD_PRESENT_TRANSFERABILITY_NOT_ESTABLISHED",
                "Lele contains B-N reactive bond parameters.",
            )

        if element_set == {"B", "H"}:
            return (
                "BH_REACTIVE_RECORD_PRESENT_TRANSFERABILITY_NOT_ESTABLISHED",
                "Lele contains B-H reactive bond parameters.",
            )

        if element_set == {"N", "H"}:
            return (
                "NH_REACTIVE_RECORD_PRESENT_TRANSFERABILITY_NOT_ESTABLISHED",
                "Lele contains N-H reactive bond parameters.",
            )

    return (
        "ELEMENT_COMBINATION_POSSIBLY_REPRESENTED_ENVIRONMENT_NOT_VALIDATED",
        "Element-level ReaxFF terms do not retain R2 environment specificity.",
    )


def risk_numeric(risk: str) -> int:
    return {
        "LOW": 1,
        "MEDIUM": 2,
        "MEDIUM_HIGH": 3,
        "HIGH": 4,
    }.get(risk, 4)


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    data = {
        name: read_csv(path)
        for name, path in INPUTS.items()
    }

    preflight_summary = data[
        "preflight_summary"
    ][0]

    if preflight_summary.get(
        "decision"
    ) != EXPECTED_PREFLIGHT_DECISION:
        raise RuntimeError(
            "Gate 3S.0 preflight decision is not accepted."
        )

    critical_counts = Counter(
        row.get("region", "UNKNOWN")
        for row in data["critical_centers"]
    )

    critical_reason_counts: dict[str, Counter[str]] = defaultdict(Counter)

    for row in data["critical_centers"]:
        region = row.get(
            "region",
            "UNKNOWN",
        )

        reasons = row.get(
            "reasons",
            "",
        )

        for reason in reasons.split("|"):
            reason = reason.strip()

            if reason:
                critical_reason_counts[
                    region
                ][reason] += 1

    critical_region_rows = []

    for region, count in sorted(
        critical_counts.items()
    ):
        reasons = critical_reason_counts[
            region
        ]

        critical_region_rows.append(
            {
                "region": region,
                "critical_center_count": count,
                "reason_counts": " | ".join(
                    f"{reason}:{number}"
                    for reason, number in sorted(
                        reasons.items()
                    )
                ),
            }
        )

    write_csv(
        CRITICAL_REGION_SUMMARY,
        critical_region_rows,
    )

    environment_rows = []

    for row in data["environment_inventory"]:
        region = row.get(
            "region",
            "",
        )

        node_type = row.get(
            "node_type",
            "",
        )

        element = row.get(
            "element",
            "",
        )

        risk, risk_basis = classify_region_risk(
            region,
            node_type,
        )

        lele_status, lele_basis = lele_environment_status(
            element,
            region,
            node_type,
        )

        fixed_status, fixed_basis = fixed_topology_environment_status(
            region,
            node_type,
        )

        qm_status, qm_basis = qm_environment_status(
            risk
        )

        environment_rows.append(
            {
                **row,
                "critical_center_count_in_region": critical_counts.get(
                    region,
                    0,
                ),
                "scientific_risk": risk,
                "risk_score": risk_numeric(
                    risk
                ),
                "risk_basis": risk_basis,
                "lele_2022_reaxff_status": lele_status,
                "lele_2022_reaxff_basis": lele_basis,
                "fixed_topology_status": fixed_status,
                "fixed_topology_basis": fixed_basis,
                "QM_reference_status": qm_status,
                "QM_reference_basis": qm_basis,
                "force_field_assignment_authorized": False,
            }
        )

    write_csv(
        ENV_MATRIX,
        environment_rows,
    )

    term_rows = []

    term_sources = [
        (
            "BOND",
            "bond_term_id",
            data["bonds"],
        ),
        (
            "ANGLE",
            "angle_term_id",
            data["angles"],
        ),
        (
            "TORSION",
            "torsion_term_id",
            data["torsions"],
        ),
    ]

    for term_kind, id_field, rows in term_sources:
        for row in rows:
            elements = term_elements(
                row,
                term_kind,
            )

            labels = term_labels(
                row,
                term_kind,
            )

            risk, risk_basis = classify_term_risk(
                labels
            )

            lele_status, lele_basis = lele_term_status(
                term_kind,
                elements,
            )

            if risk == "MEDIUM":
                fixed_status = (
                    "POTENTIALLY_TRANSFERABLE_SOURCE_EVIDENCE_PENDING"
                )
            else:
                fixed_status = (
                    "FUNCTIONALIZED_OR_CUSTOM_TERM_EVIDENCE_PENDING"
                )

            if risk in {
                "HIGH",
                "MEDIUM_HIGH",
            }:
                qm_status = (
                    "QM_REFERENCE_REQUIRED"
                )
            else:
                qm_status = (
                    "QM_REFERENCE_RECOMMENDED"
                )

            term_rows.append(
                {
                    "term_kind": term_kind,
                    "term_id": row[id_field],
                    "elements": "-".join(
                        elements
                    ),
                    "environment_labels": " | ".join(
                        labels
                    ),
                    "count": row.get(
                        "count",
                        "",
                    ),
                    "minimum": (
                        row.get("minimum_nm")
                        or row.get("minimum_deg")
                        or ""
                    ),
                    "mean": (
                        row.get("mean_nm")
                        or row.get("mean_deg")
                        or ""
                    ),
                    "maximum": (
                        row.get("maximum_nm")
                        or row.get("maximum_deg")
                        or ""
                    ),
                    "scientific_risk": risk,
                    "risk_score": risk_numeric(
                        risk
                    ),
                    "risk_basis": risk_basis,
                    "lele_2022_reaxff_status": lele_status,
                    "lele_2022_reaxff_basis": lele_basis,
                    "fixed_topology_status": fixed_status,
                    "QM_reference_status": qm_status,
                    "parameter_assignment_authorized": False,
                }
            )

    write_csv(
        TERM_MATRIX,
        term_rows,
    )

    improper_groups: dict[
        tuple[str, str, str],
        int,
    ] = Counter()

    for row in data["impropers"]:
        key = (
            row.get(
                "center_element",
                "",
            ),
            row.get(
                "center_type",
                "",
            ),
            row.get(
                "region",
                "",
            ),
        )

        improper_groups[
            key
        ] += 1

    improper_rows = []

    for (
        center_element,
        center_type,
        region,
    ), count in sorted(
        improper_groups.items()
    ):
        risk, risk_basis = classify_region_risk(
            region,
            center_type,
        )

        lele_status, lele_basis = lele_environment_status(
            center_element,
            region,
            center_type,
        )

        fixed_status, fixed_basis = fixed_topology_environment_status(
            region,
            center_type,
        )

        qm_status, qm_basis = qm_environment_status(
            risk
        )

        improper_rows.append(
            {
                "center_element": center_element,
                "center_type": center_type,
                "region": region,
                "center_count": count,
                "scientific_risk": risk,
                "risk_basis": risk_basis,
                "lele_2022_reaxff_status": lele_status,
                "lele_2022_reaxff_basis": lele_basis,
                "fixed_topology_status": fixed_status,
                "fixed_topology_basis": fixed_basis,
                "QM_reference_status": qm_status,
                "QM_reference_basis": qm_basis,
                "improper_assignment_authorized": False,
            }
        )

    write_csv(
        IMPROPER_MATRIX,
        improper_rows,
    )

    high_environment_count = sum(
        1
        for row in environment_rows
        if row["scientific_risk"] == "HIGH"
    )

    medium_high_environment_count = sum(
        1
        for row in environment_rows
        if row["scientific_risk"] == "MEDIUM_HIGH"
    )

    medium_environment_count = sum(
        1
        for row in environment_rows
        if row["scientific_risk"] == "MEDIUM"
    )

    high_term_count = sum(
        1
        for row in term_rows
        if row["scientific_risk"] == "HIGH"
    )

    medium_high_term_count = sum(
        1
        for row in term_rows
        if row["scientific_risk"] == "MEDIUM_HIGH"
    )

    medium_term_count = sum(
        1
        for row in term_rows
        if row["scientific_risk"] == "MEDIUM"
    )

    directly_validated_R2_statuses = {
        "DIRECTLY_VALIDATED_FOR_R2",
        "VALIDATED_R2_ENVIRONMENT",
    }

    lele_directly_validated_environments = sum(
        1
        for row in environment_rows
        if row["lele_2022_reaxff_status"]
        in directly_validated_R2_statuses
    )

    fixed_directly_validated_environments = 0

    qm_fragment_rows = []

    for row in data["qm_fragments"]:
        priority = row.get(
            "priority",
            "",
        ).upper()

        authorization = row.get(
            "calculation_authorized",
            "",
        )

        qm_fragment_rows.append(
            {
                **row,
                "reference_route_status": (
                    "DEFINED_NOT_AUTHORIZED"
                ),
                "parameterization_decision_dependency": (
                    "REQUIRED_BEFORE_NOVEL_TERM_ADOPTION"
                    if priority == "HIGH"
                    else
                    "RECOMMENDED_FOR_TRANSFERABILITY_VALIDATION"
                ),
                "calculation_authorized_preserved": authorization,
            }
        )

    write_csv(
        QM_MATRIX,
        qm_fragment_rows,
    )

    route_rows = [
        {
            "route": "LELE_2022_REAXFF",
            "model_type": "Reactive bond-order force field",
            "primary_demonstrated_domain": (
                "Gas-phase B/N/H chemistry and high-temperature BN nanostructure formation"
            ),
            "elemental_scope_for_R2": "B/N/H PRESENT",
            "direct_R2_environment_validation": "NONE IDENTIFIED",
            "water_validation": "NONE IDENTIFIED",
            "curved_equilibrium_BNNT_validation": "NOT ESTABLISHED",
            "four_atom_bridge_validation": "NONE",
            "main_advantage": (
                "Reactive B-N, B-H and N-H chemistry is represented."
            ),
            "main_limitation": (
                "Training and validation domain does not match equilibrium solvated R2."
            ),
            "current_risk": "HIGH",
            "current_decision": "REFERENCE_CANDIDATE_NOT_AUTHORIZED",
        },
        {
            "route": "FUNCTIONALIZED_HBN_FIXED_TOPOLOGY",
            "model_type": "Nonreactive bonded fixed-topology model",
            "primary_demonstrated_domain": (
                "Potentially equilibrium h-BN and functionalized aqueous interfaces"
            ),
            "elemental_scope_for_R2": "SOURCE AUDIT PENDING",
            "direct_R2_environment_validation": "PENDING",
            "water_validation": "PENDING PRIMARY-SOURCE AUDIT",
            "curved_equilibrium_BNNT_validation": "PENDING",
            "four_atom_bridge_validation": "UNLIKELY WITHOUT CUSTOM TERMS",
            "main_advantage": (
                "Closer conceptual match to equilibrium solvated structural MD."
            ),
            "main_limitation": (
                "Primary parameter artifacts and exact environment coverage not yet audited."
            ),
            "current_risk": "UNRESOLVED",
            "current_decision": "PRIORITY_SOURCE_AUDIT_REQUIRED",
        },
        {
            "route": "CUSTOM_QM_REFERENCED_FIXED_TOPOLOGY",
            "model_type": (
                "Fixed-topology force field with dedicated QM reference fragments"
            ),
            "primary_demonstrated_domain": (
                "Can be constructed specifically for R2 local environments"
            ),
            "elemental_scope_for_R2": "B/N/H TARGETED",
            "direct_R2_environment_validation": (
                "POSSIBLE AFTER QM DATA GENERATION"
            ),
            "water_validation": (
                "REQUIRES NONBONDED AND INTERFACE VALIDATION"
            ),
            "curved_equilibrium_BNNT_validation": (
                "CAN BE INCLUDED IN QM REFERENCE SET"
            ),
            "four_atom_bridge_validation": (
                "CAN BE DIRECTLY TARGETED"
            ),
            "main_advantage": (
                "Highest domain specificity and explicit validation path."
            ),
            "main_limitation": (
                "Highest development cost; QM calculations are not yet authorized."
            ),
            "current_risk": "MEDIUM_AFTER_VALIDATION",
            "current_decision": "CONTINGENCY_ROUTE_NOT_YET_AUTHORIZED",
        },
    ]

    write_csv(
        ROUTE_MATRIX,
        route_rows,
    )

    gates = {
        "preflight_decision_is_accepted": (
            preflight_summary.get(
                "decision"
            )
            == EXPECTED_PREFLIGHT_DECISION
        ),
        "all_40_environments_classified": (
            len(environment_rows) == 40
        ),
        "all_148_bonded_term_classes_classified": (
            len(term_rows) == 148
        ),
        "all_2112_improper_centers_grouped": (
            sum(
                as_int(row["center_count"])
                for row in improper_rows
            )
            == 2112
        ),
        "all_7_QM_fragment_classes_preserved": (
            len(qm_fragment_rows) == 7
        ),
        "three_routes_compared": (
            len(route_rows) == 3
        ),
        "lele_direct_R2_validation_not_inferred": (
            lele_directly_validated_environments == 0
        ),
        "fixed_topology_validation_not_inferred": (
            fixed_directly_validated_environments == 0
        ),
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
        "AUDIT_FUNCTIONALIZED_HBN_FIXED_TOPOLOGY_PRIMARY_SOURCES_AND_PARAMETERS"
        if not failed_gates
        else
        "REVIEW_R2_ROUTE_COVERAGE_AND_RISK_CLASSIFICATION"
    )

    summary = {
        "decision": decision,
        "environment_classes": len(
            environment_rows
        ),
        "high_risk_environment_classes": high_environment_count,
        "medium_high_risk_environment_classes": (
            medium_high_environment_count
        ),
        "medium_risk_environment_classes": medium_environment_count,
        "bonded_term_classes": len(
            term_rows
        ),
        "high_risk_bonded_term_classes": high_term_count,
        "medium_high_risk_bonded_term_classes": (
            medium_high_term_count
        ),
        "medium_risk_bonded_term_classes": medium_term_count,
        "improper_environment_groups": len(
            improper_rows
        ),
        "improper_centers": sum(
            as_int(row["center_count"])
            for row in improper_rows
        ),
        "critical_centers": sum(
            critical_counts.values()
        ),
        "QM_fragment_classes": len(
            qm_fragment_rows
        ),
        "routes_compared": len(
            route_rows
        ),
        "Lele_directly_validated_R2_environments": (
            lele_directly_validated_environments
        ),
        "fixed_topology_directly_validated_R2_environments": (
            fixed_directly_validated_environments
        ),
        "force_field_route_selected": False,
        "parameter_adoption_authorized": False,
        "topology_generation_authorized": False,
        "formal_charge_assignment_authorized": False,
        "force_field_parameterization_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "QM_calculation_authorized": False,
        "failed_gates": " | ".join(
            failed_gates
        ),
        "required_next_step": required_next_step,
    }

    write_csv(
        SUMMARY,
        [summary],
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

    JSON_OUT.write_text(
        json.dumps(
            {
                "summary": summary,
                "gates": gates,
                "route_comparison": route_rows,
                "risk_counts": {
                    "environment": {
                        "HIGH": high_environment_count,
                        "MEDIUM_HIGH": medium_high_environment_count,
                        "MEDIUM": medium_environment_count,
                    },
                    "bonded_terms": {
                        "HIGH": high_term_count,
                        "MEDIUM_HIGH": medium_high_term_count,
                        "MEDIUM": medium_term_count,
                    },
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    REPORT.write_text(
        f"""# R2 Force-Field Route Coverage and Risk Matrix — Day025

## Classified system

- Chemical environment classes: **{len(environment_rows)}**
- Bonded term classes: **{len(term_rows)}**
- Improper/planarity centers: **{sum(as_int(row['center_count']) for row in improper_rows)}**
- Parameterization-critical centers: **{sum(critical_counts.values())}**
- QM fragment classes preserved: **{len(qm_fragment_rows)}**

## Environment risk distribution

- High: **{high_environment_count}**
- Medium-high: **{medium_high_environment_count}**
- Medium: **{medium_environment_count}**

## Bonded-term risk distribution

- High: **{high_term_count}**
- Medium-high: **{medium_high_term_count}**
- Medium: **{medium_term_count}**

## Route findings

### Lele 2022 ReaxFF

B/N/H reactive records are present. No R2 environment is classified
as directly validated. The demonstrated domain remains gas-phase
reactive chemistry and high-temperature BN nanostructure formation.

### Functionalized h-BN fixed topology

This is conceptually closer to equilibrium solvated structural MD,
but its primary sources, parameter artifacts and exact R2 coverage
must still be audited.

### Custom QM-referenced fixed topology

This route can directly target the novel annulus, edge and bridge
environments, but requires dedicated QM reference data and remains
unauthorized.

## Decision

- Decision: **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Force-field route selected: **NO**
- Parameter adoption authorized: **NO**
- Topology generation authorized: **NO**
- Minimization authorized: **NO**
- MD authorized: **NO**
- QM calculations authorized: **NO**
- Required next step:
  `{required_next_step}`
""",
        encoding="utf-8",
    )

    print(
        "Day025 R2 force-field route coverage "
        "and risk matrix completed."
    )

    print(
        "Environment classes high/medium-high/medium: "
        f"{high_environment_count}/"
        f"{medium_high_environment_count}/"
        f"{medium_environment_count}"
    )

    print(
        "Bonded term classes total/high/medium-high/medium: "
        f"{len(term_rows)}/"
        f"{high_term_count}/"
        f"{medium_high_term_count}/"
        f"{medium_term_count}"
    )

    print(
        "Improper groups / centers: "
        f"{len(improper_rows)}/"
        f"{sum(as_int(row['center_count']) for row in improper_rows)}"
    )

    print(
        "Critical centers / QM fragment classes / routes: "
        f"{sum(critical_counts.values())}/"
        f"{len(qm_fragment_rows)}/"
        f"{len(route_rows)}"
    )

    print(
        "Lele directly validated R2 environments: "
        f"{lele_directly_validated_environments}"
    )

    print(
        "Fixed-topology directly validated R2 environments: "
        f"{fixed_directly_validated_environments}"
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
        "Force-field route selected: NO"
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
        ENV_MATRIX,
        TERM_MATRIX,
        IMPROPER_MATRIX,
        ROUTE_MATRIX,
        QM_MATRIX,
        CRITICAL_REGION_SUMMARY,
        SUMMARY,
        GATES,
        JSON_OUT,
        REPORT,
    ):
        print(
            f"Wrote: {path.relative_to(ROOT)}"
        )

    if failed_gates:
        raise RuntimeError(
            "R2 route coverage and risk matrix requires review."
        )


if __name__ == "__main__":
    main()
