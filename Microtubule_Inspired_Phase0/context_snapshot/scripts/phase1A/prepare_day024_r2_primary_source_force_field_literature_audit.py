#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/phase1A/day024_chemical_end_rim_design"

GATE3R = BASE / "30_r2_force_field_coverage_preflight"
GATE3R1 = BASE / "31_r2_local_force_field_asset_audit"

BOND_TERMS = GATE3R / "r2_required_bond_terms.csv"
ANGLE_TERMS = GATE3R / "r2_required_angle_terms.csv"
TORSION_TERMS = GATE3R / "r2_required_torsion_terms.csv"
IMPROPER_CENTERS = GATE3R / "r2_required_improper_centers.csv"
QM_CLASSES = GATE3R / "r2_preliminary_qm_fragment_classes.csv"

LOCAL_AUDIT_SUMMARY = (
    GATE3R1
    / "r2_local_force_field_asset_audit_summary.csv"
)

OUT = (
    BASE
    / "32_r2_primary_source_force_field_literature_audit_preparation"
)

SOURCE_LEDGER = OUT / "r2_primary_source_candidate_ledger.csv"
COVERAGE_MATRIX = OUT / "r2_force_field_requirement_coverage_matrix.csv"
SEARCH_MATRIX = OUT / "r2_primary_source_search_matrix.csv"
QM_DECISION_MATRIX = OUT / "r2_qm_fragment_decision_matrix.csv"
SUMMARY = OUT / "r2_literature_audit_preparation_summary.csv"
JSON_OUT = OUT / "r2_literature_audit_preparation.json"
REPORT = OUT / "R2_PRIMARY_SOURCE_FORCE_FIELD_LITERATURE_AUDIT_PREPARATION_DAY024.md"

EXPECTED_LOCAL_OUTCOME = (
    "NO_LOCAL_PHYSICAL_BN_H_FORCE_FIELD_CANDIDATE_IDENTIFIED"
)


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty required file: {path}")


def read_rows(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(f"No rows found in {path}")

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
        raise RuntimeError(f"No rows available for {path}")

    fields: list[str] = []

    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()

        for row in rows:
            writer.writerow(
                {
                    field: row.get(field, "")
                    for field in fields
                }
            )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    for required in (
        BOND_TERMS,
        ANGLE_TERMS,
        TORSION_TERMS,
        IMPROPER_CENTERS,
        QM_CLASSES,
        LOCAL_AUDIT_SUMMARY,
    ):
        require_file(required)

    local_summary = read_one(LOCAL_AUDIT_SUMMARY)

    if local_summary.get("audit_outcome") != EXPECTED_LOCAL_OUTCOME:
        raise RuntimeError(
            "Unexpected local force-field audit outcome: "
            f"{local_summary.get('audit_outcome')}"
        )

    bond_rows = read_rows(BOND_TERMS)
    angle_rows = read_rows(ANGLE_TERMS)
    torsion_rows = read_rows(TORSION_TERMS)
    improper_rows = read_rows(IMPROPER_CENTERS)
    qm_rows = read_rows(QM_CLASSES)

    source_ledger_rows = [
        {
            "source_id": "SRC_0001",
            "citation": "",
            "DOI_or_primary_URL": "",
            "publication_year": "",
            "source_type": "PRIMARY_RESEARCH_ARTICLE_OR_OFFICIAL_PARAMETER_SET",
            "material_system": "",
            "dimensionality": "",
            "flat_or_curved": "",
            "rigid_or_flexible_model": "",
            "water_model_used": "",
            "ion_model_used": "",
            "B_charge": "",
            "N_charge": "",
            "H_charge_types": "",
            "B_N_bonded_terms": "",
            "B_H_bonded_terms": "",
            "N_H_bonded_terms": "",
            "angle_terms": "",
            "torsion_terms": "",
            "improper_planarity_terms": "",
            "nonbonded_form": "",
            "LJ_or_Buckingham_parameters": "",
            "parameter_fitting_method": "",
            "QM_reference_level": "",
            "validation_targets": "",
            "software_format": "",
            "license_or_reuse_status": "",
            "applicability_to_parent_hBN": "UNASSESSED",
            "applicability_to_edges": "UNASSESSED",
            "applicability_to_BH_NH_passivation": "UNASSESSED",
            "applicability_to_four_atom_bridge": "UNASSESSED",
            "overall_status": "UNASSESSED",
            "review_notes": "",
        }
    ]

    write_rows(
        SOURCE_LEDGER,
        source_ledger_rows,
    )

    coverage_rows = []

    for row in bond_rows:
        coverage_rows.append(
            {
                "requirement_class": "BOND",
                "requirement_id": row["bond_term_id"],
                "term_definition": (
                    f"{row['atom_type_1']} -- {row['atom_type_2']}"
                ),
                "occurrence_count": row["count"],
                "source_id": "",
                "coverage_status": "UNASSESSED",
                "exact_or_analogous": "",
                "transferability_status": "UNASSESSED",
                "QM_reference_required": "UNASSESSED",
                "notes": "",
            }
        )

    for row in angle_rows:
        coverage_rows.append(
            {
                "requirement_class": "ANGLE",
                "requirement_id": row["angle_term_id"],
                "term_definition": (
                    f"{row['atom_type_1']} -- "
                    f"{row['center_type']} -- "
                    f"{row['atom_type_3']}"
                ),
                "occurrence_count": row["count"],
                "source_id": "",
                "coverage_status": "UNASSESSED",
                "exact_or_analogous": "",
                "transferability_status": "UNASSESSED",
                "QM_reference_required": "UNASSESSED",
                "notes": "",
            }
        )

    for row in torsion_rows:
        coverage_rows.append(
            {
                "requirement_class": "TORSION",
                "requirement_id": row["torsion_term_id"],
                "term_definition": (
                    f"{row['atom_type_1']} -- "
                    f"{row['atom_type_2']} -- "
                    f"{row['atom_type_3']} -- "
                    f"{row['atom_type_4']}"
                ),
                "occurrence_count": row["count"],
                "source_id": "",
                "coverage_status": "UNASSESSED",
                "exact_or_analogous": "",
                "transferability_status": "UNASSESSED",
                "QM_reference_required": "UNASSESSED",
                "notes": "",
            }
        )

    coverage_rows.append(
        {
            "requirement_class": "IMPROPER_OR_PLANARITY",
            "requirement_id": "IMPROPER_SCOPE",
            "term_definition": (
                "Three-coordinate B/N centers requiring "
                "planarity or pyramidalization control"
            ),
            "occurrence_count": len(improper_rows),
            "source_id": "",
            "coverage_status": "UNASSESSED",
            "exact_or_analogous": "",
            "transferability_status": "UNASSESSED",
            "QM_reference_required": "UNASSESSED",
            "notes": "",
        }
    )

    write_rows(
        COVERAGE_MATRIX,
        coverage_rows,
    )

    search_rows = [
        {
            "search_id": "SEARCH_01",
            "priority": "HIGHEST",
            "topic": "Flexible h-BN force fields",
            "query": (
                "\"hexagonal boron nitride\" flexible force field "
                "molecular dynamics B N bond angle dihedral"
            ),
            "required_evidence": (
                "Primary article plus complete parameter table or repository"
            ),
            "exclusion_criteria": (
                "Rigid sheets, dummy atoms, topology-only models"
            ),
        },
        {
            "search_id": "SEARCH_02",
            "priority": "HIGHEST",
            "topic": "Boron nitride nanotube force fields",
            "query": (
                "\"boron nitride nanotube\" force field molecular dynamics "
                "bond angle torsion parameters"
            ),
            "required_evidence": (
                "Curvature validation and explicit bonded/nonbonded terms"
            ),
            "exclusion_criteria": (
                "Models validated only for flat h-BN without curvature analysis"
            ),
        },
        {
            "search_id": "SEARCH_03",
            "priority": "HIGHEST",
            "topic": "Hydrogen-terminated BN edges",
            "query": (
                "\"hydrogen terminated\" boron nitride edge "
                "B-H N-H force field molecular dynamics"
            ),
            "required_evidence": (
                "Explicit B-H and N-H charges and bonded terms"
            ),
            "exclusion_criteria": (
                "Pure electronic-structure studies without reusable parameters"
            ),
        },
        {
            "search_id": "SEARCH_04",
            "priority": "HIGH",
            "topic": "h-BN water interfaces",
            "query": (
                "\"h-BN water interface\" molecular dynamics "
                "force field partial charges wettability"
            ),
            "required_evidence": (
                "Water compatibility, electrostatics and interfacial validation"
            ),
            "exclusion_criteria": (
                "Contact-angle fits without published parameter definitions"
            ),
        },
        {
            "search_id": "SEARCH_05",
            "priority": "HIGH",
            "topic": "Reactive BN potentials",
            "query": (
                "boron nitrogen hydrogen reactive force field "
                "ReaxFF B N H parameterization"
            ),
            "required_evidence": (
                "B/N/H training set and validation for bonds and edge chemistry"
            ),
            "exclusion_criteria": (
                "Parameter sets without B-H or N-H validation"
            ),
        },
        {
            "search_id": "SEARCH_06",
            "priority": "HIGH",
            "topic": "BN edge and defect parameterization",
            "query": (
                "boron nitride edge defect force field "
                "parameterization B-H N-H"
            ),
            "required_evidence": (
                "Edge, defect or reconstructed-site transferability"
            ),
            "exclusion_criteria": (
                "Bulk-only harmonic models"
            ),
        },
    ]

    write_rows(
        SEARCH_MATRIX,
        search_rows,
    )

    qm_decision_rows = []

    for row in qm_rows:
        qm_decision_rows.append(
            {
                **row,
                "covered_by_primary_force_field_source": "UNASSESSED",
                "coverage_source_id": "",
                "exact_environment_match": "UNASSESSED",
                "transferability_acceptable": "UNASSESSED",
                "QM_fragment_still_required": "UNASSESSED",
                "proposed_QM_method": "",
                "proposed_basis_or_pseudopotential": "",
                "charge_and_multiplicity": "",
                "geometry_protocol": "",
                "scan_coordinates": "",
                "validation_observables": "",
                "final_decision": "UNASSESSED",
            }
        )

    write_rows(
        QM_DECISION_MATRIX,
        qm_decision_rows,
    )

    summary = {
        "decision": (
            "R2_PRIMARY_SOURCE_FORCE_FIELD_LITERATURE_AUDIT_PREPARED"
        ),
        "required_bond_term_classes": len(bond_rows),
        "required_angle_term_classes": len(angle_rows),
        "required_torsion_term_classes": len(torsion_rows),
        "improper_planarity_centers": len(improper_rows),
        "coverage_matrix_rows": len(coverage_rows),
        "primary_source_search_queries": len(search_rows),
        "preliminary_QM_fragment_classes": len(qm_rows),
        "local_physical_force_field_candidates": 0,
        "literature_sources_reviewed": 0,
        "force_field_coverage_established": False,
        "topology_generation_authorized": False,
        "formal_charge_assignment_authorized": False,
        "force_field_parameterization_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "QM_calculation_authorized": False,
        "required_next_step": (
            "EXECUTE_PRIMARY_SOURCE_BN_H_FORCE_FIELD_LITERATURE_AUDIT"
        ),
    }

    write_rows(
        SUMMARY,
        [summary],
    )

    JSON_OUT.write_text(
        json.dumps(
            {
                "summary": summary,
                "search_matrix": search_rows,
                "restrictions": [
                    "No force-field source has yet been accepted.",
                    "No parameters or charges are assigned.",
                    "No topology is generated.",
                    "No minimization, MD or QM calculation is authorized.",
                    (
                        "The literature audit must use primary articles, "
                        "official documentation or original parameter repositories."
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
        f"""# R2 Primary-Source Force-Field Literature Audit Preparation

## Current evidence

- Local physical BN/H force-field candidates: **0**
- Required bond classes: **{len(bond_rows)}**
- Required angle classes: **{len(angle_rows)}**
- Required torsion classes: **{len(torsion_rows)}**
- Improper/planarity centers: **{len(improper_rows)}**
- Preliminary QM fragment classes: **{len(qm_rows)}**

## Generated audit instruments

- Primary-source candidate ledger
- Requirement-by-requirement coverage matrix
- Literature search matrix
- QM fragment decision matrix

## Restrictions

- Force-field coverage established: **NO**
- Topology generation authorized: **NO**
- Charge assignment authorized: **NO**
- Parameterization authorized: **NO**
- Energy minimization authorized: **NO**
- MD authorized: **NO**
- QM calculation authorized: **NO**

## Required next step

`{summary['required_next_step']}`
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 primary-source force-field "
        "literature-audit preparation completed."
    )

    print(
        "Required bond/angle/torsion classes: "
        f"{len(bond_rows)}/"
        f"{len(angle_rows)}/"
        f"{len(torsion_rows)}"
    )

    print(
        "Improper-planarity centers / coverage rows: "
        f"{len(improper_rows)}/"
        f"{len(coverage_rows)}"
    )

    print(
        "Search queries / QM fragment classes: "
        f"{len(search_rows)}/"
        f"{len(qm_rows)}"
    )

    print(
        "Local physical force-field candidates: 0"
    )

    print(
        "Decision: "
        f"{summary['decision']}"
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
        "Required next step: "
        f"{summary['required_next_step']}"
    )

    for path in (
        SOURCE_LEDGER,
        COVERAGE_MATRIX,
        SEARCH_MATRIX,
        QM_DECISION_MATRIX,
        SUMMARY,
        JSON_OUT,
        REPORT,
    ):
        print(
            f"Wrote: {path.relative_to(ROOT)}"
        )


if __name__ == "__main__":
    main()
