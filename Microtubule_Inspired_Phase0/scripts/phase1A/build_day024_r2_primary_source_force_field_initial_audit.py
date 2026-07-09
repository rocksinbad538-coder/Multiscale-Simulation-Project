#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/phase1A/day024_chemical_end_rim_design"

GATE3R2 = (
    BASE
    / "32_r2_primary_source_force_field_literature_audit_preparation"
)

PREPARATION_SUMMARY = (
    GATE3R2
    / "r2_literature_audit_preparation_summary.csv"
)

REQUIREMENT_MATRIX = (
    GATE3R2
    / "r2_force_field_requirement_coverage_matrix.csv"
)

QM_MATRIX = (
    GATE3R2
    / "r2_qm_fragment_decision_matrix.csv"
)

OUT = (
    BASE
    / "33_r2_primary_source_force_field_initial_audit"
)

SOURCE_LEDGER = (
    OUT
    / "r2_curated_primary_source_ledger.csv"
)

MODEL_COVERAGE = (
    OUT
    / "r2_candidate_model_domain_coverage.csv"
)

MODEL_GAPS = (
    OUT
    / "r2_candidate_model_gap_analysis.csv"
)

REQUIRED_ARTIFACTS = (
    OUT
    / "r2_required_source_artifact_retrieval.csv"
)

PRELIMINARY_DECISIONS = (
    OUT
    / "r2_preliminary_model_decisions.csv"
)

SUMMARY = (
    OUT
    / "r2_primary_source_initial_audit_summary.csv"
)

JSON_OUT = (
    OUT
    / "r2_primary_source_initial_audit.json"
)

REPORT = (
    OUT
    / "R2_PRIMARY_SOURCE_FORCE_FIELD_INITIAL_AUDIT_DAY024.md"
)

EXPECTED_PREPARATION_DECISION = (
    "R2_PRIMARY_SOURCE_FORCE_FIELD_LITERATURE_AUDIT_PREPARED"
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
        raise RuntimeError(
            f"No rows available for {path}"
        )

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


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        PREPARATION_SUMMARY,
        REQUIREMENT_MATRIX,
        QM_MATRIX,
    ):
        require_file(required)

    preparation = read_one(
        PREPARATION_SUMMARY
    )

    if preparation.get(
        "decision"
    ) != EXPECTED_PREPARATION_DECISION:
        raise RuntimeError(
            "Literature-audit preparation gate is not accepted."
        )

    requirement_rows = read_rows(
        REQUIREMENT_MATRIX
    )

    qm_rows = read_rows(
        QM_MATRIX
    )

    sources = [
        {
            "source_id": "SRC_0001",
            "year": 2018,
            "title": (
                "Ab Initio Molecular Dynamics and Lattice "
                "Dynamics-Based Force Field for Modeling "
                "Hexagonal Boron Nitride in Mechanical and "
                "Interfacial Applications"
            ),
            "journal": (
                "The Journal of Physical Chemistry Letters"
            ),
            "DOI": "10.1021/acs.jpclett.7b03443",
            "model_name": "RAJAN_HBN_CLASSICAL_FF",
            "model_family": "FIXED_TOPOLOGY_CLASSICAL",
            "primary_training_or_reference": (
                "DFT-based ab initio MD and lattice dynamics"
            ),
            "explicit_B": True,
            "explicit_N": True,
            "explicit_H": False,
            "reactive": False,
            "polarizable": False,
            "curved_BNNT_validation": "UNCONFIRMED",
            "water_validation": "PARTIAL_OR_INTERFACIAL",
            "source_status": "PRIMARY_SOURCE_IDENTIFIED",
        },
        {
            "source_id": "SRC_0002",
            "year": 2016,
            "title": (
                "Hexagonal boron nitride and water "
                "interaction parameters"
            ),
            "journal": "The Journal of Chemical Physics",
            "DOI": "10.1063/1.4947094",
            "model_name": "WU_HBN_WATER",
            "model_family": "NONBONDED_HBN_WATER",
            "primary_training_or_reference": (
                "DMC-validated RPA water-hBN interaction energies"
            ),
            "explicit_B": True,
            "explicit_N": True,
            "explicit_H": (
                "WATER_H_ONLY"
            ),
            "reactive": False,
            "polarizable": False,
            "curved_BNNT_validation": "NO",
            "water_validation": "YES_CONTACT_ANGLE",
            "source_status": "PRIMARY_SOURCE_IDENTIFIED",
        },
        {
            "source_id": "SRC_0003",
            "year": 2005,
            "title": (
                "The theoretical study on interaction of "
                "hydrogen with single-walled boron nitride nanotubes"
            ),
            "journal": "The Journal of Chemical Physics",
            "DOI": "10.1063/1.1999628",
            "model_name": "REAXFF_HBN",
            "model_family": "REACTIVE_REAXFF",
            "primary_training_or_reference": (
                "Quantum data for B-H, B-B and B-N bond "
                "dissociation and B/N/H angular strain"
            ),
            "explicit_B": True,
            "explicit_N": True,
            "explicit_H": True,
            "reactive": True,
            "polarizable": (
                "CHARGE_EQUILIBRATION"
            ),
            "curved_BNNT_validation": "PARTIAL",
            "water_validation": "NOT_ESTABLISHED",
            "source_status": "PRIMARY_SOURCE_IDENTIFIED",
        },
        {
            "source_id": "SRC_0004",
            "year": 2025,
            "title": (
                "Molecular dynamics simulations of "
                "functionalized nanoporous hexagonal boron "
                "nitride membranes: Ab initio force field and "
                "implications for water desalination"
            ),
            "journal": "The Journal of Chemical Physics",
            "DOI": "10.1063/5.0242541",
            "model_name": "GHORAI_FUNCTIONALIZED_HBN",
            "model_family": (
                "DFT_FITTED_FIXED_TOPOLOGY_AQUEOUS_EDGE_FF"
            ),
            "primary_training_or_reference": (
                "DFT/AIMD and potential-energy-surface fitting; "
                "borazine reference for H functionalization"
            ),
            "explicit_B": True,
            "explicit_N": True,
            "explicit_H": True,
            "reactive": False,
            "polarizable": "UNCONFIRMED",
            "curved_BNNT_validation": "NO",
            "water_validation": "YES",
            "source_status": "PRIMARY_SOURCE_IDENTIFIED",
        },
        {
            "source_id": "SRC_0005",
            "year": 2020,
            "title": (
                "Machine Learning Potential for Hexagonal "
                "Boron Nitride Applied to Thermally and "
                "Mechanically Induced Rippling"
            ),
            "journal": (
                "The Journal of Physical Chemistry C"
            ),
            "DOI": "10.1021/acs.jpcc.0c05831",
            "model_name": "HBN_GAP",
            "model_family": "MACHINE_LEARNING_INTERATOMIC",
            "primary_training_or_reference": (
                "DFT configuration training set"
            ),
            "explicit_B": True,
            "explicit_N": True,
            "explicit_H": False,
            "reactive": (
                "WITHIN_BN_TRAINING_DOMAIN"
            ),
            "polarizable": False,
            "curved_BNNT_validation": "YES",
            "water_validation": "NO",
            "source_status": "PRIMARY_SOURCE_IDENTIFIED",
        },
        {
            "source_id": "SRC_0006",
            "year": 2023,
            "title": (
                "Anisotropic Interfacial Force Field for "
                "Interfaces of Water with Hexagonal Boron Nitride"
            ),
            "journal": "Langmuir",
            "DOI": "10.1021/acs.langmuir.3c01612",
            "model_name": "FENG_ANISOTROPIC_HBN_WATER",
            "model_family": (
                "ANISOTROPIC_INTERFACIAL_NONBONDED"
            ),
            "primary_training_or_reference": (
                "SCAN binding curves and sliding energy surfaces"
            ),
            "explicit_B": True,
            "explicit_N": True,
            "explicit_H": (
                "WATER_H_ONLY"
            ),
            "reactive": False,
            "polarizable": False,
            "curved_BNNT_validation": "NO",
            "water_validation": "YES",
            "source_status": "PRIMARY_SOURCE_IDENTIFIED",
        },
        {
            "source_id": "SRC_0007",
            "year": 2024,
            "title": (
                "Water Electric Field Induced Modulation of "
                "the Wetting of Hexagonal Boron Nitride: "
                "Insights from Multiscale Modeling of "
                "Many-Body Polarization"
            ),
            "journal": "ACS Nano",
            "DOI": "10.1021/acsnano.3c09811",
            "model_name": "Luo_DRUDE_HBN_WATER",
            "model_family": (
                "POLARIZABLE_DRUDE_INTERFACIAL"
            ),
            "primary_training_or_reference": (
                "Anisotropic hBN polarizability and "
                "water-hBN binding reference data"
            ),
            "explicit_B": True,
            "explicit_N": True,
            "explicit_H": (
                "WATER_H_ONLY"
            ),
            "reactive": False,
            "polarizable": True,
            "curved_BNNT_validation": "NO",
            "water_validation": "YES",
            "source_status": "PRIMARY_SOURCE_IDENTIFIED",
        },
        {
            "source_id": "SRC_0008",
            "year": 2023,
            "title": (
                "Boron Nitride Nanotubes: Force Field "
                "Parameterization, Epoxy Interactions, and "
                "Comparison with Carbon Nanotubes for "
                "High-Performance Composite Materials"
            ),
            "journal": "ACS Applied Nano Materials",
            "DOI": "10.1021/acsanm.2c05285",
            "model_name": "IFF_R_BNNT",
            "model_family": (
                "INTERFACE_FORCE_FIELD_REACTIVE_EXTENSION"
            ),
            "primary_training_or_reference": (
                "BNNT/interface parameterization and "
                "composite interfacial validation"
            ),
            "explicit_B": True,
            "explicit_N": True,
            "explicit_H": "UNCONFIRMED",
            "reactive": "IFF_R_FORMALISM",
            "polarizable": (
                "VIRTUAL_ELECTRON_FORMALISM"
            ),
            "curved_BNNT_validation": "YES",
            "water_validation": "NO",
            "source_status": "PRIMARY_SOURCE_IDENTIFIED",
        },
        {
            "source_id": "SRC_0009",
            "year": 2011,
            "title": (
                "Deformation behaviors of an armchair "
                "boron-nitride nanotube under axial tensile strains"
            ),
            "journal": "Journal of Applied Physics",
            "DOI": "10.1063/1.3626065",
            "model_name": "BNNT_TERSOFF_FORCE_MATCHED",
            "model_family": "TERSOFF_BOND_ORDER",
            "primary_training_or_reference": (
                "DFT-assisted force matching for BNNT mechanics"
            ),
            "explicit_B": True,
            "explicit_N": True,
            "explicit_H": False,
            "reactive": (
                "BOND_ORDER_BN_ONLY"
            ),
            "polarizable": False,
            "curved_BNNT_validation": "YES",
            "water_validation": "NO",
            "source_status": "PRIMARY_SOURCE_IDENTIFIED",
        },
    ]

    write_rows(
        SOURCE_LEDGER,
        sources,
    )

    domains = [
        "B_N_BONDED_BULK",
        "B_N_BONDED_CURVED_BNNT",
        "B_H_BONDED",
        "N_H_BONDED",
        "EDGE_ANGLES_AND_TORSIONS",
        "FOUR_ATOM_BRIDGE",
        "PLANARITY_OR_IMPROPERS",
        "HBN_WATER_NONBONDED",
        "CURVATURE_DEPENDENT_WATER_INTERACTION",
        "POLARIZATION",
        "AQUEOUS_STABILITY",
    ]

    coverage: dict[str, dict[str, str]] = {
        "SRC_0001": {
            "B_N_BONDED_BULK": "CANDIDATE",
            "B_N_BONDED_CURVED_BNNT": "UNCONFIRMED",
            "PLANARITY_OR_IMPROPERS": "CANDIDATE",
            "HBN_WATER_NONBONDED": "PARTIAL",
        },
        "SRC_0002": {
            "HBN_WATER_NONBONDED": "STRONG_CANDIDATE",
            "AQUEOUS_STABILITY": "CONTACT_ANGLE_VALIDATED",
            "CURVATURE_DEPENDENT_WATER_INTERACTION": "NOT_COVERED",
        },
        "SRC_0003": {
            "B_N_BONDED_BULK": "CANDIDATE",
            "B_N_BONDED_CURVED_BNNT": "CANDIDATE",
            "B_H_BONDED": "STRONG_CANDIDATE",
            "N_H_BONDED": "REQUIRES_PARAMETER_FILE_AUDIT",
            "EDGE_ANGLES_AND_TORSIONS": "CANDIDATE",
            "FOUR_ATOM_BRIDGE": "REQUIRES_TARGETED_VALIDATION",
            "PLANARITY_OR_IMPROPERS": "EMERGENT_REACTIVE_TERM",
        },
        "SRC_0004": {
            "B_H_BONDED": "STRONG_CANDIDATE",
            "N_H_BONDED": "STRONG_CANDIDATE",
            "EDGE_ANGLES_AND_TORSIONS": "STRONG_CANDIDATE",
            "HBN_WATER_NONBONDED": "STRONG_CANDIDATE",
            "AQUEOUS_STABILITY": "VALIDATED_FOR_FLAT_NANOPORES",
            "FOUR_ATOM_BRIDGE": "NOT_VALIDATED",
        },
        "SRC_0005": {
            "B_N_BONDED_BULK": "STRONG_CANDIDATE",
            "B_N_BONDED_CURVED_BNNT": "STRONG_CANDIDATE",
            "PLANARITY_OR_IMPROPERS": "IMPLICIT_ML_PES",
            "B_H_BONDED": "NOT_COVERED",
            "N_H_BONDED": "NOT_COVERED",
            "HBN_WATER_NONBONDED": "NOT_COVERED",
        },
        "SRC_0006": {
            "HBN_WATER_NONBONDED": "STRONG_CANDIDATE",
            "CURVATURE_DEPENDENT_WATER_INTERACTION": "NOT_VALIDATED",
            "POLARIZATION": "NOT_EXPLICIT",
        },
        "SRC_0007": {
            "HBN_WATER_NONBONDED": "STRONG_CANDIDATE",
            "POLARIZATION": "STRONG_CANDIDATE",
            "CURVATURE_DEPENDENT_WATER_INTERACTION": "NOT_VALIDATED",
        },
        "SRC_0008": {
            "B_N_BONDED_CURVED_BNNT": "STRONG_CANDIDATE",
            "B_N_BONDED_BULK": "CANDIDATE",
            "POLARIZATION": "VIRTUAL_ELECTRON_MODEL",
            "B_H_BONDED": "UNCONFIRMED",
            "N_H_BONDED": "UNCONFIRMED",
            "HBN_WATER_NONBONDED": "NOT_VALIDATED",
        },
        "SRC_0009": {
            "B_N_BONDED_CURVED_BNNT": "STRONG_CANDIDATE",
            "B_N_BONDED_BULK": "CANDIDATE",
            "B_H_BONDED": "NOT_COVERED",
            "N_H_BONDED": "NOT_COVERED",
            "HBN_WATER_NONBONDED": "NOT_COVERED",
        },
    }

    coverage_rows = []

    for source in sources:
        source_id = source[
            "source_id"
        ]

        for domain in domains:
            coverage_rows.append(
                {
                    "source_id": source_id,
                    "model_name": source[
                        "model_name"
                    ],
                    "domain": domain,
                    "preliminary_coverage": coverage.get(
                        source_id,
                        {},
                    ).get(
                        domain,
                        "UNASSESSED",
                    ),
                    "evidence_level": (
                        "ABSTRACT_OR_ARTICLE_TEXT_CONFIRMED"
                    ),
                    "supplementary_parameter_audit_completed": False,
                    "transferability_to_R2_established": False,
                }
            )

    write_rows(
        MODEL_COVERAGE,
        coverage_rows,
    )

    gap_rows = [
        {
            "gap_id": "GAP_01",
            "gap": (
                "No single accepted model covers all R2 domains."
            ),
            "severity": "CRITICAL",
            "resolution_required": (
                "Evaluate hybrid, reactive and custom-reference routes "
                "without mixing models prematurely."
            ),
        },
        {
            "gap_id": "GAP_02",
            "gap": (
                "Four-atom alternating B-N-B-N bridges have no "
                "identified direct validation set."
            ),
            "severity": "CRITICAL",
            "resolution_required": (
                "Retain dedicated bridge and attachment QM fragments."
            ),
        },
        {
            "gap_id": "GAP_03",
            "gap": (
                "Curvature dependence of hBN-water parameters is "
                "not established for the selected R2 radius."
            ),
            "severity": "HIGH",
            "resolution_required": (
                "Compare flat-interface and curved-cluster water "
                "interaction reference energies."
            ),
        },
        {
            "gap_id": "GAP_04",
            "gap": (
                "B-H and N-H edge parameters may be environment-specific."
            ),
            "severity": "HIGH",
            "resolution_required": (
                "Audit 2025 functionalized-hBN supplementary data "
                "and preserve separate B-H/N-H fragment classes."
            ),
        },
        {
            "gap_id": "GAP_05",
            "gap": (
                "No source has yet been mapped against all "
                "22/50/76 R2 term classes."
            ),
            "severity": "CRITICAL",
            "resolution_required": (
                "Retrieve complete parameter files and perform "
                "term-by-term coverage mapping."
            ),
        },
        {
            "gap_id": "GAP_06",
            "gap": (
                "Polarizable water-hBN models do not automatically "
                "cover bonded edge chemistry."
            ),
            "severity": "HIGH",
            "resolution_required": (
                "Keep bonded and interfacial polarization coverage "
                "as separate audit dimensions."
            ),
        },
    ]

    write_rows(
        MODEL_GAPS,
        gap_rows,
    )

    artifact_rows = []

    for source in sources:
        source_id = source[
            "source_id"
        ]

        for artifact_type in (
            "FULL_ARTICLE",
            "SUPPORTING_INFORMATION",
            "PARAMETER_FILE_OR_REPOSITORY",
        ):
            artifact_rows.append(
                {
                    "source_id": source_id,
                    "DOI": source[
                        "DOI"
                    ],
                    "artifact_type": artifact_type,
                    "retrieval_status": "NOT_YET_RETRIEVED",
                    "local_file": "",
                    "sha256": "",
                    "license_or_access_notes": "",
                    "required_for_final_coverage_decision": True,
                }
            )

    write_rows(
        REQUIRED_ARTIFACTS,
        artifact_rows,
    )

    decision_rows = [
        {
            "model_name": "RAJAN_HBN_CLASSICAL_FF",
            "preliminary_decision": "RETAIN_FOR_BN_BONDED_AUDIT",
            "adoption_authorized": False,
            "reason": (
                "Relevant classical hBN bonded baseline; "
                "H-edge and R2 bridge coverage absent."
            ),
        },
        {
            "model_name": "WU_HBN_WATER",
            "preliminary_decision": "RETAIN_FOR_WATER_NONBONDED_AUDIT",
            "adoption_authorized": False,
            "reason": (
                "High-level water-hBN reference; curvature "
                "transferability unresolved."
            ),
        },
        {
            "model_name": "REAXFF_HBN",
            "preliminary_decision": "RETAIN_AS_REACTIVE_BNH_CANDIDATE",
            "adoption_authorized": False,
            "reason": (
                "Broadest B/N/H chemistry; parameter-file provenance "
                "and R2/water validation still required."
            ),
        },
        {
            "model_name": "GHORAI_FUNCTIONALIZED_HBN",
            "preliminary_decision": (
                "HIGHEST_PRIORITY_FOR_EDGE_BH_NH_AUDIT"
            ),
            "adoption_authorized": False,
            "reason": (
                "Directly targets H/OH-functionalized hBN "
                "nanopores in aqueous media."
            ),
        },
        {
            "model_name": "HBN_GAP",
            "preliminary_decision": "RETAIN_FOR_CURVED_BN_REFERENCE",
            "adoption_authorized": False,
            "reason": (
                "Strong B/N and nanotube domain, but no H or water."
            ),
        },
        {
            "model_name": "FENG_ANISOTROPIC_HBN_WATER",
            "preliminary_decision": (
                "RETAIN_FOR_ANISOTROPIC_INTERFACE_COMPARISON"
            ),
            "adoption_authorized": False,
            "reason": (
                "Interfacial nonbonded model only."
            ),
        },
        {
            "model_name": "Luo_DRUDE_HBN_WATER",
            "preliminary_decision": (
                "RETAIN_FOR_POLARIZATION_REFERENCE"
            ),
            "adoption_authorized": False,
            "reason": (
                "Relevant to water-induced polarization, "
                "not bonded edge chemistry."
            ),
        },
        {
            "model_name": "IFF_R_BNNT",
            "preliminary_decision": (
                "RETAIN_FOR_CURVED_BNNT_INTERFACE_AUDIT"
            ),
            "adoption_authorized": False,
            "reason": (
                "Curved BNNT model; H edges and water unconfirmed."
            ),
        },
        {
            "model_name": "BNNT_TERSOFF_FORCE_MATCHED",
            "preliminary_decision": (
                "REFERENCE_ONLY_FOR_BNNT_MECHANICS"
            ),
            "adoption_authorized": False,
            "reason": (
                "B/N mechanical potential with no H/water domain."
            ),
        },
    ]

    write_rows(
        PRELIMINARY_DECISIONS,
        decision_rows,
    )

    summary = {
        "decision": (
            "R2_PRIMARY_SOURCE_FORCE_FIELD_INITIAL_AUDIT_COMPLETED"
        ),
        "primary_sources_identified": len(
            sources
        ),
        "model_families_identified": len(
            {
                row[
                    "model_family"
                ]
                for row in sources
            }
        ),
        "coverage_domain_rows": len(
            coverage_rows
        ),
        "critical_or_high_gaps": len(
            gap_rows
        ),
        "required_artifacts": len(
            artifact_rows
        ),
        "required_artifacts_retrieved": 0,
        "models_authorized_for_adoption": sum(
            bool(
                row[
                    "adoption_authorized"
                ]
            )
            for row in decision_rows
        ),
        "single_complete_R2_force_field_identified": False,
        "force_field_coverage_established": False,
        "topology_generation_authorized": False,
        "formal_charge_assignment_authorized": False,
        "force_field_parameterization_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "QM_calculation_authorized": False,
        "required_next_step": (
            "RETRIEVE_AND_AUDIT_PRIMARY_SOURCE_"
            "SUPPORTING_INFORMATION_AND_PARAMETER_FILES"
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
                "sources": sources,
                "gaps": gap_rows,
                "preliminary_decisions": decision_rows,
                "limitations": [
                    (
                        "Coverage classifications are preliminary "
                        "and intentionally conservative."
                    ),
                    (
                        "No source is accepted until its full article, "
                        "supplementary information and parameter files "
                        "have been audited."
                    ),
                    (
                        "No hybrid combination of models is authorized."
                    ),
                    (
                        "No topology, charges, parameterization, "
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

    REPORT.write_text(
        f"""# R2 Primary-Source Force-Field Initial Audit

## Primary evidence identified

- Primary sources: **{len(sources)}**
- Distinct model families:
  **{summary['model_families_identified']}**
- Coverage-domain records:
  **{len(coverage_rows)}**
- Critical/high-priority gaps:
  **{len(gap_rows)}**

## Main finding

No single source has yet demonstrated simultaneous coverage of:

- flexible curved B-N scaffold chemistry;
- B-H and N-H passivation;
- reconstructed annulus and four-atom bridges;
- water interactions;
- polarization;
- the full R2 22/50/76 bonded-term requirement set.

## Highest-priority evidence

1. Ghorai et al. 2025 for functionalized B-H/N-H aqueous edges.
2. Han et al. 2005 ReaxFF for broad B/N/H chemistry.
3. Rajan et al. 2018 and hBN-GAP for B/N mechanics.
4. Wu et al. 2016, Feng et al. 2023 and Luo et al. 2024
   for water interaction and polarization.
5. IFF-R and Tersoff BNNT models for curvature/mechanical references.

## Restrictions

- Force-field coverage established: **NO**
- Model adoption authorized: **NO**
- Hybrid model construction authorized: **NO**
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
        "initial audit completed."
    )

    print(
        "Primary sources / model families: "
        f"{len(sources)}/"
        f"{summary['model_families_identified']}"
    )

    print(
        "Coverage-domain rows / identified gaps: "
        f"{len(coverage_rows)}/"
        f"{len(gap_rows)}"
    )

    print(
        "Required source artifacts / retrieved: "
        f"{len(artifact_rows)}/0"
    )

    print(
        "Models authorized for adoption: 0"
    )

    print(
        "Single complete R2 force field identified: False"
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
        MODEL_COVERAGE,
        MODEL_GAPS,
        REQUIRED_ARTIFACTS,
        PRELIMINARY_DECISIONS,
        SUMMARY,
        JSON_OUT,
        REPORT,
    ):
        print(
            f"Wrote: {path.relative_to(ROOT)}"
        )


if __name__ == "__main__":
    main()
