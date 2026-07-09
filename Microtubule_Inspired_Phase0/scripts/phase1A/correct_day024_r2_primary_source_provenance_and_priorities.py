#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/phase1A/day024_chemical_end_rim_design"

GATE3R3 = (
    BASE
    / "33_r2_primary_source_force_field_initial_audit"
)

SOURCE_LEDGER_IN = (
    GATE3R3
    / "r2_curated_primary_source_ledger.csv"
)

COVERAGE_IN = (
    GATE3R3
    / "r2_candidate_model_domain_coverage.csv"
)

ARTIFACTS_IN = (
    GATE3R3
    / "r2_required_source_artifact_retrieval.csv"
)

DECISIONS_IN = (
    GATE3R3
    / "r2_preliminary_model_decisions.csv"
)

SUMMARY_IN = (
    GATE3R3
    / "r2_primary_source_initial_audit_summary.csv"
)

OUT = (
    BASE
    / "34_r2_primary_source_provenance_correction"
)

SOURCE_LEDGER_OUT = (
    OUT
    / "r2_corrected_primary_source_ledger.csv"
)

COVERAGE_OUT = (
    OUT
    / "r2_corrected_candidate_model_domain_coverage.csv"
)

ARTIFACT_PRIORITY_OUT = (
    OUT
    / "r2_prioritized_source_artifact_retrieval.csv"
)

CORRECTIONS_OUT = (
    OUT
    / "r2_primary_source_corrections.csv"
)

DECISIONS_OUT = (
    OUT
    / "r2_corrected_preliminary_model_decisions.csv"
)

SUMMARY_OUT = (
    OUT
    / "r2_primary_source_provenance_correction_summary.csv"
)

JSON_OUT = (
    OUT
    / "r2_primary_source_provenance_correction.json"
)

REPORT = (
    OUT
    / "R2_PRIMARY_SOURCE_PROVENANCE_CORRECTION_DAY024.md"
)

EXPECTED_INITIAL_DECISION = (
    "R2_PRIMARY_SOURCE_FORCE_FIELD_INITIAL_AUDIT_COMPLETED"
)

INCORRECT_DOI = "10.1063/1.1999628"
CORRECT_REAXFF_DOI = "10.1021/acs.jpca.1c09648"

PRIORITY_BY_DOI = {
    "10.1063/5.0242541": 1,
    "10.1021/acs.jpca.1c09648": 2,
    "10.1021/acsanm.2c05285": 3,
    "10.1021/acs.jpclett.7b03443": 4,
    "10.1063/1.4947094": 5,
    "10.1021/acs.langmuir.3c01612": 6,
    "10.1021/acsnano.3c09811": 7,
    "10.1021/acs.jpcc.0c05831": 8,
    "10.1063/1.3626065": 9,
    "10.1063/1.1999628": 10,
}


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


def truthy(value: Any) -> bool:
    return str(value).strip().lower() in {
        "true",
        "yes",
        "1",
    }


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        SOURCE_LEDGER_IN,
        COVERAGE_IN,
        ARTIFACTS_IN,
        DECISIONS_IN,
        SUMMARY_IN,
    ):
        require_file(required)

    initial_summary = read_one(
        SUMMARY_IN
    )

    if initial_summary.get(
        "decision"
    ) != EXPECTED_INITIAL_DECISION:
        raise RuntimeError(
            "Gate 3R.3 initial audit is not accepted."
        )

    source_rows = read_rows(
        SOURCE_LEDGER_IN
    )

    coverage_rows = read_rows(
        COVERAGE_IN
    )

    artifact_rows = read_rows(
        ARTIFACTS_IN
    )

    decision_rows = read_rows(
        DECISIONS_IN
    )

    source_by_id = {
        row["source_id"]: row
        for row in source_rows
    }

    incorrect_sources = [
        row
        for row in source_rows
        if row.get("DOI") == INCORRECT_DOI
    ]

    if len(incorrect_sources) != 1:
        raise RuntimeError(
            "Expected exactly one incorrectly classified 2005 source; "
            f"found {len(incorrect_sources)}."
        )

    incorrect_source = incorrect_sources[0]
    incorrect_source_id = incorrect_source["source_id"]

    incorrect_source.update(
        {
            "model_name": "HYDROGEN_BNNT_THEORETICAL_STUDY",
            "model_family": "QUANTUM_THEORETICAL_REFERENCE",
            "primary_training_or_reference": (
                "Theoretical hydrogen–BNNT interaction study; "
                "not an identified ReaxFF parameter-development source"
            ),
            "explicit_B": True,
            "explicit_N": True,
            "explicit_H": True,
            "reactive": False,
            "polarizable": False,
            "curved_BNNT_validation": "THEORETICAL_REFERENCE_ONLY",
            "water_validation": "NO",
            "source_status": (
                "PRIMARY_REFERENCE_RECLASSIFIED_NOT_FORCE_FIELD"
            ),
        }
    )

    reaxff_matches = [
        row
        for row in source_rows
        if row.get("DOI") == CORRECT_REAXFF_DOI
    ]

    if reaxff_matches:
        raise RuntimeError(
            "Correct 2022 ReaxFF source already exists unexpectedly."
        )

    next_source_number = max(
        int(
            row["source_id"].split("_")[-1]
        )
        for row in source_rows
    ) + 1

    reaxff_source_id = (
        f"SRC_{next_source_number:04d}"
    )

    source_rows.append(
        {
            "source_id": reaxff_source_id,
            "year": 2022,
            "title": (
                "ReaxFF Force Field Development for Gas-Phase "
                "hBN Nanostructure Synthesis"
            ),
            "journal": (
                "The Journal of Physical Chemistry A"
            ),
            "DOI": CORRECT_REAXFF_DOI,
            "model_name": "LELE_GAS_PHASE_HBN_REAXFF",
            "model_family": "REACTIVE_REAXFF_BNH",
            "primary_training_or_reference": (
                "Quantum reference data for gas-phase B/N/H "
                "chemistry and hBN nanostructure synthesis"
            ),
            "explicit_B": True,
            "explicit_N": True,
            "explicit_H": True,
            "reactive": True,
            "polarizable": "QEQ_CHARGE_EQUILIBRATION",
            "curved_BNNT_validation": "NOT_ESTABLISHED",
            "water_validation": "NO",
            "source_status": (
                "PRIMARY_SOURCE_IDENTIFIED_PARAMETER_FILE_REPORTED_IN_SI"
            ),
        }
    )

    write_rows(
        SOURCE_LEDGER_OUT,
        source_rows,
    )

    corrected_coverage_rows = [
        row
        for row in coverage_rows
        if row["source_id"] != incorrect_source_id
    ]

    reclassified_domains = {
        "B_N_BONDED_BULK": "THEORETICAL_REFERENCE_ONLY",
        "B_N_BONDED_CURVED_BNNT": "THEORETICAL_REFERENCE_ONLY",
        "B_H_BONDED": "THEORETICAL_REFERENCE_ONLY",
        "N_H_BONDED": "THEORETICAL_REFERENCE_ONLY",
        "EDGE_ANGLES_AND_TORSIONS": "NOT_A_FORCE_FIELD_SOURCE",
        "FOUR_ATOM_BRIDGE": "NOT_COVERED",
        "PLANARITY_OR_IMPROPERS": "NOT_A_FORCE_FIELD_SOURCE",
        "HBN_WATER_NONBONDED": "NOT_COVERED",
        "CURVATURE_DEPENDENT_WATER_INTERACTION": "NOT_COVERED",
        "POLARIZATION": "NOT_COVERED",
        "AQUEOUS_STABILITY": "NOT_COVERED",
    }

    for domain, status in reclassified_domains.items():
        corrected_coverage_rows.append(
            {
                "source_id": incorrect_source_id,
                "model_name": (
                    "HYDROGEN_BNNT_THEORETICAL_STUDY"
                ),
                "domain": domain,
                "preliminary_coverage": status,
                "evidence_level": (
                    "PRIMARY_ARTICLE_CLASSIFICATION_CORRECTED"
                ),
                "supplementary_parameter_audit_completed": False,
                "transferability_to_R2_established": False,
            }
        )

    reaxff_domains = {
        "B_N_BONDED_BULK": "REACTIVE_CANDIDATE",
        "B_N_BONDED_CURVED_BNNT": "NOT_VALIDATED",
        "B_H_BONDED": "REACTIVE_CANDIDATE",
        "N_H_BONDED": "REACTIVE_CANDIDATE",
        "EDGE_ANGLES_AND_TORSIONS": "REACTIVE_CANDIDATE",
        "FOUR_ATOM_BRIDGE": "REQUIRES_TARGETED_VALIDATION",
        "PLANARITY_OR_IMPROPERS": "EMERGENT_REACTIVE_PES",
        "HBN_WATER_NONBONDED": "NOT_COVERED",
        "CURVATURE_DEPENDENT_WATER_INTERACTION": "NOT_COVERED",
        "POLARIZATION": "QEQ_ONLY_NOT_EXPLICIT_POLARIZATION",
        "AQUEOUS_STABILITY": "NOT_VALIDATED",
    }

    for domain, status in reaxff_domains.items():
        corrected_coverage_rows.append(
            {
                "source_id": reaxff_source_id,
                "model_name": (
                    "LELE_GAS_PHASE_HBN_REAXFF"
                ),
                "domain": domain,
                "preliminary_coverage": status,
                "evidence_level": (
                    "PRIMARY_ARTICLE_AND_SI_PARAMETER_FILE_REPORTED"
                ),
                "supplementary_parameter_audit_completed": False,
                "transferability_to_R2_established": False,
            }
        )

    write_rows(
        COVERAGE_OUT,
        corrected_coverage_rows,
    )

    corrected_decision_rows = []

    for row in decision_rows:
        if row.get(
            "model_name"
        ) == "REAXFF_HBN":
            corrected_decision_rows.append(
                {
                    "model_name": (
                        "HYDROGEN_BNNT_THEORETICAL_STUDY"
                    ),
                    "preliminary_decision": (
                        "RETAIN_AS_QM_REFERENCE_NOT_FORCE_FIELD"
                    ),
                    "adoption_authorized": False,
                    "reason": (
                        "The 2005 paper is not accepted as a "
                        "ReaxFF development or parameter source."
                    ),
                }
            )
        else:
            corrected_decision_rows.append(
                row
            )

    corrected_decision_rows.append(
        {
            "model_name": (
                "LELE_GAS_PHASE_HBN_REAXFF"
            ),
            "preliminary_decision": (
                "RETAIN_FOR_REACTIVE_BNH_PARAMETER_FILE_AUDIT"
            ),
            "adoption_authorized": False,
            "reason": (
                "Includes B/N/H reactive chemistry and reports "
                "an SI parameter file, but was developed for "
                "gas-phase hBN synthesis rather than aqueous R2."
            ),
        }
    )

    write_rows(
        DECISIONS_OUT,
        corrected_decision_rows,
    )

    source_by_id = {
        row["source_id"]: row
        for row in source_rows
    }

    corrected_artifacts = [
        row
        for row in artifact_rows
        if row["source_id"] != incorrect_source_id
    ]

    for artifact_type in (
        "FULL_ARTICLE",
        "SUPPORTING_INFORMATION",
        "PARAMETER_FILE_OR_REPOSITORY",
    ):
        corrected_artifacts.append(
            {
                "source_id": reaxff_source_id,
                "DOI": CORRECT_REAXFF_DOI,
                "artifact_type": artifact_type,
                "retrieval_status": "NOT_YET_RETRIEVED",
                "local_file": "",
                "sha256": "",
                "license_or_access_notes": "",
                "required_for_final_coverage_decision": True,
            }
        )

    for artifact_type in (
        "FULL_ARTICLE",
        "SUPPORTING_INFORMATION",
    ):
        corrected_artifacts.append(
            {
                "source_id": incorrect_source_id,
                "DOI": INCORRECT_DOI,
                "artifact_type": artifact_type,
                "retrieval_status": "NOT_YET_RETRIEVED",
                "local_file": "",
                "sha256": "",
                "license_or_access_notes": (
                    "Reference-only source; not a force-field parameter source"
                ),
                "required_for_final_coverage_decision": False,
            }
        )

    artifact_priority_rows = []

    for row in corrected_artifacts:
        source = source_by_id[
            row["source_id"]
        ]

        DOI = row["DOI"]

        base_priority = PRIORITY_BY_DOI.get(
            DOI,
            99,
        )

        artifact_offset = {
            "PARAMETER_FILE_OR_REPOSITORY": 0,
            "SUPPORTING_INFORMATION": 1,
            "FULL_ARTICLE": 2,
        }.get(
            row["artifact_type"],
            3,
        )

        artifact_priority_rows.append(
            {
                **row,
                "model_name": source[
                    "model_name"
                ],
                "retrieval_priority": (
                    base_priority
                    * 10
                    + artifact_offset
                ),
                "priority_group": (
                    "P1_IMMEDIATE"
                    if base_priority <= 4
                    else (
                        "P2_INTERFACE_REFERENCE"
                        if base_priority <= 7
                        else "P3_SECONDARY_REFERENCE"
                    )
                ),
                "automatic_parameter_adoption_authorized": False,
            }
        )

    artifact_priority_rows.sort(
        key=lambda row: (
            int(
                row[
                    "retrieval_priority"
                ]
            ),
            row[
                "source_id"
            ],
        )
    )

    write_rows(
        ARTIFACT_PRIORITY_OUT,
        artifact_priority_rows,
    )

    correction_rows = [
        {
            "correction_id": "CORR_001",
            "affected_DOI": INCORRECT_DOI,
            "original_classification": (
                "REAXFF_HBN reactive force-field candidate"
            ),
            "corrected_classification": (
                "Theoretical hydrogen–BNNT reference, "
                "not an accepted force-field development source"
            ),
            "reason": (
                "Primary-source verification did not support "
                "the ReaxFF attribution."
            ),
        },
        {
            "correction_id": "CORR_002",
            "affected_DOI": CORRECT_REAXFF_DOI,
            "original_classification": "ABSENT",
            "corrected_classification": (
                "B/N/H ReaxFF candidate for gas-phase "
                "hBN nanostructure synthesis"
            ),
            "reason": (
                "Primary article explicitly reports ReaxFF "
                "development and an SI parameter file."
            ),
        },
    ]

    write_rows(
        CORRECTIONS_OUT,
        correction_rows,
    )

    immediate_artifacts = [
        row
        for row in artifact_priority_rows
        if row[
            "priority_group"
        ] == "P1_IMMEDIATE"
    ]

    gates = {
        "incorrect_2005_ReaxFF_attribution_removed": (
            all(
                row.get(
                    "model_name"
                )
                != "REAXFF_HBN"
                for row in source_rows
            )
        ),
        "2022_BNH_ReaxFF_source_added": (
            any(
                row.get(
                    "DOI"
                )
                == CORRECT_REAXFF_DOI
                for row in source_rows
            )
        ),
        "2025_functionalized_hBN_is_priority1": (
            any(
                row.get(
                    "DOI"
                )
                == "10.1063/5.0242541"
                and row[
                    "priority_group"
                ]
                == "P1_IMMEDIATE"
                for row in artifact_priority_rows
            )
        ),
        "all_sources_remain_unadopted": all(
            not truthy(
                row.get(
                    "adoption_authorized",
                    False,
                )
            )
            for row in corrected_decision_rows
        ),
        "all_artifacts_prohibit_automatic_parameter_adoption": all(
            not truthy(
                row[
                    "automatic_parameter_adoption_authorized"
                ]
            )
            for row in artifact_priority_rows
        ),
        "no_topology_parameters_MD_or_QM_generated": True,
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
        "R2_PRIMARY_SOURCE_PROVENANCE_CORRECTED_"
        "ARTIFACT_RETRIEVAL_PRIORITIZED"
        if accepted
        else
        "R2_PRIMARY_SOURCE_PROVENANCE_CORRECTION_REQUIRES_REVIEW"
    )

    required_next_step = (
        "RETRIEVE_P1_PRIMARY_ARTICLES_SUPPORTING_INFORMATION_"
        "AND_PARAMETER_FILES"
        if accepted
        else
        "REVIEW_PRIMARY_SOURCE_PROVENANCE_CORRECTION"
    )

    summary = {
        "decision": decision,
        "corrected_source_count": len(
            source_rows
        ),
        "source_corrections": len(
            correction_rows
        ),
        "coverage_rows": len(
            corrected_coverage_rows
        ),
        "prioritized_artifacts": len(
            artifact_priority_rows
        ),
        "priority1_artifacts": len(
            immediate_artifacts
        ),
        "models_authorized_for_adoption": 0,
        "force_field_coverage_established": False,
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
        SUMMARY_OUT,
        [
            summary
        ],
    )

    JSON_OUT.write_text(
        json.dumps(
            {
                "summary": summary,
                "corrections": correction_rows,
                "gates": gates,
                "priority1_artifacts": immediate_artifacts,
                "restrictions": [
                    (
                        "No source is adopted solely because it "
                        "contains B, N and H."
                    ),
                    (
                        "The gas-phase hBN ReaxFF must not be "
                        "treated as validated for water or R2."
                    ),
                    (
                        "No topology, charges, parameters, "
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

    priority_lines = "\n".join(
        (
            f"- Priority {row['retrieval_priority']}: "
            f"{row['model_name']} / "
            f"{row['artifact_type']} / "
            f"DOI {row['DOI']}"
        )
        for row in immediate_artifacts
    )

    REPORT.write_text(
        f"""# R2 Primary-Source Provenance Correction

## Corrections

1. DOI `{INCORRECT_DOI}` is retained only as a theoretical
   hydrogen–BNNT reference. It is not treated as an identified
   ReaxFF parameter-development source.
2. DOI `{CORRECT_REAXFF_DOI}` is added as the B/N/H ReaxFF
   candidate for gas-phase hBN nanostructure synthesis.

## Immediate retrieval priorities

{priority_lines}

## Restrictions

- Force-field coverage established: **NO**
- Model adoption authorized: **NO**
- Topology generation authorized: **NO**
- Charge assignment authorized: **NO**
- Parameterization authorized: **NO**
- Energy minimization authorized: **NO**
- MD authorized: **NO**
- QM calculation authorized: **NO**

## Decision

- Decision: **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Required next step:
  `{required_next_step}`
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 primary-source provenance "
        "correction completed."
    )

    print(
        "Corrected sources / corrections: "
        f"{len(source_rows)}/"
        f"{len(correction_rows)}"
    )

    print(
        "Coverage rows / prioritized artifacts / P1 artifacts: "
        f"{len(corrected_coverage_rows)}/"
        f"{len(artifact_priority_rows)}/"
        f"{len(immediate_artifacts)}"
    )

    print(
        "2005 hydrogen-BNNT source classified as force field: NO"
    )

    print(
        "2022 gas-phase hBN B/N/H ReaxFF source added: YES"
    )

    print(
        "Models authorized for adoption: 0"
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
        SOURCE_LEDGER_OUT,
        COVERAGE_OUT,
        ARTIFACT_PRIORITY_OUT,
        CORRECTIONS_OUT,
        DECISIONS_OUT,
        SUMMARY_OUT,
        JSON_OUT,
        REPORT,
    ):
        print(
            f"Wrote: {path.relative_to(ROOT)}"
        )

    if not accepted:
        raise RuntimeError(
            "Primary-source provenance correction requires review."
        )


if __name__ == "__main__":
    main()
