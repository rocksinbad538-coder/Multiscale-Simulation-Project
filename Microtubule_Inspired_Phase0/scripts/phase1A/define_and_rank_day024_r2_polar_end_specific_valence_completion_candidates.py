#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

GATE3A_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "01_r2_parent_rim_chemical_audit"
)

OUTPUT_ROOT = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design/"
    "02_r2_polar_end_specific_candidate_ranking"
)

PARENT_SUMMARY = (
    GATE3A_ROOT
    / "r2_parent_rim_chemical_audit_summary.csv"
)

PARENT_GATES = (
    GATE3A_ROOT
    / "r2_parent_rim_chemical_audit_gates.csv"
)

DESIGN_CONSTRAINTS = (
    GATE3A_ROOT
    / "r2_chemical_end_rim_design_constraints.csv"
)

TERMINAL_END_SUMMARY = (
    GATE3A_ROOT
    / "r2_parent_terminal_end_summary.csv"
)

TERMINAL_ATOMS = (
    GATE3A_ROOT
    / "r2_parent_terminal_rim_atoms.csv"
)

TOPOLOGY_CLASSIFICATION = (
    GATE3A_ROOT
    / "topology_terminal_classification"
    / "hbn_terminal_coordination_classification_summary.csv"
)

VALENCE_REQUIREMENTS_CSV = (
    OUTPUT_ROOT
    / "r2_polar_end_valence_completion_requirements.csv"
)

CANDIDATE_MATRIX_CSV = (
    OUTPUT_ROOT
    / "r2_polar_end_specific_candidate_matrix.csv"
)

RANKING_CSV = (
    OUTPUT_ROOT
    / "r2_polar_end_specific_candidate_ranking.csv"
)

SELECTION_SUMMARY_CSV = (
    OUTPUT_ROOT
    / "r2_polar_end_specific_candidate_selection_summary.csv"
)

GATES_CSV = (
    OUTPUT_ROOT
    / "r2_polar_end_specific_candidate_ranking_gates.csv"
)

DESIGN_CONTRACT_JSON = (
    OUTPUT_ROOT
    / "r2_primary_graded_bn_collar_design_contract.json"
)

SOURCE_MANIFEST_CSV = (
    OUTPUT_ROOT
    / "r2_candidate_ranking_source_manifest.csv"
)

REPORT_MD = (
    OUTPUT_ROOT
    / "R2_POLAR_END_SPECIFIC_VALENCE_COMPLETION_RANKING_DAY024.md"
)

EXPECTED_PARENT_DECISION = (
    "R2_PARENT_RIM_AND_CHEMICAL_CONSTRAINTS_AUDITED"
)

EXPECTED_CLASSIFICATION_DECISION = (
    "HBN_EXPLICIT_TOPOLOGY_AND_TERMINAL_COORDINATION_CLASSIFIED"
)

PASS_DECISION = (
    "R2_POLAR_END_SPECIFIC_VALENCE_COMPLETION_CANDIDATES_RANKED"
)

PRIMARY_CANDIDATE = (
    "C2_GRADED_HETEROPOLAR_BN_COLLAR_ANNULUS"
)

CONTINGENCY_CANDIDATE = (
    "C3_RECONSTRUCTED_EDGE_PLUS_GRADED_BN_ANNULUS"
)

EXPECTED_TERMINAL_ATOMS_TOTAL = 60
EXPECTED_TERMINAL_ATOMS_PER_END = 30

EXPECTED_LOWER_B = 30
EXPECTED_LOWER_N = 0
EXPECTED_UPPER_B = 0
EXPECTED_UPPER_N = 30

CURRENT_COORDINATION = 1
TARGET_COORDINATION = 3
MISSING_BONDS_PER_TERMINAL_SITE = 2

EXPECTED_NEW_PARENT_BONDS_PER_END = (
    EXPECTED_TERMINAL_ATOMS_PER_END
    * MISSING_BONDS_PER_TERMINAL_SITE
)

EXPECTED_NEW_PARENT_BONDS_TOTAL = (
    EXPECTED_TERMINAL_ATOMS_TOTAL
    * MISSING_BONDS_PER_TERMINAL_SITE
)

SCORE_WEIGHTS = {
    "valence_completion_score": 0.20,
    "aperture_fidelity_score": 0.18,
    "end_specificity_score": 0.15,
    "neutrality_feasibility_score": 0.12,
    "heteropolar_junction_score": 0.10,
    "rigidity_score": 0.08,
    "minimal_electronic_perturbation_score": 0.07,
    "topological_constructability_score": 0.05,
    "chemical_precedent_score": 0.05,
}

RISK_PENALTIES = {
    "junction_strain_risk": 2.0,
    "charge_uncertainty_risk": 2.0,
    "flexibility_risk": 1.5,
    "forcefield_novelty_risk": 1.5,
}

PARENT_RECONSTRUCTION_PENALTY = 6.0


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
    )


def require_file(path: Path) -> None:
    if (
        not path.exists()
        or not path.is_file()
        or path.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing or empty required file: {path}"
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


def read_csv_rows(
    path: Path,
) -> list[dict[str, str]]:
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


def read_single_csv_row(
    path: Path,
) -> dict[str, str]:
    rows = read_csv_rows(path)

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


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


def parse_int(
    row: dict[str, str],
    key: str,
) -> int:
    try:
        return int(float(row[key]))
    except (
        KeyError,
        TypeError,
        ValueError,
    ) as exc:
        raise RuntimeError(
            f"Could not parse integer field {key!r}"
        ) from exc


def parse_float(
    row: dict[str, str],
    key: str,
) -> float:
    try:
        value = float(row[key])
    except (
        KeyError,
        TypeError,
        ValueError,
    ) as exc:
        raise RuntimeError(
            f"Could not parse numeric field {key!r}"
        ) from exc

    if not math.isfinite(value):
        raise RuntimeError(
            f"Non-finite field {key!r}"
        )

    return value


def weighted_score(
    candidate: dict[str, Any],
) -> tuple[
    float,
    float,
    float,
]:
    base_score = 0.0

    for field, weight in SCORE_WEIGHTS.items():
        score = float(candidate[field])

        if not 0.0 <= score <= 5.0:
            raise RuntimeError(
                f"Invalid score for {candidate['candidate_id']} "
                f"{field}: {score}"
            )

        base_score += (
            weight
            * score
            / 5.0
            * 100.0
        )

    risk_penalty = 0.0

    for field, multiplier in RISK_PENALTIES.items():
        risk = float(candidate[field])

        if not 0.0 <= risk <= 5.0:
            raise RuntimeError(
                f"Invalid risk for {candidate['candidate_id']} "
                f"{field}: {risk}"
            )

        risk_penalty += (
            multiplier
            * risk
        )

    if bool(
        candidate[
            "requires_parent_rim_reconstruction"
        ]
    ):
        risk_penalty += (
            PARENT_RECONSTRUCTION_PENALTY
        )

    final_score = (
        base_score
        - risk_penalty
    )

    return (
        base_score,
        risk_penalty,
        final_score,
    )


def main() -> None:
    OUTPUT_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        PARENT_SUMMARY,
        PARENT_GATES,
        DESIGN_CONSTRAINTS,
        TERMINAL_END_SUMMARY,
        TERMINAL_ATOMS,
        TOPOLOGY_CLASSIFICATION,
    ):
        require_file(required)

    parent = read_single_csv_row(
        PARENT_SUMMARY
    )

    classification = read_single_csv_row(
        TOPOLOGY_CLASSIFICATION
    )

    parent_gate_rows = read_csv_rows(
        PARENT_GATES
    )

    end_rows = read_csv_rows(
        TERMINAL_END_SUMMARY
    )

    terminal_atom_rows = read_csv_rows(
        TERMINAL_ATOMS
    )

    if parent.get(
        "decision"
    ) != EXPECTED_PARENT_DECISION:
        raise RuntimeError(
            "Gate 3A parent audit is not accepted."
        )

    if classification.get(
        "decision"
    ) != EXPECTED_CLASSIFICATION_DECISION:
        raise RuntimeError(
            "Explicit topology classification is not accepted."
        )

    failed_parent_gates = [
        row.get("gate", "")
        for row in parent_gate_rows
        if not parse_bool(
            row.get("pass", "false")
        )
    ]

    if failed_parent_gates:
        raise RuntimeError(
            "Gate 3A still contains failed gates: "
            + " | ".join(failed_parent_gates)
        )

    degree1_atoms = parse_int(
        parent,
        "degree1_terminal_atoms",
    )

    degree2_atoms = parse_int(
        parent,
        "degree2_atoms",
    )

    degree3_atoms = parse_int(
        parent,
        "degree3_interior_atoms",
    )

    lower_atoms = parse_int(
        parent,
        "lower_end_atoms",
    )

    lower_b = parse_int(
        parent,
        "lower_end_B_atoms",
    )

    lower_n = parse_int(
        parent,
        "lower_end_N_atoms",
    )

    upper_atoms = parse_int(
        parent,
        "upper_end_atoms",
    )

    upper_b = parse_int(
        parent,
        "upper_end_B_atoms",
    )

    upper_n = parse_int(
        parent,
        "upper_end_N_atoms",
    )

    aperture_diameter_nm = parse_float(
        parent,
        "target_aperture_diameter_nm",
    )

    aperture_radius_nm = parse_float(
        parent,
        "target_aperture_radius_nm",
    )

    open_area_fraction = parse_float(
        parent,
        "target_open_area_fraction",
    )

    parent_rim_radius_nm = parse_float(
        parent,
        "parent_rim_mean_radius_nm",
    )

    required_radial_occlusion_nm = parse_float(
        parent,
        "required_radial_occlusion_nm",
    )

    annular_area_nm2 = parse_float(
        parent,
        "annular_area_per_end_nm2",
    )

    estimated_bn_atoms_per_end = parse_float(
        parent,
        "estimated_monolayer_hBN_atoms_per_end",
    )

    terminal_deficit = parse_int(
        parent,
        "terminal_coordination_deficit_to_degree3",
    )

    lower_end_row = next(
        (
            row
            for row in end_rows
            if row.get("end") == "LOWER"
        ),
        None,
    )

    upper_end_row = next(
        (
            row
            for row in end_rows
            if row.get("end") == "UPPER"
        ),
        None,
    )

    if (
        lower_end_row is None
        or upper_end_row is None
    ):
        raise RuntimeError(
            "Could not resolve LOWER and UPPER end summaries."
        )

    valence_rows = [
        {
            "end": "LOWER",
            "parent_terminal_element": "B",
            "complementary_first_junction_element": "N",
            "forbidden_primary_same_element_junction": "B-B",
            "terminal_sites": lower_atoms,
            "current_coordination_per_site": CURRENT_COORDINATION,
            "target_coordination_per_site": TARGET_COORDINATION,
            "missing_bond_incidences_per_site": (
                MISSING_BONDS_PER_TERMINAL_SITE
            ),
            "required_new_parent_junction_bonds": (
                lower_atoms
                * MISSING_BONDS_PER_TERMINAL_SITE
            ),
            "design_requirement": (
                "Resolve every B-terminal site using two "
                "additional junction incidences without "
                "assuming B-B primary bonding."
            ),
        },
        {
            "end": "UPPER",
            "parent_terminal_element": "N",
            "complementary_first_junction_element": "B",
            "forbidden_primary_same_element_junction": "N-N",
            "terminal_sites": upper_atoms,
            "current_coordination_per_site": CURRENT_COORDINATION,
            "target_coordination_per_site": TARGET_COORDINATION,
            "missing_bond_incidences_per_site": (
                MISSING_BONDS_PER_TERMINAL_SITE
            ),
            "required_new_parent_junction_bonds": (
                upper_atoms
                * MISSING_BONDS_PER_TERMINAL_SITE
            ),
            "design_requirement": (
                "Resolve every N-terminal site using two "
                "additional junction incidences without "
                "assuming N-N primary bonding."
            ),
        },
    ]

    write_csv(
        VALENCE_REQUIREMENTS_CSV,
        valence_rows,
    )

    candidates: list[
        dict[str, Any]
    ] = [
        {
            "candidate_id": (
                "C0_SMALL_GROUP_PASSIVATION_ONLY"
            ),
            "candidate_name": (
                "Small-group termination without annulus"
            ),
            "lower_end_strategy": (
                "Conceptual two-ligand completion of each "
                "degree-1 B terminal site."
            ),
            "upper_end_strategy": (
                "Conceptual two-ligand completion of each "
                "degree-1 N terminal site."
            ),
            "structural_role": (
                "Valence-reference control only"
            ),
            "parent_bonds_per_terminal_site": 2,
            "expected_new_parent_bonds_total": 120,
            "preserves_target_aperture": False,
            "resolves_end_asymmetry": True,
            "requires_parent_rim_reconstruction": False,
            "valence_completion_score": 4,
            "aperture_fidelity_score": 0,
            "end_specificity_score": 4,
            "neutrality_feasibility_score": 2,
            "heteropolar_junction_score": 1,
            "rigidity_score": 1,
            "minimal_electronic_perturbation_score": 3,
            "topological_constructability_score": 4,
            "chemical_precedent_score": 3,
            "junction_strain_risk": 1,
            "charge_uncertainty_risk": 3,
            "flexibility_risk": 2,
            "forcefield_novelty_risk": 1,
            "status": (
                "REJECT_AS_R2_CAP_REPLACEMENT"
            ),
            "reason": (
                "May terminate dangling valence conceptually "
                "but cannot reproduce the validated annular "
                "steric surface or 0.84 nm pore."
            ),
        },
        {
            "candidate_id": (
                "C1_DIRECT_PLANAR_BN_ANNULUS_DUAL_LINK"
            ),
            "candidate_name": (
                "Direct planar BN annulus attached by two "
                "junction bonds per parent terminal atom"
            ),
            "lower_end_strategy": (
                "B parent sites connect directly to "
                "N-rich annulus-interface atoms."
            ),
            "upper_end_strategy": (
                "N parent sites connect directly to "
                "B-rich annulus-interface atoms."
            ),
            "structural_role": (
                "Direct inorganic annular candidate"
            ),
            "parent_bonds_per_terminal_site": 2,
            "expected_new_parent_bonds_total": 120,
            "preserves_target_aperture": True,
            "resolves_end_asymmetry": True,
            "requires_parent_rim_reconstruction": False,
            "valence_completion_score": 5,
            "aperture_fidelity_score": 5,
            "end_specificity_score": 5,
            "neutrality_feasibility_score": 3,
            "heteropolar_junction_score": 5,
            "rigidity_score": 5,
            "minimal_electronic_perturbation_score": 4,
            "topological_constructability_score": 2,
            "chemical_precedent_score": 1,
            "junction_strain_risk": 5,
            "charge_uncertainty_risk": 2,
            "flexibility_risk": 0,
            "forcefield_novelty_risk": 4,
            "status": (
                "DEFER_DIRECT_SEAM_HIGH_STRAIN"
            ),
            "reason": (
                "Reproduces the steric annulus but assumes "
                "a highly strained abrupt tube-to-plane junction."
            ),
        },
        {
            "candidate_id": PRIMARY_CANDIDATE,
            "candidate_name": (
                "Graded heteropolar BN collar connected "
                "to a planar BN annulus"
            ),
            "lower_end_strategy": (
                "B-terminated parent couples first to an "
                "N-rich transition collar; the collar "
                "progressively joins the annular BN plate."
            ),
            "upper_end_strategy": (
                "N-terminated parent couples first to a "
                "B-rich transition collar; the collar "
                "progressively joins the annular BN plate."
            ),
            "structural_role": (
                "Primary static topological-design hypothesis"
            ),
            "parent_bonds_per_terminal_site": 2,
            "expected_new_parent_bonds_total": 120,
            "preserves_target_aperture": True,
            "resolves_end_asymmetry": True,
            "requires_parent_rim_reconstruction": False,
            "valence_completion_score": 5,
            "aperture_fidelity_score": 5,
            "end_specificity_score": 5,
            "neutrality_feasibility_score": 4,
            "heteropolar_junction_score": 5,
            "rigidity_score": 5,
            "minimal_electronic_perturbation_score": 4,
            "topological_constructability_score": 4,
            "chemical_precedent_score": 2,
            "junction_strain_risk": 3,
            "charge_uncertainty_risk": 2,
            "flexibility_risk": 0,
            "forcefield_novelty_risk": 3,
            "status": (
                "ADVANCE_TO_CONNECTIVITY_BLUEPRINT_GATE"
            ),
            "reason": (
                "Best balance between parent valence "
                "completion, heteropolar bonding, rigidity, "
                "end specificity, and preservation of the "
                "validated annular aperture."
            ),
        },
        {
            "candidate_id": CONTINGENCY_CANDIDATE,
            "candidate_name": (
                "Parent-edge reconstruction followed by "
                "graded BN collar and annulus"
            ),
            "lower_end_strategy": (
                "Reconstruct or remove the degree-1 B row "
                "before building an N-complementary collar."
            ),
            "upper_end_strategy": (
                "Reconstruct or remove the degree-1 N row "
                "before building a B-complementary collar."
            ),
            "structural_role": (
                "Contingency if direct degree-1-site "
                "completion is topologically infeasible"
            ),
            "parent_bonds_per_terminal_site": 0,
            "expected_new_parent_bonds_total": 0,
            "preserves_target_aperture": True,
            "resolves_end_asymmetry": True,
            "requires_parent_rim_reconstruction": True,
            "valence_completion_score": 5,
            "aperture_fidelity_score": 4,
            "end_specificity_score": 5,
            "neutrality_feasibility_score": 4,
            "heteropolar_junction_score": 5,
            "rigidity_score": 5,
            "minimal_electronic_perturbation_score": 3,
            "topological_constructability_score": 3,
            "chemical_precedent_score": 3,
            "junction_strain_risk": 2,
            "charge_uncertainty_risk": 2,
            "flexibility_risk": 0,
            "forcefield_novelty_risk": 3,
            "status": (
                "RETAIN_AS_PARENT_RECONSTRUCTION_CONTINGENCY"
            ),
            "reason": (
                "Potentially cleaner chemistry but changes "
                "the accepted parent scaffold and therefore "
                "requires a separate structural-comparability gate."
            ),
        },
        {
            "candidate_id": (
                "C4_END_SPECIFIC_ORGANIC_MACROCYCLE"
            ),
            "candidate_name": (
                "Rigid organic or hybrid macrocycle with "
                "different B-end and N-end anchors"
            ),
            "lower_end_strategy": (
                "N/O-donor-rich anchors directed toward "
                "the B-terminated parent end."
            ),
            "upper_end_strategy": (
                "Electrophilic or B-containing anchors "
                "directed toward the N-terminated parent end."
            ),
            "structural_role": (
                "Organic fallback candidate"
            ),
            "parent_bonds_per_terminal_site": 2,
            "expected_new_parent_bonds_total": 120,
            "preserves_target_aperture": True,
            "resolves_end_asymmetry": True,
            "requires_parent_rim_reconstruction": False,
            "valence_completion_score": 4,
            "aperture_fidelity_score": 4,
            "end_specificity_score": 5,
            "neutrality_feasibility_score": 2,
            "heteropolar_junction_score": 3,
            "rigidity_score": 4,
            "minimal_electronic_perturbation_score": 1,
            "topological_constructability_score": 2,
            "chemical_precedent_score": 3,
            "junction_strain_risk": 3,
            "charge_uncertainty_risk": 4,
            "flexibility_risk": 1,
            "forcefield_novelty_risk": 4,
            "status": (
                "FALLBACK_AFTER_INORGANIC_CANDIDATE_FAILURE"
            ),
            "reason": (
                "Potentially modular but introduces larger "
                "chemical, electrostatic, and force-field "
                "perturbations near the confined water."
            ),
        },
        {
            "candidate_id": (
                "C5_INWARD_FLEXIBLE_TETHER_CORONA"
            ),
            "candidate_name": (
                "Flexible inward tether corona"
            ),
            "lower_end_strategy": (
                "B-specific tether anchors with inward groups."
            ),
            "upper_end_strategy": (
                "N-specific tether anchors with inward groups."
            ),
            "structural_role": (
                "Flexible confinement concept"
            ),
            "parent_bonds_per_terminal_site": 2,
            "expected_new_parent_bonds_total": 120,
            "preserves_target_aperture": False,
            "resolves_end_asymmetry": True,
            "requires_parent_rim_reconstruction": False,
            "valence_completion_score": 3,
            "aperture_fidelity_score": 2,
            "end_specificity_score": 4,
            "neutrality_feasibility_score": 2,
            "heteropolar_junction_score": 2,
            "rigidity_score": 1,
            "minimal_electronic_perturbation_score": 1,
            "topological_constructability_score": 3,
            "chemical_precedent_score": 3,
            "junction_strain_risk": 2,
            "charge_uncertainty_risk": 3,
            "flexibility_risk": 5,
            "forcefield_novelty_risk": 3,
            "status": (
                "DEFER_PORE_NOT_STRUCTURALLY_FIXED"
            ),
            "reason": (
                "Cannot guarantee a persistent 0.839 nm "
                "effective pore or rigid steric envelope."
            ),
        },
        {
            "candidate_id": (
                "C6_METAL_OR_SPIN_ACTIVE_COORDINATED_CAP"
            ),
            "candidate_name": (
                "Metal-coordinated or spin-active annular cap"
            ),
            "lower_end_strategy": (
                "Coordination network at B-terminated end."
            ),
            "upper_end_strategy": (
                "Coordination network at N-terminated end."
            ),
            "structural_role": (
                "Out-of-scope electronic architecture"
            ),
            "parent_bonds_per_terminal_site": 2,
            "expected_new_parent_bonds_total": 120,
            "preserves_target_aperture": True,
            "resolves_end_asymmetry": True,
            "requires_parent_rim_reconstruction": False,
            "valence_completion_score": 4,
            "aperture_fidelity_score": 4,
            "end_specificity_score": 5,
            "neutrality_feasibility_score": 1,
            "heteropolar_junction_score": 3,
            "rigidity_score": 4,
            "minimal_electronic_perturbation_score": 0,
            "topological_constructability_score": 2,
            "chemical_precedent_score": 2,
            "junction_strain_risk": 3,
            "charge_uncertainty_risk": 5,
            "flexibility_risk": 0,
            "forcefield_novelty_risk": 5,
            "status": (
                "REJECT_AT_CURRENT_GATE"
            ),
            "reason": (
                "Introduces unvalidated charge, spin, and "
                "electronic degrees of freedom before the "
                "confinement architecture is chemically resolved."
            ),
        },
    ]

    ranking_rows = []

    for candidate in candidates:
        (
            base_score,
            risk_penalty,
            final_score,
        ) = weighted_score(candidate)

        candidate[
            "base_weighted_score_0_100"
        ] = base_score

        candidate[
            "risk_penalty_points"
        ] = risk_penalty

        candidate[
            "final_screening_score"
        ] = final_score

        ranking_rows.append(
            dict(candidate)
        )

    ranking_rows.sort(
        key=lambda row: float(
            row[
                "final_screening_score"
            ]
        ),
        reverse=True,
    )

    for rank, row in enumerate(
        ranking_rows,
        start=1,
    ):
        row["rank"] = rank

    write_csv(
        CANDIDATE_MATRIX_CSV,
        candidates,
    )

    write_csv(
        RANKING_CSV,
        ranking_rows,
    )

    primary = next(
        row
        for row in ranking_rows
        if row[
            "candidate_id"
        ] == PRIMARY_CANDIDATE
    )

    contingency = next(
        row
        for row in ranking_rows
        if row[
            "candidate_id"
        ] == CONTINGENCY_CANDIDATE
    )

    highest_ranked = ranking_rows[0]

    gates = {
        "Gate3A_parent_audit_is_accepted": (
            parent.get(
                "decision"
            )
            == EXPECTED_PARENT_DECISION
        ),
        "Gate3A_has_no_failed_gates": (
            len(
                failed_parent_gates
            )
            == 0
        ),
        "explicit_topology_classification_is_accepted": (
            classification.get(
                "decision"
            )
            == EXPECTED_CLASSIFICATION_DECISION
        ),
        "parent_has_60_degree1_terminal_atoms": (
            degree1_atoms
            == EXPECTED_TERMINAL_ATOMS_TOTAL
        ),
        "parent_has_zero_degree2_atoms": (
            degree2_atoms == 0
        ),
        "parent_has_1620_degree3_atoms": (
            degree3_atoms == 1620
        ),
        "lower_end_is_30B_0N": (
            lower_atoms
            == EXPECTED_TERMINAL_ATOMS_PER_END
            and lower_b
            == EXPECTED_LOWER_B
            and lower_n
            == EXPECTED_LOWER_N
        ),
        "upper_end_is_0B_30N": (
            upper_atoms
            == EXPECTED_TERMINAL_ATOMS_PER_END
            and upper_b
            == EXPECTED_UPPER_B
            and upper_n
            == EXPECTED_UPPER_N
        ),
        "terminal_coordination_deficit_is_120": (
            terminal_deficit
            == EXPECTED_NEW_PARENT_BONDS_TOTAL
        ),
        "lower_end_requires_60_new_parent_bonds": (
            int(
                valence_rows[0][
                    "required_new_parent_junction_bonds"
                ]
            )
            == EXPECTED_NEW_PARENT_BONDS_PER_END
        ),
        "upper_end_requires_60_new_parent_bonds": (
            int(
                valence_rows[1][
                    "required_new_parent_junction_bonds"
                ]
            )
            == EXPECTED_NEW_PARENT_BONDS_PER_END
        ),
        "target_aperture_is_finite_and_open": (
            aperture_diameter_nm > 0.0
            and aperture_radius_nm > 0.0
            and aperture_radius_nm
            < parent_rim_radius_nm
        ),
        "open_area_fraction_is_valid": (
            0.0
            < open_area_fraction
            < 1.0
        ),
        "candidate_set_contains_at_least_six_classes": (
            len(candidates) >= 6
        ),
        "all_candidate_scores_are_finite": all(
            math.isfinite(
                float(
                    row[
                        "final_screening_score"
                    ]
                )
            )
            for row in ranking_rows
        ),
        "primary_candidate_is_highest_ranked": (
            highest_ranked[
                "candidate_id"
            ]
            == PRIMARY_CANDIDATE
        ),
        "primary_candidate_requires_two_bonds_per_site": (
            int(
                primary[
                    "parent_bonds_per_terminal_site"
                ]
            )
            == MISSING_BONDS_PER_TERMINAL_SITE
        ),
        "primary_candidate_accounts_for_120_parent_bonds": (
            int(
                primary[
                    "expected_new_parent_bonds_total"
                ]
            )
            == EXPECTED_NEW_PARENT_BONDS_TOTAL
        ),
        "primary_candidate_preserves_target_aperture": (
            bool(
                primary[
                    "preserves_target_aperture"
                ]
            )
        ),
        "primary_candidate_resolves_end_asymmetry": (
            bool(
                primary[
                    "resolves_end_asymmetry"
                ]
            )
        ),
        "primary_candidate_does_not_require_parent_reconstruction": (
            not bool(
                primary[
                    "requires_parent_rim_reconstruction"
                ]
            )
        ),
        "parent_reconstruction_is_retained_as_contingency": (
            contingency[
                "status"
            ]
            == (
                "RETAIN_AS_PARENT_RECONSTRUCTION_CONTINGENCY"
            )
        ),
        "small_group_passivation_is_not_selected_as_cap": (
            next(
                row
                for row in ranking_rows
                if row[
                    "candidate_id"
                ]
                == "C0_SMALL_GROUP_PASSIVATION_ONLY"
            )[
                "status"
            ]
            == "REJECT_AS_R2_CAP_REPLACEMENT"
        ),
        "metal_spin_active_candidate_is_rejected": (
            next(
                row
                for row in ranking_rows
                if row[
                    "candidate_id"
                ]
                == "C6_METAL_OR_SPIN_ACTIVE_COORDINATED_CAP"
            )[
                "status"
            ]
            == "REJECT_AT_CURRENT_GATE"
        ),
    }

    failed_gates = [
        name
        for name, passed
        in gates.items()
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
        else
        "R2_POLAR_END_SPECIFIC_CANDIDATE_RANKING_REQUIRES_REVIEW"
    )

    required_next_step = (
        "BUILD_AND_VALIDATE_R2_GRADED_HETEROPOLAR_COLLAR_CONNECTIVITY_BLUEPRINT"
        if accepted
        else
        "REVIEW_R2_POLAR_END_SPECIFIC_CANDIDATE_RANKING_FAILURES"
    )

    design_contract = {
        "decision": decision,
        "primary_candidate_id": (
            PRIMARY_CANDIDATE
        ),
        "primary_candidate_is_final_chemistry": False,
        "primary_candidate_status": (
            primary[
                "status"
            ]
        ),
        "contingency_candidate_id": (
            CONTINGENCY_CANDIDATE
        ),
        "lower_parent_terminal_element": "B",
        "lower_required_first_junction_element": "N",
        "upper_parent_terminal_element": "N",
        "upper_required_first_junction_element": "B",
        "terminal_sites_per_end": (
            EXPECTED_TERMINAL_ATOMS_PER_END
        ),
        "missing_bond_incidences_per_site": (
            MISSING_BONDS_PER_TERMINAL_SITE
        ),
        "required_new_parent_bonds_lower": (
            EXPECTED_NEW_PARENT_BONDS_PER_END
        ),
        "required_new_parent_bonds_upper": (
            EXPECTED_NEW_PARENT_BONDS_PER_END
        ),
        "required_new_parent_bonds_total": (
            EXPECTED_NEW_PARENT_BONDS_TOTAL
        ),
        "target_aperture_diameter_nm": (
            aperture_diameter_nm
        ),
        "target_aperture_radius_nm": (
            aperture_radius_nm
        ),
        "target_open_area_fraction": (
            open_area_fraction
        ),
        "parent_rim_radius_nm": (
            parent_rim_radius_nm
        ),
        "required_radial_occlusion_nm": (
            required_radial_occlusion_nm
        ),
        "annular_area_per_end_nm2": (
            annular_area_nm2
        ),
        "screening_estimated_BN_atoms_per_end": (
            estimated_bn_atoms_per_end
        ),
        "same_element_parent_junction_bonds_authorized": False,
        "abrupt_unstrained_90_degree_seam_assumed": False,
        "formal_charge_assignment_authorized": False,
        "explicit_coordinate_generation_authorized": False,
        "topology_generation_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "mobile_MD_authorized": False,
        "multitemperature_MD_authorized": False,
        "QM_recalculation_authorized": False,
        "next_gate_scope": (
            "Connectivity blueprint only: enumerate "
            "parent-to-collar bonds, collar-ring composition, "
            "coordination closure, atom counts, ring topology, "
            "and graph-level aperture constraints."
        ),
        "required_next_step": (
            required_next_step
        ),
    }

    DESIGN_CONTRACT_JSON.write_text(
        json.dumps(
            design_contract,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    selection_summary = {
        "decision": decision,
        "candidate_count": len(
            candidates
        ),
        "primary_candidate": (
            PRIMARY_CANDIDATE
        ),
        "primary_rank": (
            primary[
                "rank"
            ]
        ),
        "primary_base_score": (
            primary[
                "base_weighted_score_0_100"
            ]
        ),
        "primary_risk_penalty": (
            primary[
                "risk_penalty_points"
            ]
        ),
        "primary_final_score": (
            primary[
                "final_screening_score"
            ]
        ),
        "primary_is_final_chemistry": False,
        "contingency_candidate": (
            CONTINGENCY_CANDIDATE
        ),
        "contingency_rank": (
            contingency[
                "rank"
            ]
        ),
        "terminal_sites_total": (
            degree1_atoms
        ),
        "terminal_sites_per_end": (
            EXPECTED_TERMINAL_ATOMS_PER_END
        ),
        "required_new_parent_bonds_per_end": (
            EXPECTED_NEW_PARENT_BONDS_PER_END
        ),
        "required_new_parent_bonds_total": (
            EXPECTED_NEW_PARENT_BONDS_TOTAL
        ),
        "lower_parent_element": "B",
        "lower_first_junction_element": "N",
        "upper_parent_element": "N",
        "upper_first_junction_element": "B",
        "target_aperture_diameter_nm": (
            aperture_diameter_nm
        ),
        "target_open_area_fraction": (
            open_area_fraction
        ),
        "failed_gates": (
            " | ".join(
                failed_gates
            )
        ),
        "connectivity_blueprint_authorized": (
            accepted
        ),
        "explicit_geometry_generation_authorized": False,
        "topology_generation_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "QM_authorized": False,
        "required_next_step": (
            required_next_step
        ),
    }

    write_csv(
        SELECTION_SUMMARY_CSV,
        [
            selection_summary
        ],
    )

    write_csv(
        GATES_CSV,
        [
            {
                "gate": name,
                "pass": passed,
            }
            for name, passed
            in gates.items()
        ],
    )

    source_manifest = [
        {
            "role": (
                "Gate3A_parent_summary"
            ),
            "file": relative(
                PARENT_SUMMARY
            ),
            "sha256": sha256(
                PARENT_SUMMARY
            ),
        },
        {
            "role": (
                "Gate3A_parent_gates"
            ),
            "file": relative(
                PARENT_GATES
            ),
            "sha256": sha256(
                PARENT_GATES
            ),
        },
        {
            "role": (
                "Gate3A_design_constraints"
            ),
            "file": relative(
                DESIGN_CONSTRAINTS
            ),
            "sha256": sha256(
                DESIGN_CONSTRAINTS
            ),
        },
        {
            "role": (
                "Gate3A_terminal_end_summary"
            ),
            "file": relative(
                TERMINAL_END_SUMMARY
            ),
            "sha256": sha256(
                TERMINAL_END_SUMMARY
            ),
        },
        {
            "role": (
                "Gate3A_topology_classification"
            ),
            "file": relative(
                TOPOLOGY_CLASSIFICATION
            ),
            "sha256": sha256(
                TOPOLOGY_CLASSIFICATION
            ),
        },
    ]

    write_csv(
        SOURCE_MANIFEST_CSV,
        source_manifest,
    )

    rank_lines = "\n".join(
        (
            f"- Rank {row['rank']}: "
            f"`{row['candidate_id']}` — "
            f"{row['final_screening_score']:.2f}; "
            f"**{row['status']}**"
        )
        for row in ranking_rows
    )

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed in gates.items()
    )

    REPORT_MD.write_text(
        f"""# R2 Polar End-Specific Valence-Completion Candidate Ranking

## Scope

This stage ranks chemical-architecture classes for replacing the
validated neutral steric R2 cap.

No coordinates, molecular topology, force-field parameters, partial
charges, minimization, MD, or QM calculation were generated.

The ranking is a deterministic screening instrument, not a formation
energy, synthetic-yield prediction, or proof of chemical stability.

## Parent terminal requirements

- Degree-1 terminal sites:
  **{degree1_atoms}**
- Lower end:
  **{lower_atoms} B-terminated sites**
- Upper end:
  **{upper_atoms} N-terminated sites**
- Missing bond incidences per parent terminal site:
  **{MISSING_BONDS_PER_TERMINAL_SITE}**
- Required new parent–junction bonds per end:
  **{EXPECTED_NEW_PARENT_BONDS_PER_END}**
- Required new parent–junction bonds total:
  **{EXPECTED_NEW_PARENT_BONDS_TOTAL}**

The lower B-terminated end requires an N-complementary first junction
layer. The upper N-terminated end requires a B-complementary first
junction layer. Same-element primary parent-junction bonds are not
authorized by this gate.

## Geometric target

- Aperture diameter:
  **{aperture_diameter_nm:.6f} nm**
- Aperture radius:
  **{aperture_radius_nm:.6f} nm**
- Open-area fraction:
  **{open_area_fraction:.6f}**
- Parent-rim radius:
  **{parent_rim_radius_nm:.6f} nm**
- Required radial occlusion:
  **{required_radial_occlusion_nm:.6f} nm**
- Annular area per end:
  **{annular_area_nm2:.6f} nm²**
- Screening estimate:
  **{estimated_bn_atoms_per_end:.3f} BN atoms/end**

## Candidate ranking

{rank_lines}

## Primary hypothesis

`{PRIMARY_CANDIDATE}`

This candidate uses a graded, end-specific, heteropolar BN collar
between the accepted tubular parent and the planar BN annulus:

- lower B parent → N-rich first collar layer;
- upper N parent → B-rich first collar layer;
- two additional parent-junction incidences per terminal site;
- 60 new parent-junction bonds per end;
- 120 new parent-junction bonds total;
- target central aperture preserved.

The candidate is not final chemistry. The next gate must determine
whether a graph with these requirements can be constructed without
invalid coordination, same-element primary junctions, disconnected
components, or an impossible ring topology.

## Contingency

`{CONTINGENCY_CANDIDATE}`

Parent-rim reconstruction is retained only as a contingency because it
would modify the accepted scaffold and require a separate structural
comparability assessment.

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Primary candidate:
  **{PRIMARY_CANDIDATE}**
- Primary screening score:
  **{primary['final_screening_score']:.2f}**
- Primary candidate is final chemistry:
  **NO**
- Connectivity-blueprint gate authorized:
  **{'YES' if accepted else 'NO'}**
- Explicit geometry generation authorized:
  **NO**
- Topology generation authorized:
  **NO**
- Energy minimization authorized:
  **NO**
- MD authorized:
  **NO**
- QM authorized:
  **NO**
- Required next step:
  `{required_next_step}`

## Literature-status limitation

Open BNNT edges and BNNT functionalization provide chemical context,
but no cited work is treated as direct evidence that this exact
graded collar–annulus junction is stable or synthesizable. Its first
test is therefore graph-level valence and connectivity feasibility.
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 polar end-specific candidate "
        "ranking completed."
    )

    print(
        "Terminal coordination d1/d2/d3: "
        f"{degree1_atoms}/"
        f"{degree2_atoms}/"
        f"{degree3_atoms}"
    )

    print(
        "Lower end atoms / B / N: "
        f"{lower_atoms}/"
        f"{lower_b}/"
        f"{lower_n}"
    )

    print(
        "Upper end atoms / B / N: "
        f"{upper_atoms}/"
        f"{upper_b}/"
        f"{upper_n}"
    )

    print(
        "Required new parent bonds "
        "lower/upper/total: "
        f"{EXPECTED_NEW_PARENT_BONDS_PER_END}/"
        f"{EXPECTED_NEW_PARENT_BONDS_PER_END}/"
        f"{EXPECTED_NEW_PARENT_BONDS_TOTAL}"
    )

    print(
        "Target aperture diameter / open-area fraction: "
        f"{aperture_diameter_nm:.6f}/"
        f"{open_area_fraction:.6f}"
    )

    print(
        "Candidate ranking:"
    )

    for row in ranking_rows:
        print(
            f"  Rank {row['rank']} | "
            f"{row['candidate_id']} | "
            f"base={row['base_weighted_score_0_100']:.2f} | "
            f"penalty={row['risk_penalty_points']:.2f} | "
            f"final={row['final_screening_score']:.2f} | "
            f"{row['status']}"
        )

    print(
        f"Primary candidate: {PRIMARY_CANDIDATE}"
    )

    print(
        "Primary candidate is final chemistry: NO"
    )

    print(
        f"Contingency candidate: {CONTINGENCY_CANDIDATE}"
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
        "Connectivity-blueprint gate authorized: "
        f"{'YES' if accepted else 'NO'}"
    )

    print(
        "Explicit geometry generation authorized: NO"
    )

    print(
        "Topology generation authorized: NO"
    )

    print(
        "Energy minimization authorized: NO"
    )

    print(
        "MD authorized: NO"
    )

    print(
        "QM authorized: NO"
    )

    print(
        f"Required next step: {required_next_step}"
    )

    for path in (
        VALENCE_REQUIREMENTS_CSV,
        CANDIDATE_MATRIX_CSV,
        RANKING_CSV,
        SELECTION_SUMMARY_CSV,
        GATES_CSV,
        DESIGN_CONTRACT_JSON,
        SOURCE_MANIFEST_CSV,
        REPORT_MD,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not accepted:
        raise RuntimeError(
            "R2 polar end-specific candidate ranking "
            "requires review."
        )


if __name__ == "__main__":
    main()
