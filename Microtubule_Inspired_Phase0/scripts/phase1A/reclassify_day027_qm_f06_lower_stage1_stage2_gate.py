#!/usr/bin/env python3
"""
Reclassify the QM_F06 LOWER Stage-1 local-structure result for Stage 2.

Blocking conditions
-------------------
1. Any original bonded interaction outside its validated chemical range.
2. Any hard long-range contact involving an artificial cap.
3. Any previously nonbonded pair sufficiently short to indicate possible
   unintended covalent bond formation.
4. Invalid or incomplete Stage-2 ORCA input.

Non-blocking but monitored
--------------------------
Compressed long-range contacts involving only original R2 atoms, provided
they do not enter the conservative unintended-covalent-bond thresholds.

Such contacts are explicit relaxation targets for Stage 2, where the bridge
and all hydrogen atoms are mobile.

No QM calculation is executed by this script.
"""

from __future__ import annotations

import csv
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

BOND_AUDIT = ROOT / (
    "runs/phase1A/day027_qm_f06_lower_stage1_validation/"
    "local_structure_audit/QM_F06_LOWER_STAGE1_bond_audit.csv"
)

NONBONDED_AUDIT = ROOT / (
    "runs/phase1A/day027_qm_f06_lower_stage1_validation/"
    "local_structure_audit/QM_F06_LOWER_STAGE1_nonbonded_audit.csv"
)

STAGE2_INPUT = ROOT / (
    "runs/phase1A/day026_qm_reference_catalog/QM_F06/"
    "orca_workflow/QM_F06_LOWER_CAPPED_REPAIRED/stage2.inp"
)

OUTPUT_DIR = ROOT / (
    "runs/phase1A/day027_qm_f06_lower_stage1_validation/"
    "stage2_readiness_gate"
)

# Conservative screening thresholds for possible unintended covalent
# connectivity. Distances above these values are not classified as new bonds.
UNINTENDED_COVALENT_THRESHOLDS_ANGSTROM = {
    tuple(sorted(("B", "B"))): 1.90,
    tuple(sorted(("B", "N"))): 1.85,
    tuple(sorted(("N", "N"))): 1.70,
    tuple(sorted(("B", "H"))): 1.45,
    tuple(sorted(("N", "H"))): 1.30,
    tuple(sorted(("H", "H"))): 0.90,
}


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty required file: {path}")


def read_csv(path: Path) -> list[dict[str, str]]:
    require_file(path)

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise RuntimeError(f"No data rows found in: {path}")

    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "atom_1",
        "element_1",
        "atom_2",
        "element_2",
        "graph_separation",
        "distance_angstrom",
        "distance_over_vdw_sum",
        "involves_artificial_cap",
        "original_hard_clash_status",
        "unintended_covalent_threshold_angstrom",
        "possible_unintended_covalent_contact",
        "stage2_classification",
        "blocks_stage2",
        "scientific_basis",
    ]

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()

        for row in rows:
            writer.writerow(
                {field: row.get(field, "") for field in fields}
            )


def parse_stage2_input(path: Path) -> dict[str, Any]:
    require_file(path)
    text = path.read_text(encoding="utf-8")

    match = re.search(
        r"(?ms)^\s*\*\s+xyz\s+(-?\d+)\s+(\d+)\s*$"
        r"(.*?)"
        r"^\s*\*\s*$",
        text,
    )

    if not match:
        raise RuntimeError("Stage-2 XYZ block not found.")

    atom_lines = [
        line.split()
        for line in match.group(3).splitlines()
        if line.strip()
    ]

    constraints = [
        int(value)
        for value in re.findall(
            r"(?m)^\s*\{\s*C\s+(\d+)\s+C\s*\}\s*$",
            text,
        )
    ]

    checks = {
        "atom_count_22": len(atom_lines) == 22,
        "charge_zero": int(match.group(1)) == 0,
        "multiplicity_one": int(match.group(2)) == 1,
        "constraint_count_4": len(constraints) == 4,
        "constraint_indices_expected": constraints == [0, 2, 11, 12],
        "promoted_from_stage1": (
            "Coordinates promoted from Stage 1" in text
        ),
        "uses_defgrid3": "DefGrid3" in text,
        "contains_optimization": bool(
            re.search(r"(?i)(^|\s)Opt(\s|$)", text)
        ),
        "no_obsolete_grid": (
            "Grid5" not in text
            and "FinalGrid6" not in text
        ),
    }

    return {
        "checks": checks,
        "gate_pass": all(checks.values()),
        "atom_count": len(atom_lines),
        "charge": int(match.group(1)),
        "multiplicity": int(match.group(2)),
        "constraints": constraints,
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    bond_rows = read_csv(BOND_AUDIT)
    nonbonded_rows = read_csv(NONBONDED_AUDIT)
    stage2 = parse_stage2_input(STAGE2_INPUT)

    failed_bonds = [
        row
        for row in bond_rows
        if row["bond_range_gate_pass"].strip().lower()
        != "true"
    ]

    classified_contacts: list[dict[str, Any]] = []
    cap_blocking_contacts = 0
    unintended_covalent_contacts = 0
    inherited_stage2_targets = 0

    for row in nonbonded_rows:
        element_pair = tuple(
            sorted((row["element_1"], row["element_2"]))
        )

        threshold = (
            UNINTENDED_COVALENT_THRESHOLDS_ANGSTROM[
                element_pair
            ]
        )

        measured = float(
            row["optimized_distance_angstrom"]
        )

        involves_cap = (
            row["involves_artificial_cap"].strip().lower()
            == "true"
        )

        original_hard = (
            row["hard_clash"].strip().lower()
            == "true"
        )

        possible_new_bond = measured <= threshold

        if possible_new_bond:
            classification = (
                "POSSIBLE_UNINTENDED_COVALENT_CONNECTIVITY"
            )
            blocks_stage2 = True
            basis = (
                "Previously nonbonded pair is within the conservative "
                "unintended-covalent screening threshold."
            )
            unintended_covalent_contacts += 1

        elif original_hard and involves_cap:
            classification = (
                "ARTIFICIAL_CAP_INDUCED_HARD_CONTACT"
            )
            blocks_stage2 = True
            basis = (
                "Hard long-range contact involves an artificial "
                "boundary cap and must be repaired before Stage 2."
            )
            cap_blocking_contacts += 1

        elif original_hard:
            classification = (
                "INHERITED_CONSTRAINED_GEOMETRY_STAGE2_RELAXATION_TARGET"
            )
            blocks_stage2 = False
            basis = (
                "Contact involves only original R2 atoms, remains "
                "outside unintended-covalent distance, and is an "
                "explicit target of the less-constrained Stage-2 "
                "optimization."
            )
            inherited_stage2_targets += 1

        else:
            classification = (
                "NONBLOCKING_CLOSE_CONTACT_MONITOR_DURING_STAGE2"
            )
            blocks_stage2 = False
            basis = (
                "Close contact does not meet a blocking chemical "
                "criterion; monitor its evolution during Stage 2."
            )

        classified_contacts.append(
            {
                "atom_1": row["atom_1"],
                "element_1": row["element_1"],
                "atom_2": row["atom_2"],
                "element_2": row["element_2"],
                "graph_separation": row["graph_separation"],
                "distance_angstrom": (
                    row["optimized_distance_angstrom"]
                ),
                "distance_over_vdw_sum": (
                    row["distance_over_vdw_sum"]
                ),
                "involves_artificial_cap": involves_cap,
                "original_hard_clash_status": original_hard,
                "unintended_covalent_threshold_angstrom": (
                    f"{threshold:.6f}"
                ),
                "possible_unintended_covalent_contact": (
                    possible_new_bond
                ),
                "stage2_classification": classification,
                "blocks_stage2": blocks_stage2,
                "scientific_basis": basis,
            }
        )

    blocking_contacts = [
        row
        for row in classified_contacts
        if row["blocks_stage2"]
    ]

    stage2_authorized = all(
        (
            len(failed_bonds) == 0,
            len(blocking_contacts) == 0,
            stage2["gate_pass"],
        )
    )

    decision = (
        "QM_F06_LOWER_STAGE1_CHEMICAL_CONNECTIVITY_PRESERVED_"
        "STAGE2_EXECUTION_AUTHORIZED"
        if stage2_authorized
        else
        "QM_F06_LOWER_STAGE1_BLOCKING_CHEMICAL_FAILURE_"
        "STAGE2_EXECUTION_NOT_AUTHORIZED"
    )

    contact_path = (
        OUTPUT_DIR
        / "QM_F06_LOWER_STAGE1_contact_reclassification.csv"
    )
    write_csv(contact_path, classified_contacts)

    gate_rows = [
        {
            "gate": "ORIGINAL_BONDED_CONNECTIVITY",
            "required": "ZERO_BOND_RANGE_FAILURES",
            "observed": len(failed_bonds),
            "pass": len(failed_bonds) == 0,
        },
        {
            "gate": "ARTIFICIAL_CAP_HARD_CONTACTS",
            "required": "ZERO",
            "observed": cap_blocking_contacts,
            "pass": cap_blocking_contacts == 0,
        },
        {
            "gate": "UNINTENDED_COVALENT_CONNECTIVITY",
            "required": "ZERO",
            "observed": unintended_covalent_contacts,
            "pass": unintended_covalent_contacts == 0,
        },
        {
            "gate": "STAGE2_STATIC_INPUT",
            "required": "PASS",
            "observed": stage2["gate_pass"],
            "pass": stage2["gate_pass"],
        },
    ]

    gates_path = (
        OUTPUT_DIR
        / "QM_F06_LOWER_STAGE1_stage2_readiness_gates.csv"
    )

    with gates_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "gate",
                "required",
                "observed",
                "pass",
            ],
        )
        writer.writeheader()
        writer.writerows(gate_rows)

    report_path = (
        OUTPUT_DIR
        / "QM_F06_LOWER_STAGE1_STAGE2_READINESS_REPORT.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 LOWER Stage-1 → Stage-2 Readiness Gate — Day027",
                "",
                f"## Decision: **{decision}**",
                "",
                "## Corrected gate logic",
                "",
                (
                    "A van der Waals compression ratio alone is not "
                    "treated as proof of bond formation or structural "
                    "failure in a constrained covalent cluster."
                ),
                "",
                (
                    "Stage 2 is blocked only by broken original bonds, "
                    "cap-induced hard contacts, possible unintended "
                    "covalent connectivity, or an invalid Stage-2 input."
                ),
                "",
                "## Results",
                "",
                f"- Original bonded interactions: **{len(bond_rows)}**",
                f"- Bond-range failures: **{len(failed_bonds)}**",
                (
                    "- Hard contacts involving artificial caps: "
                    f"**{cap_blocking_contacts}**"
                ),
                (
                    "- Possible unintended covalent contacts: "
                    f"**{unintended_covalent_contacts}**"
                ),
                (
                    "- Inherited constrained-geometry contacts assigned "
                    f"as Stage-2 relaxation targets: "
                    f"**{inherited_stage2_targets}**"
                ),
                (
                    "- Stage-2 static input gate: "
                    f"**{'PASS' if stage2['gate_pass'] else 'FAIL'}**"
                ),
                "",
                "## Stage-2 role",
                "",
                (
                    "Stage 2 releases all hydrogen atoms while retaining "
                    "only four peripheral heavy-atom constraints. It is "
                    "therefore the appropriate controlled test of whether "
                    "the remaining local strain relaxes without changing "
                    "the intended B–N–B–N connectivity."
                ),
                "",
                "## Authorization state",
                "",
                (
                    "- Stage-2 input preparation: "
                    f"**{'AUTHORIZED' if stage2_authorized else 'NOT AUTHORIZED'}**"
                ),
                (
                    "- Stage-2 execution: "
                    f"**{'AUTHORIZED' if stage2_authorized else 'NOT AUTHORIZED'}**"
                ),
                "- Stage-2 calculation executed by this gate: **NO**",
                "- Force-field parameter adoption: **NOT AUTHORIZED**",
                "- MD execution: **NOT AUTHORIZED**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = {
        "decision": decision,
        "bond_range_failures": len(failed_bonds),
        "artificial_cap_hard_contacts": (
            cap_blocking_contacts
        ),
        "possible_unintended_covalent_contacts": (
            unintended_covalent_contacts
        ),
        "inherited_stage2_relaxation_targets": (
            inherited_stage2_targets
        ),
        "stage2_static_input_gate_pass": (
            stage2["gate_pass"]
        ),
        "stage2_execution_authorized": (
            stage2_authorized
        ),
        "stage2_executed_by_this_script": False,
        "force_field_parameter_adoption_authorized": False,
        "md_authorized": False,
        "required_next_step": (
            "BUILD_AND_RUN_STAGE2_LOWER_PREFLIGHT"
            if stage2_authorized
            else "REVIEW_BLOCKING_STAGE2_GATE"
        ),
        "reference": (
            "Cordero et al., Covalent radii revisited, "
            "Dalton Transactions, 2008, 2832-2838."
        ),
    }

    summary_path = (
        OUTPUT_DIR
        / "QM_F06_LOWER_STAGE1_stage2_readiness_summary.json"
    )
    summary_path.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print("QM_F06 LOWER Stage-2 readiness gate completed.")
    print(f"Decision: {decision}")
    print("Bond-range failures:", len(failed_bonds))
    print(
        "Artificial-cap hard contacts:",
        cap_blocking_contacts,
    )
    print(
        "Possible unintended covalent contacts:",
        unintended_covalent_contacts,
    )
    print(
        "Inherited Stage-2 relaxation targets:",
        inherited_stage2_targets,
    )
    print(
        "Stage-2 execution authorized:",
        stage2_authorized,
    )
    print("Stage-2 calculation executed: False")
    print(f"Report: {report_path}")


if __name__ == "__main__":
    main()
