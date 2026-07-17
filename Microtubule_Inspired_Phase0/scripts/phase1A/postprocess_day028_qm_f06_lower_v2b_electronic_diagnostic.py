#!/usr/bin/env python3
"""
Post-process the completed QM_F06 LOWER V2-B ORCA electronic diagnostic.

Extracts and validates:

- normal termination and SCF convergence;
- final single-point energy;
- Mayer bond orders;
- Mayer target-pair upper bound when absent below print threshold;
- Hirshfeld charges;
- MBIS charges;
- CHELPG charges;
- charge sums;
- target-atom values;
- comparison with the two intended covalent bonds involving the target atoms.

No charge model or force-field parameter is adopted.
"""

from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

DIAGNOSTIC_ROOT = ROOT / (
    "runs/phase1A/"
    "day028_qm_f06_lower_boundary_v2b_electronic_diagnostic"
)

EXECUTION_ROOT = DIAGNOSTIC_ROOT / "orca_executions"

ATOM_MAP = DIAGNOSTIC_ROOT / (
    "QM_F06_LOWER_BOUNDARY_V2B_atom_index_map.csv"
)

PREPARATION_SUMMARY = DIAGNOSTIC_ROOT / (
    "QM_F06_LOWER_BOUNDARY_V2B_"
    "electronic_diagnostic_summary.json"
)

OUTPUT_DIR = DIAGNOSTIC_ROOT / "postprocessing"

TARGET_B_ID = "BR4:LOWER:00:3"
TARGET_H_ID = "H4:LOWER:0017:0"

TARGET_B_INDEX = 5
TARGET_H_INDEX = 20

# Bond-order printing threshold specified in the ORCA input.
MAYER_PRINT_THRESHOLD = 0.01


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty file: {path}")


def latest_execution() -> Path:
    candidates = sorted(
        path
        for path in EXECUTION_ROOT.glob("electronic_diagnostic_*")
        if path.is_dir()
    )

    if not candidates:
        raise RuntimeError(
            f"No electronic-diagnostic execution found in {EXECUTION_ROOT}"
        )

    return candidates[-1]


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"Attempted to write empty table: {path}")

    fields: list[str] = []

    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(
            {
                field: row.get(field, "")
                for field in fields
            }
            for row in rows
        )


def parse_charge_block(
    text: str,
    start_pattern: str,
    end_patterns: tuple[str, ...],
    expected_atoms: int = 28,
) -> list[dict[str, Any]]:
    start = re.search(start_pattern, text, flags=re.I | re.M)

    if not start:
        raise RuntimeError(
            f"Charge section not found: {start_pattern}"
        )

    remainder = text[start.end():]

    end_positions = []

    for pattern in end_patterns:
        match = re.search(pattern, remainder, flags=re.I | re.M)
        if match:
            end_positions.append(match.start())

    block = (
        remainder[:min(end_positions)]
        if end_positions
        else remainder
    )

    rows = []

    # Supports lines of the common form:
    #  0 N  -0.238227 ...
    #  0 N : -0.238227
    line_pattern = re.compile(
        r"^\s*(\d+)\s+([A-Za-z]{1,2})"
        r"(?:\s*:)?\s+(-?\d+\.\d+(?:[Ee][+-]?\d+)?)",
        flags=re.M,
    )

    for match in line_pattern.finditer(block):
        rows.append(
            {
                "index_0based": int(match.group(1)),
                "element": match.group(2),
                "charge_e": float(match.group(3)),
            }
        )

    unique = {
        row["index_0based"]: row
        for row in rows
    }

    rows = [
        unique[index]
        for index in sorted(unique)
    ]

    if len(rows) != expected_atoms:
        raise RuntimeError(
            f"Expected {expected_atoms} charge records; "
            f"found {len(rows)} for pattern {start_pattern}"
        )

    return rows


def parse_mayer_orders(text: str) -> list[dict[str, Any]]:
    section = re.search(
        r"Mayer bond orders larger than\s+"
        r"(-?\d+\.\d+)(.*?)(?:\n\s*-{5,}|\n\s*HIRSHFELD ANALYSIS)",
        text,
        flags=re.I | re.S,
    )

    if not section:
        raise RuntimeError("Mayer bond-order table not found.")

    threshold = float(section.group(1))
    block = section.group(2)

    pattern = re.compile(
        r"B\(\s*(\d+)-([A-Za-z]{1,2})\s*,\s*"
        r"(\d+)-([A-Za-z]{1,2})\s*\)\s*:\s*"
        r"(-?\d+\.\d+)"
    )

    rows = []

    for match in pattern.finditer(block):
        index_1 = int(match.group(1))
        element_1 = match.group(2)
        index_2 = int(match.group(3))
        element_2 = match.group(4)
        value = float(match.group(5))

        rows.append(
            {
                "index_1_0based": index_1,
                "element_1": element_1,
                "index_2_0based": index_2,
                "element_2": element_2,
                "mayer_bond_order": value,
                "absolute_mayer_bond_order": abs(value),
                "printed_threshold": threshold,
            }
        )

    if not rows:
        raise RuntimeError("No Mayer bond orders parsed.")

    return rows


def find_pair(
    rows: list[dict[str, Any]],
    index_1: int,
    index_2: int,
) -> dict[str, Any] | None:
    target = {index_1, index_2}

    for row in rows:
        observed = {
            row["index_1_0based"],
            row["index_2_0based"],
        }

        if observed == target:
            return row

    return None


def merge_charge_map(
    scheme: str,
    rows: list[dict[str, Any]],
    atom_map: dict[int, dict[str, str]],
) -> list[dict[str, Any]]:
    merged = []

    for row in rows:
        index = row["index_0based"]
        atom = atom_map[index]

        merged.append(
            {
                "scheme": scheme,
                "index_0based": index,
                "atom_id": atom["atom_id"],
                "element": row["element"],
                "atom_role": atom["atom_role"],
                "node_type": atom["node_type"],
                "charge_e": f"{row['charge_e']:.10f}",
                "target_contact_atom": (
                    atom["target_contact_atom"].lower() == "true"
                ),
            }
        )

    return merged


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    execution_dir = latest_execution()
    output_path = execution_dir / "diagnostic.out"
    validation_path = execution_dir / "diagnostic_validation.json"

    for path in (
        output_path,
        validation_path,
        ATOM_MAP,
        PREPARATION_SUMMARY,
    ):
        require_file(path)

    text = output_path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    validation = json.loads(
        validation_path.read_text(encoding="utf-8")
    )

    preparation = json.loads(
        PREPARATION_SUMMARY.read_text(encoding="utf-8")
    )

    with ATOM_MAP.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        atom_rows = list(csv.DictReader(handle))

    atom_map = {
        int(row["index_0based"]): row
        for row in atom_rows
    }

    if len(atom_map) != 28:
        raise RuntimeError(
            f"Expected 28 mapped atoms; found {len(atom_map)}"
        )

    if atom_map[TARGET_B_INDEX]["atom_id"] != TARGET_B_ID:
        raise RuntimeError("Target B index-map mismatch.")

    if atom_map[TARGET_H_INDEX]["atom_id"] != TARGET_H_ID:
        raise RuntimeError("Target H index-map mismatch.")

    execution_checks = {
        "return_code_zero": validation["return_code"] == 0,
        "normal_termination": validation["normal_termination"] is True,
        "scf_converged": validation["scf_converged"] is True,
        "execution_gate_pass": (
            validation["execution_gate_pass"] is True
        ),
        "complete_diagnostic_output": (
            validation["complete_diagnostic_output"] is True
        ),
    }

    if not all(execution_checks.values()):
        raise RuntimeError(
            f"Execution validation failed: {execution_checks}"
        )

    energy_matches = re.findall(
        r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)",
        text,
    )

    if not energy_matches:
        raise RuntimeError("Final energy not found.")

    final_energy = float(energy_matches[-1])

    mayer_rows = parse_mayer_orders(text)

    target_mayer = find_pair(
        mayer_rows,
        TARGET_B_INDEX,
        TARGET_H_INDEX,
    )

    target_mayer_printed = target_mayer is not None

    target_mayer_upper_bound = (
        abs(target_mayer["mayer_bond_order"])
        if target_mayer_printed
        else MAYER_PRINT_THRESHOLD
    )

    # Intended covalent neighbors used as electronic references.
    # B5-H9 is the intended B-H bond at the target bridge B.
    # N2-H20 is the intended N-H bond containing the target hydrogen.
    reference_bh = find_pair(mayer_rows, 5, 9)
    reference_nh = find_pair(mayer_rows, 2, 20)

    if reference_bh is None:
        raise RuntimeError("Reference B5-H9 Mayer order absent.")

    if reference_nh is None:
        raise RuntimeError("Reference N2-H20 Mayer order absent.")

    hirshfeld = parse_charge_block(
        text,
        r"^\s*HIRSHFELD ANALYSIS\s*$",
        (
            r"^\s*MBIS ANALYSIS\s*$",
            r"^\s*CHELPG",
        ),
    )

    mbis = parse_charge_block(
        text,
        r"^\s*MBIS ANALYSIS\s*$",
        (
            r"^\s*MBIS VALENCE-SHELL DATA",
            r"^\s*CHELPG",
        ),
    )

    chelpg = parse_charge_block(
        text,
        r"^\s*CHELPG Charges\s*$",
        (
            r"^\s*CHELPG charges calculated",
            r"FINAL SINGLE POINT ENERGY",
        ),
    )

    all_charge_rows = []

    for scheme, rows in (
        ("HIRSHFELD", hirshfeld),
        ("MBIS", mbis),
        ("CHELPG", chelpg),
    ):
        all_charge_rows.extend(
            merge_charge_map(scheme, rows, atom_map)
        )

    write_csv(
        OUTPUT_DIR / "QM_F06_LOWER_V2B_charge_comparison.csv",
        all_charge_rows,
    )

    write_csv(
        OUTPUT_DIR / "QM_F06_LOWER_V2B_mayer_bond_orders.csv",
        mayer_rows,
    )

    charge_summaries = {}

    for scheme, rows in (
        ("HIRSHFELD", hirshfeld),
        ("MBIS", mbis),
        ("CHELPG", chelpg),
    ):
        charge_summaries[scheme] = {
            "charge_sum_e": sum(
                row["charge_e"]
                for row in rows
            ),
            "target_b_charge_e": rows[TARGET_B_INDEX]["charge_e"],
            "target_h_charge_e": rows[TARGET_H_INDEX]["charge_e"],
        }

    # Mayer interpretation:
    # The target pair is unprinted at a 0.01 cutoff, while true local
    # covalent B-H and N-H bonds have Mayer orders close to unity.
    target_no_significant_covalent_character = (
        not target_mayer_printed
        and target_mayer_upper_bound <= 0.01
        and abs(reference_bh["mayer_bond_order"]) >= 0.75
        and abs(reference_nh["mayer_bond_order"]) >= 0.75
    )

    electronic_gate_pass = all(
        (
            all(execution_checks.values()),
            target_no_significant_covalent_character,
            all(
                abs(values["charge_sum_e"]) <= 1.0e-3
                for values in charge_summaries.values()
            ),
        )
    )

    decision = (
        "QM_F06_LOWER_V2B_RESIDUAL_CONTACT_HAS_NO_"
        "SIGNIFICANT_COVALENT_CHARACTER_ELECTRONIC_GATE_PASS"
        if electronic_gate_pass
        else
        "QM_F06_LOWER_V2B_ELECTRONIC_GATE_REQUIRES_REVIEW"
    )

    target_summary_rows = [
        {
            "quantity": "TARGET_PAIR",
            "atom_1": TARGET_B_ID,
            "index_1_0based": TARGET_B_INDEX,
            "atom_2": TARGET_H_ID,
            "index_2_0based": TARGET_H_INDEX,
            "mayer_printed": target_mayer_printed,
            "mayer_value_or_upper_bound": (
                target_mayer["mayer_bond_order"]
                if target_mayer_printed
                else f"<{MAYER_PRINT_THRESHOLD:.2f}"
            ),
            "interpretation": (
                "NO_SIGNIFICANT_COVALENT_CHARACTER"
                if target_no_significant_covalent_character
                else "REQUIRES_REVIEW"
            ),
        },
        {
            "quantity": "REFERENCE_COVALENT_BH",
            "atom_1": atom_map[5]["atom_id"],
            "index_1_0based": 5,
            "atom_2": atom_map[9]["atom_id"],
            "index_2_0based": 9,
            "mayer_printed": True,
            "mayer_value_or_upper_bound": (
                reference_bh["mayer_bond_order"]
            ),
            "interpretation": "INTENDED_COVALENT_BH",
        },
        {
            "quantity": "REFERENCE_COVALENT_NH",
            "atom_1": atom_map[2]["atom_id"],
            "index_1_0based": 2,
            "atom_2": atom_map[20]["atom_id"],
            "index_2_0based": 20,
            "mayer_printed": True,
            "mayer_value_or_upper_bound": (
                reference_nh["mayer_bond_order"]
            ),
            "interpretation": "INTENDED_COVALENT_NH",
        },
    ]

    write_csv(
        OUTPUT_DIR / "QM_F06_LOWER_V2B_target_mayer_analysis.csv",
        target_summary_rows,
    )

    summary = {
        "decision": decision,
        "execution_directory": str(
            execution_dir.relative_to(ROOT)
        ),
        "output_file": str(output_path.relative_to(ROOT)),
        "final_single_point_energy_hartree": final_energy,
        "target_b_atom_id": TARGET_B_ID,
        "target_b_index_0based": TARGET_B_INDEX,
        "target_h_atom_id": TARGET_H_ID,
        "target_h_index_0based": TARGET_H_INDEX,
        "mayer_print_threshold": MAYER_PRINT_THRESHOLD,
        "target_mayer_printed": target_mayer_printed,
        "target_mayer_upper_bound": target_mayer_upper_bound,
        "reference_b5_h9_mayer": (
            reference_bh["mayer_bond_order"]
        ),
        "reference_n2_h20_mayer": (
            reference_nh["mayer_bond_order"]
        ),
        "target_no_significant_covalent_character": (
            target_no_significant_covalent_character
        ),
        "charge_summaries": charge_summaries,
        "electronic_diagnostic_gate_pass": electronic_gate_pass,
        "lower_geometry_formally_accepted": electronic_gate_pass,
        "esp_resp_protocol_definition_authorized": electronic_gate_pass,
        "esp_resp_execution_authorized": False,
        "charge_adoption_authorized": False,
        "force_field_parameter_adoption_authorized": False,
    }

    summary_path = OUTPUT_DIR / (
        "QM_F06_LOWER_V2B_electronic_validation_summary.json"
    )

    summary_path.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    report_path = OUTPUT_DIR / (
        "QM_F06_LOWER_V2B_ELECTRONIC_VALIDATION_REPORT.md"
    )

    report_path.write_text(
        "\n".join(
            [
                "# QM_F06 LOWER V2-B Electronic Validation — Day028",
                "",
                f"## Decision: **{decision}**",
                "",
                "## Execution",
                "",
                "- Return code: **0**",
                "- SCF convergence: **YES**",
                "- Normal ORCA termination: **YES**",
                (
                    "- Final single-point energy: "
                    f"**{final_energy:.12f} Eh**"
                ),
                "",
                "## Target contact",
                "",
                f"- B atom: `{TARGET_B_ID}` — ORCA index `{TARGET_B_INDEX}`",
                f"- H atom: `{TARGET_H_ID}` — ORCA index `{TARGET_H_INDEX}`",
                (
                    "- Mayer print threshold: "
                    f"**{MAYER_PRINT_THRESHOLD:.2f}**"
                ),
                (
                    "- Target Mayer bond order printed: "
                    f"**{'YES' if target_mayer_printed else 'NO'}**"
                ),
                (
                    "- Target Mayer bond-order result: "
                    f"**{'observed ' + str(target_mayer['mayer_bond_order']) if target_mayer_printed else '< 0.01'}**"
                ),
                "",
                "## Covalent-reference bonds",
                "",
                (
                    "- Intended B5–H9 Mayer bond order: "
                    f"**{reference_bh['mayer_bond_order']:.4f}**"
                ),
                (
                    "- Intended N2–H20 Mayer bond order: "
                    f"**{reference_nh['mayer_bond_order']:.4f}**"
                ),
                "",
                "## Charge analyses",
                "",
                *[
                    (
                        f"- {scheme}: sum = "
                        f"**{values['charge_sum_e']:+.6e} e**; "
                        f"B5 = **{values['target_b_charge_e']:+.6f} e**; "
                        f"H20 = **{values['target_h_charge_e']:+.6f} e**"
                    )
                    for scheme, values
                    in charge_summaries.items()
                ],
                "",
                "## Interpretation",
                "",
                (
                    "The B5···H20 target pair is absent from the Mayer "
                    "table printed at a 0.01 threshold, while the intended "
                    "local B–H and N–H bonds have Mayer bond orders close "
                    "to unity. The compressed 1–5 contact therefore has "
                    "no significant covalent-bond character at the "
                    "PBE0-D4/def2-TZVP level."
                ),
                "",
                (
                    "Hirshfeld, MBIS and CHELPG values remain diagnostic "
                    "only. They are not adopted as force-field charges by "
                    "this validation."
                ),
                "",
                "## Authorization state",
                "",
                (
                    "- LOWER V2-B geometry formally accepted: "
                    f"**{'YES' if electronic_gate_pass else 'NO'}**"
                ),
                (
                    "- ESP/RESP protocol definition: "
                    f"**{'AUTHORIZED' if electronic_gate_pass else 'NOT AUTHORIZED'}**"
                ),
                "- ESP/RESP execution: **NOT AUTHORIZED**",
                "- Charge adoption: **NOT AUTHORIZED**",
                "- Force-field parameter adoption: **NOT AUTHORIZED**",
                "",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    preparation[
        "electronic_postprocessing_completed"
    ] = True

    preparation[
        "electronic_validation_summary"
    ] = str(summary_path.relative_to(ROOT))

    preparation[
        "electronic_validation_report"
    ] = str(report_path.relative_to(ROOT))

    preparation[
        "electronic_diagnostic_gate_pass"
    ] = electronic_gate_pass

    preparation[
        "esp_resp_parameter_adoption_authorized"
    ] = False

    preparation[
        "force_field_parameter_adoption_authorized"
    ] = False

    PREPARATION_SUMMARY.write_text(
        json.dumps(preparation, indent=2) + "\n",
        encoding="utf-8",
    )

    print("Electronic diagnostic post-processing completed.")
    print("Decision:", decision)
    print("Final energy:", f"{final_energy:.12f} Eh")
    print("Target Mayer printed:", target_mayer_printed)
    print(
        "Target Mayer result:",
        (
            target_mayer["mayer_bond_order"]
            if target_mayer_printed
            else "<0.01"
        ),
    )
    print(
        "Reference B5-H9 Mayer:",
        reference_bh["mayer_bond_order"],
    )
    print(
        "Reference N2-H20 Mayer:",
        reference_nh["mayer_bond_order"],
    )
    print("Electronic gate pass:", electronic_gate_pass)
    print("Report:", report_path)


if __name__ == "__main__":
    main()
