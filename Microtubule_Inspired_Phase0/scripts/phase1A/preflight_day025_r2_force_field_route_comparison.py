#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]

DAY024 = (
    ROOT
    / "runs/phase1A/day024_chemical_end_rim_design"
)

OUT = (
    ROOT
    / "runs/phase1A/day025_force_field_route_comparison"
    / "38_r2_force_field_route_comparison_preflight"
)

INPUTS = {
    "gate29_summary": (
        DAY024
        / "29_r2_chemical_realizability_and_parameterization_scope"
        / "r2_chemical_realizability_and_parameterization_scope_summary.csv"
    ),
    "gate29_environment_inventory": (
        DAY024
        / "29_r2_chemical_realizability_and_parameterization_scope"
        / "r2_local_chemical_environment_inventory.csv"
    ),
    "gate29_assignments": (
        DAY024
        / "29_r2_chemical_realizability_and_parameterization_scope"
        / "r2_node_chemical_environment_assignments.csv"
    ),
    "gate29_critical_centers": (
        DAY024
        / "29_r2_chemical_realizability_and_parameterization_scope"
        / "r2_parameterization_critical_centers.csv"
    ),
    "gate29_scope": (
        DAY024
        / "29_r2_chemical_realizability_and_parameterization_scope"
        / "r2_parameterization_scope.csv"
    ),
    "gate30_summary": (
        DAY024
        / "30_r2_force_field_coverage_preflight"
        / "r2_force_field_coverage_preflight_summary.csv"
    ),
    "gate30_bonds": (
        DAY024
        / "30_r2_force_field_coverage_preflight"
        / "r2_required_bond_terms.csv"
    ),
    "gate30_angles": (
        DAY024
        / "30_r2_force_field_coverage_preflight"
        / "r2_required_angle_terms.csv"
    ),
    "gate30_torsions": (
        DAY024
        / "30_r2_force_field_coverage_preflight"
        / "r2_required_torsion_terms.csv"
    ),
    "gate30_impropers": (
        DAY024
        / "30_r2_force_field_coverage_preflight"
        / "r2_required_improper_centers.csv"
    ),
    "gate30_qm_fragments": (
        DAY024
        / "30_r2_force_field_coverage_preflight"
        / "r2_preliminary_qm_fragment_classes.csv"
    ),
    "gate37_summary": (
        DAY024
        / "37_r2_lele_si_parameter_file_audit"
        / "lele_2022_si_parameter_file_audit_summary.csv"
    ),
    "gate37_bnh_terms": (
        DAY024
        / "37_r2_lele_si_parameter_file_audit"
        / "lele_2022_reaxff_bnh_term_inventory.csv"
    ),
    "gate37_domain_assessment": (
        DAY024
        / "37_r2_lele_si_parameter_file_audit"
        / "lele_2022_R2_domain_assessment.csv"
    ),
}

SUMMARY = OUT / "r2_force_field_route_comparison_preflight_summary.csv"
SCHEMA = OUT / "r2_force_field_route_comparison_input_schema.csv"
JSON_OUT = OUT / "r2_force_field_route_comparison_preflight.json"
REPORT = OUT / "R2_FORCE_FIELD_ROUTE_COMPARISON_PREFLIGHT_DAY025.md"


def require_file(path: Path) -> None:
    if not path.is_file():
        raise RuntimeError(f"Missing required input: {path}")

    if path.stat().st_size == 0:
        raise RuntimeError(f"Empty required input: {path}")


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
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
        raise RuntimeError(f"No CSV header found: {path}")

    return fields, rows


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


def first_nonempty_value(
    rows: list[dict[str, str]],
    keys: list[str],
) -> str:
    for row in rows:
        for key in keys:
            value = row.get(key, "").strip()

            if value:
                return value

    return ""


def compact_preview(
    row: dict[str, str],
    limit: int = 6,
) -> str:
    parts = []

    for key, value in row.items():
        value = value.strip()

        if not value:
            continue

        if len(value) > 80:
            value = value[:77] + "..."

        parts.append(f"{key}={value}")

        if len(parts) >= limit:
            break

    return " | ".join(parts)


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    schema_rows = []
    input_data: dict[str, dict[str, Any]] = {}

    for input_name, path in INPUTS.items():
        fields, rows = read_csv(path)

        input_data[input_name] = {
            "path": path,
            "fields": fields,
            "rows": rows,
        }

        schema_rows.append(
            {
                "input_name": input_name,
                "relative_path": str(
                    path.relative_to(ROOT)
                ),
                "row_count": len(rows),
                "column_count": len(fields),
                "columns": " | ".join(fields),
                "first_row_preview": (
                    compact_preview(rows[0])
                    if rows
                    else "NO_ROWS"
                ),
            }
        )

    gate29_rows = input_data["gate29_summary"]["rows"]
    gate30_rows = input_data["gate30_summary"]["rows"]
    gate37_rows = input_data["gate37_summary"]["rows"]

    gate29_decision = first_nonempty_value(
        gate29_rows,
        [
            "decision",
            "Decision",
            "outcome",
            "status",
        ],
    )

    gate30_decision = first_nonempty_value(
        gate30_rows,
        [
            "decision",
            "Decision",
            "outcome",
            "status",
        ],
    )

    gate37_decision = first_nonempty_value(
        gate37_rows,
        [
            "decision",
            "Decision",
            "outcome",
            "status",
        ],
    )

    counts = {
        "chemical_environment_rows": len(
            input_data["gate29_environment_inventory"]["rows"]
        ),
        "chemical_assignment_rows": len(
            input_data["gate29_assignments"]["rows"]
        ),
        "critical_center_rows": len(
            input_data["gate29_critical_centers"]["rows"]
        ),
        "parameterization_scope_rows": len(
            input_data["gate29_scope"]["rows"]
        ),
        "required_bond_rows": len(
            input_data["gate30_bonds"]["rows"]
        ),
        "required_angle_rows": len(
            input_data["gate30_angles"]["rows"]
        ),
        "required_torsion_rows": len(
            input_data["gate30_torsions"]["rows"]
        ),
        "required_improper_rows": len(
            input_data["gate30_impropers"]["rows"]
        ),
        "preliminary_qm_fragment_rows": len(
            input_data["gate30_qm_fragments"]["rows"]
        ),
        "lele_bnh_parameter_rows": len(
            input_data["gate37_bnh_terms"]["rows"]
        ),
        "lele_domain_rows": len(
            input_data["gate37_domain_assessment"]["rows"]
        ),
    }

    gates = {
        "all_required_inputs_exist": True,
        "gate29_summary_has_decision": bool(
            gate29_decision
        ),
        "gate30_summary_has_decision": bool(
            gate30_decision
        ),
        "gate37_summary_has_decision": bool(
            gate37_decision
        ),
        "chemical_environment_inventory_nonempty": (
            counts["chemical_environment_rows"] > 0
        ),
        "required_bond_inventory_nonempty": (
            counts["required_bond_rows"] > 0
        ),
        "required_angle_inventory_nonempty": (
            counts["required_angle_rows"] > 0
        ),
        "required_torsion_inventory_nonempty": (
            counts["required_torsion_rows"] > 0
        ),
        "lele_parameter_inventory_nonempty": (
            counts["lele_bnh_parameter_rows"] > 0
        ),
        "no_parameter_adoption_or_simulation_performed": True,
    }

    failed_gates = [
        name
        for name, passed in gates.items()
        if not passed
    ]

    decision = (
        "R2_FORCE_FIELD_ROUTE_COMPARISON_DATA_CONTRACT_VALIDATED"
        if not failed_gates
        else
        "R2_FORCE_FIELD_ROUTE_COMPARISON_DATA_CONTRACT_REQUIRES_REVIEW"
    )

    next_step = (
        "BUILD_R2_FORCE_FIELD_ROUTE_COVERAGE_AND_RISK_MATRIX"
        if not failed_gates
        else
        "REVIEW_INPUT_SCHEMAS_AND_PREVIOUS_GATE_DECISIONS"
    )

    summary = {
        "decision": decision,
        "gate29_decision": gate29_decision,
        "gate30_decision": gate30_decision,
        "gate37_decision": gate37_decision,
        **counts,
        "failed_gates": (
            " | ".join(failed_gates)
        ),
        "parameter_adoption_authorized": False,
        "topology_generation_authorized": False,
        "formal_charge_assignment_authorized": False,
        "force_field_parameterization_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "QM_calculation_authorized": False,
        "required_next_step": next_step,
    }

    write_csv(
        SCHEMA,
        schema_rows,
    )

    write_csv(
        SUMMARY,
        [summary],
    )

    JSON_OUT.write_text(
        json.dumps(
            {
                "summary": summary,
                "gates": gates,
                "schemas": schema_rows,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    REPORT.write_text(
        f"""# R2 Force-Field Route Comparison Preflight — Day025

## Previous decisions

- Gate 29:
  `{gate29_decision}`
- Gate 30:
  `{gate30_decision}`
- Gate 37:
  `{gate37_decision}`

## Input counts

- Chemical environments:
  **{counts['chemical_environment_rows']}**
- Node/environment assignments:
  **{counts['chemical_assignment_rows']}**
- Critical centers:
  **{counts['critical_center_rows']}**
- Required bond terms:
  **{counts['required_bond_rows']}**
- Required angle terms:
  **{counts['required_angle_rows']}**
- Required torsion terms:
  **{counts['required_torsion_rows']}**
- Improper/planarity centers:
  **{counts['required_improper_rows']}**
- Preliminary QM fragment classes:
  **{counts['preliminary_qm_fragment_rows']}**
- Lele B/N/H parameter records:
  **{counts['lele_bnh_parameter_rows']}**
- Lele domain-assessment rows:
  **{counts['lele_domain_rows']}**

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Parameter adoption:
  **NOT AUTHORIZED**
- Topology, minimization, MD and QM:
  **NOT AUTHORIZED**
- Required next step:
  `{next_step}`
""",
        encoding="utf-8",
    )

    print(
        "Day025 R2 force-field route-comparison "
        "data-contract preflight completed."
    )

    print(
        f"Gate 29 decision: {gate29_decision}"
    )

    print(
        f"Gate 30 decision: {gate30_decision}"
    )

    print(
        f"Gate 37 decision: {gate37_decision}"
    )

    print(
        "Environment / assignments / critical centers: "
        f"{counts['chemical_environment_rows']}/"
        f"{counts['chemical_assignment_rows']}/"
        f"{counts['critical_center_rows']}"
    )

    print(
        "Required bonds / angles / torsions / impropers: "
        f"{counts['required_bond_rows']}/"
        f"{counts['required_angle_rows']}/"
        f"{counts['required_torsion_rows']}/"
        f"{counts['required_improper_rows']}"
    )

    print(
        "QM fragment classes / Lele BNH records / "
        "Lele domain rows: "
        f"{counts['preliminary_qm_fragment_rows']}/"
        f"{counts['lele_bnh_parameter_rows']}/"
        f"{counts['lele_domain_rows']}"
    )

    print(
        f"Decision: {decision}"
    )

    print(
        "Failed gates: "
        + (
            "NONE"
            if not failed_gates
            else " | ".join(failed_gates)
        )
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
        f"Required next step: {next_step}"
    )

    print(
        f"Wrote: {SCHEMA.relative_to(ROOT)}"
    )

    print(
        f"Wrote: {SUMMARY.relative_to(ROOT)}"
    )

    print(
        f"Wrote: {JSON_OUT.relative_to(ROOT)}"
    )

    print(
        f"Wrote: {REPORT.relative_to(ROOT)}"
    )

    print()
    print("INPUT SCHEMAS")

    for row in schema_rows:
        print(
            f"- {row['input_name']}: "
            f"rows={row['row_count']} "
            f"columns={row['column_count']}"
        )

        print(
            f"  fields={row['columns']}"
        )

        print(
            f"  first={row['first_row_preview']}"
        )

    if failed_gates:
        raise RuntimeError(
            "Force-field route-comparison preflight requires review."
        )


if __name__ == "__main__":
    main()
