#!/usr/bin/env python3
"""
DAY038 / D038-E1

Scientific audit of the QM_F06_UPPER_V7A_R1 RESP Stage 1
candidate charges.

The script compares RESP Stage 1 candidate charges against the
ORCA CHELPG charges using the preserved 52-atom order.

Outputs
-------
- Atom-by-atom CSV comparison
- Machine-readable JSON audit report

Scientific policy
-----------------
- Candidate charges are analyzed but not adopted.
- RESP Stage 2 remains blocked.
- Force-field adoption remains blocked.
"""

from __future__ import annotations

import csv
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

from resp_common import (
    load_json,
    require_file,
    sha256,
)


ROOT = Path(__file__).resolve().parents[2]

EXECUTION_PARENT = (
    ROOT
    / "runs/phase1A/day038_resp_stage1_executions"
)

LATEST_POINTER = (
    EXECUTION_PARENT
    / "LATEST_RESP_STAGE1_UPPER_V7A_R1_EXECUTION.txt"
)

DAY036_EXECUTION_POINTER = (
    ROOT
    / "runs/phase1A/day036_qm_f06_upper_v7a_r1_esp_executions"
    / "LATEST_ESP_UPPER_V7A_R1_EXECUTION.txt"
)

ATOM_CLASSES = (
    ROOT
    / "runs/phase1A/day036_qm_f06_upper_v7a_r1_resp_preparation"
    / "QM_F06_UPPER_V7A_R1_RESP_ATOM_CLASSES.csv"
)

EXPECTED_STAGE1_DECISION = (
    "D038_D1_RESP_STAGE1_EXECUTION_PASS_"
    "CANDIDATE_CHARGES_GENERATED_NOT_ADOPTED"
)

EXPECTED_ATOM_COUNT = 52
EXPECTED_COMPOSITION = {
    "B": 17,
    "H": 21,
    "N": 14,
}


def parse_numeric_tokens(path: Path) -> list[float]:
    require_file(path)

    values: list[float] = []

    for token in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).split():
        try:
            value = float(
                token.replace("D", "E").replace("d", "e")
            )
        except ValueError:
            continue

        if not math.isfinite(value):
            raise RuntimeError(
                f"Non-finite value in {path}: {token!r}"
            )

        values.append(value)

    return values


def parse_pc_chelpg(
    path: Path,
) -> tuple[list[float], list[tuple[float, float, float]]]:
    require_file(path)

    lines = path.read_text(
        encoding="utf-8"
    ).splitlines()

    declared_count = int(lines[0].strip())

    charges: list[float] = []
    coordinates: list[tuple[float, float, float]] = []

    for line_number, line in enumerate(
        lines[2:],
        start=3,
    ):
        if not line.strip():
            continue

        tokens = line.split()

        if len(tokens) != 5 or tokens[0].upper() != "Q":
            raise RuntimeError(
                "Malformed .pc_chelpg row.\n"
                f"Line: {line_number}\n"
                f"Raw: {line!r}"
            )

        charge, x, y, z = map(float, tokens[1:5])

        charges.append(charge)
        coordinates.append((x, y, z))

    if len(charges) != declared_count:
        raise RuntimeError(
            ".pc_chelpg count mismatch.\n"
            f"Declared: {declared_count}\n"
            f"Parsed: {len(charges)}"
        )

    return charges, coordinates


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def population_std(values: list[float]) -> float:
    center = mean(values)

    return math.sqrt(
        sum((value - center) ** 2 for value in values)
        / len(values)
    )


def quantile(
    sorted_values: list[float],
    probability: float,
) -> float:
    position = probability * (len(sorted_values) - 1)

    lower = math.floor(position)
    upper = math.ceil(position)

    if lower == upper:
        return sorted_values[lower]

    fraction = position - lower

    return (
        sorted_values[lower] * (1.0 - fraction)
        + sorted_values[upper] * fraction
    )


def regression(
    x_values: list[float],
    y_values: list[float],
) -> dict[str, float]:
    x_mean = mean(x_values)
    y_mean = mean(y_values)

    sxx = sum(
        (value - x_mean) ** 2
        for value in x_values
    )

    syy = sum(
        (value - y_mean) ** 2
        for value in y_values
    )

    sxy = sum(
        (x - x_mean) * (y - y_mean)
        for x, y in zip(x_values, y_values)
    )

    if sxx == 0.0 or syy == 0.0:
        raise RuntimeError(
            "Zero variance in charge regression"
        )

    slope = sxy / sxx
    intercept = y_mean - slope * x_mean
    pearson = sxy / math.sqrt(sxx * syy)

    return {
        "slope": slope,
        "intercept": intercept,
        "pearson_r": pearson,
        "r_squared": pearson * pearson,
    }


def extract_resp_metrics(
    path: Path,
) -> dict[str, float | int | None]:
    require_file(path)

    text = path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    patterns = {
        "convergence_iterations": (
            r"Convergence in\s+(\d+)\s+iterations"
        ),
        "initial_ssvpot": (
            r"The initial sum of squares \(ssvpot\)"
            r"\s+([-+0-9.DEde]+)"
        ),
        "residual_chipot": (
            r"The residual sum of squares \(chipot\)"
            r"\s+([-+0-9.DEde]+)"
        ),
        "standard_error": (
            r"The std err of estimate "
            r"\(sqrt\(chipot/N\)\)"
            r"\s+([-+0-9.DEde]+)"
        ),
        "esp_relative_rms": (
            r"ESP relative RMS "
            r"\(SQRT\(chipot/ssvpot\)\)"
            r"\s+([-+0-9.DEde]+)"
        ),
        "dipole_moment_debye": (
            r"Dipole Moment \(Debye\)="
            r"\s+([-+0-9.DEde]+)"
        ),
    }

    metrics: dict[str, float | int | None] = {}

    for name, pattern in patterns.items():
        match = re.search(pattern, text)

        if match is None:
            metrics[name] = None
            continue

        token = match.group(1)

        if name == "convergence_iterations":
            metrics[name] = int(token)
        else:
            metrics[name] = float(
                token.replace("D", "E").replace("d", "e")
            )

    return metrics


print("=" * 100)
print("DAY038 / D038-E1 — RESP STAGE 1 CANDIDATE CHARGE AUDIT")
print("=" * 100)


print("\n[1] SOURCE EXECUTIONS")

require_file(LATEST_POINTER)
require_file(DAY036_EXECUTION_POINTER)

stage1_execution = (
    ROOT
    / LATEST_POINTER.read_text(
        encoding="utf-8"
    ).strip()
)

day036_execution = (
    ROOT
    / DAY036_EXECUTION_POINTER.read_text(
        encoding="utf-8"
    ).strip()
)

stage1_report_path = (
    stage1_execution
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1.json"
)

stage1_qout = (
    stage1_execution
    / "QM_F06_UPPER_V7A_R1_resp1.qout"
)

stage1_output = (
    stage1_execution
    / "QM_F06_UPPER_V7A_R1_resp1.out"
)

chelpg_path = (
    day036_execution
    / "esp_upper_v7a_r1.pc_chelpg"
)

for path in (
    stage1_report_path,
    stage1_qout,
    stage1_output,
    chelpg_path,
    ATOM_CLASSES,
):
    require_file(path)
    print(f"FOUND  {path}")


print("\n[2] STAGE 1 AUTHORIZATION")

stage1_report = load_json(
    stage1_report_path
)

if (
    stage1_report.get("decision")
    != EXPECTED_STAGE1_DECISION
):
    raise RuntimeError(
        "Unexpected RESP Stage 1 decision.\n"
        f"Expected: {EXPECTED_STAGE1_DECISION}\n"
        f"Observed: {stage1_report.get('decision')}"
    )

authorizations = stage1_report.get(
    "authorizations",
    {},
)

if (
    authorizations.get(
        "candidate_charge_analysis_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Candidate-charge analysis is not authorized"
    )

if authorizations.get("RESP_stage2_authorized") is not False:
    raise RuntimeError(
        "Unexpected RESP Stage 2 authorization state"
    )

if authorizations.get("charge_adoption_authorized") is not False:
    raise RuntimeError(
        "Unexpected charge-adoption authorization state"
    )

print("stage1_decision_gate                 = PASS")
print("candidate_analysis_authorization_gate = PASS")
print("stage2_blocked_gate                  = PASS")
print("charge_adoption_blocked_gate         = PASS")


print("\n[3] ATOM METADATA")

with ATOM_CLASSES.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:
    atom_rows = list(
        csv.DictReader(handle)
    )

if len(atom_rows) != EXPECTED_ATOM_COUNT:
    raise RuntimeError(
        "Atom-class row count mismatch"
    )

indices = [
    int(row["atom_index_0based"])
    for row in atom_rows
]

if indices != list(range(EXPECTED_ATOM_COUNT)):
    raise RuntimeError(
        "Atom indices are not exactly 0..51"
    )

elements = [
    row["element"].strip()
    for row in atom_rows
]

composition = {
    element: elements.count(element)
    for element in sorted(set(elements))
}

if composition != EXPECTED_COMPOSITION:
    raise RuntimeError(
        "Composition mismatch.\n"
        f"Expected: {EXPECTED_COMPOSITION}\n"
        f"Observed: {composition}"
    )

print(f"atom_count  = {len(elements)}")
print(f"composition = {composition}")
print(
    f"metadata_columns = "
    f"{list(atom_rows[0].keys())}"
)
print("atom_metadata_gate = PASS")


print("\n[4] CHARGE DATASETS")

resp_charges = parse_numeric_tokens(
    stage1_qout
)

chelpg_charges, chelpg_coordinates = (
    parse_pc_chelpg(
        chelpg_path
    )
)

if len(resp_charges) != EXPECTED_ATOM_COUNT:
    raise RuntimeError(
        "RESP charge count mismatch.\n"
        f"Observed: {len(resp_charges)}"
    )

if len(chelpg_charges) != EXPECTED_ATOM_COUNT:
    raise RuntimeError(
        "CHELPG charge count mismatch"
    )

resp_sum = sum(resp_charges)
chelpg_sum = sum(chelpg_charges)

print(f"RESP_count       = {len(resp_charges)}")
print(f"RESP_sum_e       = {resp_sum:.16g}")
print(f"RESP_min_e       = {min(resp_charges):.16g}")
print(f"RESP_max_e       = {max(resp_charges):.16g}")
print(f"CHELPG_count     = {len(chelpg_charges)}")
print(f"CHELPG_sum_e     = {chelpg_sum:.16g}")
print(f"CHELPG_min_e     = {min(chelpg_charges):.16g}")
print(f"CHELPG_max_e     = {max(chelpg_charges):.16g}")

neutrality_gate = (
    abs(resp_sum) <= 5.0e-5
    and abs(chelpg_sum) <= 5.0e-5
)

print(
    "neutrality_gate = "
    + ("PASS" if neutrality_gate else "FAIL")
)


print("\n[5] RESP VERSUS CHELPG")

differences = [
    resp - chelpg
    for resp, chelpg
    in zip(resp_charges, chelpg_charges)
]

absolute_differences = [
    abs(value)
    for value in differences
]

rmse = math.sqrt(
    sum(value * value for value in differences)
    / EXPECTED_ATOM_COUNT
)

mae = mean(absolute_differences)
maximum_absolute_difference = max(
    absolute_differences
)

regression_metrics = regression(
    chelpg_charges,
    resp_charges,
)

same_sign_count = sum(
    1
    for chelpg, resp
    in zip(chelpg_charges, resp_charges)
    if (
        (chelpg > 0 and resp > 0)
        or (chelpg < 0 and resp < 0)
        or (chelpg == 0 and resp == 0)
    )
)

sign_change_indices = [
    index
    for index, (chelpg, resp)
    in enumerate(
        zip(chelpg_charges, resp_charges)
    )
    if (
        chelpg != 0.0
        and resp != 0.0
        and math.copysign(1.0, chelpg)
        != math.copysign(1.0, resp)
    )
]

print(f"RMSE_e                    = {rmse:.16g}")
print(f"MAE_e                     = {mae:.16g}")
print(
    f"maximum_absolute_difference_e = "
    f"{maximum_absolute_difference:.16g}"
)
print(
    f"same_sign_count           = "
    f"{same_sign_count}/{EXPECTED_ATOM_COUNT}"
)
print(
    f"sign_change_count         = "
    f"{len(sign_change_indices)}"
)

for name, value in regression_metrics.items():
    print(f"{name} = {value:.16g}")


print("\n[6] ABSOLUTE-DIFFERENCE QUANTILES")

sorted_absolute_differences = sorted(
    absolute_differences
)

difference_quantiles = {}

for probability in (
    0.00,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    1.00,
):
    value = quantile(
        sorted_absolute_differences,
        probability,
    )

    label = f"q{100 * probability:06.2f}"
    difference_quantiles[label] = value
    print(f"{label}_e = {value:.16g}")


print("\n[7] ELEMENT-WISE SUMMARY")

element_summary = {}

for element in sorted(set(elements)):
    selected_indices = [
        index
        for index, observed_element in enumerate(elements)
        if observed_element == element
    ]

    selected_resp = [
        resp_charges[index]
        for index in selected_indices
    ]

    selected_chelpg = [
        chelpg_charges[index]
        for index in selected_indices
    ]

    selected_differences = [
        differences[index]
        for index in selected_indices
    ]

    selected_absolute_differences = [
        abs(value)
        for value in selected_differences
    ]

    summary = {
        "count": len(selected_indices),
        "RESP_mean_e": mean(selected_resp),
        "RESP_std_e": population_std(selected_resp),
        "RESP_min_e": min(selected_resp),
        "RESP_max_e": max(selected_resp),
        "RESP_sum_e": sum(selected_resp),
        "CHELPG_mean_e": mean(selected_chelpg),
        "CHELPG_std_e": population_std(selected_chelpg),
        "CHELPG_min_e": min(selected_chelpg),
        "CHELPG_max_e": max(selected_chelpg),
        "CHELPG_sum_e": sum(selected_chelpg),
        "difference_mean_e": mean(
            selected_differences
        ),
        "difference_MAE_e": mean(
            selected_absolute_differences
        ),
        "difference_RMSE_e": math.sqrt(
            sum(
                value * value
                for value in selected_differences
            )
            / len(selected_differences)
        ),
        "difference_max_abs_e": max(
            selected_absolute_differences
        ),
    }

    element_summary[element] = summary

    print(f"\nElement {element}")

    for name, value in summary.items():
        print(f"  {name} = {value}")


print("\n[8] FIFTEEN LARGEST RESP–CHELPG DIFFERENCES")

ranked_indices = sorted(
    range(EXPECTED_ATOM_COUNT),
    key=lambda index: absolute_differences[index],
    reverse=True,
)

for rank, index in enumerate(
    ranked_indices[:15],
    start=1,
):
    print(
        f"rank={rank:>2} "
        f"atom_0based={index:>2} "
        f"atom_1based={index + 1:>2} "
        f"element={elements[index]} "
        f"CHELPG={chelpg_charges[index]: .6f} "
        f"RESP1={resp_charges[index]: .6f} "
        f"difference={differences[index]: .6f} "
        f"abs_difference={absolute_differences[index]:.6f}"
    )


print("\n[9] SIGN CHANGES")

if sign_change_indices:
    for index in sign_change_indices:
        print(
            f"atom_0based={index:>2} "
            f"atom_1based={index + 1:>2} "
            f"element={elements[index]} "
            f"CHELPG={chelpg_charges[index]: .6f} "
            f"RESP1={resp_charges[index]: .6f}"
        )
else:
    print("(none)")


print("\n[10] RESP FIT METRICS")

fit_metrics = extract_resp_metrics(
    stage1_output
)

for name, value in fit_metrics.items():
    print(f"{name} = {value}")


print("\n[11] WRITE REPRODUCIBLE OUTPUTS")

comparison_csv = (
    stage1_execution
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1_CHARGE_COMPARISON.csv"
)

audit_json = (
    stage1_execution
    / "QM_F06_UPPER_V7A_R1_RESP_STAGE1_AUDIT.json"
)

metadata_columns = [
    column
    for column in atom_rows[0].keys()
    if column not in {
        "atom_index_0based",
        "element",
    }
]

fieldnames = [
    "atom_index_0based",
    "atom_index_1based",
    "element",
    *metadata_columns,
    "x_angstrom",
    "y_angstrom",
    "z_angstrom",
    "CHELPG_charge_e",
    "RESP_stage1_charge_e",
    "RESP_minus_CHELPG_e",
    "absolute_difference_e",
    "sign_changed",
]

with comparison_csv.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
    )

    writer.writeheader()

    for index in range(EXPECTED_ATOM_COUNT):
        x, y, z = chelpg_coordinates[index]

        row = {
            "atom_index_0based": index,
            "atom_index_1based": index + 1,
            "element": elements[index],
            "x_angstrom": f"{x:.9f}",
            "y_angstrom": f"{y:.9f}",
            "z_angstrom": f"{z:.9f}",
            "CHELPG_charge_e": (
                f"{chelpg_charges[index]:.9f}"
            ),
            "RESP_stage1_charge_e": (
                f"{resp_charges[index]:.9f}"
            ),
            "RESP_minus_CHELPG_e": (
                f"{differences[index]:.9f}"
            ),
            "absolute_difference_e": (
                f"{absolute_differences[index]:.9f}"
            ),
            "sign_changed": (
                index in sign_change_indices
            ),
        }

        for column in metadata_columns:
            row[column] = atom_rows[index].get(
                column,
                "",
            )

        writer.writerow(row)


fit_metrics_complete = all(
    value is not None
    for value in fit_metrics.values()
)

charge_finiteness_gate = all(
    math.isfinite(value)
    for value in (
        resp_charges
        + chelpg_charges
        + differences
    )
)

audit_gates = {
    "stage1_decision_gate": True,
    "candidate_analysis_authorization_gate": True,
    "atom_count_gate": (
        len(resp_charges)
        == len(chelpg_charges)
        == EXPECTED_ATOM_COUNT
    ),
    "atom_order_gate": True,
    "composition_gate": (
        composition == EXPECTED_COMPOSITION
    ),
    "neutrality_gate": neutrality_gate,
    "charge_finiteness_gate": charge_finiteness_gate,
    "fit_metrics_parse_gate": fit_metrics_complete,
    "comparison_csv_created_gate": (
        comparison_csv.is_file()
        and comparison_csv.stat().st_size > 0
    ),
}

all_gates_pass = all(
    audit_gates.values()
)

decision = (
    "D038_E1_RESP_STAGE1_CANDIDATE_CHARGES_"
    "SCIENTIFIC_AUDIT_PASS_STAGE2_REMAINS_BLOCKED"
    if all_gates_pass
    else
    "D038_E1_RESP_STAGE1_CANDIDATE_CHARGE_"
    "AUDIT_REVIEW_REQUIRED"
)

audit_report = {
    "generated_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "decision": decision,
    "source_execution_directory": str(
        stage1_execution.resolve()
    ),
    "source_identity": {
        "stage1_report": str(
            stage1_report_path.resolve()
        ),
        "stage1_report_sha256": sha256(
            stage1_report_path
        ),
        "stage1_qout": str(
            stage1_qout.resolve()
        ),
        "stage1_qout_sha256": sha256(
            stage1_qout
        ),
        "stage1_output": str(
            stage1_output.resolve()
        ),
        "stage1_output_sha256": sha256(
            stage1_output
        ),
        "CHELPG_source": str(
            chelpg_path.resolve()
        ),
        "CHELPG_source_sha256": sha256(
            chelpg_path
        ),
        "atom_classes": str(
            ATOM_CLASSES.resolve()
        ),
        "atom_classes_sha256": sha256(
            ATOM_CLASSES
        ),
    },
    "composition": composition,
    "charge_summaries": {
        "RESP_stage1": {
            "count": len(resp_charges),
            "sum_e": resp_sum,
            "mean_e": mean(resp_charges),
            "std_e": population_std(resp_charges),
            "minimum_e": min(resp_charges),
            "maximum_e": max(resp_charges),
        },
        "CHELPG": {
            "count": len(chelpg_charges),
            "sum_e": chelpg_sum,
            "mean_e": mean(chelpg_charges),
            "std_e": population_std(chelpg_charges),
            "minimum_e": min(chelpg_charges),
            "maximum_e": max(chelpg_charges),
        },
    },
    "RESP_vs_CHELPG": {
        "RMSE_e": rmse,
        "MAE_e": mae,
        "maximum_absolute_difference_e": (
            maximum_absolute_difference
        ),
        "same_sign_count": same_sign_count,
        "sign_change_count": len(
            sign_change_indices
        ),
        "sign_change_atom_indices_0based": (
            sign_change_indices
        ),
        "difference_quantiles_e": (
            difference_quantiles
        ),
        "regression": regression_metrics,
    },
    "element_summary": element_summary,
    "RESP_fit_metrics": fit_metrics,
    "outputs": {
        "comparison_csv": str(
            comparison_csv.resolve()
        ),
        "comparison_csv_sha256": sha256(
            comparison_csv
        ),
    },
    "gates": audit_gates,
    "authorizations": {
        "candidate_charge_interpretation_authorized": (
            all_gates_pass
        ),
        "RESP_stage2_authorized": False,
        "charge_adoption_authorized": False,
        "force_field_adoption_authorized": False,
    },
}

audit_json.write_text(
    json.dumps(
        audit_report,
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)

print(f"comparison_csv = {comparison_csv}")
print(
    f"comparison_csv_sha256 = "
    f"{sha256(comparison_csv)}"
)
print(f"audit_json = {audit_json}")
print(f"audit_json_sha256 = {sha256(audit_json)}")


print("\n[12] GATES")

for name, value in audit_gates.items():
    print(
        f"{name} = "
        f"{'PASS' if value else 'FAIL'}"
    )


print("\n[13] DECISION")
print(f"decision = {decision}")
print(
    "candidate_charge_interpretation_authorized = "
    f"{all_gates_pass}"
)
print("RESP_stage2_authorized = False")
print("charge_adoption_authorized = False")
print("force_field_adoption_authorized = False")
print("=" * 100)
