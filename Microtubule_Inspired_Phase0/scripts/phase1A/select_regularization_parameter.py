#!/usr/bin/env python3
"""
DAY039 / D039-A6

Multi-objective scientific review of the regularization path for the
37-real-atom constrained electrostatic refit of
QM_F06_UPPER_V7A_R1.

Objectives minimized
--------------------
1. Electrostatic RMSE against the authorized quantum ESP.
2. RMS perturbation relative to the retained RESP Stage 1 charges.

Additional diagnostic criteria
------------------------------
- exact neutrality;
- maximum absolute charge;
- maximum atom-level perturbation;
- number of atomic sign changes;
- same-sign fraction of the electrostatic potential.

Scientific policy
-----------------
- Pareto dominance is evaluated only on RMSE and charge delta RMS.
- Hard admissibility filters are reported separately.
- A normalized ideal-distance recommendation is generated only for
  scientific review.
- No lambda is adopted.
- No charges are adopted.
- RESP Stage 2 is not executed.
"""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from resp_common import (
    load_json,
    require_file,
    sha256,
)


ROOT = Path(__file__).resolve().parents[2]

LATEST_POINTER = (
    ROOT
    / "runs/phase1A/day038_resp_stage1_executions"
    / "LATEST_RESP_STAGE1_UPPER_V7A_R1_EXECUTION.txt"
)

EXPECTED_A5_DECISION = (
    "D039_A5_CONSTRAINED_REAL_ATOM_REFIT_AUDIT_PASS_"
    "REGULARIZATION_SELECTION_REVIEW_AUTHORIZED"
)

NEUTRALITY_TOLERANCE_E = 1.0e-10

# These are review filters, not universal physical constants.
MAX_ABSOLUTE_CHARGE_REVIEW_LIMIT_E = 1.0
MAX_DELTA_ABSOLUTE_REVIEW_LIMIT_E = 0.85
MAX_SIGN_CHANGE_REVIEW_COUNT = 5

# Objectives used for Pareto analysis.
OBJECTIVE_NAMES = (
    "electrostatic_RMSE_au",
    "charge_delta_RMS_e",
)


def normalize_minimization(
    value: float,
    minimum: float,
    maximum: float,
) -> float:
    if maximum == minimum:
        return 0.0

    return (
        value - minimum
    ) / (
        maximum - minimum
    )


def dominates(
    first: dict,
    second: dict,
) -> bool:
    """
    Return True if first Pareto-dominates second for both
    minimization objectives.
    """

    first_values = (
        first["electrostatic_RMSE_au"],
        first["charge_delta_RMS_e"],
    )

    second_values = (
        second["electrostatic_RMSE_au"],
        second["charge_delta_RMS_e"],
    )

    no_worse = all(
        first_value <= second_value
        for first_value, second_value
        in zip(
            first_values,
            second_values,
        )
    )

    strictly_better = any(
        first_value < second_value
        for first_value, second_value
        in zip(
            first_values,
            second_values,
        )
    )

    return no_worse and strictly_better


def json_safe(value):
    try:
        import numpy as np
    except ImportError:
        np = None

    if np is not None:
        if isinstance(value, np.bool_):
            return bool(value)

        if isinstance(value, np.integer):
            return int(value)

        if isinstance(value, np.floating):
            return float(value)

        if isinstance(value, np.ndarray):
            return value.tolist()

    raise TypeError(
        f"Object of type {type(value).__name__} "
        "is not JSON serializable"
    )


print("=" * 100)
print("DAY039 / D039-A6 — MULTI-OBJECTIVE REGULARIZATION SELECTION REVIEW")
print("=" * 100)


print("\n[1] SOURCE EXECUTION")

require_file(LATEST_POINTER)

execution_dir = (
    ROOT
    / LATEST_POINTER.read_text(
        encoding="utf-8"
    ).strip()
)

a5_report_path = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_CONSTRAINED_REFIT_AUDIT.json"
)

candidate_charge_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_CONSTRAINED_REFIT_REGULARIZATION_PATH.csv"
)

require_file(a5_report_path)
require_file(candidate_charge_csv)

print(f"execution_dir       = {execution_dir}")
print(f"A5_report           = {a5_report_path}")
print(f"candidate_charge_csv = {candidate_charge_csv}")


print("\n[2] UPSTREAM AUTHORIZATION")

a5_report = load_json(
    a5_report_path
)

if (
    a5_report.get("decision")
    != EXPECTED_A5_DECISION
):
    raise RuntimeError(
        "Unexpected D039-A5 decision.\n"
        f"Expected: {EXPECTED_A5_DECISION}\n"
        f"Observed: {a5_report.get('decision')}"
    )

authorizations = a5_report.get(
    "authorizations",
    {},
)

if (
    authorizations.get(
        "regularization_selection_review_authorized"
    )
    is not True
):
    raise RuntimeError(
        "Regularization-selection review is not authorized"
    )

if (
    authorizations.get(
        "constrained_refit_candidate_adoption_authorized"
    )
    is not False
):
    raise RuntimeError(
        "Unexpected constrained-refit adoption authorization"
    )

if (
    authorizations.get(
        "RESP_stage2_execution_authorized"
    )
    is not False
):
    raise RuntimeError(
        "Unexpected RESP Stage 2 authorization"
    )

print("upstream_decision_gate                 = PASS")
print("selection_review_authorized_gate       = PASS")
print("candidate_adoption_blocked_gate        = PASS")
print("RESP_stage2_execution_blocked_gate     = PASS")


print("\n[3] LOAD REGULARIZATION PATH")

regularization_path = a5_report.get(
    "regularization_path",
    [],
)

if not regularization_path:
    raise RuntimeError(
        "A5 regularization path is empty"
    )

records: list[dict] = []

for raw_record in regularization_path:
    electrostatic = raw_record[
        "electrostatic"
    ]

    charges = raw_record[
        "charges"
    ]

    record = {
        "regularization_lambda": float(
            raw_record[
                "regularization_lambda"
            ]
        ),
        "electrostatic_RMSE_au": float(
            electrostatic["RMSE_au"]
        ),
        "electrostatic_MAE_au": float(
            electrostatic["MAE_au"]
        ),
        "electrostatic_max_abs_error_au": float(
            electrostatic[
                "maximum_absolute_error_au"
            ]
        ),
        "electrostatic_pearson_r": float(
            electrostatic["pearson_r"]
        ),
        "electrostatic_same_sign_fraction": float(
            electrostatic[
                "same_sign_fraction"
            ]
        ),
        "charge_sum_e": float(
            charges["charge_sum_e"]
        ),
        "minimum_charge_e": float(
            charges["minimum_charge_e"]
        ),
        "maximum_charge_e": float(
            charges["maximum_charge_e"]
        ),
        "maximum_absolute_charge_e": float(
            charges[
                "maximum_absolute_charge_e"
            ]
        ),
        "charge_delta_MAE_e": float(
            charges["delta_MAE_e"]
        ),
        "charge_delta_RMS_e": float(
            charges["delta_RMS_e"]
        ),
        "charge_delta_max_abs_e": float(
            charges["delta_max_abs_e"]
        ),
        "charge_sign_change_count": int(
            charges["sign_change_count"]
        ),
        "charge_sign_change_real_indices": list(
            charges[
                "sign_change_real_indices"
            ]
        ),
    }

    record["neutrality_gate"] = (
        abs(record["charge_sum_e"])
        <= NEUTRALITY_TOLERANCE_E
    )

    record["maximum_charge_review_gate"] = (
        record["maximum_absolute_charge_e"]
        <= MAX_ABSOLUTE_CHARGE_REVIEW_LIMIT_E
    )

    record["maximum_delta_review_gate"] = (
        record["charge_delta_max_abs_e"]
        <= MAX_DELTA_ABSOLUTE_REVIEW_LIMIT_E
    )

    record["sign_change_review_gate"] = (
        record["charge_sign_change_count"]
        <= MAX_SIGN_CHANGE_REVIEW_COUNT
    )

    record["review_admissible"] = all(
        (
            record["neutrality_gate"],
            record["maximum_charge_review_gate"],
            record["maximum_delta_review_gate"],
            record["sign_change_review_gate"],
        )
    )

    records.append(record)

records.sort(
    key=lambda record: (
        record["regularization_lambda"]
    )
)

print(f"candidate_count = {len(records)}")

for record in records:
    print(
        f"lambda={record['regularization_lambda']:.1e} "
        f"RMSE={record['electrostatic_RMSE_au']:.10g} "
        f"delta_RMS={record['charge_delta_RMS_e']:.10g} "
        f"max_abs_q={record['maximum_absolute_charge_e']:.10g} "
        f"delta_max={record['charge_delta_max_abs_e']:.10g} "
        f"sign_changes={record['charge_sign_change_count']:>2} "
        f"admissible={record['review_admissible']}"
    )


print("\n[4] PARETO ANALYSIS")

for record in records:
    dominators = [
        other
        for other in records
        if (
            other is not record
            and dominates(
                other,
                record,
            )
        )
    ]

    record["pareto_dominated"] = bool(
        dominators
    )

    record["dominated_by_lambdas"] = [
        other["regularization_lambda"]
        for other in dominators
    ]

pareto_records = [
    record
    for record in records
    if not record[
        "pareto_dominated"
    ]
]

pareto_records.sort(
    key=lambda record: (
        record["electrostatic_RMSE_au"]
    )
)

print(
    f"pareto_candidate_count = "
    f"{len(pareto_records)}"
)

for record in pareto_records:
    print(
        f"PARETO lambda="
        f"{record['regularization_lambda']:.1e} "
        f"RMSE="
        f"{record['electrostatic_RMSE_au']:.10g} "
        f"delta_RMS="
        f"{record['charge_delta_RMS_e']:.10g} "
        f"admissible="
        f"{record['review_admissible']}"
    )


print("\n[5] NORMALIZED IDEAL-DISTANCE ANALYSIS")

rmse_values = [
    record["electrostatic_RMSE_au"]
    for record in records
]

delta_rms_values = [
    record["charge_delta_RMS_e"]
    for record in records
]

rmse_min = min(rmse_values)
rmse_max = max(rmse_values)

delta_rms_min = min(
    delta_rms_values
)

delta_rms_max = max(
    delta_rms_values
)

for record in records:
    normalized_rmse = normalize_minimization(
        record["electrostatic_RMSE_au"],
        rmse_min,
        rmse_max,
    )

    normalized_delta_rms = normalize_minimization(
        record["charge_delta_RMS_e"],
        delta_rms_min,
        delta_rms_max,
    )

    record["normalized_electrostatic_RMSE"] = (
        normalized_rmse
    )

    record["normalized_charge_delta_RMS"] = (
        normalized_delta_rms
    )

    record["equal_weight_ideal_distance"] = math.sqrt(
        normalized_rmse ** 2
        + normalized_delta_rms ** 2
    )

admissible_pareto_records = [
    record
    for record in pareto_records
    if record["review_admissible"]
]

if admissible_pareto_records:
    recommendation_pool = (
        admissible_pareto_records
    )

    recommendation_pool_name = (
        "ADMISSIBLE_PARETO_FRONT"
    )
else:
    recommendation_pool = (
        pareto_records
    )

    recommendation_pool_name = (
        "PARETO_FRONT_WITHOUT_ADMISSIBILITY_FILTER"
    )

recommended_record = min(
    recommendation_pool,
    key=lambda record: (
        record[
            "equal_weight_ideal_distance"
        ],
        record["charge_sign_change_count"],
        record["maximum_absolute_charge_e"],
        record["regularization_lambda"],
    ),
)

for record in sorted(
    recommendation_pool,
    key=lambda item: (
        item["equal_weight_ideal_distance"]
    ),
):
    print(
        f"lambda={record['regularization_lambda']:.1e} "
        f"norm_RMSE="
        f"{record['normalized_electrostatic_RMSE']:.8f} "
        f"norm_delta_RMS="
        f"{record['normalized_charge_delta_RMS']:.8f} "
        f"ideal_distance="
        f"{record['equal_weight_ideal_distance']:.8f} "
        f"sign_changes="
        f"{record['charge_sign_change_count']}"
    )

print(
    f"\nrecommendation_pool = "
    f"{recommendation_pool_name}"
)

print(
    f"recommended_lambda_for_review = "
    f"{recommended_record['regularization_lambda']:.16g}"
)


print("\n[6] DISCRETE KNEE ANALYSIS")

# The knee is estimated on the normalized Pareto curve using
# maximum perpendicular distance from the line joining its endpoints.
knee_record = None
knee_distance = None

if len(pareto_records) >= 3:
    ordered_pareto = sorted(
        pareto_records,
        key=lambda record: (
            record[
                "normalized_electrostatic_RMSE"
            ]
        ),
    )

    first = ordered_pareto[0]
    last = ordered_pareto[-1]

    x1 = first[
        "normalized_electrostatic_RMSE"
    ]
    y1 = first[
        "normalized_charge_delta_RMS"
    ]

    x2 = last[
        "normalized_electrostatic_RMSE"
    ]
    y2 = last[
        "normalized_charge_delta_RMS"
    ]

    line_length = math.hypot(
        x2 - x1,
        y2 - y1,
    )

    if line_length > 0.0:
        knee_candidates = []

        for record in ordered_pareto[
            1:-1
        ]:
            x0 = record[
                "normalized_electrostatic_RMSE"
            ]

            y0 = record[
                "normalized_charge_delta_RMS"
            ]

            distance = abs(
                (y2 - y1) * x0
                - (x2 - x1) * y0
                + x2 * y1
                - y2 * x1
            ) / line_length

            record[
                "pareto_endpoint_line_distance"
            ] = distance

            knee_candidates.append(
                (
                    distance,
                    record,
                )
            )

        if knee_candidates:
            knee_distance, knee_record = max(
                knee_candidates,
                key=lambda item: (
                    item[0]
                ),
            )

if knee_record is None:
    print(
        "knee_status = NOT_RESOLVED"
    )
else:
    print(
        "knee_status = RESOLVED"
    )

    print(
        f"knee_lambda = "
        f"{knee_record['regularization_lambda']:.16g}"
    )

    print(
        f"knee_endpoint_line_distance = "
        f"{knee_distance:.16g}"
    )

    print(
        f"knee_review_admissible = "
        f"{knee_record['review_admissible']}"
    )


print("\n[7] RECOMMENDED CANDIDATE DETAILS")

for name, value in (
    recommended_record.items()
):
    print(f"{name} = {value}")


print("\n[8] WRITE OUTPUTS")

selection_csv = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_REGULARIZATION_SELECTION_REVIEW.csv"
)

selection_json = (
    execution_dir
    / "QM_F06_UPPER_V7A_R1_REGULARIZATION_SELECTION_REVIEW.json"
)

csv_fieldnames = [
    "regularization_lambda",
    "electrostatic_RMSE_au",
    "electrostatic_MAE_au",
    "electrostatic_max_abs_error_au",
    "electrostatic_pearson_r",
    "electrostatic_same_sign_fraction",
    "charge_sum_e",
    "minimum_charge_e",
    "maximum_charge_e",
    "maximum_absolute_charge_e",
    "charge_delta_MAE_e",
    "charge_delta_RMS_e",
    "charge_delta_max_abs_e",
    "charge_sign_change_count",
    "charge_sign_change_real_indices",
    "neutrality_gate",
    "maximum_charge_review_gate",
    "maximum_delta_review_gate",
    "sign_change_review_gate",
    "review_admissible",
    "pareto_dominated",
    "dominated_by_lambdas",
    "normalized_electrostatic_RMSE",
    "normalized_charge_delta_RMS",
    "equal_weight_ideal_distance",
    "pareto_endpoint_line_distance",
    "recommended_for_scientific_review",
    "discrete_knee_candidate",
]

with selection_csv.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=csv_fieldnames,
    )

    writer.writeheader()

    for record in records:
        output = dict(record)

        output[
            "charge_sign_change_real_indices"
        ] = json.dumps(
            output[
                "charge_sign_change_real_indices"
            ]
        )

        output[
            "dominated_by_lambdas"
        ] = json.dumps(
            output[
                "dominated_by_lambdas"
            ]
        )

        output[
            "pareto_endpoint_line_distance"
        ] = record.get(
            "pareto_endpoint_line_distance"
        )

        output[
            "recommended_for_scientific_review"
        ] = (
            record is recommended_record
        )

        output[
            "discrete_knee_candidate"
        ] = (
            knee_record is not None
            and record is knee_record
        )

        writer.writerow(
            {
                field: output.get(
                    field
                )
                for field in csv_fieldnames
            }
        )


print("\n[9] SCIENTIFIC GATES")

all_neutral_gate = all(
    record["neutrality_gate"]
    for record in records
)

pareto_nonempty_gate = (
    len(pareto_records) > 0
)

recommendation_resolved_gate = (
    recommended_record is not None
)

recommended_is_pareto_gate = (
    recommended_record
    in pareto_records
)

recommended_admissibility_gate = (
    recommended_record[
        "review_admissible"
    ]
)

gates = {
    "upstream_decision_gate": True,
    "candidate_count_gate": (
        len(records) == 12
    ),
    "all_candidate_neutrality_gate": (
        all_neutral_gate
    ),
    "pareto_front_nonempty_gate": (
        pareto_nonempty_gate
    ),
    "recommendation_resolved_gate": (
        recommendation_resolved_gate
    ),
    "recommended_is_pareto_gate": (
        recommended_is_pareto_gate
    ),
    "recommended_review_admissibility_gate": (
        recommended_admissibility_gate
    ),
    "selection_csv_created_gate": (
        selection_csv.is_file()
        and selection_csv.stat().st_size > 0
    ),
    "no_lambda_adopted_gate": True,
    "no_charge_adopted_gate": True,
    "RESP_stage2_not_executed_gate": True,
}

for name, value in gates.items():
    print(
        f"{name} = "
        f"{'PASS' if value else 'FAIL'}"
    )

all_gates_pass = all(
    gates.values()
)


print("\n[10] WRITE JSON REPORT")

decision = (
    "D039_A6_REGULARIZATION_SELECTION_REVIEW_PASS_"
    "CANDIDATE_LAMBDA_PROPOSED_NOT_ADOPTED"
    if all_gates_pass
    else
    "D039_A6_REGULARIZATION_SELECTION_REVIEW_"
    "REVIEW_REQUIRED"
)

report = {
    "generated_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "decision": decision,
    "source_execution_directory": str(
        execution_dir.resolve()
    ),
    "source_identity": {
        "A5_report": str(
            a5_report_path.resolve()
        ),
        "A5_report_sha256": sha256(
            a5_report_path
        ),
        "candidate_charge_csv": str(
            candidate_charge_csv.resolve()
        ),
        "candidate_charge_csv_sha256": sha256(
            candidate_charge_csv
        ),
    },
    "objective_definition": {
        "pareto_objectives_minimized": list(
            OBJECTIVE_NAMES
        ),
        "normalization": (
            "MIN_MAX_ACROSS_EVALUATED_REGULARIZATION_PATH"
        ),
        "ideal_distance": (
            "EQUAL_WEIGHT_EUCLIDEAN_DISTANCE_IN_"
            "NORMALIZED_OBJECTIVE_SPACE"
        ),
        "knee_method": (
            "MAXIMUM_DISTANCE_FROM_PARETO_ENDPOINT_LINE"
        ),
    },
    "review_thresholds": {
        "neutrality_tolerance_e": (
            NEUTRALITY_TOLERANCE_E
        ),
        "maximum_absolute_charge_limit_e": (
            MAX_ABSOLUTE_CHARGE_REVIEW_LIMIT_E
        ),
        "maximum_delta_absolute_limit_e": (
            MAX_DELTA_ABSOLUTE_REVIEW_LIMIT_E
        ),
        "maximum_sign_change_count": (
            MAX_SIGN_CHANGE_REVIEW_COUNT
        ),
        "threshold_policy": (
            "PROJECT_REVIEW_FILTERS_NOT_UNIVERSAL_"
            "PHYSICAL_CONSTANTS"
        ),
    },
    "candidate_records": records,
    "pareto_lambdas": [
        record[
            "regularization_lambda"
        ]
        for record in pareto_records
    ],
    "admissible_pareto_lambdas": [
        record[
            "regularization_lambda"
        ]
        for record in admissible_pareto_records
    ],
    "recommendation_pool": (
        recommendation_pool_name
    ),
    "recommended_candidate_for_review": (
        recommended_record
    ),
    "discrete_knee_candidate": (
        knee_record
    ),
    "discrete_knee_distance": (
        knee_distance
    ),
    "gates": gates,
    "authorizations": {
        "recommended_lambda_scientific_review_authorized": (
            all_gates_pass
        ),
        "regularization_lambda_adoption_authorized": False,
        "constrained_refit_charge_adoption_authorized": False,
        "RESP_stage2_protocol_design_authorized": False,
        "RESP_stage2_execution_authorized": False,
        "charge_adoption_authorized": False,
        "force_field_adoption_authorized": False,
    },
    "outputs": {
        "selection_csv": str(
            selection_csv.resolve()
        ),
        "selection_csv_sha256": sha256(
            selection_csv
        ),
    },
}

selection_json.write_text(
    json.dumps(
        report,
        indent=2,
        default=json_safe,
    )
    + "\n",
    encoding="utf-8",
)

print(f"selection_csv = {selection_csv}")
print(
    f"selection_csv_sha256 = "
    f"{sha256(selection_csv)}"
)
print(f"selection_json = {selection_json}")
print(
    f"selection_json_sha256 = "
    f"{sha256(selection_json)}"
)


print("\n[11] DECISION")

print(f"decision = {decision}")

print(
    f"recommended_lambda_for_scientific_review = "
    f"{recommended_record['regularization_lambda']:.16g}"
)

print(
    "recommended_lambda_scientific_review_authorized = "
    f"{all_gates_pass}"
)

print(
    "regularization_lambda_adoption_authorized = False"
)

print(
    "constrained_refit_charge_adoption_authorized = False"
)

print(
    "RESP_stage2_protocol_design_authorized = False"
)

print(
    "RESP_stage2_execution_authorized = False"
)

print("charge_adoption_authorized = False")
print("force_field_adoption_authorized = False")
print("=" * 100)
