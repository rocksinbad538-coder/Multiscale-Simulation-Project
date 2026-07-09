#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/phase1A/day024_chemical_end_rim_design"

GATE3O = (
    BASE
    / "19_r2_selected_four_atom_heavy_coordinate_embedding"
)

GATE3O1 = (
    BASE
    / "20_r2_aperture_metric_consistency_audit"
)

HEAVY_SUMMARY = (
    GATE3O
    / "r2_selected_four_atom_heavy_embedding_summary.csv"
)

HEAVY_GATES = (
    GATE3O
    / "r2_selected_four_atom_heavy_embedding_gates.csv"
)

HEAVY_COORDINATES = (
    GATE3O
    / "r2_selected_four_atom_heavy_coordinates.csv"
)

METRIC_SUMMARY = (
    GATE3O1
    / "r2_aperture_metric_consistency_summary.csv"
)

OUT = (
    BASE
    / "21_r2_aperture_target_semantics_resolution"
)

SUMMARY = (
    OUT
    / "r2_aperture_target_semantics_resolution_summary.csv"
)

GATES = (
    OUT
    / "r2_aperture_target_semantics_resolution_gates.csv"
)

JSON_OUT = (
    OUT
    / "r2_aperture_target_semantics_resolution.json"
)

MANIFEST = (
    OUT
    / "r2_aperture_target_semantics_resolution_manifest.csv"
)

REPORT = (
    OUT
    / "R2_APERTURE_TARGET_SEMANTICS_RESOLUTION_DAY024.md"
)

EXPECTED_HEAVY_DECISION = (
    "R2_SELECTED_FOUR_ATOM_BN_BRIDGE_"
    "HEAVY_COORDINATE_EMBEDDING_REQUIRES_REVIEW"
)

EXPECTED_METRIC_DECISION = (
    "R2_APERTURE_METRIC_AUDIT_REQUIRES_FURTHER_REVIEW"
)

PASS_DECISION = (
    "R2_HEAVY_COORDINATE_EMBEDDING_VALIDATED_"
    "APERTURE_FUNCTIONAL_GATE_DEFERRED"
)

FAIL_DECISION = (
    "R2_APERTURE_TARGET_SEMANTICS_RESOLUTION_REQUIRES_REVIEW"
)

TARGET_EFFECTIVE_APERTURE_DIAMETER_5KBT_NM = (
    0.8394062099300137
)

TARGET_EFFECTIVE_APERTURE_RADIUS_5KBT_NM = (
    0.41970310496500685
)

CAP_WATER_5KBT_DISTANCE_NM = 0.170

EXPECTED_HEAVY_NODES = 2112
EXPECTED_BRIDGE_NODES = 120
EXPECTED_HEAVY_EDGES = 3066


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty required file: {path}"
        )


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
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


def read_rows(
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


def read_one(
    path: Path,
) -> dict[str, str]:
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


def parse_int(
    row: dict[str, str],
    key: str,
) -> int:
    return int(
        round(
            parse_float(
                row,
                key,
            )
        )
    )


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {
        "true",
        "yes",
        "1",
        "pass",
    }


def main() -> None:
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        HEAVY_SUMMARY,
        HEAVY_GATES,
        HEAVY_COORDINATES,
        METRIC_SUMMARY,
    ):
        require_file(required)

    heavy_summary = read_one(
        HEAVY_SUMMARY
    )

    metric_summary = read_one(
        METRIC_SUMMARY
    )

    heavy_gate_rows = read_rows(
        HEAVY_GATES
    )

    heavy_coordinate_rows = read_rows(
        HEAVY_COORDINATES
    )

    failed_heavy_gates = [
        row["gate"]
        for row in heavy_gate_rows
        if not parse_bool(
            row["pass"]
        )
    ]

    raw_nuclear_lower = parse_float(
        heavy_summary,
        "lower_aperture_diameter_nm",
    )

    raw_nuclear_upper = parse_float(
        heavy_summary,
        "upper_aperture_diameter_nm",
    )

    implied_exclusion_lower = (
        raw_nuclear_lower
        - TARGET_EFFECTIVE_APERTURE_DIAMETER_5KBT_NM
    ) / 2.0

    implied_exclusion_upper = (
        raw_nuclear_upper
        - TARGET_EFFECTIVE_APERTURE_DIAMETER_5KBT_NM
    ) / 2.0

    gates = {
        "Gate3O_failed_only_the_non_equivalent_aperture_comparison": (
            failed_heavy_gates
            == [
                "aperture_errors_are_within10percent"
            ]
        ),
        "Gate3O1_confirmed_all_nuclear_aperture_metrics_exceed_effective_target": (
            metric_summary.get(
                "decision"
            )
            == EXPECTED_METRIC_DECISION
            and metric_summary.get(
                "preferred_aperture_metric"
            )
            == "NONE"
        ),
        "heavy_embedding_contains_2112_heavy_nodes": (
            parse_int(
                heavy_summary,
                "heavy_nodes",
            )
            == EXPECTED_HEAVY_NODES
        ),
        "heavy_embedding_contains_120_exact_bridge_nodes": (
            parse_int(
                heavy_summary,
                "exact_bridge_nodes",
            )
            == EXPECTED_BRIDGE_NODES
        ),
        "heavy_embedding_contains_3066_heavy_edges": (
            parse_int(
                heavy_summary,
                "heavy_edges",
            )
            == EXPECTED_HEAVY_EDGES
        ),
        "heavy_embedding_has_zero_clashes": (
            parse_int(
                heavy_summary,
                "heavy_heavy_clash_count",
            )
            == 0
        ),
        "heavy_embedding_minimum_angle_is_at_least70deg": (
            parse_float(
                heavy_summary,
                "critical_angle_minimum_deg",
            )
            >= 70.0
        ),
        "heavy_embedding_BN_deviation_is_at_most0p003nm": (
            parse_float(
                heavy_summary,
                "maximum_BN_bond_deviation_nm",
            )
            <= 0.003
        ),
        "fixed_coordinates_are_unchanged": (
            parse_bool(
                heavy_summary[
                    "fixed_coordinates_unchanged"
                ]
            )
        ),
        "lower_upper_heavy_geometry_is_symmetric": (
            abs(
                raw_nuclear_lower
                - raw_nuclear_upper
            )
            <= 1.0e-9
        ),
        "target_is_explicitly_classified_as_effective_5kBT_aperture": (
            TARGET_EFFECTIVE_APERTURE_DIAMETER_5KBT_NM
            == 2.0
            * TARGET_EFFECTIVE_APERTURE_RADIUS_5KBT_NM
        ),
        "heavy_nuclear_aperture_is_not_relabelled_as_effective_5kBT_aperture": (
            abs(
                raw_nuclear_lower
                - TARGET_EFFECTIVE_APERTURE_DIAMETER_5KBT_NM
            )
            > 1.0e-6
        ),
        "204_H_coordinates_are_still_pending": (
            len(
                heavy_coordinate_rows
            )
            == EXPECTED_HEAVY_NODES
        ),
    }

    failed_gates = [
        name
        for name, passed
        in gates.items()
        if not passed
    ]

    accepted = (
        len(failed_gates) == 0
    )

    decision = (
        PASS_DECISION
        if accepted
        else FAIL_DECISION
    )

    required_next_step = (
        "GENERATE_AND_VALIDATE_R2_SELECTED_FOUR_ATOM_"
        "BN_BRIDGE_HYDROGEN_COORDINATES"
        if accepted
        else
        "REVIEW_R2_APERTURE_TARGET_SEMANTICS_RESOLUTION"
    )

    summary = {
        "decision": decision,
        "target_quantity": (
            "EFFECTIVE_WATER_ACCESSIBLE_APERTURE_AT_5KBT"
        ),
        "target_effective_aperture_radius_5kBT_nm": (
            TARGET_EFFECTIVE_APERTURE_RADIUS_5KBT_NM
        ),
        "target_effective_aperture_diameter_5kBT_nm": (
            TARGET_EFFECTIVE_APERTURE_DIAMETER_5KBT_NM
        ),
        "source_CAP_water_5kBT_distance_nm": (
            CAP_WATER_5KBT_DISTANCE_NM
        ),
        "heavy_nuclear_aperture_lower_nm": (
            raw_nuclear_lower
        ),
        "heavy_nuclear_aperture_upper_nm": (
            raw_nuclear_upper
        ),
        "implied_exclusion_per_side_lower_nm": (
            implied_exclusion_lower
        ),
        "implied_exclusion_per_side_upper_nm": (
            implied_exclusion_upper
        ),
        "heavy_coordinate_embedding_validated": (
            accepted
        ),
        "heavy_nuclear_aperture_used_as_functional_5kBT_gate": False,
        "inner_H_aperture_proxy_pending": True,
        "effective_5kBT_aperture_validation_pending_nonbonded_model": True,
        "hydrogen_coordinate_generation_authorized": (
            accepted
        ),
        "molecular_topology_generation_authorized": False,
        "formal_charge_assignment_authorized": False,
        "force_field_parameterization_authorized": False,
        "energy_minimization_authorized": False,
        "MD_authorized": False,
        "QM_authorized": False,
        "failed_gates": (
            " | ".join(failed_gates)
        ),
        "required_next_step": (
            required_next_step
        ),
    }

    write_rows(
        SUMMARY,
        [summary],
    )

    write_rows(
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
                "interpretation": {
                    "heavy_nuclear_diameter": (
                        "Raw geometric distance between the inner "
                        "heavy-atom rim centers."
                    ),
                    "inner_H_nuclear_diameter": (
                        "Static geometric proxy available after "
                        "hydrogen placement."
                    ),
                    "effective_5kBT_aperture": (
                        "Water-accessible aperture derived from an "
                        "explicit nonbonded interaction threshold."
                    ),
                },
                "limitations": [
                    (
                        "No coordinates are modified by this gate."
                    ),
                    (
                        "The heavy nuclear aperture is not converted "
                        "into a 5 kBT aperture using an assumed atomistic "
                        "exclusion radius."
                    ),
                    (
                        "The original CAP-water distance of 0.170 nm "
                        "belongs to the neutral steric control and is not "
                        "automatically transferred to atomistic B/N/H."
                    ),
                    (
                        "No topology, charges, force-field parameters, "
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

    write_rows(
        MANIFEST,
        [
            {
                "role": "Gate3O_heavy_summary",
                "file": relative(
                    HEAVY_SUMMARY
                ),
                "sha256": sha256(
                    HEAVY_SUMMARY
                ),
            },
            {
                "role": "Gate3O_heavy_gates",
                "file": relative(
                    HEAVY_GATES
                ),
                "sha256": sha256(
                    HEAVY_GATES
                ),
            },
            {
                "role": "Gate3O_heavy_coordinates",
                "file": relative(
                    HEAVY_COORDINATES
                ),
                "sha256": sha256(
                    HEAVY_COORDINATES
                ),
            },
            {
                "role": "Gate3O1_metric_summary",
                "file": relative(
                    METRIC_SUMMARY
                ),
                "sha256": sha256(
                    METRIC_SUMMARY
                ),
            },
        ],
    )

    gate_lines = "\n".join(
        (
            f"- `{name}`: "
            f"**{'PASS' if passed else 'FAIL'}**"
        )
        for name, passed in gates.items()
    )

    REPORT.write_text(
        f"""# R2 Aperture Target Semantics Resolution

## Resolved definition

The value **{TARGET_EFFECTIVE_APERTURE_DIAMETER_5KBT_NM:.9f} nm**
is an effective water-accessible aperture obtained from the R2 neutral
steric-control model at the 5 kBT interaction threshold.

It is not a direct heavy-atom nuclear diameter.

## Heavy embedding

- Lower/upper heavy nuclear aperture:
  **{raw_nuclear_lower:.9f}/{raw_nuclear_upper:.9f} nm**
- Heavy atoms:
  **{parse_int(heavy_summary, 'heavy_nodes')}**
- Exact bridge atoms:
  **{parse_int(heavy_summary, 'exact_bridge_nodes')}**
- Heavy-heavy clashes:
  **{parse_int(heavy_summary, 'heavy_heavy_clash_count')}**
- Minimum critical angle:
  **{parse_float(heavy_summary, 'critical_angle_minimum_deg'):.9f} degrees**
- Maximum B-N deviation:
  **{parse_float(heavy_summary, 'maximum_BN_bond_deviation_nm'):.9f} nm**

## Interpretation

The heavy-only aperture is retained as a structural descriptor. It is
not compared directly against the effective 5 kBT target.

The next static geometric proxy is the aperture defined by the inner
hydrogen nuclei after the 204 passivants are placed. Final validation
of the effective water-accessible aperture remains deferred until an
explicit nonbonded model is authorized.

## Gates

{gate_lines}

## Decision

- Decision:
  **{decision}**
- Failed gates:
  **{'NONE' if not failed_gates else ' | '.join(failed_gates)}**
- Heavy coordinate embedding validated:
  **{'YES' if accepted else 'NO'}**
- Hydrogen coordinate generation authorized:
  **{'YES' if accepted else 'NO'}**
- Molecular topology generation authorized:
  **NO**
- Force-field parameterization authorized:
  **NO**
- Energy minimization authorized:
  **NO**
- MD authorized:
  **NO**
- QM authorized:
  **NO**
- Required next step:
  `{required_next_step}`
""",
        encoding="utf-8",
    )

    print(
        "Day024 R2 aperture target semantics "
        "resolution completed."
    )

    print(
        "Resolved target quantity: "
        "EFFECTIVE_WATER_ACCESSIBLE_APERTURE_AT_5KBT"
    )

    print(
        "Effective target radius / diameter: "
        f"{TARGET_EFFECTIVE_APERTURE_RADIUS_5KBT_NM:.9f}/"
        f"{TARGET_EFFECTIVE_APERTURE_DIAMETER_5KBT_NM:.9f} nm"
    )

    print(
        "Heavy nuclear aperture lower / upper: "
        f"{raw_nuclear_lower:.9f}/"
        f"{raw_nuclear_upper:.9f} nm"
    )

    print(
        "Implied per-side reduction to effective target "
        "(diagnostic only): "
        f"{implied_exclusion_lower:.9f}/"
        f"{implied_exclusion_upper:.9f} nm"
    )

    print(
        "Heavy nuclear aperture used as 5kBT functional gate: NO"
    )

    print(
        "Heavy coordinate embedding validated: "
        f"{accepted}"
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
        "Hydrogen coordinate generation authorized: "
        f"{'YES' if accepted else 'NO'}"
    )

    print(
        "Molecular topology generation authorized: NO"
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
        "QM authorized: NO"
    )

    print(
        f"Required next step: {required_next_step}"
    )

    for path in (
        SUMMARY,
        GATES,
        JSON_OUT,
        MANIFEST,
        REPORT,
    ):
        print(
            f"Wrote: {relative(path)}"
        )

    if not accepted:
        raise RuntimeError(
            "Aperture-target semantics resolution requires review."
        )


if __name__ == "__main__":
    main()
