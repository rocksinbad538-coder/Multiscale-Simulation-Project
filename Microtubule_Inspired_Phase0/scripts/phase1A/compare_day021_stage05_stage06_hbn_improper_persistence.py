#!/usr/bin/env python3

from __future__ import annotations

import csv
import importlib.util
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

PROTOCOL = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

MODULE_PATH = (
    ROOT
    / "scripts/phase1A/"
    "diagnose_day021_stage06_hbn_improper_phase.py"
)

GRO_PATHS = {
    "stage04": (
        PROTOCOL
        / "execution/04_nvt_k100_2ps/"
        "04_nvt_k100_2ps.gro"
    ),
    "stage05": (
        PROTOCOL
        / "execution/05_nvt_unrestrained_2ps/"
        "05_nvt_unrestrained_2ps.gro"
    ),
    "stage06": (
        PROTOCOL
        / "execution/06_nvt_unrestrained_10ps/"
        "06_nvt_unrestrained_10ps.gro"
    ),
}

HBN_ITP = (
    PROTOCOL
    / "protocol_inputs/topology/"
    "hbn_bonded_mobile_release.itp"
)

RUN06 = (
    PROTOCOL
    / "execution/06_nvt_unrestrained_10ps"
)

SUMMARY_CSV = (
    RUN06
    / "stage06_hbn_improper_persistence_summary.csv"
)

OUTLIERS_CSV = (
    RUN06
    / "stage06_hbn_improper_persistent_outliers.csv"
)

REPORT_MD = (
    RUN06
    / "STAGE06_HBN_IMPROPER_PERSISTENCE_DAY021.md"
)


def load_diagnostic_module():
    specification = (
        importlib.util.spec_from_file_location(
            "stage06_improper_module",
            MODULE_PATH,
        )
    )

    if (
        specification is None
        or specification.loader is None
    ):
        raise RuntimeError(
            "Could not load Stage06 improper module"
        )

    module = importlib.util.module_from_spec(
        specification
    )

    specification.loader.exec_module(module)

    return module


def q(values: np.ndarray, probability: float) -> float:
    return float(
        np.quantile(
            values,
            probability,
        )
    )


def metrics(values: np.ndarray) -> dict[str, float]:
    return {
        "q95": q(values, 0.95),
        "q99": q(values, 0.99),
        "maximum": float(values.max()),
        "mean": float(values.mean()),
    }


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
) -> None:
    if not rows:
        raise RuntimeError(
            f"No rows available for {path}"
        )

    fields: list[str] = []
    seen: set[str] = set()

    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
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
    required = (
        MODULE_PATH,
        HBN_ITP,
        *GRO_PATHS.values(),
    )

    missing = [
        path
        for path in required
        if (
            not path.exists()
            or path.stat().st_size == 0
        )
    ]

    if missing:
        raise RuntimeError(
            "Missing required files:\n"
            + "\n".join(
                str(path)
                for path in missing
            )
        )

    module = load_diagnostic_module()

    impropers = module.parse_impropers(
        HBN_ITP
    )

    angles: dict[str, np.ndarray] = {}
    planarity: dict[str, np.ndarray] = {}
    equilibrium: dict[str, np.ndarray] = {}
    boxes = {}

    for label, path in GRO_PATHS.items():
        positions, box = module.read_hbn_gro(
            path
        )

        boxes[label] = box

        angles[label] = np.array(
            [
                module.dihedral_angle(
                    positions,
                    box,
                    improper,
                )
                for improper in impropers
            ],
            dtype=float,
        )

        planarity[label] = (
            module.planarity_deviation(
                angles[label]
            )
        )

        transformed = module.wrap_degrees(
            -angles[label] + 180.0
        )

        equilibrium[label] = (
            module.equilibrium_residual(
                transformed,
                impropers,
            )
        )

    if not (
        np.allclose(
            boxes["stage04"],
            boxes["stage05"],
            atol=1.0e-8,
        )
        and np.allclose(
            boxes["stage05"],
            boxes["stage06"],
            atol=1.0e-8,
        )
    ):
        raise RuntimeError(
            "Simulation boxes do not match"
        )

    change_05_06 = np.abs(
        module.wrap_degrees(
            angles["stage06"]
            - angles["stage05"]
        )
    )

    metrics05 = metrics(
        planarity["stage05"]
    )

    metrics06 = metrics(
        planarity["stage06"]
    )

    equilibrium06 = metrics(
        equilibrium["stage06"]
    )

    change_metrics = metrics(
        change_05_06
    )

    count = len(impropers)
    top_count = max(
        1,
        math.ceil(
            0.01 * count
        ),
    )

    top05 = set(
        np.argsort(
            planarity["stage05"]
        )[::-1][:top_count].tolist()
    )

    top06 = set(
        np.argsort(
            planarity["stage06"]
        )[::-1][:top_count].tolist()
    )

    above20_stage05 = set(
        np.flatnonzero(
            planarity["stage05"] > 20.0
        ).tolist()
    )

    above20_stage06 = set(
        np.flatnonzero(
            planarity["stage06"] > 20.0
        ).tolist()
    )

    persistent_above20 = (
        above20_stage05
        & above20_stage06
    )

    newly_above20 = (
        above20_stage06
        - above20_stage05
    )

    recovered_below20 = (
        above20_stage05
        - above20_stage06
    )

    blocked_reasons: list[str] = []
    review_reasons: list[str] = []

    if metrics06["maximum"] > 60.0:
        blocked_reasons.append(
            "Stage06 planarity maximum exceeds 60 degrees"
        )

    if equilibrium06["maximum"] > 60.0:
        blocked_reasons.append(
            "Stage06 calibrated equilibrium maximum "
            "exceeds 60 degrees"
        )

    if change_metrics["maximum"] > 60.0:
        blocked_reasons.append(
            "Stage05-to-Stage06 improper change "
            "exceeds 60 degrees"
        )

    relative_checks = {
        "Stage06 planarity q99 <= 22.5 degrees": (
            metrics06["q99"] <= 22.5
        ),
        "Stage06 planarity maximum <= 35 degrees": (
            metrics06["maximum"] <= 35.0
        ),
        "planarity q99 increase <= 2.5 degrees": (
            metrics06["q99"]
            - metrics05["q99"]
            <= 2.5
        ),
        "planarity maximum increase <= 5 degrees": (
            metrics06["maximum"]
            - metrics05["maximum"]
            <= 5.0
        ),
        "calibrated equilibrium q99 <= 20 degrees": (
            equilibrium06["q99"] <= 20.0
        ),
        "calibrated equilibrium maximum <= 40 degrees": (
            equilibrium06["maximum"] <= 40.0
        ),
        "stage-change q99 <= 30 degrees": (
            change_metrics["q99"] <= 30.0
        ),
        "stage-change maximum <= 45 degrees": (
            change_metrics["maximum"] <= 45.0
        ),
    }

    review_reasons.extend(
        label
        for label, passed
        in relative_checks.items()
        if not passed
    )

    if blocked_reasons:
        decision = "BLOCKED"
        next_step = "REVIEW_HBN_MODEL"
    elif review_reasons:
        decision = "REVIEW"
        next_step = "ADDITIONAL_SHORT_DIAGNOSTIC"
    else:
        decision = "PASS_WITH_MONITORING"
        next_step = (
            "25PS_EXTENDED_UNRESTRAINED_VALIDATION"
        )

    summary = {
        "improper_count": count,
        "top_one_percent_count": top_count,
        "stage05_planarity_q99_deg": (
            metrics05["q99"]
        ),
        "stage05_planarity_max_deg": (
            metrics05["maximum"]
        ),
        "stage06_planarity_q99_deg": (
            metrics06["q99"]
        ),
        "stage06_planarity_max_deg": (
            metrics06["maximum"]
        ),
        "planarity_q99_increase_deg": (
            metrics06["q99"]
            - metrics05["q99"]
        ),
        "planarity_max_increase_deg": (
            metrics06["maximum"]
            - metrics05["maximum"]
        ),
        "stage06_equilibrium_q99_deg": (
            equilibrium06["q99"]
        ),
        "stage06_equilibrium_max_deg": (
            equilibrium06["maximum"]
        ),
        "stage05_stage06_change_q99_deg": (
            change_metrics["q99"]
        ),
        "stage05_stage06_change_max_deg": (
            change_metrics["maximum"]
        ),
        "top_one_percent_overlap": (
            len(top05 & top06)
        ),
        "stage05_above_20_deg": (
            len(above20_stage05)
        ),
        "stage06_above_20_deg": (
            len(above20_stage06)
        ),
        "persistent_above_20_deg": (
            len(persistent_above20)
        ),
        "newly_above_20_deg": (
            len(newly_above20)
        ),
        "recovered_below_20_deg": (
            len(recovered_below20)
        ),
        "targeted_decision": decision,
        "authorized_next_step": next_step,
        "long_mobile_production_authorized": False,
        "review_reasons": (
            " | ".join(
                review_reasons
            )
        ),
        "blocked_reasons": (
            " | ".join(
                blocked_reasons
            )
        ),
    }

    outlier_indices = sorted(
        above20_stage05
        | above20_stage06,
        key=lambda index: (
            planarity["stage06"][index]
        ),
        reverse=True,
    )

    outlier_rows = []

    for index in outlier_indices:
        improper = impropers[index]

        outlier_rows.append(
            {
                "improper_index": index + 1,
                "atom_i": improper["i"],
                "atom_j": improper["j"],
                "atom_k": improper["k"],
                "atom_l": improper["l"],
                "stage04_planarity_deg": float(
                    planarity["stage04"][index]
                ),
                "stage05_planarity_deg": float(
                    planarity["stage05"][index]
                ),
                "stage06_planarity_deg": float(
                    planarity["stage06"][index]
                ),
                "stage05_stage06_change_deg": float(
                    change_05_06[index]
                ),
                "stage06_equilibrium_residual_deg": float(
                    equilibrium["stage06"][index]
                ),
                "persistent_above_20_deg": (
                    index
                    in persistent_above20
                ),
            }
        )

    write_csv(
        SUMMARY_CSV,
        [summary],
    )

    if outlier_rows:
        write_csv(
            OUTLIERS_CSV,
            outlier_rows,
        )

    REPORT_MD.write_text(
        f"""# Day021 Stage05–Stage06 HBN Improper Persistence

- Targeted decision: **{decision}**
- Authorized next step: `{next_step}`
- Long mobile production authorized: **NO**

## Relative planarity

- Stage05 q99/max: {metrics05['q99']:.4f}/{metrics05['maximum']:.4f} degrees
- Stage06 q99/max: {metrics06['q99']:.4f}/{metrics06['maximum']:.4f} degrees
- q99/max increase: {summary['planarity_q99_increase_deg']:.4f}/{summary['planarity_max_increase_deg']:.4f} degrees
- Calibrated equilibrium q99/max: {equilibrium06['q99']:.4f}/{equilibrium06['maximum']:.4f} degrees
- Stage05–Stage06 change q99/max: {change_metrics['q99']:.4f}/{change_metrics['maximum']:.4f} degrees

## Persistence

- Top 1% overlap: {len(top05 & top06)}/{top_count}
- Stage05 impropers above 20 degrees: {len(above20_stage05)}
- Stage06 impropers above 20 degrees: {len(above20_stage06)}
- Persistent above 20 degrees: {len(persistent_above20)}
- Newly above 20 degrees: {len(newly_above20)}
- Recovered below 20 degrees: {len(recovered_below20)}
""",
        encoding="utf-8",
    )

    print(
        "Day021 Stage05-Stage06 HBN improper "
        "persistence comparison completed."
    )

    print(
        "Stage05 planarity q99/max: "
        f"{metrics05['q99']:.4f}/"
        f"{metrics05['maximum']:.4f} deg"
    )

    print(
        "Stage06 planarity q99/max: "
        f"{metrics06['q99']:.4f}/"
        f"{metrics06['maximum']:.4f} deg"
    )

    print(
        "Planarity q99/max increase: "
        f"{summary['planarity_q99_increase_deg']:.4f}/"
        f"{summary['planarity_max_increase_deg']:.4f} deg"
    )

    print(
        "Stage06 calibrated equilibrium q99/max: "
        f"{equilibrium06['q99']:.4f}/"
        f"{equilibrium06['maximum']:.4f} deg"
    )

    print(
        "Stage05-Stage06 change q99/max: "
        f"{change_metrics['q99']:.4f}/"
        f"{change_metrics['maximum']:.4f} deg"
    )

    print(
        f"Top 1% overlap: "
        f"{len(top05 & top06)}/{top_count}"
    )

    print(
        "Above 20 deg Stage05/Stage06/persistent/new/recovered: "
        f"{len(above20_stage05)}/"
        f"{len(above20_stage06)}/"
        f"{len(persistent_above20)}/"
        f"{len(newly_above20)}/"
        f"{len(recovered_below20)}"
    )

    print(
        f"Targeted decision: {decision}"
    )

    print(
        f"Authorized next step: {next_step}"
    )

    print(
        "Long mobile production authorized: NO"
    )

    if review_reasons:
        print(
            "Review reasons: "
            + " | ".join(
                review_reasons
            )
        )

    if blocked_reasons:
        print(
            "Blocking reasons: "
            + " | ".join(
                blocked_reasons
            )
        )

    print(
        f"Wrote: "
        f"{REPORT_MD.relative_to(ROOT)}"
    )


if __name__ == "__main__":
    main()
