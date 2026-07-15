#!/usr/bin/env python3

from __future__ import annotations

import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

RUN_ROOT = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol/"
    "execution/08_nvt_mobile_100ps"
)

OPERATIONAL_CSV = (
    RUN_ROOT
    / "08_nvt_mobile_100ps_summary.csv"
)

STRUCTURAL_CSV = (
    RUN_ROOT
    / "08_nvt_mobile_100ps_structural_summary.csv"
)

TEMPORAL_CSV = (
    RUN_ROOT
    / "time_resolved_stability/"
    "stage08_time_resolved_stability_summary.csv"
)

COMPARISON_CSV = (
    RUN_ROOT
    / "rmsf_tube_frozen_comparison/"
    "stage08_rmsf_tube_frozen_comparison_summary.csv"
)

XTC_CHECK_LOG = (
    RUN_ROOT
    / "run_stage08_xtc_check.log"
)

OUTPUT_CSV = (
    RUN_ROOT
    / "stage08_mobile_100ps_acceptance_summary.csv"
)

OUTPUT_MD = (
    RUN_ROOT
    / "STAGE08_MOBILE_100PS_ACCEPTANCE_DAY022.md"
)


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
    )


def read_single_row(
    path: Path,
) -> dict[str, str]:
    if (
        not path.exists()
        or path.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing or empty file: {path}"
        )

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        rows = list(csv.DictReader(handle))

    if len(rows) != 1:
        raise RuntimeError(
            f"Expected exactly one row in {path}; "
            f"found {len(rows)}"
        )

    return rows[0]


def number(
    row: dict[str, str],
    field: str,
) -> float:
    try:
        return float(row[field])
    except (KeyError, ValueError) as error:
        raise RuntimeError(
            f"Invalid numeric field: {field}"
        ) from error


def parse_xtc_check() -> tuple[int, float]:
    if (
        not XTC_CHECK_LOG.exists()
        or XTC_CHECK_LOG.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing XTC check log: {XTC_CHECK_LOG}"
        )

    text = XTC_CHECK_LOG.read_text(
        encoding="utf-8",
        errors="replace",
    )

    frame_match = re.search(
        r"Coords\s+(\d+)\s+",
        text,
    )

    time_match = re.search(
        r"Last frame\s+\d+\s+time\s+"
        r"([-+0-9.eE]+)",
        text,
    )

    if (
        frame_match is None
        or time_match is None
    ):
        raise RuntimeError(
            "Could not parse Stage08 XTC check"
        )

    return (
        int(frame_match.group(1)),
        float(time_match.group(1)),
    )


def write_csv(
    row: dict[str, object],
) -> None:
    with OUTPUT_CSV.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(row.keys()),
        )

        writer.writeheader()
        writer.writerow(row)


def main() -> None:
    operational = read_single_row(
        OPERATIONAL_CSV
    )

    structural = read_single_row(
        STRUCTURAL_CSV
    )

    temporal = read_single_row(
        TEMPORAL_CSV
    )

    comparison = read_single_row(
        COMPARISON_CSV
    )

    frames, duration_ps = parse_xtc_check()

    failures: list[str] = []

    if (
        operational.get(
            "decision",
            "",
        ).strip().upper()
        != "PASS"
    ):
        failures.append(
            "operational decision is not PASS"
        )

    if (
        structural.get(
            "structural_screen",
            "",
        ).strip().upper()
        != "STABLE_CANDIDATE"
    ):
        failures.append(
            "structural screen is not STABLE_CANDIDATE"
        )

    if (
        temporal.get(
            "time_resolved_decision",
            "",
        ).strip().upper()
        != "PASS_WITH_MONITORING"
    ):
        failures.append(
            "time-resolved decision is not "
            "PASS_WITH_MONITORING"
        )

    if (
        comparison.get(
            "comparison_decision",
            "",
        ).strip().upper()
        != "PASS_WITH_MONITORING"
    ):
        failures.append(
            "RMSF/tube comparison is not "
            "PASS_WITH_MONITORING"
        )

    if temporal.get(
        "blocked_reasons",
        "",
    ).strip():
        failures.append(
            "time-resolved analysis contains blocking reasons"
        )

    if comparison.get(
        "blocked_reasons",
        "",
    ).strip():
        failures.append(
            "RMSF/tube comparison contains blocking reasons"
        )

    if frames != 201:
        failures.append(
            f"expected 201 frames; found {frames}"
        )

    if abs(
        duration_ps - 100.0
    ) > 1.0e-6:
        failures.append(
            f"expected 100 ps; found {duration_ps}"
        )

    quantitative_checks = {
        "HBN RMSF q95 <= 0.15 nm": (
            number(
                comparison,
                "mobile_HBN_rmsf_q95_nm",
            )
            <= 0.15
        ),
        "HBN RMSF maximum <= 0.35 nm": (
            number(
                comparison,
                "mobile_HBN_rmsf_max_nm",
            )
            <= 0.35
        ),
        "radius difference <= 10 percent": (
            abs(
                number(
                    comparison,
                    "tube_radius_mobile_vs_frozen_difference_pct",
                )
            )
            <= 10.0
        ),
        "length difference <= 10 percent": (
            abs(
                number(
                    comparison,
                    "tube_length_mobile_vs_frozen_difference_pct",
                )
            )
            <= 10.0
        ),
        "mean ellipticity <= 1.20": (
            number(
                comparison,
                "mobile_tube_ellipticity_mean",
            )
            <= 1.20
        ),
        "maximum ellipticity <= 1.35": (
            number(
                comparison,
                "mobile_tube_ellipticity_max",
            )
            <= 1.35
        ),
        "PYR internal RMSF <= 0.08 nm": (
            number(
                comparison,
                "mobile_PYR_max_internal_rmsf_nm",
            )
            <= 0.08
        ),
        "no persistent impropers in >=80 percent frames": (
            int(
                number(
                    temporal,
                    "HBN_impropers_above20_persistent80_count",
                )
            )
            == 0
        ),
        "minimum intergroup contact >= 0.14 nm": (
            number(
                temporal,
                "minimum_intergroup_contact_nm",
            )
            >= 0.14
        ),
    }

    failures.extend(
        label
        for label, passed
        in quantitative_checks.items()
        if not passed
    )

    if failures:
        raise RuntimeError(
            "Stage08 scientific acceptance failed:\n"
            + "\n".join(failures)
        )

    next_step = (
        "MOBILE_VS_FROZEN_WATER_AND_"
        "SOLVENT_DISORDER_COMPARISON"
    )

    row = {
        "stage": "08_nvt_mobile_100ps",
        "ensemble": "NVT",
        "temperature_K": 300,
        "duration_ps": duration_ps,
        "trajectory_frames": frames,
        "operational_decision": "PASS",
        "structural_screen": "STABLE_CANDIDATE",
        "time_resolved_decision": (
            "PASS_WITH_MONITORING"
        ),
        "rmsf_tube_comparison_decision": (
            "PASS_WITH_MONITORING"
        ),
        "final_scientific_decision": (
            "PASS_WITH_MONITORING"
        ),
        "mobile_validation_300K_closed": True,
        "authorized_next_step": next_step,
        "longer_mobile_production_authorized": False,
        "multitemperature_production_authorized": False,
        "mobile_HBN_rmsf_mean_nm": number(
            comparison,
            "mobile_HBN_rmsf_mean_nm",
        ),
        "mobile_HBN_rmsf_q95_nm": number(
            comparison,
            "mobile_HBN_rmsf_q95_nm",
        ),
        "mobile_HBN_rmsf_max_nm": number(
            comparison,
            "mobile_HBN_rmsf_max_nm",
        ),
        "tube_radius_difference_pct": number(
            comparison,
            "tube_radius_mobile_vs_frozen_difference_pct",
        ),
        "tube_length_difference_pct": number(
            comparison,
            "tube_length_mobile_vs_frozen_difference_pct",
        ),
        "tube_ellipticity_mean": number(
            comparison,
            "mobile_tube_ellipticity_mean",
        ),
        "tube_ellipticity_max": number(
            comparison,
            "mobile_tube_ellipticity_max",
        ),
        "PYR_max_internal_rmsf_nm": number(
            comparison,
            "mobile_PYR_max_internal_rmsf_nm",
        ),
        "minimum_intergroup_contact_nm": number(
            temporal,
            "minimum_intergroup_contact_nm",
        ),
        "persistent80_improper_count": int(
            number(
                temporal,
                "HBN_impropers_above20_persistent80_count",
            )
        ),
    }

    write_csv(row)

    OUTPUT_MD.write_text(
        f"""# Stage08 Mobile 100 ps Acceptance

## Decision

- Final scientific decision: **PASS_WITH_MONITORING**
- Mobile validation at 300 K closed: **YES**
- Authorized next step: `{next_step}`
- Longer mobile production authorized: **NO**
- Multitemperature production authorized: **NO**

## Trajectory

- Ensemble: NVT
- Temperature: 300 K
- Duration: {duration_ps:.1f} ps
- Frames: {frames}
- Coordinate interval: 0.5 ps

## HBN stability

- RMSF mean/q95/max: {row['mobile_HBN_rmsf_mean_nm']:.6f}/
  {row['mobile_HBN_rmsf_q95_nm']:.6f}/
  {row['mobile_HBN_rmsf_max_nm']:.6f} nm
- Radius difference versus frozen control:
  {row['tube_radius_difference_pct']:.4f} %
- Length difference versus frozen control:
  {row['tube_length_difference_pct']:.4f} %
- Ellipticity mean/max:
  {row['tube_ellipticity_mean']:.6f}/
  {row['tube_ellipticity_max']:.6f}
- Improper angles above 20 degrees in at least
  80 percent of frames: {row['persistent80_improper_count']}

## Pyrene and contacts

- Maximum internal PYR RMSF:
  {row['PYR_max_internal_rmsf_nm']:.6f} nm
- Minimum intergroup contact:
  {row['minimum_intergroup_contact_nm']:.6f} nm

The 100 ps mobile trajectory is accepted for comparative
water-organization and solvent-disorder analysis. The frozen
trajectory remains a geometric and solvent-reference control,
not evidence of solute stability.
""",
        encoding="utf-8",
    )

    print(
        "Day022 Stage08 mobile-100ps reconciliation completed."
    )
    print(
        f"Trajectory frames / duration: "
        f"{frames} / {duration_ps:.1f} ps"
    )
    print(
        "Operational decision: PASS"
    )
    print(
        "Structural screen: STABLE_CANDIDATE"
    )
    print(
        "Time-resolved decision: PASS_WITH_MONITORING"
    )
    print(
        "RMSF/tube comparison: PASS_WITH_MONITORING"
    )
    print(
        "Final scientific decision: PASS_WITH_MONITORING"
    )
    print(
        "Mobile validation at 300 K closed: YES"
    )
    print(
        f"Authorized next step: {next_step}"
    )
    print(
        "Longer mobile production authorized: NO"
    )
    print(
        "Multitemperature production authorized: NO"
    )
    print(
        f"Wrote: {relative(OUTPUT_MD)}"
    )


if __name__ == "__main__":
    main()
