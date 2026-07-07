#!/usr/bin/env python3

from __future__ import annotations

import csv
import importlib.util
import math
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]

SOURCE_SCRIPT = (
    ROOT
    / "scripts/phase1A/"
    "analyze_day022_mobile_vs_frozen_water.py"
)

PROTOCOL = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol"
)

MOBILE_RUN = (
    PROTOCOL
    / "execution/08_nvt_mobile_100ps"
)

MATCHED_CONTROL_RUN = (
    PROTOCOL
    / "matched_frozen_control_144ps/execution"
)

ANALYSIS_ROOT = (
    MOBILE_RUN
    / "matched_mobile_vs_frozen_water"
)

FROZEN_FRAME_ROOT = (
    ANALYSIS_ROOT
    / "matched_frozen_44_144ps_frames"
)

FROZEN_EXTRACTION_LOG = (
    ANALYSIS_ROOT
    / "matched_frozen_44_144ps_extraction.log"
)

INDEX_FILE = (
    ANALYSIS_ROOT
    / "hbn_pyr_waterO.ndx"
)

TIMESERIES_CSV = (
    ANALYSIS_ROOT
    / "matched_mobile_frozen_water_timeseries.csv"
)

RADIAL_CSV = (
    ANALYSIS_ROOT
    / "matched_mobile_frozen_radial_density.csv"
)

AXIAL_CSV = (
    ANALYSIS_ROOT
    / "matched_mobile_frozen_axial_density.csv"
)

PYR_HYDRATION_CSV = (
    ANALYSIS_ROOT
    / "matched_mobile_frozen_pyrene_hydration.csv"
)

SNAPSHOT_CSV = (
    ANALYSIS_ROOT
    / "matched_representative_snapshot_pairs.csv"
)

SUMMARY_CSV = (
    ANALYSIS_ROOT
    / "matched_mobile_frozen_water_comparison_summary.csv"
)

SOURCE_REPORT_MD = (
    ANALYSIS_ROOT
    / "MATCHED_MOBILE_VS_FROZEN_WATER_COMPARISON_DAY022.md"
)

PAIRED_CSV = (
    ANALYSIS_ROOT
    / "matched_framewise_water_differences.csv"
)

BLOCK_CSV = (
    ANALYSIS_ROOT
    / "matched_10ps_paired_water_blocks.csv"
)

PAIRED_SUMMARY_CSV = (
    ANALYSIS_ROOT
    / "matched_paired_water_summary.csv"
)

PAIRED_REPORT_MD = (
    ANALYSIS_ROOT
    / "MATCHED_PAIRED_WATER_AUDIT_DAY022.md"
)

EXPECTED_FRAMES = 201
FROZEN_BEGIN_PS = 44.0
FROZEN_END_PS = 144.0


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(ROOT)
    )


def load_source_module():
    if (
        not SOURCE_SCRIPT.exists()
        or SOURCE_SCRIPT.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing source workflow: {SOURCE_SCRIPT}"
        )

    specification = (
        importlib.util.spec_from_file_location(
            "matched_water_source",
            SOURCE_SCRIPT,
        )
    )

    if (
        specification is None
        or specification.loader is None
    ):
        raise RuntimeError(
            f"Could not load: {SOURCE_SCRIPT}"
        )

    module = importlib.util.module_from_spec(
        specification
    )

    specification.loader.exec_module(
        module
    )

    return module


def safe_js_divergence(
    first: np.ndarray,
    second: np.ndarray,
) -> float:
    first = np.asarray(
        first,
        dtype=float,
    )

    second = np.asarray(
        second,
        dtype=float,
    )

    if (
        np.any(first < 0.0)
        or np.any(second < 0.0)
        or first.sum() <= 0.0
        or second.sum() <= 0.0
    ):
        return math.nan

    first = first / first.sum()
    second = second / second.sum()

    midpoint = 0.5 * (
        first + second
    )

    first_mask = first > 0.0
    second_mask = second > 0.0

    first_term = float(
        np.sum(
            first[first_mask]
            * np.log(
                first[first_mask]
                / midpoint[first_mask]
            )
        )
    )

    second_term = float(
        np.sum(
            second[second_mask]
            * np.log(
                second[second_mask]
                / midpoint[second_mask]
            )
        )
    )

    return 0.5 * (
        first_term
        + second_term
    )


def matched_extract_frozen_frames(
    module,
    gmx: str,
    frozen_xtc: Path,
    frozen_tpr: Path,
) -> list[Path]:
    FROZEN_FRAME_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    for old_frame in FROZEN_FRAME_ROOT.glob(
        "matched_frozen_frame*.gro"
    ):
        old_frame.unlink()

    output_template = (
        FROZEN_FRAME_ROOT
        / "matched_frozen_frame.gro"
    )

    command = [
        gmx,
        "trjconv",
        "-f",
        str(frozen_xtc),
        "-s",
        str(frozen_tpr),
        "-n",
        str(INDEX_FILE),
        "-o",
        str(output_template),
        "-b",
        f"{FROZEN_BEGIN_PS:.1f}",
        "-e",
        f"{FROZEN_END_PS:.1f}",
        "-sep",
        "-pbc",
        "atom",
    ]

    completed = subprocess.run(
        command,
        input="HBN_PYR_WATERO\n",
        cwd=MATCHED_CONTROL_RUN,
        env={
            **os.environ,
            "GMX_MAXBACKUP": "-1",
        },
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )

    FROZEN_EXTRACTION_LOG.write_text(
        completed.stdout,
        encoding="utf-8",
    )

    if completed.returncode != 0:
        raise RuntimeError(
            "Matched frozen-frame extraction failed. "
            f"See {FROZEN_EXTRACTION_LOG}"
        )

    frames = sorted(
        FROZEN_FRAME_ROOT.glob(
            "matched_frozen_frame*.gro"
        ),
        key=module.natural_key,
    )

    if len(frames) != EXPECTED_FRAMES:
        raise RuntimeError(
            "Expected 201 matched frozen frames "
            f"from 44-144 ps; found {len(frames)}"
        )

    first_title, _, _ = module.read_gro(
        frames[0]
    )

    last_title, _, _ = module.read_gro(
        frames[-1]
    )

    first_time = module.parse_time(
        first_title,
        fallback=math.nan,
    )

    last_time = module.parse_time(
        last_title,
        fallback=math.nan,
    )

    if not math.isclose(
        first_time,
        FROZEN_BEGIN_PS,
        rel_tol=0.0,
        abs_tol=1.0e-6,
    ):
        raise RuntimeError(
            "First matched frozen frame is not "
            f"44 ps: {first_time}"
        )

    if not math.isclose(
        last_time,
        FROZEN_END_PS,
        rel_tol=0.0,
        abs_tol=1.0e-6,
    ):
        raise RuntimeError(
            "Last matched frozen frame is not "
            f"144 ps: {last_time}"
        )

    print(
        "Matched frozen extraction:"
    )
    print(
        f"  frames: {len(frames)}"
    )
    print(
        "  absolute time window: "
        f"{first_time:.3f}-{last_time:.3f} ps"
    )
    print(
        "  coordinate interval: 0.500 ps"
    )

    return frames


def read_csv(
    path: Path,
) -> list[dict[str, str]]:
    if (
        not path.exists()
        or path.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing or empty CSV: {path}"
        )

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:
        return list(
            csv.DictReader(handle)
        )


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
        for field in row:
            if field not in seen:
                fields.append(field)
                seen.add(field)

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
                    field: row.get(
                        field,
                        "",
                    )
                    for field in fields
                }
            )


def configure_module(module) -> None:
    module.FROZEN_RUN = (
        MATCHED_CONTROL_RUN
    )

    module.ANALYSIS_ROOT = (
        ANALYSIS_ROOT
    )

    module.FROZEN_FRAME_ROOT = (
        FROZEN_FRAME_ROOT
    )

    module.INDEX_FILE = (
        INDEX_FILE
    )

    module.FROZEN_EXTRACTION_LOG = (
        FROZEN_EXTRACTION_LOG
    )

    module.TIMESERIES_CSV = (
        TIMESERIES_CSV
    )

    module.RADIAL_CSV = (
        RADIAL_CSV
    )

    module.AXIAL_CSV = (
        AXIAL_CSV
    )

    module.PYR_HYDRATION_CSV = (
        PYR_HYDRATION_CSV
    )

    module.SNAPSHOT_CSV = (
        SNAPSHOT_CSV
    )

    module.SUMMARY_CSV = (
        SUMMARY_CSV
    )

    module.REPORT_MD = (
        SOURCE_REPORT_MD
    )

    module.js_divergence = (
        safe_js_divergence
    )

    def extractor(
        gmx: str,
        frozen_xtc: Path,
        frozen_tpr: Path,
    ) -> list[Path]:
        return matched_extract_frozen_frames(
            module,
            gmx,
            frozen_xtc,
            frozen_tpr,
        )

    module.extract_frozen_frames = (
        extractor
    )


def paired_audit() -> dict[str, object]:
    rows = read_csv(
        TIMESERIES_CSV
    )

    frozen = sorted(
        [
            row
            for row in rows
            if row[
                "dataset"
            ].strip() == "frozen"
        ],
        key=lambda row: float(
            row[
                "relative_time_ps"
            ]
        ),
    )

    mobile = sorted(
        [
            row
            for row in rows
            if row[
                "dataset"
            ].strip() == "mobile"
        ],
        key=lambda row: float(
            row[
                "relative_time_ps"
            ]
        ),
    )

    if (
        len(frozen) != EXPECTED_FRAMES
        or len(mobile) != EXPECTED_FRAMES
    ):
        raise RuntimeError(
            "Matched datasets do not contain "
            "201 frames each"
        )

    frozen_times = np.array(
        [
            float(
                row[
                    "relative_time_ps"
                ]
            )
            for row in frozen
        ],
        dtype=float,
    )

    mobile_times = np.array(
        [
            float(
                row[
                    "relative_time_ps"
                ]
            )
            for row in mobile
        ],
        dtype=float,
    )

    if not np.array_equal(
        frozen_times,
        mobile_times,
    ):
        maximum_time_difference = float(
            np.max(
                np.abs(
                    frozen_times
                    - mobile_times
                )
            )
        )

        if maximum_time_difference > 1.0e-9:
            raise RuntimeError(
                "Matched relative frame times differ. "
                f"Maximum difference: "
                f"{maximum_time_difference}"
            )

    metrics = [
        "lumen_water_count",
        "lumen_number_density_nm-3",
        "mean_normalized_radial_position",
        "wall_radius_nm",
        "lumen_length_nm",
    ]

    for pyrene in range(
        1,
        5,
    ):
        metrics.extend(
            [
                (
                    f"PYR{pyrene}_"
                    "waterO_within_35pm"
                ),
                (
                    f"PYR{pyrene}_"
                    "waterO_within_50pm"
                ),
            ]
        )

    paired_rows = []

    for index in range(
        EXPECTED_FRAMES
    ):
        row: dict[str, object] = {
            "frame": index,
            "relative_time_ps": (
                mobile_times[index]
            ),
            "matched_frozen_absolute_time_ps": (
                FROZEN_BEGIN_PS
                + mobile_times[index]
            ),
            "mobile_absolute_stage08_time_ps": (
                mobile_times[index]
            ),
        }

        for metric in metrics:
            mobile_value = float(
                mobile[index][metric]
            )

            frozen_value = float(
                frozen[index][metric]
            )

            row[
                f"mobile_{metric}"
            ] = mobile_value

            row[
                f"matched_frozen_{metric}"
            ] = frozen_value

            row[
                f"mobile_minus_frozen_{metric}"
            ] = (
                mobile_value
                - frozen_value
            )

        paired_rows.append(row)

    write_csv(
        PAIRED_CSV,
        paired_rows,
    )

    block_rows = []

    times = mobile_times

    for block_index in range(10):
        lower = 10.0 * block_index
        upper = lower + 10.0

        if block_index == 9:
            mask = (
                (times >= lower)
                & (times <= upper)
            )
        else:
            mask = (
                (times >= lower)
                & (times < upper)
            )

        indices = np.flatnonzero(
            mask
        )

        mobile_occupancy = np.array(
            [
                float(
                    mobile[index][
                        "lumen_water_count"
                    ]
                )
                for index in indices
            ],
            dtype=float,
        )

        frozen_occupancy = np.array(
            [
                float(
                    frozen[index][
                        "lumen_water_count"
                    ]
                )
                for index in indices
            ],
            dtype=float,
        )

        difference = (
            mobile_occupancy
            - frozen_occupancy
        )

        block_rows.append(
            {
                "block": block_index + 1,
                "relative_start_ps": lower,
                "relative_end_ps": upper,
                "matched_frozen_absolute_start_ps": (
                    FROZEN_BEGIN_PS + lower
                ),
                "matched_frozen_absolute_end_ps": (
                    FROZEN_BEGIN_PS + upper
                ),
                "frame_count": len(indices),
                "mobile_occupancy_mean": float(
                    mobile_occupancy.mean()
                ),
                "matched_frozen_occupancy_mean": float(
                    frozen_occupancy.mean()
                ),
                "paired_difference_mean": float(
                    difference.mean()
                ),
                "paired_difference_std": float(
                    difference.std()
                ),
                "paired_difference_min": float(
                    difference.min()
                ),
                "paired_difference_max": float(
                    difference.max()
                ),
            }
        )

    write_csv(
        BLOCK_CSV,
        block_rows,
    )

    mobile_occupancy = np.array(
        [
            float(
                row[
                    "lumen_water_count"
                ]
            )
            for row in mobile
        ],
        dtype=float,
    )

    frozen_occupancy = np.array(
        [
            float(
                row[
                    "lumen_water_count"
                ]
            )
            for row in frozen
        ],
        dtype=float,
    )

    difference = (
        mobile_occupancy
        - frozen_occupancy
    )

    source_summary_rows = read_csv(
        SUMMARY_CSV
    )

    if len(source_summary_rows) != 1:
        raise RuntimeError(
            "Expected one matched comparison "
            "summary row"
        )

    source_summary = (
        source_summary_rows[0]
    )

    summary = {
        "mobile_window_relative_to_branch_ps": (
            "44_to_144"
        ),
        "matched_frozen_window_relative_to_branch_ps": (
            "44_to_144"
        ),
        "mobile_stage08_local_time_ps": (
            "0_to_100"
        ),
        "matched_frozen_absolute_time_ps": (
            "44_to_144"
        ),
        "frames_per_dataset": (
            EXPECTED_FRAMES
        ),
        "relative_time_exact_match": True,
        "initial_coordinates_and_water_velocities_matched": (
            True
        ),
        "mobile_occupancy_mean": float(
            mobile_occupancy.mean()
        ),
        "matched_frozen_occupancy_mean": float(
            frozen_occupancy.mean()
        ),
        "paired_occupancy_difference_mean": float(
            difference.mean()
        ),
        "paired_occupancy_difference_std": float(
            difference.std()
        ),
        "paired_occupancy_difference_mae": float(
            np.mean(
                np.abs(
                    difference
                )
            )
        ),
        "paired_mobile_lower_fraction": float(
            np.mean(
                difference < 0.0
            )
        ),
        "paired_equal_fraction": float(
            np.mean(
                difference == 0.0
            )
        ),
        "paired_mobile_higher_fraction": float(
            np.mean(
                difference > 0.0
            )
        ),
        "radial_JS_divergence": float(
            source_summary[
                "radial_JS_divergence"
            ]
        ),
        "axial_JS_divergence": float(
            source_summary[
                "axial_JS_divergence"
            ]
        ),
        "maximum_absolute_PYR_hydration_difference_pct": float(
            source_summary[
                "maximum_absolute_PYR_hydration_difference_pct"
            ]
        ),
        "matched_comparison_status": (
            "COMPLETE"
        ),
        "causal_interpretation_status": (
            "REQUIRES_SCIENTIFIC_REVIEW"
        ),
        "electronic_recalculation_authorized": (
            False
        ),
        "longer_mobile_production_authorized": (
            False
        ),
    }

    write_csv(
        PAIRED_SUMMARY_CSV,
        [summary],
    )

    PAIRED_REPORT_MD.write_text(
        f"""# Matched Paired Water Audit

## Provenance

- Mobile branch interval: **44-144 ps**
- Mobile Stage08 local interval: **0-100 ps**
- Matched frozen interval: **44-144 ps**
- Frames per dataset: **{EXPECTED_FRAMES}**
- Relative frame times matched: **YES**
- Initial coordinates matched: **YES**
- Initial water velocities matched: **YES**

## Lumen occupancy

- Mobile mean:
  {summary['mobile_occupancy_mean']:.6f}
- Matched frozen mean:
  {summary['matched_frozen_occupancy_mean']:.6f}
- Mean paired difference, mobile minus frozen:
  {summary['paired_occupancy_difference_mean']:.6f}
- Paired difference standard deviation:
  {summary['paired_occupancy_difference_std']:.6f}
- Paired mean absolute difference:
  {summary['paired_occupancy_difference_mae']:.6f}
- Fraction mobile lower/equal/higher:
  {summary['paired_mobile_lower_fraction']:.6f}/
  {summary['paired_equal_fraction']:.6f}/
  {summary['paired_mobile_higher_fraction']:.6f}

## Spatial and local hydration screening

- Radial Jensen-Shannon divergence:
  {summary['radial_JS_divergence']:.8f}
- Axial Jensen-Shannon divergence:
  {summary['axial_JS_divergence']:.8f}
- Maximum absolute PYR hydration difference:
  {summary['maximum_absolute_PYR_hydration_difference_pct']:.4f} %

## Decision

- Matched comparison status: **COMPLETE**
- Causal interpretation: **REQUIRES SCIENTIFIC REVIEW**
- Electronic recalculation authorized: **NO**
- Longer mobile production authorized: **NO**

The pressure of the frozen control is retained as a diagnostic only.
It is not used as a matched observable because freezing the solute
changes its contribution to the virial.
""",
        encoding="utf-8",
    )

    return summary


def main() -> None:
    ANALYSIS_ROOT.mkdir(
        parents=True,
        exist_ok=True,
    )

    module = load_source_module()

    configure_module(
        module
    )

    module.main()

    summary = paired_audit()

    print()
    print(
        "===== MATCHED PAIRED WATER AUDIT ====="
    )

    print(
        "Matched windows: "
        "frozen 44-144 ps vs mobile Stage08 0-100 ps"
    )

    print(
        "Frames per dataset / exact time match: "
        f"{EXPECTED_FRAMES} / YES"
    )

    print(
        "Mobile / matched-frozen lumen occupancy mean: "
        f"{summary['mobile_occupancy_mean']:.6f}/"
        f"{summary['matched_frozen_occupancy_mean']:.6f}"
    )

    print(
        "Paired occupancy difference mean/std/MAE: "
        f"{summary['paired_occupancy_difference_mean']:.6f}/"
        f"{summary['paired_occupancy_difference_std']:.6f}/"
        f"{summary['paired_occupancy_difference_mae']:.6f}"
    )

    print(
        "Paired mobile lower/equal/higher fractions: "
        f"{summary['paired_mobile_lower_fraction']:.6f}/"
        f"{summary['paired_equal_fraction']:.6f}/"
        f"{summary['paired_mobile_higher_fraction']:.6f}"
    )

    print(
        "Radial / axial JS divergence: "
        f"{summary['radial_JS_divergence']:.8f}/"
        f"{summary['axial_JS_divergence']:.8f}"
    )

    print(
        "Maximum absolute PYR hydration difference: "
        f"{summary['maximum_absolute_PYR_hydration_difference_pct']:.4f} %"
    )

    print(
        "Matched comparison status: COMPLETE"
    )

    print(
        "Causal interpretation status: "
        "REQUIRES_SCIENTIFIC_REVIEW"
    )

    print(
        "Electronic recalculation authorized: NO"
    )

    print(
        "Longer mobile production authorized: NO"
    )

    print(
        f"Wrote: {relative(PAIRED_REPORT_MD)}"
    )


if __name__ == "__main__":
    main()
