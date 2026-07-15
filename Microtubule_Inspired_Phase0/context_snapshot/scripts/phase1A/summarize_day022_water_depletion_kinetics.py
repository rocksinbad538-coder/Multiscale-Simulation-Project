#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

ANALYSIS_ROOT = (
    ROOT
    / "runs/phase1A/day021_mobile_restraint_protocol/"
    "execution/08_nvt_mobile_100ps/"
    "mobile_vs_frozen_water"
)

ENDPOINT_CSV = (
    ANALYSIS_ROOT
    / "water_depletion_endpoint_audit.csv"
)

BLOCK_CSV = (
    ANALYSIS_ROOT
    / "water_depletion_10ps_block_audit.csv"
)

REPORT_MD = (
    ANALYSIS_ROOT
    / "WATER_DEPLETION_KINETICS_SUMMARY_DAY022.md"
)


def read_csv(path: Path) -> list[dict[str, str]]:
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
        return list(
            csv.DictReader(handle)
        )


def as_float(
    row: dict[str, str],
    field: str,
) -> float:
    value = row.get(
        field,
        "",
    ).strip()

    if not value:
        raise RuntimeError(
            f"Missing numeric field {field} "
            f"for row {row}"
        )

    return float(value)


def main() -> None:
    endpoint_rows = read_csv(
        ENDPOINT_CSV
    )

    block_rows = read_csv(
        BLOCK_CSV
    )

    analyzed = [
        row
        for row in endpoint_rows
        if row.get(
            "status",
            "",
        ).strip() == "ANALYZED"
    ]

    analyzed.sort(
        key=lambda row: as_float(
            row,
            "elapsed_mobile_ps",
        )
    )

    if not analyzed:
        raise RuntimeError(
            "No analyzed endpoint rows were found"
        )

    transitions = []

    for previous, current in zip(
        analyzed[:-1],
        analyzed[1:],
    ):
        previous_count = as_float(
            previous,
            "lumen_water_count",
        )

        current_count = as_float(
            current,
            "lumen_water_count",
        )

        previous_time = as_float(
            previous,
            "elapsed_mobile_ps",
        )

        current_time = as_float(
            current,
            "elapsed_mobile_ps",
        )

        duration = (
            current_time
            - previous_time
        )

        change = (
            current_count
            - previous_count
        )

        transitions.append(
            {
                "from_stage": previous["stage"],
                "to_stage": current["stage"],
                "start_ps": previous_time,
                "end_ps": current_time,
                "duration_ps": duration,
                "start_count": previous_count,
                "end_count": current_count,
                "change": change,
                "rate_water_per_ps": (
                    change / duration
                    if duration > 0.0
                    else 0.0
                ),
            }
        )

    largest_loss = min(
        transitions,
        key=lambda row: row["change"],
    )

    zero_rows = [
        row
        for row in analyzed
        if as_float(
            row,
            "lumen_water_count",
        ) <= 0.5
    ]

    first_zero = (
        zero_rows[0]
        if zero_rows
        else None
    )

    frozen_blocks = sorted(
        (
            row
            for row in block_rows
            if row.get(
                "dataset",
                "",
            ).strip() == "frozen"
        ),
        key=lambda row: int(
            float(
                row["block"]
            )
        ),
    )

    mobile_blocks = sorted(
        (
            row
            for row in block_rows
            if row.get(
                "dataset",
                "",
            ).strip() == "mobile"
        ),
        key=lambda row: int(
            float(
                row["block"]
            )
        ),
    )

    if (
        len(frozen_blocks) != 10
        or len(mobile_blocks) != 10
    ):
        raise RuntimeError(
            "Expected ten frozen and ten mobile blocks"
        )

    lines = [
        "# Water-Depletion Kinetics Summary",
        "",
        "## Endpoint progression",
        "",
        "| Stage | Elapsed mobile time (ps) | Lumen water | Density (nm^-3) | Radius (nm) | Length (nm) |",
        "|---|---:|---:|---:|---:|---:|",
    ]

    print(
        "===== ENDPOINT WATER OCCUPANCY ====="
    )

    print(
        f"{'STAGE':34s} "
        f"{'TIME_PS':>8s} "
        f"{'WATER':>8s} "
        f"{'DENSITY':>12s} "
        f"{'RADIUS':>10s} "
        f"{'LENGTH':>10s}"
    )

    for row in analyzed:
        stage = row["stage"]
        time_ps = as_float(
            row,
            "elapsed_mobile_ps",
        )

        count = as_float(
            row,
            "lumen_water_count",
        )

        density = as_float(
            row,
            "lumen_number_density_nm-3",
        )

        radius = as_float(
            row,
            "wall_radius_nm",
        )

        length = as_float(
            row,
            "lumen_length_nm",
        )

        print(
            f"{stage:34s} "
            f"{time_ps:8.1f} "
            f"{count:8.0f} "
            f"{density:12.6f} "
            f"{radius:10.6f} "
            f"{length:10.6f}"
        )

        lines.append(
            f"| {stage} | {time_ps:.1f} | "
            f"{count:.0f} | {density:.6f} | "
            f"{radius:.6f} | {length:.6f} |"
        )

    lines.extend(
        [
            "",
            "## Endpoint transitions",
            "",
            "| From | To | Interval (ps) | Water change | Rate (water/ps) |",
            "|---|---|---:|---:|---:|",
        ]
    )

    print()
    print(
        "===== ENDPOINT TRANSITIONS ====="
    )

    print(
        f"{'FROM -> TO':64s} "
        f"{'INTERVAL':>12s} "
        f"{'CHANGE':>10s} "
        f"{'RATE':>12s}"
    )

    for row in transitions:
        label = (
            f"{row['from_stage']} -> "
            f"{row['to_stage']}"
        )

        interval = (
            f"{row['start_ps']:.1f}-"
            f"{row['end_ps']:.1f}"
        )

        print(
            f"{label:64s} "
            f"{interval:>12s} "
            f"{row['change']:10.0f} "
            f"{row['rate_water_per_ps']:12.6f}"
        )

        lines.append(
            f"| {row['from_stage']} | "
            f"{row['to_stage']} | "
            f"{interval} | "
            f"{row['change']:.0f} | "
            f"{row['rate_water_per_ps']:.6f} |"
        )

    lines.extend(
        [
            "",
            "## Ten-picosecond blocks",
            "",
            "| Dataset | Block | Interval (ps) | Mean | Std | Min | Max |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )

    print()
    print(
        "===== TEN-PS BLOCK OCCUPANCY ====="
    )

    print(
        f"{'DATASET':10s} "
        f"{'BLOCK':>6s} "
        f"{'INTERVAL':>12s} "
        f"{'MEAN':>10s} "
        f"{'STD':>10s} "
        f"{'MIN':>8s} "
        f"{'MAX':>8s}"
    )

    for dataset, rows in (
        ("frozen", frozen_blocks),
        ("mobile", mobile_blocks),
    ):
        for row in rows:
            block = int(
                float(
                    row["block"]
                )
            )

            start = as_float(
                row,
                "start_ps",
            )

            end = as_float(
                row,
                "end_ps",
            )

            mean = as_float(
                row,
                "occupancy_mean",
            )

            std = as_float(
                row,
                "occupancy_std",
            )

            minimum = as_float(
                row,
                "occupancy_min",
            )

            maximum = as_float(
                row,
                "occupancy_max",
            )

            interval = (
                f"{start:.0f}-{end:.0f}"
            )

            print(
                f"{dataset:10s} "
                f"{block:6d} "
                f"{interval:>12s} "
                f"{mean:10.4f} "
                f"{std:10.4f} "
                f"{minimum:8.0f} "
                f"{maximum:8.0f}"
            )

            lines.append(
                f"| {dataset} | {block} | "
                f"{interval} | {mean:.4f} | "
                f"{std:.4f} | {minimum:.0f} | "
                f"{maximum:.0f} |"
            )

    first_zero_description = (
        (
            f"{first_zero['stage']} at "
            f"{as_float(first_zero, 'elapsed_mobile_ps'):.1f} ps"
        )
        if first_zero is not None
        else "not reached"
    )

    lines.extend(
        [
            "",
            "## Diagnostic findings",
            "",
            (
                "- Largest endpoint loss: "
                f"`{largest_loss['from_stage']}` → "
                f"`{largest_loss['to_stage']}`, "
                f"{largest_loss['change']:.0f} waters over "
                f"{largest_loss['duration_ps']:.1f} ps."
            ),
            (
                "- First endpoint with zero lumen waters: "
                f"**{first_zero_description}**."
            ),
            (
                "- Existing frozen and mobile windows remain "
                "temporally unmatched."
            ),
            (
                "- No electronic recalculation is authorized "
                "from the unmatched comparison."
            ),
        ]
    )

    REPORT_MD.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )

    print()
    print(
        "===== KINETIC SUMMARY ====="
    )

    print(
        "Largest endpoint loss: "
        f"{largest_loss['from_stage']} -> "
        f"{largest_loss['to_stage']}: "
        f"{largest_loss['change']:.0f} waters over "
        f"{largest_loss['duration_ps']:.1f} ps"
    )

    print(
        "First zero-occupancy endpoint: "
        f"{first_zero_description}"
    )

    print(
        "Matched frozen continuation required: YES"
    )

    print(
        "Electronic recalculation authorized: NO"
    )

    print(
        "Wrote: "
        + str(
            REPORT_MD.relative_to(ROOT)
        )
    )


if __name__ == "__main__":
    main()
